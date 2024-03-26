package bdmmprime.flow;

import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;

public class FlowTest {

    static void assertMatrixEquality(RealMatrix matrix1, RealMatrix matrix2) {
        for (int i = 0; i < matrix1.getRowDimension(); i++) {
            for (int j = 0; j < matrix1.getColumnDimension(); j++) {
                assertEquals(matrix1.getEntry(i, j), matrix2.getEntry(i, j), 1e-5);
            }
        }
    }

    @Test
    public void testFlow() {
        Tree tree = new TreeParser("((3[&type=2]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(3),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("6 2 5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.45 0.55")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.333333 0.45")),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1 0.2 0.15 0.12 0.12 0.15")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 3));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter((1.0/3.0) + " " + (1.0/3.0) + " " + (1.0/3.0)),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        // setup extinction ODE

        double[] initialState = new double[parameterization.getNTypes()];
        Arrays.fill(initialState, 1.0);

        IntervalODESystem extinctionSystem = new ExtinctionODESystem(parameterization);
        ContinuousOutputModel extinctionProbabilities = extinctionSystem.integrateBackwardsOverIntegrals(
                initialState.clone(), 1e-100, 1e-20
        );

        // set up flow system

        FlowODESystem system = new FlowODESystem(
                parameterization,
                new ExtinctionProbabilities(new ContinuousOutputModel[] {extinctionProbabilities}, new double[] {parameterization.getTotalProcessLength()})
        );

        // set up normal flow

        initialState = new double[parameterization.getNTypes() * parameterization.getNTypes()];
        for (int j = 0; j < parameterization.getNTypes(); j++) {
            // fill diagonal entries with 1 to get the identity matrix
            initialState[j * parameterization.getNTypes() + j] = 1;
        }
        ContinuousOutputModel output = system.integrateOverIntegrals(
                initialState, 1e-100, 1e-10
        );
        Flow normalFlow = new Flow(
                new ContinuousOutputModel[]{output}, new double[]{parameterization.getTotalProcessLength()}, parameterization.getNTypes()
        );

        // set up interval flow

        int numIntervals = 3;

        double intervalSize = parameterization.getTotalProcessLength() / numIntervals;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[numIntervals];
        double[] endTimes = new double[numIntervals];

        for (int i = 0; i < numIntervals; i++) {
            initialState = new double[parameterization.getNTypes() * parameterization.getNTypes()];
            for (int j = 0; j < parameterization.getNTypes(); j++) {
                // fill diagonal entries with 1 to get the identity matrix
                initialState[j * parameterization.getNTypes() + j] = 1;
            }

            double startTime = i * intervalSize;
            double endTime = (i+1) * intervalSize;

            output = system.integrateOverIntegrals(
                    startTime, endTime, initialState, 1e-100, 1e-10
            );
            outputModels[i] = output;
            endTimes[i] = endTime;
        }

        Flow intervalFlow = new Flow(outputModels, endTimes, parameterization.getNTypes());

        assertMatrixEquality(normalFlow.getFlow(0.0), intervalFlow.getFlow(0.0));
        assertMatrixEquality(normalFlow.getFlow(1.0), intervalFlow.getFlow(1.0));
        assertMatrixEquality(normalFlow.getFlow(2.0), intervalFlow.getFlow(2.0));
        assertMatrixEquality(normalFlow.getFlow(3.0), intervalFlow.getFlow(3.0));
        assertMatrixEquality(normalFlow.getFlow(4.0), intervalFlow.getFlow(4.0));
    }
}
