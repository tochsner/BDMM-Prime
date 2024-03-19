package bdmmprime.flow;


import bdmmprime.distribution.*;
import bdmmprime.parameterization.*;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;
import static org.junit.jupiter.api.Assertions.*;

public class FlowODESystemTest {

    static class NormalFlowODESystem extends FlowODESystem {
        public NormalFlowODESystem(Parameterization parameterization, ContinuousOutputModel extinctionProbabilities) {
            super(parameterization, extinctionProbabilities);
        }

        @Override
        public int getDimension() {
            return param.getNTypes();
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {
            RealMatrix systemMatrix = this.buildSystemMatrix(t);
            RealVector yDotMatrix = systemMatrix.operate(Utils.toVector(y));
            Utils.fillVector(yDotMatrix, yDot);
        }
    }

    private void testFlowODESolver(
            Parameterization parameterization,
            double startTime,
            double endTime,
            double[] initialState
    ) {
        // setup extinction ODE

        IntervalODESystem extinctionSystem = new ExtinctionODESystem(parameterization);

        double[] initialExtinctionState = new double[parameterization.getNTypes()];
        Arrays.fill(initialExtinctionState, 1.0);

        ContinuousOutputModel extinctionProbabilities = extinctionSystem.integrateBackwardsOverIntegrals(
                initialExtinctionState, 1e-100, 1e-20
        );

        // integrate using flow

        FlowODESystem flowSystem = new FlowODESystem(parameterization, extinctionProbabilities);

        double[] initialFlowState = new double[parameterization.getNTypes() * parameterization.getNTypes()];
        for (int i = 0; i < parameterization.getNTypes(); i++) {
            // fill diagonal entries with 1 to get the identity matrix
            initialFlowState[i * parameterization.getNTypes() + i] = 1.0;
        }

        ContinuousOutputModel flow = flowSystem.integrateOverIntegrals(
                initialFlowState, 1e-100, 1e-20
        );
        double[] flowIntegral = FlowODESystem.integrateUsingFlow(startTime, endTime, initialState.clone(), flow);

        // integrate using normal integration method

        NormalFlowODESystem normalSystem = new NormalFlowODESystem(parameterization, extinctionProbabilities);
        ContinuousOutputModel normalOutput = normalSystem.integrateOverIntegrals(
                startTime, endTime, initialState.clone(), 1e-100, 1e-20
        );

        normalOutput.setInterpolatedTime(endTime);
        double[] normalIntegral = normalOutput.getInterpolatedState();

        // integrate using the old non-flow system

        P0GeSystem oldSystem = new P0GeSystem(
                parameterization, 1e-100, 1e-20
        );

        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                parameterization.getTotalProcessLength() * 1e-100,
                parameterization.getTotalProcessLength() / 10, 1e-100, 1e-10
        );
        extinctionProbabilities.setInterpolatedTime(startTime);
        double[] oldIntegral = new double[]{
                extinctionProbabilities.getInterpolatedState()[0],
                extinctionProbabilities.getInterpolatedState()[1],
                initialState[0],
                initialState[1]
        };
        integrator.integrate(oldSystem, startTime, oldIntegral.clone(), endTime, oldIntegral);
        oldIntegral = Arrays.copyOfRange(oldIntegral, 2, 4);

        System.out.println(Arrays.toString(oldIntegral));
        System.out.println(Arrays.toString(normalIntegral));
        System.out.println(Arrays.toString(flowIntegral));

        for (int i = 0; i < flowIntegral.length; i++) {
            assertEquals(oldIntegral[i], normalIntegral[i], 1e-10);
            assertEquals(oldIntegral[i], flowIntegral[i], 1e-10);
        }
    }

    void testDerivatives(
            Parameterization parameterization,
            double t,
            double[] initialState
    ) {
        // setup extinction ODE

        IntervalODESystem extinctionSystem = new ExtinctionODESystem(parameterization);

        double[] initialExtinctionState = new double[parameterization.getNTypes()];
        Arrays.fill(initialExtinctionState, 1.0);

        ContinuousOutputModel extinctionProbabilities = extinctionSystem.integrateBackwardsOverIntegrals(
                initialExtinctionState, 1e-100, 1e-20
        );

        // get flow ODE derivatives

        FlowODESystem flowSystem = new NormalFlowODESystem(parameterization, extinctionProbabilities);
        double[] flowDerivatives = new double[initialState.length];
        flowSystem.computeDerivatives(t, initialState.clone(), flowDerivatives);

        // get old ODE derivatives

        P0GeSystem oldSystem = new P0GeSystem(parameterization, 1e-100, 1e-20);

        extinctionProbabilities.setInterpolatedTime(t);
        double[] startExtinction = extinctionProbabilities.getInterpolatedState();
        double[] oldInitialState = new double[]{
                startExtinction[0],
                startExtinction[1],
                initialState[0],
                initialState[1]
        };
        double[] oldDerivatives = new double[2*initialState.length];
        oldSystem.computeDerivatives(t, oldInitialState, oldDerivatives);

        for (int i = 0; i < initialState.length; i++) {
            assertEquals(oldDerivatives[i+2], flowDerivatives[i], 1e-7);
        }
    }

    @Test
    public void testFlowODESolver2() {
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealParameter("6"),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(4.0 / 3.0) + " " + String.valueOf(4.0 / 3.0))),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0"),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.00001 0.00001"),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(1.0 / 3.0) + " " + String.valueOf(1.0 / 3.0))),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3 0.7")));

        double[] initialState = new double[]{6.52173913043E-1, 0.0};
        this.testFlowODESolver(
                parameterization, 4.5, 6, initialState
        );
    }

    @Test
    public void testFlowODESolver3() {
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("2.5"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(40.0 / 3.0) + " " + Double.toString(0.0 / 3.0)), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(1.0 / 3.0)), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1"),
                        2)
        );

        double[] initialState = new double[]{0.5, 0};

        this.testDerivatives(
                parameterization, 1, new double[]{0.5, 0.}
        );
        this.testDerivatives(
                parameterization, 0.25, new double[]{0.18, 0.01}
        );
        this.testDerivatives(
                parameterization, 0.5, new double[]{0.1792880448169992, 0.01028044974753326}
        );

        this.testFlowODESolver(
                parameterization,
                0, 2.5, initialState
        );
    }

    @Test
    public void testFlowODESolver5() {
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("2.5"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(40.0 / 3.0) + " " + Double.toString(1.0 / 3.0)), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(1.0 / 3.0)), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1"),
                        2)
        );

        testDerivatives(
                parameterization,
                0.0062384240274003236,
                new double[]{0.4425952958374468, 2.922586221786918E-4}
        );
    }

}