package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.ODEIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.util.Arrays;

public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The equilibrium frequencies for each type",
            new RealParameter("1.0"));

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Attribute key used to specify sample trait values in tree.");

    public Input<Double> relativeToleranceInput = new Input<>("relTolerance",
            "Relative tolerance for numerical integration.",
            1e-7);

    public Input<Double> absoluteToleranceInput = new Input<>("absTolerance",
            "Absolute tolerance for numerical integration.",
            1e-100 /*Double.MIN_VALUE*/);

    private Parameterization parameterization;
    private double finalSampleOffset;
    private TreeInterface tree;
    double absoluteTolerance;
    double relativeTolerance;

    @Override
    public void initAndValidate() {
        // unpack input values

        this.parameterization = this.parameterizationInput.get();
        this.finalSampleOffset = this.finalSampleOffsetInput.get().getArrayValue();
        this.tree = this.treeInput.get();
        this.absoluteTolerance = this.absoluteToleranceInput.get();
        this.relativeTolerance = this.relativeToleranceInput.get();

        // validate typeLabel

        if (this.parameterization.getNTypes() != 1 && this.typeLabelInput.get() == null)
            throw new RuntimeException("Error: For models with >1 type, typeLabel must be specified.");

        // validate frequenciesInput

        if (this.frequenciesInput.get().getDimension() != this.parameterization.getNTypes())
            throw new RuntimeException("Error: dimension of equilibrium frequencies " +
                    "parameter must match number of types.");

        double freqSum = 0;
        for (double f : this.frequenciesInput.get().getValues()) freqSum += f;
        if (Math.abs(1.0 - freqSum) > 1e-10)
            throw new RuntimeException(
                    "Error: equilibrium frequencies must add " + "up to 1 but currently add to " + freqSum + "."
            );
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {
        ContinuousOutputModel flow = this.calculateFlow();

        Node root = this.tree.getRoot();
        double[] rootLikelihoodPerState = this.calculateSubTreeLikelihood(root, flow);

        double treeLikelihood = rootLikelihoodPerState[0];

        return treeLikelihood;
    }

    private ContinuousOutputModel calculateFlow() {
        // set up ODE

        FlowODESystem flowODE = new FlowODESystem(this.parameterization);

        double integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = parameterization.getTotalProcessLength() / 10;
        FirstOrderIntegrator flowIntegrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, this.absoluteTolerance, this.relativeTolerance
        );

        ContinuousOutputModel flow = new ContinuousOutputModel();
        flowIntegrator.addStepHandler(flow);

        // run integration over the entire timespan (all the intervals) to calculate the flow

        double[] intervalEndTimes = this.parameterization.getIntervalEndTimes();

        // TODO: initial state?
        double[] state = new double[]{};

        for (int interval = 0; interval < this.parameterization.getTotalIntervalCount(); interval++) {
            flowODE.setInterval(interval);

            double startTime = interval == 0 ? 0 : intervalEndTimes[interval - 1];
            double endTime = intervalEndTimes[interval];

            flowIntegrator.integrate(flowODE, startTime, state, endTime, state);
        }

        return flow;
    }

    private double[] calculateSubTreeLikelihood(Node root, double timeEdgeStart, double timeEdgeEnd, ContinuousOutputModel flow) {
        if (root.isFake()) {
            return new double[]{};
        }

        if (root.isLeaf()) {
            return new double[]{};
        }

        // normal edge

        Node child1 = root.getChild(0);
        double[] likelihoodChild1 = this.calculateSubTreeLikelihood(
                child1,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child1, this.finalSampleOffset),
                flow
        );

        Node child2 = root.getChild(1);
        double[] likelihoodChild2 = this.calculateSubTreeLikelihood(
                child2,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child1, this.finalSampleOffset),
                flow
        );

        // combine the child likelihoods to get the likelihood at the edge end

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);
        double[] birthRatesEdgeEnd = this.parameterization.getBirthRates()[intervalEdgeEnd];

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            likelihoodEdgeEnd[i] = likelihoodChild1[i];

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                likelihoodEdgeEnd[i + this.parameterization.getNTypes()] += (
                        birthRatesEdgeEnd[i]
                                * likelihoodChild1[i + this.parameterization.getNTypes()]
                                * likelihoodChild2[j + this.parameterization.getNTypes()]
                );
                likelihoodEdgeEnd[i + this.parameterization.getNTypes()] += (
                        birthRatesEdgeEnd[i]
                                * likelihoodChild1[j + this.parameterization.getNTypes()]
                                * likelihoodChild2[i + this.parameterization.getNTypes()]
                );
            }
        }

        double[] likelihoodEdgeStart;

        DecompositionSolver linearSolver = new QRDecomposition()

        flow.setInterpolatedTime(timeEdgeStart);
        double[] flowEdgeStart = flow.getInterpolatedState();

        flow.setInterpolatedTime(timeEdgeEnd);
        double[] flowEdgeEnd = flow.getInterpolatedState();


        return new double[]{};
    }
}
