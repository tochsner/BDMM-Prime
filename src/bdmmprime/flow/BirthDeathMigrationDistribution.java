package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

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
        ContinuousOutputModel extinctionProbabilities = this.calculateExtinctionProbabilities();
        ContinuousOutputModel flow = this.calculateFlow(extinctionProbabilities);

        Node root = this.tree.getRoot();
        double[] rootLikelihoodPerState = this.calculateSubTreeLikelihood(
                root,
                0,
                this.parameterization.getNodeTime(tree.getRoot(), this.finalSampleOffset),
                flow
        );

        double treeLikelihood = rootLikelihoodPerState[0];

        return treeLikelihood;
    }

    private ContinuousOutputModel calculateExtinctionProbabilities() {
        IntervalODESystem system = new ExtinctionODESystem(this.parameterization);
        return system.integrateOverIntegrals(this.absoluteTolerance, this.relativeTolerance);
    }

    private ContinuousOutputModel calculateFlow(ContinuousOutputModel extinctionProbabilities) {
        FlowODESystem system = new FlowODESystem(this.parameterization, extinctionProbabilities);
        return system.integrateOverIntegrals(this.absoluteTolerance, this.relativeTolerance);
    }

    private double[] calculateSubTreeLikelihood(Node root, double timeEdgeStart, double timeEdgeEnd, ContinuousOutputModel flow) {
        if (root.isFake()) {
            throw new UnsupportedOperationException();
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
                this.parameterization.getNodeTime(child2, this.finalSampleOffset),
                flow
        );

        // combine the child likelihoods to get the likelihood at the edge end

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);
        double[][] birthRatesEdgeEnd = this.parameterization.getCrossBirthRates()[intervalEdgeEnd];

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            for (int j = 0; j < parameterization.getNTypes(); j++) {
                likelihoodEdgeEnd[i] += birthRatesEdgeEnd[i][j] * (
                        likelihoodChild1[i] * likelihoodChild2[j] + likelihoodChild1[j] * likelihoodChild2[i]
                );
            }
        }

        // solve a linear system instead of integrating

        RealVector likelihoodVectorEdgeEnd = Utils.toVector(likelihoodEdgeEnd);

        flow.setInterpolatedTime(timeEdgeStart);
        double[] flowEdgeStart = flow.getInterpolatedState();
        RealMatrix flowMatrixEdgeStart = Utils.toMatrix(flowEdgeStart, this.parameterization.getNTypes());

        DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEdgeStart).getSolver();
        RealVector solution = linearSolver.solve(likelihoodVectorEdgeEnd);

        RealVector likelihoodVectorEdgeStart = flowMatrixEdgeStart.operate(solution);
        return likelihoodVectorEdgeStart.toArray();
    }
}
