package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import static jeigen.Shortcuts.*;
import jeigen.DenseMatrix;

public class FlowODESystem extends IntervalODESystem {
    private final ExtinctionProbabilities extinctionProbabilities;
    private DenseMatrix constantMatrix;
    private RealMatrix realConstantMatrix;

    public FlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
            double absoluteTolerance,
            double relativeTolerance
    ) {
        super(parameterization, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;
    }

    @Override
    public int getDimension() {
        return param.getNTypes() * param.getNTypes();
    }

    protected RealMatrix buildTimeInvariantSystemMatrix(double t) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        int interval = this.param.getIntervalIndex(t);

        // fill transitions

        for (int i = 0; i < param.getNTypes(); i++) {
            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        j,
                        -this.param.getMigRates()[interval][i][j]
                );
            }
        }

        // fill diagonals

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    this.param.getDeathRates()[interval][i] + this.param.getSamplingRates()[interval][i]
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.param.getMigRates()[interval][i][j]
                );
            }
        }

        return system;
    }

    protected RealMatrix buildTimeVaryingSystemMatrix(double t) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        int interval = this.param.getIntervalIndex(t);
        double[] extinctProbabilities = this.extinctionProbabilities.getProbabilities(t, interval);

        // fill transitions

        for (int i = 0; i < param.getNTypes(); i++) {
            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        j,
                        -this.param.getCrossBirthRates()[interval][i][j] * extinctProbabilities[i]
                );
            }
        }

        // fill diagonals

        for (int i = 0; i < param.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    -this.param.getBirthRates()[interval][i] * (2 * extinctProbabilities[i] - 1)
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        -this.param.getCrossBirthRates()[interval][i][j] * (extinctProbabilities[j] - 1)
                );
            }
        }

        return system;
    }

    protected RealMatrix buildSystemMatrix(double t) {
        RealMatrix timeInvariantMatrix = this.buildTimeInvariantSystemMatrix(t);
        RealMatrix timeVaryingMatrix = this.buildTimeVaryingSystemMatrix(t);

        double[] values = new double[this.param.getNTypes() * this.param.getNTypes()];
        Utils.fillArray(timeVaryingMatrix, values);

        DenseMatrix A = this.constantMatrix.mul(this.param.getTotalProcessLength() - t).mexp();
        DenseMatrix B = this.constantMatrix.mul(-this.param.getTotalProcessLength() + t).mexp();

        return Utils.toMatrix(A).multiply(timeVaryingMatrix.add(timeInvariantMatrix).subtract(this.realConstantMatrix)).multiply(Utils.toMatrix(B));
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] state) {
        double[] v = new double[this.param.getNTypes()*this.param.getNTypes()];

        RealMatrix constantMatrix = new BlockRealMatrix(this.param.getNTypes(), this.param.getNTypes()); //  this.buildTimeInvariantSystemMatrix(this.param.getTotalProcessLength());
        Utils.fillArray(constantMatrix, v);
        this.realConstantMatrix = constantMatrix;

        this.constantMatrix = new DenseMatrix(this.param.getNTypes(), this.param.getNTypes());
        Utils.fillMatrix(this.constantMatrix, v);

        return super.integrateBackwardsOverIntegrals(state);
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        int numTypes = this.param.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);

        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    public static double[] integrateUsingFlow(
            double timeStart,
            int intervalStart,
            double timeEnd,
            int intervalEnd,
            double[] initialState,
            Flow flow
    ) {
        RealMatrix flowMatrixStart = flow.getFlow(timeStart, intervalStart);
        RealMatrix flowMatrixEnd = flow.getFlow(timeEnd, intervalEnd);

        RealVector likelihoodVectorEnd = Utils.toVector(initialState);

        try {
            DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEnd).getSolver();
            RealVector solution = linearSolver.solve(likelihoodVectorEnd);

            RealVector likelihoodVectorStart = flowMatrixStart.operate(solution);

            return likelihoodVectorStart.toArray();
        } catch (SingularMatrixException e) {
            throw e;
        }
    }

    @Override
    protected void handleIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.param.getNTypes(); i++) {
            for (int j = 0; j < this.param.getNTypes(); j++) {
                state[i*this.param.getNTypes() + j] *= (1 - this.param.getRhoValues()[newInterval][i]);
            }
        }
    }
}
