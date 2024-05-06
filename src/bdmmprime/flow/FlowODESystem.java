package bdmmprime.flow;

import bdmmprime.flow.benchmark.MatrixExponential;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jblas.Eigen;

public class FlowODESystem extends IntervalODESystem {
    private final ExtinctionProbabilities extinctionProbabilities;
    public RealMatrix constantMatrix;
    public DoubleMatrix jblasConstantMatrix;

    public FlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
            double absoluteTolerance,
            double relativeTolerance
    ) {
        super(parameterization, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;
        this.setUpExponentials(absoluteTolerance, relativeTolerance);
    }

    protected void setUpExponentials(
            double absoluteTolerance,
            double relativeTolerance
    ) {
        this.constantMatrix = this.buildTimeInvariantSystemMatrix(this.param.getTotalProcessLength());
        this.jblasConstantMatrix = Utils.toMatrix(this.constantMatrix);
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

        RealMatrix A = Utils.toMatrix(MatrixFunctions.expm(this.jblasConstantMatrix.mul(this.param.getTotalProcessLength() - t)));
        RealMatrix B = Utils.toMatrix(MatrixFunctions.expm(this.jblasConstantMatrix.mul(t - this.param.getTotalProcessLength())));

        return A.multiply(timeVaryingMatrix.add(timeInvariantMatrix).subtract(constantMatrix)).multiply(B);
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
                state[i * this.param.getNTypes() + j] *= (1 - this.param.getRhoValues()[newInterval][i]);
            }
        }
    }
}
