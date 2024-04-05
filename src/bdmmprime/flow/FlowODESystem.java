package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.util.Arrays;

public class FlowODESystem extends IntervalODESystem {
    private final ContinuousOutputModel[] extinctionProbabilities;

    public FlowODESystem(
            Parameterization parameterization,
            ContinuousOutputModel[] extinctionProbabilities,
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

    protected RealMatrix buildSystemMatrix(double t) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        int interval = this.param.getIntervalIndex(t);

        this.extinctionProbabilities[interval].setInterpolatedTime(t);
        double[] extinctProbabilities = this.extinctionProbabilities[interval].getInterpolatedState();

        // fill transitions

        for (int i = 0; i < param.getNTypes(); i++) {
            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        j,
                        -this.param.getMigRates()[interval][i][j]
                                - this.param.getCrossBirthRates()[interval][i][j] * extinctProbabilities[i]
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
            system.addToEntry(
                    i,
                    i,
                    -this.param.getBirthRates()[interval][i] * (2 * extinctProbabilities[i] - 1)
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.param.getMigRates()[interval][i][j]
                );
                system.addToEntry(
                        i,
                        i,
                        -this.param.getCrossBirthRates()[interval][i][j] * (extinctProbabilities[j] - 1)
                );
            }
        }

        return system;
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
            ContinuousOutputModel[] flow
    ) {
        int n = initialState.length;

        flow[intervalStart].setInterpolatedTime(timeStart);
        double[] flowStart = flow[intervalStart].getInterpolatedState();
        RealMatrix flowMatrixStart = Utils.toMatrix(flowStart, n);

        flow[intervalEnd].setInterpolatedTime(timeEnd);
        double[] flowEnd = flow[intervalEnd].getInterpolatedState();
        RealMatrix flowMatrixEnd = Utils.toMatrix(flowEnd, n);

        RealVector likelihoodVectorEnd = Utils.toVector(initialState);

        DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEnd).getSolver();
        RealVector solution = linearSolver.solve(likelihoodVectorEnd);

        RealVector likelihoodVectorStart = flowMatrixStart.operate(solution);

        return likelihoodVectorStart.toArray();
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
