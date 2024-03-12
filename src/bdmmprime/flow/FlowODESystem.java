package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class FlowODESystem extends IntervalODESystem {
    private final ContinuousOutputModel extinctionProbabilities;

    public FlowODESystem(Parameterization parameterization, ContinuousOutputModel extinctionProbabilities) {
        super(parameterization);
        this.extinctionProbabilities = extinctionProbabilities;
    }

    @Override
    public int getDimension() {
        return param.getNTypes() * param.getNTypes();
    }

    private RealMatrix buildSystemMatrix(double t) {
        RealMatrix system = new BlockRealMatrix(param.getNTypes(), param.getNTypes());

        int interval = this.param.getIntervalIndex(t);

        this.extinctionProbabilities.setInterpolatedTime(t);
        double[] extinctProbabilities = this.extinctionProbabilities.getInterpolatedState();

        // fill transitions
        for (int i = 0; i < param.getNTypes(); i++) {
            for (int j = 0; j < param.getNTypes(); j++) {
                if (i == j) continue;

                system.setEntry(
                        i,
                        j,
                        this.param.getMigRates()[interval][i][j]
                                + this.param.getCrossBirthRates()[interval][i][j] * extinctProbabilities[i]
                );
            }
        }

        // fill diagonals
        for (int i = 0; i < param.getNTypes(); i++) {
            system.setEntry(
                    i,
                    i,
                    -this.param.getDeathRates()[interval][i] - this.param.getSamplingRates()[interval][i]
            );
            system.addToEntry(
                    i,
                    i,
                    this.param.getBirthRates()[interval][i]
                            * (2*extinctProbabilities[i] - 1)
            );

            for (int j = 0; j < param.getNTypes(); j++) {
                if (i == j) continue;

                system.addToEntry(
                        i,
                        i,
                        -this.param.getMigRates()[interval][i][j]
                );
                system.addToEntry(
                        i,
                        i,
                        this.param.getCrossBirthRates()[interval][i][j]
                                * (extinctProbabilities[j] - 1)
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

        RealMatrix yDotMatrix = yMatrix.multiply(systemMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }
}
