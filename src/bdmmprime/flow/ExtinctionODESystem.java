package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;


public class ExtinctionODESystem extends IntervalODESystem {
    public ExtinctionODESystem(Parameterization parameterization) {
        super(parameterization);
    }

    @Override
    public int getDimension() {
        return this.param.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        for (int i = 0; i < this.param.getNTypes(); i++) {
            yDot[i] = param.getBirthRates()[interval][i] * y[i];
            yDot[i] += param.getDeathRates()[interval][i] * y[i];
            yDot[i] += param.getSamplingRates()[interval][i] * y[i];

            yDot[i] -= param.getBirthRates()[interval][i] * y[i] * y[i];
            yDot[i] -= param.getDeathRates()[interval][i];

            for (int j = 0; j < this.param.getNTypes(); j++) {
                if (i == j) continue;

                yDot[i] += param.getCrossBirthRates()[interval][i][j] * (y[i] - y[i] * y[j]);
                yDot[i] += param.getMigRates()[interval][i][j] * (y[i] - y[j]);
            }
        }
    }
}
