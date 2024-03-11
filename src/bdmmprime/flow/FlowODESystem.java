package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class FlowODESystem implements FirstOrderDifferentialEquations {

    private Parameterization param;
    private int interval;

    public FlowODESystem(Parameterization parameterization) {
        this.param = parameterization;
    }

    public void setInterval(int interval) {
        this.interval = interval;
    }


    /* FirstOrderDifferentialEquations implementation */

    @Override
    public int getDimension() {
        return param.getNTypes()*2;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {

        int nTypes = param.getNTypes();

        for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			yDot[i] = (param.getBirthRates()[interval][i]
                    + param.getDeathRates()[interval][i]
                    + param.getSamplingRates()[interval][i]) * y[i]
					- param.getBirthRates()[interval][i] * y[i] * y[i]
					- param.getDeathRates()[interval][i];

			for (int j = 0; j < nTypes; j++){

			    if (i==j)
			        continue;

                yDot[i] += param.getCrossBirthRates()[interval][i][j] * (y[i] - y[i]*y[j]);
                yDot[i] += param.getMigRates()[interval][i][j] * (y[i] - y[j]);
			}

			/*  ge equations: (dim .. 2*dim-1) */

			yDot[nTypes + i] = (param.getBirthRates()[interval][i]
                    + param.getDeathRates()[interval][i]
                    + param.getSamplingRates()[interval][i])*y[nTypes+i]
                    - 2*param.getBirthRates()[interval][i]*y[nTypes+i]*y[i];

			for (int j = 0; j< nTypes; j++) {

                if (i==j)
			        continue;

                yDot[nTypes + i] += param.getCrossBirthRates()[interval][i][j]
                        * (y[nTypes+i] - (y[nTypes+i]*y[j] + y[nTypes+j]*y[i]));

                yDot[nTypes + i] += param.getMigRates()[interval][i][j]
                        * (y[nTypes+i] - y[nTypes+j]);
			}
		}

    }
}
