package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


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
                yDot[i] += param.getCrossBirthRates()[interval][i][j] * (y[i] - y[i] * y[j]);
                yDot[i] += param.getMigRates()[interval][i][j] * (y[i] - y[j]);
            }
        }
    }

    @Override
    public ContinuousOutputModel integrateBackwardsOverIntegrals(
            double timeStart,
            double timeEnd,
            double[] state,
            double absoluteTolerance,
            double relativeTolerance
    ) {
        if (timeStart > timeEnd) {
            throw new IllegalArgumentException("timeStart has to be smaller than timeEnd.");
        }

        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                this.integrationMinStep, this.integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        ContinuousOutputModel result = new ContinuousOutputModel();

        int startInterval = this.param.getIntervalIndex(timeStart);
        int endInterval = this.param.getIntervalIndex(timeEnd);

        double[] intervalEndTimes = this.param.getIntervalEndTimes();

        for (int interval = endInterval; interval >= startInterval; interval--) {
            this.setInterval(interval);

            double endTime = interval == endInterval ? timeEnd : intervalEndTimes[interval];
            double startTime = interval == startInterval ? timeStart : intervalEndTimes[interval - 1];

            if (Utils.equalWithPrecision(endTime, intervalEndTimes[interval])) {
                // we have to adjust the boundary condition for rho sampling
                for (int i = 0; i < this.param.getNTypes(); i++) {
                    state[i] *= (1 - this.param.getRhoValues()[interval][i]);
                }
            }

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();
            integrator.addStepHandler(intervalResult);

            integrator.integrate(this, endTime, state, startTime, state);

            result.append(intervalResult);
            integrator.clearStepHandlers();
        }

        return result;
    }
}
