package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected int interval;

    protected final double integrationMinStep;
    protected final double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization) {
        this.param = parameterization;
        this.integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        this.integrationMaxStep = this.param.getTotalProcessLength() / 20;
    }

    protected void setInterval(int interval) {
        this.interval = interval;
    }

    public ContinuousOutputModel integrateBackwardsOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        return this.integrateBackwardsOverIntegrals(
                0, this.param.getTotalProcessLength(), state, absoluteTolerance, relativeTolerance
        );
    }

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

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();
            integrator.addStepHandler(intervalResult);

            integrator.integrate(this, endTime, state, startTime, state);

            result.append(intervalResult);
            integrator.clearStepHandlers();
        }

        return result;
    }

}
