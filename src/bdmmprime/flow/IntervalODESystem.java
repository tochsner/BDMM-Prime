package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected int interval;

    private final double integrationMinStep;
    private final double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization) {
        this.param = parameterization;
        this.integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        this.integrationMaxStep = this.param.getTotalProcessLength() / 10;
    }

    protected void setInterval(int interval) {
        this.interval = interval;
    }

    public ContinuousOutputModel integrateOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                this.integrationMinStep, this.integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        ContinuousOutputModel result = new ContinuousOutputModel();

        double[] intervalEndTimes = this.param.getIntervalEndTimes();

        for (int interval = 0; interval < this.param.getTotalIntervalCount(); interval++) {
            this.setInterval(interval);

            double startTime = interval == 0 ? 0 : intervalEndTimes[interval - 1];
            double endTime = intervalEndTimes[interval];

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();
            integrator.addStepHandler(intervalResult);

            integrator.integrate(this, startTime, state, endTime, state);

            result.append(intervalResult);
            integrator.clearStepHandlers();
        }

        return result;
    }

    public ContinuousOutputModel integrateBackwardsOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                this.integrationMinStep, this.integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        ContinuousOutputModel result = new ContinuousOutputModel();

        double[] intervalEndTimes = this.param.getIntervalEndTimes();

        for (int interval = this.param.getTotalIntervalCount() - 1; 0 <= interval; interval--) {
            this.setInterval(interval);

            double startTime = intervalEndTimes[interval];
            double endTime = interval == 0 ? 0 : intervalEndTimes[interval - 1];

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();
            integrator.addStepHandler(intervalResult);

            integrator.integrate(this, startTime, state, endTime, state);

            result.append(intervalResult);
            integrator.clearStepHandlers();
        }

        return result;
    }

}
