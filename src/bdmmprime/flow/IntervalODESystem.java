package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected int interval;

    public IntervalODESystem(Parameterization parameterization) {
        this.param = parameterization;
    }

    protected void setInterval(int interval) {
        this.interval = interval;
    }

    public ContinuousOutputModel integrateOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 10;
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        ContinuousOutputModel result = new ContinuousOutputModel();
        integrator.addStepHandler(result);

        // run integration over the entire timespan (all the intervals) to calculate the flow

        double[] intervalEndTimes = this.param.getIntervalEndTimes();

        for (int interval = 0; interval < this.param.getTotalIntervalCount(); interval++) {
            this.setInterval(interval);

            double startTime = interval == 0 ? 0 : intervalEndTimes[interval - 1];
            double endTime = intervalEndTimes[interval];

            integrator.integrate(this, startTime, state, endTime, state);
        }

        return result;
    }

    public ContinuousOutputModel integrateBackwardsOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 10;
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        ContinuousOutputModel result = new ContinuousOutputModel();
        integrator.addStepHandler(result);

        // run integration over the entire timespan (all the intervals) to calculate the flow

        double[] intervalEndTimes = this.param.getIntervalEndTimes();

        for (int interval = this.param.getTotalIntervalCount() - 1; 0 <= interval; interval--) {
            this.setInterval(interval);

            double startTime = intervalEndTimes[interval];
            double endTime = interval == 0 ? 0 : intervalEndTimes[interval - 1];

            integrator.integrate(this, startTime, state, endTime, state);
        }

        return result;
    }


}
