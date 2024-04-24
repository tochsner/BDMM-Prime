package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected Parameterization param;
    protected int interval;

    protected FirstOrderIntegrator integrator;

    public IntervalODESystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        this.param = parameterization;
        double integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        double integrationMaxStep = this.param.getTotalProcessLength() / 20;
        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep, absoluteTolerance, relativeTolerance
        );
    }

    public ContinuousOutputModel[] integrateBackwardsOverIntegrals(double[] state) {
        this.interval = this.param.getTotalIntervalCount() - 1;

        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[this.param.getTotalIntervalCount()];

        for (int interval = this.param.getTotalIntervalCount() - 1; interval >= 0; interval--) {
            double endTime = this.param.getIntervalEndTimes()[interval];
            double startTime = 0 < interval ? this.param.getIntervalEndTimes()[interval - 1] : 0;

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();
            integrator.addStepHandler(intervalResult);

            integrator.integrate(this, endTime, state, startTime, state);

            integrator.clearStepHandlers();

            if (0 < interval) {
                this.handleIntervalBoundary(startTime, interval, interval - 1, state);
            }

            outputModels[interval] = (intervalResult);
        }

        return outputModels;
    }

    protected void handleIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        this.interval = newInterval;
    }

}
