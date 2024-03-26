package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.*;


public abstract class IntervalODESystem implements FirstOrderDifferentialEquations, EventHandler {

    protected Parameterization param;
    protected int interval;

    private final double integrationMinStep;
    private final double integrationMaxStep;

    private final int MIN_RATE_CHANGE_NUM_CHECKS = 10;
    private final double RATE_CHANGE_CHECK_CONVERGENCE = 1e-5;
    private final int RATE_CHANGE_MAX_ITERATIONS = 1000;

    public IntervalODESystem(Parameterization parameterization) {
        this.param = parameterization;
        this.integrationMinStep = this.param.getTotalProcessLength() * 1e-100;
        this.integrationMaxStep = this.param.getTotalProcessLength() / 20;
    }

    public ContinuousOutputModel integrateOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        return this.integrateOverIntegrals(
                0, this.param.getTotalProcessLength(), state, absoluteTolerance, relativeTolerance
        );
    }

    public ContinuousOutputModel integrateBackwardsOverIntegrals(double[] state, double absoluteTolerance, double relativeTolerance) {
        return this.integrateOverIntegrals(
                this.param.getTotalProcessLength(), 0, state, absoluteTolerance, relativeTolerance
        );
    }


    public ContinuousOutputModel integrateOverIntegrals(
            double timeStart,
            double timeEnd,
            double[] state,
            double absoluteTolerance,
            double relativeTolerance
    ) {
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                this.integrationMinStep, this.integrationMaxStep, absoluteTolerance, relativeTolerance
        );

        integrator.addEventHandler(
                this,
                (timeEnd-timeStart)/ MIN_RATE_CHANGE_NUM_CHECKS,
                RATE_CHANGE_CHECK_CONVERGENCE, RATE_CHANGE_MAX_ITERATIONS
        );

        ContinuousOutputModel result = new ContinuousOutputModel();
        integrator.addStepHandler(result);

        this.interval = this.param.getIntervalIndex(timeStart);
        integrator.integrate(this, timeStart, state, timeEnd, state);

        return result;
    }

    /* EventHandler Implementation */

    @Override
    public void init(double v, double[] doubles, double v1) { }

    @Override
    public double g(double v, double[] doubles) {
        double result = 1.0;
        for (double boundary : param.getIntervalEndTimes())
            result *= v - boundary;

        return result;
    }

    @Override
    public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
        return EventHandler.Action.RESET_STATE;
    }

    @Override
    public void resetState(double t, double[] y) {
        this.interval = this.param.getIntervalIndex(t);
    }

}
