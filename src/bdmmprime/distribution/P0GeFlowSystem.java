package bdmmprime.distribution;


import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.analysis.UnivariateFunction;

public class P0GeFlowSystem extends P0GeSystem {

    UnivariateFunction flow;

	public P0GeFlowSystem(Parameterization parameterization,
                          double absoluteTolerance,
                          double relativeTolerance) {

	    super(parameterization, absoluteTolerance, relativeTolerance);
        this.calculateFlow();
	}

    private void calculateFlow() {
        // precalculate flow

        // times, trajectories = integrate(0, T)
        // this.flow = LinearInterpolator().interpolate(times, trajectories)
    }

    @Override
    public ScaledNumbers safeIntegrate(ScaledNumbers pgScaled, double tStart, double tEnd) {
        // integrate using the precalculated flow

        // flow = this.flow.value(tStart)
        // y = linearSolve(flow * y = state)
        // state = flow(tEnd) * y
    }

    @Override
    public void integrate(P0State state, double tStart, double tEnd) {
        // integrate using the precalculated flow

        // flow = this.flow.value(tStart)
        // y = linearSolve(flow * y = state)
        // state = flow(tEnd) * y
    }
}

