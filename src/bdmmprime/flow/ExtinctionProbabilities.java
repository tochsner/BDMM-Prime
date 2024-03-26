package bdmmprime.flow;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class Flow {
    ContinuousOutputModel[] outputModels;
    double[] endTimes;
    RealMatrix[] initialFlows;

    int n;

    public Flow(ContinuousOutputModel[] flows, double[] endTimes, int n) {
        // test if array length match

        if (flows.length != endTimes.length) {
            throw new IllegalArgumentException("The number of flows must match the number of provided end times.");
        }

        // test if endTimes are sorted

        for (int i = 0; i < endTimes.length - 1; i++) {
            if (endTimes[i] > endTimes[i + 1]) {
                throw new IllegalArgumentException("The provided end times must be sorted.");
            }
        }


        this.n = n;
        this.outputModels = flows;
        this.endTimes = endTimes;

        this.initialFlows = new RealMatrix[this.outputModels.length];

        this.initialFlows[0] = MatrixUtils.createRealIdentityMatrix(this.n);

        for (int i = 1; i < this.initialFlows.length; i++) {
            RealMatrix flowEnd = this.getFlow(this.outputModels[i - 1], this.endTimes[i - 1]);
            this.initialFlows[i] = flowEnd.multiply(this.initialFlows[i - 1]);
        }
    }

    protected RealMatrix getFlow(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        double[] flow = output.getInterpolatedState();
        return Utils.toMatrix(flow, this.n);
    }

    public RealMatrix getFlow(double time) {
        for (int i = 0; i < this.outputModels.length; i++) {
            if (time <= this.endTimes[i]) {
                return this.getFlow(this.outputModels[i], time).multiply(this.initialFlows[i]);
            }
        }

        return this.initialFlows[this.initialFlows.length - 1].multiply(
                this.getFlow(this.outputModels[this.initialFlows.length - 1], time)
        );
    }
}
