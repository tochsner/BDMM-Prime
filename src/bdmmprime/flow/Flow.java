package bdmmprime.flow;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class Flow {
    ContinuousOutputModel[] outputModels;
    int n;

    public Flow(ContinuousOutputModel[] outputModels, int n) {
        this.outputModels = outputModels;
        this.n = n;
    }

    public RealMatrix getFlow(double time, int intervalIndex) {
        outputModels[intervalIndex].setInterpolatedTime(time);
        double[] flow = outputModels[intervalIndex].getInterpolatedState();
        return Utils.toMatrix(flow, n);
    }
}
