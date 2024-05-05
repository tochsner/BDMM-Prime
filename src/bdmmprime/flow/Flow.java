package bdmmprime.flow;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class Flow {
    ContinuousOutputModel[] outputModels;
    ContinuousOutputModel posExponentials;
    int n;
    double totalTime;

    public Flow(ContinuousOutputModel[] outputModels, int n) {
        this.outputModels = outputModels;
        this.n = n;
    }

    public Flow(ContinuousOutputModel[] outputModels, int n, ContinuousOutputModel posExponentials, double totalTime) {
        this(outputModels, n);
        this.posExponentials = posExponentials;
        this.totalTime = totalTime;
    }

    public RealMatrix getFlow(double time, int intervalIndex) {
        outputModels[intervalIndex].setInterpolatedTime(time);
        double[] flow = outputModels[intervalIndex].getInterpolatedState();
        RealMatrix flowMatrix = Utils.toMatrix(flow, n);

        if (this.posExponentials == null) return flowMatrix;

        posExponentials.setInterpolatedTime(time - totalTime);
        double[] exp = posExponentials.getInterpolatedState();
        RealMatrix expMatrix = Utils.toMatrix(exp, n);

        return expMatrix.multiply(flowMatrix);
    }
}
