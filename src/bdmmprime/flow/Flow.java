package bdmmprime.flow;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class Flow {
    ContinuousOutputModel[] outputModels;
    DoubleMatrix constantMatrix;
    int n;
    double totalTime;

    public Flow(ContinuousOutputModel[] outputModels, int n) {
        this.outputModels = outputModels;
        this.n = n;
    }

    public Flow(ContinuousOutputModel[] outputModels, int n, DoubleMatrix constantMatrix, double totalTime) {
        this(outputModels, n);
        this.constantMatrix = constantMatrix;
        this.totalTime = totalTime;
    }

    public RealMatrix getFlow(double time, int intervalIndex) {
        outputModels[intervalIndex].setInterpolatedTime(time);
        double[] flow = outputModels[intervalIndex].getInterpolatedState();
        RealMatrix flowMatrix = Utils.toMatrix(flow, n);

        if (this.constantMatrix == null) return flowMatrix;

        RealMatrix expMatrix = Utils.toMatrix(MatrixFunctions.expm(this.constantMatrix.mul(time - this.totalTime)));

        return expMatrix.multiply(flowMatrix);
    }
}
