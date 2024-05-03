package bdmmprime.flow;

import org.apache.commons.math3.ode.ContinuousOutputModel;

public class ExtinctionProbabilities {
    ContinuousOutputModel[] outputModels;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels) {
        this.outputModels = outputModels;
    }

    public double[] getProbabilities(double time, int intervalIndex) {
        outputModels[intervalIndex].setInterpolatedTime(time);
        return outputModels[intervalIndex].getInterpolatedState();
    }
}
