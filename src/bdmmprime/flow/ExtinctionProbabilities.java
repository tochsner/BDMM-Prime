package bdmmprime.flow;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.ContinuousOutputModel;

public class ExtinctionProbabilities {
    ContinuousOutputModel[] outputModels;
    double[] endTimes;

    public ExtinctionProbabilities(ContinuousOutputModel[] flows, double[] endTimes) {
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

        this.outputModels = flows;
        this.endTimes = endTimes;
    }

    protected double[] getExtinctionProbability(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        return output.getInterpolatedState();
    }

    public double[] getExtinctionProbability(double time) {
        for (int i = 0; i < this.outputModels.length; i++) {
            if (time <= this.endTimes[i]) {
                return this.getExtinctionProbability(this.outputModels[i], time);
            }
        }

        return this.getExtinctionProbability(this.outputModels[this.outputModels.length - 1], time);
    }
}
