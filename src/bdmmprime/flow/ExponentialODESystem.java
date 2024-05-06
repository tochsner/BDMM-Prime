package bdmmprime.flow;

import bdmmprime.parameterization.Parameterization;
import cern.colt.matrix.impl.RCDoubleMatrix2D;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


public class ExponentialODESystem extends IntervalODESystem {
    RealMatrix constantMatrix;

    public ExponentialODESystem(RealMatrix constantMatrix, Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, absoluteTolerance, relativeTolerance);
        this.constantMatrix = constantMatrix;
    }

    @Override
    public int getDimension() {
        return this.param.getNTypes() * this.param.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
        RealMatrix yMatrix = Utils.toMatrix(y, this.param.getNTypes());
        RealMatrix yDotMatrix = this.constantMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }
}
