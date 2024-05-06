package bdmmprime.flow;

import jeigen.DenseMatrix;
import jeigen.SparseMatrixLil;
import org.jblas.MatrixFunctions;
import org.junit.Test;

import org.jblas.DoubleMatrix;

public class JblasTest {

    @Test
    public void testJeigen() {
        DoubleMatrix a = DoubleMatrix.rand(10, 10);
        MatrixFunctions.expm(a);
    }

}
