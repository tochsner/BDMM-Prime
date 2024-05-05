package bdmmprime.flow;

import org.junit.Test;

import static jeigen.Shortcuts.*;
import jeigen.DenseMatrix;
import jeigen.SparseMatrixLil;

public class JeigenTest {

    @Test
    public void testJeigen() {
        DenseMatrix dm1;
        DenseMatrix dm2;
        dm1 = new DenseMatrix("1 2; 3 4"); // create new matrix
        // with rows {1,2} and {3,4}
        dm1.mexp();
        dm1 = new DenseMatrix(new double[][]{{1, 2}, {3, 4}}); // create new matrix
        // with rows {1,2} and {3,4}
        dm1 = zeros(5, 3);  // creates a dense matrix with 5 rows, and 3 columns
        dm1 = rand(5, 3); // create a 5*3 dense matrix filled with random numbers
        dm1 = ones(5, 3);  // 5*3 matrix filled with '1's
        dm1 = diag(rand(5, 1)); // creates a 5*5 diagonal matrix of random numbers
        dm1 = eye(5); // creates a 5*5 identity matrix

        SparseMatrixLil sm1;
        sm1 = spzeros(5, 3); // creates an empty 5*3 sparse matrix
        sm1 = spdiag(rand(5, 1)); // creates a sparse 5*5 diagonal matrix of random
        // numbers
        sm1 = speye(5); // creates a 5*5 identity matrix, sparse
    }

}
