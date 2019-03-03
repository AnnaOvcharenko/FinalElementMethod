package matrixEquation;

public class MatrixEquation {
    private static final double EPSILON = 1e-10;
    public MatrixEquation() {
    }

    public double[] solveEquationGauss(double[][] A, double[] b) {
        int n = b.length;

        for (int p = 0; p < n; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double   t    = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (Math.abs(A[p][p]) <= EPSILON) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;

    }

    public double[] solveEquationKramer(double[][] A, double[] b) {
        double[] resultVector = new double[A.length];
        double[][] tmp = new double[A.length][A.length];
        double aDet = determinant(A, A.length);

        for (int i = 0; i < A.length; ++i) {
            //tmp = a;
            for(int j=0;j<A[i].length;++j){
                System.arraycopy(A[j],0,tmp[j],0,A.length);
            }
            for (int j=0;j<A.length;++j){
                tmp[j][i]=b[j];
            }
            resultVector[i]=determinant(tmp,tmp.length)/aDet;
        }


        return resultVector;
    }

    public static double determinant(double A[][], int n) {
        double det = 0;
        if (n == 1) {
            det = A[0][0];
        } else if (n == 2) {
            det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
        } else if (n == 3) {
            det = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[2][1] * A[0][2] - A[0][2] * A[1][1] * A[2][0] - A[0][0] * A[1][2] * A[2][1] - A[2][2] * A[1][0] * A[0][1];
//            00 01 02
//            10 11 12
//            20 21 22
        } else {
            det = 0;
            for (int j1 = 0; j1 < n; j1++) {
                double[][] m = new double[n - 1][];
                for (int k = 0; k < (n - 1); k++) {
                    m[k] = new double[n - 1];
                }
                for (int i = 1; i < n; i++) {
                    int j2 = 0;
                    for (int j = 0; j < n; j++) {
                        if (j == j1)
                            continue;
                        m[i - 1][j2] = A[i][j];
                        j2++;
                    }
                }
                det += Math.pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * determinant(m, n - 1);
            }
        }
        return det;
    }
}
