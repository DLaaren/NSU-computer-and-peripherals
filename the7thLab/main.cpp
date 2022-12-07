#include "tests.h"

int main() {
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> A = {0};
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> backMatrix1 = {0}; //A^(-1)
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> backMatrix2 = {0};
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> backMatrix3 = {0};

    A = generateMatrix();

    int type = NO_OPT;
    double res;
    res = matrixInverse(A, backMatrix1, type); //returns time 
    std::cout << "no optimization: " << res << " ms\n";

    type = VECTOR_OPT;
    res = matrixInverse(A, backMatrix2, type); //returns time 
    std::cout << "vector optimization: " << res << " ms\n";

    type = BLAS_OPT;
    res = matrixInverse(A, backMatrix3, type); //returns time 
    std::cout << "blas optimization: " << res << " ms\n";

    //mainTest(A, backMatrix1, backMatrix2, backMatrix3);

    return 0;
}

