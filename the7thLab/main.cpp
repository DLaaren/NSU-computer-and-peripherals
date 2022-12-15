#define NO_OPT 0
#define VECTOR_OPT 1
#define BLAS_OPT 2

#include "matrixInverse.h"
#include "tests.h"

int main() {
    auto A = new float [SQUARED_SIZE];
    auto backMatrix = new float [SQUARED_SIZE]; //A^(-1)
    memset(backMatrix, 0, sizeof(float) * SQUARED_SIZE);

    generateMatrix(A);

    //printMatrix(A);

    int type = NO_OPT;
    double res;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "no optimization: " << res << " ms\n";
    mainTest(A, backMatrix);
    // printMatrix(backMatrix);
    memset(backMatrix, 0, sizeof(float) * SQUARED_SIZE);    

    type = VECTOR_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "vector optimization: " << res << " ms\n";
    mainTest(A, backMatrix);
    // printMatrix(backMatrix);
    memset(backMatrix, 0, sizeof(float) * SQUARED_SIZE);

    type = BLAS_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "blas optimization: " << res << " ms\n";
    mainTest(A, backMatrix);
    // printMatrix(backMatrix);

    delete [] A;
    delete [] backMatrix;
 
    return 0;
}
