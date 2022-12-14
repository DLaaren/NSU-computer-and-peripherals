#pragma once

#include <array>
#include <cassert>
#include "matrixInverse.h"

void isInvertable(const float *backMatrix) { //inverse matrix's elements are equal if its determinant is zero 
    for (int i = 0; i < SQUARED_SIZE - 1; i++) {
        assert((backMatrix[i] != backMatrix[i+1]) && "Matrix is not invertable (the determinant is zero)\n");
    }
}

void isInverseCorrect(const float *A, const float *M) {
    auto I = new float [SQUARED_SIZE];
    multMatrices_BLAS_OPT(A, M, I);
    std::cout << "I\n";
    printMatrix(I);
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            if (i == j) {
                assert(((float)1 - (float)fabs(I[i*MATRIX_SIZE + j]) <= 0.4f) && "Incorrect inverse matrix\n");
                continue;
            }
            assert(((float)1 - (float)fabs(I[i*MATRIX_SIZE + j]) >= 0.01f) && "Incorrect inverse matrix\n");
        }
    }

    delete [] I;
}

void mainTest(const float *A, const float *M) {
    isInvertable(M);
    std::cout << "Matrix is invertable\n";
    isInverseCorrect(A, M);
    std::cout << "Inverse matrix is correct\n";
}