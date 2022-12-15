#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include "matrixInverse.h"

void isInvertable(const float *backMatrix) { //inverse matrix's elements are equal if its determinant is zero 
    for (int i = 0; i < SQUARED_SIZE - 1; i++) {
        assert((backMatrix[i] != backMatrix[i+1]) && "Matrix is not invertable (the determinant is zero)\n");
    }
}

void isInverseCorrect(const float *A, const float *M) {
    auto I = new float [SQUARED_SIZE];
    multMatrices_BLAS_OPT(A, M, I);
    // std::cout << "I\n"; printMatrix(I);
    float module = 0.0f;
    for (int i = 0; i < SQUARED_SIZE; i++) {
        module += M[i] * M[i];
    }
    module = sqrt(module);
    float inaccuracy = (float)fabs(1 - module);
    std::cout << "inaccuracy: " << inaccuracy << "\n";
    assert((inaccuracy < 0.1f && "too big inaccuracy\n"));
    delete [] I;
}

void mainTest(const float *A, const float *M) {
    isInvertable(M);
    std::cout << "Matrix is invertable\n";
    isInverseCorrect(A, M);
    std::cout << "Inverse matrix is correct\n\n";
}