#pragma once

#include <array>
#include <cassert>
#include "matrixInverse.h"

void isMatricesEqual(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M2,
                     const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M3) {
    assert((M1 == M2) && "Matrices M1 and M2 are not equal\n");
    assert((M1 == M3) && "Matrices M1 and M3 are not equal\n");
    assert((M2 == M3) && "Matrices M2 and M3 are not equal\n");
}

void isInvertable(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &backMatrix) { //inverse matrix's elements are equal if its determinant is zero 
    for (int i = 0; i < MATRIX_SIZE*MATRIX_SIZE - 1; i++) {
        assert((backMatrix.at(i) != backMatrix.at(i+1)) && "Matrix is not invertable (the determinant is zero)\n");
    }
}

void testResult(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &A, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &backMatrix) {
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> I = {0};
    multMatrices_BLAS_OPT(A, backMatrix, I);
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            
        }
    }
}

void mainTest(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &A, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M1, 
              const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M2, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M3) {
    isMatricesEqual(M1, M2, M3);
    std::cout << "Matrices are equal\n";
    isInvertable(M1);
    std::cout << "Matrix is invertable\n";
    testResult(A, M1);
    std::cout << "All tests are passed\n";
}