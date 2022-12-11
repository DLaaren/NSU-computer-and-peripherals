#pragma once

#include <iostream>
#include <cstring>
#include <cmath>
#include <xmmintrin.h>

#define MATRIX_SIZE 9
#define SQUARED_SIZE ((MATRIX_SIZE)*(MATRIX_SIZE)) 
#define MATRIX_SIZE_VECTOR_OPT ((MATRIX_SIZE) - ((MATRIX_SIZE) % (4)))
#define iterations 10000
//for tests it should be 10'000

#include "matrixInverseNOOPT.h"
#include "matrixInverseVECTOROPT.h"
#include "matrixInverseBLASOPT.h"
#include "generator.h"

void printMatrix(const float *M) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            std::cout << M[i*MATRIX_SIZE + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


void generateMatrix(float *A) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i*MATRIX_SIZE + j] = std::rand() % 10;
        }
    }
}

void createIdentityMatrix(float *M);
void findMatrixB(const float *A, float *B);
void multMatrices(const float *M1, const float *M2, float *res, int type);
void sumMatrices(const float *M1, const float *M2, bool sign, float *res, int type);


double matrixInverse(const float *A, float *resMatrix, int type) {
    clock_t start = clock();
    
    auto I = new float [SQUARED_SIZE];
    auto R = new float [SQUARED_SIZE]; //R = I - B*A
    auto B = new float [SQUARED_SIZE];
    auto tmp = new float [SQUARED_SIZE];

    memset(I, 0, sizeof(float) * SQUARED_SIZE);
    memset(R, 0, sizeof(float) * SQUARED_SIZE);
    memset(B, 0, sizeof(float) * SQUARED_SIZE);
    memset(tmp, 0, sizeof(float) * SQUARED_SIZE);

    createIdentityMatrix(I);
    findMatrixB(A, B);   

    multMatrices(B, A, tmp, type); //B*A = tmp
    sumMatrices(I, tmp, false, R, type); // I - tmp = R
    sumMatrices(I, R, true, resMatrix, type); //I + R = resMatrix
    multMatrices(R, R, tmp, type); //tmp = R^2

    auto R_k = new float[SQUARED_SIZE];
    memset(R_k, 0, sizeof(float) * SQUARED_SIZE);

    for (int k = 2; k < iterations - 1; k++) {
        sumMatrices(resMatrix, tmp, true, R_k, type); //resMatrix + R^k
        memcpy(resMatrix, R_k, sizeof(float) * SQUARED_SIZE);
        memset(R_k, 0, sizeof(float) * SQUARED_SIZE);
        multMatrices(tmp, R, R_k, type); //R^(k+1)
        memcpy(tmp, R_k, sizeof(float) * SQUARED_SIZE);
    }
    memset(R_k, 0, sizeof(float) * SQUARED_SIZE);
    sumMatrices(resMatrix, tmp, true, R_k, type); //last addition (resMatrix + R^9)

    multMatrices(R_k, B, resMatrix, type); // (I + R + ... + R^9) * B
    
    delete [] R_k;
    delete [] B;
    delete [] I;
    delete [] tmp;
    delete [] R;

    clock_t end = clock();
    double time = 1000.0 * ((double)(end - start)/CLOCKS_PER_SEC);
    return time;
}

void createIdentityMatrix(float *I) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        I[i*MATRIX_SIZE + i] = 1;
    }
}

void findMatrixB(const float *A, float *B) {
    //||A||1 = максимальный по сумме элементов столбец -> maxSumJ
    //||A||inf = максимальная по сумме элементов строчка -> maxSumI
    int maxSumI = 0, maxSumJ = 0;
    int tmpSumI = 0, tmpSumJ = 0;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            tmpSumI += (float)fabs(A[i*MATRIX_SIZE + j]);
            tmpSumJ += (float)fabs(A[j*MATRIX_SIZE + i]);
        }
        if (tmpSumI > maxSumI) maxSumI = tmpSumI;
        if (tmpSumJ > maxSumJ) maxSumJ = tmpSumJ;
        tmpSumI = 0;
        tmpSumJ = 0;
    }

    float divider = maxSumI * maxSumJ;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            B[i*MATRIX_SIZE +j] = A[j*MATRIX_SIZE + i] / divider;
        }
    }
}

void multMatrices(const float *M1, const float *M2, float *res, int type) {
    switch (type) {
        case NO_OPT:
            multMatrices_NO_OPT(M1, M2, res);
            break;
        case VECTOR_OPT:
            multMatrices_VECTOR_OPT(M1, M2, res);
            break;
        case BLAS_OPT:
            multMatrices_BLAS_OPT(M1, M2, res);
            break;
    }
}

void sumMatrices(const float *M1, const float *M2, bool sign, float *res, int type) {
    if (sign == true) {
        switch (type) {
            case NO_OPT:
                sumMatrices_NO__OPT(M1, M2, res);
                break;
            case VECTOR_OPT:
                sumMatrices_VECTOR__OPT(M1, M2, res);
                break;
            case BLAS_OPT:
                sumMatrices_BLAS__OPT(M1, M2, res);
                break;
        }
    }
    else if (sign == false) {
        switch (type) {
            case NO_OPT:
                subMatrices_NO__OPT(M1, M2, res);
                break;
            case VECTOR_OPT:
                subMatrices_VECTOR__OPT(M1, M2, res);
                break;
            case BLAS_OPT:
                subMatrices_BLAS__OPT(M1, M2, res);
                break;
        }
    }
}