#pragma once

#include <iostream>
#include <array>
#include <cstring>
#include <cmath>
#include <xmmintrin.h>

#include <assert.h>

#ifdef __cplusplus
extern "C" {
    #endif
    #include <cblas.h>
    #ifdef __cplusplus
}
#endif

#define MATRIX_SIZE 400
#define MATRIX_SIZE_VECTOR_OPT ((MATRIX_SIZE) - ((MATRIX_SIZE) % (4)))
#define iterations 10

#define NO_OPT 0
#define VECTOR_OPT 1
#define BLAS_OPT 2

#include "generator.h"

void printMatrix(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            std::cout << M[i*MATRIX_SIZE + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


std::array<float, MATRIX_SIZE*MATRIX_SIZE> generateMatrix() {
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> A = {0};
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i*MATRIX_SIZE + j] = std::rand() % 10;
        }
    }

    return A;
}

void createIdentityMatrix(std::array<float, MATRIX_SIZE*MATRIX_SIZE> &M);

void findMatrixB(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &A ,std::array<float, MATRIX_SIZE*MATRIX_SIZE> &B);

void multMatrices(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                  std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res, int type);

void sumMatrices(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2, bool sign, 
                 std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res, int type);


double matrixInverse(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &A,
                     std::array<float, MATRIX_SIZE*MATRIX_SIZE> &resMatrix, int type) {
    clock_t start = clock();
    
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> I = {0};
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> R = {0}; //R = I - B*A
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> B = {0};
    std::array<float, MATRIX_SIZE*MATRIX_SIZE> tmp = {0};

    createIdentityMatrix(I);
    findMatrixB(A, B);    
    multMatrices(B, A, tmp, type); //B*A = tmp

    sumMatrices(I, tmp, false, R, type); // I - tmp = R
    sumMatrices(I, R, true, resMatrix, type); //I + R = resMatrix

    multMatrices(R, R, tmp, type); //tmp = R^2

    for (int k = 2; k < iterations - 1; k++) {
        sumMatrices(resMatrix, tmp, true, resMatrix, type); //resMatrix + R^k
        multMatrices(tmp, R, tmp, type); //R^(k+1)
    }
    sumMatrices(resMatrix, tmp, true, resMatrix, type); //last addition (resMatrix + R^9)

    multMatrices(resMatrix, B, resMatrix, type); // (I + R + ... + R^9) * B

    clock_t end = clock();
    double time = 1000.0 * ((double)(end - start)/CLOCKS_PER_SEC);
    return time;
}

void createIdentityMatrix(std::array<float, MATRIX_SIZE*MATRIX_SIZE> &I) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        I[i*MATRIX_SIZE + i] = 1;
    }
}

void findMatrixB(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> &A, std::array<float, MATRIX_SIZE*MATRIX_SIZE> &B) {
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

void multMatrices_NO_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                         std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res);

void multMatrices_VECTOR_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, 
                             const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2, std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res);

void multMatrices_BLAS_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                           std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res);

void multMatrices(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                  std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res, int type) {
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
void multMatrices_NO_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                         std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
           float tmpSum = 0;
           for (int k = 0; k < MATRIX_SIZE; k++) {
                tmpSum = tmpSum + M1[i*MATRIX_SIZE + k] * M2[k*MATRIX_SIZE + j];
           }
           res[i*MATRIX_SIZE + j] = tmpSum;
        }
    }
}

void multMatrices_VECTOR_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                             std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res) {
    static std::array<float, MATRIX_SIZE*MATRIX_SIZE> tmp;
    std::fill(tmp.begin(), tmp.end(), 0);
    float tmpVector[4];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int k = 0; k < MATRIX_SIZE; k++) {
            __m128 M1_scalar = _mm_load_ps1(&M1[i*MATRIX_SIZE + k]);
            for (int j = 0; j < MATRIX_SIZE_VECTOR_OPT; j += 4) {
                __m128 M2_row = _mm_loadu_ps(&M2[k*MATRIX_SIZE + j]);
                __m128 result;
                result = _mm_mul_ps(M1_scalar, M2_row);
                _mm_storeu_ps(tmpVector, result);
                for (int m = 0, n = j; m < 4; m++, n++) {
                    tmp[i*MATRIX_SIZE + n] += tmpVector[m];
                }
            }
        }
    }
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int k = 0; k < MATRIX_SIZE; k++) {
            float tmpSum = 0;
            for (int j = MATRIX_SIZE_VECTOR_OPT; j < MATRIX_SIZE; j++) {
                float scalar = M1[i*MATRIX_SIZE + k];
                tmpSum = M2[k*MATRIX_SIZE + j] * scalar;
                tmp[i*MATRIX_SIZE + j] += tmpSum;
            }
        }
    }
    res = tmp;
}

void multMatrices_BLAS_OPT(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2,
                           std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res) {
    //https://www.ibm.com/docs/en/essl/6.3?topic=mos-sgemm-dgemm-cgemm-zgemm-combined-matrix-multiplication-addition-general-matrices-their-transposes-conjugate-transposes
    const float *M1_p = &M1[0];
    const float *M2_p = &M2[0];
    float *res_p = &res[0];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE, 1.0, 
                M1_p, MATRIX_SIZE, M2_p, MATRIX_SIZE, 0.0, res_p, MATRIX_SIZE);
}

void sumMatrices(const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M1, const std::array<float, MATRIX_SIZE*MATRIX_SIZE> M2, bool sign,
                 std::array<float, MATRIX_SIZE*MATRIX_SIZE> &res, int type) {
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