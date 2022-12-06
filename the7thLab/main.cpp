#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
    #endif
    #include <cblas.h>
    #ifdef __cplusplus
}
#endif

#define MATRIX_SIZE 4
#define iterations 10

#define NO_OPT 0
#define VECTOR_OPT 1
#define BLAS_OPT 2

/*tasks:
1) Убрать if'ы из for'ов 
2) Найти blas сложение векторов  //done//
3) Убрать транспонирование матрицы в векторизации и векторизовать то, что будет //done//
4) Сделать тесты на корректность разложения
5) Сделать проверку на существование обратной матрицы
*/

#include "generator.h"

double matrixInverse(std::vector<float>, std::vector<float>&, int type);
void printMatrix(const std::vector<float>&);
std::vector<float> generateMatrix();

int main() {
    std::vector<float> A(MATRIX_SIZE*MATRIX_SIZE);
    std::vector<float> backMatrix(MATRIX_SIZE*MATRIX_SIZE, 0); //A^(-1)

    A = generateMatrix();

    int type = NO_OPT;
    double res;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "no optimization: " << res << " ms\n";

    printMatrix(backMatrix);

    std::fill(backMatrix.begin(), backMatrix.end(), 0);

    type = VECTOR_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "vector optimization: " << res << " ms\n";

    printMatrix(backMatrix);

    std::fill(backMatrix.begin(), backMatrix.end(), 0);

    type = BLAS_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "blas optimization: " << res << " ms\n";

    printMatrix(backMatrix);

    return 0;
}

void printMatrix(const std::vector<float> &M) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            std::cout << M[i*MATRIX_SIZE + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


std::vector<float> generateMatrix() {
    std::vector<float> A(MATRIX_SIZE*MATRIX_SIZE);
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            A[i*MATRIX_SIZE + j] = std::rand() % 10;
        }
    }

    return A;
}

void createIdentityMatrix(std::vector<float> &M);
void findMatrixB(const std::vector<float> A ,std::vector<float> &B);
void multMatrices(const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res, int type);
void sumMatrices(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res, int type);

double matrixInverse(std::vector<float> A, std::vector<float> &resMatrix, int type) {
    clock_t start = clock();
    
    std::vector<float> I(MATRIX_SIZE*MATRIX_SIZE, 0);
    std::vector<float> R(MATRIX_SIZE*MATRIX_SIZE, 0); //R = I - B*A
    std::vector<float> B(MATRIX_SIZE*MATRIX_SIZE, 0);
    std::vector<float> tmp(MATRIX_SIZE*MATRIX_SIZE, 0);

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

void createIdentityMatrix(std::vector<float> &I) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        I[i*MATRIX_SIZE + i] = 1;
    }
}

void findMatrixB(const std::vector<float> A, std::vector<float> &B) {
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

void multMatrices_NO_OPT(const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res);
void multMatrices_VECTOR_OPT(std::vector<float> M1, std::vector<float> M2, std::vector<float> &res);
void multMatrices_BLAS_OPT(std::vector<float> M1, std::vector<float> M2, std::vector<float> &res);

void multMatrices(const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res, int type) {
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
void multMatrices_NO_OPT(const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res) {
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

void multMatrices_VECTOR_OPT(std::vector<float> M1, std::vector<float> M2, std::vector<float> &res) {
    std::vector<float> tmp(MATRIX_SIZE*MATRIX_SIZE);
    float tmpVector[4];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int k = 0; k < MATRIX_SIZE; k++) {
            __m128 M1_scalar = _mm_load_ps1(&M1[i*MATRIX_SIZE + k]);
            for (int j = 0; j < MATRIX_SIZE; j += 4) {
                __m128 M2_row = _mm_load_ps(&M2[k*MATRIX_SIZE + j]);
                __m128 result;
                result = _mm_mul_ps(M1_scalar, M2_row);
                _mm_store_ps(tmpVector, result);
                for (int m = 0, n = j; m < 4; m++, n++) {
                    tmp[i*MATRIX_SIZE + n] += tmpVector[m];
                }
            }
        }
    }
    res = tmp;
}

void multMatrices_BLAS_OPT(std::vector<float> M1, std::vector<float> M2, std::vector<float> &res) {
    //https://www.ibm.com/docs/en/essl/6.3?topic=mos-sgemm-dgemm-cgemm-zgemm-combined-matrix-multiplication-addition-general-matrices-their-transposes-conjugate-transposes
    float *M1_p = &M1[0];
    float *M2_p = &M2[0];
    float *res_p = &res[0];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE, 1.0, M1_p, MATRIX_SIZE, M2_p, MATRIX_SIZE, 0.0, res_p, MATRIX_SIZE);
}

void sumMatrices(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res, int type) {
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