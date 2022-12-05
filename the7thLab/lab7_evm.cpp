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

/*tasks:
1) Убрать if'ы из for'ов 
2) Найти blas сложение векторов  //done//
3) Убрать транспонирование матрицы в векторизации и векторизовать то, что будет //done//
4) Сделать тесты на корректность разложения 
*/


#define N 1024
//for correct work it should be multiple of 4
#define iterations 10
//it should be 500 and over to be accurate as much as possible

#define NO_OPT 0
#define VECTOR_OPT 1
#define BLAS_OPT 2

double matrixInverse(std::vector<float>, std::vector<float>&, int type);
void printMatrix(const std::vector<float>&);

int main() {
    std::vector<float> A(N*N);
    std::vector<float> backMatrix(N*N, 0); //A^(-1)

    //creates random matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i*N + j] = std::rand() % 10;
        }
    }

    int type = NO_OPT;
    double res;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "no optimization: " << res << " ms\n";

    //printMatrix(backMatrix);

    std::fill(backMatrix.begin(), backMatrix.end(), 0);

    type = VECTOR_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "vector optimization: " << res << " ms\n";

    //printMatrix(backMatrix);

    std::fill(backMatrix.begin(), backMatrix.end(), 0);

    type = BLAS_OPT;
    res = matrixInverse(A, backMatrix, type); //returns time 
    std::cout << "blas optimization: " << res << " ms\n";

    //printMatrix(backMatrix);

    return 0;
}

void printMatrix(const std::vector<float> &M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << M[i*N + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void createIdentityMatrix(std::vector<float> &M);
void findMatrixB(const std::vector<float> A ,std::vector<float> &B);
void multMatrices(const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res, int type);
void sumMatrices(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res, int type);

double matrixInverse(std::vector<float> A, std::vector<float> &resMatrix, int type) {
    clock_t start = clock();
    
    std::vector<float> I(N*N, 0);
    std::vector<float> R(N*N, 0); //R = I - B*A
    std::vector<float> B(N*N, 0);
    std::vector<float> tmp(N*N, 0);

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
    for (int i = 0; i < N; i++) {
        I[i*N + i] = 1;
    }
}

void findMatrixB(const std::vector<float> A, std::vector<float> &B) {
    //||A||1 = максимальный по сумме элементов столбец -> maxSumJ
    //||A||inf = максимальная по сумме элементов строчка -> maxSumI
    int maxSumI = 0, maxSumJ = 0;
    int tmpSumI = 0, tmpSumJ = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            tmpSumI += (float)fabs(A[i*N + j]);
            tmpSumJ += (float)fabs(A[j*N + i]);
        }
        if (tmpSumI > maxSumI) maxSumI = tmpSumI;
        if (tmpSumJ > maxSumJ) maxSumJ = tmpSumJ;
        tmpSumI = 0;
        tmpSumJ = 0;
    }

    float divider = maxSumI * maxSumJ;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i*N +j] = A[j*N + i] / divider;
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
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
           float tmpSum = 0;
           for (int k = 0; k < N; k++) {
                tmpSum = tmpSum + M1[i*N + k] * M2[k*N + j];
           }
           res[i*N + j] = tmpSum;
        }
    }
}

void multMatrices_VECTOR_OPT(std::vector<float> M1, std::vector<float> M2, std::vector<float> &res) {
    std::vector<float> tmp(N*N);
    float tmpVector[4];
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            __m128 M1_scalar = _mm_load_ps1(&M1[i*N + k]);
            for (int j = 0; j < N; j += 4) {
                __m128 M2_row = _mm_load_ps(&M2[k*N + j]);
                __m128 result;
                result = _mm_mul_ps(M1_scalar, M2_row);
                _mm_store_ps(tmpVector, result);
                for (int m = 0, n = j; m < 4; m++, n++) {
                    tmp[i*N + n] += tmpVector[m];
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
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, M1_p, N, M2_p, N, 0.0, res_p, N);
}

void sumMatrices_NO_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res);
void sumMatrices_VECTOR_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res);
void sumMatrices_BLAS_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res);

void sumMatrices(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res, int type) {
    switch (type) {
        case NO_OPT:
            sumMatrices_NO_OPT(M1, M2, sign, res);
            break;
        case VECTOR_OPT:
            sumMatrices_VECTOR_OPT(M1, M2, sign, res);
            break;
        case BLAS_OPT:
            sumMatrices_BLAS_OPT(M1, M2, sign, res);
            break;
    }
}

void sumMatrices_NO_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res) {
    if (sign == true) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                res[i*N + j] = M1[i*N + j] + M2[i*N + j];
            }
        }
    }
    else if (sign == false) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                res[i*N + j] = M1[i*N + j] - M2[i*N + j];
            }
        }
    }
}

void sumMatrices_VECTOR_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res) {
    for (int i = 0; i < N; i++) {
        __m128 result;
        for (int j = 0; j < N; j += 4) {
            __m128 M1_row = _mm_load_ps(&M1[i*N + j]);
            __m128 M2_row = _mm_load_ps(&M2[i*N + j]);
            if (sign == true) {
                result = _mm_add_ps(M1_row, M2_row);
            }
            else if (sign == false) {
                result = _mm_sub_ps(M1_row, M2_row);
            }
            _mm_store_ps(&res[i*N + j], result);
        }
    }
}

void sumMatrices_BLAS_OPT(const std::vector<float> M1, const std::vector<float> M2, bool sign, std::vector<float> &res) {
    for (int i = 0; i < N; i++) {
        const float *M1_p = &M1[i*N];
        const float *M2_p = &M2[i*N];
        float *res_p = &res[i*N];
        cblas_scopy(N, M1_p, 1, res_p, 1);
        if (sign) cblas_saxpy(N, 1, M2_p, 1, res_p, 1);
        if (!sign) cblas_saxpy(N, -1, M2_p, 1, res_p, 1);
    }
}
