#pragma once

#ifdef __cplusplus
extern "C" {
    #endif
    #include <cblas.h>
    #ifdef __cplusplus
}
#endif

void sumMatrices_BLAS__OPT(const float *M1, const float *M2, float *res);
void subMatrices_BLAS__OPT(const float *M1, const float *M2, float *res);

void multMatrices_BLAS_OPT(const float *M1, const float *M2, float *res) {
    //https://www.ibm.com/docs/en/essl/6.3?topic=mos-sgemm-dgemm-cgemm-zgemm-combined-matrix-multiplication-addition-general-matrices-their-transposes-conjugate-transposes
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE, 1.0, 
                M1, MATRIX_SIZE, M2, MATRIX_SIZE, 0.0, res, MATRIX_SIZE);
}