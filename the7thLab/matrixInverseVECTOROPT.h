#pragma once

void sumMatrices_VECTOR__OPT(const float *M1, const float *M2, float *res);
void subMatrices_VECTOR__OPT(const float *M1, const float *M2, float *res);

void multMatrices_VECTOR_OPT(const float *M1, const float *M2, float *res) {

    auto tmp = new float [SQUARED_SIZE];
    memset(tmp, 0, sizeof(float) * SQUARED_SIZE);

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int k = 0; k < MATRIX_SIZE; k++) {
            __m128 M1_scalar = _mm_load_ps1(&M1[i*MATRIX_SIZE + k]);
            for (int j = 0; j < MATRIX_SIZE_VECTOR_OPT; j += 4) {
                __m128 M2_row = _mm_loadu_ps(&M2[k*MATRIX_SIZE + j]);
                __m128 result = _mm_mul_ps(M1_scalar, M2_row);
                tmp[i*MATRIX_SIZE + j] += result[0];
                tmp[i*MATRIX_SIZE + j+1] += result[1];
                tmp[i*MATRIX_SIZE + j+2] += result[2];
                tmp[i*MATRIX_SIZE + j+3] += result[3];
            }
        }
    }

    if (MATRIX_SIZE != MATRIX_SIZE_VECTOR_OPT) {
        for (int i = 0; i < MATRIX_SIZE; i++) {
            for (int k = 0; k < MATRIX_SIZE; k++) {
                for (int j = MATRIX_SIZE_VECTOR_OPT; j < MATRIX_SIZE; j++) {
                    tmp[i*MATRIX_SIZE + j] += M2[k*MATRIX_SIZE + j] * M1[i*MATRIX_SIZE + k];
                }
            }
        }
    }
    memcpy(res, tmp, sizeof(float) * SQUARED_SIZE);
    delete [] tmp;
}