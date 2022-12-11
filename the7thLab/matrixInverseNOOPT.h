#pragma once

void sumMatrices_NO__OPT(const float *M1, const float *M2, float *res);
void subMatrices_NO__OPT(const float *M1, const float *M2, float *res);

void multMatrices_NO_OPT(const float *M1, const float *M2, float *res) {
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