#define TCAT2(a, b) a##_##b
#define CAT2(a, b) TCAT2(a, b)

void CAT2(SUMMATRIX_FUNCTION_NAME, NO__OPT) (const float *M1, const float *M2, float *res) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            #ifndef SUBSTRACT
                res[i*MATRIX_SIZE + j] = M1[i*MATRIX_SIZE + j] + M2[i*MATRIX_SIZE + j];
            #else
                res[i*MATRIX_SIZE + j] = M1[i*MATRIX_SIZE + j] - M2[i*MATRIX_SIZE + j];
            #endif
        }
    }
}

void CAT2(SUMMATRIX_FUNCTION_NAME, VECTOR__OPT) (const float *M1, const float *M2, float *res) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        __m128 result;
        for (int j = 0; j < MATRIX_SIZE_VECTOR_OPT; j += 4) {
            __m128 M1_row = _mm_loadu_ps(&M1[i*MATRIX_SIZE + j]);
            __m128 M2_row = _mm_loadu_ps(&M2[i*MATRIX_SIZE + j]);
            # ifndef SUBSTRACT
                result = _mm_add_ps(M1_row, M2_row);
            #else
                result = _mm_sub_ps(M1_row, M2_row);
            #endif
            _mm_storeu_ps(&res[i*MATRIX_SIZE+ j], result);
        }
    }
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = MATRIX_SIZE_VECTOR_OPT; j < MATRIX_SIZE; j++) {
            #ifndef SUBSTRACT
                res[i*MATRIX_SIZE + j] = M1[i*MATRIX_SIZE + j] + M2[i*MATRIX_SIZE + j];
            #else
                res[i*MATRIX_SIZE + j] = M1[i*MATRIX_SIZE + j] - M2[i*MATRIX_SIZE + j];
            #endif
        }
    }
}

void CAT2(SUMMATRIX_FUNCTION_NAME, BLAS__OPT) (const float *M1, const float *M2, float *res) {
        cblas_scopy(SQUARED_SIZE, M1, 1, res, 1);

        #ifndef SUBSTRACT
            cblas_saxpy(SQUARED_SIZE, 1, M2, 1, res, 1);
        #else
            cblas_saxpy(SQUARED_SIZE, -1, M2, 1, res, 1);
        #endif
}