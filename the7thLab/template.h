
#define CONVERT(c) c

#define TCAT2(a, b) a##_##b
#define CAT2(a, b) TCAT2(a, b)

void CAT2(SUMMATRIX_FUNCTION_NAME, NO__OPT) (const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res) {
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

void CAT2(SUMMATRIX_FUNCTION_NAME, VECTOR__OPT) (const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        __m128 result;
        for (int j = 0; j < MATRIX_SIZE; j += 4) {
            __m128 M1_row = _mm_load_ps(&M1[i*MATRIX_SIZE + j]);
            __m128 M2_row = _mm_load_ps(&M2[i*MATRIX_SIZE + j]);
            #ifndef SUBSTRACT
                result = _mm_add_ps(M1_row, M2_row);
            #else
                result = _mm_sub_ps(M1_row, M2_row);
            #endif
            _mm_store_ps(&res[i*MATRIX_SIZE + j], result);
        }
    }
}

void CAT2(SUMMATRIX_FUNCTION_NAME, BLAS__OPT) (const std::vector<float> M1, const std::vector<float> M2, std::vector<float> &res) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        const float *M1_p = &M1[i*MATRIX_SIZE];
        const float *M2_p = &M2[i*MATRIX_SIZE];
        float *res_p = &res[i*MATRIX_SIZE];
        cblas_scopy(MATRIX_SIZE, M1_p, 1, res_p, 1);

        #ifndef SUBSTRACT
            cblas_saxpy(MATRIX_SIZE, 1, M2_p, 1, res_p, 1);
        #else
            cblas_saxpy(MATRIX_SIZE, -1, M2_p, 1, res_p, 1);
        #endif
    }
}

#undef CONVERT