#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iomanip>
//#pragma GCC target("avx2")

void simple_multiply(double *mat1,
              double *mat2, std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O,
              double *res) {
    for (std::ptrdiff_t i = 0; i < N; i++) {
        for (std::ptrdiff_t j = 0; j < M; j++) {
            for (std::ptrdiff_t k = 0; k < O; k++)
                res[i * M + j] += mat1[i * O + k] * mat2[k * M + j];
        }
    }
}

void sum(const double *A, std::ptrdiff_t A_row, const double *B, std::ptrdiff_t B_row, double *C, std::ptrdiff_t C_row,
         std::ptrdiff_t N, std::ptrdiff_t M) {
    for (std::ptrdiff_t i = 0; i < N; i++) {
        for (std::ptrdiff_t j = 0; j < M; j++) {
            C [i * C_row + j] = A [i * A_row + j] + B [i * B_row + j];
        }
    }
}

void sub(const double *A, std::ptrdiff_t A_row, const double *B, std::ptrdiff_t B_row, double *C, std::ptrdiff_t C_row,
         std::ptrdiff_t N, std::ptrdiff_t M) {
    for (std::ptrdiff_t i = 0; i < N; i++) {
        for (std::ptrdiff_t j = 0; j < M; j++) {
            C [i * C_row + j] = A [i * A_row + j] - B [i * B_row + j];
        }
    }
}

void GeneralMult(std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O, const double *A, std::ptrdiff_t A_row, const double *B, std::ptrdiff_t B_row, double *C, std::ptrdiff_t C_row) {
    for (std::ptrdiff_t i = 0; i < N; i++) {
        for (std::ptrdiff_t j = 0; j < M; j++) {
            C [i * C_row + j] = 0.0;  
            for (std::ptrdiff_t k = 0; k < O; k++)  
                C [i * C_row + j] += (double)(A [i * A_row + k] * B [k * B_row + j]);  
        } 
    }
}

#include <iomanip>

void printMatrix(const double *Mat, std::ptrdiff_t N, std::ptrdiff_t M) {
    for (std::ptrdiff_t i = 0; i < N; i++) {
        for (std::ptrdiff_t j = 0; j < M; j++) {
            std::cout << std::setprecision(17) << Mat[i * M + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

//  Strassen algorithm:  https://habr.com/ru/post/313258
 
//  S1  = A11 + A22
//  S2  = B11 + B22
//  P1  = S1 * S2 = (A11 + A22) * (B11 + B22)
//  S3  = A21 + A22
//  P2  = S3 * B11 = (A21 + A22) * B11
//  S4  = B12 - B22
//  P3  = A11 * S4 = A11 * (B12 - B22)
//  S5  = B21 - B11
//  P4  = A22 * S5 = A22 * (B21 - B11)
//  S6  = A11 + A12
//  P5  = S6 * B22 = (A11+A12) * B22
//  S7  = A21 - A11
//  S8  = B11 + B12
//  P6  = S7 * S8 = (A21 - A11) * (B11 + B12)
//  S9  = A12 - A22
//  S10 = B21 + B22
//  P7  = S9 * S10 = (A12 - A22) * (B21 + B22)
//  C11 = P1 + P4 - P5 + P7
//  C12 = P3 + P5
//  C21 = P2 + P4
//  C22 = P1 - P2 + P3 + P6

void StrassenMult(std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O, const double *A, const std::ptrdiff_t A_row,
                  const double *B, std::ptrdiff_t B_row, double *C, std::ptrdiff_t C_row) {

    if ( N == 1 || M == 1 || O == 1 /*|| ((unsigned long long)N * M * O < 512 * 512 * 512)*/ )
        GeneralMult(N, M, O, A, A_row, B, B_row, C, C_row); //if matrices are too small
    else {
        std::ptrdiff_t N_2 = N / 2;
        std::ptrdiff_t M_2 = M / 2;
        std::ptrdiff_t O_2 = O / 2;

        //A = N * O
        const double *A11 = &A [0]; 
        const double *A12 = &A [O_2];
        const double *A21 = &A [N_2 * A_row]; 
        const double *A22 = &A [N_2 * A_row + O_2];

        //B = O * M
        const double *B11 = &B [0];
        const double *B12 = &B [M_2];
        const double *B21 = &B [O_2 * B_row];
        const double *B22 = &B [O_2 * B_row + M_2];

        //C = N * M
        double *C11 = &C [0];
        double *C12 = &C [M_2];
        double *C21 = &C [N_2 * C_row];
        double *C22 = &C [N_2 * C_row + M_2];

        double *S1 = (double *)malloc(sizeof(double) * N_2 * O_2); //A
        double *S2 = (double *)malloc(sizeof(double) * O_2 * M_2); //B
        double *S3 = (double *)malloc(sizeof(double) * N_2 * O_2); //A
        double *S4 = (double *)malloc(sizeof(double) * O_2 * M_2); // B
        double *S5 = (double *)malloc(sizeof(double) * O_2 * M_2); // B
        double *S6 = (double *)malloc(sizeof(double) * N_2 * O_2); // A
        double *S7 = (double *)malloc(sizeof(double) * N_2 * O_2); // A
        double *S8 = (double *)malloc(sizeof(double) * O_2 * M_2); // B
        double *S9 = (double *)malloc(sizeof(double) * N_2 * O_2); // A
        double *S10 = (double *)malloc(sizeof(double) * O_2 * M_2); // B

        double *P1 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P2 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P3 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P4 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P5 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P6 = (double *)malloc(sizeof(double) * N_2 * M_2);
        double *P7 = (double *)malloc(sizeof(double) * N_2 * M_2);
        
        sum(A11, A_row, A22, A_row, S1, O_2, N_2, O_2); //S1
        sum(B11, B_row, B22, B_row, S2, M_2, O_2, M_2); //S2
        StrassenMult(N_2, M_2, O_2, S1, O_2, S2, M_2, P1, M_2); //P1  = S1 * S2

        sum(A21, A_row, A22, A_row, S3, O_2, N_2, O_2); //S3
        StrassenMult(N_2, M_2, O_2, S3, O_2, B11, B_row, P2, M_2); //P2  = S3 * B11

        sub(B12, B_row, B22, B_row, S4, M_2, O_2, M_2); //S4
        StrassenMult(N_2, M_2, O_2, A11, A_row, S4, M_2, P3, M_2); //P3  = A11 * S4

        sub(B21, B_row, B11, B_row, S5, M_2, O_2, M_2); //S5
        StrassenMult(N_2, M_2, O_2, A22, A_row, S5, M_2, P4, M_2); //P4  = A22 * S5

        sum(A11, A_row, A12, A_row, S6, O_2, N_2, O_2); //S6
        StrassenMult(N_2, M_2, O_2, S6, O_2, B22, B_row, P5, M_2); //P5  = S6 * B22

        sub(A21, A_row, A11, A_row, S7, O_2, N_2, O_2);
        sum(B11, B_row, B12, B_row, S8, M_2, O_2, M_2);
        StrassenMult(N_2, M_2, O_2, S7, O_2, S8, M_2, P6, M_2); //P6  = S7 * S8

        sub(A12, A_row, A22, A_row, S9, O_2, N_2, O_2); //S9
        sum(B21, B_row, B22, B_row, S10, M_2, O_2, M_2); //S10
        StrassenMult(N_2, M_2, O_2, S9, O_2, S10, M_2, P7, M_2); //P7  = S9 * S10

        for (std::ptrdiff_t i = 0; i < N_2; i++) {
            for (std::ptrdiff_t j = 0; j < M_2; j++) {
                C11 [i * C_row + j] = P1 [i * M_2 + j] + P4 [i * M_2 + j] - P5 [i * M_2 + j] + P7 [i * M_2 + j];
                C12 [i * C_row + j] = P3 [i * M_2 + j] + P5 [i * M_2 + j];
                C21 [i * C_row + j] = P2 [i * M_2 + j] + P4 [i * M_2 + j];
                C22 [i * C_row + j] = P1 [i * M_2 + j] - P2 [i * M_2 + j] + P3 [i * M_2 + j] + P6 [i * M_2 + j];
            }
        }

        free(S1); free(S2); free(S3); free(S4); free(S5); free(S6); free(S7); free(S8); free(S9); free(S10);
        free(P1); free(P2); free(P3); free(P4); free(P5); free(P6); free(P7);
    }
}

void multiply(double *mat1,
              double *mat2, std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O,
              double *res) {
    StrassenMult(N, M, O, mat1, O, mat2, M, res, M); //this algorithm is good only for big matrices
}

