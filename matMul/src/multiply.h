#include <cstddef>

void simple_multiply(double *mat1,
              double *mat2, std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O,
              double *res);

void multiply(double *mat1,
              double *mat2, std::ptrdiff_t N, std::ptrdiff_t M, std::ptrdiff_t O,
              double *res);
