#include <stdio.h>
#include <math.h>
#include <time.h> //-lrt

#define Pi 3.1415926536

int main() {
    struct timespec start, end; //time checkpoints

    long long int degree; //input
    scanf("%lld", &degree);

    degree %= 360;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double x;
    x = degree * (Pi / 180); //degree to radians

    double numerator = x; //числитель
    double denominator = 1; //знаменатель

    int i = 2; //for denominator
    double term, sum = 0; //term - слагаемое
    do {
        term = numerator / denominator;
        numerator = -numerator * x * x;
        denominator = denominator * i * (i+1);
        sum += term;
        i += 2;
    } while (fabs(term) >= 0.00001);

    /*for (int j = 0; j < 100000; j++) {
        for (int k = 0; k < 10000; k++) {}
    }*/

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    printf("answer: %lf\n", sum);
    printf("Time taken: %lld nsec.\n", 1000000000*(end.tv_sec-start.tv_sec) + (end.tv_nsec-start.tv_nsec));
    return 0;
}
