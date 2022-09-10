#include <stdio.h>
#include <math.h>
#include <x86intrin.h> // -msse2

#define Pi 3.1415926536

int main() {

    unsigned int start, end;

    long long int cpu_Hz = 2300000000;     //3085.922


    long long int degree; //input
    scanf("%lld", &degree);

    start = __rdtsc();

    degree %= 360;

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

    end = __rdtsc();

    printf("answer: %lf\n", sum);
    printf("Time taken: %lld sec.\n", (end - start) / cpu_Hz);
    printf("ticks: %u\n", end - start);
    return 0;
}