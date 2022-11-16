#include <stdio.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>

#define Pi 3.1415926536

int main() {

    struct tms start, end;
    long clocks_per_sec = sysconf(_SC_CLK_TCK);    
    long int clocks;

    long long int degree; //input
    scanf("%lld", &degree);

    times(&start);

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

    times(&end);

    printf("answer: %lf\n", sum);
    clocks = (long int)end.tms_utime - (long int)start.tms_utime;
    printf("Time taken: %lf sec.\n", (double)clocks / clocks_per_sec);
    return 0;
}
