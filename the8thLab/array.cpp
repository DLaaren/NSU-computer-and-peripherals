#pragma O1
#include <x86intrin.h>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <limits.h>

#define ITERATIONS 10

void generateArray(unsigned int *array, size_t arraySize) {
    for (size_t i = 0; i < arraySize; i++) {
        array[i] = i+1;
    }
    array[arraySize - 1] = 0;
}

void generateReverseArray(unsigned int *array, size_t arraySize) {
    for (size_t i = arraySize - 1; i > 0; i--) {
        array[i] = i-1;
    }
    array[0] = arraySize - 1;
}

void test(unsigned int *array, size_t arraySize) {
    bool *tmp = new bool [arraySize];
    int m = 0;
    for (int i = 0; i < arraySize; i++) {
        tmp [ array[i] ] = true;
    }
    for (int i = 0; i < arraySize; i++) {
        assert((tmp[i] == true) && "incorrect random array\n");
    }  

    delete [] tmp;
}

void generateRandomArray(unsigned int *array, size_t arraySize) {

    for (long long i = 0; i < arraySize; i++) {
        array[i] = i;
    }
    srand(time(NULL));

    for (long long i = arraySize - 1; i > 0; i--) {
        unsigned int j = 0;
        j = rand() % i;
        std::swap(array[i], array[j]);
    }

    test(array, arraySize);
}

void run(unsigned int *array, size_t arraySize) {
    int m = 0;
    for (int i = 0; i < arraySize * ITERATIONS; i++) {
        m = array[m];
    }
}

unsigned int nextSize(unsigned int currentSize) {
    if (currentSize < ((25 * 1024) / 4) ) {
        currentSize *= 1.5;
        return currentSize;
    }
    else if (currentSize >= ((25 * 1024) / 4) && currentSize <= ((50 * 1024) / 4)) { //25...40 Kib
        currentSize += 200;
        return currentSize;
    }
    else if (currentSize > ((50 * 1024) / 4) && currentSize < ((100 * 1024) / 4) ) {
        currentSize *= 1.5;
        return currentSize;
    }   
    else if (currentSize >= ((100 * 1024) / 4) && currentSize <= ((280 * 1024) / 4)) { //100...280 Kib
        currentSize += 1500;
        return currentSize;
    }
    else if (currentSize > ((280 * 1024) / 4) && currentSize < ((7 * 1024 * 1024) / 4) ) {
        currentSize *= 1.1;
        return currentSize;
    }
    else if (currentSize >= ((7 * 1024 * 1024) / 4) && currentSize <= ((8.4 * 1024 * 1024) / 4)) { //7 ... 8.4 Mib
        currentSize += 1024 * 10;
        return currentSize;
    }
    else if (currentSize > ((8.4 * 1024 * 1024) / 4)) {
        currentSize *= 1.2;
        return currentSize;
    }
    return currentSize;
}

int main() {

    //lscpu 
    //L1 32 Kib (per core)
    //L2 256 Kib (per core)
    //L3 8 Mib (for all cores)

    unsigned int minSize = 256; //256 * int = 1024 = 1Kib
    unsigned int maxSize = 8 * 1024 * 1024; // 8Mib * int = 32 Mib
    unsigned int currentSize;

    currentSize = minSize;
    while (currentSize < maxSize) {
        unsigned int *array = new unsigned int [currentSize] ;
        generateArray(array, currentSize);
        run(array, currentSize);
        unsigned long long start = __rdtsc();
        run(array, currentSize);
        unsigned long long end = __rdtsc();
        unsigned long long time = end - start;
        delete [] array;
        std::cout << "forward: " << (currentSize * 4) / 1024 << "Kib; " << time / (ITERATIONS * currentSize) << " ticks\n";

        currentSize = nextSize(currentSize);
    }

    currentSize = minSize;
    while (currentSize < maxSize) {
        unsigned int *array = new unsigned int [currentSize] ;
        generateReverseArray(array, currentSize);
        run(array, currentSize);
        unsigned long long start = __rdtsc();
        run(array, currentSize);
        unsigned long long end = __rdtsc();
        unsigned long long time = end - start;
        delete [] array;
        std::cout << "reverse: " << (currentSize * 4) / 1024 << "Kib; " << time / (ITERATIONS * currentSize) << " ticks\n";

        currentSize = nextSize(currentSize);
    }

    currentSize = minSize;
    while (currentSize < maxSize) {
        unsigned int *array = new unsigned int [currentSize] ;
        generateRandomArray(array, currentSize);
        run(array, currentSize);
        unsigned long long start = __rdtsc();
        run(array, currentSize);
        unsigned long long end = __rdtsc();
        unsigned long long time = end - start;
        delete [] array;
        std::cout << "random: " << (currentSize * 4) / 1024 << "Kib; " << time / (ITERATIONS * currentSize) << " ticks\n";

        currentSize = nextSize(currentSize);
    }
    return 0;
}