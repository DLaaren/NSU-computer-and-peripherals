#include <iostream>
#include <climits>
#include <x86intrin.h>

#define ARRAY_SIZE 200000000
#define RUN_TIMES 4

void initArray(unsigned int *array, unsigned int fragments, size_t offset, size_t size) {
    size_t j = 1;
    size_t i = 0;

    while (i < size) {
        for(j = 1; j < fragments; j++) {
            array[i + (j - 1) * offset] = i + j * offset;
        }
        array[i + (j - 1) * offset] = i + 1;
        i++;
    }

    array[i - 1 + (j - 1) * offset] = 0;
}

unsigned long long runArray(unsigned int const *array) {
    unsigned long long startTime, endTime;
    unsigned long long minTime = ULLONG_MAX;

    for (size_t j = 0; j < RUN_TIMES; j++) {
        startTime = __rdtsc();

        for(size_t k = 0, i = 0; i < ARRAY_SIZE; i++){
            k = array[k];
        }

        endTime = __rdtsc();
        if (minTime > endTime - startTime)
            minTime = endTime - startTime;
    }
    return minTime;
}

void countTime(unsigned int *array, unsigned int fragments, int offset, int size) {

    initArray(array, fragments, offset, size);
    std::cout << fragments << " fragments  " << runArray(array) / ARRAY_SIZE << "  tacts" << std::endl;

}

// getconf -a | grep CACHE 
int main() {
    unsigned int *array = new unsigned int [ARRAY_SIZE];

    if (array == nullptr) {
        return 1;
    }

    unsigned int offset = 8*1024*1024; // 1MB
    unsigned int BlockSize = 32*8*1024; // 32KB ( L1 cache size )

    for(int fragments = 1; fragments <= 32; fragments++) {
        countTime(array, fragments, offset / sizeof(int), BlockSize / sizeof(int));
    }

    delete [] array;
    return 0;
}