/*
*  This file is part of Christian's OpenMP software lab 
*
*  Copyright (C) 2016 by Christian Terboven <terboven@itc.rwth-aachen.de>
*  Copyright (C) 2016 by Jonas Hahnfeld <hahnfeld@itc.rwth-aachen.de>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>
#include <omp.h>

#include <iostream>
#include <algorithm>

#include <cstdlib>
#include <cstdio>

#include <cmath>
#include <ctime>
#include <cstring>

#define CUT_OFF 1000000

/**
  * helper routine: check if array is sorted correctly
  */
bool isSorted(int ref[], int data[], const size_t size) {
    std::sort(ref, ref + size);
    for (size_t idx = 0; idx < size; ++idx){
        if (ref[idx] != data[idx]) {
            return false;
        }
    }
    return true;
}


/*
 * helper routine: "lower bound" binary search
 * */
long binarySearch(int *a, int el, long begin, long end) {
    long s = begin;
    long e = end;

    while(s <= (e-1)) {
        long pivot = (s+e)/2;

        if(el < a[pivot]) {
            e = pivot;
        } else {
            s = pivot+1;
        }
    }
    
    return s;
}

/**
  * sequential merge step (straight-forward implementation)
  */
// TODO: cut-off could also apply here (extra parameter?)
// TODO: optional: we can also break merge in two halves
void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
    long left = begin1;
    long right = begin2;

    long idx = outBegin;

    while (left < end1 && right < end2) {
        if (in[left] <= in[right]) {
            out[idx] = in[left];
            left++;
        } else {
            out[idx] = in[right];
            right++;
        }
        idx++;
    }

    while (left < end1) {
        out[idx] = in[left];
        left++, idx++;
    }

    while (right < end2) {
        out[idx] = in[right];
        right++, idx++;
    }
}

void MsMergeParallel(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
    if((end1-begin1) + (end2-begin2) <= CUT_OFF) {
        MsMergeSequential(out, in, begin1, end1, begin2, end2, outBegin);
    } else {
        long pivot1 = begin1+(end1-begin1)/2;
        long pivot2 = binarySearch(in, in[pivot1], begin2, end2);
        long pivot3 = outBegin+(pivot1-begin1)+(pivot2-begin2);

        #pragma omp task firstprivate(begin1, pivot1, begin2, pivot2, outBegin)
        MsMergeParallel(out, in, begin1, pivot1, begin2, pivot2, outBegin);

        MsMergeParallel(out, in, pivot1, end1, pivot2, end2, pivot3);

        #pragma omp taskwait
    }
}

/**
  * Sequential MergeSort
  */
// OK: remember one additional parameter (depth)
// OK: recursive calls could be taskyfied
// OK: task synchronization also is required
void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
    if (begin < (end - 1)) {
        const long half = (begin + end) / 2;
        MsSequential(array, tmp, !inplace, begin, half);
        MsSequential(array, tmp, !inplace, half, end);
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }
    } else if (!inplace) {
        tmp[begin] = array[begin];
    }
}

/*
 * Parallel MergeSort
 * */
void MsParallel_(int *array, int *tmp, bool inplace, long begin, long end, long depth) {
    // When the sub-arrays are "small enough" don't create more tasks,
    //  instead, use the sequential mergesort
    //  => less overhead
    if(depth <= 0) {
        MsSequential(array, tmp, inplace, begin, end);
    } else {
        if (begin < (end - 1)) {
            const long half = (begin + end) / 2;
            // Create a task to be executed by an available thread
            #pragma omp task firstprivate(inplace, begin, half, depth)
            MsParallel_(array, tmp, !inplace, begin, half, depth-1);

            // Start executing another MsParellel_ instance as soon as
            //  the previous task is created.
            MsParallel_(array, tmp, !inplace, half, end, depth-1);

            // Wait for the previous task to complete execution
            #pragma omp taskwait

            // Merge the 2 ordered sub-arrays into a single array
            if (inplace) {
                // MsMergeSequential(array, tmp, begin, half, half, end, begin);
                MsMergeParallel(array, tmp, begin, half, half, end, begin);
            } else {
                // MsMergeSequential(tmp, array, begin, half, half, end, begin);
                MsMergeParallel(tmp, array, begin, half, half, end, begin);
            }
        } else if (!inplace) {
            tmp[begin] = array[begin];
        }
    }
}



/**
  * Serial MergeSort
  */
// OK: this function should create the parallel region
// OK: good point to compute a good depth level (cut-off)
void MsSerial(int *array, int *tmp, const size_t size) {

    // OK: parallel version of MsSequential will receive one more parameter: 'depth' (used as cut-off)
    MsSequential(array, tmp, true, 0, size);
}

/*
 * Parallel MergeSort
 * */
void MsParallel(int *array, int *tmp, const size_t size) {
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Compute a depth level based on the number of available threads
            //  and the dimension of the array
            long depth = std::log2(omp_get_num_threads());

            MsParallel_(array, tmp, true, 0, size, depth);
        }
    }
}



/** 
  * @brief program entry point
  */
int main(int argc, char* argv[]) {
    // variables to measure the elapsed time
    struct timeval t1, t2;
    double etime;

    // expect one command line arguments: array size
    if (argc != 2) {
        printf("Usage: MergeSort.exe <array size> \n");
        printf("\n");
        return EXIT_FAILURE;
    }
    else {
        const size_t stSize = strtol(argv[1], NULL, 10);
        int *data = (int*) malloc(stSize * sizeof(int));
        int *tmp = (int*) malloc(stSize * sizeof(int));
        int *ref = (int*) malloc(stSize * sizeof(int));

        printf("Initialization...\n");

        // Initialize "data" with random values
        srand(95);
        for (size_t idx = 0; idx < stSize; ++idx){
            data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
        }
        std::copy(data, data + stSize, ref);

        double dSize = (stSize * sizeof(int)) / 1024 / 1024;
        printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

        gettimeofday(&t1, NULL);
        // MsSerial(data, tmp, stSize);
        MsParallel(data, tmp, stSize);
        gettimeofday(&t2, NULL);

        etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
        etime = etime / 1000;

        printf("Time: %f\n", etime);
        printf("Verification... ");
        if (isSorted(ref, data, stSize)) {
            printf(" successful.\n");
        }
        else {
            printf(" FAILED.\n");
        }

        free(data);
        free(tmp);
        free(ref);
    }

    return EXIT_SUCCESS;
}
