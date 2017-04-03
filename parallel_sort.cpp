/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include <time.h>
#include <math.h>
#include "parallel_sort.h"

// implementation of your parallel sorting
void parallel_sort(int* begin, int* end, MPI_Comm comm) {

    int worldsize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);

    // Terminating condition
	int commsize;
	int arrSize = sizeof(begin)/sizeof(begin[0]);
	MPI_Comm_size(comm, &commsize);
	if(commsize > 1){
		qsort(begin, arrSize, sizeof(int), cmpfunc);
	}
    // Call seeding helper function
    seed_rand(commsize, worldsize);

    // Generate a pivot
    int index = floor(rand()/RAND_MAX*arrSize);
    int p, rank, pivot;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    // Check for pivot in local array
    // If you have the pivot, broadcast it, otherwise receive it
    int source = floor(index/p);
    if(rank == source){
    	pivot = begin[index%p];
    }
    MPI_Bcast(&pivot,1,MPI_INT,source,comm);


    // Split local array based on pivot
    //  - allocate second array
    //  - keep track of size of both arrays
    //  - realloc
    int num, le_size = 0, g_size = 0;
    int *greater = (int*) malloc(end - begin);
    for (int i = 0; i < arrSize; i++) {
        num = begin[i];
        if (num <= pivot) {
            begin[le_size++] = num;
        } else {
            greater[g_size++] = num;
        }
    }
    begin = (int*) realloc(begin, sizeof(int) * le_size);
    greater = (int*) realloc(greater, sizeof(int) * g_size);
    

    // Allgather to find total # of elements < and > pivot
    int* small = (int*) malloc(sizeof(int) * p);
    int* big = (int*) malloc(sizeof(int) * p);
    MPI_Allgather(&le_size, 1, MPI_INT, small, 1, MPI_INT,comm);
    MPI_Allgather(&g_size, 1, MPI_INT, big, 1, MPI_INT,comm);
    int smallsum = 0, bigsum = 0;
    for(int i = 0; i < p; i++){
    	smallsum += small[i];
    	bigsum += big[i];
    }

    // Decide # of processors for < and > pivot
    int l_proc_num = floor(worldsize * smallsum / (smallsum + bigsum));
    int g_proc_num = commsize - l_proc_num;

    // Send < and > arrays to appropriate processors (using alltoall)
    // Receive array for next recursion


    // Dealloc greater
    if (rank < l_proc_num) {
        free(greater);
    } else {
        free(begin);
    }

    // Create two new communicators
    // MPI_Comm_split
    // Dealloc old comm


    // Call self recursively



}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

// Function to seed RNG once with the same seed on each processor
void seed_rand(int commsize, int worldsize) {
    // Only seed once
    if (commsize == worldsize) {
        // Each processor should have the same number for the current minute.
        // There's a small chance the seed will be different between processors
        // if the program runs in a very small window around when the minute
        // turns over. This is acceptably unlikely for now.
        srand(time(NULL) / 60);
    }
}

