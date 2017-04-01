/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include <time.h>
#include "parallel_sort.h"

// implementation of your parallel sorting
void parallel_sort(int* begin, int* end, MPI_Comm comm) {

    // Terminating condition
	int commsize;
	int arrSize = sizeof(begin)/sizeof(begin[0]);
	MPI_Comm_size(comm, &commsize);
	if(commsize > 1){
		qsort(begin, arrSize, sizeof(int), cmpfunc);
	}
    // Call seeding helper function
    seed_rand(commsize);

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
    int le_size = 0, g_size = 0;
    int *greater = (int*) malloc(end - begin);
    

    // Allgather to find total # of elements < and > pivot
    MPI_Allgather(void* sbuf, int scount, MPI_INT,void* rbuf, int rcount, MPI_INT,comm);


    // Decide # of processors for < and > pivot


    // Send < and > arrays to appropriate processors (using alltoall)


    // Create two new communicators
    // MPI_Comm_split


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
void seed_rand(int commsize) {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Only seed once
    if (commsize == world_size) {
        // Each processor should have the same number for the current minute.
        // There's a small chance the seed will be different between processors
        // if the program runs in a very small window around when the minute
        // turns over. This is acceptably unlikely for now.
        srand(time(NULL) / 60);
    }
}

