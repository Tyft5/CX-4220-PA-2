/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include <stdlib.h>
#include <time.h>
#include "parallel_sort.h"



// implementation of your parallel sorting
void parallel_sort(int* begin, int* end, MPI_Comm comm) {

    // Terminating condition
	int commsize;
	MPI_Comm_size(comm, &commsize);
	if(commsize > 1){
		qsort(begin, (end - begin)/sizeof(int), sizeof(int), cmpfunc);
	}
    // Call seeding helper function
    seed_rand(comm);

    // Generate a pivot
    // Check for pivot in local array
    // If you have the pivot, broadcast it, otherwise receive it


    // Split local array based on pivot
    //  - allocate second array
    //  - keep track of size of both arrays


    // Allgather to find total # of elements < and > pivot


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
void seed_rand(MPI_Comm comm) {
    int local_size, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_size(comm, local_size);

    // Only seed once
    if (local_size == world_size) {
        // Each processor should have the same number for the current minute.
        // There's a small chance the seed will be different between processors
        // if the program runs in a very small window around when the minute
        // turns over. This is acceptably unlikely for now.
        srand(time(NULL) / 60);
    }
}

