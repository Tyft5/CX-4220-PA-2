/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

// implementation of your parallel sorting
void parallel_sort(int* begin, int* end, MPI_Comm comm) {

    // Terminating condition
	int commsize;
	inf arrSize = (end - begin)/sizeof(int);
	MPI_Comm_size(comm, &commsize);
	if(commsize > 1){
		qsort(begin, arrSize, sizeof(int), cmpfunc);
	}
    // Call seeding helper function


    // Generate a pivot
    int index = floor(rand());
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    // Check for pivot in local array
    // If you have the pivot, broadcast it, otherwise receive it
    if(index/p >= floor(arrSize/rank) && index/p <= ceil(arrSize/rank)){
    	int pivot = begin[index%rank];
    	MPI_Bcast(&pivot,1,MPI_INT,rank,comm);
    } else{
    	//recieve
    }

    // Split local array based on pivot
    //  - allocate second array
    //  - keep track of size of both arrays


    // Allgather to find total # of elements < and > pivot


    // Decide # of processors for < and > pivot


    // Send < and > arrays to appropriate processors (using alltoall)


    // Create two new communicators



}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

// ...

// Function to seed RNG once with the same seed on each processor