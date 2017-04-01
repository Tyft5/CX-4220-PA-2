/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"

// implementation of your parallel sorting
void parallel_sort(int * begin, int* end, MPI_Comm comm) {

    // Terminating condition


    // Call seeding helper function


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



}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

// ...

// Function to seed RNG once with the same seed on each processor