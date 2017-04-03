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
    int p, rank, pivot;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(comm, &rank);

    // Terminating condition
	int commsize;
	int arrSize = (end - begin);
	MPI_Comm_size(comm, &commsize);
	if(commsize == 1){
		qsort(begin, arrSize, sizeof(int), cmpfunc);
        return;
	}
    // Call seeding helper function
    seed_rand(commsize, worldsize);

    // Generate a pivot
    int index = (int)floor(((float)rand()/RAND_MAX) * arrSize * p);
    // Check for pivot in local array
    // If you have the pivot, broadcast it, otherwise receive it
    int source = floor(index/arrSize);
    if(rank == source){
    	pivot = begin[index%arrSize];
    }
    MPI_Bcast(&pivot, 1, MPI_INT, source, comm);
    printf("rank %d is working\n", rank);

    // Split local array based on pivot
    //  - allocate second array
    //  - keep track of size of both arrays
    //  - realloc
    int num, le_size = 0, g_size = 0;
    int *greater = (int*) malloc(arrSize * sizeof(int));
    int *lesser = (int*) malloc(arrSize * sizeof(int));
    for (int i = 0; i < arrSize; i++) {
        num = begin[i];
        if (num <= pivot) {
            lesser[le_size++] = num;
        } else {
            greater[g_size++] = num;
        }
    }
    // begin = (int*) realloc(begin, sizeof(int) * le_size);
    // greater = (int*) realloc(greater, sizeof(int) * g_size);
    

    // Allgather to find total # of elements < and > pivot
    int* small = (int*) malloc(sizeof(int) * p);
    int* big = (int*) malloc(sizeof(int) * p);
    MPI_Allgather(&le_size, 1, MPI_INT, small, p, MPI_INT, comm);
    MPI_Allgather(&g_size, 1, MPI_INT, big, p, MPI_INT, comm);
    int smallsum = 0, bigsum = 0;
    for(int i = 0; i < p; i++){
    	smallsum += small[i];
    	bigsum += big[i];
    }

    // Decide # of processors for < and > pivot
    int l_proc_num = floor(worldsize * smallsum / (smallsum + bigsum));
    //int g_proc_num = commsize - l_proc_num;

    // Send < and > arrays to appropriate processors (using alltoall)
    // int** sendarr = (int**)calloc(p, sizeof(int));
    int* receivearr = (int*)calloc(arrSize, sizeof(int));
    // if(rank < l_proc_num){
    // 	sendarr[rank+l_proc_num] = greater;	
    // } else {
    // 	sendarr[(rank-l_proc_num)%l_proc_num] = begin;    
    // }

    int *space = (int*) malloc(sizeof(int) * p);
    for (int i = 0; i < p; i++) {
        if (i < l_proc_num) {
            space[i] = big[i];
        } else {
            space[i] = small[i];
        }
    }

    int b_send_size = big[rank];
    int s_send_size = small[rank];
    int myspace = space[rank];
    int *send_disp = (int*) calloc(p, sizeof(int));
    int *send_count = (int*) malloc(p * sizeof(int));
    int *rec_disp = (int*) calloc(p, sizeof(int));
    int *rec_count = (int*) malloc(p * sizeof(int));
    for (int i = 0; i < p; i++) {

        if (i > 0) {
            send_disp[i] = send_disp[i-1] + send_count[i-1];
            rec_disp[i] = rec_disp[i-1] + rec_count[i-1];
        }

        if (rank < l_proc_num) {
            if (i < l_proc_num) {
                send_count[i] = 0;
                rec_count[i] = 0;
            } else {
                // Sending counts
                if (space[i] >= b_send_size) {
                    send_count[i] = b_send_size;
                    b_send_size = 0;
                    space[i] -= b_send_size;
                } else {
                    send_count[i] = space[i];
                    b_send_size -= space[i];
                    space[i] = 0;
                }

                // Receiving counts
                if (myspace >= small[i]) {
                    rec_count[i] = small[i];
                    myspace -= small[i];
                } else {
                    rec_count[i] = myspace;
                    myspace = 0;
                }
            }
        } else {
            if (i >= l_proc_num) {
                send_count[i] = 0;
                rec_count[i] = 0;
            } else {
                // Sending counts
                if (space[i] >= s_send_size) {
                    send_count[i] = s_send_size;
                    s_send_size = 0;
                    space[i] -= s_send_size;
                } else {
                    send_count[i] = space[i];
                    s_send_size -= space[i];
                    space[i] = 0;
                }

                // Receiving counts
                if (myspace >= big[i]) {
                    rec_count[i] = big[i];
                    myspace -= big[i];
                } else {
                    rec_count[i] = myspace;
                    myspace = 0;
                }
            }
        }
    }

    if (rank < l_proc_num) {
        MPI_Alltoallv(greater, send_count, send_disp, MPI_INT,
            receivearr, rec_count, rec_disp, MPI_INT, comm);
    } else {
        MPI_Alltoallv(lesser, send_count, send_disp, MPI_INT,
            receivearr, rec_count, rec_disp, MPI_INT, comm);
    }

    // free(space);
    // free(send_count);
    // free(send_disp);
    // free(rec_count);
    // free(rec_disp);

    // if (rank >= l_proc_num) {
    //     // begin = (int*) realloc(begin, g_size * sizeof(int));
    //     for (int i = 0; i < g_size; i++) {
    //         begin[i] = greater[i];
    //     }
    // }

    for (int i = 0; i < arrSize; i++) {
        if (rank < l_proc_num) {
            if (i < le_size) {
                begin[i] = lesser[i];
            } else {
                begin[i] = receivearr[i - le_size];
            }
        } else {
            if (i < g_size) {
                begin[i] = greater[i];
            } else {
                begin[i] = receivearr[i - g_size];
            }
        }
    }

    // for(int i = 0; i < p; i++){
    // 	if(receivearr[i] != 0){
    // 		if (rank < l_proc_num) {
    // 			begin = (int *) realloc(begin, sizeof(int) * (le_size + small[i]));
    // 			for (int j = 0; j < small[i]; j++) {
    // 				begin[le_size + j] = receivearr[i][j];
    // 			}
    // 			le_size += small[i];
    // 		} else {
    // 			begin = (int *) realloc(begin, sizeof(int) * (g_size + big[i]));
    // 			for (int j = 0; j < big[i]; j++) {
    // 				begin[g_size + j] = receivearr[i][j];
    // 			}
    // 			g_size += big[i];
    // 		}
    // 	}
    // }

    // Dealloc greater
    // free(greater);
    // free(lesser);

    // Create two new communicators
    // MPI_Comm_split
    // Dealloc old comm
    MPI_Comm newcomm;
    int color;
    if(rank < l_proc_num){
    	color = 0;
    } else {
    	color = 1;
    }
    MPI_Comm_split(comm, color, rank, &newcomm);
    // MPI_Comm_free(&comm);

    // Call self recursively
    parallel_sort(begin, begin + (arrSize * sizeof(int)), newcomm);
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

