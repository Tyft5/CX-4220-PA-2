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

    // Each processor should have the same number for the current minute.
    // There's a small chance the seed will be different between processors
    // if the program runs in a very small window around when the minute
    // turns over. This is acceptably unlikely for now.
    srand(time(NULL) / 60);

    int p;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(comm, &rank);

    // Perform sort
    int *output;
    printf("1, rank %d\n", rank);
    MPI_Barrier(comm);
    int *temp = (int *) malloc((end - begin) * sizeof(int));
    for (int i = 0; i < end - begin; i++) {
        temp[i] = begin[i];
    }
    int arrSize = recursive_sort(temp, temp + (end - begin), &output, comm);

    // Communicate local array size
    int *sizes = (int *) malloc(p * sizeof(int));
    MPI_Allgather(&arrSize, 1, MPI_INT, sizes, 1, MPI_INT, comm);

    // Communicate local arrays
    int** arrays = (int **) malloc(p * sizeof(int*));
    MPI_Allgather(&output, 1, MPI_INT, arrays, 1, MPI_INT, comm);

    int whole_arr[(end - begin + 1) * p];
    int index = 0;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < sizes[i]; j++) {
            whole_arr[index++] = arrays[i][j];
        }
    }

    int total_size = index + 1;
    int rem = total_size % p;
    int start;
    if (rank < rem) {
        start = rank * (end - begin);
    } else {
        start = rem * (1 + end - begin) + (rank - rem) * (end - begin);
    }

    for (int i = 0; i < end - begin; i++) {
        begin[i] = whole_arr[start + i];
    }

    // Rearrange numbers so that all processors have equal size arrays
    // int sum = 0, corr_size = end - begin;
    // for (int i = 0; i <= rank; i++) {
    //     sum += sizes[i];
    // }
    // int to_send = (rank + 1) * (corr_size) - sum;

    // int *send_count = (int *) calloc(p, sizeof(int));
    // int *send_disp = (int *) calloc(p, sizeof(int));
    // for (int i = 0; i < p; i++) {
    //     if (i == rank - 1) {
    //         for (int j = 0; j < rank; j++) {
    //             if (sizes[j] < corr_size) {
    //                 to_send += corr_size - sizes[j];
    //                 sizes[rank] -= to_send;
    //                 send_count[j] = to_send;
    //             }
    //         }
    //         if (sizes[rank] > corr_size) {

    //         }
    //     } else {
    //         if (sizes[i] < corr_size) {
    //             sizes[i+1] -= corr_size - sizes[i];
    //         } else {
    //             sizes[i+1] += sizes[i] - corr_size;
    //         }
    //     }
    // }


    /////////////////////////////////////////////////////////////////

 //    int p, rank, pivot;
 //    MPI_Comm_size(MPI_COMM_WORLD, &p);
 //    MPI_Comm_rank(comm, &rank);

 //    // Terminating condition
	// int commsize;
	// int arrSize = (end - begin);
	// MPI_Comm_size(comm, &commsize);
	// if(commsize == 1){
	// 	qsort(begin, arrSize, sizeof(int), cmpfunc);
 //        return;
	// }
    // // Call seeding helper function
    // seed_rand(commsize, p);

 //    // Generate a pivot
 //    int index = (int)floor(((float)rand()/RAND_MAX) * arrSize * p);
 //    // Check for pivot in local array
 //    // If you have the pivot, broadcast it, otherwise receive it
 //    int source = floor(index/arrSize);
 //    if(rank == source){
 //    	pivot = begin[index%arrSize];
 //    }
 //    MPI_Bcast(&pivot, 1, MPI_INT, source, comm);

 //    // Split local array based on pivot
 //    //  - allocate second array
 //    //  - keep track of size of both arrays
 //    //  - realloc
 //    int num, le_size = 0, g_size = 0;
 //    int *greater = (int*) malloc(arrSize * sizeof(int));
 //    int *lesser = (int*) malloc(arrSize * sizeof(int));
 //    for (int i = 0; i < arrSize; i++) {
 //        num = begin[i];
 //        if (num <= pivot) {
 //            lesser[le_size++] = num;
 //        } else {
 //            greater[g_size++] = num;
 //        }
 //    }
 //    // begin = (int*) realloc(begin, sizeof(int) * le_size);
 //    // greater = (int*) realloc(greater, sizeof(int) * g_size);
    

 //    // Allgather to find total # of elements < and > pivot
 //    int* small = (int*) calloc(p,sizeof(int));
 //    int* big = (int*) calloc(p,sizeof(int));
 //    MPI_Allgather(&le_size, 1, MPI_INT, small, 1, MPI_INT, comm);
 //    MPI_Barrier(comm);
 //    MPI_Allgather(&g_size, 1, MPI_INT, big, 1, MPI_INT, comm);
 //    int smallsum = 0, bigsum = 0;
 //    printf("Rank %d sizes: %d %d\n", rank, le_size, g_size);
 //    for(int i = 0; i < p; i++){
 //    	smallsum += small[i];
 //    	bigsum += big[i];
 //    }

 //    // Decide # of processors for < and > pivot
 //    int l_proc_num = ceil(commsize * smallsum / (smallsum + bigsum));
 //    //int g_proc_num = commsize - l_proc_num;

 //    // Send < and > arrays to appropriate processors (using alltoall)
 //    // int** sendarr = (int**)calloc(p, sizeof(int));
 //    int* receivearr = (int*)calloc(arrSize, sizeof(int));
 //    // if(rank < l_proc_num){
 //    // 	sendarr[rank+l_proc_num] = greater;	
 //    // } else {
 //    // 	sendarr[(rank-l_proc_num)%l_proc_num] = begin;    
 //    // }

    // int *space = (int*) malloc(sizeof(int) * p);
    // for (int i = 0; i < p; i++) {
    //     if (i < l_proc_num) {
    //         space[i] = big[i];
    //     } else {
    //         space[i] = small[i];
    //     }
    // }

    // int b_send_size = big[rank];
    // int s_send_size = small[rank];
    // int myspace = space[rank];
    // int *send_disp = (int*) calloc(p, sizeof(int));
    // int *send_count = (int*) malloc(p * sizeof(int));
    // int *rec_disp = (int*) calloc(p, sizeof(int));
    // int *rec_count = (int*) malloc(p * sizeof(int));
    // for (int i = 0; i < p; i++) {

    //     if (i > 0) {
    //         send_disp[i] = send_disp[i-1] + send_count[i-1];
    //         rec_disp[i] = rec_disp[i-1] + rec_count[i-1];
    //     }

    //     if (rank < l_proc_num) {
    //         if (i < l_proc_num) {
    //             send_count[i] = 0;
    //             rec_count[i] = 0;
    //         } else {
    //             // Sending counts
    //             if (space[i] >= b_send_size) {
    //                 send_count[i] = b_send_size;
    //                 b_send_size = 0;
    //                 space[i] -= b_send_size;
    //             } else {
    //                 send_count[i] = space[i];
    //                 b_send_size -= space[i];
    //                 space[i] = 0;
    //             }

    //             // Receiving counts
    //             if (myspace >= small[i]) {
    //                 rec_count[i] = small[i];
    //                 myspace -= small[i];
    //             } else {
    //                 rec_count[i] = myspace;
    //                 myspace = 0;
    //             }
    //         }
    //     } else {
    //         if (i >= l_proc_num) {
    //             send_count[i] = 0;
    //             rec_count[i] = 0;
    //         } else {
    //             // Sending counts
    //             if (space[i] >= s_send_size) {
    //                 send_count[i] = s_send_size;
    //                 s_send_size = 0;
    //                 space[i] -= s_send_size;
    //             } else {
    //                 send_count[i] = space[i];
    //                 s_send_size -= space[i];
    //                 space[i] = 0;
    //             }

    //             // Receiving counts
    //             if (myspace >= big[i]) {
    //                 rec_count[i] = big[i];
    //                 myspace -= big[i];
    //             } else {
    //                 rec_count[i] = myspace;
    //                 myspace = 0;
    //             }
    //         }
    //     }
    // }
    
    // if (rank < l_proc_num) {
    //     MPI_Barrier(comm);
    //     MPI_Alltoallv(greater, send_count, send_disp, MPI_INT,
    //         receivearr, rec_count, rec_disp, MPI_INT, comm);
    // } else {
    //     MPI_Barrier(comm);
    //     MPI_Alltoallv(lesser, send_count, send_disp, MPI_INT,
    //         receivearr, rec_count, rec_disp, MPI_INT, comm);
    // }

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

    // for (int i = 0; i < arrSize; i++) {
    //     if (rank < l_proc_num) {
    //         if (i < le_size) {
    //             begin[i] = lesser[i];
    //         } else {
    //             begin[i] = receivearr[i - le_size];
    //         }
    //     } else {
    //         if (i < g_size) {
    //             begin[i] = greater[i];
    //         } else {
    //             begin[i] = receivearr[i - g_size];
    //         }
    //     }
    // }

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

    // // Create two new communicators
    // // MPI_Comm_split
    // // Dealloc old comm
    // MPI_Comm newcomm;
    // int color;
    // if(rank < l_proc_num){
    // 	color = 0;
    // } else {
    // 	color = 1;
    // }
    // MPI_Comm_split(comm, color, rank, &newcomm);
    // // MPI_Comm_free(&comm);

    // Call self recursively
    // parallel_sort(begin, end, newcomm);
}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

// Function to seed RNG once with the same seed on each processor
// void seed_rand(int commsize, int worldsize) {
//     // Only seed once
//     if (commsize == worldsize) {
//         // Each processor should have the same number for the current minute.
//         // There's a small chance the seed will be different between processors
//         // if the program runs in a very small window around when the minute
//         // turns over. This is acceptably unlikely for now.
//         srand(time(NULL) / 60);
//     }
// }

int recursive_sort(int *begin, int *end, int** out, MPI_Comm comm) {
    int p, rank, pivot, commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    unsigned int the_time;
    int arrSize = (end - begin);

    // Terminating condition
    if(commsize == 1){
        qsort(begin, arrSize, sizeof(int), cmpfunc);
        out = &begin;
        return end - begin;
    }

    // Generate a pivot
    int index = (int)floor(((float)rand()/RAND_MAX) * arrSize * p);

    // Check for pivot in local array
    // If you have the pivot, broadcast it, otherwise receive it
    int source = floor(index/arrSize);
    if(rank == source){
        pivot = begin[index%arrSize];
    }
    MPI_Bcast(&pivot, 1, MPI_INT, source, comm);

    // Split local array based on pivot
    //  - allocate second array
    //  - keep track of size of both arrays
    //  - realloc
    int num, le_size = 0, g_size = 0;
    printf("2, rank %d\n", rank);
    MPI_Barrier(comm);
    int *greater = (int*) malloc(arrSize * sizeof(int));
    // int *lesser = (int*) malloc(arrSize * sizeof(int));
    for (int i = 0; i < arrSize; i++) {
        num = begin[i];
        if (num <= pivot) {
            begin[le_size++] = num;
        } else {
            greater[g_size++] = num;
        }
    }
    // begin = (int*) realloc(begin, sizeof(int) * le_size);
    // greater = (int*) realloc(greater, sizeof(int) * g_size);
    

    // Allgather to find total # of elements < and > pivot
    printf("3, rank %d\n", rank);
    MPI_Barrier(comm);
    int* small = (int*) malloc(p * sizeof(int));
    printf("4, rank %d\n", rank);
    MPI_Barrier(comm);
    int* big = (int*) malloc(p * sizeof(int));
    MPI_Allgather(&le_size, 1, MPI_INT, small, 1, MPI_INT, comm);
    MPI_Barrier(comm);
    MPI_Allgather(&g_size, 1, MPI_INT, big, 1, MPI_INT, comm);
    int smallsum = 0, bigsum = 0;
    printf("Rank %d sizes: %d %d\n", rank, le_size, g_size);
    MPI_Barrier(comm);
    for(int i = 0; i < p; i++){
        smallsum += small[i];
        bigsum += big[i];
    }

    // Decide # of processors for < and > pivot
    printf("Rank %d sums: %d %d\n", rank, smallsum, bigsum);
    MPI_Barrier(comm);
    int l_proc_num = ceil(commsize * smallsum / (smallsum + bigsum));
    int g_proc_num = commsize - l_proc_num;

    // Send < and > arrays to appropriate processors (using alltoall)
    printf("5, rank %d\n", rank);
    MPI_Barrier(comm);
    // int** sendarr = (int**) malloc(p, sizeof(int*));

    printf("6, rank %d\n", rank);
    the_time = time(0);
    while (time(0) < (the_time + 2));
    MPI_Barrier(comm);
    
    // int** receivearr = (int**) calloc(p, sizeof(int*));
    // if(rank < l_proc_num){
    //     sendarr[l_proc_num + (rank % g_proc_num)] = greater; 
    // } else if (l_proc_num != 0) {
    //     sendarr[(rank-l_proc_num) % l_proc_num] = begin;    
    // } else {}

    // MPI_Alltoall(sendarr, p, MPI_INT, receivearr, p, MPI_INT, comm);

    int *send_counts = (int *) calloc(p, sizeof(int));
    int *send_disp = (int *) calloc(p, sizeof(int));
    int *rec_count = (int *) malloc(p * sizeof(int));
    int *rec_disp = (int *) calloc(p, sizeof(int));
    int *rec_arr = (int *) calloc(smallsum + bigsum, sizeof(int));
    for (int i = 0; i < p; i++) {
        if (rank < l_proc_num) {
            rec_count[i] = small[i];
            if (i > 0) {
                rec_disp[i] = rec_disp[i - 1] + small[i];
            }
            // else {
            //     rec_disp[i] = small[i];
            // }
        } else {
            rec_count[i] = big[i];
            if (i > 0) {
                rec_disp[i] = rec_disp[i - 1] + big[i];
            }
            // else {
            //     rec_disp[i] = big[i];
            // }
        }
    }

    if (rank < l_proc_num) {
        send_counts[l_proc_num + (rank % g_proc_num)] = g_size;
        MPI_Alltoallv(greater, send_counts, send_disp, MPI_INT,
            rec_arr, rec_count, rec_disp, MPI_INT, comm);
    } else {
        send_counts[(rank-l_proc_num) % l_proc_num] = le_size;
        MPI_Alltoallv(begin, send_counts, send_disp, MPI_INT,
            rec_arr, rec_count, rec_disp, MPI_INT, comm);
    }

    printf("Before realloc, rank %d\n", rank);
    MPI_Barrier(comm);
    // Move greater into begin so all processors have numbers that weren't
    // send in begin.
    // begin = (int*) realloc(begin, (smallsum + bigsum) * sizeof(int));
    int *newbegin = (int *) malloc((smallsum + bigsum) * sizeof(int));
    if (rank >= l_proc_num) {
        // begin = (int*) realloc(begin, g_size * sizeof(int));
        for (int i = 0; i < g_size; i++) {
            newbegin[i] = greater[i];
        }
    } else {
        for (int i = 0; i < le_size; i++) {
            newbegin[i] = begin[i];
        }
    }

    // Add the received numbers onto the existing numbers
    printf("Before concat, rank %d\n", rank);
    MPI_Barrier(comm);
    the_time = time(0);
    while (time(0) < (the_time + 2));
    // int *dummy;
    int sum = 0, c = 1, d = 0;
    if (rank < l_proc_num) {
        for (int i = 0; i < l_proc_num + rank; i++) { sum += small[i]; }
    } else {
        for (int i = 0; i < rank - l_proc_num; i++) { sum += big[i]; }
    }
    for(int i = 0; i < p; i++){
        // if(receivearr[i] != 0){
            if (rank < l_proc_num) {
                // dummy = (int *) realloc(begin, sizeof(int) * (le_size + small[i]));
                // if (dummy) {
                //     begin = dummy;
                // } else {
                //     printf("Failed to realloc\n");
                //     exit(1);
                // }
                // begin = (int *) realloc(begin, sizeof(int) * (le_size + small[i]));

                // for (int j = 0; j < small[i]; j++) {
                //     dummy = receivearr[i];
                //     // newbegin[le_size + j] = receivearr[i][j];
                //     newbegin[le_size + j] = dummy[j];
                // }
                // le_size += small[i];

                while ((c * l_proc_num) + rank < p) {
                    for (int k = 0; k < small[c * l_proc_num + rank]; k++) {
                        newbegin[le_size++] = rec_arr[sum + k];
                    }
                    for (int l = c * l_proc_num + rank; l < (c+1) * l_proc_num + rank; l++) {
                        sum += small[l];
                    }
                    c++;
                }

            } else {
                // dummy = (int *) realloc(begin, sizeof(int) * (g_size + big[i]));
                // if (dummy) {
                //     begin = dummy;
                // } else {
                //     printf("Failed to realloc\n");
                //     exit(1);
                // }
                // begin = (int *) realloc(begin, sizeof(int) * (g_size + big[i]));

                // for (int j = 0; j < big[i]; j++) {
                //     // newbegin[g_size + j] = receivearr[i][j];
                //     dummy = receivearr[i];
                //     newbegin[g_size + j] = dummy[j];
                // }
                // g_size += big[i];

                while (rank - l_proc_num + (d * g_proc_num) < l_proc_num) {
                    for (int k = 0; k < big[rank - l_proc_num + (d * g_proc_num)]; k++) {
                        newbegin[g_size++] = rec_arr[sum + k];
                    }
                    for (int l = rank - l_proc_num + (d * g_proc_num);
                            l < rank - l_proc_num + ((d+1) * g_proc_num); l++) {
                        sum += big[l];
                    }
                    d++;
                }
            }
        //}
    }
    printf("After concat, rank %d\n", rank);

    // Create new communicator
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

    if (rank < l_proc_num) {
        newbegin = (int *) realloc(newbegin, le_size * sizeof(int));
    } else {
        newbegin = (int *) realloc(newbegin, g_size * sizeof(int));
    }

    // Call self recursively and pass final output back down the stack
    int *output;
    int fin_count;
    if (rank < l_proc_num) {
        fin_count = recursive_sort(newbegin, newbegin + (le_size * sizeof(int)),
            &output, newcomm);
    } else {
        fin_count = recursive_sort(newbegin, newbegin + (g_size * sizeof(int)),
            &output, newcomm);
    }
    out = &output;
    return fin_count;
}

