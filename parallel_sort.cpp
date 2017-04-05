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
    unsigned int the_time;

    // Perform sort
    int *output;
    // printf("1, rank %d\n", rank);
    MPI_Barrier(comm);
    int *temp = (int *) malloc((end - begin) * sizeof(int));
    for (int i = 0; i < end - begin; i++) {
        temp[i] = begin[i];
    }

    MPI_Barrier(comm);
    printf("0, rank %d\n", rank);

    int arrSize = recursive_sort(temp, temp + (end - begin), &output, (end-begin-1)*p, comm);
    MPI_Barrier(comm);

    MPI_Barrier(comm);
    printf("Rec done, rank %d\n", rank);
    the_time = time(0);
    while (time(0) < the_time + 0.5) ;

    // Communicate local array size
    int *sizes = (int *) malloc(p * sizeof(int));
    MPI_Allgather(&arrSize, 1, MPI_INT, sizes, 1, MPI_INT, comm);

    MPI_Barrier(comm);
    printf("6, rank %d\n", rank);

    // Communicate local arrays
    int *send_disp = (int *) calloc(p, sizeof(int));
    int *send_count = (int *) malloc(p * sizeof(int));
    int *rec_disp = (int *) calloc(p, sizeof(int));
    int *rec_arr = (int *) malloc((end - begin + 1) * p * sizeof(int));
    int total_size = 0;
    for (int i = 0; i < p; i++) {
        send_count[i] = arrSize;
        if (i > 0) {
            rec_disp[i] = rec_disp[i - 1] + sizes[i - 1];
        }
        total_size += sizes[i];
    }

    MPI_Barrier(comm);
    printf("7, rank %d\n", rank);
    MPI_Alltoallv(output, send_count, send_disp, MPI_INT,
        rec_arr, sizes, rec_disp, MPI_INT, comm);

    int rem = total_size % p;
    int start;
    if (rank < rem) {
        start = rank * (end - begin);
    } else {
        start = rem * (1 + end - begin) + (rank - rem) * (end - begin);
    }

    MPI_Barrier(comm);
    printf("8, rank %d\n", rank);

    for (int i = 0; i < end - begin; i++) {
        begin[i] = rec_arr[start + i];
    }


}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int recursive_sort(int *begin, int *end, int** out, int comm_arr_size, MPI_Comm comm) {
    int p, rank, pivot, commsize, g_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

    // unsigned int the_time;
    int arrSize = (end - begin);

    unsigned int the_time;

    // Terminating condition
    if(commsize == 1){
    //     MPI_Barrier(comm);
    //     printf("Term, rank %d\n", g_rank);
    //     the_time = time(0);
    //     while (time(0) < the_time + 0.5) ;

    //     printf("Rank %d, arrSize %d\n", g_rank, arrSize);

        if (arrSize > 1) {
            qsort(begin, arrSize, sizeof(int), cmpfunc);
        }

        // MPI_Barrier(comm);
        // printf("Sorted, rank %d\n", g_rank);
        // the_time = time(0);
        // while (time(0) < the_time + 0.5) ;
        *out = begin;
        // printf("Done, rank %d %d\n", rank, g_rank);
        return end - begin;
    }

    // MPI_Barrier(comm);
    // printf("1, rank %d\n", g_rank);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;

    int source, l_proc_num, g_proc_num = 0;
    int num, le_size = 0, g_size = 0;
    int *small = (int*) malloc(commsize * sizeof(int));
    int *big = (int*) malloc(commsize * sizeof(int));
    int *greater = (int*) malloc(arrSize * sizeof(int));
    int smallsum = 0, bigsum = 0;
    double temp_rand;
    while (g_proc_num == 0) {
        le_size = 0;
        g_size = 0;
        // Generate a pivot
        rand();
        // Check for pivot in local array
        // If you have the pivot, broadcast it, otherwise receive it
        //source = floor(index/arrSize);
        temp_rand = rand();
        source = (int)(temp_rand / ((double)RAND_MAX +1)*commsize);
        if(rank == source){
            pivot = begin[(int)((double)temp_rand / ((double)RAND_MAX +1)* arrSize)];
            //pivot = begin[index%arrSize];
        }
        MPI_Barrier(comm);
        MPI_Bcast(&pivot, 1, MPI_INT, source, comm);

        // Split local array based on pivot
        //  - allocate second array
        //  - keep track of size of both arrays
        //  - realloc
        for (int i = 0; i < arrSize; i++) {
            num = begin[i];
            if (num <= pivot) {
                begin[le_size++] = num;
            } else {
                greater[g_size++] = num;
            }
        }

        MPI_Allgather(&le_size, 1, MPI_INT, small, 1, MPI_INT, comm);
        MPI_Barrier(comm);
        MPI_Allgather(&g_size, 1, MPI_INT, big, 1, MPI_INT, comm);

        MPI_Barrier(comm);
        smallsum = 0;
        for(int i = 0; i < commsize; i++){
            smallsum += small[i];
            bigsum += big[i];
        }

        // Decide # of processors for < and > pivot
        MPI_Barrier(comm);
        l_proc_num = (int) ceil(float(commsize * smallsum) / (smallsum + bigsum));
        // if (l_proc_num == commsize) l_proc_num--;
        g_proc_num = commsize - l_proc_num;

        if (commsize == 2) {
            l_proc_num = 1;
            g_proc_num = 1;
        }

        if (l_proc_num == 0) {
            g_proc_num = 0;
        }
    }

    // MPI_Barrier(comm);
    // printf("2, rank %d\n", g_rank);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;
    // MPI_Alltoall(sendarr, p, MPI_INT, receivearr, p, MPI_INT, comm);

    int *send_counts = (int *) calloc(commsize, sizeof(int));
    int *send_disp = (int *) calloc(commsize, sizeof(int));
    int *rec_count = (int *) calloc(commsize, sizeof(int));
    int *rec_disp = (int *) calloc(commsize, sizeof(int));
    int *rec_arr = (int *) calloc(smallsum + bigsum, sizeof(int));
    int ind, rec_sum = 0;
    for (int i = 0; i < commsize; i++) {
        if (rank < l_proc_num) {
            ind = ((i+1) * l_proc_num) + rank;
            if (ind < commsize) {
                rec_count[ind] = small[ind];

            }

            if (i > 0) {
                rec_disp[i] = rec_disp[i - 1] + rec_count[i - 1];
            }
 
        } else {
            ind = rank - l_proc_num + (i * g_proc_num);
            if (ind < l_proc_num) {
                rec_count[ind] = big[ind];

            }

            if (i > 0) {
                rec_disp[i] = rec_disp[i - 1] + rec_count[i - 1];
            }
        }
        rec_sum += rec_count[i];
    }

    // MPI_Barrier(comm);
    // printf("3, rank %d\n", g_rank);
    // printf("Rank %d, l/g_proc_num %d %d\n", g_rank, l_proc_num, g_proc_num);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;

    if (rank < l_proc_num) {
        send_counts[l_proc_num + (rank % g_proc_num)] = g_size;
        for (int i = l_proc_num + (rank % g_proc_num) + 1; i < commsize; i++) {
            send_disp[i] = g_size;
        }

        MPI_Barrier(comm);
        MPI_Alltoallv(greater, send_counts, send_disp, MPI_INT,
            rec_arr, rec_count, rec_disp, MPI_INT, comm);
    } else {
        send_counts[(rank-l_proc_num) % l_proc_num] = le_size;
        for (int i = (rank-l_proc_num) % l_proc_num + 1; i < commsize; i++) {
            send_disp[i] = le_size;
        }

        MPI_Barrier(comm);
        MPI_Alltoallv(begin, send_counts, send_disp, MPI_INT,
            rec_arr, rec_count, rec_disp, MPI_INT, comm);
    }

    // MPI_Barrier(comm);
    // printf("4, rank %d\n", g_rank);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;

    MPI_Barrier(comm);
    // Move greater into begin so all processors have numbers that weren't
    // send in begin.
    int *newbegin = (int *) malloc((smallsum + bigsum) * sizeof(int));
    if (rank < l_proc_num) {
        for (int i = 0; i < le_size; i++) {
            newbegin[i] = begin[i];
        }
        for (int i = 0; i < rec_sum; i++) {
            newbegin[le_size + i] = rec_arr[i];
        }
        le_size += rec_sum;
    } else {
        for (int i = 0; i < g_size; i++) {
            newbegin[i] = greater[i];
        }
        for (int i = 0; i < rec_sum; i++) {
            newbegin[g_size + i] = rec_arr[i];
        }
        g_size += rec_sum;
    }

    // Create new communicator
    // Dealloc old comm
    MPI_Comm newcomm;
    int color;
    if(rank < l_proc_num){
        color = 0;
    } else {
        color = 1;
    }
    MPI_Barrier(comm);
    MPI_Comm_split(comm, color, rank, &newcomm);
    //MPI_Comm_free(&comm);

    // MPI_Barrier(comm);
    // printf("5, rank %d\n", g_rank);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;

    if (rank < l_proc_num) {
        if (le_size > 0) {
            newbegin = (int *) realloc(newbegin, le_size * sizeof(int));
        }
    } else {
        if (g_size > 0) {
            newbegin = (int *) realloc(newbegin, g_size * sizeof(int));
        }
    }

    // MPI_Barrier(comm);
    // printf("10, rank %d\n", g_rank);
    // the_time = time(0);
    // while (time(0) < the_time + 0.5) ;

    // Call self recursively and pass final output back down the stack
    int *output;
    int fin_count;
    if (rank < l_proc_num) {
        fin_count = recursive_sort(newbegin, newbegin + le_size,
            &output, smallsum, newcomm);
    } else {
        fin_count = recursive_sort(newbegin, newbegin + g_size,
            &output, bigsum, newcomm);
    }
    *out = output;
    return fin_count;
}

