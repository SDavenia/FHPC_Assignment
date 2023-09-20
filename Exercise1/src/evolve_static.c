#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "evolve_static.h"
#include "read_write_parallel.h"


void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s){
  if(s == 0)
    s = n_steps; // To ensure if we call s = 0, only last one is printed. 
  if(size == 1)
    evolve_static_OMP(current, next, k, n_steps, s);
  else
    evolve_static_MPI(current, next,k , n_steps, rank, size, rows_read, s);
}

void evolve_static_OMP(unsigned char* current, unsigned char* next, int k, int n_steps, int s){
    for(int n_step=0; n_step < n_steps; n_step++){
      int nthreads;
      #pragma omp parallel
      {
        int myid = omp_get_thread_num();
        #pragma omp master
          nthreads = omp_get_num_threads();
      
        #pragma omp for
        for(int i=1;i<k+1;i++){
          for(int j=0; j<k;j++){
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                  current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                  current[(i-1)*k + j] + current[(i+1)*k + j];
              next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
        // Here there is an implicit barrier of the for loop

        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[i] = next[(k)*k + i];
        }

        #pragma omp single nowait
        for (int i = 0; i < k; i++){
            next[(k+1)*k+i]=next[k+i];
        }
      }
      // Here there is an implicit barrier for the end of the parallel region.

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

      if((n_step+1) % s == 0){
        char file_path[45] = "images/evolve_static/";
        char filename[20];
        
        snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
        strcat(file_path, filename);
        write_pgm_parallel(current+k, 255, k, k, file_path, 0, 1, k+2);
      }
    
    }

}

void evolve_static_MPI_blocking(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s){
  MPI_Request request[2]; // Array of MPI_Request, necessary to ensure that the non blocking receive is complete

  for(int n_step=0; n_step < n_steps; n_step++){
      int nthreads;
      #pragma omp parallel
      {
        int myid = omp_get_thread_num();
        #pragma omp master
          nthreads = omp_get_num_threads();
      
        #pragma omp for
        for(int i=1;i<rows_read+1;i++){
          for(int j=0; j<k;j++){
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                  current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                  current[(i-1)*k + j] + current[(i+1)*k + j];
              next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
        
      }
      // Here there is an implcit barrier for the end of the parallel region.

      if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
        // Send message to last process
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
        // Send message to the next process
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // Blocking receive message
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Upper row receive from final process
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }else if(rank == size-1){
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process 0
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
        
        // BLOCKING RECEIVE
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }else{
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process after
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // BLOCKING RECEIVE
        // Upper row receive
        MPI_Recv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Lower row receive
        MPI_Recv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      MPI_Waitall(2, request, MPI_STATUS_IGNORE); // To ensure they all do not alter the buffer by moving to next iteration.

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

    if((n_step+1) % s == 0){
      char file_path[45] = "images/evolve_static/";
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
    }
  }

}

void evolve_static_MPI(unsigned char* current, unsigned char* next, int k, int n_steps, int rank, int size, int rows_read, int s){
  MPI_Request request[4]; // Array of MPI_Request, necessary to ensure that the non blocking receive is complete

  for(int n_step=0; n_step < n_steps; n_step++){
      int nthreads;
      #pragma omp parallel
      {
        int myid = omp_get_thread_num();
        #pragma omp master
          nthreads = omp_get_num_threads();
      
        #pragma omp for
        for(int i=1;i<rows_read+1;i++){
          for(int j=0; j<k;j++){
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                  current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                  current[(i-1)*k + j] + current[(i+1)*k + j];
              next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
      }
      // Here there is an implcit barrier for the end of the parallel region.

      if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
        // Send message to last process
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
        // Send message to the next process
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // Blocking receive message
        // Lower row receive
        MPI_Irecv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, &request[2]);
        // Upper row receive from final process
        MPI_Irecv(next, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, &request[3]);

      }else if(rank == size-1){
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process 0
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
        
        // Lower row receive
        MPI_Irecv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, &request[2]);
        // Upper row receive
        MPI_Irecv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, &request[3]);
      }else{
        // Send message to process before
        MPI_Isend(next+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
        // Send message to process after
        MPI_Isend(next + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
        
        // Upper row receive
        MPI_Irecv(next, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, &request[2]);
        // Lower row receive
        MPI_Irecv(next + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+n_step, MPI_COMM_WORLD, &request[3]);
      }
      MPI_Waitall(4, request, MPI_STATUS_IGNORE); // To ensure they all do not alter the buffer by moving to next iteration.

      unsigned char* tmp;
      tmp = next;
      next = current;
      current=tmp;

    if((n_step+1) % s == 0){
      char file_path[45] = "images/evolve_static/";
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
    }
  }

}

