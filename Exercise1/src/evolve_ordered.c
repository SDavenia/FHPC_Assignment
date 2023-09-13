#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "evolve_ordered.h"
#include "read_write_parallel.h"


void evolve_ordered(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s){
  if(s == 0)
    s = n_steps; // To ensure if we call s = 0, only last one is printed. 
  if(size == 1)
    evolve_ordered_OMP(current, k, n_steps, s);
  else
    evolve_ordered_MPI(current, k, n_steps, rank, size, rows_read, s);
}

void evolve_ordered_OMP(unsigned char* current, int k, int n_steps, int s){
  for(int n_step=0; n_step < n_steps; n_step++){
    /*
    FILE* prova_file;
    char nome_file[] = "prova_file.txt";
    prova_file=fopen(nome_file, "w");
    */
    int nthreads;
    #pragma omp parallel
    {
      int myid = omp_get_thread_num();
      #pragma omp master
        nthreads = omp_get_num_threads();
        
      #pragma omp for ordered
      for(int i=1;i<k;i++){
        for(int j=0; j<k;j++){
          #pragma omp ordered
          {
            // fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid,i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
            current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
      }

      // Update bottom row of frame: since we are only copying no need for ordered
      #pragma omp for
      for(int j=0;j<k;j++){
        current[(k+1)*k+j] = current[k+j];
      }

      // Update last inner row of current and upper row of frame
      #pragma omp for ordered
      for(int j = 0; j< k; j++){
        #pragma omp ordered
        {
          int i = k;
          int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                current[(i-1)*k + j] + current[(i+1)*k + j];
          current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
          // Frame
          current[j] = current[i*k+j];
        }
      }
    }
    
    if((n_step+1) % s == 0){
      char file_path[45] = "images/evolve_ordered/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, 0, 1, k+2);
    }
  }
}

void evolve_ordered_MPI(unsigned char* current, int k, int n_steps, int rank, int size, int rows_read, int s){
  MPI_Request request[2];
  for(int n_step=0; n_step < n_steps; n_step++){

    // All processes wait until they have the current status of the matrix to update.

    if(rank == size -1){   // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, 0+n_step+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(rank != 0){ // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, rank-1, rank-1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    

    // Evolve parallel using openMP
    int nthreads;
    #pragma omp parallel
    {
      int myid = omp_get_thread_num();
      #pragma omp master
        nthreads = omp_get_num_threads();
        
      #pragma omp for ordered
      for(int i=1;i<rows_read+1;i++){
        for(int j=0; j<k;j++){
          #pragma omp ordered
          {
            // fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid,i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
            current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
          }
        }
      }
    } // Implicit barrier at the end of parallel region

    if(rank == 0){ // TAG IS ALWAYS THE RANK OF THE SENDING PROCESS + n_step
      // Send message to the next process
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      // Send message to last process
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, size-1, rank + n_step + 1, MPI_COMM_WORLD, &request[0]);
      
      // Blocking receive message
      // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1 + n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Upper row receive
      MPI_Recv(current, k, MPI_UNSIGNED_CHAR, size-1, size-1 + n_step + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }else if(rank == size-1){
      // Send message to process 0
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, 0, rank + n_step +1 , MPI_COMM_WORLD, &request[1]);
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
    }else{
      // Send message to process after
      MPI_Isend(current + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank + n_step, MPI_COMM_WORLD, &request[1]);
      // Send message to process before
      MPI_Isend(current+k, k, MPI_UNSIGNED_CHAR, rank-1, rank + n_step, MPI_COMM_WORLD, &request[0]);
      
      // BLOCKING RECEIVE
      // Lower row receive
      MPI_Recv(current + k + rows_read*k, k, MPI_UNSIGNED_CHAR, rank+1, rank+1+n_step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    /*
    if(n_step == ){
      FILE* prova_file;
      char nome_file[] = "prova_read3.txt";
      if(n_step==0){
        if(rank==0){
          prova_file = fopen(nome_file, "w");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          
          fprintf(prova_file, "\n");
          fclose(prova_file);
        }
      }else{
        if(rank==0){
          prova_file = fopen(nome_file, "w");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);
          
          fprintf(prova_file, "\n");
          fclose(prova_file);
        } 
      }

      MPI_Barrier(MPI_COMM_WORLD);
      
      if(rank==1){
          prova_file = fopen(nome_file, "a");
          fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
          for(int i = 0; i < (rows_read+2) * k; i++)
              fprintf(prova_file, "%u ",current[i]);

        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==2){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==3){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if(rank==4){
        prova_file = fopen(nome_file, "a");
        fprintf(prova_file,"I am process %d of step %d\n", rank, n_step);
        for(int i = 0; i < (rows_read+2) * k; i++)
            fprintf(prova_file, "%u ",current[i]);
        fprintf(prova_file, "\n");
        fclose(prova_file);
      }
    }
    */
    
    // MPI_Waitall(2, request, MPI_STATUS_IGNORE);
    // MPI_Barrier(MPI_COMM_WORLD);
    if((n_step+1) % s == 0){
      char file_path[45] = "images/evolve_ordered/"; // Sufficiently large
      char filename[20];
      
      snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
      strcat(file_path, filename);
      write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
    }
  }
}
