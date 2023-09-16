#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "black_white.h"


void evolve_black_white_even(unsigned char *current, int k, int n_steps, int s){
  int nthreads;
  for(int n_step=0; n_step < n_steps; n_step++){
      #pragma omp parallel
      {
          int myid = omp_get_thread_num();
          // Update white inner positions (consider that the first cell is white)
          for(int i=1;i<k;i++){
              int t = (i+1)%2;
              #pragma omp for
              for(int j = t; j<k;j+=2){
                  int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
                  current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              }
          }
          // Update bottom row of frame (only white positions)
          #pragma omp for
          for(int j=0;j<k;j+=2){
              current[(k+1)*k+j] = current[k+j];
          }

          // Update last inner row of current and upper row of frame
          // Update white positions (consider that the first cell is white)
          #pragma omp for
          for(int j = 0; j< k; j+=2){
              int i = k;
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
              current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              // Frame
              current[j] = current[i*k+j];
          }
          
          // Update black positions (consider that the first black cell is the second one of the row)
          for(int i=1;i<k;i++){
              int t = i%2;
              #pragma omp for
              for(int j=t; j<k;j+=2){
                  int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
                  current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              }
          }
          // Update bottom row of frame (only black positions)
          #pragma omp for
          for(int j=1;j<k;j+=2){
              current[(k+1)*k+j] = current[k+j];
          }

          // Update last inner row of current and upper row of frame
          // Update black positions (consider that the first black cell is the second one of the row)
          #pragma omp for
          for(int j = 1; j< k; j+=2){
              int i = k;
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
              current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              // Frame
              current[j] = current[i*k+j];
          }
      }
        if((n_step+1) % s == 0){
            char file_path[45] = "images/evolve_black_white/"; // Sufficiently large
            char filename[20];
            
            snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
            strcat(file_path, filename);
            write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
        }
  }
}

void evolve_black_white_odd(unsigned char *current, int k, int n_steps, int s){
    int nthreads;
    /*FILE* prova_file;
    char nome_file[] = "prova_file.txt";
    prova_file=fopen(nome_file, "w");
    */
    for(int n_step=0; n_step < n_steps; n_step++){
        #pragma omp parallel
        {
            int myid = omp_get_thread_num();
            // Update white positions (consider that the first cell is white)
            for(int i=1;i<k;i++){
                int t = (i+1)%2;
                #pragma omp for
                for(int j = t; j<k-1;j+=2){
                    //fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid, i,j);
                    int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
                    current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                }
                #pragma omp single
                {
                  // If number of columns is odd, update last column of odd rows
                  if(t == 0){
                    int j = k-1;
                    int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
                    current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                  }
                }
            }
            // Update bottom row of frame (only white positions)
            #pragma omp for
            for(int j=0;j<k;j+=2){
                current[(k+1)*k+j] = current[k+j];
            }

            // Update last inner row of current and upper row of frame
            // Update white positions (consider that the first cell is white)
            #pragma omp for
            for(int j = 0; j < k-1; j+=2){
              int i = k;
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
              current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              // Frame
              current[j] = current[i*k+j];
            }
            {
              // Update last element of last innser row
              int i = k;
              int j = k-1;
              int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
              current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              // Frame
              current[j] = current[i*k+j];
            }
            // Update black positions (consider that the first black cell is the second one of the row)
            for(int i=1;i<k;i++){
                int t = i%2;
                #pragma omp for
                for(int j=t; j<k-1;j+=2){
                    //fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid, i,j);
                    int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
                    current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                }
                #pragma omp single
                {
                  // If number of columns is odd, update last column of odd rows
                  if(t == 0){
                    int j = k-1;
                    int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
                    current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                  }
                }
                
            }
            // Update bottom row of frame (only black positions)
            #pragma omp for
            for(int j=1;j<k;j+=2){
                current[(k+1)*k+j] = current[k+j];
            }

            // Update last inner row of current and upper row of frame
            // Update black positions (consider that the first black cell is the second one of the row)
            #pragma omp for
            for(int j = 1; j< k; j+=2){
                int i = k;
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                        current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                        current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                // Frame
                current[j] = current[i*k+j];
            }
        }
        if((n_step+1) % s == 0){
            char file_path[45] = "images/evolve_black_white/"; // Sufficiently large
            char filename[20];
            
            snprintf(filename, 20, "snapshot_%05d.pgm", n_step+1);
            strcat(file_path, filename);
            write_pgm_parallel(current+k, 255, k, k, file_path, rank, size, rows_read);
        }
    }
}

void evolve_black_white(unsigned char *current, int k, int n_steps, int s){
  if(k%2 == 0){
    evolve_black_white_even(current, k, n_steps, s);
  }else{
    evolve_black_white_odd(current, k, n_steps, s);
  }

}

