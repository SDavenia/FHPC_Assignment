#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "black_white.h"


void evolve_black_white(unsigned char *current, int k, int n_steps){
  int nthreads;
  for(int n_step=0; n_step < n_steps; n_step++){
      #pragma omp parallel
      {
          int myid = omp_get_thread_num();
          // Update white positions (consider that the first cell is white)
          /*
          AGGIUNGERE CONTROLLO: se il numero di colonne k della matrice è pari,
          prima di aggiornare l'ultima colonna è necessario che la colonna 0 sia stata aggiornata
          */
          for(int i=1;i<k;i++){
              int t = (i+1)%2;
              #pragma omp for
              for(int j = t; j<k;j+=2){
                  int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                      current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                      current[(i-1)*k + j] + current[(i+1)*k + j];
                  current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
              }
              #pragma omp master
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
  }
}

