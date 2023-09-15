#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
/*#include <mpi.h>
#include <omp.h>
*/
// To run: srun gcc -fopenmp black_white.c -o black_white.exe
// To execute: srun ./black_white.exe 5


void print_image(unsigned char* ptr, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
            //printf("%d ", ptr[i*ncol + j]);
        }
        printf("\n");
    }
}

void initialize_matrix(unsigned char *current, int k){
  // generate a random matrix of 0 and 255 and create its rows frame

  /*unsigned int seed = clock();
  for (int i = 0; i < k*k; i++) {
    unsigned char rand_num = (unsigned char) rand_r(&seed) % 2;
    current[i+k] = rand_num==1 ? 255 : 0;
  }*/
  int values[] = {
        255, 255, 0, 255, 255,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 255,
        0, 0, 0, 0, 0,
        255, 255, 255, 0, 255
    };
    for(int i = 0; i<k*k;i++){
        current[i+k]=values[i];
    }

    // Last row of current -> First row of current
  for (int j = 0; j < k; j++){
    current[(k+1)*k+j] = current[j+k];
    //current[i] = current[(k-1)*k + i];
  }
  // First row of current -> Last row of current
  for (int j = 0; j < k; j++){
    current[j]=current[k*k+j];
  }  
}

/*
void evolve_black_white_parallel(unsigned char *current,  int k, int n_steps){
    int nthreads;
    FILE* prova_file;
    char nome_file[] = "prova_file.txt";
    prova_file=fopen(nome_file, "w");
    
    for(int n_step=0; n_step < n_steps; n_step++){
        #pragma omp parallel
        {
            int myid = omp_get_thread_num();
            // Update black positions (consider that the first cell is black)
            for(int i=1;i<k;i++){
                #pragma omp for
                for(int j=0; j<k;j++){
                    fprintf(prova_file,"I am thread %d and I am updating element (%d,%d)\n", myid, i,j);
                    if((i-1+j)%2 == 0){
                        int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                            current[(i-1)*k + j] + current[(i+1)*k + j];
                        current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
                    }
                }
            }
            // Update white positions (consider that the first white cell is the second one of the row)
            for(int i=1;i<k;i++){
                #pragma omp for
                for(int j=0; j<k;j++){
                    if((i-1+j)%2 == 1){
                        int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                            current[(i-1)*k + j] + current[(i+1)*k + j];
                        current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
                    }
                }
            }
            // Update bottom row of frame
            #pragma omp for
            for(int j=0;j<k;j++){
                current[(k+1)*k+j] = current[k+j];
            }

            // Update last inner row of current and upper row of frame
            #pragma omp for
            for(int j = 0; j< k; j++){
                if((i-1+j)%2 == 0){
                    int i = k;
                    int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                            current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                            current[(i-1)*k + j] + current[(i+1)*k + j];
                    current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
                    // Frame
                    current[j] = current[i*k+j];
                }
            }
            #pragma omp for
            for(int j = 0; j< k; j++){
                if((i-1+j)%2 == 1){
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
}
*/

void evolve_black_white_serial(unsigned char *current, int k, int n_steps){
    for(int n_step=0; n_step < n_steps; n_step++){
            // Update black positions (consider that the first cell is black)
        printf("Black inner positions except last row:\n");
        for(int i=1;i<k;i++){
            int t = (i+1)%2;
            for(int j = t; j<k;j+=2){
                printf("(%d, %d) ", i,j);
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
            }
        }
        printf("\n");
        // Update bottom row of frame (only black positions)
        printf("Update bottom row of frame (only black positions)\n");
        for(int j=0;j<k;j+=2){
            int u = k+j;
            printf("(%d,%d) ",u,j);
            current[(k+1)*k+j] = current[k+j];
        }
        printf("\n");
        printf("White inner positions except last row:\n");
        // Update white positions (consider that the first white cell is the second one of the row)
        for(int i=1;i<k;i++){
            int t = i%2;
            for(int j=t; j<k;j+=2){
                printf("(%d,%d) ", i,j);
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
            }
        }
        printf("\n");
        // Update bottom row of frame (only white positions)
        printf("Update bottom row of frame (only white positions)\n");
        for(int j=1;j<k;j++){
            int u = k+j;
            printf("(%d,%d) ", u,j);
            current[(k+1)*k+j] = current[k+j];
        }
        printf("\n");

        // Update last inner row of current and upper row of frame
        // Update black positions (consider that the first cell is black)
        printf("Last inner row (black positions)\n");
        for(int j = 0; j< k; j+=2){
            int i = k;
            printf("(%d,%d) ", i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            // Frame
            current[j] = current[i*k+j];
        }
        printf("\n");
        // Update white positions (consider that the first white cell is the second one of the row)
        printf("Last inner row (white positions)\n");
        for(int j = 1; j< k; j+=2){
            int i = k;
            printf("(%d,%d) ", i,j);
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            // Frame
            current[j] = current[i*k+j];
        }
        printf("\n");

    }
}


int main(int argc, char** argv){
    int k;
    if (argc > 1)
        k = atoi(argv[1]);
    else
        k = 5;

    unsigned char* current = (unsigned char*)malloc((k+2)*k*sizeof(unsigned char));
    
    initialize_matrix(current, k);

    print_image(current,k+2,k);

    evolve_black_white_serial(current, k, 1);
    printf("Matrix after 1 step of update:\n");
    print_image(current,k+2,k);

    free(current);
    return 0;
}