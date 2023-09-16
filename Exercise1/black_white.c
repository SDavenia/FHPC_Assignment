#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
//#include <mpi.h>
#include <omp.h>

// To run: srun gcc -fopenmp black_white.c -o black_white.exe
// To execute: srun ./black_white.exe 5


void read_pgm_image( unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name)
/*
 * image        : a pointer to the pointer that will contain the image
 * maxval       : a pointer to the int that will store the maximum intensity in the image
 * xsize, ysize : pointers to the x and y sizes
 * image_name   : the name of the file to be read
 *
 */
/*
"image" is a pointer to a pointer so:
- with *image I access the address of the first element of the string (image)
- with **image I access the value of the first element of the string (image)
*/
{
  FILE* image_file; 
  image_file = fopen(image_name, "r");

  *image = NULL; //address of the first element of image
  *xsize = *ysize = *maxval = 0; // set to 0 the value of xsize, ysize and maxval
  
  char    MagicN[2]; // define a string of 2 elements
  char   *line = NULL; //define a pointer "line" to NULL
  size_t  k, n = 0;
  
  // get the Magic Number
  k = fscanf(image_file, "%2s%*c", MagicN ); // This one reads P5

  // skip all the comments
  k = getline( &line, &n, image_file); // Here we read all the lines starting with #, i.e. all the comments.
  while ( (k > 0) && (line[0]=='#') )
    k = getline( &line, &n, image_file);

  if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);  // This one reads the number
      if ( k < 3 )
	fscanf(image_file, "%d%*c", maxval);
    }
  else
    {
      *maxval = -1;         // this is the signal that there was an I/O error
			    // while reading the image header
      free( line );
      return;
    }
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
  unsigned int size = *xsize * *ysize * color_depth;
  
  if ( (*image = (unsigned char*)malloc( size )) == NULL )
    {
      fclose(image_file);
      *maxval = -2;         // this is the signal that memory was insufficient
      *xsize  = 0;
      *ysize  = 0;
      return;
    }
  
  if ( fread( *image, 1, size, image_file) != size )
    {
      free( image );
      image   = NULL;
      *maxval = -3;         // this is the signal that there was an i/o error
      *xsize  = 0;
      *ysize  = 0;
    }  
/*
  for (int u = 0; u < 9; u++) {
    //printf("%u\t", (*image)[u]);
    printf("%u\t", *(*image+u));
  }
*/
  
  fclose(image_file);
  return;
}

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
  /*int values[] = {
        255, 255, 0, 255, //255,
        0, 0, 0, 0, //0,
        0, 0, 0, 0, //255,
        0, 0, 0, 0, //0,
        255, 255, 255, 0, //255
    };*/
    int values[] = {
        0, 255, 0, 0, 255,
        255, 0, 255, 255, 0,
        255, 0, 0, 255, 0,
        0, 0, 0, 0, 0,
        255, 255, 0, 255, 255
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


void evolve_black_white_parallel_even(unsigned char *current, int k, int n_steps){
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

void evolve_black_white_parallel_odd(unsigned char *current, int k, int n_steps){
    int nthreads;
    /*FILE* prova_file;
    char nome_file[] = "prova_file.txt";
    prova_file=fopen(nome_file, "w");
    */
    for(int n_step=0; n_step < n_steps; n_step++){
        #pragma omp parallel
        {
            int myid = omp_get_thread_num();
            // Update black positions (consider that the first cell is black)
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
            // Update bottom row of frame (only black positions)
            #pragma omp for
            for(int j=0;j<k;j+=2){
                current[(k+1)*k+j] = current[k+j];
            }

            // Update last inner row of current and upper row of frame
            // Update black positions (consider that the first cell is black)
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
            // Update white positions (consider that the first white cell is the second one of the row)
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
            // Update bottom row of frame (only white positions)
            #pragma omp for
            for(int j=1;j<k;j+=2){
                current[(k+1)*k+j] = current[k+j];
            }

            // Update last inner row of current and upper row of frame
            // Update white positions (consider that the first white cell is the second one of the row)
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

      printf("Step %d\n", n_step);
      print_image(current,k+2,k);
    }
}

void evolve_black_white_parallel(unsigned char *current, int k, int n_steps){
  if(k%2 == 0){
    evolve_black_white_parallel_even(current, k, 5);
  }else{
    evolve_black_white_parallel_odd(current, k, 5);
  }

}

void evolve_black_white_serial(unsigned char *current, int k, int n_steps){
    for(int n_step=0; n_step < n_steps; n_step++){
        // Update black positions (consider that the first cell is black)
        for(int i=1;i<k;i++){
            int t = (i+1)%2;
            for(int j = t; j<k;j+=2){
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            }
        }             

  // Update bottom row of frame (only black positions)
        for(int j=0;j<k;j+=2){
            current[(k+1)*k+j] = current[k+j];
        }

        // Update last inner row of current and upper row of frame
        // Update black positions (consider that the first cell is black)

        for(int j = 0; j< k; j+=2){
            int i = k;
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            // Frame
            current[j] = current[i*k+j];
        }
        
        // Update white positions (consider that the first white cell is the second one of the row)
        for(int i=1;i<k;i++){
            int t = i%2;
            for(int j=t; j<k;j+=2){
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            }
        }
        // Update bottom row of frame (only white positions)
        for(int j=1;j<k;j+=2){
            current[(k+1)*k+j] = current[k+j];
        }

        // Update last inner row of current and upper row of frame
        // Update white positions (consider that the first white cell is the second one of the row)
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

void initialize_current(unsigned char* input, unsigned char* current, int k){

  // First row of current -> Last row of current
  for (int j = 0; j < k; j++){
    current[(k+1)*k+j] = current[j+k];
  }
  // Last row of current -> First row of current
  for (int j = 0; j < k; j++){
    current[j]=current[k*k+j];
    }
}


int main(int argc, char** argv){
    int k;
    if (argc > 1)
        k = atoi(argv[1]);
    else
        k = 5;
    int xsize;
    int ysize;
    int maxval;  
    
    /*unsigned char* input;
    char *file_path = (char*)malloc( sizeof(optarg)+ 1 +30);
    strcpy(file_path,"images/initial_matrices/init_20000.pgm");
    read_pgm_image(&input, &maxval, &xsize, &ysize, file_path);
    

    unsigned char* current = (unsigned char*)malloc((k+2)*k*sizeof(unsigned char));
    initialize_current(input, current, k);
    free(input);*/

    
    unsigned char* current = (unsigned char*)malloc((k+2)*k*sizeof(unsigned char));
    initialize_matrix(current, k);
    print_image(current,k+2,k);
    evolve_black_white_parallel(current, k, 10);
    /*printf("Matrix after 1 step of update:\n");
    print_image(current,k+2,k);*/
    

    /*double Tstart_bw;
    Tstart_bw = omp_get_wtime();
    evolve_black_white_parallel(current, k, 1);
    double Time_bw = omp_get_wtime() - Tstart_bw;
    printf("Black and white time: %lf\n", Time_bw);*/
    

    free(current);
    //free(file_path);
    return 0;
}