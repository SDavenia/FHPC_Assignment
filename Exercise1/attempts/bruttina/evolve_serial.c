#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
/*
THIS IS THE CORRECT VERSION FOR THE SERIAL CODE, with frame only above and below.
Use this as reference to check whether the Conway evolution was correct.
*/

void read_pgm_image( unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name);
void initialize_current(unsigned char* input, unsigned char* current, int k);
void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps);
void evolve_dynamic(unsigned char* current, int k, int n_steps);

void print_image(unsigned char* ptr, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]){
    printf("Run\n");
    int xsize;
    int ysize;
    int maxval;
    char fname[] = "init.pgm";
    unsigned char* input;
    read_pgm_image(&input, &maxval, &xsize, &ysize, fname);
    int k;
    if(argc > 1)
      k = atoi(argv[1]);
    else
      k = 5;
    printf("Input:\n");
    print_image(input,k,k);

    unsigned char* current = (unsigned char*)malloc((k+2)*k*sizeof(unsigned char));
    initialize_current(input, current, k);
    printf("\nCurrent:\n");
    print_image(current,k+2,k);


    //if(e == 0){ // Ordered
    printf("ORDERED EXECUTION\n");
    evolve_dynamic(current, k, 1);
    //}else{ // Static
    //printf("STATIC EXECUTION\n");
    //unsigned char* next = (unsigned char*)malloc((k+2)*k*sizeof(unsigned char));
    //printf("Allocated memory for next using malloc\n");

    //evolve_static(current, next, k, 100);
    //printf("Finished static evolution\n");
    free(current);
    //free(next);
    //}

    return 0;
}




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
  printf("\n");
  
  fclose(image_file);
  return;
}

void initialize_current(unsigned char* input, unsigned char* current, int k){

  // Last row of input -> First row of current
  for (int i = 0; i < k; i++){
    current[i] = input[(k-1)*k + i];
  }
  
  // Initialize the inner values
  for(int i = 0; i < k; i++){
    for(int j = 0; j < k; j++){
        current[(i+1)*k+j] = input[i*k+j];
    }
  }
  // Initialize the corners and the frame rows and columns
  // First row of input -> Last row of current
  for (int i = 0; i < k; i++){
    current[(k+1)*k+i]=input[i];
  }
}

void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps){
    for(int n_step=0; n_step < n_steps; n_step++){

        for(int i=1;i<k+1;i++){
            for(int j=0; j<k;j++){
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                next[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
            }
        }
        for (int i = 0; i < k; i++){
            next[i] = next[(k)*k + i];
        }
        for (int i = 0; i < k; i++){
            next[(k+1)*k+i]=next[k+i];
        }
        unsigned char* tmp;
        tmp = next;
        next = current;
        current=tmp;

        if(n_step == 50){
          printf("Step %d:\n", n_step);
          print_image(current, k+2,k);
        }
    }

}

void evolve_dynamic(unsigned char* current, int k, int n_steps){
    for(int n_step=0; n_step < n_steps; n_step++){
        
        for(int i=1;i<k;i++){
            for(int j=0; j<k;j++){
                int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
                current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255; 
            }
        }
        // Update bottom row of frame
        for(int j=0;j<k;j++){
            current[(k+1)*k+j] = current[k+j];
        }

        // Update last inner row of current and upper row of frame
        for(int j = 0; j< k; j++){
            int i = k;
            int n_neigh = current[(j-1 + k)%k + i*k] + current[(j+1 + k)%k + i*k] + current[(j-1 + k)%k + (i-1)*k] +
                    current[(j+1 + k)%k + (i-1)*k] + current[(j-1 + k)%k + (i+1)*k] + current[(j+1 + k)%k + (i+1)*k]+
                    current[(i-1)*k + j] + current[(i+1)*k + j];
            current[i*k+j] = (n_neigh > 765 || n_neigh < 510) ? 0 : 255;
            // Frame
            current[j] = current[i*k+j];
        }
        if(n_step==0){
          printf("Step %d of dynamic:\n", n_step);
          print_image(current, k+2,k);
        }
    }
}