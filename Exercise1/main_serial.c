// FINAL VERSION OF THE OPTIMIZED MAIN (ONLY SERIAL FUNCTIONS)
// To compile: gcc main.c -o main.exe
// To run executable to generate playground: ./main.exe -i -k 5 -f init.pgm
// To run execubtable to play on playground: ./main.exe -r -k 5 -f init.pgm -n 3
// OPTIMIZED and with the errors fixed
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

void write_pgm_image( unsigned char *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name);

void initialize_current(unsigned char* input, unsigned char* current, int k);
void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps);
void evolve_ordered(unsigned char* current, int k, int n_steps);

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9)


char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;


void print_image(unsigned char* ptr, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
        }
        printf("\n");
    }
}

void random_playground(int k, char *fname){
    /*
    This function creates the playground and saves it to the pgm file.
    Always read the matrix from the pgm file to use it.
    */
  unsigned char* ptr = (unsigned char*)calloc(k*k, sizeof(unsigned char)); // creates a k*k array of unsigned char
  // generate a random matrix of 0 and 255
  unsigned int seed = clock();

  for (int i = 0; i < k*k; i++) {
    unsigned char rand_num = (unsigned char) rand_r(&seed) % 2;
    ptr[i] = rand_num==1 ? 255 : 0;
    /*if(rand_num==1){
      ptr[i] = 255;
    }else{
      ptr[i]=rand_num;
    }*/
  }
  write_pgm_image(ptr, 255, k, k,fname);
  free(ptr); 
}

int main ( int argc, char **argv )
{
  struct timespec ts;
  int action = 0;
  char *optstring = "irk:e:f:n:s:"; //optstring is a list of characters, each representing a single character option

  int c;
  /*
  Generally, the getopt() function is called from inside of a loop’s conditional statement.
  The loop terminates when the getopt() function returns -1.
  A switch statement is then executed with the value returned by getopt() function.
  */
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch(c) {
      //c takes the current option (ex. k)
      // optarg takes the value of the current option (ex. 400 for k)
      
    case 'i':
      action = INIT; break;
      
    case 'r':
      action = RUN; break;
      
    case 'k':
      k = atoi(optarg); break;

    case 'e':
      e = atoi(optarg); break;

    case 'f':
      fname = (char*)malloc( sizeof(optarg)+1 ); //fname now is an "array of char" variable (a string)
      sprintf(fname, "%s", optarg );
      printf("%s \n",fname);
      break;

    case 'n':
      n = atoi(optarg); break;

    case 's':
      s = atoi(optarg); break;

    default :
      printf("argument -%c not known\n", c ); break;
    }
  }

  if(action == INIT){
    // create initial conditions
    printf("Initialize\n");
    random_playground(k,fname);
  }else{ 
    // Read and run a playground
    printf("Run\n");
    int xsize;
    int ysize;
    int maxval;  
    
    unsigned char* input;
    read_pgm_image(&input, &maxval, &xsize, &ysize, fname);
 
    printf("Initial image\n");
    print_image(input, k, k);
    printf("INITALIZING THE FRAME\n");
    unsigned char* current = (unsigned char*)malloc((k+2)*(k+2)*sizeof(unsigned char));
    printf("Allocated memory for current using malloc\n");

    double Tstart_init = CPU_TIME;
    initialize_current(input, current, k);
    double Time_init = CPU_TIME - Tstart_init;
    printf("Initializing time: %lf\n", Time_init);

    printf("FREE USELESS STUFF\n");
    free(input); // Since we do not need it anymore
    if (fname != NULL)
        free(fname);

    double Tstart_exec = CPU_TIME;
    if(e == 0){ // Ordered
        printf("ORDERED EXECUTION\n");
        evolve_ordered(current, k, n);
    }else{ // Static
        printf("STATIC EXECUTION\n");
        unsigned char* next = (unsigned char*)malloc((k+2)*(k+2)*sizeof(unsigned char));
        printf("Allocated memory for next using malloc\n");

        evolve_static(current, next, k, n);
        printf("Finished static evolution");
        free(next);
    }
    double Time_exec = CPU_TIME - Tstart_exec;
    printf("Execution time: %lf\n", Time_exec);
    free(current);
  }

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

void write_pgm_image( unsigned char *image, int maxval, int xsize, int ysize, const char *image_name)
/*
 * image        : a pointer to the memory region that contains the image
 * maxval       : either 255 or 65536
 * xsize, ysize : x and y dimensions of the image
 * image_name   : the name of the file to be written
 *
 */
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1 + ( maxval > 255 ); // 1+0 if maxval <= 255; 1+1 if maxval > 255

  fprintf(image_file, "P5\n# generated by\n# Elena Rivaroli and Samuele D'Avenia\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, 1, xsize*ysize*color_depth, image_file);
  /*
  size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
    - ptr − This is the pointer to the array of elements to be written.
    - size − This is the size in bytes of each element to be written.
    - nmemb − This is the number of elements, each one with a size of size bytes.
    - stream − This is the pointer to a FILE object that specifies an output stream.
  */  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
           ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
  ------------------------------------------------------------------ */
}

void initialize_current(unsigned char* input, unsigned char* current, int k){

  current[0] = input[k*k-1]; // Top left of current
  // Last row of input -> First row of current
  for (int i = 0; i < k; i++){
    current[i+1] = input[(k-1)*k + i];
  }
  current[(k+2) - 1] = input[(k-1)*k]; // Top right of current


  // Initialize the inner values
  for(int i = 0; i < k; i++){
    // First column of that row
    current[(i+1)*(k+2)] = input[(i+1)*k - 1];
    for(int j = 0; j < k; j++){
      current[(i+1)*(k+2) + (j+1)] = input[i*k + j];
    }
    // Last column of that row
    current[(i+2)*(k+2)-1] = input[i*k];
  }

  current[(k+1)*(k+2)] = input[k-1]; // Bottom left of current

  // Initialize the corners and the frame rows and columns
  // First row of input -> Last row of current
  for (int i = 0; i < k; i++){
    current[(k+1)*(k+2) + i + 1]=input[i];
  }
  current[(k+2)*(k+2)-1] = input[0]; // Bottom right of current
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

void evolve_ordered(unsigned char* current, int k, int n_steps){
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
          printf("Step %d of ordered:\n", n_step);
          print_image(current, k+2,k);
        }
    }
}