// To compile: gcc main.c -o main.exe
// To run executable to generate playground: ./main.exe -i -k 5 -f init.pgm
// 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

void write_pgm_image( unsigned char *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name);

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1


char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;

void random_playground(int k, char *fname){
    /*
    This function creates the playground and saves it to the pgm file.
    Always read the matrix from the pgm file to use it.
    */
  unsigned char* ptr = (unsigned char*)calloc(k*k, sizeof(unsigned char)); // creates a k*k array of unsigned char

  // generate a random matrix of 0 and 255
  for (int i = 0; i < k*k; i++) {
    unsigned char rand_num = (unsigned char) rand() % 2;
    if(rand_num==1){
      ptr[i] = 255;
    }else{
      ptr[i]=rand_num;
    }
  }
  write_pgm_image(ptr, 255, k, k,fname);
  free(ptr); 
}

int main ( int argc, char **argv )
{
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
    // run a play ground
    printf("Run\n");
  }

  int xsize;
  int ysize;
  int maxval;
  unsigned char* myimage; 
  //void* myimage;
  // myimage is a pointer so with &myimage I access the address of the pointer => pointer of a pointer
  read_pgm_image(&myimage, &maxval, &xsize, &ysize, fname);
  
  if ( fname != NULL ){
    free ( fname );
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
  k = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, image_file);
  while ( (k > 0) && (line[0]=='#') )
    k = getline( &line, &n, image_file);

  if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
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
      printf("Here\n");
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

