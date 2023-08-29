
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#ifndef READ_WRITE_PGM
#define READ_WRITE_PGM
void write_pgm_image( unsigned char *image, int maxval, int xsize, int ysize, const char *image_name);
#endif


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

void init_matrix(int k, char *fname){
  unsigned char* ptr = (unsigned char*)calloc(k*k, sizeof(unsigned char)); // creates a k*k array of unsigned char

  for (int i = 0; i < k*k; i++) {
    ptr[i] = (unsigned char) rand() % 2;
  }
  // Just to check if the matrix is ok
  for (int j = 0; j < k*k; j++) {
    printf("%u", ptr[j]);
  }

  write_pgm_image(ptr, 255, k, k,fname);

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
      fname = (char*)malloc( sizeof(optarg)+1 );
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
    init_matrix(k,fname);
  }else{
    // run a play ground
    printf("Run\n");
  }

    if ( fname != NULL )
      free ( fname );

  return 0;
}
