#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*
Only problem is that if we run it consecutively (run two times ./a.out) we get basically the same random generated numbers.
*/

int main()
{
    unsigned int seed = time(NULL);    // seed taken from the clock for rand_r
    double p = 0.8;                    // probability of having 1

    double random_01;                  // store the uniformly generated number
    unsigned char pixel;
    for(int i = 0; i < 10; i++){
        random_01 = (double)rand_r(&seed) / (double)RAND_MAX ; 
        printf("generated %f\n", random_01);
        pixel = random_01 <= p ? 1 : 0;
        printf("Coincides with %u\n", pixel); 
    }        
    
}
