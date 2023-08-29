#include <stdio.h>
#include <stdlib.h>
/*
In all these functions we have that k is the dimension of the matrices, which are assumed to be square.
*/

void print_image(unsigned char* ptr, int ncol){
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]);
        }
        printf("\n");
    }
}
// Ora faccio quella con la cornice totale

void initialize_current(unsigned char* input, unsigned char* current, int k){

    // Initialize the innter values
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            current[(i+1)*(k+2) + (j+1)] = input[i*k + j];
        }
    }

    // Initialize the corners and the frame rows and columns
    // First row of input -> Last row of current
    for (int i = 0; i < k; i++){
        current[(k+1)*(k+2) + i + 1]=input[i];
    }

    // Last row of input -> First row of current
    for (int i = 0; i < k; i++){
        current[i+1] = input[(k-1)*k + i];
    }

    // Last column of input -> First column of current
    for (int i = 0; i < k; i++){
        current[(i+1)*(k+2)] = input[(i+1)*k - 1];
    }

    // First column of input -> Last column of current
    for(int i = 0; i < k; i++){
        current[(i+2)*(k+2)-1] = input[i*k];
    }

    // Now add the corners
    current[0] = input[k*k-1]; // Top left of current
    current[(k+2) - 1] = input[k*(k-1)]; // Top right of current
    current[(k+1)*(k+2)] = input[k-1]; // Bottom left of current
    current[(k+2)*(k+2)-1] = input[0]; // Bottom right of current
}

void evolve_static(unsigned char* current, unsigned char* next, int k){
    // We work with i and j referring to the extended matrix unlike before

    // Start by updating the internal values
    for(int i = 1; i < k+2; i++){
        for(int j = 1; j < k+2; j++){
            int alive_neighbours = current[(i+1)*k + j] + current[(i-1)*k + j] + current[(i)*k + (j+1)] + 
                                   current[(i)*k + (j-1)] + current[(i+1)*k + (j+1)] + current[(i+1)*k + (j-1)] + 
                                   current[(i-1)*k + (j-1)] + current[(i-1)*k + (j+1)];
            if (alive_neighbours > 3 | alive_neighbours < 2)
                next[i*k + j] = 0; // Dead
            else
                next[i*k + j] = 1; // Alive
        }
    }

    // Now update the frame accordingly
    // First row
    for(int i = 1; i < k+2; i++){
        next[i] = next[(k)*(k+2) + i]; // row k is the last inner row
    }

    // Last row
    for(int i = 1; i < k+2; i++){
        next[(k+1)*(k+2) + i] = next[(k+2) + i];
    }

    // First column
    for(int i = 1; i < k+2; i++){
        next[(k+2)*i] = next[(k+2)*(i+1) - 2];
    }

    // Last column
    for(int i = 1; i < k+2; i++){
        next[(k+2)*(i+1) - 1] = next[(i)*(k+2) + 1];
    } 

    // Update the corners
    next[0] = next[(k+1)*(k+2)-2];          // Top left corner
    next[(k+2)-1] = next[(k)*(k+2)+1];      // Top right corner
    next[(k+1)*(k+2)] = next[2*(k+2)-2];    // Bottom left corner
    next[(k*2)*(k*2)-1] = next[(k+2)+1];    // Bottom right corner
}


int main(int argc, char* argv[]){
    int n_steps = 5;
    int k = 3;

    printf("HI\n");
    // Take the input
    unsigned char* input = (unsigned char*)calloc(k*k, sizeof(unsigned char));
    input[8] = 1;
    input[1] = 1;
    input[5] = 1;
    input[7] = 1;
    input[4] = 1;

    // Allocate memory for new ones, and create frame
    unsigned char* current = (unsigned char*)calloc((k+2)*(k+2), sizeof(unsigned char)); 
    unsigned char* next = (unsigned char*)malloc((k+2)*(k+2)*sizeof(unsigned char));
    
    print_image(input, k);

    printf("\nNow we initialize by adding the frame\n\n");

    initialize_current(input, current, k);
    print_image(current, k+2);
    printf("\nNow we evolve\n\n");

    evolve_static(current, next, k);
    print_image(next, k+2);


    free(input);
    free(current);
    free(next);
}

