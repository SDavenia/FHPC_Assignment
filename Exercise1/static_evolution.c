#include <stdio.h>
#include <stdlib.h>

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
    unsigned char* current = (unsigned char*)calloc((k+1)*(k+1), sizeof(unsigned char)); 
    unsigned char* next = (unsigned char*)malloc((k+1)*(k+1)*sizeof(unsigned char));
    
    print_image(input, k);

    printf("\nNow we initialize\n\n");

    initialize_current(input, current, k);
    print_image(current, k+2);

    free(input);
    free(current);
    free(next);

}

