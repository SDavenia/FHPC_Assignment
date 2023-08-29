#include <stdio.h>
#include <stdlib.h>
/*
In all these functions we have that k is the dimension of the matrices, which are assumed to be square.
*/

void print_image(unsigned char* ptr, int ncol){
    for(int i = 0; i < ncol; i++){
        for(int j = 0; j < ncol; j++){
            printf("%d ", ptr[i*ncol + j]/255);
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


void evolve_static(unsigned char* current, unsigned char* next, int k, int n_steps){
    for (int n = 0; n < n_steps; n++){

        // Start by updating the internal values
        for(int i = 1; i < k+1; i++){
            for(int j = 1; j < k+1; j++){
                int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                    current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                    current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
                alive_neighbours /= 255;
                //printf("I am (%d, %d) and I have %d\n", i, j,  alive_neighbours);
                if (alive_neighbours > 3 || alive_neighbours < 2)
                    next[i*(k+2) + j] = 0; // Dead
                else
                    next[i*(k+2) + j] = 255; // Alive
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


        // Swap the pointers so that you have the right one
        unsigned char* tmp;
        // printf("next is at %p\n", next);
        // printf("current is at %p\n", current);
        tmp = next;
        next = current;
        current = tmp;
        // printf("next is at %p\n", next);
        // printf("current is at %p\n", current);

        // printf("PRINTING CURRENT\n");
        printf("Result after iteration %d:\n", n+1);
        print_image(current, k+2);
    }
    
}



void evolve_dynamic(unsigned char* current, int k, int n_steps){

    for (int n = 0 ; n < n_steps; n++){
        // Update first inner row and consequently update last row

        // First entry should also update last column
        {int i = 1;
        int j = 1;
        int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
        
        alive_neighbours /= 255;
        if (alive_neighbours > 3 || alive_neighbours < 2){
            current[i*(k+2) + j] = 0; // Dead
            current[(k+1)*(k+2)+j] = 0; // Update last row as well
            current[2*(k+2)-1] = 0; // Update last column as well
        }else{
            current[i*(k+2) + j] = 255; // Alive
            current[(k+1)*(k+2)+j] = 255; // Update last row as well
            current[2*(k+2)-1] = 255; // Update last column as well
        }
        }

        // Now all other entries
        for(int j = 2; j < k+1; j++){
            int i = 1;

            int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
            
            alive_neighbours /= 255;
            if (alive_neighbours > 3 || alive_neighbours < 2){
                current[i*(k+2) + j] = 0; // Dead
                current[(k+1)*(k+2)+j] = 0; // Update last row as well
            }else{
                current[i*(k+2) + j] = 255; // Alive
                current[(k+1)*(k+2)+j] = 255; // Update last row as well
            }
        }
        // You can now also update the first row entry (not needed to be done before since not used immediately)
        current[(k+2)] = current[2*(k+2)-2];

        // Now that they are correct, update the bottom corners using the values you just calculated 
        // No need to update them before since they are only needed later
        current[(k+1)*(k+2)] = current[(k+2)*2 - 2];   // Bottom left corner
        current[(k+2)*(k+2)-1] = current[(k+2) + 1]; // Bottom right corner
        

        // Update all the other rows (from second to penultimate), keep in mind that whenever first or last value of the row get updated, 
        //  they should also be propagated to last and first column respectively
        for(int i = 2; i < k; i++){
            
            // First inner entry so we update last column
            {int j = 1;
            int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
            alive_neighbours /= 255;
            if (alive_neighbours > 3 || alive_neighbours < 2){
                current[i*(k+2) + j] = 0;
                current[(i+1)*(k+2) - 1] = 0;  // Update last column as well
            }else{
                current[i*(k+2) + j] = 255; // Alive
                current[(i+1)*(k+2) - 1] = 255; // Update last column as well
            }                 
            }
            for (int j = 2; j < k; j++){
                int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
                alive_neighbours /= 255;
                if (alive_neighbours > 3 || alive_neighbours < 2)
                    current[i*(k+2) + j] = 0;
                else
                    current[i*(k+2) + j] = 255; // Alive  
            }
            // Finally last inner entry so we update last first column as well
            {int j = k;
            int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
            alive_neighbours /= 255;
            if (alive_neighbours > 3 || alive_neighbours < 2){
                current[i*(k+2) + j] = 0;
                current[i*(k+2)] = 0;  // Update first column as well
            }else{
                current[i*(k+2) + j] = 255; // Alive
                current[i*(k+2)] = 255; // Update first column as well
            }  
            }
        }

        
        // Update last inner row and consequently update first row
        // First entry should also update last column
        {int i = k;
        int j = 1;
        int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
        
        alive_neighbours /= 255;
        if (alive_neighbours > 3 || alive_neighbours < 2){
            current[i*(k+2) + j] = 0;
            current[j] = 0;  // Update first row as well
            current[(k+1)*(k+2)-1] = 0;// Update last column as well
        }else{
            current[i*(k+2) + j] = 255; // Alive
            current[j] = 255; // Update first row as well
            current[(k+1)*(k+2)-1] = 255;
        }  
        }

        for(int j = 2; j < k+1; j++){
            int i = k;

            int alive_neighbours = current[(i+1)*(k+2) + j] + current[(i-1)*(k+2) + j] + current[(i)*(k+2) + (j+1)] + 
                                        current[(i)*(k+2) + (j-1)] + current[(i+1)*(k+2) + (j+1)] + current[(i+1)*(k+2) + (j-1)] + 
                                        current[(i-1)*(k+2) + (j-1)] + current[(i-1)*(k+2) + (j+1)];
            alive_neighbours /= 255;
            if (alive_neighbours > 3 || alive_neighbours < 2){
                current[i*(k+2) + j] = 0;
                current[j] = 0;  // Update first row as well
            }else{
                current[i*(k+2) + j] = 255; // Alive
                current[j] = 255; // Update first row as well
            }  
        }
        // You can now update the first row entry 
        current[(k)*(k+2)] = current[(k+1)*(k+2)-2];
        
        // Now that they are correct, update the top corners using the values you just calculated 
        // No need to update them before since they are only needed later
        current[0] = current[(k+1)*(k+2)-2];   // Top left corner
        current[(k+2)-1] = current[(k*(k+2))+1]; // Top right corner

        printf("Result after %d steps:\n", n+1);
        print_image(current, k+2);
    }
}


int main(int argc, char* argv[]){
    int n_steps = 2;
    int k = 4;

    printf("HI\n");
    // Take the input
    unsigned char* input = (unsigned char*)calloc(k*k, sizeof(unsigned char));
    input[8] = 255;
    input[1] = 255;
    input[5] = 255;
    input[7] = 255;
    input[4] = 255;

    // Allocate memory for new ones, and create frame
    unsigned char* current = (unsigned char*)calloc((k+2)*(k+2), sizeof(unsigned char)); 
    unsigned char* next = (unsigned char*)malloc((k+2)*(k+2)*sizeof(unsigned char));
    
    //print_image(input, k); // Initial one without frame
    initialize_current(input, current, k);
    
    printf("Below is the initial playground with the extended frame:\n");
    print_image(current, k+2); // Initial one with frame

    // STATIC EVOLUTION
    // evolve_static(current, next, k, n_steps);

    // DYNAMIC EVOLUION
    evolve_dynamic(current, k, n_steps);


    free(input);
    free(current);
    free(next);
}

