// File usato solo per riordinare le idee sui puntatori in C, poi lo cancello
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main ( int argc, char **argv ){
    int* a;
    int b = 5;
    a = &b;

    *a = 7;
    printf("b: %p\n", &b);
    printf("a: %p\n", &(*a));
    printf("value of b: %d\n",b);
    printf("value of b: %d\n",*a);

    char* image;
    image="ciao";
    printf("%c\n", *image);
    printf("%p\n", image);
    printf("%p\n", &(*image));

    //DOUBLE POINTERS
    int var = 789;
 
    // pointer for var
    int* ptr1;
 
    // double pointer for ptr1
    int** ptr2;
 
    // storing address of var in ptr2
    ptr1 = &var;
 
    // Storing address of ptr2 in ptr1
    ptr2 = &ptr1;
 
    printf("Value of var = %d\n", var);
    printf("Value of var using single pointer = %d\n", *ptr1);
    printf("Value of var using double pointer = %d\n", **ptr2);
    printf("Address of var: %p\n", &var);
    printf("Address of var: %p\n", ptr1);
    printf("Address of var: %p\n", *ptr2);
    printf("Address of ptr1: %p\n", &ptr1);
    printf("Address of ptr1: %p\n", ptr2);


    return 0;
}