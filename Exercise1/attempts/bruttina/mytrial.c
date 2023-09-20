#include <stdio.h>

int main(int argc, char** argv){

    int i=1;
    int t = i%2;
    int k = 5;
    printf("%d\n",t);
    for(int j = t; j<k;j+=2){
        printf("%d\n",j);
    }

}