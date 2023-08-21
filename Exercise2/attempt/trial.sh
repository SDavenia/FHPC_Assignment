#!/bin/bash
gcc OB.c -o OB.x
gcc MKL.c -o MKL.x

rm -rf *.csv # Remove previous csv file.

for i in 1000 2000 3000
do

    # Run everything with openBLAS
    for j in 1 2 3 4 5 # Take multiple measurements
    do
        ./OB.x $i > output.txt #just a temporary file

        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt | cut -d' ' -f2)
        time=$(grep -o 'Time: [0-9]*' output.txt | cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9]*' output.txt | cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="OBLAS_$size.csv"

        if [ ! -e $filename ]; then
            echo "Size,Time,GFLOPS" > $filename
        fi
        echo "$size,$time,$gflops" >> $filename
    done

    # Now repeat for MKL library
    for j in 1 2 3 4 5
    do
        ./MKL.x $i > output.txt #just a temporary file

        # Extract information using grep and regular expressions
        size=$(grep -o 'Size: [0-9]*' output.txt | cut -d' ' -f2)
        time=$(grep -o 'Time: [0-9]*' output.txt | cut -d' ' -f2)
        gflops=$(grep -o 'GFLOPS: [0-9]*' output.txt | cut -d' ' -f2)
        # Store the extracted information in a CSV file
        filename="MKL_$size.csv"

        if [ ! -e $filename ]; then
            echo "Size,Time,GFLOPS" > $filename
        fi
        echo "$size,$time,$gflops" >> $filename
    done

done