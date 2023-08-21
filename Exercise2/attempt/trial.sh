#!/bin/bash
gcc trial.c -o trial.x
rm output.csv # Remove previous file.
for i in 1000 2000 3000
do
    echo $i
    ./trial.x $i > output.txt

    # Extract information using grep and regular expressions
    size=$(grep -o 'Size: [0-9]*' output.txt | cut -d' ' -f2)
    time=$(grep -o 'Time: [0-9]*' output.txt | cut -d' ' -f2)
    gflops=$(grep -o 'GFLOPS: [0-9]*' output.txt | cut -d' ' -f2)

    # Store the extracted information in a CSV file
    #!/bin/bash
    if [ ! -e "output.csv" ]; then
        echo "Size,Time,GFLOPS" > output.csv
    fi
    echo "$size,$time,$gflops" >> output.csv

    echo "Extracted information has been stored in output.csv"
done