#!/bin/bash

# Set a maximum number of iterations (to prevent infinite loops)
max_retries=5000000
retry_count=0

# Specify the output file to check
output_file="Output/Branch/H1/mlradBH1_0.4.out"

while [ $retry_count -lt $max_retries ]; do
    # Run the first command
    reEvoRadical ctl/evoRadical.dat

    # Check if the output file exists
    if [ ! -f "$output_file" ]; then
        echo "$output_file not found, retrying..."
        retry_count=$((retry_count + 1))
    else
        echo "$output_file found. Execution completed."
        break
    fi
done

# If the maximum number of repetitions is reached
if [ $retry_count -eq $max_retries ]; then
    echo "Reached maximum retries. The output file $output_file was not found after $max_retries attempts."
fi
