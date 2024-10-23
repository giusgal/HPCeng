#!/bin/bash
# Executes the file passed as argument for different input sizes
#  and saves the results in a csv file.

# Check if the output file is provided as an argument
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <executable> <csv_file>"
  exit 1
fi
executable=$1
output_file=$2

# Initialize CSV file with headers
echo "Input Size,Time (s)" > $output_file

for input_size in $((10**4)) $((10**5)) $((10**6)) $((10**7)) $((10**8))
do
  # Extract the time part using grep and awk
  time_taken=$(./$executable $input_size | grep "Time:" | awk '{print $2}')

  # Save to the CSV file
  echo $input_size,$time_taken >> $output_file

  echo "Input Size: $input_size, Time: $time_taken"
done
