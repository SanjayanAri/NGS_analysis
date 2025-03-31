#!/bin/bash
echo "Genome Assembly Script"
echo "Run this script in conda enviorment with the following packages installed: UNICYCLER, QUAST, BANDAGE"
# Input prompt
prompt_for_input () {
  echo "Input read file 1: "
  read read1
  if [ ! -f $read1 ]; then
    echo "File not found!"
    prompt_for_input
  fi

  echo "Input read file 2: "
  read read2
  if [ ! -f $read2 ]; then
    echo "File not found!"
    prompt_for_input
  fi

  echo "Enter where do you want the output files to be saved: "
  read output_dir 
}

# Main Script
while true; do
  # Prompt for input
    prompt_for_input

  # Run Unicycler
  unicycler -1 "$read1" -2 "$read2" -o "$output_dir" -t 16

  # Check if the assembly is successful
  if [[ $? -eq 0 ]]; then
    echo "Assembly successful!"
    break
  else
    echo "Assembly failed!"
    echo "Do you want to try again? (y/n)"
    read response
    if [[ $response == "n" ]]; then
      break
    fi
  fi
done
