#!/bin/bash

#Function to ask user for input
prompt_for_input() {
    echo "Enter the first read file: "
    read read1
    if [[ ! -f $read1 ]]; then
        echo "File not found. Please try again."
        prompt_for_input
    fi

    echo "Enter the second read file: "
    read read2
    if [[ ! -f $read2 ]]; then
        echo "File not found. Please try again."
        prompt_for_input
    fi

    echo "Enter the desired output file name(Fw read): "
    read Trim_read1

    echo "Enter the desired output file name(Rv read): "
    read Trim_read2

    echo "Enter the html output file name: "
    read Trim_reporthtml

}

#Script
while true; do
    #Prompt user for input
    prompt_for_input

    #Run fastp command
    fastp --detect_adapter_for_pe --overrepresentation_analysis --cut_right --html "$Trim_reporthtml" --thread 8 -i "$read1" -I "$read2" -o "$Trim_read1" -O "$Trim_read2"

    #check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Trimming successful. Files saved as $Trim_read1 and $Trim_read2 and report saved as $Trim_reporthtml"
        break
    else
        echo "Trimming failed. Please try again."
    fi
done