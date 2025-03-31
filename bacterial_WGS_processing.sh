#!/bin/bash

#File path: /home/parasmals/Scripts/bacterial_WGS_processing.sh
echo "Thank you for running this script!"
echo
echo " You are the universe in ecstatic motion. ~ Rumi "
echo
echo "Bacterial WGS Processing Script"
echo
echo "This script is for processing bacterial WGS data."
echo
echo "This script does trimming of reads, genome assembly, and variant calling."
echo
echo " Run this script where the following tools are installed: Fastp (for trimming), Unicycler (for assembly), and Snippy (for variant calling)."
echo
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

    echo "Input reference genome (in fasta format): "
    read ref_genome
    if [ ! -f $ref_genome ]; then
        echo "File not found!"
        prompt_for_input
    fi

    echo "Enter the prefix for VCF files: "
    read vcf_prefix

    echo "Enter the name of the result directory:(this used to create a new directory)"
    read user_name
    parent_dir="/home/parasmals/processing/bacterial/${user_name}_results"
    mkdir -p "$parent_dir"
    if [[ $? -ne 0 ]]; then
        echo "Error creating directory! check for permission"
        exit 1
    fi
}

echo "The result will be stored in the directory: $parent_dir"

# Main Script
prompt_for_input

# Run Fastp

echo "Running Fastp for quality trimming of reads..."
fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html "$parent_dir/fastp_report.html" --thread 16 -i "$read1" -I "$read2" -o "$parent_dir/trimmed_read1.fastq" -O "$parent_dir/trimmed_read2.fastq"

if [[ $? -eq 0 ]]; then
    echo "Trimming successful!"
else
    echo "Trimming failed! exiting script"
    exit 1
fi

echo "Trimmed reads are saved as $parent_dir/trimmed_read1.fastq and $parent_dir/trimmed_read2.fastq"
echo
# Run Snippy

echo "Running Snippy for variant calling..."

snippy --cpus 16 --outdir "$parent_dir/snippy_results" --ref "$ref_genome" --R1 "$parent_dir/trimmed_read1.fastq" --R2 "$parent_dir/trimmed_read2.fastq" --prefix "$vcf_prefix"

if [[ $? -eq 0 ]]; then
    echo "Variant calling successful!"
else
    echo "Variant calling failed! exiting script"
    exit 1
fi
echo
echo "VCF files are saved in $parent_dir/snippy_results"
echo


# Run Unicycler

echo "Running Unicycler for genome assembly..."

unicycler -1 "$parent_dir/trimmed_read1.fastq" -2 "$parent_dir/trimmed_read2.fastq" -o "$parent_dir/genome_assembly" -t 16

if [[ $? -eq 0 ]]; then
    echo "Assembly successful!"
else
    echo "Assembly failed! exiting script"
    exit 1
fi

echo "Assembly files are saved in $parent_dir/genome_assembly"
echo
echo "run the following command to visualize the assembly: bandage image $parent_dir/genome_assembly/assembly.gfa $parent_dir/genome_assembly/assembly.png"


# End of script
echo "All processes completed successfully!"
echo
echo "Results are saved in $parent_dir"
echo
echo "Won't you come and dance in the dark with me? "
echo "Diamonds in the trees, pentagrams in the night sky" 
echo "~ assensionism, sleep token"
echo

echo "now, ugh go away! bye bye"

exit 0