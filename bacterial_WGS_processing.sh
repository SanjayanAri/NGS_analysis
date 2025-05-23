#!/bin/bash

#File path: /home/parasmals/Scripts/bacterial_WGS_processing.sh

echo
echo " You are the universe in ecstatic motion. ~ Rumi "
echo
echo "Script made by Sanjayan. A"
echo
echo "Contact me at: a.sanjayan2002@gmail.com"
echo
echo "This script is for processing bacterial WGS data."
echo
echo "This script does trimming of reads, genome assembly, annotation, plasmid detection and variant calling."
echo
echo " Run this script where the following tools are installed: Fastp (for trimming), Unicycler (for assembly), prokka (for annotation), Mob_suite (for plasmid detection) and Snippy (for variant calling)."
echo 
echo "I hope you are running this script in the conda environment 'auto' "
echo "if not then please activate the conda environment 'auto' before running this script."
echo "By ctrl+c to stop this script and then type ' conda activate auto ' and then run this script again."
echo
echo "In terminal ctrl+shift+v to paste the command"
echo "In terminal ctrl+shift+c to copy the command"
echo " Doing ctrl+c and ctrl+v in terminal will stop the script so know the difference "

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

    echo "Iput genus name: (This is used for prokka annotation first letter should be capitalized)"
    echo
    echo "For example: Escherichia, Salmonella, Klebsiella"
    read genus_name
    if [[ -z "$genus_name" ]]; then
        echo "Genus name cannot be empty! Please provide a valid genus name."
        prompt_for_input
    fi

    echo "Input species name: (This is used for prokka annotation)"
    read species_name
    if [[ -z "$species_name" ]]; then
        echo "Species name cannot be empty! Please provide a valid species name."
        prompt_for_input
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

# Run Prokka
echo "Running Prokka for annotation..."

prokka --outdir "$parent_dir/prokka_results" --prefix "$user_name" --cpus 16 --addgenes --addmrna --genus "$genus_name" --species "$species_name" --usegenus $parent_dir/genome_assembly/assembly.fasta
if [[ $? -eq 0 ]]; then
    echo "Annotation successful!"
else
    echo "Annotation failed! exiting script"
    exit 1
fi
echo "Annotation files are saved in $parent_dir/prokka_results"
echo

# Run Mob_recon
echo "Running Mob_recon for plasmid typing..."

mob_recon --infile "$parent_dir/genome_assembly/assembly.fasta" --outdir "$parent_dir/mob_recon_results" -u -n 8
if [[ $? -eq 0 ]]; then
    echo "Plasmid typing successful!"
else
    echo "Plasmid typing failed! exiting script"
    exit 1
fi
echo "Plasmid typing files are saved in $parent_dir/mob_recon_results"
echo


# End of script
echo "All processes completed successfully!"
echo
echo "Results are saved in $parent_dir"
echo
echo "Won't you come and dance in the dark with me? "
echo  
echo "~ assensionism, sleep token"
echo
echo
echo "Thank you for using this script!"

exit 0
