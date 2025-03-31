#!/bin/bash


# File path: /home/parasmals/Scripts/human_genome_processing.sh
echo " Thank you for running this script!"
echo
echo " Ensure you are running this in the environment where the following tools are installed: BWA (for alignment), Samtools (for sorting and indexing), and GATK (for variant calling)."
echo
echo " This script is for processing human genome data, version hg38 (GRCh38) "
echo
echo " This is a huge long script that will take long time to run. Please be patient."
echo 
echo " Ensure you have enough space in your disk to store the output files."
echo 
echo " Make sure you have mounted the disk where before you run this script, or this will fail."
echo
echo " If you have not mounted the disk use the following command to mount the disk: " 
echo
echo " sudo mount -t ntfs-3g /dev/sdb2 /mnt/  copy paste this command in fish terminal to mount the disk and enter the password."
echo 
echo " If all is ready press Y to continue "

# prompt user to continue or exit
read -p "Continue? (Y/N): " response
if [[ "$response" == "Y" || "$response" == "y" ]]; then
    echo "Great! Let's get started."
else
    echo "Exiting script, good bye."
    exit 1
fi
# Input prompt

echo "Enter the sequence name for prefix for the run: "
read prefix

# Create base directory
base_dir="/mnt/NGS_work/vc_somatic/${prefix}"
mkdir -p "$base_dir"

# Define directories
align_dir="${base_dir}/align"
somatic_dir="${base_dir}/somatic"
filt_and_annotated_dir="${base_dir}/filtered_and_annotated"
final_result_dir="/home/parasmals/processing/human/${prefix}"

# Create directories
mkdir -p "$align_dir" "$somatic_dir" "$filt_and_annotated_dir" "$final_result_dir"

# Define paths
ref="/home/parasmals/Swathi/Reference/hg38.fa"
known_sites="/mnt/NGS_work/BQSR/Homo_sapiens_assembly38.dbsnp138.vcf"
germline_resource="/mnt/NGS_work/somatic_db/af-only-gnomad.hg38.vcf.gz"
pon="/mnt/NGS_work/somatic_db/1000g_pon.hg38.vcf.gz"
exome_calling_regions="/mnt/NGS_work/somatic_db/exome_calling_regions.v1.1.interval_list"
funcotator_data_sources="/mnt/NGS_work/Funcotator/funcotator_dataSources.v1.8.hg38.20230908s"

# Input for reads
echo "Enter the trimmed read  1"
read read1
echo "Enter the trimmed read  2"
read read2

# Step 1: Alignment with BWA
echo " Running BWA for alignment "
echo
bwa mem -t 16 -R "@RG\tID:${prefix}\tPL:ILLUMINA\tSM:${prefix}" "$ref" "$read1" "$read2" > "$align_dir/${prefix}.sam"

if [[ $? -ne 0 ]]; then
    echo "Error: BWA failed exiting script"
    exit 1
fi
echo "Alignment completed successfully. output: $align_dir/${prefix}.sam"
echo 
echo "Running MarkduplicatesSpark"

# Step 2: Mark dups

gatk MarkDuplicatesSpark -I "$align_dir/${prefix}.sam" -O "$align_dir/${prefix}_dedup.bam"

if [[ $? -ne 0 ]]; then
    echo "Error: MarkDuplicatesSpark failed. Exiting script."
    exit 1
fi
echo "MarkDuplicatesSpark completed successfully. Output: $align_dir/${prefix}_dedup.bam"

# Removing SAM file
rm "$align_dir/${prefix}.sam"

# Step #: BQSR
echo "Running BaseRecalibrator"
gatk BaseRecalibrator -R "$ref" -I "$align_dir/${prefix}_dedup.bam" --known-sites "$known_sites" -O "$align_dir/${prefix}_recal_data.table"

if [[ $? -ne 0 ]]; then
    echo "Error: BaseRecalibrator failed. Exiting script."
    exit 1
fi
echo "BaseRecalibrator completed successfully. Output: $align_dir/${prefix}_recal_data.table"
echo
echo "Running ApplyBQSR"
gatk ApplyBQSR -R "$ref" -I "$align_dir/${prefix}_dedup.bam" --bqsr-recal-file "$align_dir/${prefix}_recal_data.table" -O "$align_dir/${prefix}_BQSR.bam"

if [[ $? -ne 0 ]]; then
    echo "Error: ApplyBQSR failed. Exiting script."
    exit 1
fi
echo "ApplyBQSR completed successfully. Output: $align_dir/${prefix}_BQSR.bam"

# Step 4: Variant calling
echo " Running GATK Mutect2 for variant calling"
echo " Note since there is no tumor WES data, This will not be accurate"
echo
gatk Mutect2 -R "$ref" -I "$align_dir/${prefix}_BQSR.bam" --germline-resource "$germline_resource" --panel-of-normals "$pon" -O "$somatic_dir/${prefix}_som_raw.vcf" --f1r2-tar-gz "$somatic_dir/${prefix}_f1r2.tar.gz"

if [[ $? -ne 0 ]]; then
    echo "Error: Mutect2 failed. Exiting script."
    exit 1
fi
echo "Mutect2 completed successfully. Output: $somatic_dir/${prefix}_som_raw.vcf"
echo
echo "Filtering the variants"
echo
# Step 5: Filter variants
echo "getting pile up summary"
gatk GetPileupSummaries -I "$align_dir/${prefix}_BQSR.bam" -V "$germline_resource" -L "$exome_calling_regions" -O "$somatic_dir/${prefix}_pileups.table"

if [[ $? -ne 0 ]]; then
    echo "Error: GetPileupSummaries failed. Exiting script."
    exit 1
fi
echo "GetPileupSummaries completed successfully. Output: $somatic_dir/${prefix}_pileups.table"
echo
echo "Running CalculateContamination"
gatk CalculateContamination -I "$somatic_dir/${prefix}_pileups.table" -O "$somatic_dir/${prefix}_contamination.table"

if [[ $? -ne 0 ]]; then
    echo "Error: CalculateContamination failed. Exiting script."
    exit 1
fi
echo "CalculateContamination completed successfully. Output: $somatic_dir/${prefix}_contamination.table"
echo
echo " Running ReadOrientationModel"
gatk LearnReadOrientationModel -I "$somatic_dir/${prefix}_f1r2.tar.gz" -O "$somatic_dir/${prefix}_read_orientation_model.tar.gz"

if [[ $? -ne 0 ]]; then
    echo "Error: LearnReadOrientationModel failed. Exiting script."
    exit 1
fi
echo "LearnReadOrientationModel completed successfully. Output: $somatic_dir/${prefix}_read_orientation_model.tar.gz"
echo
echo "Running FilterMutectCalls"

gatk FilterMutectCalls -V "$somatic_dir/${prefix}_som_raw.vcf" -R "$ref" --contamination-table "$somatic_dir/${prefix}_contamination.table" --ob-priors "$somatic_dir/${prefix}_read_orientation_model.tar.gz" -O "$filt_and_annotated_dir/${prefix}_filt_som.vcf"

if [[ $? -ne 0 ]]; then
    echo "Error: FilterMutectCalls failed. Exiting script."
    exit 1
fi  
echo "FilterMutectCalls completed successfully. Output: $filt_and_annotated_dir/${prefix}_filt_som.vcf"
echo
echo "Running Funcotator"
gatk Funcotator -V "$filt_and_annotated_dir/${prefix}_filt_som.vcf" -R "$ref" --ref-version hg38 --data-sources-path "$funcotator_data_sources" -O "$filt_and_annotated_dir/${prefix}_annotated_som.vcf" --output-file-format VCF

if [[ $? -ne 0 ]]; then
    echo "Error: Funcotator failed. Exiting script."
    exit 1
fi
echo "Funcotator completed successfully. Output: $filt_and_annotated_dir/${prefix}_annotated_som.vcf"
echo
echo "All processes completed successfully."
echo
cp "$filt_and_annotated_dir/${prefix}_annotated_som.vcf" "$final_result_dir/"
cp "$filt_and_annotated_dir/${prefix}_filt_som.vcf" "$final_result_dir/"
echo
echo "The final filtered somatic variant file is stored in: $final_result_dir/${prefix}_filt_som.vcf"
echo
echo "The final annotated somatic variant file is stored in: $final_result_dir/${prefix}_annotated_som.vcf"
echo
echo "Maybe you've been down too long ~ Slipknot"
echo
echo "Thank you for using this script. Good bye!"


exit 0

