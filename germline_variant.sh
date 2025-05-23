#!/bin/bash

# File path: /home/parasmals/Scripts/germline_variant.sh

echo " Thank you for running this script!"
echo
echo " Script made by Sanjayan. A"
echo
echo "make sure you have enough space in your disk to store the output files."
echo
echo "only trimmed reads are taken as input"
echo
# set variables
ref="/home/parasmals/Swathi/Reference/hg38.fa"
funcotator_data_sources="/mnt/NGS_work/Funcotator/germline/funcotator_dataSources.v1.8.hg38.20230908g/gnomAD_exome/hg38"
known_sites="/mnt/NGS_work/BQSR/Homo_sapiens_assembly38.dbsnp138.vcf"

echo "Enter the sequence name for prefix for the run: "
read prefix

# Create base directory
base_dir="/mnt/NGS_work/Germline${prefix}"
mkdir -p "$base_dir"

# defining directories
align_dir="${base_dir}/align"
germline_dir="${base_dir}/germline"
filt_and_annotated_dir="${base_dir}/filtered_and_annotated"
final_result_dir="/home/parasmals/processing/human/germline/${prefix}"

# Create directory
mkdir -p "$align_dir" "$germline_dir" "$filt_and_annotated_dir" "$final_result_dir"

#input 
echo "Enter the trim read 1"
read read1
echo "Enter the trim read 2"
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

# Step 3: BQSR
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
echo
echo
echo "Will you halt this eclipse in me? ~ Sleeptoken"
echo 
# Step 4: Germline variant calling
echo "Running HaplotypeCaller"
gatk HaplotypeCaller -R "$ref" -I "$align_dir/${prefix}_BQSR.bam" -O "$germline_dir/${prefix}.g.vcf.gz" 

if [[ $? -ne 0 ]]; then
    echo "Error: HaplotypeCaller failed. Exiting script."
    exit 1
fi
echo "HaplotypeCaller completed successfully. Output: $germline_dir/${prefix}.g.vcf.gz"
echo

# filtering and annotation
echo "spliting, filtering and annotating the germline vcf file"

gatk SelectVariants -R "$ref" -V "$germline_dir/${prefix}.g.vcf.gz" -select-type SNP -O "$germline_dir/${prefix}_snps.vcf.gz"
gatk SelectVariants -R "$ref" -V "$germline_dir/${prefix}.g.vcf.gz" -select-type INDEL -O "$germline_dir/${prefix}_indels.vcf.gz"

if [[ $? -ne 0 ]]; then
    echo "Error: SelectVariants failed. Exiting script."
    exit 1
fi
echo "SelectVariants completed successfully. Output: $germline_dir/${prefix}_snps.vcf.gz and $germline_dir/${prefix}_indels.vcf.gz" 
echo
echo "Filtering SNPs"
gatk VariantFiltration -R "$ref" -V "$germline_dir/${prefix}_snps.vcf.gz" --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 20" -genotype-filter-name "GQ_filter" -O "$germline_dir/${prefix}_snps_filtered.vcf.gz"
gatk VariantFiltration -R "$ref" -V "$germline_dir/${prefix}_indels.vcf.gz" --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0" --filter-name "InbreedingCoeff_filter" -filter "InbreedingCoeff < -0.8" -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" -genotype-filter-expression "GQ < 20" -genotype-filter-name "GQ_filter" -O "$germline_dir/${prefix}_indels_filtered.vcf.gz"

if [[ $? -ne 0 ]]; then
    echo "Error: VariantFiltration failed. Exiting script."
    exit 1
fi  


echo "Filtering completed successfully. Output: $germline_dir/${prefix}_snps_filtered.vcf.gz and $germline_dir/${prefix}_indels_filtered.vcf.gz"
echo
echo "Preping for annotation"
echo
gatk SelectVariants -R "$ref" -V "$germline_dir/${prefix}_snps_filtered.vcf.gz" -O "$germline_dir/${prefix}_snps_ready.vcf.gz" --exclude-filtered
gatk SelectVariants -R "$ref" -V "$germline_dir/${prefix}_indels_filtered.vcf.gz" -O "$germline_dir/${prefix}_indels_ready.vcf.gz" --exclude-filtered

echo "Annotating SNPs"
gatk Funcotator -R "$ref" -V "$germline_dir/${prefix}_snps_ready.vcf.gz" --output "$filt_and_annotated_dir/${prefix}_snps_annotated.vcf.gz" --data-sources-path "$funcotator_data_sources" --output-file-format VCF --ref-version hg38

echo "Annotating Indels"
gatk Funcotator -R "$ref" -V "$germline_dir/${prefix}_indels_ready.vcf.gz" --output "$filt_and_annotated_dir/${prefix}_indels_annotated.vcf.gz" --data-sources-path "$funcotator_data_sources" --output-file-format VCF --ref-version hg38
echo "Annotation completed successfully. Output: $filt_and_annotated_dir/${prefix}_snps_annotated.vcf.gz and $filt_and_annotated_dir/${prefix}_indels_annotated.vcf.gz"
echo
echo "Copying the final results to $final_result_dir"
cp "$filt_and_annotated_dir/${prefix}_snps_annotated.vcf.gz" "$final_result_dir"
cp "$filt_and_annotated_dir/${prefix}_indels_annotated.vcf.gz" "$final_result_dir"  
cp "$filt_and_annotated_dir/${prefix}_snps_filtered.vcf.gz" "$final_result_dir"
cp "$filt_and_annotated_dir/${prefix}_indels_filtered.vcf.gz" "$final_result_dir"
echo 
if 
[[ $? -ne 0 ]]; then
    echo "Error: Copying files failed. Exiting script."
    exit 1
fi
echo "Copying completed successfully. Output: $final_result_dir"

echo "The final filtered and annotated germline variant files are stored in: $final_result_dir"
echo

echo "All processes completed successfully."
echo
echo "Thank you for using this script!"
echo

exit 0

