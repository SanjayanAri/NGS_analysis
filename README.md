# NGS_analysis
Shell scripts for NGS analysis

Contains scripts for trimming, genome assembly, somatic and germline variant calling, and bacterial wgs processing.

## Bacterial wgs processing 
Contains:
Trimming through fastp 
Variant calling through snippy
Genome assembly through Unicycler
mob_recon for plasmid 
prokka for annotation

## Human genome processing (somatic variants)
Note: Since no tumor WES is available the script is just using wes from the patient (normal,germline)
which can be possible as research suggests that we can obtain somatic variants form pheriphral blood samples.

using gatk workflow

## Germline
germline variant calling

## New addition of heatmaps scripts and HLA typing workflow. 
python scripts for heatmaps and HLA typing is done manual workflow using Kourami tool.

