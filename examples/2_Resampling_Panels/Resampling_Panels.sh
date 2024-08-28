#!/bin/bash

# This file implements an example of panel resampling.

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh\", you can it copy it from under scripts/ directory.."
    exit 1
fi

##################################################################################################
# Import the VCF file and save it to a genotype panel matrix named "KG_chr20". We already saw this in previous example.
# Download 1000 Genomes Project's chromosome 20.
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# We should have the vcf file here:
if [[ ! -f "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ]]
then
    echo "Failed to download the VCF file: ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    exit 1
fi

ref_VCF="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

all_REF_panel_ID=KG_chr20
is_panel_phased=1
update_variant_ID=1

# Clean all the files for the panel.
rm -f ${all_REF_panel_ID}_*
rm -f -r genetic_maps
rm -f -r plink.*

./ProxyTyper.sh -import_VCF ${ref_VCF} ${is_panel_phased} ${update_variant_ID} ${all_REF_panel_ID}
if [[ $? -ne 0 ]]
then
    echo "VCF importing failed"
    exit 1
fi

##################################################################################################

./ProxyTyper.sh -uniquefy_panel_variants ${all_REF_panel_ID} ${all_REF_panel_ID}_unique
if [[ $? -ne 0 ]]
then
    echo "Unique panel extraction failed"
    exit 1
fi

all_REF_panel_ID=${all_REF_panel_ID}_unique

##################################################################################################
# We first download and setup the genetic maps:
mkdir genetic_maps
./ProxyTyper.sh -setup_BEAGLE_genetic_maps https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
if [[ $? -ne 0 ]]
then
    echo "VCF importing failed"
    exit 1
fi

# ProxyTyper uses the 4 column plink-formatted genetic maps. Above command should download and copy the 
# plink formatted genetic maps to the directory specific in PROXYTYPER.ini as following:
# PROXYTYPER_GENETIC_MAPS_DIR=genetic_maps

##################################################################################################

# We will now resample the panel up to 5,000 subjects. The resampling configuration is defined by following lines 
# in PROXYTYPER.ini:
#
# haplocoded_flag=1
# pre_mosaic_Ne=0.125 # This is the effective population size per # of reference haplotypes. Overall, this corresponds to 12.5 times the original reference size.
# pre_mosaic_allele_eps=0.0 # We don't add any random errors. If user wishes to add errors, this parameter can be set.
# pre_mosaic_n_threads=${n_threads} # Number of threads that resampling uses.
# min_cM_delta_per_anchoring_vars=0.001 # This is the minimum cM distance between anchor points used for selecting haplotype states.
# max_l_seg_n_bps=10000000 # Maximum segment length in base pairs to ensure that we don't sample very long haplotypes by chance. This is a hard cutoff.
# max_l_seg_cM=0 # Maximum segment length in cMs, 0 means ignore.
# max_l_seg_nvars=0 # Maximum segment length in number of variants, 0 means ignore.
#

# To re-sampling the panel, we run re-sample option of ProxyTyper:
# Below we resample the panel to 5,000 subjects.
n_resample_size=5000
./ProxyTyper.sh -resample_panel ${all_REF_panel_ID} ${n_resample_size} ${all_REF_panel_ID}_resample
if [[ $? -ne 0 ]]
then
    echo "Resampling failed"
    exit 1
fi

# Update the panel id so we can use the resampled panel identifier.
all_REF_panel_ID=${all_REF_panel_ID}_resample

##################################################################################################

# We now calculate the alternate allele frequencies in the resampled panel.
# Note that the resampled is expected to have variations in allele frequencies due to random haplotype selection.
./ProxyTyper.sh -calculate_panel_AAF ${all_REF_panel_ID} ${all_REF_panel_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

################################

# We can save the resampled panel as a VCF formatted file, if it is to be used as is:
./ProxyTyper.sh -export_VCF ${all_REF_panel_ID} hg19 1 ${all_REF_panel_ID}.vcf.gz
if [[ $? -ne 0 ]]
then
    echo "Variant subset failed."
    exit 1
fi

echo "The exported VCF file can be found @ ${all_REF_panel_ID}.vcf.gz"..
echo "Completed.."

