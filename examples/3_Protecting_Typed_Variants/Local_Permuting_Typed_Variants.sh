#!/bin/bash

# This file implements an example of mechanisms that are used for protecting typed variants.

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh}\", you can it copy it from under scripts/ directory.."
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

##################################################################################################
# Import.
./ProxyTyper.sh -import_VCF ${ref_VCF} ${is_panel_phased} ${update_variant_ID} ${all_REF_panel_ID}
if [[ $? -ne 0 ]]
then
    echo "VCF importing failed"
    exit 1
fi

##################################################################################################
# Take the unique variants.
./ProxyTyper.sh -uniquefy_panel_variants ${all_REF_panel_ID} ${all_REF_panel_ID}_unique
if [[ $? -ne 0 ]]
then
    echo "Unique panel extraction failed"
    exit 1
fi

all_REF_panel_ID=${all_REF_panel_ID}_unique

# Store the panel with all variants in a variable since we will use this later on.
TYPED_UNTYPED_VARIANTS_panel_ID=${all_REF_panel_ID}

##################################################################################################
# Subselect 10,000 variants as tag variants; in this case, we also subselect the variants with vertain AAF cutoff:
awk {'split($4, arr, "_");if(arr[4]>0.05 && arr[4]<0.95){print $0}'} ${all_REF_panel_ID}_variants.bed > ALL_VARIANTS.bed
n_all_variants=`wc -l ALL_VARIANTS.bed | awk {'print $1'}`
n_typed=10000
shuf ALL_VARIANTS.bed | head -n ${n_typed} | sort -n -k2,2 > TYPED_VARIANTS.bed

# Get the genotype panel for the typed variants.
./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} TYPED_VARIANTS.bed ${all_REF_panel_ID}_typed
if [[ $? -ne 0 ]]
then
    echo "Variant subsetting failed."
    exit 1
fi

TYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}_typed

##################################################################################################
# 1) Build the permutation-based proxy model:
./ProxyTyper.sh -generate_tag_permuter_model TYPED_VARIANTS.bed 
if [[ $? -ne 0 ]]
then
    echo "Tag proxization model generation failed."
    exit 1
fi

# this model is also stored under ${MODELS_DIR} (defined in PROXYTYPER.ini file). Permutations are stored in an extended BED file (PROXIZATION_MODELS/permute_proxization_regions.bed) that define a mapping between the regions in original and permuted panels.
# The other parameters related to permuting typed variants in PROXYTYPER.ini:
# perm_proxization_vicinity_size=2    : The vicinity around each variant where we perform permutation of variants (while keeping genotypes constant)
# perm_proxy_geno_inversion_prob=0.5  : Probability of inverting the alleles (i.e., bias probability in permutation)
# perm_proxy_perm_prob=0.1            : Probability of performing permutation at each locus.
# Remarks:
# Permutation is a very strong protection technique but it may adversely impact accuracy. When other mechanisms are use, it is useful to set perm_proxy_perm_prob low (e.g., 0.1) is beneficial.
##################################################################################################
# Now we can use the model to calculate the hashed alleles for each tag variant:
./ProxyTyper.sh -permute_proxize_tag_variants ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_proxy.permute
if [[ $? -ne 0 ]]
then
    echo "Tag permutations failed in typed variants panel."
    exit 1
fi

# Now we update the current proxized panel to the proxy tag panel.
TYPED_VARIANTS_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_proxy.permute

# We can also proxize the typed variants in the original panel, which contains the typed+untyped variants. In this case, the untyped variants are copied as they are:
./ProxyTyper.sh -permute_proxize_tag_variants ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_proxy.permute
if [[ $? -ne 0 ]]
then
    echo "Tag permutations failed in typed+untyped variants panel."
    exit 1
fi

TYPED_UNTYPED_VARIANTS_panel_ID=${TYPED_UNTYPED_VARIANTS_panel_ID}_proxy.permute

# Calculate the alternate allele frequencies in the variants. Note that these don't look like the original values any more.
./ProxyTyper.sh -calculate_panel_AAF ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

./ProxyTyper.sh -calculate_panel_AAF ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

# Report the number of times we permuted the typed variants.
awk '{if($7!=$8){n_permuted+=1}}END{print n_permuted" variants among "NR" variants are permuted in positions."}' PROXIZATION_MODELS/permute_proxization_regions.bed