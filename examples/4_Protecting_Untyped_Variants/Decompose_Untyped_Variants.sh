#!/bin/bash

# This file implements an example of variant decomposition mechanism.

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

TYPED_UNTYPED_VARIANTS_panel_ID=${all_REF_panel_ID}

# Store the original panel id with all variants in a variable since we will use this later on.
ORIGINAL_TYPED_UNTYPED_VARIANTS_panel_ID=${all_REF_panel_ID}

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
# Decomposition is performed with the following command. This command takes the typed variants, and reference panel as input.
# Note that decomposition is only performed for protecting the untyped variants at the reference site (query does not have untyped variants to begin with.)
./ProxyTyper.sh -decompose_variants TYPED_VARIANTS.bed ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_decomposed
if [[ $? -ne 0 ]]
then
    echo "Decomposition failed"
    exit 1
fi

# Update the typed+untyped reference panel id.
DECOMPOSING_INTERVAL_MODEL_FILE=${TYPED_UNTYPED_VARIANTS_panel_ID}.interval
TYPED_UNTYPED_VARIANTS_panel_ID=${TYPED_UNTYPED_VARIANTS_panel_ID}_decomposed
# There are 3 parameters that are used to tune untyped variant decomposition.
# variant_decomp_min_AAF=0.1          : The minimum alternate allele frequency at which ProxyTyper performs partitioning. Decreasing this too much will increase decomposed panel size a lot.
# shuffle_decomp_vars=1               : This flag indicates ProxyTyper to use a 50% random bias on each decomposed variant's alleles (independently)
# USE_MT_RECOMPOSITION=1              : Asks ProxyTyper to use multithreaded option for decomposition/recomposition.
# Remarks:
# Unlike previous mechanisms that protect typed variants (which did not affect untyped variants), decomposition has no affect on the typed variants. The typed variants are returned as they are.
# Note that unlike previous mechanisms that we saw, this option generates the decomposition model (i.e., the interval file) and decomposes the reference panel simultaneously. This is by design since this mechanism is run only at the reference site. The decomposition model file should not be immediately shared with the query site.
# As with previous mechanisms, decomposition is randomly performed, i.e., if you run it on the same panel twice, the results will be different for the decomposed variants.
# The model for variant decomposition is an interval file that described the partitioning of untyped variants. This file is not stored with other model files because it is specific to the reference site.
##################################################################################################
# Calculate the alternate allele frequencies in the variants. Note that these don't look like the original values any more.
./ProxyTyper.sh -calculate_panel_AAF ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

./ProxyTyper.sh -calculate_panel_AAF ${ORIGINAL_TYPED_UNTYPED_VARIANTS_panel_ID} ${ORIGINAL_TYPED_UNTYPED_VARIANTS_panel_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

N_VARIANTS_IN_DECOMP_PANEL=`wc -l ${TYPED_UNTYPED_VARIANTS_panel_ID}_variants.bed | awk {'print $1'}`
N_VARIANTS_IN_ORIGINAL_PANEL=`wc -l ${ORIGINAL_TYPED_UNTYPED_VARIANTS_panel_ID}_variants.bed | awk {'print $1'}`

# Report the total number of variants.
echo "Decomposed panel contains ${N_VARIANTS_IN_DECOMP_PANEL} variants in augmented typed variants panel)."
echo "Proginal panel contains ${N_VARIANTS_IN_ORIGINAL_PANEL} variants)."

echo "The decomposition model interval file can be viewed in ${DECOMPOSING_INTERVAL_MODEL_FILE}"

# We can continue decomposing the decomposed panel to further expand the MAF spectrum into lower frequency but this would increase the number of variants quadratically.