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
# Setup genetic maps: Proxyization model uses the genetic maps, we initialize them first.
mkdir genetic_maps
./ProxyTyper.sh -setup_BEAGLE_genetic_maps https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

# 1) Build the tag hashing proxy model:
./ProxyTyper.sh -generate_tag_proxizer_model TYPED_VARIANTS.bed 
if [[ $? -ne 0 ]]
then
    echo "Tag proxization model generation failed."
    exit 1
fi

# this model is also stored under ${MODELS_DIR} (defined in PROXYTYPER.ini file). Note that the hashing model file is a binary file and cannot be viewed with a text editor.
# The other parameters related to vicinty proxization in PROXYTYPER.ini:
# filter_proxization_vicinity_size=5        : This is the size of left/right vicinity for the variants that are used to calculate the hash. Note that a subset of these variants are used.
# coding_modulus=2                          : This is always 2.
# allele_err_eps=0.00                       : This is an independent error term but is always set to 0.
# normalized_N_e=600                        : This parameter describes how much we should take genetic distance while choosing the surrounding variants in hash model. Higher values of this value forces model to select closer variants to calculate the hash for the center variant.

# var_weight_prob=0.3                       : First order parameter selection probability; i.e., w_i
# var2var_interaction_prob=0.8              : Second order interaction parameter selection probability, i.e., w_ij
# var2var2var_interaction_prob=0.6          : Third order interaction parameter selection probability, i.e., w_ijk
# filter_weight_inversion_prob=0.5          : The probability of inverting the 
# filter_proxization_min_n_params_per_var=2 : The minimum number of surrounding variants used to calculate the hash for a typed variant. This weight selection is repeated until ProyTyper chooses at least this many parameters.
# Remarks:
# Most of these parameters do not make significant difference in results but may impact privacy considerations. Following parameters are important to ensure model is effective in mixing alleles:
# 1) normalized_N_e>1500 may constrain the model to select very close variants only, which causes only an inversion of the alleles without any mixing with surrounding variants. 
# 2) filter_weight_inversion_prob=0.5 is an important parameter to include the bias terms with 50% probability. This value should not be changed in any case.
# 3) filter_proxization_vicinity_size=0 is a special case that forces the model to use only the center variant with no mixing with surrounding variants. This is used in unphased genotype imputation. In this case, filter_proxization_min_n_params_per_var=1 must be set. In this case, users must rely on other mechanisms such as local permutations of typed variants.
##################################################################################################
# Now we can use the model to calculate the hashed alleles for each tag variant:
./ProxyTyper.sh -proxize_tag_variants ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_proxy.hashing
if [[ $? -ne 0 ]]
then
    echo "Tag proxy allele calculation failed in typed variants panel."
    exit 1
fi

# Now we update the current proxized panel to the proxy tag panel.
TYPED_VARIANTS_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_proxy.hashing

# We can also prozie the typed variants in the original panel, which contains the typed+untyped variants. In this case, the untyped variants are copied as they are:
./ProxyTyper.sh -proxize_tag_variants ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_proxy.hashing
if [[ $? -ne 0 ]]
then
    echo "Tag proxy allele calculation failed in typed+untyped variants panel."
    exit 1
fi

TYPED_UNTYPED_VARIANTS_panel_ID=${TYPED_UNTYPED_VARIANTS_panel_ID}_proxy.hashing

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