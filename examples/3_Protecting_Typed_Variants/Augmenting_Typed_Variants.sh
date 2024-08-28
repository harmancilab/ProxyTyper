#!/bin/bash

# This file implements an example of mechanisms that are used for protecting typed variants.

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
# Setup genetic maps.
mkdir genetic_maps
#./ProxyTyper.sh -setup_BEAGLE_genetic_maps https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

##################################################################################################
# Subselect 10,000 variants as tag variants.
cp ${all_REF_panel_ID}_variants.bed ALL_VARIANTS.bed
n_all_variants=`wc -l ALL_VARIANTS.bed | awk {'print $1'}`
n_typed=10000
shuf ${all_REF_panel_ID}_variants.bed | head -n ${n_typed} | sort -n -k2,2 > TYPED_VARIANTS.bed

# Get the genotype panel for the typed variants.
./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} TYPED_VARIANTS.bed ${all_REF_panel_ID}_typed
if [[ $? -ne 0 ]]
then
    echo "Variant subsetting failed."
    exit 1
fi

TYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}_typed

#######
# 1) Augment the variants: We need to first generate the augmentation "model", which is basically a mapping between the typed regions and their augmented counterparts:
./ProxyTyper.sh -generate_tag_augmenter_model TYPED_VARIANTS.bed ALL_VARIANTS.bed
if [[ $? -ne 0 ]]
then
    echo "Augment model generation failed."
    exit 1
fi

# this model is stored under ${MODELS_DIR} (defined in PROXYTYPER.ini file).
# The parameters related to augmentation mechanism in PROXYTYPER.ini are as following:
# n_tag_augments=3                  : The number of times we will augment typed variants, this is cumulative, meaning the each augmentation starts from previous step's augmented typed variants.
# tag_augmenter_probability=0.99    : Probability of augmenting any typed variant
# n_tags_per_augment_vicinity=0     : The number of tags in the vicinity to augment each typed variant. For typed variant i; a new proxy typed variant is augmented within [i-n_vic, i+n_vic] block.
# tag_genotype_augmentation_type=1  : The type of genotype copy option. 0: Genotype is initialized as "0". 1: Genotype of augmented variants are copied from their "mates".
##################################################################################################
# Now we can augment the tag variants:
./ProxyTyper.sh -augment_tag_variants ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_tag_augmented
if [[ $? -ne 0 ]]
then
    echo "Tag augmentation failed in typed variants panel."
    exit 1
fi

# Now we update the current proxized panel to the tag augmented panel.
TYPED_VARIANTS_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_tag_augmented

# We now have a panel with tag variants augmented at random positions in the vicinity of the original variants.
# With the current setup, we get observe around 70,000 new variants.
N_AUGMENTED_TAGS_IN_TYPED_PANEL=`wc -l ${TYPED_VARIANTS_PANEL_ID}_variants.bed | awk {'print $1'}`

# We can also augment the typed variants in the original panel, which contains the typed+untyped variants. In this case, the untyped variants are copied as they are:
./ProxyTyper.sh -augment_tag_variants ${TYPED_UNTYPED_VARIANTS_panel_ID} ${TYPED_UNTYPED_VARIANTS_panel_ID}_tag_augmented
if [[ $? -ne 0 ]]
then
    echo "Tag augmentation failed in typed+untyped variants panel."
    exit 1
fi

TYPED_UNTYPED_VARIANTS_panel_ID=${TYPED_UNTYPED_VARIANTS_panel_ID}_tag_augmented

N_AUGMENTED_TAGS_IN_TYPED_UNTYPED_PANEL=`wc -l ${TYPED_UNTYPED_VARIANTS_panel_ID}_variants.bed | awk {'print $1'}`

# Report the total number of variants.
echo "Augmented to ${N_AUGMENTED_TAGS_IN_TYPED_PANEL} variants in augmented typed variants panel (Originally was ${n_typed})."
echo "Augmented to ${N_AUGMENTED_TAGS_IN_TYPED_UNTYPED_PANEL} variants in augmented typed+untyped variants panel (Originally was ${n_all_variants})."

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