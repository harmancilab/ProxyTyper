#!/bin/bash

# This file implements an example of coordinate/genetic map anonymization.

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

# Use panel with unique variants from here on.
all_REF_panel_ID=${all_REF_panel_ID}_unique

# Keep the panel identifier for the all variants separately since we may need it later.
TYPED_UNTYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}

# Store the original panel id with all variants in a variable since we will use this later on.
ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}

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
ORIGINAL_TYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}_typed

##################################################################################################
# Initializet the genetic maps, this is needed to anonymize the genetic maps.
mkdir genetic_maps
BEAGLE_MAPS_URL=https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
./ProxyTyper.sh -setup_BEAGLE_genetic_maps ${BEAGLE_MAPS_URL}
if [[ $? -ne 0 ]]
then
    echo "Genetic map setup failed, make sure directory exists."
    exit 1
fi

# Build the tag/map anonymization and save it. This step will take place at the reference site.
./ProxyTyper.sh -generate_coordinate_anonymizer_model ${TYPED_VARIANTS_PANEL_ID}_variants.bed ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed
if [[ $? -ne 0 ]]
then
    echo "Coordinate/Map anonymization failed."
    exit 1
fi

# List the files that we built
echo "The anoymization model files are written under PROXIZATION_MODELS/anon_genetic_maps directory: Below are the first 10 lines from each file"
head -n 10 PROXIZATION_MODELS/anon_genetic_maps/*

# Note that the coodinate anonymization is performed by mapping the original coordinates to a stretch of predefined chromosomal length. 
# The parameters used for building the anonymizer are defined under PROXYTYPER.ini file:
# Coordinate anonymization parameter:
# coord_anon_cM_noise_SD=0.05   : This is the standard deviation
# l_anonym_chrom=100000000      : This is the length of chromosome length that ProxyTyper maps all variant coordinates uniformly. 100 megabase is set by default and works for all chromosomes. Making this too large or too small may cause BEAGLE to fail or following mechanisms to fail to run.

# Now, anonymize the coordinates of query and reference panels:
# It important to make sure that the coordinate systems for the anonymization to match, otherwise ProxyTyper will give an error and exit.
# The basic requirement is that the panel that is being anonymized must be a strict subset of the coordinate anonymizing model's source coordinates.
./ProxyTyper.sh -anonymize_coordinates ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_anon_coords
if [[ $? -ne 0 ]]
then
    echo "Coordinate/Map anonymization failed."
    exit 1
fi
TYPED_VARIANTS_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_anon_coords

./ProxyTyper.sh -anonymize_coordinates ${TYPED_UNTYPED_VARIANTS_PANEL_ID} ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_anon_coords
if [[ $? -ne 0 ]]
then
    echo "Coordinate/Map anonymization failed."
    exit 1
fi
TYPED_UNTYPED_VARIANTS_PANEL_ID=${TYPED_UNTYPED_VARIANTS_PANEL_ID}_anon_coords

#########################################################################################################
# At this point, we have the coordiantes anonymized and concordant with the anonymized genetic map.
#########################################################################################################

# We can just de-anonymize these coordiantes using a similar option.
./ProxyTyper.sh -deanonymize_coordinates ${TYPED_VARIANTS_PANEL_ID} ${TYPED_VARIANTS_PANEL_ID}_deanonymized
if [[ $? -ne 0 ]]
then
    echo "Coordinate/Map anonymization failed."
    exit 1
fi
TYPED_VARIANTS_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_deanonymized

./ProxyTyper.sh -deanonymize_coordinates ${TYPED_UNTYPED_VARIANTS_PANEL_ID} ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_deanonymized
if [[ $? -ne 0 ]]
then
    echo "Coordinate/Map anonymization failed."
    exit 1
fi
TYPED_UNTYPED_VARIANTS_PANEL_ID=${TYPED_UNTYPED_VARIANTS_PANEL_ID}_deanonymized

# At this point, the variant coordinates should match exactly the original panels:
echo "Comparing original anonymized-then-deanonymized coordinates with the original coordinates:"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>(Following must be empty)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
diff ${TYPED_VARIANTS_PANEL_ID}_variants.bed ${ORIGINAL_TYPED_VARIANTS_PANEL_ID}_variants.bed
diff ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed ${ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"