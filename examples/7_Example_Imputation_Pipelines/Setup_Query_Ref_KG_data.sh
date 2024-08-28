#!/bin/bash

############################################################################################################
# This file implements data extraction step for the example protocols in this directory.
# Run this script to generate the query-typed and reference-{typed+untyped} panels. 
# This script also generates query-untyped variant panel to calculate accuracy of imputation protocols.
############################################################################################################

chr_id=20

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh\", you can it copy it from under scripts/ directory.."
    exit 1
fi

# We include this to be able to download the BEAGLE's correct version in the config files for the phased protocol. Note that unphased protocol also uses the same version.
if [[ ! -f "PHASED_PROXYTYPER.ini" ]]
then
    echo "Could not find the proxytyper config file named \"PHASED_PROXYTYPER.ini\", this script is distributed with this config file, make sure to copy it here."
    exit 1
fi

source PHASED_PROXYTYPER.ini

# Download the BEAGLE jar file; this must match the version in PROXYTYPER.ini.
wget -c https://faculty.washington.edu/browning/beagle/${BEAGLE_JAR}

if [[ ! -f ${BEAGLE_JAR} ]]
then
    echo "Could not find BEAGLE's JAR file @ ${BEAGLE_JAR}"
    exit 1
fi

##################################################################################################
# Import the VCF file and save it to a genotype panel matrix named "KG_chr${chr_id}". We already saw this in previous example.
# Download 1000 Genomes Project's chromosome 20.
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${chr_id}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# We should have the vcf file here:
if [[ ! -f "ALL.chr${chr_id}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ]]
then
    echo "Failed to download the VCF file: ALL.chr${chr_id}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    exit 1
fi

# This is the downloaded VCF file.
ref_VCF="ALL.chr${chr_id}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

all_REF_panel_ID=KG_chr${chr_id}

# Clean all the files for the panel.
rm -f ${all_REF_panel_ID}_*
rm -f -r genetic_maps
rm -f -r plink.*

##################################################################################################
# Import.
is_panel_phased=1
update_variant_ID=1
./ProxyTyper.sh -import_VCF ${ref_VCF} ${is_panel_phased} ${update_variant_ID} ${all_REF_panel_ID}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

##################################################################################################
# Take the unique variants.
./ProxyTyper.sh -uniquefy_panel_variants ${all_REF_panel_ID} ${all_REF_panel_ID}_unique
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Use panel with unique variants from here on.
all_REF_panel_ID=${all_REF_panel_ID}_unique
##################################################################################################
# To cut down computation time, we subselect 10-20meg region on the chromosome.
awk {'if($2>10000000 && $3<20000000){print $0}'} ${all_REF_panel_ID}_variants.bed > focus_reg_vars.bed

./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} focus_reg_vars.bed ${all_REF_panel_ID}_focus
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

all_REF_panel_ID=${all_REF_panel_ID}_focus

##################################################################################################
# Set this identifier aside to be used later.
TYPED_UNTYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}

##################################################################################################
# Subselect 3,000 variants as tag variants; in this case, we also subselect the variants with vertain AAF cutoff:
awk {'split($4, arr, "_");if(arr[4]>0.05 && arr[4]<0.95){print $0}'} ${all_REF_panel_ID}_variants.bed > ALL_TYPED_VARIANT_CANDIDATES.bed
n_all_variants=`wc -l ALL_TYPED_VARIANT_CANDIDATES.bed | awk {'print $1'}`
n_typed=3000
shuf ALL_TYPED_VARIANT_CANDIDATES.bed | head -n ${n_typed} | sort -n -k2,2 > TYPED_VARIANTS.bed

# Get the genotype panel for the typed variants.
./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} TYPED_VARIANTS.bed ${all_REF_panel_ID}_typed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

TYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}_typed

##################################################################################################
# Divide the panels into 1,000 query subjects and 1,504 reference subjects.
n_query_size=1000
shuf ${all_REF_panel_ID}_subjects.list | head -n ${n_query_size} > QUERY_subjects.list

# Exclude the query subjects from all subjects. Query and Reference does not have any overlaps after this step.
grep -v -w -f QUERY_subjects.list ${all_REF_panel_ID}_subjects.list > REF_subjects.list

# Subset the typed variants panel to select the query subjects only.
./ProxyTyper.sh -subset_panel_subjects ${TYPED_VARIANTS_PANEL_ID} QUERY_subjects.list ${TYPED_VARIANTS_PANEL_ID}_QUERY
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Query panel only contains the typed variants.
QUERY_PANEL_ID=${TYPED_VARIANTS_PANEL_ID}_QUERY

# Subset the typed+untyped panel to select reference subjects only. This is our reference panel.
./ProxyTyper.sh -subset_panel_subjects ${TYPED_UNTYPED_VARIANTS_PANEL_ID} REF_subjects.list ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_REFERENCE
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Set this as our reference panel identifier.
REFERENCE_PANEL_ID=${TYPED_UNTYPED_VARIANTS_PANEL_ID}_REFERENCE

# Extract the untyped variants for the query panel.
./ProxyTyper.sh -subset_panel_subjects ${TYPED_UNTYPED_VARIANTS_PANEL_ID} QUERY_subjects.list ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_QUERY
./ProxyTyper.sh -exclude ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_QUERY_variants.bed ${QUERY_PANEL_ID}_variants.bed UNTYPED_VARIANTS.bed
./ProxyTyper.sh -subset_panel_variants ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_QUERY UNTYPED_VARIANTS.bed ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_QUERY_TARGETS
QUERY_TARGET_PANEL_ID=${TYPED_UNTYPED_VARIANTS_PANEL_ID}_QUERY_TARGETS

##################################################################################################
# Download and setup the genetics maps, which is needed by BEAGLE.
rm -f -r genetic_maps
mkdir genetic_maps
BEAGLE_MAPS_URL=https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
./ProxyTyper.sh -setup_BEAGLE_genetic_maps ${BEAGLE_MAPS_URL}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

##################################################################################################
# At this point, we setup the query panel with 1,000 query subjects and only with typed variants; and reference panels with 1,504 subjects with typed+untyped variants.
##################################################################################################
# We first run BEAGLE as if it was being run without any protections:
# Note that the panels are not encode/hashed or anonymized, we just use the plaintext (non-anonymized) genetic maps and run BEAGLE.
genetic_map=${PROXYTYPER_GENETIC_MAPS_DIR}/20.map
do_impute=true
phased_query=0 # Set this to 1 to input phased panel from query.
./ProxyTyper.sh -run_BEAGLE ${QUERY_PANEL_ID} ${REFERENCE_PANEL_ID} ${genetic_map} ${do_impute} ${phased_query} PLAIN_IMPUTED_BEAGLE_DIR
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Import the imputed vcf file that contains typed+untyped variants for the query panel so we can use it with ProxyTyper.
./ProxyTyper.sh -import_VCF PLAIN_IMPUTED_BEAGLE_DIR/imputed.op.vcf.gz 1 1 PLAIN_IMPUTED_BEAGLE
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Extract the untyped variants only. We need this because we will calculate accuracy only on the untyped variants.
./ProxyTyper.sh -exclude ${TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed TYPED_VARIANTS.bed UNTYPED_VARIANTS.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Finally, subset the original typed+untyped panel to extract only query subjects' variant genotypes. This includes the ground truth for the untyped genotypes of query panel.
./ProxyTyper.sh -subset_panel_subjects ${TYPED_UNTYPED_VARIANTS_PANEL_ID} QUERY_subjects.list QUERY_PANEL_TYPED_UNTYPED_VARIANTS
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Calculate accuracy statistics.
./ProxyTyper.sh -genotype_R2 UNTYPED_VARIANTS.bed PLAIN_IMPUTED_BEAGLE QUERY_PANEL_TYPED_UNTYPED_VARIANTS ${chr_id}_plaintext_BEAGLE_accuracy.txt
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
echo "BEAGLE accuracy:"
./ProxyTyper.sh -summarize_accuracy ${chr_id}_plaintext_BEAGLE_accuracy.txt 0 0.001
./ProxyTyper.sh -summarize_accuracy ${chr_id}_plaintext_BEAGLE_accuracy.txt 0.001 0.005
./ProxyTyper.sh -summarize_accuracy ${chr_id}_plaintext_BEAGLE_accuracy.txt 0.005 0.01
./ProxyTyper.sh -summarize_accuracy ${chr_id}_plaintext_BEAGLE_accuracy.txt 0.01 0.05
./ProxyTyper.sh -summarize_accuracy ${chr_id}_plaintext_BEAGLE_accuracy.txt 0.05 1
#################################################################################################
# We are finished with plain BEAGLE, which is our baseline; copy the panel. We need to do this so we can run protocol scripts.
echo "COPYING THE PANELS TO PROTOCOL'S EXPECTED PANELS.."
PROTOCOL_REF_ID=${chr_id}_REF_tags_targets
PROTOCOL_QUERY_ID=${chr_id}_QUERY_tags
PROTOCOL_QUERY_targets=${chr_id}_QUERY_targets
./ProxyTyper.sh -copy_panel ${QUERY_PANEL_ID} ${PROTOCOL_QUERY_ID}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

./ProxyTyper.sh -copy_panel ${QUERY_TARGET_PANEL_ID} ${PROTOCOL_QUERY_targets}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

./ProxyTyper.sh -copy_panel ${REFERENCE_PANEL_ID} ${PROTOCOL_REF_ID}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

echo "COMPLETED.. You can, or should be able to, run the protocols in this directory.. (with some luck..)"