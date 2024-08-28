#!/bin/bash

# This file implements an example of run BEAGLE using ProxyTyper.sh. 
# We demonstrate coordinate/genetic map anonymization and variant decomposition mechanisms on the variants and compare accuracy with the plain run of BEAGLE.
# Note that this is a custom pipeline that is used for demonstration purposes. Also, a lot of data wrangling and setup code is at the beginning, the actual protocol takes up a small portion of the code.
# Also, we run the steps of the protocol sequentially that would normally run in parallel on query and reference (or data owner's site) site.
# We have streamlined the error checks into single lines to make the code more easily readable.

################################################################################################################
# There are two main options in PROXYTYPER.ini file that are used to correctly run BEAGLE:
#BEAGLE_JAR=beagle.01Mar24.d36.jar  : This is the main jar file for BEAGLE.
#JAVA_MEM_FLAG=Xmx400000m           : This is the memory flag for java, this must be set to as high as possible depending on the system to have enough memory.
################################################################################################################

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh\", you can it copy it from under scripts/ directory.."
    exit 1
fi

source PROXYTYPER.ini

# Download the BEAGLE jar file; this must match the version in PROXYTYPER.ini.
wget -c https://faculty.washington.edu/browning/beagle/${BEAGLE_JAR}

if [[ ! -f ${BEAGLE_JAR} ]]
then
    echo "Could not find BEAGLE's JAR file @ ${BEAGLE_JAR}"
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

# This is the downloaded VCF file.
ref_VCF="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

all_REF_panel_ID=KG_chr20

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

# Store the original panel id with all variants in a variable since we will use this later on.
ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}

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
ORIGINAL_TYPED_VARIANTS_PANEL_ID=${all_REF_panel_ID}_typed

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
./ProxyTyper.sh -exclude ${ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed TYPED_VARIANTS.bed UNTYPED_VARIANTS.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Finally, subset the original typed+untyped panel to extract only query subjects' variant genotypes. This includes the ground truth for the untyped genotypes of query panel.
./ProxyTyper.sh -subset_panel_subjects ${ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID} QUERY_subjects.list QUERY_PANEL_TYPED_UNTYPED_VARIANTS
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Calculate accuracy statistics.
./ProxyTyper.sh -genotype_R2 UNTYPED_VARIANTS.bed PLAIN_IMPUTED_BEAGLE QUERY_PANEL_TYPED_UNTYPED_VARIANTS plaintext_BEAGLE_accuracy.txt
#################################################################################################
# We are finished with plain BEAGLE, which is our baseline.
# Below, we do anonymization and decomposition and run BEAGLE as if this is a simple protocol.
# This protocol does not protect typed variant alleles via hashing/permutations. It also does not include typed variant augmentations.
# We present more detailed 

# First generate the coordinate anonymizer model.
./ProxyTyper.sh -generate_coordinate_anonymizer_model ${QUERY_PANEL_ID}_variants.bed ${REFERENCE_PANEL_ID}_variants.bed 
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Anonymize query.
./ProxyTyper.sh -anonymize_coordinates ${QUERY_PANEL_ID} ${QUERY_PANEL_ID}_anon_coords
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Anonymize reference.
./ProxyTyper.sh -anonymize_coordinates ${REFERENCE_PANEL_ID} ${REFERENCE_PANEL_ID}_anon_coords
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Update the panel identifiers for both panels.
QUERY_PANEL_ID=${QUERY_PANEL_ID}_anon_coords
REFERENCE_PANEL_ID=${REFERENCE_PANEL_ID}_anon_coords

# We have the anonymized coordinates for both query and reference panel.
##################################################################################################
# We will now apply decomposition to the untyped variants in reference panel:
# Note that we have already learned about this mechanism in one of the previous examples.
./ProxyTyper.sh -decompose_variants ${QUERY_PANEL_ID}_variants.bed ${REFERENCE_PANEL_ID} ${REFERENCE_PANEL_ID}_decomposed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Update the reference panel's identifier.
REFERENCE_PANEL_ID=${REFERENCE_PANEL_ID}_decomposed

##################################################################################################
# At this point we have the decomposed reference panel on anonymized coordinates at the reference site and query panel on anonymized coordinates. 
# We will run BEAGLE to impute the untyped variants on the query panel.
# It is important to use anonymized genetic maps in BEAGLE since we are in anonymized coordinates.
anonymized_genetic_map=${ANONYMIZED_COORDINATES_DIR}/20_BEAGLE.map
do_impute=true
phased_query=0
./ProxyTyper.sh -run_BEAGLE ${QUERY_PANEL_ID} ${REFERENCE_PANEL_ID} ${anonymized_genetic_map} ${do_impute} ${phased_query} PROXY_IMPUTED_BEAGLE_DIR
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# The resulting VCF file is sent to the query site. We first re-compose, assuming data owner sends the decomposition intervals as they are..
./ProxyTyper.sh -recompose_BEAGLE_VCF_variants PROXY_IMPUTED_BEAGLE_DIR/imputed.op.vcf.gz ${QUERY_PANEL_ID}_subjects.list ${REFERENCE_PANEL_ID}.interval ${QUERY_PANEL_ID}_imputed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Update the query panel identifier. Note that this panel contains typed+untyped variants from BEAGLE (but they are in anonymized coordinates).
QUERY_PANEL_ID=${QUERY_PANEL_ID}_imputed

# Note that the variants are in anonymized coordinates, de-anonymize them:
./ProxyTyper.sh -deanonymize_coordinates ${QUERY_PANEL_ID} ${QUERY_PANEL_ID}_deanonymized
QUERY_PANEL_ID=${QUERY_PANEL_ID}_deanonymized

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# We are now finished and have a fully decoded panel of typed+untyped variants at the query site.
# In the following we extract the known untyped variant genotypes for the query panel, then we estimate accuracy.
# We can calculate accuracy on untyped variants, which were imputed.
./ProxyTyper.sh -exclude ${ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID}_variants.bed TYPED_VARIANTS.bed UNTYPED_VARIANTS.bed
./ProxyTyper.sh -subset_panel_subjects ${ORIGINAL_TYPED_UNTYPED_VARIANTS_PANEL_ID} QUERY_subjects.list QUERY_PANEL_TYPED_UNTYPED_VARIANTS

./ProxyTyper.sh -genotype_R2 UNTYPED_VARIANTS.bed ${QUERY_PANEL_ID} QUERY_PANEL_TYPED_UNTYPED_VARIANTS proxy_protocol_accuracy.txt

# Summarize the accuracy.
awk -v MAF=0.01 '{split($4, arr, "_");if(arr[4]>MAF && arr[4]<(1-MAF)){tot+=$5;n_tot+=1;}}END{print tot/n_tot"\t"tot"\t"n_tot}' proxy_protocol_accuracy.txt
awk -v MAF=0.01 '{split($4, arr, "_");if(arr[4]>MAF && arr[4]<(1-MAF)){tot+=$5;n_tot+=1;}}END{print tot/n_tot"\t"tot"\t"n_tot}' plaintext_BEAGLE_accuracy.txt
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# End of Protocol.