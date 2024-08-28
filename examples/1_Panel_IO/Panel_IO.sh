#!/bin/bash

# This file implements the VCF importing, some simple subsetting operations, and exporting a VCF file.
# We demonstrate usage of the I/O options.
# There are error checks after each command as an example of using error returned by the ProxyTyper.sh script.

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh\", you can it copy it from under scripts/ directory.."
    exit 1
fi

# Download 1000 Genomes Project's chromosome 20.
wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# We should have the vcf file here:
if [[ ! -f "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ]]
then
    echo "Failed to download the VCF file: ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    exit 1
fi

ref_VCF="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# Import the VCF file and save it to a genotype panel matrix named "KG_chr20"
# update_variant_ID: This option is set to 1 to update the variant identifiers with the ref/alt alleles and alternate allele frequencies.
# is_panel_phased: This option specifies whether this is a phased or unphased panel. Since 1000 Genomes data is phased, we set it to 1.
# If the option does not match the genotypes in the file, this script will still process the file but will write a lot of warnings.
# Also, this command uses multithreading to import the VCF file. The number of threads is set in PROXYTYPER.ini file, "n_threads=40" option.
# You can change the number of threads by modifying this line in PROXYTYPER.ini file.
all_REF_panel_ID=KG_chr20
is_panel_phased=1
update_variant_ID=1

# Clean all the files for the panel.
rm -f ${all_REF_panel_ID}_*

./ProxyTyper.sh -import_VCF ${ref_VCF} ${is_panel_phased} ${update_variant_ID} ${all_REF_panel_ID}
if [[ $? -ne 0 ]]
then
    echo "VCF importing failed"
    exit 1
fi

# The import generates these files:
# KG_chr20_genotypes.matrix.gz: This is a binary file that contains the genotype matrix.
# KG_chr20_variants.bed: This is a bed file that contains the variants in a bed formatted file.
# KG_chr20_subjects.list: This is a text file that contains the list of subjects in the VCF file.
# The text files can be viewed with a text editor.
# The names of these files are formed using the panel id we used, i.e., "all_REF_panel_ID=KG_chr20".
# All mechanisms in ProxyTyper use this format to process and save panel genotypes. 
# As we process the panels, we update the panel identifiers to indicate which mechanisms have operated
# on the panels.
# By default, ProxyTyper removes indels and multiallelic variants 

##################################################################################################

# Below, we demonstrate several commands to perform simple operations in the panel for exploring it.
# Note that each panel processing step consists of taking a panel as input and updating the panel id with 
# the processed output panel. This simplifies the genotype panel processing pipeline building.

##################################################################################################

./ProxyTyper.sh -uniquefy_panel_variants ${all_REF_panel_ID} ${all_REF_panel_ID}_unique
if [[ $? -ne 0 ]]
then
    echo "Unique panel extraction failed"
    exit 1
fi

all_REF_panel_ID=${all_REF_panel_ID}_unique

##################################################################################################

# It is a good idea to sort the panel after importing the panel.
./ProxyTyper.sh -sort_panel ${all_REF_panel_ID} ${all_REF_panel_ID}_sorted
if [[ $? -ne 0 ]]
then
    echo "Panel sorting failed"
    exit 1
fi

all_REF_panel_ID=${all_REF_panel_ID}_sorted

################################
# We can extract the variants between positions 20:20,000,000-30,000,000
awk {'if($2>20000000 && $3<30000000){print $0}'} ${all_REF_panel_ID}_variants.bed > subset_variants.bed

# Extract the subset of variants.
./ProxyTyper.sh -subset_panel_variants ${all_REF_panel_ID} subset_variants.bed ${all_REF_panel_ID}_variant_subset
if [[ $? -ne 0 ]]
then
    echo "Variant subset failed."
    exit 1
fi

# Update the panel identifier in the protocol.
all_REF_panel_ID=${all_REF_panel_ID}_variant_subset
################################
# We now extract a random sample of 1,000 subjects:
shuf ${all_REF_panel_ID}_subjects.list | head -n 1000 > subset_subjects.list

# extract the subset panel.
./ProxyTyper.sh -subset_panel_subjects ${all_REF_panel_ID} subset_subjects.list ${all_REF_panel_ID}_subject_subset
if [[ $? -ne 0 ]]
then
    echo "Subject subset failed."
    exit 1
fi

# Update the panel identifier in the protocol.
all_REF_panel_ID=${all_REF_panel_ID}_subject_subset
################################

# Finally, we can estimate the alternate allele frequencies to just do something with the panel.
./ProxyTyper.sh -calculate_panel_AAF ${all_REF_panel_ID} ${all_REF_panel_ID}_AAF.bed
if [[ $? -ne 0 ]]
then
    echo "AF calculation failed."
    exit 1
fi

# ${all_REF_panel_ID}_AAF.bed is a BED formatted file, containing the alternate allele frequencies @ 5th column.
# It can be viewed in a text editor.
# We don't update the panel identifier here since we did not update the panel..

################################
# We can now save the panel file as a VCF formatted file:
./ProxyTyper.sh -export_VCF ${all_REF_panel_ID} hg19 1 ${all_REF_panel_ID}.vcf.gz
if [[ $? -ne 0 ]]
then
    echo "VCF exporting failed."
    exit 1
fi

echo "The exported VCF file can be found @ KG_chr20_unique_sorted_variant_subset_subject_subset.vcf.gz"..
echo "Completed.."

