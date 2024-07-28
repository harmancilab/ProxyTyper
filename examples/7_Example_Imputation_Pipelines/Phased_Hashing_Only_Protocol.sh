#!/bin/bash

#
# This script implements the protocol that uses only resampling, hashing, and coordinate anonymization. Also, reference site uses untyped variant decompositions.
# This is a simplification of the Phased_Panel_Imputation_Protocol.sh script by commenting out augmentation, and permutation.
# This script is provided for demonstration purposes to explain how we can build new custom pipelines.
# This protocol can be modified arbirarily by removing blocks of mechanisms. This should, however, be done in encoding and decoding stages consistently and care must be taken that coordinates systems match betweeen modified encoding/decoding stages of the custom protocol.
# In the protocol below, we intentionally did not remove the commented out sections to demonstrate that they can be simply uncommented to include them in new protocols.
# Note that the ordering of the mechanisms could be changed to build new protocols, this is a matter of future research. If you made this far into the examples (you are awesome!) and have some ideas for new mecanisms or protocols, please let us know!!
# 

if [[ $# -ne 3 ]]
then
    echo "USAGE: $0 [chromosome id] [Reference resample size] [Query resample size]"
    exit 1
fi

# Make sure the script is here.
if [[ ! -f "ProxyTyper.sh" ]]
then
    echo "Could not find ProxyTyper script @ \"ProxyTyper.sh}\", you can it copy it from under scripts/ directory.."
    exit 1
fi

chr_id=$1
n_REF_resample=$2
n_QUERY_resample=$3

echo "Starting $0 on chromosome=${chr_id}; Ref. resampling=${n_REF_resample}; Query resampling=${n_QUERY_resample}"

############################################################################################################
# This protocol uses tag augmentation.
# This protocol is the new version.
# This protocol performs final local imputation on de-anonymized plaintext coordinates.
############################################################################################################
# Setup beagle if necessary.
echo "Downloading BEAGLE.."

# Download the BEAGLE jar file; this must match the version in PROXYTYPER.ini.
source UNPHASED_PROXYTYPER.ini
wget -c https://faculty.washington.edu/browning/beagle/${BEAGLE_JAR}
if [[ ! -f ${BEAGLE_JAR} ]]
then
    echo "Could not find BEAGLE's JAR file @ ${BEAGLE_JAR}"
    exit 1
fi

############################################################################################################
# Setup the genetic maps: genetic_maps/ directory name is important since it is used in PROXYTYPER.ini to setup the genetic maps.
BEAGLE_MAPS_URL=https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
./ProxyTyper.sh -setup_BEAGLE_genetic_maps ${BEAGLE_MAPS_URL}
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
############################################################################################################
# Clean the models.
rm -f WALLTIME.LOG
rm -f -r PROXIZATION_MODELS/*
rm -f temp_*

# Copy the appropriate config file.
PROTOCOL_INI=PHASED_PROXYTYPER.ini
if [[ ! -f "${PROTOCOL_INI}" ]]
then
    echo "Could not find ${PROTOCOL_INI}."
    exit 1
fi
cp ${PROTOCOL_INI} PROXYTYPER.ini

script_name=$(basename "$0")
INTERMEDIATE_DATA_DIR=INTERMEDIATE
rm -f -r ${INTERMEDIATE_DATA_DIR}
mkdir ${INTERMEDIATE_DATA_DIR}

cp PROXYTYPER.ini ${INTERMEDIATE_DATA_DIR}

# We need to have Reference, Query-typed, and Query-untyped (for accuracy estimation at the end) panels.
REF_ID=${chr_id}_REF_tags_targets
QUERY_ID=${chr_id}_QUERY_tags
QUERY_targets=${chr_id}_QUERY_targets

./ProxyTyper.sh -check_files "${REF_ID}_genotypes.matrix.gz" "${REF_ID}_variants.bed" "${REF_ID}_subjects.list"
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

./ProxyTyper.sh -check_files "${QUERY_ID}_genotypes.matrix.gz" "${QUERY_ID}_variants.bed" "${QUERY_ID}_subjects.list"
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

./ProxyTyper.sh -check_files "${QUERY_targets}_genotypes.matrix.gz" "${QUERY_targets}_variants.bed" "${QUERY_targets}_subjects.list"
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Copy all the data to INTERMEDIATE
cp ${REF_ID}_subjects.list INTERMEDIATE
cp ${REF_ID}_genotypes.matrix.gz INTERMEDIATE
cp ${REF_ID}_variants.bed INTERMEDIATE

cp ${QUERY_ID}_subjects.list INTERMEDIATE
cp ${QUERY_ID}_genotypes.matrix.gz INTERMEDIATE
cp ${QUERY_ID}_variants.bed INTERMEDIATE

cp ${QUERY_targets}_subjects.list INTERMEDIATE
cp ${QUERY_targets}_genotypes.matrix.gz INTERMEDIATE
cp ${QUERY_targets}_variants.bed INTERMEDIATE

# Update the identifiers.
QUERY_ID=${INTERMEDIATE_DATA_DIR}/${QUERY_ID}
REF_ID=${INTERMEDIATE_DATA_DIR}/${REF_ID}
QUERY_targets=${INTERMEDIATE_DATA_DIR}/${QUERY_targets}

chr_id=`cut -f1,1 ${QUERY_ID}_variants.bed | sort -u | head -n 1`
echo "Imputing ${chr_id}"

TIME_MEM_LOG=${script_name}_TIME_MEM_LOG.txt
time_exec=/usr/bin/time
rm -f -r ${TIME_MEM_LOG}

################################################################################################
# We are finished with setting up the datasets. Below is where the actual protocol starts.
################################################################################################
# We don't augment, we just comment it out.
# # # Augment the tags first. This can be done later as well. We augment for a predefined number of times.
# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate tag augmenter model\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_augmenter_model ${QUERY_ID}_variants.bed ${REF_ID}_variants.bed
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Query tag augment\t"%e"\t"%M ./ProxyTyper.sh -augment_tag_variants ${QUERY_ID} ${QUERY_ID}_tag_augmented
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
# QUERY_ID=${QUERY_ID}_tag_augmented

# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Reference tag augment\t"%e"\t"%M ./ProxyTyper.sh -augment_tag_variants ${REF_ID} ${REF_ID}_tag_augmented
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
# REF_ID=${REF_ID}_tag_augmented

################################################################################################
# we store this for using later while doing the local imputation.
LOCAL_QUERY_PANEL_ID=${QUERY_ID}

################################################################################################
# Resampling is the first step. We resample both query and reference and update panel id.
${time_exec} -o ${TIME_MEM_LOG} --append -f "REFERENCE Resampling\t"%e"\t"%M ./ProxyTyper.sh -resample_panel ${REF_ID} ${n_REF_resample} ${REF_ID}_resampled
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_resampled

${time_exec} -o ${TIME_MEM_LOG} --append -f "QUERY Resampling\t"%e"\t"%M ./ProxyTyper.sh -resample_panel ${QUERY_ID} ${n_QUERY_resample} ${QUERY_ID}_resampled
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_resampled

################################################################################################
# Query site generates the tag proxizer model regions; this is file with hashing parameters.
# This step does not generate new panels.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate hashing weights\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_proxizer_model ${QUERY_ID}_variants.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Proxize both 
${time_exec} -o ${TIME_MEM_LOG} --append -f "Hash reference tag variants\t"%e"\t"%M ./ProxyTyper.sh -proxize_tag_variants ${REF_ID} ${REF_ID}_proxy.hash
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_proxy.hash

${time_exec} -o ${TIME_MEM_LOG} --append -f "Hash query tag variants\t"%e"\t"%M ./ProxyTyper.sh -proxize_tag_variants ${QUERY_ID} ${QUERY_ID}_proxy.hash
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_proxy.hash

################################################################################################
# In this demo, we do not use permutation-based proxizaton, so we comment this out.
# # Now we will generate the tag permuter model; this writes one tag permutation BED file.
# # We are still in the same coordinate system, we are just permuting the tag variants.
# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate tag permuter model\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_permuter_model ${QUERY_ID}_variants.bed
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# # Proxize both panels.
# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Locally permute query tags\t"%e"\t"%M ./ProxyTyper.sh -permute_proxize_tag_variants ${QUERY_ID} ${QUERY_ID}_proxy.permute
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
# QUERY_ID=${QUERY_ID}_proxy.permute

# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Locally permute reference tags\t"%e"\t"%M ./ProxyTyper.sh -permute_proxize_tag_variants ${REF_ID} ${REF_ID}_proxy.permute
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
# REF_ID=${REF_ID}_proxy.permute

################################################################################################
# Server does this -- server cannot send target variants to client.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate anonymizer model\t"%e"\t"%M ./ProxyTyper.sh -generate_coordinate_anonymizer_model ${QUERY_ID}_variants.bed ${REF_ID}_variants.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

${time_exec} -o ${TIME_MEM_LOG} --append -f "Anonymize query coordinates\t"%e"\t"%M ./ProxyTyper.sh -anonymize_coordinates ${QUERY_ID} ${QUERY_ID}_anon_coords
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_anon_coords

${time_exec} -o ${TIME_MEM_LOG} --append -f "Anonymize reference coordinates\t"%e"\t"%M ./ProxyTyper.sh -anonymize_coordinates ${REF_ID} ${REF_ID}_anon_coords
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_anon_coords

################################################################################################

# At this point, we are in the new coordinate system.
# Client still does not know anything about the targets. Only tag information is sent to them.
# Decompose the untyped variants. 
# We get a new proxy panel for reference dataset so we update the identifier.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Decompose reference untyped variants\t"%e"\t"%M ./ProxyTyper.sh -decompose_variants ${QUERY_ID}_variants.bed ${REF_ID} ${REF_ID}_decomp
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_decomp
DECOMPOSING_INTERVAL=${REF_ID}.interval

################################################################################################
# Central server:
# Run BEAGLE with decomposed reference; we need the anonymized genetic map for BEAGLE.
anon_beagle_map=PROXIZATION_MODELS/anon_genetic_maps/${chr_id}.map
${time_exec} -o ${TIME_MEM_LOG} --append -f "BEAGLE (Resampled/Proxy Query)\t"%e"\t"%M ./ProxyTyper.sh -run_BEAGLE ${QUERY_ID} ${REF_ID} ${anon_beagle_map} true 0 ${INTERMEDIATE_DATA_DIR}/QUERY_RESAMP_IMPUTE_beagle_dir
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

QUERY_RESAMP_IMPUTE_VCF=${INTERMEDIATE_DATA_DIR}/QUERY_RESAMP_IMPUTE_beagle_dir/imputed.op.vcf.gz

if [[ ! -f ${QUERY_RESAMP_IMPUTE_VCF} ]]
then
    echo "${QUERY_RESAMP_IMPUTE_VCF} does not exist, imputation of the resampled/proxized query seems to have failed.."
    exit 1
fi
################################################################################################
# We need to first recompose the target variants on the resampled query subjects that we just imputed.
# This ensures that targets we have are mapped back to anonymized coordinates system.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Recompose reference untyped variants\t"%e"\t"%M ./ProxyTyper.sh -recompose_BEAGLE_VCF_variants ${QUERY_RESAMP_IMPUTE_VCF} \
${QUERY_ID}_subjects.list \
${REF_ID}.interval \
${INTERMEDIATE_DATA_DIR}/QUERY_RESAMP_IMPUTED
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_LOCAL_IMPUTER_REF_PANEL_ID=${INTERMEDIATE_DATA_DIR}/QUERY_RESAMP_IMPUTED

# Now we have recomposed variants for the resampled query panel. We will now use this to do one more local imputation at client's side.

# We have the imputation panel for re-imputing tag variants. First deanonymize_coordinates the tag/target coordinates:
${time_exec} -o ${TIME_MEM_LOG} --append -f "De-anonymize imputed resampled-query panel\t"%e"\t"%M ./ProxyTyper.sh -deanonymize_coordinates ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID} ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID}_coords_deanonym
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_LOCAL_IMPUTER_REF_PANEL_ID=${QUERY_LOCAL_IMPUTER_REF_PANEL_ID}_coords_deanonym

# Now, we have the local imputation panel whose tags are proxized. We will proxize tags of the 
# original tag augmented Query panel and then use the proxized original panel tags to do imputation with this panel.
# We set LOCAL_QUERY_PANEL_ID to our very original typed variant panel (you can see this before the first step of the protocol above) for the query (before resampling) so we use it here to do the local imputation.
LOCAL_QUERY_ID=${LOCAL_QUERY_PANEL_ID}
./ProxyTyper.sh -check_files ${LOCAL_QUERY_ID}_variants.bed ${LOCAL_QUERY_ID}_subjects.list ${LOCAL_QUERY_ID}_variants.bed
if [[ $? -ne 0 ]]
then
    echo "Sanity check failed: Could not find the panel for ${LOCAL_QUERY_ID} (${LINENO})"
    exit 1
fi

# Make sure we are not overwriting anything. This should hold since the protocol starts with resampling (not hashing).
# We should not be overwriting this file otherwise sth went wrong with the protocol. 
if [[ -f ${LOCAL_QUERY_ID}_proxy.hash_subjects.list ]]
then
    echo "SANITY CHECK ERROR: We will overwrite a file, exiting (${LINENO})."
    exit 1
fi

${time_exec} -o ${TIME_MEM_LOG} --append -f "Hash-proxize query tags\t"%e"\t"%M ./ProxyTyper.sh -proxize_tag_variants ${LOCAL_QUERY_ID} ${LOCAL_QUERY_ID}_proxy.hash
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
LOCAL_QUERY_ID=${LOCAL_QUERY_ID}_proxy.hash

# # We do not use permute proxization in this demo, so we don't need to permute proxize tags. If you uncomment permute proxy step above, make sure to uncomment it here..
# ${time_exec} -o ${TIME_MEM_LOG} --append -f "Permute-proxize hash-proxized query tags\t"%e"\t"%M ./ProxyTyper.sh -permute_proxize_tag_variants ${LOCAL_QUERY_ID} ${LOCAL_QUERY_ID}_proxy.permute
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
# LOCAL_QUERY_ID=${LOCAL_QUERY_ID}_proxy.permute

# The local panel and the resampled imputed panel (panel id is LOCAL_QUERY_PANEL_ID, as received from central server's imputation above) are now comparable.
# We need to do a local imputation at the query site... 

# WAIT!!! WHY DO WE NEED A SECOND IMPUTATION??? (In case ProxyTyper paper did not make any sense, which it probably did not..)
# REASON FOR A LOCAL IMPUTATION STEP is that we resampled our query panel. Because of resampling, the imputed untyped variant genotypes from the central server's first step imputation above correspond to a mosaic of the query panel subjects.
# To retrieve the imputed genotypes for the original (non-resampled) query subjects, query site uses BEAGLE again in this local imputation step to impute the genotypes for untyped variants. 
# Query site use the imputed typed+untyped variants from central server as a reference in this local imputation step.
# The input to this imptuation step is, obviously, the original query panel.

# We must use original beagle genetic map here since we have deanonymized the coordinates. We can skip deanonymization step above and do imputation immediately on the anonymized coordinates but then we would need to deanonymize the imputed targets. In that case, we would need to use anonymized genetic maps.
beagle_maps_dir=genetic_maps
beagle_genetic_map_file=${beagle_maps_dir}/${chr_id}.map
${time_exec} -o ${TIME_MEM_LOG} --append -f "BEAGLE (local query)\t"%e"\t"%M ./ProxyTyper.sh -run_BEAGLE ${LOCAL_QUERY_ID} ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID} ${beagle_genetic_map_file} true 0 ${INTERMEDIATE_DATA_DIR}/QUERY_PLAIN_REIMPUTE_beagle_dir
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Import this VCF; this VCF has the final targets (which were recomposed and deanonymized)
add_AF_info_2_id=1  # Make sure to include this flag since ProxyTyper uses it to calculate the minor allele accuracy metrics (not related to R2).
is_VCF_phased=1     # This flag indicates that the VCF is phased, which it is.
${time_exec} -o ${TIME_MEM_LOG} --append -f "VCF Import\t"%e"\t"%M ./ProxyTyper.sh -import_VCF ${INTERMEDIATE_DATA_DIR}/QUERY_PLAIN_REIMPUTE_beagle_dir/imputed.op.vcf.gz ${is_VCF_phased} ${add_AF_info_2_id} ${INTERMEDIATE_DATA_DIR}/QUERY_PLAIN_IMPUTE
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

################################################################################################
# We are finished with protocol here. Below is the accuracy estimation.
################################################################################################
# Get accuracy. We use query targets only here.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Genotype R2\t"%e"\t"%M ./ProxyTyper.sh -genotype_R2 ${QUERY_targets}_variants.bed ${INTERMEDIATE_DATA_DIR}/QUERY_PLAIN_IMPUTE ${QUERY_targets} ${chr_id}_CUSTOM_PHASED_PROTOCOL1_accs.txt

# Report the overall accuracy.
./ProxyTyper.sh -summarize_accuracy ${chr_id}_CUSTOM_PHASED_PROTOCOL1_accs.txt 0.001 0.005
./ProxyTyper.sh -summarize_accuracy ${chr_id}_CUSTOM_PHASED_PROTOCOL1_accs.txt 0.005 0.01
./ProxyTyper.sh -summarize_accuracy ${chr_id}_CUSTOM_PHASED_PROTOCOL1_accs.txt 0.01 0.05
./ProxyTyper.sh -summarize_accuracy ${chr_id}_CUSTOM_PHASED_PROTOCOL1_accs.txt 0.05 1

exit 0

# The End.