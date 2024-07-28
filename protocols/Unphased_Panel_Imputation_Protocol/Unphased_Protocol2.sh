#!/bin/bash

if [[ $# -ne 2 ]]
then
    echo "USAGE: $0 [chromosome id] [Reference resample size]"
    exit 1
fi

chr_id=$1
n_REF_resample=$2

#
# This is unphased query imputation protocol. It does include typed variant aggregation step. This protocol is more flexible than other unphased protocol.
# 

echo "Starting $0 on chromosome=${chr_id}; Ref. resampling=${n_REF_resample}"

# Setup beagle if necessary.
echo "Downloading BEAGLE.."

# Download the BEAGLE jar file; this must match the version in PROXYTYPER.ini.
# wget -c https://faculty.washington.edu/browning/beagle/beagle.01Mar24.d36.jar
wget -c https://faculty.washington.edu/browning/beagle/beagle.13Mar20.38e.jar

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
PROTOCOL_INI=UNPHASED_PROXYTYPER.ini
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

REF_ID=${chr_id}_REF_tags_targets
QUERY_ID=${chr_id}_QUERY_tags
QUERY_targets=${chr_id}_QUERY_targets

./ProxyTyper.sh -check_files "${REF_ID}_genotypes.matrix.gz" "${REF_ID}_variants.bed" "${REF_ID}_subjects.list"
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

./ProxyTyper.sh -check_files "${QUERY_ID}_genotypes.matrix.gz" "${QUERY_ID}_variants.bed" "${QUERY_ID}_subjects.list"
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

###############################################################################################
# Since we are assuming we take unphased genotypes, we need to start from unphased query panel.
# Start from genocoded QUERY panel.
./ProxyTyper.sh -haplo2geno ${QUERY_ID} ${QUERY_ID}_unphased
QUERY_ID=${QUERY_ID}_unphased

################################################################################################
# We are finished with setting up the datasets. Below is where the actual protocol starts.
################################################################################################
################################################################################################
# We will first randomly phase the query panel and then proxize it.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Random phase query tags\t"%e"\t"%M ./ProxyTyper.sh -random_phase_panel ${QUERY_ID} ${QUERY_ID}_phased.random
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_phased.random

# ################################################################################################
# Augment the typed variants.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate tag augmenter model\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_augmenter_model ${QUERY_ID}_variants.bed ${REF_ID}_variants.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

${time_exec} -o ${TIME_MEM_LOG} --append -f "Query tag augment\t"%e"\t"%M ./ProxyTyper.sh -augment_tag_variants ${QUERY_ID} ${QUERY_ID}_tag_augmented
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_tag_augmented

${time_exec} -o ${TIME_MEM_LOG} --append -f "Reference tag augment\t"%e"\t"%M ./ProxyTyper.sh -augment_tag_variants ${REF_ID} ${REF_ID}_tag_augmented
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_tag_augmented

################################################################################################
# Resampling is the first step. We resample the reference and update panel id.
${time_exec} -o ${TIME_MEM_LOG} --append -f "REFERENCE Resampling\t"%e"\t"%M ./ProxyTyper.sh -resample_panel ${REF_ID} ${n_REF_resample} ${REF_ID}_resampled
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_resampled

# We are skipping the resampling for query sample in unphased query protocol (i.e., genocoded query).
# We currently don't use this because it does not work well. 
# ./ProxyTyper.sh -resample_unphased_panel ${QUERY_ID} ${n_QUERY_resample} ${QUERY_ID}_resampled
# QUERY_ID=${QUERY_ID}_resampled
# if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

################################################################################################
# Query site generates the tag proxizer model regions; this is file with hashing parameters.
# It may look like we are using hashing on unphased genotypes, this is not the case since in UNPHASED_PROXYTYPER.ini file (which is used for this protocol), we set the vicinity length to 0, which enforces a 1 parameter hashing model (only bias term). While doing this, we set the maximum number of weights to 1.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate hashing weights\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_proxizer_model ${QUERY_ID}_variants.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Now we have the random phased query panel; we hypothesize that this should give us good accuracy using a 0-vicinity hashing model that only inverts alleles.
# Continue with proxizing the query and ref panels.
# Proxize query panel.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Hash reference tag variants\t"%e"\t"%M ./ProxyTyper.sh -proxize_tag_variants ${REF_ID} ${REF_ID}_proxy.hash
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_proxy.hash

# Now proxize it with first order model.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Hash query tag variants\t"%e"\t"%M ./ProxyTyper.sh -proxize_tag_variants ${QUERY_ID} ${QUERY_ID}_proxy.hash
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_proxy.hash

################################################################################################
# Now we will generate the tag permuter model; this writes one tag permutation BED file.
# We are still in the same coordinate system, we are just permuting the tag variants.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Generate tag permuter model\t"%e"\t"%M ./ProxyTyper.sh -generate_tag_permuter_model ${QUERY_ID}_variants.bed
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi

# Proxize both panels.
${time_exec} -o ${TIME_MEM_LOG} --append -f "Locally permute query tags\t"%e"\t"%M ./ProxyTyper.sh -permute_proxize_tag_variants ${QUERY_ID} ${QUERY_ID}_proxy.permute
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
QUERY_ID=${QUERY_ID}_proxy.permute

${time_exec} -o ${TIME_MEM_LOG} --append -f "Locally permute reference tags\t"%e"\t"%M ./ProxyTyper.sh -permute_proxize_tag_variants ${REF_ID} ${REF_ID}_proxy.permute
if [[ $? -ne 0 ]]; then echo "FAILED @ ${LINENO}"; exit 1; fi
REF_ID=${REF_ID}_proxy.permute

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

# We are now testing just deanonymization:
# De-anonymize the panel to move back to the correct coordinate system.
${time_exec} -o ${TIME_MEM_LOG} --append -f "De-anonymize imputed resampled-query panel\t"%e"\t"%M ./ProxyTyper.sh -deanonymize_coordinates ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID} ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID}_coords_deanonym

# We have the target variants for the query panel. Test accuracy:
${time_exec} -o ${TIME_MEM_LOG} --append -f "Genotype R2\t"%e"\t"%M ./ProxyTyper.sh -genotype_R2 ${QUERY_targets}_variants.bed ${QUERY_LOCAL_IMPUTER_REF_PANEL_ID}_coords_deanonym ${QUERY_targets} ${chr_id}_UNPHASED_PROTOCOL2_accs.txt
./ProxyTyper.sh -summarize_accuracy ${chr_id}_UNPHASED_PROTOCOL2_accs.txt 0.001 0.005
./ProxyTyper.sh -summarize_accuracy ${chr_id}_UNPHASED_PROTOCOL2_accs.txt 0.005 0.01
./ProxyTyper.sh -summarize_accuracy ${chr_id}_UNPHASED_PROTOCOL2_accs.txt 0.01 0.05
./ProxyTyper.sh -summarize_accuracy ${chr_id}_UNPHASED_PROTOCOL2_accs.txt 0.05 1

exit 0
