#!/bin/bash

# Check to make sure config file exists.
config_file=PROXYTYPER.ini

if [[ ! -f ${config_file} ]]
then
    echo "Could not find ${config_file}"
    exit 1
fi

if [[ $# -lt 1 ]]
then
    echo "USAGE: $0 [Options] [Arguments]
Options:
    Data Input:
        -check_files
        -check_panel_sortedness
        -import_VCF
        -get_VCF_subjects
        -export_VCF
    Panel Processing:
        -copy_panel
        -sort_panel
        -delete_panel
        -uniquefy_panel_variants
        -subset_panel_variants
        -subset_panel_subjects
        -calculate_panel_AAF
    Resampling:
        -setup_BEAGLE_genetic_maps
        -resample_panel
        -resample_unphased_panel
    Variant Protection:
        -generate_tag_augmenter_model
        -augment_tag_variants
        -generate_tag_proxizer_model 
        -proxize_tag_variants 
        -generate_tag_permuter_model
        -permute_proxize_tag_variants
        -random_phase_panel
    Coordinate/Genetic Map Anonymization:
        -generate_coordinate_anonymizer_model
        -anonymize_coordinates
        -deanonymize_coordinates
    Variant Decomposition:
        -decompose_variants 
        -recompose_BEAGLE_VCF_variants
    Imputation/Phasing (BEAGLE only):
        -run_BEAGLE
    BED Processing:
        -intersect
        -exclude
    Accuracy:
        -genotype_R2
        -summarize_accuracy"
    exit 1
fi

# Load the config.
source ${config_file}

exec_check=`type -P ${PROXYTYPER_EXEC}`
if [[ ${exec_check} == "" ]]
then
    echo "Could not find ProxyTyper @ ${PROXYTYPER_EXEC}"
    exit 1
fi

cmd_option=$1

######################################################################
# Validate system before running anything.
${PROXYTYPER_EXEC} -validate_system ${n_threads} >& /dev/null
if [[ $? -ne 0 ]]
then
    echo "System check failed (${LINENO})"
    exit 1
fi
######################################################################

if [[ ${cmd_option} == "-cat2stdout" ]]
then
	if [[ $# -ne 2 ]]
	then
		echo "USAGE: $0 $1 [Input file]" >&2
		exit 1
	fi

	fullfile=$2

	filename=$(basename "$fullfile")

	#echo $filename

	extension="${filename##*.}"

	#echo $extension

	#filename="${filename%.*}"

	if [[ "${extension}" == "vcf" ]]
	then
		cat $fullfile
	elif [[ "${extension}" == "txt" ]]
	then
		cat $fullfile
	elif [[ "${extension}" == "gz" ]]
	then
		gzip -cd $fullfile
	elif [[ "${extension}" == "zip" ]]
	then
		unzip -c $fullfile
	else
		cat $fullfile
	fi

	exit 0
fi

if [[ ${cmd_option} == "-get_now_str" ]]
then
    now_str=`date +"["%m/%d/%y::%H:%M:%S"]>"`
    echo ${now_str}

    exit 0
fi

######################################################################################################################################################
# Above options may report text that is used by other tools, don't put any unoptioned text there that output messages.
######################################################################################################################################################
if [[ ! -d "${MODELS_DIR}" ]]
then
    echo "Creating \"${MODELS_DIR}\""
    mkdir "${MODELS_DIR}"
fi

start_time=$(date +%s)

if [[ "${cmd_option}" == "-check_directory" ]]
then
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [Directory Path]"
        exit 1
    fi

    directory_path=$2
    if [[ -d "${directory_path}" ]]
    then
        exit 0
    fi

    echo "Could not find directory \"${directory_path}\".."    
    exit 1
fi

###################################################################################################################

if [[ "${cmd_option}" == "-check_file" ]]
then
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [File path]"
        exit 1
    fi

    file_path=$2
    if [[ -f "${file_path}" ]]
    then
        exit 0
    fi

    echo "Could not find file \"${file_path}\".."    
    exit 1
fi

###################################################################################################################

if [[ "${cmd_option}" == "-check_files" ]]
then
    if [[ $# -lt 2 ]]
    then
        echo "USAGE: $0 $1 [File paths]"
        exit 1
    fi

    n_args=$#
    file_i=`seq 2 ${n_args}`
    for i_f in ${file_i[@]}
    do
        echo "Checking ${!i_f}"
        if [[ ! -f "${!i_f}" ]]
        then
            echo "Could not find file \"${!i_f}\".."
            exit 1
        fi
    done

    exit 0
fi

###################################################################################################################

# Checks strict sortedness in a BED file, including overlaps.
if [[ ${cmd_option} == "-check_BED_sortedness" ]]
then
	if [[ $# -lt 2 ]]
	then
		echo "$0 $1 [BED file]" >&2
		exit 1
	fi

    BED_file=$2

    $0 -check_files "${BED_file}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    echo "Checking sortedness in ${BED_file}.."

    # Check if the variants are sorted.
    panel_is_unsorted=`awk 'BEGIN{FS="\t";cur_pos=0;panel_is_unsorted=0;}{if(NR>1){if(cur_pos >= $2){panel_is_unsorted=1}};cur_pos=$2}END{print panel_is_unsorted}' ${BED_file}`

    if [[ ${panel_is_unsorted} -eq 1 ]]
    then
        awk 'BEGIN{FS="\t";cur_pos=0;panel_is_unsorted=0;}{if(NR>1){if(cur_pos >= $2){panel_is_unsorted=1;print "MISSORTED LINE: "$0}};cur_pos=$2}END{print panel_is_unsorted}' ${BED_file}
    fi

    exit ${panel_is_unsorted}
fi

###################################################################################################################

if [[ ${cmd_option} == "-check_panel_sortedness" ]]
then
	if [[ $# -lt 2 ]]
	then
		echo "$0 $1 [Panel prefix]" >&2
		exit 1
	fi

    panel_prefix=$2

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    echo "Checking sortedness in ${panel_prefix}.."

    # Check if the variants are sorted.
    panel_is_unsorted=`awk 'BEGIN{FS="\t";cur_pos=0;panel_is_unsorted=0;}{if(NR>1){if(cur_pos >= $2){panel_is_unsorted=1}};cur_pos=$2}END{print panel_is_unsorted}' ${panel_prefix}_variants.bed`

    if [[ ${panel_is_unsorted} -eq 1 ]]
    then
        awk 'BEGIN{FS="\t";cur_pos=0;panel_is_unsorted=0;}{if(NR>1){if(cur_pos >= $2){panel_is_unsorted=1;print "MISSORTED LINE: "$0}};cur_pos=$2}END{print panel_is_unsorted}' ${panel_prefix}_variants.bed
    fi

    exit ${panel_is_unsorted}
fi

###################################################################################################################

if [[ ${cmd_option} == "-get_VCF_subjects" ]]
then
	if [[ $# -lt 3 ]]
	then
		echo "$0 $1 [VCF file path] [VCF subject list file]" >&2
		exit 1
	fi

	vcf_file=$2
    subjects_list_file=$3

    $0 -check_files ${vcf_file}
    if [[ $? -ne 0 ]]
    then
        echo "Could not find VCF file @ ${vcf_file} (${LINENO})"
        exit 1
    fi

    echo "Extracting VCF subject identifiers from ${vcf_file}."
	$0 -cat2stdout $vcf_file | head -n 1000 | awk {'if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i}}'} > ${subjects_list_file}

    if [[ ! -f ${subjects_list_file} ]]
    then
        echo "Failed to extract subjects list from ${vcf_file} (${LINENO})"
        exit 1
    fi

    exit $?
fi

###################################################################################################################

if [[ ${cmd_option} == "-import_VCF" ]]
then
	if [[ $# -lt 5 ]]
	then
		echo "$0 $1 [VCF file path] [Is VCF phased? (0/1)] [Add Ref/Alt allele-AF Info to name (0/1)] [Output prefix]" >&2
		exit 1
	fi

	vcf_file=$2
    panel_is_phased=$3
    add_AF_info_2_id=$4
    output_prefix=$5

    # [add_AF_info_2_id] adds ref/alt alleles and AAF information to the variant identifiers. 
    # This information is needed while calculating genotype-R2 (for assigning minor alleles). 
    # This can be changed to keep the imported variant ids original as in the VCF file.

    file_path=$2
    $0 -check_file "${vcf_file}"
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

	GENOME_SEQ_DIR="none"

    rm -f ${output_prefix}_subjects.list ${output_prefix}_variants.bed ${output_prefix}_genotypes.matrix.gz

    echo "Extracting VCF subject identifiers."
	$0 -cat2stdout $vcf_file | head -n 1000 | awk {'if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i}}'} > vcf_subjects.list

    ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF_no_buffer_multithreaded ${vcf_file} vcf_subjects.list ${panel_is_phased} ${add_AF_info_2_id} ${n_threads} ${output_prefix}

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Report time
    end_time=$(date +%s)

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

	exit $?
fi

########################################################################################################################
if [[ "${cmd_option}" == "-setup_BEAGLE_genetic_maps" ]]
then
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [4-column PLINK formatted (chr rsid cM posn) genetic map file URL, e.g., https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip]"
        exit 1
    fi
    
    BEAGLE_maps_URL=$2

    beagle_map_archive=beagle_maps.zip

    echo "Downloading BEAGLE maps from ${BEAGLE_maps_URL}"
    curl -s "${BEAGLE_maps_URL}" -o "${beagle_map_archive}"
    unzip -o ${beagle_map_archive}

    rm -f -r ${PROXYTYPER_GENETIC_MAPS_DIR}
    if [[ ! -d "${PROXYTYPER_GENETIC_MAPS_DIR}" ]]
    then
        #echo "Make sure genetic maps directory exists under: ${PROXYTYPER_GENETIC_MAPS_DIR}"
        #exit 1
        echo "Creating genetic maps directory @ ${PROXYTYPER_GENETIC_MAPS_DIR}"
        mkdir "${PROXYTYPER_GENETIC_MAPS_DIR}"
    fi

    if [[ ! -d "${PROXYTYPER_GENETIC_MAPS_DIR}" ]]
    then
        echo "Could not create genetic maps directory @ ${PROXYTYPER_GENETIC_MAPS_DIR}"
        exit 1
    fi

    # Convert the autosomes, only.
    for cur_chr in $(seq 1 22)
    do
        n_files=`ls *.chr${cur_chr}.*.map | wc -l | awk {'print $1'}`
        if [[ ${n_files} -ne 1 ]]
        then
            echo "Could not find the map file for ${cur_chr}, multiple files exist.."
            continue
        fi
        cur_map_file=`ls *.chr${cur_chr}.*.map`
        echo "Copying ${cur_chr} (${cur_map_file}).."
        cp ${cur_map_file} ${PROXYTYPER_GENETIC_MAPS_DIR}/${cur_chr}.map
        if [[ $? -ne 0 ]]
        then
            echo "Could not copy ${cur_map_file}, most likely there is a permission issue.. (${LINENO})"
            exit 1
        fi
    done

    exit $?
fi

###################################################################################################################

if [[ "${cmd_option}" == "-resample_panel" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Resample size] [Output prefix]"
        exit 1
    fi
    
    panel_prefix=$2
    n_resample_size=$3
    output_prefix=$4

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -d "${PROXYTYPER_GENETIC_MAPS_DIR}" ]]
    then
        echo "Could not find ProxyTyper's source of genetic maps @ \"${PROXYTYPER_GENETIC_MAPS_DIR}\""
        exit 1
    fi

    # Check resample size:
    panel_size=`wc -l ${panel_prefix}_subjects.list | awk {'print $1'}`

    if [[ ${panel_size} -gt ${n_resample_size} ]]
    then
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
    fi

    echo "Resampling ${panel_prefix} and saving to ${output_prefix}"

    rm -f ${output_prefix}_subjects.list ${output_prefix}_variants.bed ${output_prefix}_genotypes.matrix.gz
    ${PROXYTYPER_EXEC} -resample_phased_haplotypes_per_state_only_sampling ${panel_prefix} ${panel_prefix}_subjects.list \
            ${PROXYTYPER_GENETIC_MAPS_DIR} ${n_resample_size} ${pre_mosaic_Ne} ${min_cM_delta_per_anchoring_vars} ${pre_mosaic_allele_eps} \
            ${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
            ${haplocoded_flag} \
            ${n_threads} \
            1 250000000 \
            ${save_recomb_patterns} \
            ${output_prefix} >& ${panel_prefix}_RESAMPLING.op

    ####################################################################################
    # We can summarize the resampling information, if it is saved.
    #${PROXYTYPER_EXEC} -summarize_sampled_segments_per_resampled_haplotype_info INTERMEDIATE/REF_tags_targets_resampled_resampling_pattern_signal_var_regs_20.sigbed.gz INTERMEDIATE/REF_tags_targets_resampled_subjects.list INTERMEDIATE/REF_tags_targets_subjects.list summarized.txt
    ####################################################################################
    
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

#########################################################################################################


if [[ "${cmd_option}" == "-calculate_panel_AAF" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Alternate allele frequency BED file]"
        exit 1
    fi
    
    panel_prefix=$2
    AAF_BED_file=$3

        $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi
    
    ${PROXYTYPER_EXEC} -compute_AF_per_geno_signal_regions ${panel_prefix} ${panel_prefix}_subjects.list ${AAF_BED_file}

    $0 -check_files "${AAF_BED_file}" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    exit $?
fi

#########################################################################################################

if [[ "${cmd_option}" == "-uniquefy_panel_variants" ]]
then
if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    unique_panel_prefix=$3
    
    ${PROXYTYPER_EXEC} -uniquefy_genotype_signal_regions ${panel_prefix} ${panel_prefix}_subjects.list ${unique_panel_prefix}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    $0 -check_files "${unique_panel_prefix}_subjects.list" \
"${unique_panel_prefix}_variants.bed" \
"${unique_panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"
    exit $?
fi

#########################################################################################################

if [[ "${cmd_option}" == "-subset_panel_variants" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Variants BED file] [Output prefix]"
        exit 1
    fi
    
    panel_prefix=$2
    variants_BED_file=$3
    output_prefix=$4

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" \
"${variants_BED_file}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -extract_genotype_signals_per_region_list ${panel_prefix} ${panel_prefix}_subjects.list ${variants_BED_file} ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi


if [[ "${cmd_option}" == "-subset_panel_subjects" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Subject list file] [Output prefix]"
        exit 1
    fi
    
    panel_prefix=$2
    subject_list_file=$3
    output_prefix=$4

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" \
"${subject_list_file}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -extract_genotype_signals_per_subsample_list ${panel_prefix} ${panel_prefix}_subjects.list ${subject_list_file} ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

#########################################################################################################

if [[ "${cmd_option}" == "-resample_unphased_panel" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Resample size] [Output prefix]"
        exit 1
    fi
    
    panel_prefix=$2
    n_resample_size=$3
    output_prefix=$4

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -d "${PROXYTYPER_GENETIC_MAPS_DIR}" ]]
    then
        echo "Could not find ProxyTyper's source of genetic maps @ \"${PROXYTYPER_GENETIC_MAPS_DIR}\""
        exit 1
    fi

    # Check resample size:
    panel_size=`wc -l ${panel_prefix}_subjects.list | awk {'print $1'}`

    if [[ ${panel_size} -gt ${n_resample_size} ]]
    then
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
        echo "***WARNING::It looks like sampling will decrease the size of the panel (panel_size=${panel_size}, resample-size=${n_resample_size}), this is not recommended..***"
    fi

    echo "Resampling ${panel_prefix} and saving to ${output_prefix}"

    rm -f ${output_prefix}_subjects.list ${output_prefix}_variants.bed ${output_prefix}_genotypes.matrix.gz
    ${PROXYTYPER_EXEC} -resample_unphased_genotypes_per_state_only_sampling ${panel_prefix} ${panel_prefix}_subjects.list \
            ${PROXYTYPER_GENETIC_MAPS_DIR} ${n_resample_size} ${pre_mosaic_Ne} ${min_cM_delta_per_anchoring_vars} ${pre_mosaic_allele_eps} \
            ${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
            ${n_threads} \
            1 250000000 \
            ${save_recomb_patterns} \
            ${output_prefix} >& ${panel_prefix}_UNPHASED_RESAMPLING.op

    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################
if [[ "${cmd_option}" == "-generate_tag_proxizer_model" ]]
then
    #-resample_panel [Panel prefix]
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [Tag variants BED file]"
        exit 1
    fi
    
    tag_vars_BED=$2

    $0 -check_file "${tag_vars_BED}"
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

    rm -f ${HASHING_MODEL_FILE}
    ${PROXYTYPER_EXEC} -generate_save_per_site_mixing_parameters_LD_aware ${tag_vars_BED} \
    ${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
    ${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${PROXYTYPER_GENETIC_MAPS_DIR} ${HASHING_MODEL_FILE}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -f "${HASHING_MODEL_FILE}" ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

if [[ "${cmd_option}" == "-random_phase_panel" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -random_phase_genotypes ${panel_prefix} ${panel_prefix}_subjects.list ${output_prefix}

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################
if [[ "${cmd_option}" == "-proxize_tag_variants" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    allele_err_eps=0.0
    # ${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters \
    # ${panel_prefix} ${panel_prefix}_subjects.list \
    # ${HASHING_MODEL_FILE} ${allele_err_eps} ${n_threads} \
    # ${output_prefix}

    # This should run faster???
    rm -f ${output_prefix}_subjects.list ${output_prefix}_variants.bed ${output_prefix}_genotypes.matrix.gz
    ${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer \
    ${panel_prefix} ${panel_prefix}_subjects.list \
    ${HASHING_MODEL_FILE} ${allele_err_eps} ${n_threads} \
    ${output_prefix}
    res=$?
    if [[ ${res} -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################

if [[ "${cmd_option}" == "-generate_tag_permuter_model" ]]
then
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [Tag variants BED]"
        exit 1
    fi

    tag_vars_BED=$2

    $0 -check_file "${tag_vars_BED}"
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

    rm -f ${TAG_PERMUTING_MODEL_FILE}
    ${PROXYTYPER_EXEC} -generate_permute_proxizing_parameters ${tag_vars_BED} ${perm_proxization_vicinity_size} \
    ${perm_proxy_perm_prob} ${perm_proxy_geno_inversion_prob} ${TAG_PERMUTING_MODEL_FILE}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -f "${TAG_PERMUTING_MODEL_FILE}" ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Validate the formatting of the file:
    n_unique_n_cols=`awk 'BEGIN{FS="\t"}{print NF}' ${TAG_PERMUTING_MODEL_FILE} | sort -u | wc -l | awk {'print $1'}`
    if [[ ${n_unique_n_cols} -ne 1 ]]
    then
        echo "Sanity check failed (${LINENO}): There are ${n_unique_n_cols} in ${TAG_PERMUTING_MODEL_FILE}"
        exit 1
    fi

    # Validate the number of columns:
    n_cols=`awk 'BEGIN{FS="\t"}{print NF}' ${TAG_PERMUTING_MODEL_FILE} | sort -u`
    if [[ ${n_cols} -ne ${N_COLS_PER_PERMUTE_MODEL_FILE} ]]
    then
        echo "Sanity check failed (${LINENO}): Number of columns is not as expected in ${TAG_PERMUTING_MODEL_FILE}: ${N_COLS_PER_PERMUTE_MODEL_FILE}!=${n_cols}"
        exit 1
    fi

    # Write mismatching statistics:
    n_permuted_vars=`awk 'BEGIN{FS="\t"}{if($7!=$8){n_permed+=1}}END{print n_permed}' ${TAG_PERMUTING_MODEL_FILE}`
    echo "Found ${n_permuted_vars} permuted variant positions in ${TAG_PERMUTING_MODEL_FILE}"

    # Calculate the elapsed time, return.
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################

if [[ "${cmd_option}" == "-permute_proxize_tag_variants" ]]
then
    #-resample_panel [Panel prefix]
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3
 
    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi
 
    # Permute proxize the QUERY data.
    rm -f ${output_prefix}_subjects.list ${output_prefix}_variants.bed ${output_prefix}_genotypes.matrix.gz
    ${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions ${panel_prefix} ${panel_prefix}_subjects.list \
    ${TAG_PERMUTING_MODEL_FILE} ${n_threads} ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################
if [[ "${cmd_option}" == "-generate_coordinate_anonymizer_model" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [tag variants BED (query site variants)] [tag+target variants BED (reference site variants)]"
        exit 1
    fi

    tag_vars_BED=$2
    tag_target_vars_BED=$3

        $0 -check_files "${tag_vars_BED}" "${tag_target_vars_BED}" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_file "${tag_vars_BED}"
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

    $0 -check_file "${tag_target_vars_BED}"
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

    rm -f -r ${ANONYMIZED_COORDINATES_DIR}
    mkdir ${ANONYMIZED_COORDINATES_DIR}

    chr_id=`sort -n -k2,2 ${tag_vars_BED} | head -n 1 | cut -f1`
    echo "Using chr_id=${chr_id}.."
    genetic_map_file=${PROXYTYPER_GENETIC_MAPS_DIR}/${chr_id}.map

    if [[ ! -f ${genetic_map_file} ]]
    then
        echo "Could not find the genetic map file \"${genetic_map_file}\""
        exit 1
    fi

    ${PROXYTYPER_EXEC} -anonymize_tag_target_genetic_map_coords_per_query_ref_variants \
${tag_vars_BED} ${tag_target_vars_BED} \
${genetic_map_file} ${coord_anon_cM_noise_SD} ${l_anonym_chrom} \
${ANONYMIZED_COORDINATES_DIR}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${ANONYMIZED_COORDINATES_DIR}/mapping_coords.bed" \
"${ANONYMIZED_COORDINATES_DIR}/original_coords.bed" \
""${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map"" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ###########################################################################################
    # TODO:: Check cMs per megabase; this is important while running BEAGLE?
    ##`awk 'BEGIN{FS=" "}{print $3}' ${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map | tail -n 1`
    ###########################################################################################

    # This file is written in case BEAGLE is ran.    
    # Note that right now BEAGLE and ProxyTyper use the same type of genetic map format (4-column PLINK maps) and we simply copy it, which creates a redundant file, which is ok for now.
    #awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map > ${ANONYMIZED_COORDINATES_DIR}/${chr_id}_BEAGLE.map
    cp ${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map ${ANONYMIZED_COORDINATES_DIR}/${chr_id}_BEAGLE.map

    # Validate the we have 4 columns.
    map_check=`awk '{if(NF != 4){print 1;exit}}END{print 0}' ${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map`
    if [[ ${map_check} -ne 0 ]]
    then
        echo "Could not validate the formatting of ${ANONYMIZED_COORDINATES_DIR}/${chr_id}.map.."
        exit 1
    fi

    if [[ ! -f "${ANONYMIZED_COORDINATES_DIR}/${chr_id}_BEAGLE.map" ]]
    then
        echo "Could not find anonymized BEAGLE map @ \"${ANONYMIZED_COORDINATES_DIR}/${chr_id}_BEAGLE.map\" (${LINE})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################

if [[ "${cmd_option}" == "-anonymize_coordinates" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3
 
    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi
 
    ${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs \
${panel_prefix} ${panel_prefix}_subjects.list \
${ANONYMIZED_COORDINATES_DIR}/original_coords.bed \
${ANONYMIZED_COORDINATES_DIR}/mapping_coords.bed \
${n_threads} \
${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

####################################################################################
if [[ "${cmd_option}" == "-deanonymize_coordinates" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3
 
    $0 -check_files "${panel_prefix}_genotypes.matrix.gz" \
"${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed"

    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi
 
    ${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs \
${panel_prefix} ${panel_prefix}_subjects.list \
${ANONYMIZED_COORDINATES_DIR}/mapping_coords.bed \
${ANONYMIZED_COORDINATES_DIR}/original_coords.bed \
${n_threads} \
${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

######################################################################

if [[ "${cmd_option}" == "-generate_tag_augmenter_model" ]]
then
    #-resample_panel [Panel prefix]
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Tag variants BED] [All panel variants BED]"
        exit 1
    fi

    tag_vars_BED=$2
    all_panel_vars_BED=$3

    $0 -check_files ${tag_vars_BED} ${all_panel_vars_BED}
    if [[ $? -ne 0 ]]
    then
        exit 1
    fi

    if [[ ${n_tag_augments} -eq 0 ]]
    then
        echo "Make sure n_tag_augments parameter is defined in the config file (i.e., Define, for example, \"n_tag_augments=2\" in PROXYTYPER.ini)."
        exit 1
    fi
    
    rm -f ${TAG_AUGMENTER_PREFIX}*

    tag_augment_iters=`seq 1 ${n_tag_augments}`
    for aug_i in ${tag_augment_iters[@]}
    do
        echo "Creating the tag augmentation mapper @ iteration ${aug_i}"
        # wc ${tag_vars_BED} ${all_panel_vars_BED}
        
        cur_iter_prefix=${TAG_AUGMENTER_PREFIX}_${aug_i}
        ${PROXYTYPER_EXEC} -build_augmenting_variants_mapper ${tag_vars_BED} ${all_panel_vars_BED} ${tag_augmenter_probability} ${cur_iter_prefix}
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        # echo "Before check_files"
        #         wc ${tag_vars_BED} ${all_panel_vars_BED}

        $0 -check_files ${cur_iter_prefix}_mapped.bed ${cur_iter_prefix}_original.bed
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        # echo "Before concat"
        #         wc ${tag_vars_BED} ${all_panel_vars_BED}

        # This is the new set of tags that will be used in the next iteration.
        # The pooled set of variants must be unique and sorted.
        # The new augmented set of variants is basically the concated list of original variants and the randomly positioned mapped variants.
        # We need to keep track of this in both whole set and the tags.
        # Note that some of the focus variants may not have augmented variants (due to randomness of augmentation).
        cat ${tag_vars_BED} ${cur_iter_prefix}_mapped.bed | sort -n -k2,2 > temp_tags_${aug_i}.bed
        cat ${all_panel_vars_BED} ${cur_iter_prefix}_mapped.bed | sort -n -k2,2 > temp_all_vars_${aug_i}.bed

        # echo "Before reassigning"
        #         wc ${tag_vars_BED} ${all_panel_vars_BED}

        tag_vars_BED=temp_tags_${aug_i}.bed
        all_panel_vars_BED=temp_all_vars_${aug_i}.bed

        $0 -check_BED_sortedness ${tag_vars_BED}
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        $0 -check_BED_sortedness ${all_panel_vars_BED}
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        # echo "After reassigning"
        #         wc ${tag_vars_BED} ${all_panel_vars_BED}

    done # aug_i loop.

    # Clean the temp files.
    rm -f temp_tags_* temp_all_vars_*

    # exit 1

    # Calculate the elapsed time, return.
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

##############################################################################################################################
if [[ "${cmd_option}" == "-copy_panel" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3

    echo "Copying panel: ${panel_prefix} >>==>> ${output_prefix}"

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    echo "Copying.."
    cp "${panel_prefix}_subjects.list" "${output_prefix}_subjects.list"
    cp "${panel_prefix}_variants.bed" "${output_prefix}_variants.bed"
    cp "${panel_prefix}_genotypes.matrix.gz" "${output_prefix}_genotypes.matrix.gz"

    echo "Checking copied panel.."
    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    exit $?
fi

##############################################################################################################################
if [[ "${cmd_option}" == "-intersect" ]]
then
    if [[ $# -ne 5 ]]
    then
        echo "USAGE: $0 $1 [Source BED] [Destination BED] [Region selector (Reg1/Reg2/Reg12)] [Overlap BED]"
        exit 1
    fi

    src_BED=$2
    dest_BED=$3
    output_selector=$4
    output_BED=$5

    if [[ ${output_selector} != "Reg1" && ${output_selector} != "Reg2" && ${output_selector} != "Reg12" ]]
    then
        echo "Output selector must be one of \"Reg1\", \"Reg2\", \"Reg12\" (${LINENO})"
        exit 1
    fi

    $0 -check_files "${src_BED}" "${dest_BED}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -intersect ${src_BED} ${dest_BED} no no ${output_selector}
    if [[ $? -ne 0 ]]
    then
        echo "Intersect failed (${LINENO})"
        exit 1
    fi

    # Copy the intersected file.
    cp intersected.bed ${output_BED}

    exit $?
fi

##############################################################################################################################
if [[ "${cmd_option}" == "-exclude" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Source BED] [Exclusion regions BED] [Output BED]"
        exit 1
    fi

    src_BED=$2
    dest_BED=$3
    output_BED=$4

    $0 -check_files "${src_BED}" "${dest_BED}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -exclude ${src_BED} ${dest_BED} n
    if [[ $? -ne 0 ]]
    then
        echo "Intersect failed (${LINENO})"
        exit 1
    fi

    # Copy the intersected file.
    cp excluded.bed ${output_BED}

    exit $?
fi 

##############################################################################################################################
if [[ "${cmd_option}" == "-augment_tag_variants" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    output_prefix=$3

    $0 -check_files "${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed" \
"${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    latest_output_prefix=${output_prefix}
    tag_augment_iters=`seq 1 ${n_tag_augments}`
    for aug_i in ${tag_augment_iters[@]}
    do
        cur_iter_augmenter_prefix=${TAG_AUGMENTER_PREFIX}_${aug_i}
        cur_iter_output_prefix=${output_prefix}_${aug_i}

        # Make sure the mapping coordinates exist.
        $0 -check_files ${cur_iter_augmenter_prefix}_mapped.bed ${cur_iter_augmenter_prefix}_original.bed
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        ${PROXYTYPER_EXEC} -augment_tag_variants ${cur_iter_augmenter_prefix} ${panel_prefix} ${panel_prefix}_subjects.list ${tag_genotype_augmentation_type} ${cur_iter_output_prefix}

        $0 -check_files "${cur_iter_output_prefix}_subjects.list" \
"${cur_iter_output_prefix}_variants.bed" \
"${cur_iter_output_prefix}_genotypes.matrix.gz" 
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi

        $0 -check_panel_sortedness ${cur_iter_output_prefix}
        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO}): ${cur_iter_output_prefix} is not sorted"
            exit 1
        fi

        # Update the panel prefix.
        panel_prefix=${cur_iter_output_prefix}

        # Update the latest output prefix.
        latest_output_prefix=${cur_iter_output_prefix}
    done

    # Make sure output prefix does not exist.
    if [[ -f ${output_prefix}_genotypes.matrix.gz ]]
    then
        echo "Cannot overwrite augmented output panel to: ${output_prefix}"
        exit 1
    fi

    if [[ -f ${output_prefix}_subjects.list ]]
    then
        echo "Cannot overwrite augmented output panel to: ${output_prefix}"
        exit 1
    fi

    if [[ -f ${output_prefix}_variants.bed ]]
    then
        echo "Cannot overwrite augmented output panel to: ${output_prefix}"
        exit 1
    fi

    # Copy the latest augmented panel.
    $0 -copy_panel ${latest_output_prefix} ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO}): Could not copy ${latest_output_prefix} ${output_prefix}"
        exit 1
    fi

    # Calculate the elapsed time, return.
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

############################################################################################################

if [[ "${cmd_option}" == "-decompose_variants" ]]
then
    if [[ $# -ne 4 ]]
    then
        echo "USAGE: $0 $1 [Tag variants BED] [Tag/Target Reference Panel prefix] [Output prefix]"
        exit 1
    fi

    tag_vars_BED=$2
    ref_tag_target_panel_prefix=$3
    output_prefix=$4
 
    $0 -check_files "${tag_vars_BED}" \
"${ref_tag_target_panel_prefix}_genotypes.matrix.gz" \
"${ref_tag_target_panel_prefix}_subjects.list" \
"${ref_tag_target_panel_prefix}_variants.bed"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    # TODO::Ensure full overlap between tags and ref panel; this is tested in decompose option.
    ${PROXYTYPER_EXEC} -intersect ${ref_tag_target_panel_prefix}_variants.bed ${tag_vars_BED} no no Reg2
    n_tags=`wc -l ${tag_vars_BED} | awk {'print $1'}`
    n_overlaps=`wc -l intersected.bed | awk {'print $1'}`

    if [[ "${n_tags}" -ne "${n_overlaps}" ]]
    then
        echo "# of tags in the reference panel is not same as the tag variant region. Make sure reference panel contains all tags.."
        exit 1
    fi

    ${PROXYTYPER_EXEC} -simple_decompose_untyped_variants_multithreaded ${tag_vars_BED} \
    ${ref_tag_target_panel_prefix} \
    ${ref_tag_target_panel_prefix}_subjects.list \
    ${variant_decomp_min_AAF} \
    ${shuffle_decomp_vars} \
    ${n_threads} \
    ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz" \
"${output_prefix}.interval" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

##########################################################################################################################################
if [[ "${cmd_option}" == "-recompose_BEAGLE_VCF_variants" ]]
then
    if [[ $# -ne 5 ]]
    then
        echo "USAGE: $0 $1 [BEAGLE VCF file] [BEAGLE imputed tag sample subject id's] [Untyped variant decomposition interval file] [Output prefix]"
        exit 1
    fi

    BEAGLE_VCF=$2
    beagle_imputed_panel_subject_list=$3
    untyped_decomposing_interval=$4
    output_prefix=$5

    $0 -check_files "${untyped_decomposing_interval}" \
"${BEAGLE_VCF}" \
"${beagle_imputed_panel_subject_list}"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # This is low memory, multithreaded slightly lower accuracy version:
    if [[ ${USE_MT_RECOMPOSITION} == 1 ]]
    then
        n_combine_threads=40
        echo "Using ${n_combine_threads}-threaded recomposition of the resampled query panel untyped variant genotypes."
        ${PROXYTYPER_EXEC} -combine_BEAGLE_imputed_decomposed_genotype_probabilities_multithreaded \
${BEAGLE_VCF} \
${untyped_decomposing_interval} ${beagle_imputed_panel_subject_list} ${n_combine_threads} ${output_prefix}

        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi
    fi

    if [[ ${USE_MT_RECOMPOSITION} == 0 ]]
    then
        echo "Using SINGLE-threaded recomposition of the resampled query panel untyped variant genotypes."
        ${PROXYTYPER_EXEC} -combine_BEAGLE_imputed_decomposed_genotype_probabilities \
${BEAGLE_VCF} \
${untyped_decomposing_interval} ${beagle_imputed_panel_subject_list} ${output_prefix}

        if [[ $? -ne 0 ]]
        then
            echo "Sanity check failed (${LINENO})"
            exit 1
        fi
    fi

    $0 -check_files "${output_prefix}_subjects.list" \
"${output_prefix}_variants.bed" \
"${output_prefix}_genotypes.matrix.gz"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    $0 -check_panel_sortedness ${output_prefix}
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

##########################################################################################################################################
if [[ "${cmd_option}" == "-run_BEAGLE" ]]
then
    if [[ $# -ne 7 ]]
    then
        echo "USAGE: $0 $1 [Query typed variant panel prefix] [Reference typed/untyped variant panel prefix] [BEAGLE map file, e.g., ${ANONYMIZED_COORDINATES_DIR}/20_BEAGLE.map] [Impute targets? (true/false)] [Phased query tag? (0/1)] [BEAGLE directory]" >&2
        exit 1
    fi

    query_panel_prefix=$2
    ref_panel_prefix=$3
    beagle_genetic_map_file=$4
    do_impute=$5
    USE_PHASED_SITE1_GENO=$6
    BEAGLE_dir=$7

    # Check parameter validity.
    if [[ "${do_impute}" != "true" ]]
    then
        if [[ "${do_impute}" != "false" ]]
        then
            echo "do_impute parameter must be true or false (${LINENO})"
            exit 1
        fi
    fi

    # If phasing-only is requested, check if options are inconsistent, first.
    if [[ "${do_impute}" == "false" ]]
    then
        echo "PERFORMING PHASING ONLY WITH BEAGLE @ test=${query_panel_prefix}, reference=${ref_panel_prefix}"
        if [[ ${USE_PHASED_SITE1_GENO} -eq 1 ]]
        then
            echo "Cannot save phased genotypes for test panel in a phasing-only BEAGLE run ($LINENO)"
            exit 1
        fi
    fi

    $0 -check_files "${beagle_genetic_map_file}" \
"${query_panel_prefix}_genotypes.matrix.gz" \
"${query_panel_prefix}_subjects.list" \
"${query_panel_prefix}_variants.bed" \
"${ref_panel_prefix}_genotypes.matrix.gz" \
"${ref_panel_prefix}_subjects.list" \
"${ref_panel_prefix}_variants.bed"
    res=$?
    if [[ ${res} -ne 0 ]]
    then
        exit 1
    fi

    #USE_PHASED_SITE1_GENO=0
    rm -f -r ${BEAGLE_dir}
    mkdir ${BEAGLE_dir}
    $0 -setup_BEAGLE_files ${query_panel_prefix} ${query_panel_prefix}_subjects.list \
${ref_panel_prefix} ${ref_panel_prefix}_subjects.list ${USE_PHASED_SITE1_GENO} ${BEAGLE_dir}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "BEAGLE run failed (${LINENO})."
        exit 1
    fi

    $0 -check_files "${beagle_genetic_map_file}" \
"${BEAGLE_dir}/ref_option.gens.vcf.gz" \
"${BEAGLE_dir}/gt_option.gens.vcf.gz"
    res=$?
    if [[ ${res} -ne 0 ]]
    then
        exit 1
    fi

    $0 -do_run_BEAGLE ${beagle_genetic_map_file} ${do_impute} ${BEAGLE_dir}

    if [[ ! -f "${BEAGLE_dir}/imputed.op.vcf.gz" ]]
    then
        echo "BEAGLE running error!"
        exit 1
    fi

    # Calculate the elapsed time: Note that this option includes previous times.
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

    exit $?
fi

##########################################################################################################################################
if [[ ${cmd_option} == "-haplo2geno" ]]
then
	if [[ $# -lt 3 ]]
	then
		echo "$0 $1 [HC Panel prefix] [GC Output prefix]" >&2
		exit 1
	fi

    HC_panel_prefix=$2
    GC_panel_prefix=$3

    $0 -check_files "${HC_panel_prefix}_genotypes.matrix.gz" \
"${HC_panel_prefix}_subjects.list" \
"${HC_panel_prefix}_variants.bed" 
    res=$?
    if [[ $res -ne 0 ]]
    then
        exit 1
    fi

    ${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${HC_panel_prefix} ${HC_panel_prefix}_subjects.list ${GC_panel_prefix}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Extraction failed (${LINENO})."
        exit 1
    fi

    $0 -check_files "${GC_panel_prefix}_subjects.list" \
"${GC_panel_prefix}_variants.bed" \
"${GC_panel_prefix}_genotypes.matrix.gz" 
    res=$?
    if [[ $res -ne 0 ]]
    then
        exit 1
    fi 

    exit $?
fi

##########################################################################################################################################
if [[ ${cmd_option} == "-summarize_accuracy" ]]
then
    if [[ $# -ne 4 ]]
	then
		echo "$0 $1 [prxytypr accuracy list file] [Min MAF] [Max MAF]"
		exit 1
	fi

    prxytypr_accs_list_file=$2
    min_MAF=$3
    max_MAF=$4

    $0 -check_files ${prxytypr_accs_list_file}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "File not found: ${prxytypr_accs_list_file} (${LINENO})."
        exit 1
    fi

    # Check the column number:
    n_cols=`awk 'BEGIN{FS="\t"}{print NF}' ${prxytypr_accs_list_file} | sort -n | tail -n 1`
    if [[ ${n_cols} -ne 9 ]]
    then
        echo "${prxytypr_accs_list_file} is not as expected (9 columns)"
        exit 1
    fi

    awk -v min_MAF=${min_MAF} -v max_MAF=${max_MAF} 'BEGIN{if(min_MAF >= max_MAF){print min_MAF" >= "max_MAF", cannot process this range..";exit 1}}{split($4, arr, "_");var_MAF=arr[4];if(var_MAF>0.5){var_MAF=1-arr[4]};if(var_MAF>min_MAF && var_MAF<max_MAF){tot+=$5;n_tot+=1;}}END{if(n_tot>0){print "R2 @ "max_MAF">MAF>"min_MAF"="tot/n_tot" ["n_tot" variants]"}}' ${prxytypr_accs_list_file}

    exit 0
fi

##########################################################################################################################################
if [[ ${cmd_option} == "-genotype_R2" ]]
then
	if [[ $# -ne 5 ]]
	then
		echo "$0 $1 [Target variants BED file] [Imputed genotypes matbed] [Known genotypes matbed] [Output file]"
		exit 1
	fi

    target_vars_BED_fp=$2
    imputed_panel_prefix=$3
    known_panel_prefix=$4
    R2_file=$5

	exec_check=`type -P ${PROXYTYPER_EXEC}`
	if [[ ${exec_check} == "" ]]
	then
		echo "Could not find ProxyTyper @ ${PROXYTYPER_EXEC}"
		exit 1
	fi

    $0 -check_files ${target_vars_BED_fp} 
    if [[ $? -ne 0 ]]
    then
        echo "Target variants file not found @ ${target_vars_BED_fp} (${LINENO})."
        exit 1
    fi

    ########################################################################################################################
    rm -f temp1_*
    $0 -subset_panel_variants ${imputed_panel_prefix} ${target_vars_BED_fp} ${imputed_panel_prefix}_targets_only
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Extraction failed (${LINENO})."
        exit 1
    fi
    imputed_panel_prefix=${imputed_panel_prefix}_targets_only

    rm -f temp2_*
    $0 -subset_panel_variants ${known_panel_prefix} ${target_vars_BED_fp} ${known_panel_prefix}_targets_only
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Extraction failed (${LINENO})."
        exit 1
    fi
    known_panel_prefix=${known_panel_prefix}_targets_only

    ########################################################################################################################
    $0 -haplo2geno ${imputed_panel_prefix} ${imputed_panel_prefix}_imputed_GC
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi
    imputed_panel_prefix=${imputed_panel_prefix}_imputed_GC

    $0 -haplo2geno ${known_panel_prefix} ${known_panel_prefix}_imputed_GC
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi
    known_panel_prefix=${known_panel_prefix}_imputed_GC

    ########################################################################################################################
	${PROXYTYPER_EXEC} -get_R2_per_imputed_genotypes ${imputed_panel_prefix} ${imputed_panel_prefix}_subjects.list ${known_panel_prefix} ${known_panel_prefix}_subjects.list
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

	if [[ ! -f "R2_stats.txt" ]]
	then
		echo "Could not find the R2 stats @ R2_stats.txt"
		exit 1
	fi

	mv R2_stats.txt ${R2_file}
	exit $?
fi

if [[ "${cmd_option}" == "-delete_panel" ]]
then
    if [[ $# -ne 2 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix]"
        exit 1
    fi

    panel_prefix=$2

    $0 -check_files "${panel_prefix}_genotypes.matrix.gz" \
"${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    # Delete panel files, most likely for saving space..
    echo "Deleting ${panel_prefix}.."
    rm -f "${panel_prefix}_genotypes.matrix.gz" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    rm -f "${panel_prefix}_subjects.list" 
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    rm -f "${panel_prefix}_variants.bed"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    exit $?
fi

if [[ "${cmd_option}" == "-sort_panel" ]]
then
    if [[ $# -ne 3 ]]
    then
        echo "USAGE: $0 $1 [Panel prefix] [Output prefix]"
        exit 1
    fi

    panel_prefix=$2
    sorted_panel_prefix=$3

    $0 -check_files "${panel_prefix}_genotypes.matrix.gz" \
"${panel_prefix}_subjects.list" \
"${panel_prefix}_variants.bed"
    if [[ $? -ne 0 ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi
    
    ${PROXYTYPER_EXEC} -sort_genosignal_matrix ${panel_prefix} ${panel_prefix}_subjects.list ${sorted_panel_prefix}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    if [[ ! -f "${sorted_panel_prefix}_subjects.list" ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -f "${sorted_panel_prefix}_variants.bed" ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    if [[ ! -f "${sorted_panel_prefix}_genotypes.matrix.gz" ]]
    then
        echo "Sanity check failed (${LINENO})"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"
    exit $?
fi

if [[ ${cmd_option} == "-setup_BEAGLE_files" ]]
then
	if [[ $# -lt 7 ]]
	then
		echo "$0 $1 [Haplocoded query tag genotypes matbed] [Query sample ids] [Haplocoded reference tag+target genotype matbed] [Reference sample ids] [Use phased(1)/unphased(0) genotypes] [Output directory]"
		exit 1
	fi
	
	query_tag_haplo_matbed=$2
	query_sample_ids=$3
	ref_tag_target_haplo_matbed=$4
	ref_sample_ids=$5
	SAVE_PHASED_FLAG=$6
	op_dir=$7
	
	BEAGLE_INPUT_DATA_DIR=${op_dir}
	
	if [[ ! -d ${BEAGLE_INPUT_DATA_DIR} ]]
	then
		echo "Could not find the output directory for BEAGLE @ \"${BEAGLE_INPUT_DATA_DIR}\""
		exit 1
	fi

    $0 -check_files ${ref_tag_target_haplo_matbed}_subjects.list \
${ref_tag_target_haplo_matbed}_genotypes.matrix.gz \
${ref_tag_target_haplo_matbed}_variants.bed \
${query_tag_haplo_matbed}_subjects.list \
${query_tag_haplo_matbed}_genotypes.matrix.gz \
${query_tag_haplo_matbed}_variants.bed 
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    ######
    # Validate the strict subsetting of typed in reference typed/untyped.
    $0 -intersect ${ref_tag_target_haplo_matbed}_variants.bed ${query_tag_haplo_matbed}_variants.bed Reg12 temp_query_ref_overlap.bed
    query_ref_overlap_size=`wc -l temp_query_ref_overlap.bed | awk {'print $1'}`
    query_size=`wc -l ${query_tag_haplo_matbed}_variants.bed | awk {'print $1'}`
    if [[ ${query_size} -ne ${query_ref_overlap_size} ]]
    then
        echo "Query does not perfectly overlap with reference.."
        exit 1
    fi
    ID_CHECK_PASS=`awk 'BEGIN{ID_CHECK=1}{if($5 != $11){ID_CHECK=0}}END{print ID_CHECK}' temp_query_ref_overlap.bed`
    if [[ ${ID_CHECK_PASS} -ne 1 ]]
    then
        echo "ID's do not match between typed and untyped panels. This could be due to mixing/merging panels or deanonymization inconsistency. This will create inconsistent VCFs.. (${LINENO}: ${ref_tag_target_haplo_matbed}, ${query_tag_haplo_matbed})"
        exit 1
    fi
    ######

    echo "Writing Beagle input files."
    rm -f gt_option_*.vcf
    rm -f input_option_*.vcf
    rm -f gt_option_*.vcf.gz
    rm -f input_option_*.vcf.gz
    ${PROXYTYPER_EXEC} -extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_multithreaded ${ref_tag_target_haplo_matbed} ${ref_sample_ids} ${query_tag_haplo_matbed} ${query_sample_ids} ${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf.gz ${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf.gz ${SAVE_PHASED_FLAG} ${n_threads}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    $0 -check_files ${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf.gz \
${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf.gz
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

	exit 0
fi


if [[ ${cmd_option} == "-do_run_BEAGLE" ]]
then
	if [[ $# -lt 4 ]]
	then
		echo "USAGE: $0 $1 [Genetic map file] [Impute (true/false)] [BEAGLE data directory]"
		exit 1
	fi

	genetic_map_file=$2
    do_impute=$3
	BEAGLE_INPUT_DATA_DIR=$4

    # Check parameter validity.
    if [[ "${do_impute}" != "true" ]]
    then
        if [[ "${do_impute}" != "false" ]]
        then
            echo "do_impute parameter must be true or false (${LINENO})"
            exit 1
        fi
    fi

    if [[ "${do_impute}" == "false" ]]
    then
        echo "**Doing phasing only with BEAGLE @ ${BEAGLE_INPUT_DATA_DIR}**"
    fi
	
    $0 -check_files ${genetic_map_file} 
    if [[ $res -ne 0 ]]
    then
        echo "Could not find PLINK formatted genetic map file @ ${genetic_map_file} -- You can download it, for example, from https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip (${LINENO})."
        exit 1
    fi

	$0 -check_files ${BEAGLE_JAR} 
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Could not find BEAGLE JAR file @ ${BEAGLE_JAR} -- You can download it from https://faculty.washington.edu/browning/beagle/${BEAGLE_JAR} (${LINENO})."
        echo "Make sure that java version is compliant with the BEAGLE!!!"        
        exit 1
    fi

	# Run beagle.
	java -${JAVA_MEM_FLAG} -jar ${BEAGLE_JAR} nthreads=${n_threads} map=${genetic_map_file} impute=${do_impute} ref=${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf.gz gt=${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf.gz gp=true ap=true out=${BEAGLE_INPUT_DATA_DIR}/imputed.op
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "BEAGLE run failed (${LINENO})."
        exit 1
    fi

    if [[ ! -f "${BEAGLE_INPUT_DATA_DIR}/imputed.op.vcf.gz" ]]
    then
        echo "Sanity check failed @ ${LINENO}"
        exit 1
    fi

    # Calculate the elapsed time
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    date_str=`date +%M/%D/%Y-%H:%M:%S`
    echo "[${date_str}]: \"$0 $@\" (${elapsed_time} seconds)" | tee -a "${WALLTIME_LOG}"

	exit $?
fi

# This is the main option to export a VCF from ProxyTyper's matrix formatted datasets.
if [[ ${cmd_option} == "-export_VCF" ]]
then
	if [[ $# -lt 5 ]]
	then
		echo "USAGE: $0 $1 [Panel identifier] [Assembly ID (e.g., hg19)] [Phased output? (0/1)] [Output VCF.gz file path]"
		exit 1
	fi

	panel_id=$2
    assembly_id=$3
    use_phased_output=$4
	gzip_op_VCF_file=$5
	
	if [[ ! -f ${panel_id}_subjects.list ]]
	then
		echo "Could not find ${panel_id} (${LINENO})"
		exit 1
    fi
	
    if [[ ! -f ${panel_id}_genotypes.matrix.gz ]]
	then
        echo "Could not find ${panel_id} (${LINENO})"
		exit 1
	fi 

    if [[ ! -f ${panel_id}_variants.bed ]]
	then
        echo "Could not find ${panel_id} (${LINENO})"
		exit 1
	fi 

    fetchChromSizes_EXEC=fetchChromSizes
    exec_check=`type -P ${fetchChromSizes_EXEC}`
    if [[ ${exec_check} == "" ]]
    then
        echo "Could not find fetchChromSizes, attempting to download from UCSC Genome Browser @ https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes ..."
        wget -c https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
        chmod 755 fetchChromSizes
        fetchChromSizes_EXEC=./fetchChromSizes_EXEC

        echo "Checking again.."
        exec_check=`type -P ${fetchChromSizes_EXEC}`
        if [[ ${exec_check} == "" ]]
        then
            echo "Could not find fetchChromSizes, make sure it is downloaded from UCSC Genome Browser @ https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes"
            exit 1
        fi

        echo "Success at installing fetchChromSizes.."
    fi

    ${fetchChromSizes_EXEC} ${assembly_id} > temp_${assembly_id}.list

    ${PROXYTYPER_EXEC} -get_me_VCF_header_per_chrom_lengths temp_${assembly_id}.list ${assembly_id} temp_VCF_header.list

    $0 -check_files temp_${assembly_id}.list
    if [[ $? -ne 0 ]]
    then
        echo "Failed to retrieve chromosome sizes for ${assembly_id}"
        exit 1
    fi

    ${PROXYTYPER_EXEC} -convert_genotype_signal_regions_2_VCF_w_header_multithreaded ${panel_id} ${panel_id}_subjects.list temp_VCF_header.list ${use_phased_output} ${n_threads} ${gzip_op_VCF_file}
    res=$?
    if [[ $res -ne 0 ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    if [[ ! -f "${gzip_op_VCF_file}" ]]
    then
        echo "Conversion failed (${LINENO})."
        exit 1
    fi

    exit $?
fi

# Give an error when we dont understand the option.
echo "UNKNOWN OPTION: ${cmd_option}"
exit 1