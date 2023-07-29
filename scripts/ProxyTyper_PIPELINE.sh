#!/bin/bash

#########################################################################################################
# This script uses backend script/code to run experiments on multiple chromosomes, etc.
#########################################################################################################

if [[ $# -lt 1 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
	Data Extraction:
		-extract_SITE12_genotypes_from_ALL_panel
	Resampling:
		-do_Resample_Site1_tag_Genotypes
		-do_Resample_Site2_tag_target_Genotypes
	Imputation:
		-do_central_imputes
		-do_impute_site1_per_proxy_protocol
	Leakage:
		-do_get_model1_proxy_2_model2_proxy_hapfreq_corrs
		-do_HapFreq_HMM_Attack
		-compare_Proxy_2_Original_HapFreq_stats
		-calculate_ProxyVar_2_OriginalVar_Correlations
		-compute_proxy2cleartext_linking_stats
		-calculate_windowed_Homer_t_attack_statistics
		-calculate_Homer_t_attack_statistics
		-calculate_Sankararaman_LRT_attack_statistics
		-hapfreq_signature_attack
	Diagnostics:
		-calculate_proxy_LD_patterns
		-get_untyped_variant_LD_stats
		-generate_PCA_proxy_original_data
	Parameter Testing:
		-grid_impute_per_proxized_tags_original_panel"

	exit 1
fi

cmd_option=$1

ASSEMBLY_ID=hg19

#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1
#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/KG_Synthetic_Array/1kG_GENOTYPE/AFR_ALL
#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/KG_Synthetic_Array/1kG_GENOTYPE/1kG_2_per_pop_Site1

chr_ids=(19 20 21 22)

########################################################################################################################
# Following is the placeholder for the data per chromosome.
SITE1_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE1_tag_haplocoded.matbed.gz
SITE1_haplocoded_target_matbed=../../DATA/Pop_Specific_Sites/${chr_id}_SITE1_target_haplocoded.matbed.gz
SITE1_sample_list=${PER_CHROM_GENO_DATA_DIR}/site1_sample_ids.list

SITE2_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_tag_haplocoded.matbed.gz
SITE2_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_target_haplocoded.matbed.gz
SITE2_sample_list=${PER_CHROM_GENO_DATA_DIR}/site2_sample_ids.list

per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
#per_chrom_beagle_maps_dir=../../beagle/genetic_maps/plink.chr22.GRCh37.map
per_chrom_beagle_maps_dir=../../beagle/genetic_maps/
########################################################################################################################

# Use this option to randomly partition site1 and site2 panels from an ALL panel.
# This is useful when we did not preselect the panels or when more complex sampling is used.
if [[ ${cmd_option} == "-extract_SITE12_genotypes_from_ALL_panel" ]]
then
	SITE1_SAMPLE_SIZE=100

	## Extract site1 and site2 genotypes: Use AFR ALL data, single pop.
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/KG_Synthetic_Array/1kG_GENOTYPE/AFR_ALL
	ALL_sample_list=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps

	if [[ ! -d ${PER_CHROM_GENO_DATA_DIR} ]]
	then
		echo "Could not find \"${PER_CHROM_GENO_DATA_DIR}\""
		exit 1
	fi

	if [[ ! -f ${ALL_sample_list} ]]
	then
		echo "Could not find \"${ALL_sample_list}\""
		exit 1
	fi

	if [[ ! -d ${per_chrom_beagle_maps_dir} ]]
	then
		echo "Could not find \"${per_chrom_beagle_maps_dir}\""
		exit 1
	fi

	echo "Extracting "

	RUN_CENTRAL_IMPUTES=0
	GET_SITE12_SAMPLES=1
	if [[ ${GET_SITE12_SAMPLES} == 1 ]]
	then
		shuf ${ALL_sample_list} | head -n ${SITE1_SAMPLE_SIZE} > SITE1_sample.list
		grep -w -v -f SITE1_sample.list ${ALL_sample_list} > SITE2_sample.list
	fi

	SITE1_sample_list=SITE1_sample.list
	SITE2_sample_list=SITE2_sample.list

	chr_ids=(19)

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		ASSEMBLY_ID=hg19

		ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
		ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz
		ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

		########################################################################
		echo "Extracting the genotypes for the site1 and site 2 genotypes."
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_target_haplocoded.matbed.gz

		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_target_haplocoded.matbed.gz

		SITE1_haplocoded_tag_matbed=${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${chr_id}_SITE1_target_haplocoded.matbed.gz

		SITE2_haplocoded_tag_matbed=${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${chr_id}_SITE2_target_haplocoded.matbed.gz

		# Do central imputation on the current subjects.
		if [[ ${RUN_CENTRAL_IMPUTES} == 1 ]]
		then
            ASSEMBLY_ID=hg19

            beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

            echo "Performing site2 imputation"
            # Calculate the pure site2-based imputation accuracy.
            ./ProxyTyper.sh -centralized_beagle_impute ${chr_id} ${ASSEMBLY_ID} ${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} ${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} ${beagle_genetic_map_file}
            cp -r central_impute_beagle_dir ${chr_id}_central_impute_beagle_dir

			exit 0
		fi
	done

	exit 0
fi

if [[ ${cmd_option} == "-get_untyped_variant_LD_stats" ]]
then
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=20
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_UNTYPED_LD_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	l_half_hap_win=6

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		#grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list
		#shuf site1_part1_samples.list > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_target_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_target_haplocoded.matbed.gz
		
			SITE1_haplocoded_tag_matbed=site1_part1_tag_haplocoded.matbed.gz
			SITE1_haplocoded_target_matbed=site1_part1_target_haplocoded.matbed.gz
			SITE1_sample_list=site1_part1_samples.list

			resampling_op_prefix=PART1_${chr_id}
			resampled_size=${n_site1_part_samples}
			N_e_frac=0.05
			allele_eps=0
			max_l_seg_n_bps=10000000
			max_l_seg_cM=0
			max_l_seg_nvars=0
			n_threads=20
			start_posn=0
			end_posn=250000000
			./ProxyTyper.sh -resample_tag_target_genotypes_w_length_cutoff ${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${resampling_op_prefix} 

			n_blocks_2_process=100
			ProxyTyper -get_untyped_variant_LD_statistics ${resampling_op_prefix}_resampled_tags.matbed.gz \
${resampling_op_prefix}_resampled_targets.matbed.gz \
${resampling_op_prefix}_resampled_sample_ids.list ${n_blocks_2_process} ${LEAKAGE_OP_DIR}/resampled_untyped_LD

			ProxyTyper -get_untyped_variant_LD_statistics site1_part2_tag_haplocoded.matbed.gz \
site1_part2_target_haplocoded.matbed.gz \
site1_part2_samples.list \
${n_blocks_2_process} ${LEAKAGE_OP_DIR}/part2_untyped_LD

			exit
		done
	done

	exit 0
fi

if [[ ${cmd_option} == "-do_Resample_Site1_tag_Genotypes" ]]
then
	PER_CHROM_GENO_DATA_DIR=../../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1/
	#POP_INFO_FILE=../../../DATA/1kG/population_info.list
	per_chrom_maps_dir=../../../DATA/beagle_genetic_maps/genetic_maps
	SITE1_sample_list=${PER_CHROM_GENO_DATA_DIR}/site1_sample_ids.list

	if [[ ! -d ${PER_CHROM_GENO_DATA_DIR} ]]
	then
		echo "Could not find \"${PER_CHROM_GENO_DATA_DIR}\""
		exit 1
	fi

	if [[ ! -f ${SITE1_sample_list} ]]
	then
		echo "Could not find \"${SITE1_sample_list}\""
		exit 1
	fi

	if [[ ! -d ${per_chrom_maps_dir} ]]
	then
		echo "Could not find \"${per_chrom_maps_dir}\""
		exit 1
	fi

	#tail -n +2 ${POP_INFO_FILE} | cut -f2,2 | sort -u > ALL_POPs.list

	#pops=`cat ALL_POPs.list`

	#chr_ids=(19 20 21 22)
	chr_ids=(19)

	postfix=`date +%d_%m_%y_%H_%M_%S`
	DATA_DIR=RESAMPLED_DATA_${postfix}
	rm -f -r ${DATA_DIR}
	mkdir ${DATA_DIR}

	for chr_id in ${chr_ids[@]}
	do
		echo "Resampling ${chr_id}"

		SITE1_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE1_target_haplocoded.matbed.gz
		SITE1_sample_list=${PER_CHROM_GENO_DATA_DIR}/site1_sample_ids.list

		op_prefix=${DATA_DIR}/SITE1_${chr_id}
		resampled_size=1000
		N_e_frac=0.125
		allele_eps=0
		max_l_seg_n_bps=10000000
		max_l_seg_cM=0
		max_l_seg_nvars=0
		n_threads=10
		start_posn=0
		end_posn=250000000
		./ProxyTyper.sh -resample_query_genotypes_w_length_cutoff ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${op_prefix} 

		cp ${op_prefix}_resampled_tags.matbed.gz ${chr_id}_SITE1_resampled_tag_haplocoded.matbed.gz

		cp ${op_prefix}_resampled_sample_ids.list ${chr_id}_SITE1_resampled_samples.list
	done #chr_id loop.

	exit 0
fi

if [[ ${cmd_option} == "-do_Resample_Site2_tag_target_Genotypes" ]]
then
	PER_CHROM_GENO_DATA_DIR=../../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1/
	#POP_INFO_FILE=../../../DATA/1kG/population_info.list
	per_chrom_maps_dir=../../../DATA/beagle_genetic_maps/genetic_maps
	SITE2_sample_list=${PER_CHROM_GENO_DATA_DIR}/site2_sample_ids.list

	if [[ ! -d ${PER_CHROM_GENO_DATA_DIR} ]]
	then
		echo "Could not find \"${PER_CHROM_GENO_DATA_DIR}\""
		exit 1
	fi

	if [[ ! -f ${SITE2_sample_list} ]]
	then
		echo "Could not find \"${SITE1_sample_list}\""
		exit 1
	fi

	if [[ ! -d ${per_chrom_maps_dir} ]]
	then
		echo "Could not find \"${per_chrom_maps_dir}\""
		exit 1
	fi

	#tail -n +2 ${POP_INFO_FILE} | cut -f2,2 | sort -u > ALL_POPs.list
	#pops=`cat ALL_POPs.list`

	#chr_ids=(19 20 21 22)
	chr_ids=(19)

	postfix=`date +%d_%m_%y_%H_%M_%S`
	DATA_DIR=RESAMPLED_DATA_${postfix}
	rm -f -r ${DATA_DIR}
	mkdir ${DATA_DIR}

	for chr_id in ${chr_ids[@]}
	do
		echo "Resampling ${chr_id}"

		SITE2_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_target_haplocoded.matbed.gz
		SITE2_sample_list=${PER_CHROM_GENO_DATA_DIR}/site2_sample_ids.list

		op_prefix=${DATA_DIR}/SITE2_${chr_id}
		resampled_size=2000
		N_e_frac=0.125
		allele_eps=0
		max_l_seg_n_bps=10000000
		max_l_seg_cM=0
		max_l_seg_nvars=0
		n_threads=10
		start_posn=0
		end_posn=250000000

		./ProxyTyper.sh -resample_tag_target_genotypes_w_length_cutoff ${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${op_prefix} 

		cp ${op_prefix}_resampled_tags.matbed.gz ${chr_id}_SITE2_resampled_tag_haplocoded.matbed.gz
		cp ${op_prefix}_resampled_targets.matbed.gz ${chr_id}_SITE2_resampled_target_haplocoded.matbed.gz

		cp ${op_prefix}_resampled_sample_ids.list ${chr_id}_SITE2_resampled_samples.list
	done #chr_id loop.

	exit 0
fi

if [[ ${cmd_option} == "-do_get_model1_proxy_2_model2_proxy_hapfreq_corrs" ]]
then
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=20
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_PROXY1_2_PROXY2_HAPFREQ_CORRELATION_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	l_half_hap_win=6

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list
		#shuf site1_part1_samples.list > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models for part1."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} model1_${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} model1_${permute_proxying_mapping_BED}

			########################################################################
			echo "Generating the proxizing models for part2 data."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} model2_${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} model2_${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz

			########################################################################
			echo "Proxizing part1 with model1."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list model1_${permute_proxying_mapping_BED} site1_part1_perm_prox_model1.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox_model1.matbed.gz site1_part1_samples.list model1_${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_part1_model1_proxized_haplocoded.matbed.gz

			echo "Proxizing part1 with model2."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list model2_${permute_proxying_mapping_BED} site1_part1_perm_prox_model2.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox_model2.matbed.gz site1_part1_samples.list model2_${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_part1_model2_proxized_haplocoded.matbed.gz

			# Calculate the correlations.
			ProxyTyper -get_proxization1_2_proxization_mapped_haplotype_frequency_correlation site1_part1_haplocoded.matbed.gz site1_part1_model1_proxized_haplocoded.matbed.gz site1_part1_model2_proxized_haplocoded.matbed.gz site1_part1_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_site1_part1_model1_2_model2_proxy_hapfreq_corrs_${cur_iter_i}

            ProxyTyper -get_proxization1_2_proxization_mapped_haplotype_frequency_correlation site1_part1_haplocoded.matbed.gz site1_part1_model1_proxized_haplocoded.matbed.gz site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_site1_part1_model1_2_orig_hapfreq_corrs_${cur_iter_i}

		done
	done

	exit 0
fi


if [[ ${cmd_option} == "-do_HapFreq_HMM_Attack" ]]
then
	filter_proxization_vicinity_size=6
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=15

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=1
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	#SITE1_SAMPLE_SIZE=2500
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_HAPFREQ_DP_ATTACK_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	l_proxy_vicinity_array=(3 4 6 8)
	l_attack_vicinity_array=(4 6 8 10)

	# This sets the matching state of the variants.
	VARIANT_POSITIONS_DECODED=0

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} > ${SITE1_sample_list}

		# Extract half of the samples in SITE1's list.
		#n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/6)'}`
		n_site1_part_samples=100

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		##################################################################
		# part1 is the target panel that is proxytyped.
		# part2 is the reference panel used for decoding.
		# part3 is the validation/control panel.
		##################################################################
		rm -f site1_part1_samples.list
		rm -f site1_part2_samples.list
		rm -f site1_part3_samples.list

		# Generate part1.
		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list

		# Generate part2.
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part3_samples.list

		# Generate part3.
		grep -w -v -f site1_part1_samples.list ${SITE1_sample_list} | grep -w -v -f site1_part3_samples.list > site1_part2_samples.list

		##############################################################################
		# Start tracing all window sizes.
		for cur_l_proxy_win in ${l_proxy_vicinity_array[@]}
		do
			for cur_l_attack_win in ${l_attack_vicinity_array[@]}
			do
				#if [[ ${cur_l_proxy_win} -gt ${cur_l_attack_win} ]]
				#then
				#	echo "Skipping: ${cur_l_proxy_win}, ${cur_l_attack_win}"
				#	continue
				#fi

				echo "Processing: Proxy_win_l: ${cur_l_proxy_win}, Attack_win_l: ${cur_l_attack_win}"

				filter_proxization_vicinity_size=${cur_l_proxy_win}

				for chr_id in ${chr_ids[@]}
				do
					# 100 meg window.
					l_win_in_bps=100000000

					# Process 10meg-110meg.
					var_posn_start_array=(`seq 10000000 ${l_win_in_bps} 10000000`)

					for var_start_posn in ${var_posn_start_array[@]}
					do
						var_end_posn=`awk -v l_win_in_bps=${l_win_in_bps} -v var_start_posn=${var_start_posn} 'BEGIN{print var_start_posn+l_win_in_bps}'`
						echo "Processing [${var_start_posn}-${var_end_posn}]"

						# Init output prefix and clean.
						op_prefix=${LEAKAGE_OP_DIR}/LiStephens_Viterbi_stats_${cur_iter_i}_${chr_id}_${var_start_posn}_${var_end_posn}_${cur_l_proxy_win}_${cur_l_attack_win}
						rm -f ${op_prefix}_*

						cp site1_part1_samples.list ${op_prefix}_site1_part1_samples.list
						cp site1_part2_samples.list ${op_prefix}_site1_part2_samples.list
						cp site1_part3_samples.list ${op_prefix}_site1_part3_samples.list

						########################################################################################################################	
						# Use the all tag/target data.
						ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
						ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

						echo "Processing chromosome ${chr_id}"

						variation_tools -dump_plain_geno_signal_regions ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} 1 ALL_var_regs.bed

						n_all_even_vars=`awk 'END{print NR-NR%2}' ALL_var_regs.bed`

						#shared_frac=0.1
						#n_shared_vars=`wc -l ALL_var_regs.bed | awk -v shared_frac=${shared_frac} {'print int($1*shared_frac/2)'}`

						if [[ ${VARIANT_POSITIONS_DECODED} == 0 ]]
						then
							echo "Separating ${n_all_even_vars} variants"
							sort -n -k2,2 ALL_var_regs.bed | head -n ${n_all_even_vars} | awk {'if(NR % 2 == 0){print $0}'} > site1_part1_regs.bed
							sort -n -k2,2 ALL_var_regs.bed | head -n ${n_all_even_vars} | awk {'if(NR % 2 == 1){print $0}'} > site1_part23_regs.bed

							cp site1_part1_regs.bed ${op_prefix}_site1_part1_regs.bed
							cp site1_part23_regs.bed ${op_prefix}_site1_part23_regs.bed
						fi

						if [[ ${VARIANT_POSITIONS_DECODED} == 1 ]]
						then
							echo "Using all variants"
							cp ALL_var_regs.bed site1_part1_regs.bed 
							cp ALL_var_regs.bed site1_part23_regs.bed 

							cp ALL_var_regs.bed ${op_prefix}_site1_part1_regs.bed
							cp ALL_var_regs.bed ${op_prefix}_site1_part23_regs.bed
						fi

						########################################################################
						echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
						variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

						# Resample the genotypes before extracting the target regions; we need this to estimate the allele accuracy and LRT/t statistics on the decoded genotypes.
						# This is necessary because when subregs are used for resampling, we cannot compare them with attacker's variant list or for doing LRT/t attack.
						upsampled_sample_size=${n_site1_part_samples}
						N_e=0.05
						eps_allele=0
						l_seg_in_bps=10000000
						l_seg_in_cMs=0
						l_seg_in_n_vars=0
						n_threads=${N_THREADS}
						start_pos=0
						end_pos=0
						ProxyTyper -resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb site1_part1_haplocoded.matbed.gz site1_part1_samples.list \
${per_chrom_maps_dir} ${upsampled_sample_size} ${N_e} ${eps_allele} \
${l_seg_in_bps} ${l_seg_in_cMs} ${l_seg_in_n_vars} \
${n_threads} ${start_pos} ${end_pos} ${op_prefix} >& ${op_prefix}_RESAMPLING.OP

						# Summarize the resampled segments.
						ProxyTyper -summarize_sampled_segments_per_resampled_haplotype_info ${op_prefix}_resampling_pattern_signal_var_regs.sigbed.gz \
${op_prefix}_resampled_sample_ids.list site1_part1_samples.list \
${op_prefix}_resampled_segments.bed

						# After this point, all decoded data are from resampled part 1 panel.
						resampled_site1_part1_matbed=${op_prefix}_resampled_tags.matbed.gz
						resampled_site1_part1_samples=${op_prefix}_resampled_sample_ids.list
				
						# Extract the subregions from the resampled panel.
						resampled_site1_part1_subreg_matbed=${op_prefix}_resampled_subreg_tags.matbed.gz
						variation_tools -extract_genotype_signals_per_region_list ${resampled_site1_part1_matbed} ${resampled_site1_part1_samples} site1_part1_regs.bed ${resampled_site1_part1_subreg_matbed}

						# Extract the subregion genotypes for part1.
						variation_tools -extract_genotype_signals_per_region_list site1_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part1_regs.bed site1_part1_haplocoded_subregs.matbed.gz

						# Extract the part2 genotypes.
						variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz
						variation_tools -extract_genotype_signals_per_region_list site1_part2_haplocoded.matbed.gz site1_part2_samples.list site1_part23_regs.bed site1_part2_haplocoded_subregs.matbed.gz

						# Extract the part3 genotypes.
						variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part3_samples.list site1_part3_haplocoded.matbed.gz
						variation_tools -extract_genotype_signals_per_region_list site1_part3_haplocoded.matbed.gz site1_part3_samples.list site1_part23_regs.bed site1_part3_haplocoded_subregs.matbed.gz

						cp site1_part1_haplocoded_subregs.matbed.gz ${op_prefix}_site1_part1_haplocoded_subregs.matbed.gz
						cp site1_part2_haplocoded_subregs.matbed.gz ${op_prefix}_site1_part2_haplocoded_subregs.matbed.gz
						cp site1_part3_haplocoded_subregs.matbed.gz ${op_prefix}_site1_part3_haplocoded_subregs.matbed.gz

						########################################################################
						echo "Generating the proxizing models for part1."
						# Generate a new model for the current chromosome.
						ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${resampled_site1_part1_subreg_matbed} ${resampled_site1_part1_samples} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${op_prefix}_${weight_params_file}

						# Generate the permutation proxization mapping regions.
						ProxyTyper -generate_permute_proxizing_parameters ${resampled_site1_part1_subreg_matbed} ${resampled_site1_part1_samples} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${op_prefix}_${permute_proxying_mapping_BED}

						########################################################################
						echo "Proxizing part1 genotypes."
						# Proxize the genotypes.
						#ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded_subregs.matbed.gz site1_part1_samples.list ${op_prefix}_${permute_proxying_mapping_BED} ${op_prefix}_site1_part1_perm_prox.matbed.gz
						#ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters ${op_prefix}_site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${op_prefix}_${weight_params_file} ${allele_err_eps} ${N_THREADS} ${op_prefix}_site1_proxized_part1_haplocoded_subregs.matbed.gz

						ProxyTyper -permute_proxize_genotype_signal_regions ${resampled_site1_part1_subreg_matbed} ${resampled_site1_part1_samples} ${op_prefix}_${permute_proxying_mapping_BED} ${op_prefix}_site1_part1_perm_prox.matbed.gz
						ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters ${op_prefix}_site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${op_prefix}_${weight_params_file} ${allele_err_eps} ${N_THREADS} ${op_prefix}_site1_proxized_part1_haplocoded_subregs.matbed.gz

						######################################################################
						# Extract variant start/end indices.
						#var_start_posn=40000000
						#var_end_posn=60000000

						var_start_i=`sort -n -k2,2 site1_part1_regs.bed | awk -v posn=${var_start_posn} '{if($2>=posn){print NR;exit}}'`
						var_end_i=`sort -n -k2,2 site1_part1_regs.bed | awk -v posn=${var_end_posn} '{if($2>=posn){print NR;exit}}'`

						echo "Decoding proxy haplotypes on variants [${var_start_i}-${var_end_i}] at [${var_start_posn}-${var_end_posn}]"
						######################################################################

						######################################################################
						# Decode the variants using histogram matching Viterbi.
						n_vic=${cur_l_attack_win}
						kmer_concordance_weight=`echo ${cur_l_attack_win} | awk {'print (2*$1+1)*8.8'}`
						N_e=1000000
						n_subj_2_decode=${n_site1_part_samples}
						n_threads=40
						ProxyTyper -decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_MT \
${resampled_site1_part1_subreg_matbed} ${op_prefix}_site1_proxized_part1_haplocoded_subregs.matbed.gz ${resampled_site1_part1_samples} \
site1_part2_haplocoded_subregs.matbed.gz site1_part2_samples.list \
${per_chrom_maps_dir} \
${var_start_i} ${var_end_i} ${n_vic} ${kmer_concordance_weight} ${N_e} ${n_subj_2_decode} ${n_threads} ${op_prefix}

						# At this point, we have the decoded genotypes for the resampled panel.

						# Extract the ref regions with imputed genotypes.
						n_max_ref_kmers=100000
						cat ${op_prefix}_per_per_sample_decoding_stats_thread_*.txt | awk -v n_max_ref_kmers=${n_max_ref_kmers} '{if($5<n_max_ref_kmers){print $0;}}' | cut -f3,4 | awk {'print $1"\t"$2-1"\t"$2'} | sort -u | sort -n -k2,2 > ${op_prefix}_decoded_focus_regs.bed

						n_dec_haplotypes=`awk -v n_subj_2_decode=${n_subj_2_decode} 'BEGIN{print 2*n_subj_2_decode}'`

						# Extract the genotypes for part1, part2, part3 for the decoded regions. Note that query may not originally have genotypes in these regions.
						variation_tools -extract_genotype_signals_per_region_list ${op_prefix}_decoded_ref_genotypes.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_decoded_focus_regs.bed ${op_prefix}_decoded_ref_genotypes_decoded_regs.matbed.gz

						# Extract the genotypes from among all resampled genotypes the target variant regions, which will be used for doing LRT and t attacks.
						variation_tools -extract_genotype_signals_per_region_list ${resampled_site1_part1_matbed} ${resampled_site1_part1_samples} ${op_prefix}_decoded_focus_regs.bed ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz

						# Extract the cleartext genotypes on the focus regions.
						variation_tools -extract_genotype_signals_per_region_list site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${op_prefix}_decoded_focus_regs.bed ${op_prefix}_site1_part1_decoded_regs.matbed.gz
						variation_tools -extract_genotype_signals_per_region_list site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${op_prefix}_decoded_focus_regs.bed ${op_prefix}_site1_part2_decoded_regs.matbed.gz
						variation_tools -extract_genotype_signals_per_region_list site1_part3_haplocoded.matbed.gz site1_part3_samples.list ${op_prefix}_decoded_focus_regs.bed ${op_prefix}_site1_part3_decoded_regs.matbed.gz

						# Here, we need to first extract the individuals whose genotypes are decoded, this is necessary to build the database AF panel correctly from these.
						head -n ${n_subj_2_decode} ${resampled_site1_part1_samples} > ${op_prefix}_decoded_samples.list
						variation_tools -extract_genotype_signals_per_subsample_list ${op_prefix}_decoded_ref_genotypes_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_decoded_samples.list ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz

						# At this stage, we are using the decoded resampled genotypes and subselected subject id's of resampled subjects.

						# Link the ${n_subj_2_decode} subjects to the matching and control references.
						# We are linking the decoded resampled-proxy panel to cleartext part1 and part2.
						l_corr_win=0
						ProxyTyper -calculate_proxy2clear_pairwise_distance_stats ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list 0 0 ${l_corr_win} ${op_prefix}_decoded_2_part1_linking
						ProxyTyper -calculate_proxy2clear_pairwise_distance_stats ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list 0 0 ${l_corr_win} ${op_prefix}_decoded_2_part2_linking

						# Get the allele errors for the matching known panel: the sample id's need to match for this to work: We compare the decoded genotypes on focus regions with the original resampled genotypes.
						# Note that this is used only for reporting the error rates, resampled cleartext panel is not available to the adversary.
						ProxyTyper -get_allele_error_per_decoded_panel_known_panel ${op_prefix}_decoded_ref_genotypes_decoded_regs.matbed.gz ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_decoded_errors.txt

						# Now, shuffle the subjects of resampled cleartext genotypes and calculate the allele error.
						shuf ${resampled_site1_part1_samples} > ${op_prefix}_rand_site1_part1_samples.list
						variation_tools -extract_genotype_signals_per_subsample_list ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_rand_site1_part1_samples.list ${op_prefix}_rand_site1_part1_decoded_regs.matbed.gz

						# This is done on all decoded subjects since we need them to match with original resampled subjects.
						ProxyTyper -get_allele_error_per_decoded_panel_known_panel ${op_prefix}_decoded_ref_genotypes_decoded_regs.matbed.gz ${op_prefix}_rand_site1_part1_decoded_regs.matbed.gz ${op_prefix}_rand_site1_part1_samples.list ${op_prefix}_shuf_decoded_errors.txt

						# Get the linking statistics.
						part2_linking_stats=`head -n ${n_dec_haplotypes} ${op_prefix}_decoded_2_part2_linking_per_hap_corr_stats.txt | awk '{tot+=$6}END{print tot/NR}'`
						part1_linking_stats=`head -n ${n_dec_haplotypes} ${op_prefix}_decoded_2_part1_linking_per_hap_corr_stats.txt | awk '{tot+=$6}END{print tot/NR}'`

						decoded_err_stat=`head -n ${n_dec_haplotypes} ${op_prefix}_decoded_errors.txt | awk '{tot+=$3/$4}END{print tot/NR}'`
						decoded_shuff_err_stat=`head -n ${n_dec_haplotypes} ${op_prefix}_shuf_decoded_errors.txt | awk '{tot+=$3/$4}END{print tot/NR}'`

						####################################################################################################################################
						####################################################################################################################################
						# Start doing LRT reidentifications.
						min_LRT_MAF=0.05
						min_LRT_var2var_dist=1000

						# Do LRT linking using decoded genotypes as query, part2 as reference, part1 and part3 as the target AF database.
						# note that this uses the whole decoded genotypes, if we generated a proxy-panel, this would not be applicable.
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_match_decoded_query
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_nomatch_decoded_query

						# The first haplotypes are the ones that are decoded.
						LRT_match_decoded=`head -n ${n_dec_haplotypes} ${op_prefix}_match_decoded_query_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`
						LRT_nomatch_decoded=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_decoded_query_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`

						# Do LRT linking using part1 and part3 as query, part2 as reference, decoded genotypes as target AF database.
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_match_clear_query
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_nomatch_clear_query	

						# The first haplotypes are the ones that are positive.
						LRT_match_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_match_clear_query_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`
						LRT_nomatch_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_clear_query_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`

						# Do LRT linking using part1 and part3 as query, part2 as reference, clear genotypes as target AF database.
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_match_all_clear
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_nomatch_all_clear

						# The first haplotypes are the ones that are positive.
						LRT_match_all_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_match_all_clear_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`
						LRT_nomatch_all_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_all_clear_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`

						# Do LRT linking using resampled panel as the mixture.
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_match_all_clear_resampled_part1
						ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${min_LRT_MAF} ${min_LRT_var2var_dist} ${op_prefix}_nomatch_all_clear_resampled_part1

						# The first haplotypes are the ones that are positive.
						LRT_match_all_clear_resampled_part1=`head -n ${n_dec_haplotypes} ${op_prefix}_match_all_clear_resampled_part1_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`
						LRT_nomatch_all_clear_resampled_part1=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_all_clear_resampled_part1_per_query_Sankararaman_LRT_stats.txt | awk '{tot+=$3}END{print tot/NR}'`

						####################################################################################################################################
						####################################################################################################################################
						# Start doing Homer reidentifications.

						# Do Homer linking using decoded genotypes as query, part2 as reference, part1 and part3 as the target AF database.
						# note that this uses the whole decoded genotypes, if we generated a proxy-panel, this would not be applicable.
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_match_decoded_query
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_nomatch_decoded_query

						# The first haplotypes are the ones that are decoded.
						Homer_match_decoded=`head -n ${n_dec_haplotypes} ${op_prefix}_match_decoded_query_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`
						Homer_nomatch_decoded=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_decoded_query_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`

						# Do LRT linking using part1 and part3 as query, part2 as reference, decoded genotypes as target AF database.
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_match_clear_query
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_decoded_ref_genotypes_decoded_regs_decoded_subjects.matbed.gz ${op_prefix}_decoded_samples.list ${op_prefix}_nomatch_clear_query	

						# The first haplotypes are the ones that are positive.
						Homer_match_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_match_clear_query_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`
						Homer_nomatch_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_clear_query_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`

						# Do Homer-t linking using part1 and part3 as query, part2 as reference, clear genotypes as target AF database.
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_match_all_clear
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_nomatch_all_clear

						# The first haplotypes are the ones that are positive.
						Homer_match_all_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_match_all_clear_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`
						Homer_nomatch_all_clear=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_all_clear_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`

						# Do Homer-t linking using resampled panel as the mixture.
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part1_decoded_regs.matbed.gz site1_part1_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_match_all_clear_resampled_part1
						ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${op_prefix}_site1_part3_decoded_regs.matbed.gz site1_part3_samples.list ${op_prefix}_site1_part2_decoded_regs.matbed.gz site1_part2_samples.list ${op_prefix}_site1_part1_resampled_decoded_regs.matbed.gz ${resampled_site1_part1_samples} ${op_prefix}_nomatch_all_clear_resampled_part1

						# The first haplotypes are the ones that are positive.
						Homer_match_all_clear_resampled_part1=`head -n ${n_dec_haplotypes} ${op_prefix}_match_all_clear_resampled_part1_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`
						Homer_nomatch_all_clear_resampled_part1=`head -n ${n_dec_haplotypes} ${op_prefix}_nomatch_all_clear_resampled_part1_per_query_t_stats.txt | awk '{tot+=$2}END{print tot/NR}'`

						####################################################################################################################################
						####################################################################################################################################

						# Add the leakage summary statistics.
						echo -e "${op_prefix}\t\
${decoded_err_stat}\t${decoded_shuff_err_stat}\t\
${part1_linking_stats}\t${part2_linking_stats}\t\
${LRT_match_clear}\t${LRT_nomatch_clear}\t\
${LRT_match_decoded}\t${LRT_nomatch_decoded}\t\
${LRT_match_all_clear}\t${LRT_nomatch_all_clear}\t\
${LRT_match_all_clear_resampled_part1}\t${LRT_nomatch_all_clear_resampled_part1}\t\
${Homer_match_clear}\t${Homer_nomatch_clear}\t\
${Homer_match_decoded}\t${Homer_nomatch_decoded}\t\
${Homer_match_all_clear}\t${Homer_nomatch_all_clear}\t\
${Homer_match_all_clear_resampled_part1}\t${Homer_nomatch_all_clear_resampled_part1}" >> ${LEAKAGE_OP_DIR}/SUMMARIZED_ERROR_LINKING_STATS.txt

						####################################################################################################################################
					done # var_start_posn loop.
				done # chr_id loop.
			done # cur_iter_i loop.
		done
	done

	exit 0
fi



if [[ ${cmd_option} == "-compare_Proxy_2_Original_HapFreq_stats" ]]
then

	filter_proxization_vicinity_size=6
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=3
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_PROXY2CLEAR_HAPFREQ_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	l_half_hap_win=${filter_proxization_vicinity_size}

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list
		#shuf site1_part1_samples.list > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models for part1."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} part1_${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} part1_${permute_proxying_mapping_BED}

			########################################################################
			echo "Generating the proxizing models for part2 data."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} part2_${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} part2_${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz


			########################################################################
			echo "Proxizing part1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list part1_${permute_proxying_mapping_BED} site1_part1_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox.matbed.gz site1_part1_samples.list part1_${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part1_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_cleartext_AFs_${l_half_hap_win}_${cur_iter_i}.bed
			variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_proxized_part1_AFs_${l_half_hap_win}_${cur_iter_i}.bed

			echo "Proxizing part2 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list part2_${permute_proxying_mapping_BED} site1_part2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part2_perm_prox.matbed.gz site1_part2_samples.list part2_${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part2_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_part2_proxized_AFs_${l_half_hap_win}_${cur_iter_i}.bed

			ProxyTyper -get_per_win_haplotype_frequencies_per_reference site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_matching_proxy_hapfreq_${l_half_hap_win}_${cur_iter_i}.bed &

			ProxyTyper -get_per_win_haplotype_frequencies_per_reference site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_matching_orig_hapfreq_${l_half_hap_win}_${cur_iter_i}.bed &

			ProxyTyper -get_per_win_haplotype_frequencies_per_reference site1_proxized_part2_haplocoded.matbed.gz site1_part2_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_nonmatching_proxy_hapfreq_${l_half_hap_win}_${cur_iter_i}.bed &

			ProxyTyper -get_per_win_haplotype_frequencies_per_reference site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${l_half_hap_win} ${LEAKAGE_OP_DIR}/${chr_id}_nonmatching_orig_hapfreq_${l_half_hap_win}_${cur_iter_i}.bed &

			wait
		done
	done

	exit 0
fi

if [[ ${cmd_option} == "-calculate_ProxyVar_2_OriginalVar_Correlations" ]]
then

	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=1
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_PROXYVAR2CLEARVAR_CORRELATION_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		#grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list
		shuf site1_part1_samples.list > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz

			#variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} SITE2_tag_haplocoded.matbed.gz

			# Set the site2 genotype matbed file for the current iteration and current chromosome.
			SITE2_haplocoded_tag_matbed=SITE2_tag_haplocoded.matbed.gz


			########################################################################
			echo "Resampling the genotypes for the site1."
			
			# Resample the genotypes before extracting the target regions; we need this to estimate the allele accuracy and LRT/t statistics on the decoded genotypes.
			# This is necessary because when subregs are used for resampling, we cannot compare them with attacker's variant list or for doing LRT/t attack.
			op_prefix=${cur_iter_i}_${chr_id}
			upsampled_sample_size=${n_site1_part_samples}
			N_e=0.05
			eps_allele=0
			l_seg_in_bps=10000000
			l_seg_in_cMs=0
			l_seg_in_n_vars=0
			n_threads=${N_THREADS}
			start_pos=0
			end_pos=0
			ProxyTyper -resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb site1_part1_haplocoded.matbed.gz site1_part1_samples.list \
${per_chrom_maps_dir} ${upsampled_sample_size} ${N_e} ${eps_allele} \
${l_seg_in_bps} ${l_seg_in_cMs} ${l_seg_in_n_vars} \
${n_threads} ${start_pos} ${end_pos} ${op_prefix} >& ${op_prefix}_RESAMPLING.OP

			# Summarize the resampled segments.
			ProxyTyper -summarize_sampled_segments_per_resampled_haplotype_info ${op_prefix}_resampling_pattern_signal_var_regs.sigbed.gz \
${op_prefix}_resampled_sample_ids.list site1_part1_samples.list \
${op_prefix}_resampled_segments.bed

			# After this point, all decoded data are from resampled part 1 panel.
			resampled_site1_part1_matbed=${op_prefix}_resampled_tags.matbed.gz
			resampled_site1_part1_samples=${op_prefix}_resampled_sample_ids.list

			########################################################################
			echo "Proxizing site1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${permute_proxying_mapping_BED} site1_part1_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part1_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed
			variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_proxized_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			########################################################################
			echo "Proxizing resampled site1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions ${resampled_site1_part1_matbed} ${resampled_site1_part1_samples} ${permute_proxying_mapping_BED} site1_part1_resampled_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_resampled_perm_prox.matbed.gz ${resampled_site1_part1_samples} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_resampled_proxized_part1_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_resampled_proxized_part1_haplocoded.matbed.gz ${resampled_site1_part1_samples} ${LEAKAGE_OP_DIR}/${chr_id}_site1_resampled_proxized_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			########################################################################
			echo "Proxizing site2 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${permute_proxying_mapping_BED} site1_part2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part2_perm_prox.matbed.gz site1_part2_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part2_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part2_haplocoded.matbed.gz site1_part1_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_part2_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed
			variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part2_haplocoded.matbed.gz site1_part2_samples.list ${LEAKAGE_OP_DIR}/${chr_id}_site1_proxized_part2_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			########################################################################
			# Link resampled proxized panel to original panel.
			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_resampled_proxized_part1_haplocoded.matbed.gz ${resampled_site1_part1_samples}
			mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_resampledproxy2matchingorig_iter_${cur_iter_i}.bed

			# Link resampled original panel to original panel.
			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${site1_resampled_panel_tag_haplocoded} ${resampled_site1_part1_samples}
			mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_resampledoriginal2matchingorig_iter_${cur_iter_i}.bed

			# Link proxized panel to original panel.
			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part1_haplocoded.matbed.gz site1_part1_samples.list
			mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_proxy2matchingorig_iter_${cur_iter_i}.bed

			########################################################################
			# Link to non-matching panel (part2): Link the original proxized panel to part2.
			# Link proxized panel to original part2 panel.
			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part2_haplocoded.matbed.gz site1_part2_samples.list
			mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_proxy2nonmatchingorig_iter_${cur_iter_i}.bed

			# Link resampled proxized panel to original part2 panel.
			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_resampled_proxized_part1_haplocoded.matbed.gz ${resampled_site1_part1_samples} site1_part2_haplocoded.matbed.gz site1_part2_samples.list
			mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_resampledproxy2nonmatchingorig_iter_${cur_iter_i}.bed

			## Link proxized panel to proxized part2 panel.
			#ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_proxized_part2_haplocoded.matbed.gz site1_part2_samples.list
			#mv var2var_correlations.bed ${LEAKAGE_OP_DIR}/${chr_id}_var2var_correlations_proxy2proxy_iter_${cur_iter_i}.bed
		done
	done

	exit 0
fi


if [[ ${cmd_option} == "-hapfreq_signature_attack" ]]
then
        filter_proxization_vicinity_size=6
        perm_proxization_vicinity_size=2
        var_weight=1
        coding_modulus=2
        allele_err_eps=0.00
        normalized_N_e=10

        var_weight_prob=0.3
        var2var_interaction_prob=0.8
        var2var2var_interaction_prob=0.6
        filter_weight_inversion_prob=0.5

        filter_proxization_min_n_params_per_var=2

        perm_proxy_geno_inversion_prob=0

        recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

        # This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
        # Note that for proxization requires site1 genotypes to be phased.
        USE_PHASED_SITE1_GENO=0

	chr_id=19

        N_THREADS=40

        weight_params_file=per_var_proxy.params
        permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
	ALL_sample_list=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	iters=(1)
	for cur_iter_i in ${iters[@]}
	do
		echo "Processing ${cur_iter_i}. iteration"
		
		if [[ 1 == 0 ]]
		then
			ref_sample_size=1500
			query_sample_size=500
			proxy_sample_size=500

			shuf ${ALL_sample_list} | head -n ${ref_sample_size} > ref_sample_ids.list
			grep -w -v -f ref_sample_ids.list ${ALL_sample_list} | shuf | head -n ${query_sample_size} > query_sample_ids.list
			grep -w -v -f ref_sample_ids.list ${ALL_sample_list} | grep -w -v -f query_sample_ids.list | head -n ${proxy_sample_size} > proxy_sample_ids.list
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} ref_sample_ids.list ref_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} query_sample_ids.list query_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} proxy_sample_ids.list proxy_tag_haplocoded.matbed.gz
			exit 1
		fi

		####################################################
		# Extract AFR sample id's.
		n_proxy_subject=200
		n_query_subject=200
		if [[ 1 == 1 ]]
		then
			shuf ${ALL_sample_list} | head -n ${n_proxy_subject} > proxy_sample_ids.list

			grep -w -v -f proxy_sample_ids.list ${ALL_sample_list} | shuf | head -n ${n_proxy_subject} > query_sample_ids.list
			grep -w -v -f proxy_sample_ids.list ${ALL_sample_list} | grep -w -v -f query_sample_ids.list > ref_sample_ids.list
			
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} ref_sample_ids.list ref_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} query_sample_ids.list query_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_sample_list} proxy_sample_ids.list proxy_tag_haplocoded.matbed.gz
		fi

		variation_tools -dump_plain_geno_signal_regions ref_tag_haplocoded.matbed.gz ref_sample_ids.list 1 tag_vars.bed

		sort -n -k2,2 tag_vars.bed | awk '{if(NR%2==0){print $0}}' > sub_ref_tag_vars.bed
		sort -n -k2,2 tag_vars.bed | awk '{if(NR%2==1){print $0}}' > sub_proxy_tag_vars.bed
		cp tag_vars.bed sub_ref_tag_vars.bed
		cp tag_vars.bed sub_proxy_tag_vars.bed
		#cp sub_proxy_tag_vars.bed sub_ref_tag_vars.bed

		variation_tools -extract_genotype_signals_per_region_list ref_tag_haplocoded.matbed.gz ref_sample_ids.list sub_ref_tag_vars.bed ref_tag_haplocoded_ref_sub.matbed.gz

		variation_tools -extract_genotype_signals_per_region_list query_tag_haplocoded.matbed.gz query_sample_ids.list sub_ref_tag_vars.bed query_tag_haplocoded_ref_sub.matbed.gz

		variation_tools -extract_genotype_signals_per_region_list proxy_tag_haplocoded.matbed.gz proxy_sample_ids.list sub_ref_tag_vars.bed proxy_tag_haplocoded_ref_sub.matbed.gz

		# We are assuming variants are exactly matched.
		variation_tools -extract_genotype_signals_per_region_list proxy_tag_haplocoded.matbed.gz proxy_sample_ids.list sub_ref_tag_vars.bed proxy_tag_haplocoded_proxy_sub.matbed.gz
		####################################################

		####################################################
		# Re-sample.
		echo -e "RESAMPLE\t${n_proxy_subject}\t0.125\t0.0000\t0\t0\t0" > PROXY_RESAMPLING.params
		echo -e "RECOMB\t${n_proxy_subject}\t0.125\t0.0\t0\t0\t0" >> PROXY_RESAMPLING.params
		rm -f PROXY_RESAMPLING_*
		./ProxyTyper.sh -resample_recomb_generating_panel_per_config ${chr_id} PROXY proxy_tag_haplocoded_proxy_sub.matbed.gz proxy_sample_ids.list PROXY_RESAMPLING.params ${per_chrom_maps_dir} PROXY_RESAMPLING
		if [[ ! -f PROXY_RESAMPLING_samples.list ]]
		then
			echo "Could not find site1 resampled sample list: PROXY_RESAMPLING_samples.list"
			exit 1
		fi

		if [[ ! -f PROXY_RESAMPLING_haplocoded.matbed.gz ]]
		then
			echo "Could not find site1 resampled genotypes: PROXY_RESAMPLING_haplocoded.matbed.gz"
			exit 1
		fi

		awk {'print "SITE1_"$0'} PROXY_RESAMPLING_samples.list > proxy_resampled_panel_sample_ids.list
		cp PROXY_RESAMPLING_haplocoded.matbed.gz proxy_resampled_panel_tag_haplocoded.matbed.gz
		####################################################

		####################################################
		# Set the file names.	
        resampled_proxy_haplocoded_tag_matbed=proxy_resampled_panel_tag_haplocoded.matbed.gz
        resampled_proxy_sample_list=proxy_resampled_panel_sample_ids.list

		original_proxy_haplocoded_tag_matbed=proxy_tag_haplocoded_proxy_sub.matbed.gz
        original_proxy_sample_list=proxy_sample_ids.list

        beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

        per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
        per_chrom_beagle_maps_dir=../../beagle/genetic_maps
		####################################################

		####################################################
		# Generate a new model for the current chromosome.
        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${resampled_proxy_haplocoded_tag_matbed} ${resampled_proxy_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

        # Generate the permutation proxization mapping regions.
        ProxyTyper -generate_permute_proxizing_parameters ${resampled_proxy_haplocoded_tag_matbed} ${resampled_proxy_sample_list} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}
		####################################################

		####################################################
		# Proxize the resampled proxy panel.
		ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters ${resampled_proxy_haplocoded_tag_matbed} ${resampled_proxy_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} resampled_proxy_perm_prox.matbed.gz

		ProxyTyper -permute_proxize_genotype_signal_regions resampled_proxy_perm_prox.matbed.gz ${resampled_proxy_sample_list} ${permute_proxying_mapping_BED} resampled_proxy_proxized_haplocoded_sub.matbed.gz
		####################################################

		####################################################
		# Proxize the original proxy panel.
		ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters ${original_proxy_haplocoded_tag_matbed} ${original_proxy_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} original_proxy_perm_prox.matbed.gz

		ProxyTyper -permute_proxize_genotype_signal_regions original_proxy_perm_prox.matbed.gz ${original_proxy_sample_list} ${permute_proxying_mapping_BED} original_proxy_proxized_haplocoded_sub.matbed.gz
		####################################################


		####################################################
		# Attack.
		ProxyTyper -linking_attack_per_haplotype_frequency_signatures_per_ref_panel query_tag_haplocoded_ref_sub.matbed.gz query_sample_ids.list ref_tag_haplocoded_ref_sub.matbed.gz ref_sample_ids.list resampled_proxy_proxized_haplocoded_sub.matbed.gz ${resampled_proxy_sample_list} 6 1 query_2_resampled_proxy_linking_stats.txt

		ProxyTyper -linking_attack_per_haplotype_frequency_signatures_per_ref_panel proxy_tag_haplocoded_ref_sub.matbed.gz proxy_sample_ids.list ref_tag_haplocoded_ref_sub.matbed.gz ref_sample_ids.list resampled_proxy_proxized_haplocoded_sub.matbed.gz ${resampled_proxy_sample_list} 6 1 proxy_2_resampled_proxy_linking_stats.txt

		ProxyTyper -linking_attack_per_haplotype_frequency_signatures_per_ref_panel query_tag_haplocoded_ref_sub.matbed.gz query_sample_ids.list ref_tag_haplocoded_ref_sub.matbed.gz ref_sample_ids.list original_proxy_proxized_haplocoded_sub.matbed.gz ${original_proxy_sample_list} 6 1 query_2_original_proxy_linking_stats.txt

		ProxyTyper -linking_attack_per_haplotype_frequency_signatures_per_ref_panel proxy_tag_haplocoded_ref_sub.matbed.gz proxy_sample_ids.list ref_tag_haplocoded_ref_sub.matbed.gz ref_sample_ids.list original_proxy_proxized_haplocoded_sub.matbed.gz ${original_proxy_sample_list} 6 1 proxy_2_original_proxy_linking_stats.txt
		####################################################
exit
	done

	exit 0
fi



if [[ ${cmd_option} == "-generate_PCA_proxy_original_data" ]]
then
	KG_DIR=../../DATA/1kG
	SNVs_BED_FILE=${KG_DIR}/snv_region_0.005.bed
	GENOME_SEQ_DIR=../../../../genomes/hg19
	genetic_maps_dir=../../DATA/genetic_maps

	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL/
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	sort -u -k2,2 ${KG_DIR}/super_pop_info.txt | grep AFR | cut -f1,1 | shuf > AFR_EUR_EAS_pops.list
	sort -u -k2,2 ${KG_DIR}/super_pop_info.txt | grep EUR | cut -f1,1 | shuf >> AFR_EUR_EAS_pops.list
	sort -u -k2,2 ${KG_DIR}/super_pop_info.txt | grep EAS | cut -f1,1 | shuf >> AFR_EUR_EAS_pops.list

	shuf AFR_EUR_EAS_pops.list > temp.txt
	mv temp.txt AFR_EUR_EAS_pops.list

	site1_pops=`cat AFR_EUR_EAS_pops.list`

	rm -f all_sample_id_sample_pops.list
	for cur_pop in ${site1_pops[@]}
	do
		echo "Extracting Pop.: ${cur_pop}"
		grep -w ${cur_pop} ${KG_DIR}/population_info.list | awk {'print $1"\t"$2'} >> all_sample_id_sample_pops.list
	done

	grep -w -f ${KG_DIR}/1kg_sample_ids.list all_sample_id_sample_pops.list > subj_super_pops.list

	cut -f1,1 subj_super_pops.list > site1_sample_ids.list
	cut -f2,2 subj_super_pops.list > site1_pops.list

	multicolumn_processing_tools -extract_rows_per_query_column_preserve_query_order site1_pops.list ../../DATA/1kG/super_pop_info.txt 0 1 site1_superpop_info.list
	awk 'BEGIN{FS="\t"}{print $3}' site1_superpop_info.list > site1_superpops.list

	SITE1_sample_list=site1_sample_ids.list

	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=1

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	coord_anon_cM_noise_SD=0.05

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	l_block=300

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed
	

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	rm -f *.OP
	rm -f *.bed
	rm -f *_meta_accs.txt
	rm -f *.matbed.gz

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		########################################################################################################################	
		# Use the all tag/target data.
		ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
		ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

		########################################################################
		echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_target_haplocoded.matbed.gz

		#variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_tag_haplocoded.matbed.gz
		#variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_target_haplocoded.matbed.gz

		########################################################################################################################
		# Set the Site1 and Site2 genotype names.
		SITE1_haplocoded_tag_matbed=${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${chr_id}_SITE1_target_haplocoded.matbed.gz
		########################################################################################################################	

		########################################################################################################
		# Generate the parameters using LD parameters. 
		ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

		# Generate the permutation proxization mapping regions.
		ProxyTyper -generate_permute_proxizing_parameters ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

		####################################################################################################
		# Site1 setup: Proxize tags.

		# Permute proxize the site1 data.
		ProxyTyper -permute_proxize_genotype_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${permute_proxying_mapping_BED} site1_perm_prox.matbed.gz

		# Filter proxize site1 variants.
		allele_err_eps=0
		ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_perm_prox.matbed.gz ${SITE1_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_tag_haplocoded.matbed.gz

		#####################################################################
		# Coordinate anonymization.
		anon_genetic_maps_dir=anon_coords_data
		beagle_genetic_map_file=${anon_genetic_maps_dir}/anon_beagle_map.map

		# Dump the original tag and target variants.
		variation_tools -dump_plain_geno_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} 1 tag_vars.bed
		variation_tools -dump_plain_geno_signal_regions ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} 1 target_vars.bed

		#10      .       0.000000        72765
		rm -f -r ${anon_genetic_maps_dir}
		mkdir ${anon_genetic_maps_dir}
		ProxyTyper -anonymize_tag_target_genetic_map_coords tag_vars.bed target_vars.bed ${per_chrom_maps_dir}/${chr_id}.map ${coord_anon_cM_noise_SD} tag_target_mapping ${anon_genetic_maps_dir}

		# Writet the beagle formatted genetic map file.
		awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${anon_genetic_maps_dir}/${chr_id}.map > ${beagle_genetic_map_file} 

		# Map the coordinates for both site1 and site2.
		ProxyTyper -generic_map_genotype_regs_per_src_dest_regs site1_proxized_tag_haplocoded.matbed.gz ${SITE1_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_tag_anon_coord_haplocoded.matbed.gz

		variation_tools -convert_haplocoded_2_genocoded site1_proxized_tag_anon_coord_haplocoded.matbed.gz ${SITE1_sample_list} site1_proxized_tag_anon_coord_genocoded.matbed.gz

		variation_tools -dump_plain_geno_signal_regions site1_proxized_tag_anon_coord_genocoded.matbed.gz ${SITE1_sample_list} 0 GENOCODED_PROXY.txt

		variation_tools -dump_plain_geno_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} 0 GENOCODED_ORIGINAL.txt

		#######################################################################################################################################
		# Now do re-sampling, proxize genotypes (no coordinate anonymization) and perform 
		#./ProxyTyper.sh -resample_recomb_generating_panel_per_config ${chr_id} SITE1 ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} PCA_RESAMPLING.params ${per_chrom_maps_dir} PCA_RESAMPLING
		if [[ ! -f ${SITE1_haplocoded_tag_matbed} ]]
		then
			echo "Could not find ${SITE1_haplocoded_tag_matbed}"
			exit 1
		fi

		if [[ ! -f ${SITE1_sample_list} ]]
		then
		echo "Could not find ${SITE1_sample_list}"
			exit 1
		fi

		########################################################################################################################
		# 
		op_prefix=SITE1_${chr_id}
		resampled_size=2000
		N_e_frac=0.125
		allele_eps=0
		max_l_seg_n_bps=10000000
		max_l_seg_cM=0
		max_l_seg_nvars=0
		n_threads=40
		start_posn=0
		end_posn=250000000
		./ProxyTyper.sh -resample_query_genotypes_w_length_cutoff ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${op_prefix}

		#cp ${op_prefix}_resampled_tags.matbed.gz ${chr_id}_SITE1_resampled_tag_haplocoded.matbed.gz
		#cp ${op_prefix}_resampled_sample_ids.list ${chr_id}_SITE1_resampled_samples.list

		if [[ ! -f ${op_prefix}_resampled_sample_ids.list ]]
		then
			echo "Could not find site1 resampled sample list: ${op_prefix}_resampled_sample_ids.list"
			exit 1
		fi

		if [[ ! -f ${op_prefix}_resampled_tags.matbed.gz ]]
		then
			echo "Could not find site1 resampled genotypes: ${op_prefix}_resampled_tags.matbed.gz"
			exit 1
		fi

		# Copy the proxy-panel and sample list name for site1.
		awk {'print "SITE1_"$0'} ${op_prefix}_resampled_sample_ids.list > site1_resampled_panel_sample_ids.list
		cp ${op_prefix}_resampled_tags.matbed.gz site1_resampled_panel_tag_haplocoded.matbed.gz

		# Proxize the site1 data.
		allele_err_eps=0

		# Permute proxize the site1 data.
		ProxyTyper -permute_proxize_genotype_signal_regions site1_resampled_panel_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${permute_proxying_mapping_BED} site1_resampled_panel_perm_prox.matbed.gz

		# Filter proxize site1 variants.
		allele_err_eps=0
		ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_resampled_panel_perm_prox.matbed.gz site1_resampled_panel_sample_ids.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_resampled_panel_tag_haplocoded.matbed.gz

		# Map the coordinates for both site1.
		ProxyTyper -generic_map_genotype_regs_per_src_dest_regs site1_proxized_resampled_panel_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz

		variation_tools -convert_haplocoded_2_genocoded site1_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list site1_proxized_resampled_panel_tag_anon_coord_genocoded.matbed.gz

		variation_tools -dump_plain_geno_signal_regions site1_proxized_resampled_panel_tag_anon_coord_genocoded.matbed.gz site1_resampled_panel_sample_ids.list 0 GENOCODED_RESAMPLED_PROXY.txt
		exit 0
	done
fi


if [[ ${cmd_option} == "-calculate_proxy_LD_patterns" ]]
then
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=1

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	coord_anon_cM_noise_SD=0.05

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	l_block=500
	l_step=5

	chr_ids=(19)

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed
	
	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	rm -f *.OP
	rm -f *.bed
	rm -f *_meta_accs.txt
	rm -f *.matbed.gz

	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=330
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list

	# Each iteration updates site1 and site2 subject lists.
	shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
	grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		########################################################################################################################	
		# Use the all tag/target data.
		ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
		ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

		########################################################################
		echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_target_haplocoded.matbed.gz

		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_target_haplocoded.matbed.gz

		########################################################################################################################
		# Set the Site1 and Site2 genotype names.
		SITE1_haplocoded_tag_matbed=${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${chr_id}_SITE1_target_haplocoded.matbed.gz

		SITE2_haplocoded_tag_matbed=${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${chr_id}_SITE2_target_haplocoded.matbed.gz
		########################################################################################################################	

		########################################################################################################
		# Generate the parameters using LD parameters. 
		ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

		# Generate the permutation proxization mapping regions.
		ProxyTyper -generate_permute_proxizing_parameters ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

		####################################################################################################
		# Site1 setup: Proxize tags.

		# First resample the site1 tag panel.
		op_prefix=SITE1_${chr_id}
		resampled_size=1000
		N_e_frac=0.125
		allele_eps=0
		max_l_seg_n_bps=10000000
		max_l_seg_cM=0
		max_l_seg_nvars=0
		n_threads=40
		start_posn=0
		end_posn=250000000
		./ProxyTyper.sh -resample_query_genotypes_w_length_cutoff ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${op_prefix}

		#cp ${op_prefix}_resampled_tags.matbed.gz ${chr_id}_SITE1_resampled_tag_haplocoded.matbed.gz
		#cp ${op_prefix}_resampled_sample_ids.list ${chr_id}_SITE1_resampled_samples.list

		if [[ ! -f ${op_prefix}_resampled_sample_ids.list ]]
		then
			echo "Could not find site1 resampled sample list: ${op_prefix}_resampled_sample_ids.list"
			exit 1
		fi

		if [[ ! -f ${op_prefix}_resampled_tags.matbed.gz ]]
		then
			echo "Could not find site1 resampled genotypes: ${op_prefix}_resampled_tags.matbed.gz"
			exit 1
		fi

		# Copy the proxy-panel and sample list name for site1.
		awk {'print "SITE1_"$0'} ${op_prefix}_resampled_sample_ids.list > site1_resampled_panel_sample_ids.list
		cp ${op_prefix}_resampled_tags.matbed.gz site1_resampled_panel_tag_haplocoded.matbed.gz
		#####################################################################

		# Permute proxize the site1 data.
		#ProxyTyper -permute_proxize_genotype_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${permute_proxying_mapping_BED} site1_perm_prox.matbed.gz
		ProxyTyper -permute_proxize_genotype_signal_regions site1_resampled_panel_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${permute_proxying_mapping_BED} site1_perm_prox.matbed.gz

		# Filter proxize site1 variants.
		allele_err_eps=0
		ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_perm_prox.matbed.gz site1_resampled_panel_sample_ids.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_tag_haplocoded.matbed.gz

		#####################################################################
		# Coordinate anonymization.
		anon_genetic_maps_dir=anon_coords_data
		beagle_genetic_map_file=${anon_genetic_maps_dir}/anon_beagle_map.map

		# Dump the original tag and target variants.
		variation_tools -dump_plain_geno_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} 1 tag_vars.bed
		variation_tools -dump_plain_geno_signal_regions ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} 1 target_vars.bed

		#10      .       0.000000        72765
		rm -f -r ${anon_genetic_maps_dir}
		mkdir ${anon_genetic_maps_dir}
		ProxyTyper -anonymize_tag_target_genetic_map_coords tag_vars.bed target_vars.bed ${per_chrom_maps_dir}/${chr_id}.map ${coord_anon_cM_noise_SD} tag_target_mapping ${anon_genetic_maps_dir}

		# Writet the beagle formatted genetic map file.
		awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${anon_genetic_maps_dir}/${chr_id}.map > ${beagle_genetic_map_file} 

		# Map the coordinates for both site1 and site2.
		ProxyTyper -generic_map_genotype_regs_per_src_dest_regs site1_proxized_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_tag_anon_coord_haplocoded.matbed.gz

		# Start calculating the var-var correlations.
		ProxyTyper -get_consecutive_block_variant_correlations site1_proxized_tag_anon_coord_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${l_block} ${l_step} ${chr_id}_proxy_LD_R2_corrs_${l_block}_${l_step}.txt &

		ProxyTyper -get_consecutive_block_variant_correlations ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${l_block} ${l_step} ${chr_id}_original_LD_R2_corrs_${l_block}_${l_step}.txt &

		ProxyTyper -get_consecutive_block_variant_correlations ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${l_block} ${l_step} ${chr_id}_site2_original_LD_R2_corrs_${l_block}_${l_step}.txt &

		echo "Waiting for processes.."
		wait

		exit 0
	done
fi

if [[ ${cmd_option} == "-compute_proxy2cleartext_linking_stats" ]]
then
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=20
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=(1)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_CLEAR2PROXY_LINKAGE_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2,3."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list

		# Use all variants.
		l_corr_win=20000

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz

			#variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} SITE2_tag_haplocoded.matbed.gz

			# Set the site2 genotype matbed file for the current iteration and current chromosome.
			SITE2_haplocoded_tag_matbed=SITE2_tag_haplocoded.matbed.gz

			########################################################################
			echo "Proxizing site1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${permute_proxying_mapping_BED} site1_part1_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part1_haplocoded.matbed.gz
			#variation_tools -compute_AF_per_geno_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed
			#variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_proxized_part1_haplocoded.matbed.gz_iter_${cur_iter_i}.bed

			ProxyTyper -permute_proxize_genotype_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${permute_proxying_mapping_BED} site1_part2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part2_perm_prox.matbed.gz site1_part2_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part2_haplocoded.matbed.gz
			#variation_tools -compute_AF_per_geno_signal_regions site1_part2_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part2_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			# We don't need site2 for linking tests.
			#echo "Proxizing site2 genotypes."
			#ProxyTyper -permute_proxize_genotype_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${permute_proxying_mapping_BED} site2_perm_prox.matbed.gz
			#ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site2_perm_prox.matbed.gz ${SITE2_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site2_proxized_tag_haplocoded.matbed.gz
			#variation_tools -compute_AF_per_geno_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${chr_id}_SITE2_tag_AFs.bed

			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part1_haplocoded.matbed.gz site1_part1_samples.list
			mv var2var_correlations.bed ${chr_id}_var2var_correlations_part1_iter_${cur_iter_i}.bed

			ProxyTyper -calculate_proxy2clear_var2var_correlation_stats site1_proxized_part2_haplocoded.matbed.gz site1_part2_samples.list site1_part2_haplocoded.matbed.gz site1_part2_samples.list
			mv var2var_correlations.bed ${chr_id}_var2var_correlations_part2_iter_${cur_iter_i}.bed

			ProxyTyper -calculate_proxy2clear_pairwise_distance_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${l_cleartext_windowizing_win} ${l_proxy_windowizing_win} ${l_corr_win} ${LEAKAGE_OP_DIR}/${chr_id}_proxy2site1_part1_cleartext_linking_stats_${cur_iter_i}

			ProxyTyper -calculate_proxy2clear_pairwise_distance_stats site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${l_cleartext_windowizing_win} ${l_proxy_windowizing_win} ${l_corr_win} ${LEAKAGE_OP_DIR}/${chr_id}_proxy2site1_part2_cleartext_linking_stats_${cur_iter_i}
		done # Chromosome loop.
	done # iteration loop.
	exit 0
fi

if [[ ${cmd_option} == "-calculate_Sankararaman_LRT_attack_statistics" ]]
then
	# Clean the Dyij statistics files.
	rm -f *_per_query_Sankararaman_LRT_stats.txt

	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=20
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=`seq 1 22`

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=300
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_LRT_LEAKAGES_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2,3."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} SITE2_tag_haplocoded.matbed.gz

			# Set the site2 genotype matbed file for the current iteration and current chromosome.
			SITE2_haplocoded_tag_matbed=SITE2_tag_haplocoded.matbed.gz

			########################################################################
			echo "Proxizing site1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${permute_proxying_mapping_BED} site1_part1_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part1_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed
			variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_proxized_part1_haplocoded.matbed.gz_iter_${cur_iter_i}.bed

			ProxyTyper -permute_proxize_genotype_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${permute_proxying_mapping_BED} site1_part2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part2_perm_prox.matbed.gz site1_part2_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part2_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part2_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part2_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			echo "Proxizing site2 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${permute_proxying_mapping_BED} site2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site2_perm_prox.matbed.gz ${SITE2_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site2_proxized_tag_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${chr_id}_SITE2_tag_AFs.bed

			#####################################################################
			## FOLLOWING IS THE REFERENCE FOR SETTING PANELS:
			## This is the target panel that adversary is attacking.
			#target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			##target_query_panel_matbed=site1_proxized_part1_haplocoded.matbed.gz
			#target_query_panel_sample_list_fp=site1_part1_samples.list

			## This is the panel that is used as the reference.
			#ref_panel_matbed=site1_proxized_part2_haplocoded.matbed.gz
			##ref_panel_matbed=site1_part2_haplocoded.matbed.gz
			#ref_panel_sample_list_fp=site1_part2_samples.list

			##ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			##ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			##ref_panel_sample_list_fp=${SITE2_sample_list}

			## This is the database from within which we are searching for individuals.
			#target_database_panel_matbed=site1_proxized_part3_haplocoded.matbed.gz
			##target_database_panel_matbed=site1_part3_haplocoded.matbed.gz
			#target_database_panel_sample_list_fp=site1_part3_samples.list

			#op_prefix=query_part1_target_part3

			#ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}

			####################################################################
			# Run the in the database statistics.
 			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_proxized_part1_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part1_samples.list

			op_prefix=${chr_id}_query_part1_target_proxy_part1_iter_${cur_iter_i}

			ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${MAF_cutoff} ${min_var2var_dist} ${op_prefix}

			####################################################################
			# Run the not in the database:
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_proxized_part2_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part2_samples.list

			op_prefix=${chr_id}_query_part1_target_proxy_part2_iter_${cur_iter_i}

			ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${MAF_cutoff} ${min_var2var_dist} ${op_prefix}

			####################################################################
			# Run not in the database with noproxy.
	       	# This is the target panel that adversary is attacking, it is sequenced locally.
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference; we are assuming only has access to this.
			ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part1_samples.list

			op_prefix=${chr_id}_query_part1_target_noproxy_part1_iter_${cur_iter_i}

			ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${MAF_cutoff} ${min_var2var_dist} ${op_prefix}

			####################################################################
			# Run in the database with noproxy.
			# This is the target panel that adversary is attacking.
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_part2_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part2_samples.list
	
			op_prefix=${chr_id}_query_part1_target_noproxy_part2_iter_${cur_iter_i}

			ProxyTyper -calculate_Sankararaman_LRT_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${MAF_cutoff} ${min_var2var_dist} ${op_prefix}

			####################################################################
			#awk '{tot+=$2}END{print tot/NR}' ${op_prefix}_per_query_LRT_stats.txt
		done # chromosome index.

		echo "POOLING AND SUMMARIZING RESULTS FROM ALL CHROMOSOMES."

		# Pool and summarize the results for this iteration.
		op_prefix=query_part1_target_proxy_part1_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_Sankararaman_LRT_stats.txt > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Sankararaman_LRT_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Sankararaman_LRT_stats.txt

		op_prefix=query_part1_target_proxy_part2_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_Sankararaman_LRT_stats.txt > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Sankararaman_LRT_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Sankararaman_LRT_stats.txt

		op_prefix=query_part1_target_noproxy_part1_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_Sankararaman_LRT_stats.txt > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Sankararaman_LRT_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Sankararaman_LRT_stats.txt

		op_prefix=query_part1_target_noproxy_part2_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_Sankararaman_LRT_stats.txt > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Sankararaman_LRT_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Sankararaman_LRT_stats.txt
	done # iteration index.
fi


if [[ ${cmd_option} == "-calculate_Homer_t_attack_statistics" ]]
then
	# Clean the Dyij statistics files.
	rm -f *_per_query_per_var_t_stat_Dyij.txt.gz
	
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	n_iters=20
	iters_i=`seq 1 ${n_iters}`

	#chr_ids=(19 20 21 22)
	chr_ids=`seq 1 22`

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=300
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_t_LEAKAGES_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	for cur_iter_i in ${iters_i[@]}
	do
		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2,3."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/2)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

	        echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Generating the proxizing models."
	        # Generate a new model for the current chromosome.
	        ProxyTyper -generate_save_per_site_mixing_parameters_LD_aware ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

        	# Generate the permutation proxization mapping regions.
	        ProxyTyper -generate_permute_proxizing_parameters ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} SITE2_tag_haplocoded.matbed.gz

			# Set the site2 genotype matbed file for the current iteration and current chromosome.
			SITE2_haplocoded_tag_matbed=SITE2_tag_haplocoded.matbed.gz

			########################################################################
			echo "Proxizing site1 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${permute_proxying_mapping_BED} site1_part1_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part1_perm_prox.matbed.gz site1_part1_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part1_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part1_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed
			variation_tools -compute_AF_per_geno_signal_regions site1_proxized_part1_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_proxized_part1_haplocoded.matbed.gz_iter_${cur_iter_i}.bed

			ProxyTyper -permute_proxize_genotype_signal_regions site1_part2_haplocoded.matbed.gz site1_part2_samples.list ${permute_proxying_mapping_BED} site1_part2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_part2_perm_prox.matbed.gz site1_part2_samples.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_part2_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions site1_part2_haplocoded.matbed.gz site1_part1_samples.list ${chr_id}_site1_part2_haplocoded.matbed.gz_AFs_iter_${cur_iter_i}.bed

			echo "Proxizing site2 genotypes."
			ProxyTyper -permute_proxize_genotype_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${permute_proxying_mapping_BED} site2_perm_prox.matbed.gz
			ProxyTyper -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site2_perm_prox.matbed.gz ${SITE2_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site2_proxized_tag_haplocoded.matbed.gz
			variation_tools -compute_AF_per_geno_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${chr_id}_SITE2_tag_AFs.bed

			#####################################################################
			## FOLLOWING IS THE REFERENCE FOR SETTING PANELS:
			## This is the target panel that adversary is attacking.
			#target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			##target_query_panel_matbed=site1_proxized_part1_haplocoded.matbed.gz
			#target_query_panel_sample_list_fp=site1_part1_samples.list

			## This is the panel that is used as the reference.
			#ref_panel_matbed=site1_proxized_part2_haplocoded.matbed.gz
			##ref_panel_matbed=site1_part2_haplocoded.matbed.gz
			#ref_panel_sample_list_fp=site1_part2_samples.list

			##ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			##ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			##ref_panel_sample_list_fp=${SITE2_sample_list}

			## This is the database from within which we are searching for individuals.
			#target_database_panel_matbed=site1_proxized_part3_haplocoded.matbed.gz
			##target_database_panel_matbed=site1_part3_haplocoded.matbed.gz
			#target_database_panel_sample_list_fp=site1_part3_samples.list

			#op_prefix=query_part1_target_part3

			#ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}

			####################################################################
			# Run the in the database statistics.
 			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_proxized_part1_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part1_samples.list

			op_prefix=${chr_id}_query_part1_target_proxy_part1_iter_${cur_iter_i}

			ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}

			####################################################################
			# Run the not in the database:
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=site2_proxized_tag_haplocoded.matbed.gz
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_proxized_part2_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part2_samples.list

			op_prefix=${chr_id}_query_part1_target_proxy_part2_iter_${cur_iter_i}

			ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}
			####################################################################
			# Run not in the database with noproxy.
	       	# This is the target panel that adversary is attacking, it is sequenced locally.
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference; we are assuming only has access to this.
			ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part1_samples.list

			op_prefix=${chr_id}_query_part1_target_noproxy_part1_iter_${cur_iter_i}

			ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}

			####################################################################
			# Run in the database with noproxy.
			# This is the target panel that adversary is attacking.
			target_query_panel_matbed=site1_part1_haplocoded.matbed.gz
			target_query_panel_sample_list_fp=site1_part1_samples.list

			# This is the panel that is used as the reference.
			ref_panel_matbed=${SITE2_haplocoded_tag_matbed}
			ref_panel_sample_list_fp=${SITE2_sample_list}

			# This is the database from within which we are searching for individuals.
			target_database_panel_matbed=site1_part2_haplocoded.matbed.gz
			target_database_panel_sample_list_fp=site1_part2_samples.list
	
			op_prefix=${chr_id}_query_part1_target_noproxy_part2_iter_${cur_iter_i}

			ProxyTyper -calculate_Homer_t_statistics_on_proxized_panels ${target_query_panel_matbed} ${target_query_panel_sample_list_fp} ${ref_panel_matbed} ${ref_panel_sample_list_fp} ${target_database_panel_matbed} ${target_database_panel_sample_list_fp} ${op_prefix}

			####################################################################
			#awk '{tot+=$2}END{print tot/NR}' ${op_prefix}_per_query_LRT_stats.txt
		done # chromosome index.

		echo "POOLING AND SUMMARIZING RESULTS FROM ALL CHROMOSOMES."

		# Pool and summarize the results for this iteration.
		op_prefix=query_part1_target_proxy_part1_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_per_var_t_stat_Dyij.txt.gz > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Homer_t_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Homer_t_stats.txt

		op_prefix=query_part1_target_proxy_part2_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_per_var_t_stat_Dyij.txt.gz > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Homer_t_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Homer_t_stats.txt

		op_prefix=query_part1_target_noproxy_part1_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_per_var_t_stat_Dyij.txt.gz > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Homer_t_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Homer_t_stats.txt

		op_prefix=query_part1_target_noproxy_part2_iter_${cur_iter_i}
		ls *_${op_prefix}_per_query_per_var_t_stat_Dyij.txt.gz > ${op_prefix}_files.list
		ProxyTyper -pool_summarize_Homer_t_statistics_per_query ${op_prefix}_files.list site1_part1_samples.list ${LEAKAGE_OP_DIR}/${op_prefix}_summarized_Homer_t_stats.txt
	done # iteration index.
fi

if [[ ${cmd_option} == "-decomposed_impute_untyped_vars" ]]
then
	filter_proxization_vicinity_size=5
	perm_proxization_vicinity_size=0
	var_weight=1
	coding_modulus=2
	allele_err_eps=0.00
	normalized_N_e=25

	var_weight_prob=0.3
	var2var_interaction_prob=0.8
	var2var2var_interaction_prob=0.6
	filter_weight_inversion_prob=0.5

	filter_proxization_min_n_params_per_var=2

	perm_proxy_geno_inversion_prob=0

	recomb_rate_dir=../../DATA/beagle_genetic_maps/genetic_maps

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	N_THREADS=40

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	#chr_ids=(19 20 21 22)
	chr_ids=(22)
	
	l_cleartext_windowizing_win=0
	l_proxy_windowizing_win=0

	#PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_ALL
	#ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list
	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=660
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list
	MAF_cutoff=0.05
	min_var2var_dist=10000

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	LEAKAGE_OP_DIR="SUMMARIZED_UNTYPED_LD_STATS_${date_str}"

	if [[ -d ${LEAKAGE_OP_DIR} ]]
	then
		echo "${LEAKAGE_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${LEAKAGE_OP_DIR}
	mkdir ${LEAKAGE_OP_DIR}

	l_half_hap_win=6

		echo "********************************************************************"
		echo "PROCESSING ITERATION ${cur_iter_i}; subsampling site1 parts 1,2."
		echo "This makes sure we use the same sample set for all chromosomes."

		# Each iteration updates site1 and site2 subject lists.
		shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
		#grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

		# Extract half of the samples in SITE1's list.
		n_site1_part_samples=`wc -l ${SITE1_sample_list} | awk {'print int($1/6)'}`

		echo "Extracting ${n_site1_part_samples} subjects for each part of site1."

		shuf ${SITE1_sample_list} | head -n ${n_site1_part_samples} > site1_part1_samples.list
		grep -v -w -f site1_part1_samples.list ${SITE1_sample_list} > site1_part2_samples.list
		#shuf site1_part1_samples.list > site1_part2_samples.list

		for chr_id in ${chr_ids[@]}
		do
			########################################################################################################################	
			# Use the all tag/target data.
			ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
			ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

			echo "Processing chromosome ${chr_id}"

			########################################################################
			echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} site1_part1_samples.list site1_part1_target_haplocoded.matbed.gz

			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_tag_haplocoded.matbed.gz
			variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} site1_part2_samples.list site1_part2_target_haplocoded.matbed.gz
		
			SITE1_haplocoded_tag_matbed=site1_part1_tag_haplocoded.matbed.gz
			SITE1_haplocoded_target_matbed=site1_part1_target_haplocoded.matbed.gz
			SITE1_sample_list=site1_part1_samples.list

			SITE2_haplocoded_tag_matbed=site1_part2_tag_haplocoded.matbed.gz
                        SITE2_haplocoded_target_matbed=site1_part2_target_haplocoded.matbed.gz
                        SITE2_sample_list=site1_part2_samples.list

			min_AAF_per_decomp_var=0.999
			ProxyTyper -simple_decompose_untyped_variants site1_part2_tag_haplocoded.matbed.gz \
site1_part2_target_haplocoded.matbed.gz \
site1_part2_samples.list \
${min_AAF_per_decomp_var} PART2

			variation_tools -compute_AF_per_geno_signal_regions site1_part2_target_haplocoded.matbed.gz site1_part2_samples.list site1_part2_AFs.bed
			variation_tools -compute_AF_per_geno_signal_regions PART2_decomposed.matbed.gz site1_part2_samples.list site1_part2_decomp_AFs.bed

			USE_PHASED_SITE1_GENO=0
			beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

			rm -f -r decomposed_impute_beagle_dir
			mkdir decomposed_impute_beagle_dir			
			./ProxyTyper.sh -setup_BEAGLE_files site1_part1_tag_haplocoded.matbed.gz site1_part1_samples.list site1_part2_tag_haplocoded.matbed.gz PART2_decomposed.matbed.gz site1_part2_samples.list ${USE_PHASED_SITE1_GENO} decomposed_impute_beagle_dir
			./ProxyTyper.sh -run_BEAGLE ${beagle_genetic_map_file} decomposed_impute_beagle_dir

			rm -f -r clear_impute_beagle_dir
			mkdir clear_impute_beagle_dir
			./ProxyTyper.sh -setup_BEAGLE_files site1_part1_tag_haplocoded.matbed.gz site1_part1_samples.list site1_part2_tag_haplocoded.matbed.gz site1_part2_target_haplocoded.matbed.gz site1_part2_samples.list ${USE_PHASED_SITE1_GENO} clear_impute_beagle_dir
		        ./ProxyTyper.sh -run_BEAGLE ${beagle_genetic_map_file} clear_impute_beagle_dir


			variation_tools -dump_plain_geno_signal_regions site1_part1_target_haplocoded.matbed.gz site1_part1_samples.list 1 targets.bed

                        annot_region_tools -BED_2_Interval targets.bed targets.bed.interval
                        awk {'print $0"\t0\t1.0"'} targets.bed.interval > targets.bed.interval_ext.int

                        ProxyTyper -combine_BEAGLE_imputed_decomposed_genotype_probabilities clear_impute_beagle_dir/imputed.op.vcf.gz targets.bed.interval_ext.int site1_part1_samples.list cleartext

			ProxyTyper -combine_BEAGLE_imputed_decomposed_genotype_probabilities decomposed_impute_beagle_dir/imputed.op.vcf.gz PART2_untyped_decomposing.interval site1_part1_samples.list decomp


			# Convert all to genocoded, including known genotypes.
			variation_tools -convert_haplocoded_2_genocoded cleartext_untyped_undecomp.matbed.gz site1_part1_samples.list cleartext_untyped_undecomp_genocoded.matbed.gz

			variation_tools -convert_haplocoded_2_genocoded decomp_untyped_undecomp.matbed.gz site1_part1_samples.list decomp_untyped_undecomp_genocoded.matbed.gz

			variation_tools -convert_haplocoded_2_genocoded site1_part1_target_haplocoded.matbed.gz site1_part1_samples.list site1_part1_target_genocoded.matbed.gz

			# Get R2.
			ProxyTyper -get_R2_per_imputed_genotypes cleartext_untyped_undecomp_genocoded.matbed.gz site1_part1_samples.list site1_part1_target_genocoded.matbed.gz site1_part1_samples.list
			mv R2_stats.txt cleartext_R2_stats.txt
			awk -f ../GET_ACC_SUMMARY.awk cleartext_R2_stats.txt

			ProxyTyper -get_R2_per_imputed_genotypes decomp_untyped_undecomp_genocoded.matbed.gz site1_part1_samples.list site1_part1_target_genocoded.matbed.gz site1_part1_samples.list
			mv R2_stats.txt decomp_R2_stats.txt
			awk -f ../GET_ACC_SUMMARY.awk decomp_R2_stats.txt

			exit
		done

	exit 0
fi


if [[ ${cmd_option} == "-do_impute_site1_per_proxy_protocol" ]]
then
	#BEAGLE_JAR=/internal/aharmanci1/dir/PAPER_DATA/ProxyGIMP_5.2023/beagle/beagle.13Mar20.38e.jar
	##CHR22_MAP=/internal/aharmanci1/dir/PAPER_DATA/MetaGIMS_5.2023/beagle/genetic_maps/plink.chr22.GRCh37.map
	#LOHAMMER_EXEC=/internal/aharmanci1/dir/codebase/genomics-codebase/Genomics/LoHaMMer_Deployment/bin/LoHaMMer
	#N_THREADS=80

	## Extract site1 and site2 genotypes: Use AFR ALL data, single pop.
    PER_CHROM_SITE1_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1/
	PER_CHROM_SITE2_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1/
	PER_CHROM_PRESAMPLED_SITE1_GENO_DATA_DIR=PRESAMPLED_SITE1
	PER_CHROM_PRESAMPLED_SITE2_GENO_DATA_DIR=PRESAMPLED_SITE2

	chr_ids=(19)

	RUN_CENTRAL_IMPUTATION=0

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		ASSEMBLY_ID=hg19
		filter_proxization_vicinity_size=6
		perm_proxization_vicinity_size=0
		var_weight=1
		coding_modulus=2
		allele_err_eps=0.00

		# Variant selector probabilities.
		var_weight_prob=0.3
		var2var_interaction_prob=0.8
		var2var2var_interaction_prob=0.4
		filter_proxization_min_n_params_per_var=2

		weight_inversion_prob=0.5
		permute_proxy_geno_inv_prob=0.0

		coord_anon_cM_noise_SD=0.05

		untyped_var_perm_n_vicinity=20
		untyped_var_allele_switch_prob=0.5

		variant_decomp_min_AAF=0.1

		weight_params_file=per_var_proxy.params

		normalized_N_e=15

		SITE1_haplocoded_tag_matbed=${PER_CHROM_SITE1_GENO_DATA_DIR}/${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${PER_CHROM_SITE1_GENO_DATA_DIR}/${chr_id}_SITE1_target_haplocoded.matbed.gz
		SITE1_sample_list=${PER_CHROM_SITE1_GENO_DATA_DIR}/site1_sample_ids.list

		SITE1_presampled_tag_matbed=${PER_CHROM_PRESAMPLED_SITE1_GENO_DATA_DIR}/${chr_id}_SITE1_resampled_tag_haplocoded.matbed.gz
		SITE1_presampled_sample_list=${PER_CHROM_PRESAMPLED_SITE1_GENO_DATA_DIR}/${chr_id}_SITE1_resampled_samples.list

		SITE2_haplocoded_tag_matbed=${PER_CHROM_SITE2_GENO_DATA_DIR}/${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${PER_CHROM_SITE2_GENO_DATA_DIR}/${chr_id}_SITE2_target_haplocoded.matbed.gz
		SITE2_sample_list=${PER_CHROM_SITE2_GENO_DATA_DIR}/site2_sample_ids.list

		SITE2_presampled_tag_matbed=${PER_CHROM_PRESAMPLED_SITE2_GENO_DATA_DIR}/${chr_id}_SITE2_resampled_tag_haplocoded.matbed.gz
		SITE2_presampled_target_matbed=${PER_CHROM_PRESAMPLED_SITE2_GENO_DATA_DIR}/${chr_id}_SITE2_resampled_target_haplocoded.matbed.gz
		SITE2_presampled_sample_list=${PER_CHROM_PRESAMPLED_SITE2_GENO_DATA_DIR}/${chr_id}_SITE2_resampled_samples.list

		# Do central imputation on the current subjects.
		if [[ ${RUN_CENTRAL_IMPUTATION} == 1 ]]
		then
            ASSEMBLY_ID=hg19

            beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

            echo "Performing site2 imputation"
            # Calculate the pure site2-based imputation accuracy.
            ./ProxyTyper.sh -centralized_beagle_impute ${chr_id} ${ASSEMBLY_ID} ${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} ${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} ${beagle_genetic_map_file}
            cp -r central_impute_beagle_dir ${chr_id}_central_impute_beagle_dir
		fi

		per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
		#per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

		./ProxyTyper.sh -impute_site1_per_proxy_protocol ${chr_id} ${ASSEMBLY_ID} \
${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} \
${SITE1_presampled_tag_matbed} ${SITE1_presampled_sample_list} \
${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} \
${SITE2_presampled_tag_matbed} ${SITE2_presampled_target_matbed} ${SITE2_presampled_sample_list} \
${filter_proxization_vicinity_size} ${perm_proxization_vicinity_size} ${coding_modulus} ${allele_err_eps} \
${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_proxization_min_n_params_per_var} \
${normalized_N_e} ${weight_inversion_prob} \
${permute_proxy_geno_inv_prob} \
${coord_anon_cM_noise_SD} \
${untyped_var_perm_n_vicinity} ${untyped_var_allele_switch_prob} \
${variant_decomp_min_AAF} \
${per_chrom_maps_dir}
	done
	exit 0
fi


if [[ ${cmd_option} == "-do_central_imputes" ]]
then
    #PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/1kG_2_per_pop_Site1/
    PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/KG_Synthetic_Array/1kG_GENOTYPE/1kG_2_per_pop_Site1/

	rm -f -r central_impute_beagle_dir meta_gen_imputed_beagle_dir

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		ASSEMBLY_ID=hg19
		
		SITE1_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE1_target_haplocoded.matbed.gz
		SITE1_sample_list=${PER_CHROM_GENO_DATA_DIR}/site1_sample_ids.list

		SITE2_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_SITE2_target_haplocoded.matbed.gz
		SITE2_sample_list=${PER_CHROM_GENO_DATA_DIR}/site2_sample_ids.list

		beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

		echo "Performing site2 imputation"
		# Calculate the pure site2-based imputation accuracy.
		./ProxyTyper.sh -centralized_beagle_impute ${chr_id} ${ASSEMBLY_ID} ${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} ${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} ${beagle_genetic_map_file}
		cp -r central_impute_beagle_dir ${chr_id}_central_impute_beagle_dir
	done
	exit 0
fi

if [[ ${cmd_option} == "-grid_impute_per_proxized_tags_original_panel" ]]
then
	rm -f *.OP
	rm -f *.bed
	rm -f *_meta_accs.txt
	rm -f *.matbed.gz

	PER_CHROM_GENO_DATA_DIR=../../DATA/Array_Genotype_Data/Illumina_Duo1M/1kG_GENOTYPE/AFR_ALL
	ALL_SAMPLE_LIST=${PER_CHROM_GENO_DATA_DIR}/ALL_sample_ids.list

	per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps
	per_chrom_beagle_maps_dir=../../beagle/genetic_maps/

	SITE1_SAMPLE_SIZE=330
	SITE1_sample_list=site1_sample_ids.list
	SITE2_sample_list=site2_sample_ids.list

	# Setup the output directory for summarized statistics.
	date_str=`date +%d_%m_%y_%H_%M_%S`
	STATS_OP_DIR="SUMMARIZED_GRID_ACCURACY_STATS_${date_str}"

	if [[ -d ${STATS_OP_DIR} ]]
	then
		echo "${STATS_OP_DIR} exists."
		exit 1
	fi

	rm -f -r ${STATS_OP_DIR}
	mkdir ${STATS_OP_DIR}

	# Each iteration updates site1 and site2 subject lists.
	shuf ${ALL_SAMPLE_LIST} | head -n ${SITE1_SAMPLE_SIZE} > ${SITE1_sample_list}
	grep -w -v -f ${SITE1_sample_list} ${ALL_SAMPLE_LIST} > ${SITE2_sample_list}

	for chr_id in ${chr_ids[@]}
	do
		echo "--------------------------------------"
		echo "Processing chromosome ${chr_id}"
		echo "--------------------------------------"

		########################################################################################################################	
		# Use the all tag/target data.
		ALL_haplocoded_tag_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_tag_haplocoded.matbed.gz
		ALL_haplocoded_target_matbed=${PER_CHROM_GENO_DATA_DIR}/${chr_id}_ALL_target_haplocoded.matbed.gz

		########################################################################
		echo "Extracting the genotypes for the site1 parts and site 2 genotypes."
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE1_sample_list} ${chr_id}_SITE1_target_haplocoded.matbed.gz

		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_tag_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_tag_haplocoded.matbed.gz
		variation_tools -extract_genotype_signals_per_subsample_list ${ALL_haplocoded_target_matbed} ${ALL_SAMPLE_LIST} ${SITE2_sample_list} ${chr_id}_SITE2_target_haplocoded.matbed.gz

		########################################################################################################################
		# Set the Site1 and Site2 genotype names.
		SITE1_haplocoded_tag_matbed=${chr_id}_SITE1_tag_haplocoded.matbed.gz
		SITE1_haplocoded_target_matbed=${chr_id}_SITE1_target_haplocoded.matbed.gz

		SITE2_haplocoded_tag_matbed=${chr_id}_SITE2_tag_haplocoded.matbed.gz
		SITE2_haplocoded_target_matbed=${chr_id}_SITE2_target_haplocoded.matbed.gz
		########################################################################################################################	

		ASSEMBLY_ID=hg19
		filter_proxization_vicinity_size=6
		perm_proxization_vicinity_size=1
		var_weight=1
		coding_modulus=2
		allele_err_eps=0.00

		# Variant selector probabilities.
		var_weight_prob=0.3
		var2var_interaction_prob=0.8
		var2var2var_interaction_prob=0.4
		filter_proxization_min_n_params_per_var=1

		weight_inversion_prob=0.5
		permute_proxy_geno_inv_prob=0.0

		weight_params_file=per_var_proxy.params

		normalized_N_e=25
		
		coord_anon_cM_noise_SD=0.05

		untyped_var_perm_n_vicinity=20
		untyped_var_allele_switch_prob=0.5

		# Enumerate parameters, then loop over all of them,
		echo "null" | awk -f enumerate_parameters.awk

		param_combs=`cat param_enums.txt`

		Ne_par_i=`awk -v par_2_select=normalized_N_e 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`

		var2var_par_i=`awk -v par_2_select=var2var_interaction_prob 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		var_weight_par_i=`awk -v par_2_select=var_weight_prob 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		filt_proxy_vicinity_par_i=`awk -v par_2_select=filter_proxization_vicinity_size 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		var2var2var_par_i=`awk -v par_2_select=var2var2var_interaction_prob 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		filt_proxy_min_n_params_par_i=`awk -v par_2_select=filter_proxization_min_n_params_per_var 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		perm_proxy_vicinity_par_i=`awk -v par_2_select=perm_proxization_vicinity_size 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		coord_anon_cM_noise_SD_par_i=`awk -v par_2_select=coord_anon_cM_noise_SD 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`
		untyped_var_perm_n_vicinity_par_i=`awk -v par_2_select=untyped_var_perm_n_vicinity 'BEGIN{par_select_i=0}{if($1==par_2_select){par_select_i=NR}}END{print par_select_i+1}' param_keys.list`

		echo "Column indices for the parameters:"
		echo "Ne_par_i: ${Ne_par_i}
var2var_par_i: ${var2var_par_i}
var_weight_par_i: ${var_weight_par_i}
filt_proxy_vicinity_par_i: ${filt_proxy_vicinity_par_i}
var2var2var_par_i: ${var2var2var_par_i}
filt_proxy_min_n_params_par_i: ${filt_proxy_min_n_params_par_i}
perm_proxy_vicinity_par_i: ${perm_proxy_vicinity_par_i}
coord_anon_cM_noise_SD_par_i: ${coord_anon_cM_noise_SD_par_i}
untyped_var_perm_n_vicinity_par_i: ${untyped_var_perm_n_vicinity_par_i}"

		for cur_params_line in ${param_combs[@]}
		do
			echo "----------------------------------------------"

			enum_par_id=`echo $cur_params_line | awk -v par_select_i=1 {'split($1,arr, "*");print arr[par_select_i]'}`

			normalized_N_e=`echo $cur_params_line | awk -v par_select_i=${Ne_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			var2var_interaction_prob=`echo $cur_params_line | awk -v par_select_i=${var2var_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			var_weight_prob=`echo $cur_params_line | awk -v par_select_i=${var_weight_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			filter_proxization_vicinity_size=`echo $cur_params_line | awk -v par_select_i=${filt_proxy_vicinity_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			var2var2var_interaction_prob=`echo $cur_params_line | awk -v par_select_i=${var2var2var_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			filter_proxization_min_n_params_per_var=`echo $cur_params_line | awk -v par_select_i=${filt_proxy_min_n_params_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			perm_proxization_vicinity_size=`echo $cur_params_line | awk -v par_select_i=${perm_proxy_vicinity_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			coord_anon_cM_noise_SD=`echo $cur_params_line | awk -v par_select_i=${coord_anon_cM_noise_SD_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`
			untyped_var_perm_n_vicinity=`echo $cur_params_line | awk -v par_select_i=${untyped_var_perm_n_vicinity_par_i} {'split($1,arr, "*");print arr[par_select_i]'}`

			per_chrom_maps_dir=../../DATA/beagle_genetic_maps/genetic_maps

			# Following are used by BEAGLE.
			per_chrom_beagle_maps_dir=../../beagle/genetic_maps/
			#beagle_genetic_map_file=${per_chrom_beagle_maps_dir}/plink.chr${chr_id}.GRCh37.map

			echo "${enum_par_id} Parameters: 
ASSEMBLY_ID=${ASSEMBLY_ID}
filter_proxization_vicinity_size=${filter_proxization_vicinity_size}
perm_proxization_vicinity_size=${perm_proxization_vicinity_size}
var_weight=${var_weight}
coding_modulus=${coding_modulus}
allele_err_eps=${allele_err_eps}

var_weight_prob=${var_weight_prob}
var2var_interaction_prob=${var2var_interaction_prob}
var2var2var_interaction_prob=${var2var2var_interaction_prob}
filter_proxization_min_n_params_per_var=${filter_proxization_min_n_params_per_var}

weight_inversion_prob=${weight_inversion_prob}
permute_proxy_geno_inv_prob=${permute_proxy_geno_inv_prob}

normalized_N_e=${normalized_N_e}
coord_anon_cM_noise_SD=${coord_anon_cM_noise_SD}
untyped_var_perm_n_vicinity=${untyped_var_perm_n_vicinity}"

			prefix=${chr_id}_${enum_par_id}_${filter_proxization_vicinity_size}_${perm_proxization_vicinity_size}_${var_weight}_${coding_modulus}_${allele_err_eps}_${var_weight_prob}_${var2var_interaction_prob}_${var2var2var_interaction_prob}_${filter_proxization_min_n_params_per_var}_${weight_inversion_prob}_${permute_proxy_geno_inv_prob}_${normalized_N_e}_${coord_anon_cM_noise_SD}_${untyped_var_perm_n_vicinity}
	
			./ProxyTyper.sh -impute_site1_per_proxized_variants ${chr_id} ${ASSEMBLY_ID} \
${SITE1_haplocoded_tag_matbed} ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} \
${SITE2_haplocoded_tag_matbed} ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} \
${filter_proxization_vicinity_size} ${perm_proxization_vicinity_size} ${coding_modulus} ${allele_err_eps} \
${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_proxization_min_n_params_per_var} \
${normalized_N_e} \
${weight_inversion_prob} ${permute_proxy_geno_inv_prob} \
${coord_anon_cM_noise_SD} \
${untyped_var_perm_n_vicinity} ${untyped_var_allele_switch_prob} \
${per_chrom_maps_dir} >& ${STATS_OP_DIR}/ALL_OP_${prefix}.OP

			cp ${chr_id}_meta_accs.txt ${STATS_OP_DIR}/${prefix}_meta_accs.txt
			cp per_var_vicinity_stats.bed ${STATS_OP_DIR}/${prefix}_per_var_vicinity_stats.bed
			cp proxizing_permute_reg_pairs.bed ${STATS_OP_DIR}/${prefix}_proxizing_permute_reg_pairs.bed
		done
	done
	exit 0
fi

