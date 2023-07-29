#!/bin/bash

########################################################################################################################
# Move to data_config.params ?
BEAGLE_JAR=/internal/aharmanci1/dir/PAPER_DATA/ProxyTyper_5.2023/beagle/beagle.13Mar20.38e.jar
PROXYTYPER_EXEC=ProxyTyper_Release
N_THREADS=40
#ASSEMBLY_ID=hg19
#GENOME_SEQ_DIR=../../../../genomes/hg19/
GENOME_SEQ_DIR=none
########################################################################################################################

if [[ ! -f ${BEAGLE_JAR} ]]
then
	echo "Could not find BEAGLE jar file @ ${BEAGLE_JAR}"
	exit 1
fi

proxytyper_exec_check=`type -P ${PROXYTYPER_EXEC}`
if [[ "${proxytyper_exec_check}" == "" ]]
then
	echo "Could not find ProxyTyper executable @ ${PROXYTYPER_EXEC}"
	exit 1
fi

exec_check=`type -P fetchChromSizes`
if [[ "${exec_check}" == "" ]]
then
	echo "Could not find fetchChromSizes executable"
	exit 1
fi

if [[ $# -lt 2 ]]
then
	echo "USAGE: $0 [Options] [Arguments]
Options:
	-clean_directory [Directory]
	-import_VCF
	-extract_tag_target_genotypes
BEAGLE Imputation:
	-setup_BEAGLE_files
	-run_BEAGLE
Imputation:
	Central Imputation:
		-centralized_beagle_impute
	Proxy Imputation:
		-impute_site1_per_proxy_protocol
Resampling:
	-resample_query_genotypes_w_length_cutoff: This option is used for re-sampling tag variants (Site1).
	-resample_tag_target_genotypes_w_length_cutoff: This option is used for re-sampling tags and targets together (Site2).
Accuracy Estimation:
	-get_R2_per_imputed_genotypes"
	
	exit 1
fi

cmd_option=$1

if [[ ${cmd_option} == "-clean_directory" ]]
then
	echo "Cleaning current working directory..."

	# Clean all genotype data.
	rm -f *.matbed.gz *.matbed.gz.txt *.list *.OP *.matbed *.sigbed.gz *.params *.bed

	# Delete txt files conditionally; dont remove accuracies.
	ls *.txt | grep -v accs.txt | grep -v LEAKAGES_PER_R.txt | xargs -Ifiles rm -f files

	exit 0
fi

if [[ ${cmd_option} == "-resample_tag_target_genotypes_w_length_cutoff" ]]
then
	if [[ $# -lt 13 ]]
	then
		echo "$0 $1 [Haplocoded query tag genotypes matbed] [Haplocoded query target genotypes matbed] [Query sample ids] 
[Upsampled size] [Genetic map directory] 
[N_e (Eff. pop. size)] [Allelic error rate] 
[Max segment length in bps] [Max segment length in cMs] [Max segment length in # variants] 
[# threads]
[Output prefix]"
		exit 1
	fi	
	
	haplo_tag_geno_matbed=$2
	haplo_target_geno_matbed=$3
	sample_ids=$4
	upsampled_sample_size=$5
	genetic_maps_dir=$6
	N_e=$7
	eps_allele=$8
	l_seg_in_bps=$9
	l_seg_in_cMs=${10}
	l_seg_in_n_vars=${11}
	n_threads=${12}
	op_prefix=${13}
	
	# This is fixed for the chromosome.
	start_pos=0
	end_pos=250000000
	
	if [[ ! -f ${haplo_tag_geno_matbed} ]]
	then
		echo "Could not find ${haplo_tag_geno_matbed}"
		exit 1
	fi

	if [[ ! -f ${haplo_target_geno_matbed} ]]
	then
		echo "Could not find ${haplo_target_geno_matbed}"
		exit 1
	fi
	
	if [[ ! -f ${sample_ids} ]]
	then
		echo "Could not find ${sample_ids}"
		exit 1
	fi
	
	if [[ ! -d ${genetic_maps_dir} ]]
	then
		echo "Could not find ${genetic_maps_dir}"
		exit 1
	fi
	
	${PROXYTYPER_EXEC} -resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling \
${haplo_tag_geno_matbed} ${haplo_target_geno_matbed} ${sample_ids} \
${genetic_maps_dir} ${upsampled_sample_size} ${N_e} ${eps_allele} \
${l_seg_in_bps} ${l_seg_in_cMs} ${l_seg_in_n_vars} \
${n_threads} ${start_pos} ${end_pos} ${op_prefix} >& ${op_prefix}_RESAMPLING.OP
	ret_res=$?

	# Summarize resampling patterns.
	${PROXYTYPER_EXEC} -summarize_sampled_segments_per_resampled_haplotype_info ${op_prefix}_resampling_pattern_signal_var_regs.sigbed.gz \
${op_prefix}_resampled_sample_ids.list ${sample_ids} \
${op_prefix}_resampled_segments.bed
	
	exit ${ret_res}
fi

if [[ ${cmd_option} == "-resample_query_genotypes_w_length_cutoff" ]]
then
	if [[ $# -lt 12 ]]
	then
		echo "$0 $1 [Haplocoded query tag genotypes matbed] [Query sample ids] 
[Upsampled size] [Genetic map directory] 
[N_e (Eff. pop. size)] [Allelic error rate] 
[Max segment length in bps] [Max segment length in cMs] [Max segment length in # variants] 
[# threads]
[Output prefix]"
		exit 1
	fi	
	
	haplo_geno_matbed=$2
	sample_ids=$3
	upsampled_sample_size=$4
	genetic_maps_dir=$5
	N_e=$6
	eps_allele=$7
	l_seg_in_bps=$8
	l_seg_in_cMs=$9
	l_seg_in_n_vars=${10}
	n_threads=${11}
	op_prefix=${12}
	
	# This is fixed for the chromosome.
	start_pos=0
	end_pos=250000000
	
	if [[ ! -f ${haplo_geno_matbed} ]]
	then
		echo "Could not find ${haplo_geno_matbed}"
		exit 1
	fi
	
	if [[ ! -f ${sample_ids} ]]
	then
		echo "Could not find ${sample_ids}"
		exit 1
	fi
	
	if [[ ! -d ${genetic_maps_dir} ]]
	then
		echo "Could not find ${genetic_maps_dir}"
		exit 1
	fi
	
	${PROXYTYPER_EXEC} -resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb ${haplo_geno_matbed} ${sample_ids} \
${genetic_maps_dir} ${upsampled_sample_size} ${N_e} ${eps_allele} \
${l_seg_in_bps} ${l_seg_in_cMs} ${l_seg_in_n_vars} \
${n_threads} ${start_pos} ${end_pos} ${op_prefix} >& ${op_prefix}_RESAMPLING.OP

	# Summarize resampling patterns.
	${PROXYTYPER_EXEC} -summarize_sampled_segments_per_resampled_haplotype_info ${op_prefix}_resampling_pattern_signal_var_regs.sigbed.gz \
${op_prefix}_resampled_sample_ids.list ${sample_ids} \
${op_prefix}_resampled_segments.bed

	ret_res=$?
	
	exit ${ret_res}
fi

#
# -impute_site1_per_proxy_protocol : This simulates the imputation among a query and a reference site using proxy panels.
#
if [[ ${cmd_option} == "-impute_site1_per_proxy_protocol" ]]
then
	if [[ $# != 30 ]]
	then
		echo "USAGE: $0 $1 [Chromosome ID] [Assembly ID (e.g., hg19)] 
[Site1 haplocoded tag geno matbed] [Site1 haplocoded target geno matbed] [Site1 sample list] 
[Site1 presampled tag geno matbed file] [Site1 presampled sample list]
[Site2 haplocoded tag geno matbed] [Site2 haplocoded target geno matbed] [Site2 sample list]
[Site1 presampled tag geno matbed file] [Site1 presampled target geno matbed file] [Site2 presampled sample list]
[Filter Proxy. Vicinity size] [Permutation Proxy. Vicinity size] [Coding modulus] [Allele error epsilon] 
[First order parameter prob.] [Second order parameter prob.] [Third order parameter prob.] [Filter Proxy. minimum # of weights per variant]
[Normalized N_e for mixing model (higher: More local models)] [Filter weight Inversion Probability] \
[Permute proxy genotype inversion probability]
[Coordinate anonymization cM-noise SD] \
[Untyped variant permutation window length (in variants)] [Untyped variant genotype switch prob.] \
[Minimum alternate AF for partitioned variants] \
[Genetic maps directory]"
		exit 1
	fi

	chr_id=$2
	ASSEMBLY_ID=$3

	SITE1_haplocoded_tag_matbed=$4
	SITE1_haplocoded_target_matbed=$5
	SITE1_sample_list=$6

	SITE1_presampled_haplocoded_tag_matbed=$7
	SITE1_presampled_sample_list=$8

	SITE2_haplocoded_tag_matbed=$9
	SITE2_haplocoded_target_matbed=${10}
	SITE2_sample_list=${11}

	SITE2_presampled_haplocoded_tag_matbed=${12}
	SITE2_presampled_haplocoded_target_matbed=${13}
	SITE2_presampled_sample_list=${14}

	filter_proxization_vicinity_size=${15}
	perm_proxization_vicinity_size=${16}
	coding_modulus=${17}
	allele_err_eps=${18}

	var_weight_prob=${19}
	var2var_interaction_prob=${20}
	var2var2var_interaction_prob=${21}

	filter_proxization_min_n_params_per_var=${22}

	normalized_N_e=${23}

	filter_weight_inversion_prob=${24}
	perm_proxy_geno_inversion_prob=${25}

	# Amount of noise we add to the "anonymized" genetic map.
	coord_anon_cM_noise_SD=${26}

	untyped_var_perm_n_vicinity=${27}
	untyped_var_allele_switch_prob=${28}
	variant_decomp_min_AAF=${29}

	per_chrom_maps_dir=${30}

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	# Proxies are generated by:
	# resample(G) -> G1, G2 where G1 is the typed variant haplotypes ; G2 -> untyped variant genotypes (Only for Site2, i.e., reference site).
	# G1 -> f(G1; w)	-> anon_coord_map(f(G1; w); C, D) 
	# G2				-> anon_coord_map(G2)     ; C, D)	-> permute_per_target_blocks((anon_coord_map(G2); C, D), P) -> partition_genotypes_per_target_blocks(...)

	########################################################################################################
	# Generate the parameters using LD parameters. 
	${PROXYTYPER_EXEC} -generate_save_per_site_mixing_parameters_LD_aware ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

	if [[ ! -f ${weight_params_file} ]]
	then
		echo "Could not generate mixing proxy model file @ ${weight_params_file}"
		exit 1
	fi
	
	# Generate the permutation proxization mapping regions.
	${PROXYTYPER_EXEC} -generate_permute_proxizing_parameters ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${perm_proxization_vicinity_size} \
${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

	if [[ ! -f ${permute_proxying_mapping_BED} ]]
	then
		echo "Could not generate permutation proxy mapper file @ ${permute_proxying_mapping_BED}"
		exit 1
	fi

	####################################################################################################
	# Site1 setup: Copy resampled tags to correct naming.
	if [[ ! -f ${SITE1_presampled_haplocoded_tag_matbed} ]]
	then
		echo "Could not find ${SITE1_presampled_haplocoded_tag_matbed}"
		exit 1
	fi

	if [[ ! -f ${SITE1_presampled_sample_list} ]]
	then
		echo "Could not find ${SITE1_presampled_sample_list}"
		exit 1
	fi

	cp ${SITE1_presampled_sample_list} site1_resampled_panel_sample_ids.list
	cp ${SITE1_presampled_haplocoded_tag_matbed} site1_resampled_panel_tag_haplocoded.matbed.gz
	####################################################################################################

	# Proxize the site1 data.
	allele_err_eps=0

	# Permute proxize the site1 data.
	${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions site1_resampled_panel_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${permute_proxying_mapping_BED} site1_resampled_panel_perm_prox.matbed.gz

	# Filter proxize site1 variants.
	allele_err_eps=0
	${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_resampled_panel_perm_prox.matbed.gz site1_resampled_panel_sample_ids.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_resampled_panel_tag_haplocoded.matbed.gz

	####################################################################################################
	# Site2 setup: Copy the presampled tag/targets to correct naming.
	if [[ ! -f ${SITE2_presampled_haplocoded_tag_matbed} ]]
	then
		echo "Could not find ${SITE2_presampled_haplocoded_tag_matbed}"
		exit 1
	fi

	if [[ ! -f ${SITE2_presampled_haplocoded_target_matbed} ]]
	then
		echo "Could not find ${SITE2_presampled_haplocoded_target_matbed}"
		exit 1
	fi

	if [[ ! -f ${SITE2_presampled_sample_list} ]]
	then
		echo "Could not find ${SITE2_presampled_sample_list}"
		exit 1
	fi

	cp ${SITE2_presampled_haplocoded_tag_matbed} SITE2_resampled_tag_target_panel_tags.matbed.gz
	cp ${SITE2_presampled_haplocoded_target_matbed} SITE2_resampled_tag_target_panel_targets.matbed.gz
	cp ${SITE2_presampled_sample_list} site2_resampled_panel_sample_ids.list

	SITE2_haplocoded_resampled_panel_tag_matbed=SITE2_resampled_tag_target_panel_tags.matbed.gz
	SITE2_haplocoded_resampled_panel_target_matbed=SITE2_resampled_tag_target_panel_targets.matbed.gz
	#SITE2_haplocoded_resampled_panel_tag_matbed=${SITE2_presampled_haplocoded_tag_matbed}
	#SITE2_haplocoded_resampled_panel_target_matbed=${SITE2_presampled_haplocoded_target_matbed}
	#cp ${SITE2_presampled_sample_list} site2_resampled_panel_sample_ids.list
	####################################################################################################

	# Do imputation using site1 sample at site 2.
	allele_err_eps=0
	# Permute proxize the site1 data.
	${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions ${SITE2_haplocoded_resampled_panel_tag_matbed} site2_resampled_panel_sample_ids.list ${permute_proxying_mapping_BED} site2_resampled_panel_perm_prox.matbed.gz

	# Filter proxize site1 variants.
	allele_err_eps=0
	${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site2_resampled_panel_perm_prox.matbed.gz site2_resampled_panel_sample_ids.list ${weight_params_file} ${allele_err_eps} ${N_THREADS} site2_proxized_resampled_panel_tag_haplocoded.matbed.gz

	#####################################################################
	# Coordinate anonymization: This applies to both tag and target variants. Note that this protects the variant positions and the genetic map.
	anon_genetic_maps_dir=anon_coords_data
	beagle_genetic_map_file=${anon_genetic_maps_dir}/anon_beagle_map.map

	# Dump the original tag and target variants.
    ${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} 1 tag_vars.bed
    ${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} 1 target_vars.bed

    #10      .       0.000000        72765
    rm -f -r ${anon_genetic_maps_dir}
    mkdir ${anon_genetic_maps_dir}
    ${PROXYTYPER_EXEC} -anonymize_tag_target_genetic_map_coords tag_vars.bed target_vars.bed ${per_chrom_maps_dir}/${chr_id}.map ${coord_anon_cM_noise_SD} tag_target_mapping ${anon_genetic_maps_dir}

	# Writet the beagle formatted genetic map file.
    awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${anon_genetic_maps_dir}/${chr_id}.map > ${beagle_genetic_map_file} 

	# Map the coordinates for both site1 and site2.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs site1_proxized_resampled_panel_tag_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz

	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs site2_proxized_resampled_panel_tag_haplocoded.matbed.gz site2_resampled_panel_sample_ids.list ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site2_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz

	# Anonymize coordinates of target variants. Note that these are further permuted below.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs ${SITE2_haplocoded_resampled_panel_target_matbed} site2_resampled_panel_sample_ids.list ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site2_resampled_target_anon_coord_haplocoded.matbed.gz

	#####################################################################
	# Proxize the target variants by permutation: Note that this is done only at SITE2.
    #untyped_recoding_op_prefix="untyped_proxy"
	untyped_recoding_op_prefix="site2_resampled_panel_target_anon_coord"
	untyped_proxy_2_target_mapping_BED=${untyped_recoding_op_prefix}_target_proxy_mapping.bed
	untyped_proxized_haplocoded_target_matbed=${untyped_recoding_op_prefix}_proxized_targets.matbed.gz

	# Recode the targets for the site2: Use the anon-coord for this since we are operating on anonymized coordinates.
	${PROXYTYPER_EXEC} -recode_untyped_variant_reference_panel_per_target_permutation \
site2_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz site2_resampled_target_anon_coord_haplocoded.matbed.gz \
site2_resampled_panel_sample_ids.list \
${untyped_var_perm_n_vicinity} \
${untyped_var_allele_switch_prob} \
${untyped_recoding_op_prefix}

	# Check the files.
	if [[ ! -f ${untyped_proxy_2_target_mapping_BED} ]]
	then
		echo "Could not find encoding coordinates @ \"${untyped_proxy_2_target_mapping_BED}\""
		exit 1
	fi

	if [[ ! -f ${untyped_proxized_haplocoded_target_matbed} ]]
	then
		echo "Could not find encoded target genotypes @ \"${untyped_proxized_haplocoded_target_matbed}\""
		exit 1
	fi

	##########################################################################################
	# Do partitioning of the untyped variants.
	# Shuffling is always hardcoded to increase privacy.
	shuffle_decomp_vars=1

	site2_untyped_decomp_prefix="site2_resampled_panel_target_anon_coord"

	# Decompose the current untyped variants, this does not change the typed variants, which are proxied above.
	${PROXYTYPER_EXEC} -simple_decompose_untyped_variants site2_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz \
${untyped_proxized_haplocoded_target_matbed} \
site2_resampled_panel_sample_ids.list \
${variant_decomp_min_AAF} \
${shuffle_decomp_vars} ${site2_untyped_decomp_prefix}

	untyped_decomposing_interval=${site2_untyped_decomp_prefix}_untyped_decomposing.interval
	untyped_decomposed_target_matbed=${site2_untyped_decomp_prefix}_decomposed.matbed.gz

	if [[ ! -f ${untyped_decomposing_interval} ]]
	then
		echo "Could not find the decomposing interval file @ ${untyped_decomposing_interval}"
		exit 1
	fi

	if [[ ! -f ${untyped_decomposed_target_matbed} ]]
	then
		echo "Could not find the decompose untyped genotypes file @ ${untyped_decomposed_target_matbed}"
		exit 1
	fi

	##########################################################################################
	# Following takes place at the central server.
	##########################################################################################

	# Do the first phase of imputation at the central server.
	rm -f -r resampled_proxized_impute_beagle_dir
	mkdir resampled_proxized_impute_beagle_dir

	# We send the coordinate anonymized tags from both sites and permutation-recoded target list from site2 as the reference.
	$0 -setup_BEAGLE_files site1_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list \
site2_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz ${untyped_decomposed_target_matbed} site2_resampled_panel_sample_ids.list \
${USE_PHASED_SITE1_GENO} resampled_proxized_impute_beagle_dir

	# This performs the imputation using the resampled-proxized-anon-coord-panels and the anonymized genetic map.
	$0 -run_BEAGLE ${beagle_genetic_map_file} resampled_proxized_impute_beagle_dir

	#####################################################################
	#####################################################################
	# Site 1 downloads the panel, and first combines the decomposed untyped variants.
	${PROXYTYPER_EXEC} -combine_BEAGLE_imputed_decomposed_genotype_probabilities resampled_proxized_impute_beagle_dir/imputed.op.vcf.gz \
${untyped_decomposing_interval} site1_resampled_panel_sample_ids.list decomp
	#decomp_untyped_undecomp.matbed.gz
	cp decomp_untyped_undecomp.matbed.gz SITE1_imputed_resampled_proxy_panel_target_haplocoded.matbed.gz

	## Now, we need to extract the tag variants, but we already have these in a previous file.
	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions site1_proxized_resampled_panel_tag_anon_coord_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list 1 anon_coord_tag_vars.bed
	gzip -cd resampled_proxized_impute_beagle_dir/imputed.op.vcf.gz | awk -v chr_id=${chr_id} {'if($1==chr_id && length($4)==1 && length($5)==1)print $0'} | ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF stdin site1_resampled_panel_sample_ids.list anon_coord_tag_vars.bed ${ASSEMBLY_ID}.list ${GENOME_SEQ_DIR} 0 0 1 SITE1_imputed_resampled_proxy_panel_tag_haplocoded.matbed.gz

	## Extract the tag and target variants: These are the coordinate proxized tags and untyped permuted coordinates.
	# These are used below to extract variants from imputed site1 variants.
	cut -f1,2,3,6 ${untyped_proxy_2_target_mapping_BED} > proxized_anon_coord_target_vars.bed
	cat anon_coord_tag_vars.bed proxized_anon_coord_target_vars.bed > tag_proxy_target_vars.bed

	# This is the new imputation panel for site1 now.

	#####################################################################
	# Moving to 2nd imputation step: Site1 proxizes variants only and does imputation.

	#####################
	# Now do imputation at site 1 using the imputed meta-panel from site2..
	# First, recode the site1 genotype data.
	${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${permute_proxying_mapping_BED} site1_panel_perm_prox.matbed.gz

	# Filter proxize site1 variants.
	allele_err_eps=0
	${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_panel_perm_prox.matbed.gz ${SITE1_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_tag_haplocoded.matbed.gz

	# Map the genetic variants to the coordinate anonymized coordinates.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs site1_proxized_tag_haplocoded.matbed.gz ${SITE1_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_tag_anon_coord_haplocoded.matbed.gz
	#######################################################################################

	rm -f -r site1_impute_per_imputed_resampled_proxy_panel
	mkdir site1_impute_per_imputed_resampled_proxy_panel
	$0 -setup_BEAGLE_files site1_proxized_tag_anon_coord_haplocoded.matbed.gz ${SITE1_sample_list} SITE1_imputed_resampled_proxy_panel_tag_haplocoded.matbed.gz SITE1_imputed_resampled_proxy_panel_target_haplocoded.matbed.gz site1_resampled_panel_sample_ids.list ${USE_PHASED_SITE1_GENO} site1_impute_per_imputed_resampled_proxy_panel

	$0 -run_BEAGLE ${beagle_genetic_map_file} site1_impute_per_imputed_resampled_proxy_panel

	fetchChromSizes ${ASSEMBLY_ID} > ${ASSEMBLY_ID}.list
	gzip -cd site1_impute_per_imputed_resampled_proxy_panel/imputed.op.vcf.gz | awk -v chr_id=${chr_id} {'if($1==chr_id && length($4)==1 && length($5)==1)print $0'} | ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF stdin ${SITE1_sample_list} tag_proxy_target_vars.bed ${ASSEMBLY_ID}.list ${GENOME_SEQ_DIR} 0 0 1 site1_imputed_proxy_panel_tag_target_haplocoded.matbed.gz

	# Finally, extract the targets from the imputed values.
	${PROXYTYPER_EXEC} -extract_genotype_signals_per_region_list site1_imputed_proxy_panel_tag_target_haplocoded.matbed.gz ${SITE1_sample_list} proxized_anon_coord_target_vars.bed site1_imputed_proxy_panel_target_haplocoded.matbed.gz

	#######################################################################################
	# Decode the target variants from permuted ordering to original ordering.
	# Decode the untyped variants back to original coordinates: These are now in anonymized coordinates.
	${PROXYTYPER_EXEC} -decode_untyped_variant_reference_panel_per_target_permutation site1_imputed_proxy_panel_target_haplocoded.matbed.gz ${untyped_proxy_2_target_mapping_BED} ${SITE1_sample_list} decoded_site1_imputed_proxy_panel_target_haplocoded.matbed.gz

	#######################################################################################
	# Now, de-anonymize the coordinates of all variants; after this, we should be back on the original target coordinates.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs decoded_site1_imputed_proxy_panel_target_haplocoded.matbed.gz ${SITE1_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed SITE1_impute_cleartext_target_haplocoded.matbed.gz

	#######################################################################################
	# At this point, we should be done with the original coordinates for SITE1's untyped variants.
	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded SITE1_impute_cleartext_target_haplocoded.matbed.gz ${SITE1_sample_list} SITE1_impute_cleartext_target_genocoded.matbed.gz

	# Extract the known genotypes for SITE1.
	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} SITE1_target_genocoded.matbed.gz

	$0 -get_R2_per_imputed_genotypes SITE1_impute_cleartext_target_genocoded.matbed.gz ${SITE1_sample_list} SITE1_target_genocoded.matbed.gz ${SITE1_sample_list} ${chr_id}_meta_accs.txt

	exit 0
fi


if [[ ${cmd_option} == "-impute_site1_per_proxized_variants" ]]
then
	if [[ $# != 24 ]]
	then
		echo "USAGE: $0 $1 [Chromosome ID] [Assembly ID (e.g., hg19)] [Site1 haplocoded tag geno matbed] [Site1 haplocoded target geno matbed] [Site1 sample list] \
[Site2 haplocoded tag geno matbed] [Site2 haplocoded target geno matbed] [Site2 sample list] \
[Filter Proxy Vicinity size] [Permutation Proxy. Vicinity size] \
[Coding modulus] [Allele error epsilon] \
[First order parameter prob.] [Second order parameter prob.] [Third order parameter prob.] \
[Filter Proxy minimum # of weights per variant] \
[Normalized N_e for mixing model (higher: More local models)] [Filter weight Inversion Probability] [Permute proxy genotype inversion probability] \
[Coordinate anonymization cM-noise SD] \
[Untyped variant permutation window length (in variants)] \
[Untyped variant genotype switch prob.] \
[Genetic maps directory]"
		exit 1
	fi

	chr_id=$2
	ASSEMBLY_ID=$3
	SITE1_haplocoded_tag_matbed=$4
	SITE1_haplocoded_target_matbed=$5
	SITE1_sample_list=$6

	SITE2_haplocoded_tag_matbed=$7
	SITE2_haplocoded_target_matbed=$8
	SITE2_sample_list=$9

	filter_proxization_vicinity_size=${10}
	perm_proxization_vicinity_size=${11}
	coding_modulus=${12}
	allele_err_eps=${13}

	var_weight_prob=${14}
	var2var_interaction_prob=${15}
	var2var2var_interaction_prob=${16}

	filter_proxization_min_n_params_per_var=${17}

	normalized_N_e=${18}

	filter_weight_inversion_prob=${19}
	perm_proxy_geno_inversion_prob=${20}

	# Amount of noise we add to the "anonymized" genetic map.
	coord_anon_cM_noise_SD=${21}

	untyped_var_perm_n_vicinity=${22}
	untyped_var_allele_switch_prob=${23}

	per_chrom_maps_dir=${24}

	# This selects whether site1 sends phased (haplocoded) or unphased (genocoded) genotypes.
	# Note that for proxization requires site1 genotypes to be phased.
	USE_PHASED_SITE1_GENO=0

	weight_params_file=per_var_proxy.params
	permute_proxying_mapping_BED=proxizing_permute_reg_pairs.bed

	########################################################################################################
	# Generate the parameters using LD parameters. 
	${PROXYTYPER_EXEC} -generate_save_per_site_mixing_parameters_LD_aware ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

	# Generate the permutation proxization mapping regions.
	${PROXYTYPER_EXEC} -generate_permute_proxizing_parameters ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${perm_proxization_vicinity_size} ${perm_proxy_geno_inversion_prob} ${permute_proxying_mapping_BED}

	####################################################################################################
	# Site1 setup: Proxize tags.

	# Permute proxize the site1 data.
	${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} ${permute_proxying_mapping_BED} site1_perm_prox.matbed.gz

	# Filter proxize site1 variants.
	allele_err_eps=0
	${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site1_perm_prox.matbed.gz ${SITE1_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site1_proxized_tag_haplocoded.matbed.gz

	#####################################################################
	# Site2 setup: Proxize tags, proxize targets.
	# Permute proxize site2 variants.
	${PROXYTYPER_EXEC} -permute_proxize_genotype_signal_regions ${SITE2_haplocoded_tag_matbed} ${SITE2_sample_list} ${permute_proxying_mapping_BED} site2_perm_prox.matbed.gz
	
	# Filter proxize site2 variants.
	allele_err_eps=0
	${PROXYTYPER_EXEC} -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters site2_perm_prox.matbed.gz ${SITE2_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} site2_proxized_tag_haplocoded.matbed.gz

	#####################################################################
	# Coordinate anonymization.
	anon_genetic_maps_dir=anon_coords_data
	beagle_genetic_map_file=${anon_genetic_maps_dir}/anon_beagle_map.map

	# Dump the original tag and target variants.
    ${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${SITE1_haplocoded_tag_matbed} ${SITE1_sample_list} 1 tag_vars.bed
    ${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} 1 target_vars.bed

    #10      .       0.000000        72765
    rm -f -r ${anon_genetic_maps_dir}
    mkdir ${anon_genetic_maps_dir}
    ${PROXYTYPER_EXEC} -anonymize_tag_target_genetic_map_coords tag_vars.bed target_vars.bed ${per_chrom_maps_dir}/${chr_id}.map ${coord_anon_cM_noise_SD} tag_target_mapping ${anon_genetic_maps_dir}

	# Writet the beagle formatted genetic map file.
    awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${anon_genetic_maps_dir}/${chr_id}.map > ${beagle_genetic_map_file} 

	# Map the coordinates for both site1 and site2.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs site1_proxized_tag_haplocoded.matbed.gz ${SITE1_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site1_proxized_tag_anon_coord_haplocoded.matbed.gz

	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs site2_proxized_tag_haplocoded.matbed.gz ${SITE2_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site2_proxized_tag_anon_coord_haplocoded.matbed.gz

	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs ${SITE2_haplocoded_target_matbed} ${SITE2_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed site2_target_anon_coord_haplocoded.matbed.gz

	#####################################################################
	# Proxize the target variants by recoding: Note that this is done only at SITE2.
    untyped_recoding_op_prefix="untyped_proxy"
	untyped_proxy_2_target_mapping_BED=${untyped_recoding_op_prefix}_target_proxy_mapping.bed
	untyped_proxized_haplocoded_target_matbed=${untyped_recoding_op_prefix}_proxized_targets.matbed.gz

	# Recode the targets for the site2: Use the anon-coord for this since we are operating on anonymized coordinates.
	${PROXYTYPER_EXEC} -recode_untyped_variant_reference_panel_per_target_permutation site2_proxized_tag_anon_coord_haplocoded.matbed.gz \
site2_target_anon_coord_haplocoded.matbed.gz \
${SITE2_sample_list} \
${untyped_var_perm_n_vicinity} \
${untyped_var_allele_switch_prob} \
${untyped_recoding_op_prefix}

	# Check the files.
	if [[ ! -f ${untyped_proxy_2_target_mapping_BED} ]]
	then
		echo "Could not find encoding coordinates @ \"${untyped_proxy_2_target_mapping_BED}\""
		exit 1
	fi

	if [[ ! -f ${untyped_proxized_haplocoded_target_matbed} ]]
	then
		echo "Could not find encoded target genotypes @ \"${untyped_proxized_haplocoded_target_matbed}\""
		exit 1
	fi

	##########################################################################################
	# Following checks the untyped variant encoding.
	ENCDEC_CHECK=0
	if [[ $ENCDEC_CHECK == 1 ]]
	then
		${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${untyped_proxized_haplocoded_target_matbed} ${SITE2_sample_list} 0 ${untyped_proxized_haplocoded_target_matbed}.txt
		${PROXYTYPER_EXEC} -decode_untyped_variant_reference_panel_per_target_permutation ${untyped_proxized_haplocoded_target_matbed} ${untyped_proxy_2_target_mapping_BED} ${SITE2_sample_list} decoded_proxy.matbed.gz
		${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions decoded_proxy.matbed.gz ${SITE2_sample_list} 0 decoded_proxy.matbed.gz.txt
		${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions site2_target_anon_coord_haplocoded.matbed.gz ${SITE2_sample_list} 0 original_genotypes.txt
		sort -n -k2,2 decoded_proxy.matbed.gz.txt > sorted_decoded_proxy.matbed.gz.txt

		echo "##################################################################################"
		echo "Difference between original and decoded target variants: (Must be empty below)"
		echo "##################################################################################"
		diff sorted_decoded_proxy.matbed.gz.txt original_genotypes.txt
		echo "##################################################################################"

		exit 0
	fi
	##########################################################################################


	#####################################################################

	# Impute the site1's proxy panel.
	rm -f -r proxized_impute_beagle_dir
	mkdir proxized_impute_beagle_dir

	# We send the coordinate anonymized tags from both sites and permutation-recoded target list from site2 as the reference.
	$0 -setup_BEAGLE_files site1_proxized_tag_anon_coord_haplocoded.matbed.gz ${SITE1_sample_list} site2_proxized_tag_anon_coord_haplocoded.matbed.gz ${untyped_proxized_haplocoded_target_matbed} ${SITE2_sample_list} ${USE_PHASED_SITE1_GENO} proxized_impute_beagle_dir

	$0 -run_BEAGLE ${beagle_genetic_map_file} proxized_impute_beagle_dir

	#####################################################################
	# Extract the tag and target coordinates: These are the final recoded coordinates, not the anon-coords from first stage.
	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions site2_proxized_tag_anon_coord_haplocoded.matbed.gz ${SITE2_sample_list} 1 anon_coord_tag_vars.bed
	#variation_tools -dump_plain_geno_signal_regions site2_target_anon_coord_haplocoded.matbed.gz ${SITE2_sample_list} 1 proxized_anon_coord_target_vars.bed
	cut -f1,2,3,6 ${untyped_proxy_2_target_mapping_BED} > proxized_anon_coord_target_vars.bed
	#####################################################################

	cat anon_coord_tag_vars.bed proxized_anon_coord_target_vars.bed > tag_proxy_target_vars.bed

	fetchChromSizes ${ASSEMBLY_ID} > ${ASSEMBLY_ID}.list
	gzip -cd proxized_impute_beagle_dir/imputed.op.vcf.gz | awk -v chr_id=${chr_id} {'if($1==chr_id && length($4)==1 && length($5)==1)print $0'} | ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF stdin ${SITE1_sample_list} tag_proxy_target_vars.bed ${ASSEMBLY_ID}.list ${GENOME_SEQ_DIR} 0 0 1 SITE1_proxized_imputed_haplocoded.matbed.gz

	# Extract the proxized target variant coordinates genotypes.
	${PROXYTYPER_EXEC} -extract_genotype_signals_per_region_list SITE1_proxized_imputed_haplocoded.matbed.gz ${SITE1_sample_list} proxized_anon_coord_target_vars.bed SITE1_imputed_proxy_target_haplocoded.matbed.gz

	# Decode the untyped variants back to original coordinates: These are now in anonymized coordinates.
	${PROXYTYPER_EXEC} -decode_untyped_variant_reference_panel_per_target_permutation SITE1_imputed_proxy_target_haplocoded.matbed.gz ${untyped_proxy_2_target_mapping_BED} ${SITE1_sample_list} decoded_SITE1_imputed_target.matbed.gz

	#######################################################################################
	# Now, de-anonymize the coordinates of all variants; we should be in original target coordinates.
	${PROXYTYPER_EXEC} -generic_map_genotype_regs_per_src_dest_regs decoded_SITE1_imputed_target.matbed.gz ${SITE1_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed SITE1_impute_cleartext_target_haplocoded.matbed.gz

	# Convert to genocoded coordinates.
	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded SITE1_impute_cleartext_target_haplocoded.matbed.gz ${SITE1_sample_list} SITE1_impute_cleartext_target_genocoded.matbed.gz

	# Extract the known genotypes for SITE1 for calculating accuracy.
	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${SITE1_haplocoded_target_matbed} ${SITE1_sample_list} SITE1_target_genocoded.matbed.gz

	# Get accuracy statistics.
	$0 -get_R2_per_imputed_genotypes SITE1_impute_cleartext_target_genocoded.matbed.gz ${SITE1_sample_list} SITE1_target_genocoded.matbed.gz ${SITE1_sample_list} ${chr_id}_meta_accs.txt

	# We should be done at this point.
	exit 0
fi

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

if [[ ${cmd_option} == "-import_VCF" ]]
then
	if [[ $# -lt 3 ]]
	then
		echo "$0 $1 [VCF file path] [Assembly ID] [Output file]" >&2
		exit 1
	fi

	vcf_file=$2
	ASSEMBLY_ID=$3
	output_file=$4

	GENOME_SEQ_DIR="none"

	fetchChromSizes ${ASSEMBLY_ID} > ${ASSEMBLY_ID}.list

	echo "Extracting VCF subject identifiers."
	$0 -cat2stdout $vcf_file | head -n 1000 | awk {'if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i}}'} > ${vcf_file}_sample_ids.list

	echo "Extracting VCF file variants."
	$0 -cat2stdout $vcf_file | cut -f1-10 | awk {'ref_all=$4;alt_all=$5;if(length(ref_all)!=1 || length(alt_all)!=1){next;};if(NR%10000==0){printf("@ %d. variant..\r", NR) > "/dev/stderr";};chr_id=$1;start_pos=$2-1;end_pos=$2;snp_id=$3"_"ref_all"_"alt_all;print chr_id"\t"start_pos"\t"end_pos"\t"snp_id"\t.\t+"'} > ${vcf_file}_var_regs.bed

	echo "Importing VCF file genotypes."
	$0 -cat2stdout $vcf_file | awk {'ref_all=$4;alt_all=$5;if(length(ref_all)!=1 || length(alt_all)!=1){next;};print $0'} | ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF stdin ${vcf_file}_sample_ids.list ${vcf_file}_var_regs.bed ${ASSEMBLY_ID}.list ${GENOME_SEQ_DIR} 0 0 1 ${output_file}

	exit 0
fi

if [[ ${cmd_option} == "-extract_tag_target_genotypes" ]]
then
	if [[ $# -lt 6 ]]
	then
		echo "$0 $1 [Panel haplocoded tag genotypes matbed] [Panel haplocoded target genotypes matbed] [Panel sample ids] [Sub. sample ids] [Output prefix]" >&2
		exit 1
	fi

	panel_haplo_tag_matbed=$2
	panel_haplo_target_matbed=$3
	panel_sample_ids=$4
	sub_sample_ids=$5
	op_prefix=$6

	rm -f ${op_prefix}_*

	${PROXYTYPER_EXEC} -extract_genotype_signals_per_subsample_list ${panel_haplo_tag_matbed} ${panel_sample_ids} ${sub_sample_ids} ${op_prefix}_tag_haplocoded.matbed.gz

	if [[ ! -f ${op_prefix}_tag_haplocoded.matbed.gz ]]
	then
		echo "Could not generate ${op_prefix}_tag_haplocoded.matbed.gz"
		exit 1
	fi

	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${op_prefix}_tag_haplocoded.matbed.gz ${sub_sample_ids} ${op_prefix}_tag_genocoded.matbed.gz
	if [[ ! -f ${op_prefix}_tag_genocoded.matbed.gz ]]
	then
		echo "Could not generate ${op_prefix}_tag_genocoded.matbed.gz"
		exit 1
	fi

	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${op_prefix}_tag_genocoded.matbed.gz ${sub_sample_ids} 0 ${op_prefix}_tag_genocoded.matbed.gz.txt
	if [[ ! -f ${op_prefix}_tag_genocoded.matbed.gz.txt ]]
	then
		echo "Could not generate ${op_prefix}_tag_genocoded.matbed.gz.txt"
		exit 1
	fi

	sort -n -k2,2 ${op_prefix}_tag_genocoded.matbed.gz.txt > ${op_prefix}_sorted_tag_genocoded.matbed.gz.txt
	if [[ ! -f ${op_prefix}_sorted_tag_genocoded.matbed.gz.txt ]]
	then
		echo "Could not generate ${op_prefix}_sorted_tag_genocoded.matbed.gz.txt"
		exit 1
	fi

	${PROXYTYPER_EXEC} -extract_genotype_signals_per_subsample_list ${panel_haplo_target_matbed} ${panel_sample_ids} ${sub_sample_ids} ${op_prefix}_target_haplocoded.matbed.gz
	if [[ ! -f ${op_prefix}_target_haplocoded.matbed.gz ]]
	then
		echo "Could not generate ${op_prefix}_target_haplocoded.matbed.gz"
		exit 1
	fi

	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${op_prefix}_target_haplocoded.matbed.gz ${sub_sample_ids} ${op_prefix}_target_genocoded.matbed.gz
	if [[ ! -f ${op_prefix}_target_genocoded.matbed.gz ]]
	then
		echo "Could not generate ${op_prefix}_target_genocoded.matbed.gz"
		exit 1
	fi

	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${op_prefix}_target_genocoded.matbed.gz ${sub_sample_ids} 0 ${op_prefix}_target_genocoded.matbed.gz.txt
	if [[ ! -f ${op_prefix}_target_genocoded.matbed.gz.txt ]]
	then
		echo "Could not generate ${op_prefix}_target_genocoded.matbed.gz.txt"
		exit 1
	fi

	sort -n -k2,2 ${op_prefix}_target_genocoded.matbed.gz.txt > ${op_prefix}_sorted_target_genocoded.matbed.gz.txt
	if [[ ! -f ${op_prefix}_sorted_target_genocoded.matbed.gz.txt ]]
	then
		echo "Could not generate ${op_prefix}_sorted_target_genocoded.matbed.gz.txt"
		exit 1
	fi

	exit 0
fi


if [[ ${cmd_option} == "-get_R2_per_imputed_genotypes" ]]
then
	if [[ $# -lt 6 ]]
	then
		echo "$0 $1 [Imputed genotypes matbed] [Imputed sample ids] [Known genotypes matbed] [known sample ids] [Output file]" >&2
		exit 1
	fi

	imputed_geno_matbed=$2
	imputed_sample_ids=$3
	known_geno_matbed=$4
	known_sample_ids=$5
	op_file=$6

	exec_check=`type -P ${PROXYTYPER_EXEC}`
	if [[ ${exec_check} == "" ]]
	then
		echo "Could not find ${PROXYTYPER_EXEC} executable @ ${PROXYTYPER_EXEC}" >&2
		exit 1
	fi

	${PROXYTYPER_EXEC} -get_R2_per_imputed_genotypes ${imputed_geno_matbed} ${imputed_sample_ids} ${known_geno_matbed} ${known_sample_ids}

	if [[ ! -f "R2_stats.txt" ]]
	then
		echo "Could not find the R2 stats @ R2_stats.txt"
		exit 1
	fi

	mv R2_stats.txt ${op_file}
	exit
fi

if [[ ${cmd_option} == "-setup_BEAGLE_files" ]]
then
	if [[ $# -lt 8 ]]
	then
		echo "$0 $1 [Haplocoded query tag genotypes matbed] [Query sample ids] [Haplocoded reference tag genotype matbed] [Haplocoded reference target genotype matbed] [Reference sample ids] [Use phased(1)/unphased(0) genotypes] [Output directory]"
		exit 1
	fi
	
	query_tag_haplo_matbed=$2
	query_sample_ids=$3
	ref_tag_haplo_matbed=$4
	ref_target_haplo_matbed=$5
	ref_sample_ids=$6
	SAVE_PHASED_FLAG=$7
	op_dir=$8
	
	BEAGLE_INPUT_DATA_DIR=${op_dir}
	
	if [[ ! -d ${BEAGLE_INPUT_DATA_DIR} ]]
	then
		echo "Could not find the output directory for BEAGLE @ \"${BEAGLE_INPUT_DATA_DIR}\""
		exit 1
	fi

	if [[ ! -f ${ref_tag_haplo_matbed} ]]
	then
		echo "Could not find ${ref_tag_haplo_matbed}"
		exit	
	fi

	if [[ ! -f ${ref_target_haplo_matbed} ]]
	then
		echo "Could not find ${ref_target_haplo_matbed}"
		exit
	fi

	if [[ ! -f ${ref_sample_ids} ]]
	then
		echo "Could not find ${ref_sample_ids}"
		exit
	fi

	if [[ ! -f ${query_tag_haplo_matbed} ]]
	then
		echo "Could not find ${query_tag_haplo_matbed}"
		exit
	fi

	if [[ ! -f ${query_sample_ids} ]]
	then
		echo "Could not find ${query_sample_ids}"
		exit
	fi


	 echo "Concatenating tag and target SNVs in the reference panel data."
    #variation_tools -dump_plain_geno_signal_regions ${ref_tag_haplo_matbed} ${ref_sample_ids} 0 ${ref_tag_haplo_matbed}.txt

    #variation_tools -dump_plain_geno_signal_regions ${ref_target_haplo_matbed} ${ref_sample_ids} 0 ${ref_target_haplo_matbed}.txt
    #cat ${ref_tag_haplo_matbed}.txt ${ref_target_haplo_matbed}.txt | sort -u -k2,2 | sort -n -k2,2 > ${BEAGLE_INPUT_DATA_DIR}/sorted_tag_target_ref_haplo.txt

	remove_dup_vars_per_coord=1
    ${PROXYTYPER_EXEC} -concat_variant_wide_genotype_signal_regions ${ref_tag_haplo_matbed} ${ref_sample_ids} ${ref_target_haplo_matbed} ${ref_sample_ids} ${remove_dup_vars_per_coord} ${BEAGLE_INPUT_DATA_DIR}/sorted_tag_target_ref_haplo.matbed.gz

    echo "Writing Beagle input files."
    #[Reference Panel matbed file path] [reference panel sample ids list file path] [Input Study Panel matbed file path] [Input Study panel sample ids list file path] [h_option file path] [l_option file path] [g_option file path] [strand_g_option file path]
    #variation_tools -dump_plain_geno_signal_regions ${query_tag_haplo_matbed} ${query_sample_ids} 0 ${query_tag_haplo_matbed}.txt

    #variation_tools -extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix ${BEAGLE_INPUT_DATA_DIR}/sorted_tag_target_ref_haplo.txt ${ref_sample_ids} ${query_tag_haplo_matbed}.txt ${query_sample_ids} ${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf ${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf ${SAVE_PHASED_FLAG}
    ${PROXYTYPER_EXEC} -extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix ${BEAGLE_INPUT_DATA_DIR}/sorted_tag_target_ref_haplo.matbed.gz ${ref_sample_ids} ${query_tag_haplo_matbed} ${query_sample_ids} ${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf ${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf ${SAVE_PHASED_FLAG}
	
	exit 0
fi

if [[ ${cmd_option} == "-run_BEAGLE" ]]
then
	if [[ $# -lt 3 ]]
	then
		echo "USAGE: $0 $1 [Genetic map file] [BEAGLE data directory]"
		exit 1
	fi

	genetic_map_file=$2
	BEAGLE_INPUT_DATA_DIR=$3
	
	if [[ ! -f ${BEAGLE_JAR} ]]
	then
		echo "Could not find BEAGLE jar file @ \"${BEAGLE_JAR}\""
		exit 1
	fi
	
	if [[ ! -f ${genetic_map_file} ]]
	then
		echo "Could not find genetic map @ ${genetic_map_file}"
		exit 1
	fi

	# Run beagle.
	java -Xmx182044m -jar ${BEAGLE_JAR} nthreads=${N_THREADS} map=${genetic_map_file} ref=${BEAGLE_INPUT_DATA_DIR}/ref_option.gens.vcf gt=${BEAGLE_INPUT_DATA_DIR}/gt_option.gens.vcf gp=true ap=true out=${BEAGLE_INPUT_DATA_DIR}/imputed.op
	ret_rest=$?

	exit ${ret_res}
fi

if [[ ${cmd_option} == "-centralized_beagle_impute" ]]
then
	# Get site1 and site2 as parameters.
	if [[ $# != 10 ]]
	then
		echo "USAGE: $0 $1 [Chromosome identifier] [Assembly ID (e.g., hg19)] [Site1 haplocoded tag geno matbed] [Site1 haplocoded target geno matbed] [Site1 sample list] [Site2 haplocoded tag geno matbed] [Site2 haplocoded target geno matbed] [Site2 sample list] [BEAGLE genetic map file]"
		exit 1
	fi

	chr_id=$2
	ASSEMBLY_ID=$3
	site1_tag_matbed=$4
	site1_target_matbed=$5
	site1_sample_list=$6
	site2_tag_matbed=$7
	site2_target_matbed=$8
	site2_sample_list=$9
	beagle_genetic_map_file=${10}

	rm -f -r central_impute_beagle_dir
	mkdir central_impute_beagle_dir
	$0 -setup_BEAGLE_files ${site1_tag_matbed} ${site1_sample_list} ${site2_tag_matbed} ${site2_target_matbed} ${site2_sample_list} 0 central_impute_beagle_dir
	$0 -run_BEAGLE ${beagle_genetic_map_file} central_impute_beagle_dir

	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${site2_tag_matbed} ${site2_sample_list} 1 tag_vars.bed
	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${site2_target_matbed} ${site2_sample_list} 1 target_vars.bed

	cat tag_vars.bed target_vars.bed > tag_target_vars.bed

	fetchChromSizes ${ASSEMBLY_ID} > ${ASSEMBLY_ID}.list
	gzip -cd central_impute_beagle_dir/imputed.op.vcf.gz | awk -v chr_id=${chr_id} {'if($1==chr_id && length($4)==1 && length($5)==1)print $0'} | ${PROXYTYPER_EXEC} -extract_genotype_signals_per_VCF stdin ${site1_sample_list} tag_target_vars.bed ${ASSEMBLY_ID}.list ${GENOME_SEQ_DIR} 0 0 1 SITE1_central_imputed_haplocoded.matbed.gz

    ${PROXYTYPER_EXEC} -extract_genotype_signals_per_region_list SITE1_central_imputed_haplocoded.matbed.gz ${site1_sample_list} target_vars.bed SITE1_central_imputed_target_haplocoded.matbed.gz

    ${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded SITE1_central_imputed_target_haplocoded.matbed.gz ${site1_sample_list} SITE1_central_imputed_target_genocoded.matbed.gz
	${PROXYTYPER_EXEC} -convert_haplocoded_2_genocoded ${site1_target_matbed} ${site1_sample_list} SITE1_target_genocoded.matbed.gz

	$0 -get_R2_per_imputed_genotypes SITE1_central_imputed_target_genocoded.matbed.gz ${site1_sample_list} SITE1_target_genocoded.matbed.gz ${site1_sample_list} ${chr_id}_central_accs.txt

	exit 0
fi



exit 0
