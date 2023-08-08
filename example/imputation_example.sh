#!/bin/bash

if [[ $# -lt 1 ]]
then
	echo "$0 [Option] [Arguments]
	-clean_directory
	-setup_reference_query_panels
	-resample_query_reference_panels
	-do_imputation"

	exit 1
fi

cmd_option=$1

N_TYPED=16000
N_REF_PANEL=2000
PER_CHROM_MAPS_DIR=../../DATA/beagle_genetic_maps/genetic_maps
PER_CHROM_BEAGLE_MAPS_DIR=../../beagle/genetic_maps/
CHR_ID=22
ASSEMBLY_ID=hg19
PROXYTYPER_EXEC=ProxyTyper_Release

if [[ ! -d ${PER_CHROM_MAPS_DIR} ]]
then
        echo "Could not find genetic maps directory"
        exit 1
fi


if [[ ${cmd_option} == "-clean_directory" ]]
then
	./ProxyTyper.sh -clean_directory ${PWD}

	git clone https://github.com/harmancilab/ProxyTyper.git

	if [[ ! -d "ProxyTyper" ]]
	then
		echo "Could not find the ProxyTyper directory."
		exit 1
	fi

	cp ProxyTyper/scripts/*.sh .
	cp ProxyTyper/example/*.sh .
	dos2unix.sh *.sh
	chmod 755 *.sh

	exit 0
fi

if [[ "${cmd_option}" == "-setup_reference_query_panels" ]]
then
	wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	#./ProxyTyper.sh -import_VCF ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz hg19 ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz

	# Calculate the allele frequencies.
	ProxyTyper_Release -compute_AF_per_geno_signal_regions ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz_AFs.bed

	awk 'BEGIN{FS="\t";OFS="\t"}{$4=$4"_"$5;print $0}'  ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz_AFs.bed > ALL_VAR_REGS.bed

	awk {'if($5>0.05){print $0}'} ALL_VAR_REGS.bed | shuf | head -n ${N_TYPED} > tags.bed

	ProxyTyper_Release -exclude ALL_VAR_REGS.bed tags.bed n
	mv excluded.bed all_untyped.bed

	awk {'if($5>=0.005){print $0}'} all_untyped.bed | shuf | head -n ${N_UNTYPED} > targets.bed

	ProxyTyper_Release -extract_genotype_signals_per_region_list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list tags.bed ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz

	ProxyTyper_Release -extract_genotype_signals_per_region_list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list targets.bed ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz

	shuf ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list | head -n ${N_REF_PANEL} > ref_panel_sample_ids.list

	grep -w -v -f ref_panel_sample_ids.list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list > query_panel_sample_ids.list

	./ProxyTyper.sh -extract_tag_target_genotypes ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list ref_panel_sample_ids.list chr${CHR_ID}_ref_panel

	./ProxyTyper.sh -extract_tag_target_genotypes ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list query_panel_sample_ids.list chr${CHR_ID}_query_panel

	exit 0
fi

if [[ "${cmd_option}" == "-resample_query_reference_panels" ]]
then
	############################################################################################################
	echo "Resampling reference panel typed+untyped variants..."
	Ref_Panel_haplocoded_tag_matbed=chr${CHR_ID}_ref_panel_tag_haplocoded.matbed.gz
	Ref_Panel_haplocoded_target_matbed=chr${CHR_ID}_ref_panel_target_haplocoded.matbed.gz
	ref_panel_sample_list=ref_panel_sample_ids.list

	if [[ ! -f "${Ref_Panel_haplocoded_tag_matbed}" ]]
	then
		echo "Could not find ${Ref_Panel_haplocoded_tag_matbed}"
		exit 1
	fi

	if [[ ! -f "${Ref_Panel_haplocoded_target_matbed}" ]]
	then
		echo "Could not find ${Ref_Panel_haplocoded_target_matbed}"
		exit 1
	fi

	if [[ ! -f "${ref_panel_sample_list}" ]]
	then
		echo "Could not find ${ref_panel_sample_list}"
		exit 1
	fi

	# First resample the reference panel typed+untyped variant genotypes.
	per_chrom_maps_dir=${PER_CHROM_MAPS_DIR}
	resampling_op_prefix=chr${CHR_ID}_ref_panel
	resampled_size=`wc -l ${ref_panel_sample_list} | awk {'print 2*$1'}`
	N_e_frac=0.05
	allele_eps=0
	max_l_seg_n_bps=10000000
	max_l_seg_cM=0
	max_l_seg_nvars=0
	n_threads=20
	start_posn=0
	end_posn=250000000

	echo "Resampling reference typed+untyped panel to ${resampled_size} subjects.."

	./ProxyTyper.sh -resample_tag_target_genotypes_w_length_cutoff ${Ref_Panel_haplocoded_tag_matbed} ${Ref_Panel_haplocoded_target_matbed} ${ref_panel_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${resampling_op_prefix} 

	############################################################################################################
	echo "Resampling query panel typed variants..."

	Query_Panel_haplocoded_tag_matbed=chr${CHR_ID}_query_panel_tag_haplocoded.matbed.gz
	Query_Panel_haplocoded_target_matbed=chr${CHR_ID}_query_panel_target_haplocoded.matbed.gz
	Query_panel_sample_list=query_panel_sample_ids.list

	if [[ ! -f "${Query_Panel_haplocoded_tag_matbed}" ]]
	then
		echo "Could not find ${Query_Panel_haplocoded_tag_matbed}"
		exit 1
	fi

	if [[ ! -f "${Query_Panel_haplocoded_target_matbed}" ]]
	then
		echo "Could not find ${Query_Panel_haplocoded_target_matbed}"
		exit 1
	fi

	if [[ ! -f "${Query_panel_sample_list}" ]]
	then
		echo "Could not find ${Query_panel_sample_list}"
		exit 1
	fi

	per_chrom_maps_dir=${PER_CHROM_MAPS_DIR}
	resampling_op_prefix=chr${CHR_ID}_query_panel
	resampled_size=`wc -l ${Query_panel_sample_list} | awk {'print 2*$1'}`
	N_e_frac=0.05
	allele_eps=0
	max_l_seg_n_bps=10000000
	max_l_seg_cM=0
	max_l_seg_nvars=0
	n_threads=20
	start_posn=0
	end_posn=250000000

	echo "Resampling query typed panel to ${resampled_size} subjects.."

	./ProxyTyper.sh -resample_query_genotypes_w_length_cutoff ${Query_Panel_haplocoded_tag_matbed} ${Query_panel_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${resampling_op_prefix} 

	############################################################################################################

	exit 0 
fi 

if [[ "${cmd_option}" == "-do_imputation" ]]
then

	chr_id=${CHR_ID}

	echo "--------------------------------------"
	echo "Processing chromosome ${chr_id}"
	echo "--------------------------------------"

	per_chrom_maps_dir=${PER_CHROM_MAPS_DIR}

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

	# Weight inversion orobability for the typed variant hashing functions.
	weight_inversion_prob=0.5
	permute_proxy_geno_inv_prob=0.0

	# Genetic map noise standard deviation.
	coord_anon_cM_noise_SD=0.05

	# Untyped variant permutation vicinity size and switching probability.
	untyped_var_perm_n_vicinity=20
	untyped_var_allele_switch_prob=0.5

	# This is the minimum alternate allele frequency at which we will partition the variants among haplotypes.
	variant_decomp_min_AAF=0.1

	weight_params_file=per_var_proxy.params

	normalized_N_e=15

	QUERY_haplocoded_tag_matbed=chr${CHR_ID}_query_panel_tag_haplocoded.matbed.gz
	QUERY_haplocoded_target_matbed=chr${CHR_ID}_query_panel_target_haplocoded.matbed.gz
	QUERY_sample_list=query_panel_sample_ids.list

	QUERY_presampled_tag_matbed=chr${CHR_ID}_query_panel_resampled_tags.matbed.gz
	QUERY_presampled_sample_list=chr${CHR_ID}_query_panel_resampled_sample_ids.list

	REFERENCE_haplocoded_tag_matbed=chr${CHR_ID}_ref_panel_tag_haplocoded.matbed.gz
	REFERENCE_haplocoded_target_matbed=chr${CHR_ID}_ref_panel_target_haplocoded.matbed.gz
	REFERENCE_sample_list=ref_panel_sample_ids.list

	REFERENCE_presampled_tag_matbed=chr${CHR_ID}_ref_panel_resampled_tags.matbed.gz
	REFERENCE_presampled_target_matbed=chr${CHR_ID}_ref_panel_resampled_targets.matbed.gz
	REFERENCE_presampled_sample_list=chr${CHR_ID}_ref_panel_resampled_sample_ids.list

	# Do central imputation on the current subjects.
	RUN_CENTRAL_IMPUTATION=1
	if [[ ${RUN_CENTRAL_IMPUTATION} == 1 ]]
	then
		beagle_genetic_map_file=${PER_CHROM_BEAGLE_MAPS_DIR}/plink.chr${chr_id}.GRCh37.map

		if [[ ! -f ${beagle_genetic_map_file} ]]
		then
			echo "Could not find the genetic map @ \"${beagle_genetic_map_file}\""
			exit 1
		fi

		echo "Performing central imputation without the proxy panels."
		# Calculate the pure site2-based imputation accuracy.
		./ProxyTyper.sh -centralized_beagle_impute ${chr_id} ${ASSEMBLY_ID} ${QUERY_haplocoded_tag_matbed} ${QUERY_haplocoded_target_matbed} ${QUERY_sample_list} \
${REFERENCE_haplocoded_tag_matbed} ${REFERENCE_haplocoded_target_matbed} ${REFERENCE_sample_list} ${beagle_genetic_map_file}
		cp -r central_impute_beagle_dir ${chr_id}_central_impute_beagle_dir
	fi

	./ProxyTyper.sh -impute_site1_per_proxy_protocol ${chr_id} ${ASSEMBLY_ID} \
${QUERY_haplocoded_tag_matbed} ${QUERY_haplocoded_target_matbed} ${QUERY_sample_list} \
${QUERY_presampled_tag_matbed} ${QUERY_presampled_sample_list} \
${REFERENCE_haplocoded_tag_matbed} ${REFERENCE_haplocoded_target_matbed} ${REFERENCE_sample_list} \
${REFERENCE_presampled_tag_matbed} ${REFERENCE_presampled_target_matbed} ${REFERENCE_presampled_sample_list} \
${filter_proxization_vicinity_size} ${perm_proxization_vicinity_size} ${coding_modulus} ${allele_err_eps} \
${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_proxization_min_n_params_per_var} \
${normalized_N_e} ${weight_inversion_prob} \
${permute_proxy_geno_inv_prob} \
${coord_anon_cM_noise_SD} \
${untyped_var_perm_n_vicinity} ${untyped_var_allele_switch_prob} \
${variant_decomp_min_AAF} \
${per_chrom_maps_dir}

	exit 0
fi







