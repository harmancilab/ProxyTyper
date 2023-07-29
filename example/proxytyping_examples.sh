#!/bin/bash

# Setup the parameters that we will use below to select tags/targets, sample sizes, etc.
N_TYPED=16000
N_REF_PANEL=2000
PER_CHROM_MAPS_DIR=../../DATA/beagle_genetic_maps/genetic_maps
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
	dos2unix.sh *.sh
	chmod 755 *.sh

	exit 0
fi

if [[ "${cmd_option}" == "-setup_reference_query_panels" ]]
then
	wget -c https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	./ProxyTyper.sh -import_VCF ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz hg19 ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz

	# Calculate the allele frequencies.
	ProxyTyper_Release -compute_AF_per_geno_signal_regions ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz_AFs.bed

	awk {'if($5>0.05){print $0}'} ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz_AFs.bed | shuf | head -n ${N_TYPED} > tags.bed

	ProxyTyper_Release -exclude ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_var_regs.bed tags.bed n
mv excluded.bed targets.bed

	ProxyTyper_Release -extract_genotype_signals_per_region_list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list tags.bed ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz

	ProxyTyper_Release -extract_genotype_signals_per_region_list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list targets.bed ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz

	shuf ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list | head -n ${N_REF_PANEL} > ref_panel_sample_ids.list

	grep -w -v -f ref_panel_sample_ids.list ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list > query_panel_sample_ids.list

	./ProxyTyper.sh -extract_tag_target_genotypes ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list ref_panel_sample_ids.list chr${CHR_ID}_ref_panel

	./ProxyTyper.sh -extract_tag_target_genotypes ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_tags.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded_targets.matbed.gz ALL.chr${CHR_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_sample_ids.list query_panel_sample_ids.list chr${CHR_ID}_query_panel

	exit 0
fi

# Anonymize the coordinates and genetic maps.
Query_Panel_haplocoded_tag_matbed=chr${CHR_ID}_Query_panel_tag_haplocoded.matbed.gz
Query_Panel_haplocoded_target_matbed=chr${CHR_ID}_Query_panel_target_haplocoded.matbed.gz
Query_panel_sample_list=query_panel_sample_ids.list
per_chrom_maps_dir=${PER_CHROM_MAPS_DIR}
resampling_op_prefix=chr${CHR_ID}_Query_panel
resampled_size=`wc -l ${Query_panel_sample_list} | awk {'print 2*$1'}`
N_e_frac=0.05
allele_eps=0
max_l_seg_n_bps=10000000
max_l_seg_cM=0
max_l_seg_nvars=0
n_threads=20
start_posn=0
end_posn=250000000

if [[ ! -f "${Query_Panel_haplocoded_tag_matbed}" ]]
then
	echo "Could not find extracted tag genotypes."
	exit 1
fi

if [[ ! -f "${Query_Panel_haplocoded_target_matbed}" ]]
then
	echo "Could not find extracted target genotypes."
	exit 1
fi

if [[ ! -f "${Query_panel_sample_list}" ]]
then
	echo "Could not find extracted sample identifiers."
	exit 1
fi

echo "Writing ${resampled_size}"

./ProxyTyper.sh -resample_tag_target_genotypes_w_length_cutoff ${Query_Panel_haplocoded_tag_matbed} ${Query_Panel_haplocoded_target_matbed} ${Query_panel_sample_list} \
${resampled_size} ${per_chrom_maps_dir} ${N_e_frac} ${allele_eps} \
${max_l_seg_n_bps} ${max_l_seg_cM} ${max_l_seg_nvars} \
${n_threads} ${resampling_op_prefix} 

haplocoded_resampled_target_matbed=chr${CHR_ID}_Query_panel_resampled_targets.matbed.gz
haplocoded_resampled_tag_matbed=chr${CHR_ID}_Query_panel_resampled_tags.matbed.gz
resampled_sample_list=chr${CHR_ID}_Query_panel_resampled_sample_ids.list
filter_proxization_vicinity_size=6
var_weight_prob=0.5
var2var_interaction_prob=0.8
var2var2var_interaction_prob=0.3
filter_weight_inversion_prob=0.5
coding_modulus=2
normalized_N_e=15
filter_proxization_min_n_params_per_var=2
per_chrom_maps_dir=${PER_CHROM_MAPS_DIR}
weight_params_file=proxy_generation_parameters.model
ProxyTyper_Release -generate_save_per_site_mixing_parameters_LD_aware ${haplocoded_resampled_tag_matbed} ${resampled_sample_list} \
${filter_proxization_vicinity_size} ${var_weight_prob} ${var2var_interaction_prob} ${var2var2var_interaction_prob} ${filter_weight_inversion_prob} \
${coding_modulus} ${normalized_N_e} ${filter_proxization_min_n_params_per_var} ${per_chrom_maps_dir} ${weight_params_file}

	if [[ ! -f ${weight_params_file} ]]
	then
		echo "Could not generate mixing proxy model file @ ${weight_params_file}"
		exit 1
	fi

# Filter proxize typed variants.
haplocoded_resampled_proxy_tag_matbed=chr${CHR_ID}_Query_panel_haplocoded_resampled_proxy_tags.matbed.gz
allele_err_eps=0
N_THREADS=40
ProxyTyper_Release -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters ${haplocoded_resampled_tag_matbed} ${resampled_sample_list} ${weight_params_file} ${allele_err_eps} ${N_THREADS} ${haplocoded_proxy_resampled_tag_matbed}

	#####################################################################
	# Coordinate anonymization: This applies to both tag and target variants. Note that this protects the variant positions and the genetic map.
	anon_genetic_maps_dir=anon_coords_data
	beagle_genetic_map_file=${anon_genetic_maps_dir}/anon_beagle_map.map

	# Dump the original tag and target variants.
	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${haplocoded_resampled_proxy_tag_matbed} ${resampled_sample_list} 1 tag_vars.bed
	${PROXYTYPER_EXEC} -dump_plain_geno_signal_regions ${haplocoded_resampled_target_matbed} ${resampled_sample_list} 1 target_vars.bed

	#10      .       0.000000        72765
	rm -f -r ${anon_genetic_maps_dir}
	mkdir ${anon_genetic_maps_dir}
	ProxyTyper -anonymize_tag_target_genetic_map_coords tag_vars.bed target_vars.bed ${per_chrom_maps_dir}/${chr_id}.map ${coord_anon_cM_noise_SD} tag_target_mapping ${anon_genetic_maps_dir}

	# Writet the beagle formatted genetic map file.
	awk -v chr_id=${chr_id} {'if(NR>1){print chr_id"\t.\t"$3"\t"$1}'} ${anon_genetic_maps_dir}/${chr_id}.map > ${beagle_genetic_map_file} 

	# Map the coordinates for both sites.
	haplocoded_resampled_proxy_anon_coord_tag_matbed=chr${CHR_ID}_Query_panel_haplocoded_resampled_proxy_anon_coord_tags.matbed.gz
	haplocoded_resampled_proxy_anon_coord_target_matbed=chr${CHR_ID}_Query_panel_haplocoded_resampled_proxy_anon_coord_targets.matbed.gz

	ProxyTyper -generic_map_genotype_regs_per_src_dest_regs ${haplocoded_resampled_proxy_tag_matbed} ${resampled_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed ${haplocoded_resampled_proxy_anon_coord_tag_matbed}

	# Anonymize coordinates of target variants. Note that these are further permuted below.
	ProxyTyper -generic_map_genotype_regs_per_src_dest_regs ${haplocoded_resampled_target_matbed} ${resampled_sample_list} ${anon_genetic_maps_dir}/tag_target_mapping_original_coords.bed ${anon_genetic_maps_dir}/tag_target_mapping_mapping_coords.bed ${haplocoded_resampled_proxy_anon_coord_target_matbed}

	#####################################################################
	# Proxize the target variants by permutation: Note that this is done only at SITE2.
    #untyped_recoding_op_prefix="untyped_proxy"
	untyped_recoding_op_prefix="haplocoded_resampled_anon_coord"
	untyped_proxy_2_target_mapping_BED=${untyped_recoding_op_prefix}_target_proxy_mapping.bed
	haplocoded_resampled_anon_coord_permuted_target_matbed=${untyped_recoding_op_prefix}_proxized_targets.matbed.gz

	# Recode the targets for the site2: Use the anon-coord for this since we are operating on anonymized coordinates.
	ProxyTyper -recode_untyped_variant_reference_panel_per_target_permutation \
${haplocoded_resampled_proxy_anon_coord_tag_matbed} ${haplocoded_resampled_proxy_anon_coord_target_matbed} \
${resampled_sample_list} \
${untyped_var_perm_n_vicinity} \
${untyped_var_allele_switch_prob} \
${untyped_recoding_op_prefix}

	# Check the files.
	if [[ ! -f ${untyped_proxy_2_target_mapping_BED} ]]
	then
		echo "Could not find encoding coordinates @ \"${untyped_proxy_2_target_mapping_BED}\""
		exit 1
	fi

	if [[ ! -f ${haplocoded_resampled_anon_coord_permuted_target_matbed} ]]
	then
		echo "Could not find encoded target genotypes @ \"${haplocoded_resampled_anon_coord_permuted_target_matbed}\""
		exit 1
	fi

	##########################################################################################
	# Do partitioning of the untyped variants.
	# Shuffling is always hardcoded to increase privacy.
	shuffle_decomp_vars=1

	untyped_decomp_prefix="resampled_panel_target_anon_coord"

	# Decompose the current untyped variants, this does not change the typed variants, which are proxied above.
	ProxyTyper -simple_decompose_untyped_variants ${haplocoded_resampled_proxy_anon_coord_tag_matbed} 
${haplocoded_resampled_anon_coord_permuted_target_matbed} \
${resampled_sample_list} \
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



