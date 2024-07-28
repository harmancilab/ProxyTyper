#ifndef __UNPHASED_PANEL_RESAMPLING_UTILS__
#define __UNPHASED_PANEL_RESAMPLING_UTILS__

int get_kmer_distance(char* kmer1, char* kmer2, int l_win);

void sample_unphased_genotype_phasings(char* unphased_geno_matrix_matbed, char* unphased_geno_matrix_sample_list,
	char* ref_haplocoded_geno_matrix_matbed, char* ref_haplocoded_geno_matrix_sample_list,
	int n_resample_per_subject,
	int l_win, double dist_weighter,
	char* op_fp);

void resample_unphased_genotypes_per_state_only_sampling(char* genocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	int n_resampled_sample_size,
	double N_e_2_n_ref_haplotypes,
	double var_selection_segment_length_min_cM, // This is the minimum distance between the variants to switch states.
	double allele_error_prob,
	double segment_length_cutoff_bp,
	double segment_length_cutoff_cM,
	double segment_length_cutoff_in_var_number,
	int n_threads,
	int start_pos, int end_pos,
	int save_recombination_patterns,
	char* op_prefix);

#endif // __UNPHASED_PANEL_RESAMPLING_UTILS__