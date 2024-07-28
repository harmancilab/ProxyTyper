#ifndef __HAPLOTYPE_RESAMPLING_UTILS__
#define __HAPLOTYPE_RESAMPLING_UTILS__

#include <vector>
#include <algorithm>

using namespace std;

struct t_annot_region;

// Following write the sampling patterns into a parseable text file from all variants.
void summarize_sampled_segments_per_resampled_haplotype_info(char* resampling_hap_info_sigbed_fp, char* resampling_sample_list_fp, char* generating_sample_list_fp, char* op_fp);

vector<t_annot_region*>* load_resampling_pattern_signal_var_regs(char* resampling_pattern_sig_BED_fp, double& rel_resampling_N_e, int& n_sampling_sample_size, int& n_resampled_sample_size);
void save_resampling_pattern_signal_var_regs(vector<t_annot_region*>* haplocoded_geno_regs, int n_sampling_sample_size, int n_resampled_sample_size, double rel_resampling_N_e, const char* resampling_pattern_sig_BED_fp);

void resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb(char* haplocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	int n_resampled_sample_size,
	double N_e_2_n_ref_haplotypes,
	double allele_error_prob,
	double segment_length_cutoff_bp,
	double segment_length_cutoff_cM,
	double segment_length_cutoff_in_var_number,
	vector<int>* all_haplotype_indices_2_sample,
	int n_threads,
	int start_pos, int end_pos,
	char* op_prefix);

void resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling(char* haplocoded_genotype_matrix_fp,
	char* haplocoded_target_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	int n_resampled_sample_size,
	double N_e_2_n_ref_haplotypes,
	double min_cM_delta_per_recomb,
	double allele_error_prob,
	double segment_length_cutoff_bp,
	double segment_length_cutoff_cM,
	double segment_length_cutoff_in_var_number,
	vector<int>* all_haplotype_indices_2_sample,
	int n_threads,
	int start_pos, int end_pos,
	char* op_prefix);

double*** score_resampled_haplotypes_on_target_regions_per_ProxyTyper_sampling_patterns(vector<t_annot_region*>* target_regions,
	vector<t_annot_region*>* all_per_var_resampling_patterns,
	double rel_N_e,
	int n_resampled_sample_size,
	int n_original_haps);

void save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns(char* haplocoded_tag_genotype_matrix_fp, char* haplocoded_target_genotype_matrix_fp, char* sample_ids_list_fp,
	char* resampling_hap_info_sigbed_fp,
	char* recombination_rate_dir,
	int n_vic_tag_vars,
	int start_target_region,
	int end_target_region,
	char* op_dir);

void save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns_multithreaded(char* haplocoded_tag_genotype_matrix_fp, char* haplocoded_target_genotype_matrix_fp, char* sample_ids_list_fp,
	char* resampling_hap_info_sigbed_fp,
	char* recombination_rate_dir,
	int n_vic_tag_vars,
	int start_target_region,
	int end_target_region,
	int n_threads,
	char* op_dir);

#endif // __HAPLOTYPE_RESAMPLING_UTILS__