#ifndef __PROXYGIMP_UTILS__
#define __PROXYGIMP_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;
class t_rng;

struct t_var_block
{
	int start_var_i;
	int end_var_i;
	vector<vector<int>*>* haplogroup_haplo_indices;
};


struct t_kmer
{
	int kmer_length;
	char* kmer;
};

bool compare_kmers(t_kmer* kmer1, t_kmer* kmer2);

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp);

void get_allele_error_per_decoded_panel_known_panel(char* decoded_panel_matbed_fp, char* known_panel_matbed_fp, char* sample_ids_list_fp, char* f_op);

void resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling(char* haplocoded_genotype_matrix_fp,
	char* haplocoded_target_genotype_matrix_fp,
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
	char* op_fp);

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
	char* op_fp);

void summarize_sampled_segments_per_resampled_haplotype_info(char* resampling_hap_info_sigbed_fp, char* resampling_sample_list_fp, char* generating_sample_list_fp, char* op_fp);
vector<t_annot_region*>* load_resampling_pattern_signal_var_regs(char* resampling_pattern_sig_BED_fp);
void save_resampling_pattern_signal_var_regs(vector<t_annot_region*>* haplocoded_geno_regs, int n_resampled_sample_size, const char* resampling_pattern_sig_BED_fp);

void random_per_sample_hap_switch(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp, double switch_prob, char* op_matbed_fp);

void linking_attack_per_haplotype_frequency_signatures(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vicinity,
	int n_step,
	char* op_fp);

void linking_attack_per_haplotype_frequency_signatures_per_ref_panel(char* target_query_panel_matbed_fp, char* target_query_sample_list_fp,
	char* ref_panel_matbed_fp, char* ref_panel_sample_list_fp,
	char* target_proxy_panel_matbed_fp, char* target_proxy_panel_sample_list_fp,
	int n_vicinity,
	int n_step,
	char* op_fp);

void combine_BEAGLE_imputed_decomposed_genotype_probabilities(char* beagle_imputed_vcf_fp,
	char* untyped_variant_decomposed_intervals_fp,
	char* panel_sample_list_fp,
	char* op_prefix);

vector<t_annot_region*>* load_decomposed_variant_Intervals(char* interval_fp);

void dump_decomposed_variant_Intervals(char* interval_fp, vector<t_annot_region*>* annot_regions);

void simple_decompose_untyped_variants(char* typed_var_geno_sig_regs_fp, char* untyped_var_geno_sig_regs_fp, char* panel_sample_list_fp, double min_AAF_per_decomp_var,
	bool shuffle_decomp_variants, char* op_prefix);

void get_untyped_variant_LD_statistics(char* typed_var_geno_sig_regs_fp, char* untyped_var_geno_sig_regs_fp, char* panel_sample_list_fp, int n_blocks_2_process, char* op_prefix);

void get_consecutive_block_variant_correlations(char* panel_target_geno_sig_regs_fp, char* panel_sample_list_fp, int l_block, int l_step, char* op_fp);

void decode_untyped_variant_reference_panel_per_target_permutation(char* recoded_ref_panel_target_haplocoded_matbed_fp,
	char* proxy_2_target_mapping_BED_fp,
	char* panel_sample_list_fp,
	char* op_matbed_fp);

void recode_untyped_variant_reference_panel_per_target_permutation(char* ref_panel_tag_haplocoded_matbed_fp,
	char* ref_panel_target_haplocoded_matbed_fp,
	char* panel_sample_list_fp,
	int untyped_var_perm_n_vicinity,
	double allele_switch_prob,
	char* op_prefix);

vector<int>* locally_permute_indices(int n_elements, int n_vicinity, t_rng* rng);

//int* matrix_get_independent_row_indices(double** matrix, int numRows, int numCols, int& numIndependentRows);
vector<int>* matrix_get_independent_row_indices(double** matrix, int numRows, int numCols);

bool matrix_is_identity(double** ident_check_mat, int nrow, int ncol, double ident_check_eps);

double** get_random_allele1_2_hap_matrix(t_rng* rng, int n_block_target_vars, int n_panel_haplotypes, double allele1_prob);

void calculate_proxy2clear_pairwise_distance_stats(char* proxized_panel_matbed, char* proxized_panel_sample_list_fp,
	char* cleartext_panel_matbed, char* cleartext_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	int l_cleartext_windowizing_win,
	int l_proxy_windowizing_win,
	int l_corr_win,
	char* op_prefix);

void generic_map_regions(vector<t_annot_region*>* target_regs, char* source_coords_bed_fp, char* dest_coords_bed_fp);

vector<t_annot_region*>* anonymize_genetic_distances_ret_regs(char* variants_BED_fp, char* genetic_map_file, double cM_noise_SD);

void anonymize_genetic_distances(char* variants_BED_fp, char* genetic_map_file, double cM_noise_SD, char* op_fp);

void anonymize_tag_target_genetic_map_coordinates(char* tag_variants_BED_fp, char* target_variants_BED_fp,
	char* genetic_map_file, double cM_noise_SD, char* mapper_op_prefix, char* op_dir);

bool sort_haplotypes(char* hap1, char* hap2);

void haplotype_frequency_attack(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int win_length,
	char* op_fp);

bool compare_haplotypes(char* hap1, char* hap2, int n_vars_per_win);

void generate_untyped_variant_recoded_reference_panel(char* ref_panel_tag_haplocoded_matbed_fp, char* ref_panel_target_haplocoded_matbed_fp, char* panel_sample_list_fp,
	char* op_dir, char* op_prefix);

void calculate_proxy2clear_var2var_correlation_stats(char* proxized_panel_matbed, char* proxized_panel_sample_list_fp,
	char* cleartext_panel_matbed, char* cleartext_panel_sample_list_fp);

void get_consecutive_variant_correlations(char* panel_target_geno_sig_regs_fp, char* panel_sample_list_fp);

void calculate_Sankararaman_LRT_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	double maf_cutoff,
	int min_var2var_dist,
	char* op_prefix);

void write_securegenome_input_files(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* op_prefix);

void calculate_windowed_Homer_t_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	int n_vicinity,
	char* op_prefix);

void calculate_Homer_t_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	char* op_prefix);

void decode_site_alleles_per_proxized_reference(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_proxized_haplocoded_geno_fp, char* ref_sample_ids_fp,
	int n_vicinity);

void decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	int var_start_i, int var_end_i,
	int n_vicinity,
	char* op_prefix);

void get_proxization1_2_proxization_mapped_haplotype_frequency_correlation(char* haplocoded_geno_signal_fp,
	char* proxy1_haplocoded_geno_signal_fp,
	char* proxy2_haplocoded_geno_signal_fp,
	char* sample_ids_list_fp,
	int n_vicinity,
	char* op_prefix);

vector<t_kmer*>* extract_kmers_per_haplotype(vector<char**>* per_ind_haplotypes, int i_win_start, int i_win_end);

void get_unique_kmers_w_counts(vector<t_kmer*>* par_kmers_list, vector<t_kmer*>* unique_kmers, vector<int>* unique_kmer_counts);

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	char* op_prefix);

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_MT(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	int n_threads,
	char* op_prefix);

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	char* op_prefix);

void decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	int n_vicinity);

void save_per_variant_site_mixing_parameters(vector<t_annot_region*>* geno_sig_regs, int n_vicinity, int avg_geno_mod, char* parameter_op_fp);

void proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
	char* proxized_haplocoded_var_regs_op_fp);

void MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
	int n_threads,
	char* proxized_haplocoded_var_regs_op_fp);

vector<t_annot_region*>* load_per_variant_site_mixing_parameters(char* parameters_fp,
	int& loaded_n_vicinity, int& loaded_avg_geno_mod);

void generate_save_per_site_mixing_parameters(vector<t_annot_region*>* geno_sig_regs,
	t_rng* rng,
	int n_vicinity,
	double per_var_weight_prob,
	double per_var2var_interaction_prob,
	int avg_geno_mod,
	char* parameter_op_fp);

void proxize_variants_per_vicinity_modular_average(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	int n_vicinity, int per_var_weight,
	int avg_geno_mod,
	double allele_err_eps,
	char* proxized_haplocoded_var_regs_op_fp);

void proxize_variants_per_vicinity_non_linear_modular_average_uniform(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	int n_vicinity,
	char* weight_params_fp,
	double per_var_weight_prob,
	double per_var2var_interaction_prob,
	int avg_geno_mod,
	double allele_err_eps,
	char* proxized_haplocoded_var_regs_op_fp);

void get_hapfreq_confusion_matrix(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vicinity,
	double min_HF_ref1_freq_2_quantify,
	double max_LF_ref2_freq_2_quantify,
	char* op_fp);

void get_per_win_haplotype_frequencies_per_reference(char* ref_panel_matbed_fp, char* ref_sample_list_fp, int n_vars_per_win, char* op_fp);

void convert_BEAGLE_imputed_meta_panel_2_matbed(char* vcf_fp, char* sample_ids_list_fp, double imp_geno_rand_weight, char* op_fp);

vector<char**>* get_per_subject_haplotypes_per_haplocoded_var_regs(vector<t_annot_region*>* haplocoded_panel_var_regs, vector<char*>* panel_sample_ids);

void assign_variant_AFs_to_var_regs(vector<t_annot_region*>* panel_var_regs, vector<char*>* sample_ids);

void generate_permute_proxizing_parameters(char* haplocoded_geno_sig_regs_fp, char* sample_ids_list_fp, int n_vicinity, double geno_inversion_prob, char* op_BED_fp);

void proxize_variants_per_locality_permutation(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* permute_proxy_mappings_BED_fp,
	char* proxized_haplocoded_var_regs_op_fp);

void pool_summarize_Homer_t_statistics_per_query(char* per_query_per_variant_files_list_fp, char* sample_ids_list_fp, char* summarized_results_op_fp);

void pool_summarize_Sankararaman_LRT_statistics_per_query(char* per_query_per_variant_files_list_fp, char* sample_ids_list_fp, char* summarized_results_op_fp);

void generate_save_per_site_mixing_parameters_LD_aware(vector<t_annot_region*>* all_geno_sig_regs,
	t_rng* rng,
	int n_vicinity,
	double per_var_weight_prob,
	double weight_inversion_prob,
	double per_var2var_interaction_prob,
	double per_var2var2var_interaction_prob,
	int avg_geno_mod,
	double normalized_N_e,
	int min_n_params_per_var,
	char* recombination_rate_dir,
	char* parameter_op_fp);

void compare_per_win_haplotypes(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vars_per_win,
	char* op_fp);

void get_query_haplotype_frequency_per_reference(char* query_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* op_fp);

#endif // __PROXYGIMP_UTILS__