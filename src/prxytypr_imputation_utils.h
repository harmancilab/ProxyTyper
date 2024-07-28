#ifndef __IMPUTATION_UTILS__
#define __IMPUTATION_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;

enum { MATCH_BY_NAME, MATCH_BY_START_POSN, N_MATCH_BY_IDS };
char** get_match_by_identifiers();


void random_phase_het_genotypes(char* genotype_regs_BED_fp, char* sample_ids_list_fp, char* haplocoded_genotype_matrix_op_fp);


vector<double*>* get_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes);

void count_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes);

vector<t_annot_region*>* load_recombination_rates(char* recombination_rate_fp);


double get_cumulative_recomb_rate_per_variant(t_annot_region* var_reg, vector<t_annot_region*>* recomb_regs);
double get_cumulative_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs);
double get_avg_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs);

void get_unique_haplotype_indices(vector<double*>* haplotype_alleles, vector<vector<int>*>* per_uniq_haplo_indices);

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option);

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_w_header_multithreaded(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	char* header_fp,
	bool save_phased_gt_option,
	int n_threads);

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_multithreaded(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option,
	int n_threads);

#endif // __IMPUTATION_UTILS__
