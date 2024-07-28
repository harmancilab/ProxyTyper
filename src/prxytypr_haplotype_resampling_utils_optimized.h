#ifndef __OPTIM_RESAMPLING__
#define __OPTIM_RESAMPLING__

#include <vector>
//using namespace std;

void resample_phased_haplotypes_per_state_only_sampling(char* haplocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	int n_resampled_sample_size,
	double N_e_2_n_ref_haplotypes,
	double var_selection_segment_length_min_cM, // This is the minimum distance between the variants to switch states.
	double allele_error_prob,
	double segment_length_cutoff_bp,
	double segment_length_cutoff_cM,
	double segment_length_cutoff_in_var_number,
	vector<int>* all_haplotype_indices_2_sample,
	bool use_haplocoded_genotypes,
	int n_threads,
	int start_pos, int end_pos,
	int save_recombination_patterns,
	char* op_prefix);


#endif