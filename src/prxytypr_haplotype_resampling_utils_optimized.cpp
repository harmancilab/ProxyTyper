#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ctime>
#include "prxytypr_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_variation_tools.h"
#include "prxytypr_imputation_utils.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_nomenclature.h"
#include "prxytypr_histogram.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_vector_macros.h"
#include "prxytypr_haplotype_resampling_utils.h"
#include "prxytypr_proxytyper.h"


using namespace std;

const bool __DUMP_OPTIM_HAPRESAMPLING_MSGS__ = false;

static void* thread_callback_save_resampled_genotypes_per_recombination_patterns(void* __thread_ptr)
{
	void** thread_info_ptr = (void**)(__thread_ptr);
	int* int_vals = (int*)(thread_info_ptr[0]);
	int thread_i = int_vals[0];
	//int n_threads = int_vals[1];
	int n_resampled_sample_size = int_vals[2];
	int var_i_start = int_vals[3];
	int var_i_end = int_vals[4];
	int use_haplocoded_genotypes_val = int_vals[5];
	int seconds_since_epoch = int_vals[6];

	// This is the sorted list of regions.
	vector<t_annot_region*>* haplocoded_geno_regs = (vector<t_annot_region*>*)(thread_info_ptr[1]);
	//vector<char*>* sample_ids = (vector<char*>*)(thread_info_ptr[2]);
	char* op_prefix = (char*)(thread_info_ptr[3]);

	//// Divide the regions into blocks.
	//int var_i_start = thread_i * vecsize(haplocoded_geno_regs) / n_threads;
	//int var_i_end = MIN(var_i_start + vecsize(haplocoded_geno_regs) / n_threads, vecsize(haplocoded_geno_regs));

	char cur_thread_resampled_subject_ids_fp[1000];
	sprintf(cur_thread_resampled_subject_ids_fp, "%s_%d_subjects.list", op_prefix, thread_i);
	FILE* f_cur_thread_resampled_subject_ids = open_f(cur_thread_resampled_subject_ids_fp, "w");
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		fprintf(f_cur_thread_resampled_subject_ids, "RES_%d_%d\n", seconds_since_epoch, i_s);
	} // i_s loop.
	close_f(f_cur_thread_resampled_subject_ids, NULL);

	char cur_thread_resampled_vars_BED_fp[1000];
	sprintf(cur_thread_resampled_vars_BED_fp, "%s_%d_variants.bed", op_prefix, thread_i);
	FILE* f_cur_thread_resampled_vars_BED = open_f(cur_thread_resampled_vars_BED_fp, "w");

	char cur_thread_resampled_geno_mat_fp[1000];
	sprintf(cur_thread_resampled_geno_mat_fp, "%s_%d_genotypes.matrix.gz", op_prefix, thread_i);
	FILE* f_geno_matrix = open_f(cur_thread_resampled_geno_mat_fp, "wb");

	// Note that these arrays hold the state over all variants for all subjects while we are tracing the variants.
	char* cur_var_per_subj_resampled_geno = new char[n_resampled_sample_size + 1];
	int* cur_var_per_subj_resampled_hap0_states = new int[n_resampled_sample_size + 2];
	int* cur_var_per_subj_resampled_hap1_states = new int[n_resampled_sample_size + 2];
	memset(cur_var_per_subj_resampled_hap0_states, 0xff, sizeof(int) * (n_resampled_sample_size));
	memset(cur_var_per_subj_resampled_hap1_states, 0xff, sizeof(int) * (n_resampled_sample_size));

	// Sanity check to make sure we are at a recombable variant:
	if ((signed int)haplocoded_geno_regs->at(var_i_start)->score == -1)
	{
		fprintf(stderr, "Thread-%d block is not at a recombable variant @ %s(%d)\n", thread_i, __FILE__, __LINE__);
		exit(1);
	}

	for (int var_i = var_i_start; var_i < var_i_end; var_i++)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////
		// First, write the signal region to the bed file.
		if (haplocoded_geno_regs->at(var_i)->strand == '-' ||
			haplocoded_geno_regs->at(var_i)->strand == '+')
		{
			if (haplocoded_geno_regs->at(var_i)->name == NULL)
			{
				haplocoded_geno_regs->at(var_i)->name = new char[5];
				strcpy(haplocoded_geno_regs->at(var_i)->name, ".");
			}

			// Translate the start and end to CODEBASE's start and end.
			fprintf(f_cur_thread_resampled_vars_BED, "%s\t%d\t%d\t%s\t.\t%c\n", haplocoded_geno_regs->at(var_i)->chrom,
				translate_coord(haplocoded_geno_regs->at(var_i)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(haplocoded_geno_regs->at(var_i)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				haplocoded_geno_regs->at(var_i)->name,
				haplocoded_geno_regs->at(var_i)->strand);
		}
		else
		{
			fprintf(f_cur_thread_resampled_vars_BED, "%s\t%d\t%d\n", haplocoded_geno_regs->at(var_i)->chrom,
				translate_coord(haplocoded_geno_regs->at(var_i)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(haplocoded_geno_regs->at(var_i)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		// Start writing the genotype signal.
		void** cur_reg_info = (void**)(haplocoded_geno_regs->at(var_i)->data);
		char* var_geno_sig = (char*)(cur_reg_info[0]);
		int** per_hap_resampled_hap_i = NULL;
		if ((signed int)haplocoded_geno_regs->at(var_i)->score != -1)
		{
			per_hap_resampled_hap_i = (int**)(cur_reg_info[2]);
		}

		// Reset the genotype.
		memset(cur_var_per_subj_resampled_geno, 0, sizeof(char) * n_resampled_sample_size);

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{		
			// Update the current states using the states.
			if (per_hap_resampled_hap_i != NULL && per_hap_resampled_hap_i[0][i_s] != -1)
			{
				cur_var_per_subj_resampled_hap0_states[i_s] = per_hap_resampled_hap_i[0][i_s];
			}

			if (per_hap_resampled_hap_i != NULL && per_hap_resampled_hap_i[1][i_s] != -1)
			{
				cur_var_per_subj_resampled_hap1_states[i_s] = per_hap_resampled_hap_i[1][i_s];
			}

			// This is the critical check to make sure we have an assigned state.
			// The basic idea is that first variant initializes the states so we should always have a state for every subject at any variant.
			if (cur_var_per_subj_resampled_hap0_states[i_s] == -1 ||
				cur_var_per_subj_resampled_hap1_states[i_s] == -1)
			{
				fprintf(stderr, "Sanity check failed: Found non-tracked site @ subject-%d; var_i=%d\n", i_s, var_i);
				exit(1);
			}

			int hap0_ref_hap_i = cur_var_per_subj_resampled_hap0_states[i_s];
			int hap0_ref_i_s = hap0_ref_hap_i >> 1;
			int hap0_ref_all_i = hap0_ref_hap_i % 2;

			int hap1_ref_hap_i = cur_var_per_subj_resampled_hap1_states[i_s];
			int hap1_ref_i_s = hap1_ref_hap_i >> 1;
			int hap1_ref_all_i = hap1_ref_hap_i % 2;

			if (use_haplocoded_genotypes_val > 0)
			{
				cur_var_per_subj_resampled_geno[i_s] = get_allele_per_haplotype(var_geno_sig[hap0_ref_i_s], hap0_ref_all_i)+
														get_allele_per_haplotype(var_geno_sig[hap1_ref_i_s], hap1_ref_all_i) * 2;
			}
			else
			{
				cur_var_per_subj_resampled_geno[i_s] = get_allele_per_haplotype(var_geno_sig[hap0_ref_i_s], hap0_ref_all_i) +
														get_allele_per_haplotype(var_geno_sig[hap1_ref_i_s], hap1_ref_all_i);
			}
		} // i_s loop.

		// Write the current variant.
		fwrite(cur_var_per_subj_resampled_geno, sizeof(char), n_resampled_sample_size, f_geno_matrix);
	} // var_i loop.
	close_f(f_geno_matrix, cur_thread_resampled_geno_mat_fp);
	close_f(f_cur_thread_resampled_vars_BED, NULL);
	
	return NULL;
} // thread_callback_save_resampled_genotypes_per_recombination_patterns option.

// This is the final step in saving the genotypes using the recombination patterns.
void save_resampled_genotypes_per_recombination_patterns_multithreaded(vector<t_annot_region*>* haplocoded_geno_regs, vector<t_annot_region*>* recombable_geno_regs, vector<char*>* sample_ids,
	int n_resampled_sample_size,
	bool use_haplocoded_genotypes,
	int n_threads,
	char* op_prefix)
{
	auto genotype_saving_start_chrono = std::chrono::high_resolution_clock::now();
	vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();

	std::time_t seconds_since_epoch = std::time(nullptr);

	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		void** thread_info = new void* [20];

		// Divide the regions into blocks of recombable variants so that we know we are aligned on recombination spots.
		int n_recombable_vars_per_thread = 1 + (vecsize(recombable_geno_regs) / n_threads);
		int recombable_var_i_start = MIN(thread_i * n_recombable_vars_per_thread, vecsize(recombable_geno_regs) - 1);
		int recombable_var_i_end = MIN(recombable_var_i_start + n_recombable_vars_per_thread, vecsize(recombable_geno_regs)-1);

		if (recombable_var_i_start == recombable_var_i_end)
		{
			fprintf(stderr, "Early stopping @ thread %d\n", thread_i);
			break;
		}

		// Translate the recombable variant indices.
		int var_i_start = recombable_geno_regs->at(recombable_var_i_start)->score;
		int var_i_end = recombable_geno_regs->at(recombable_var_i_end)->score;

		if (recombable_geno_regs->at(recombable_var_i_end)->start != haplocoded_geno_regs->at(var_i_end)->start ||
			recombable_geno_regs->at(recombable_var_i_start)->start != haplocoded_geno_regs->at(var_i_start)->start)
		{
			fprintf(stderr, "Sanity check failed for recombable regions @ %s(%d)..\n", __FILE__, __LINE__);
			exit(1);
		}

		// First block always starts from first variant.
		if (recombable_var_i_start == 0)
		{
			var_i_start = 0;
		}

		// If we are at the end of recombable blocks, set this to the last block.
		if (recombable_var_i_end == (vecsize(recombable_geno_regs) - 1))
		{
			var_i_end = vecsize(haplocoded_geno_regs);
		}

		if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
		{ 
			t_string::print_padded_string(stderr, '\r', 100, "Thread-%d: Variants @ All-vars:[%d-%d]; Recomb-vars:[%d-%d]/%d", thread_i,
				var_i_start, var_i_end,
				recombable_var_i_start, recombable_var_i_end,
				vecsize(recombable_geno_regs));
		}

		if (var_i_start > var_i_end)
		{
			fprintf(stderr, "Sanity check failed @ %s(%d)\n", __FILE__, __LINE__);
			exit(1);
		}
		
		int* int_vals = new int[20];
		int_vals[0] = thread_i;
		int_vals[1] = n_threads;
		int_vals[2] = n_resampled_sample_size;
		int_vals[3] = var_i_start;
		int_vals[4] = var_i_end;
		int_vals[5] = (int)use_haplocoded_genotypes;
		int_vals[6] = (int)seconds_since_epoch;
		thread_info[0] = int_vals;

		thread_info[1] = haplocoded_geno_regs;
		thread_info[2] = sample_ids;
		thread_info[3] = op_prefix;

		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_save_resampled_genotypes_per_recombination_patterns, thread_info);
		threads->push_back(cur_thread);
		cur_thread->run_thread();
		t_string::print_padded_string(stderr, '\r', 100, "Started Thread-%d..", thread_i);
	} // thread_i loop.

	t_string::print_padded_string(stderr, '\n', 100, "Waiting for threads..");

	for (int thread_i = 0; thread_i < vecsize(threads); thread_i++)
	{
		threads->at(thread_i)->wait_thread();
		t_string::print_padded_string(stderr, '\r', 100, "Thread %d finished..", thread_i);
	} // thread_i loop.

	t_string::print_padded_string(stderr, '\n', 100, "Pooling, saving, and cleaning up..");

	// Concatenating results:
	vector<t_annot_region*>* pooled_var_regs = new vector<t_annot_region*>();
	vector<char*>* var_BED_files = new vector<char*>();
	vector<char*>* matrix_gzip_files = new vector<char*>();
	vector<char*>* sample_ids_list_files = new vector<char*>();
	for (int thread_i = 0; thread_i < vecsize(threads); thread_i++)
	{
		char cur_var_regs_BED[1000];
		sprintf(cur_var_regs_BED, "%s_%d_variants.bed", op_prefix, thread_i);
		vector<t_annot_region*>* cur_vars = load_BED(cur_var_regs_BED);
		pooled_var_regs->insert(pooled_var_regs->end(), cur_vars->begin(), cur_vars->end());
		var_BED_files->push_back(t_string::copy_me_str(cur_var_regs_BED));

		char cur_matrix_gz[1000];
		sprintf(cur_matrix_gz, "%s_%d_genotypes.matrix.gz", op_prefix, thread_i);
		matrix_gzip_files->push_back(t_string::copy_me_str(cur_matrix_gz));

		char cur_sample_list_fp[1000];
		sprintf(cur_sample_list_fp, "%s_%d_subjects.list", op_prefix, thread_i);
		sample_ids_list_files->push_back(t_string::copy_me_str(cur_sample_list_fp));
	} // thread_i loop.

	// Save the genotypes matrix (M).
	char pooled_geno_matrix_fp[1000];
	sprintf(pooled_geno_matrix_fp, "%s_genotypes.matrix.gz", op_prefix);
	concatenateGzipFiles(pooled_geno_matrix_fp, matrix_gzip_files);

	// Use one of the files as the pooled sample identifiers (S).
	vector<char*>* pooled_sample_ids = buffer_file(sample_ids_list_files->at(0));
	char subject_ids_list_fp[1000];
	sprintf(subject_ids_list_fp, "%s_subjects.list", op_prefix);
	save_lines(pooled_sample_ids, subject_ids_list_fp);

	// Save the regions (R).
	char pooled_var_BED_fp[1000];
	sprintf(pooled_var_BED_fp, "%s_variants.bed", op_prefix);
	dump_BED(pooled_var_BED_fp, pooled_var_regs);

	// Clean up the files.
	fprintf(stderr, "Cleaning up...                \n");
	for (int i_f = 0; i_f < vecsize(matrix_gzip_files); i_f++)
	{
		delete_file(matrix_gzip_files->at(i_f));
		delete_file(sample_ids_list_files->at(i_f));
		delete_file(var_BED_files->at(i_f));
	} // i_f loop.

	vector<t_annot_region*>* validation_regs = load_BED(pooled_var_BED_fp);
	if (vecsize(validation_regs) != vecsize(haplocoded_geno_regs))
	{
		fprintf(stderr, "Sanity check failed while validating the saved regions: %d/%d @ %s(%d)\n",
				vecsize(validation_regs), vecsize(haplocoded_geno_regs), __FILE__, __LINE__);

		exit(1);
	}

	auto genotype_saving_end_chrono = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> genotype_saving_duration = genotype_saving_end_chrono - genotype_saving_start_chrono;
	fprintf(stderr, "Genotype saving finished in %.3f seconds..\n", genotype_saving_duration.count());
}

// THIS IS THE FUNCTION USED FOR RESAMPLING TAGS WHILE STORING RECOMB INFORMATION.
static void* thread_callback_haplotype_state_only_resampling(void* thread_info_ptr)
{
	void** thread_ptrs_list = (void**)(thread_info_ptr);

	int* int_vals = (int*)(thread_ptrs_list[0]);
	//int_vals[0] = thread_i;
	int thread_i = int_vals[0];
//int_vals[1] = n_threads;
	int n_threads = int_vals[1];
//int_vals[2] = start_pos;
	int start_coord = int_vals[2];
//int_vals[3] = end_pos;
	int end_coord = int_vals[3];
//int_vals[4] = n_original_haps;
	int n_original_haps = int_vals[4];
//int_vals[5] = n_resampled_sample_size;
	int n_resampled_sample_size = int_vals[5];

	double* dbl_vals = (double*)(thread_ptrs_list[1]);
	//dbl_vals[2] = N_e;
	double N_e = dbl_vals[0];
	//dbl_vals[3] = allele_error_prob;
	double allele_error = dbl_vals[1];
	//dbl_vals[4] = segment_length_cutoff_bp;
	double length_cutoff_in_bps = dbl_vals[2];
	//dbl_vals[5] = segment_length_cutoff_cM;
	double length_cutoff_in_cM = dbl_vals[3];
	//dbl_vals[6] = segment_length_cutoff_in_var_number;
	double length_cutoff_in_var_number = dbl_vals[4];

	t_restr_annot_region_list* restr_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[2]);

	//char* recombination_rate_dir = (char*)(thread_ptrs_list[3]);

	vector<int>* haplotype_indices_2_sample = (vector<int>*)(thread_ptrs_list[4]);

	if ((int)haplotype_indices_2_sample->size() != n_original_haps)
	{
		fprintf(stderr, "Could not match the number of haplotypes: %d, %d\n", (int)haplotype_indices_2_sample->size(), (int)n_original_haps);
		exit(1);
	}

	if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
	{
		fprintf(stderr, "Length cutoffs:\n\
Max segment in bps: %.1f\n\
Max segment in cMs: %.1f\n\
Max segment in n vars.: %.1f\n", length_cutoff_in_bps, length_cutoff_in_cM, length_cutoff_in_var_number);
	}

	// Load the recombination rates.
	int cur_thread_seed = time(NULL) + thread_i;
	t_rng* rng = new t_rng(cur_thread_seed);
	//t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
	{
		fprintf(stderr, "Setting thread %d/%d: N_e: %d, error_prob: %.5f, n_ref_haps: %d, n_resampled_size: %d; Interval: [%d-%d]; random_seed: %d\n",
			thread_i,
			n_threads,
			(int)N_e,
			allele_error,
			n_original_haps,
			n_resampled_sample_size,
			start_coord, end_coord,
			cur_thread_seed);
	}

	if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
	{
		fprintf(stderr, "First 5 random numbers: ");
		for (int r_i = 0; r_i < 5; r_i++)
		{
			fprintf(stderr, "%.4f, ", rng->random_double_ran3());
		} // r_i loop.
		fprintf(stderr, "\n");
	}

	for (int i_chr = 0; i_chr < (int)restr_geno_regs->chr_ids->size(); i_chr++)
	{
		//fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];

		//vector<char**>* cur_chr_per_subj_haplotypes = per_chr_per_subj_haplotypes[i_chr];

		double* per_var_self_probs = new double[vecsize(cur_chrom_var_regs) + 2];
		double* per_var_other_probs = new double[vecsize(cur_chrom_var_regs) + 2];

		per_var_self_probs[0] = 1.0;
		per_var_other_probs[0] = 0.0;

		// Get the recombination rate between these positions.
		for (int i_reg = 1; i_reg < vecsize(cur_chrom_var_regs); i_reg++)
		{
			double r_m = fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_chrom_var_regs->at(i_reg - 1)->dbl_score);
			double rho_m = 4 * N_e * r_m;

			double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

			double other_prob = tau_m / n_original_haps;
			double self_prob = (1 - tau_m) + (tau_m / n_original_haps);

			per_var_self_probs[i_reg] = self_prob;
			per_var_other_probs[i_reg] = other_prob;
		} // i_reg loop.

		// Start resampling.
		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			if (i_s % n_threads != thread_i)
			{
				continue;
			}

			//if (thread_i == 0)
			//{
			//	fprintf(stderr, "Sampling Subject %d..\n", i_s);
			//}

			// Sample the two haplotypes.
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				int TRACKED_HAPLO_STATE = -1;
				int n_recombs = 0;
				int n_erroneous_alleles = 0;
				int cur_segment_start_var_i = -1;
				int cur_segment_start_bps = -1;
				double cur_segment_start_cM = -1;

				// Following loop tracks all variants in order from left to right while updating haplo state.
				for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
				{
					if (start_coord <= 0 || end_coord <= 0 ||
						(cur_chrom_var_regs->at(i_reg)->start > start_coord &&
							cur_chrom_var_regs->at(i_reg)->end < end_coord))
					{

					}
					else
					{
						continue;
					}

					void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
					int* cur_reg_n_recombs = (int*)(cur_reg_info[1]);
					int** cur_reg_per_hap_resampled_hap_i = (int**)(cur_reg_info[2]);

					// This check may impact performance but it is used to make sure initializations are correctly done.
					if (cur_reg_per_hap_resampled_hap_i[i_hap][i_s] != -1)
					{
						fprintf(stderr, "Sanity check failed: Hap state Initialization error!!!\n");
						exit(1);
					}

					if (TRACKED_HAPLO_STATE == -1)
					{
						// Following selects the next haplotype as the initial haplotype, this does a uniform sampling of the initial states.
						int cur_haplo_state_i = MIN((n_original_haps - 1), floor(rng->random_double_ran3() * n_original_haps));
						//int cur_haplo_state_i = (i_s * 2 + i_hap) % n_original_haps;
						TRACKED_HAPLO_STATE = haplotype_indices_2_sample->at(cur_haplo_state_i);

						// Update the new segment's coordinates.
						cur_segment_start_bps = cur_chrom_var_regs->at(i_reg)->start;
						cur_segment_start_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
						cur_segment_start_var_i = i_reg;
					}
					else if ((signed int)cur_chrom_var_regs->at(i_reg)->score != -1) // Make sure that tihs is a recombable variant. If not, we just continue on the current tracked state.
					{
						double other_prob = per_var_other_probs[i_reg];
						double self_prob = per_var_self_probs[i_reg];

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// This calculates the self and other probs, use for sanity checks.
						//// Get the recombination rate between these positions.
						//double r_m = fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_chrom_var_regs->at(i_reg - 1)->dbl_score);
						//double rho_m = 4 * N_e * r_m;

						//double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

						//double other_prob_check = tau_m / n_original_haps;
						//double self_prob_check = (1 - tau_m) + (tau_m / n_original_haps);
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// Keep the distance constraints in.
						// Check the segment length so far and assign uniform if this segment is too long.
						bool reset_trans_probs = false;
						if (length_cutoff_in_bps > 0 &&
							(cur_chrom_var_regs->at(i_reg)->start - cur_segment_start_bps) > length_cutoff_in_bps)
						{
							reset_trans_probs = true;
						}

						if (length_cutoff_in_cM > 0 &&
							(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_segment_start_cM) > length_cutoff_in_cM)
						{
							reset_trans_probs = true;
						}

						if (length_cutoff_in_var_number > 0 &&
							(i_reg - cur_segment_start_var_i) > length_cutoff_in_var_number)
						{
							reset_trans_probs = true;
						}

						if (reset_trans_probs)
						{
							other_prob = (double)1.0 / n_original_haps;
							self_prob = (double)1.0 / n_original_haps;
						}
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#define __TOTAL_LOOP__
#define __SELF_OTHER_SAMPLER__

#ifdef __SELF_OTHER_SAMPLER__
						// This performs haplotype sampling faster than full loop.
						double total_other_prob = (n_original_haps - 1) * other_prob;
						double total_self_prob = self_prob;

						if (fabs(total_other_prob + total_self_prob - 1.0) > 0.001)
						{
							fprintf(stderr, "Sanity check failed: Total probability is not as expected: %.5f vs %.5f\n", total_other_prob, total_self_prob);
							exit(1);
						}

						// Use this random value for sampling next haplotype state.
						double rand_cumul_val_self_other = rng->random_double_ran3();

						if (rand_cumul_val_self_other < total_self_prob)
						{
							TRACKED_HAPLO_STATE = TRACKED_HAPLO_STATE;
						} // self sampling check.
						else
						{
							// Select a non-self haplotype by rejection sampling.
							int cur_sampled_haplo_i = MIN(n_original_haps - 1, (int)(rng->random_double_ran3() * n_original_haps));
							while (cur_sampled_haplo_i == TRACKED_HAPLO_STATE)
							{
								cur_sampled_haplo_i = MIN(n_original_haps - 1, (int)(rng->random_double_ran3() * n_original_haps));
							} // resampler for non-self haplotype.

							// Update the recombinations.
							cur_reg_n_recombs[i_s]++;
							n_recombs++;

							// Update the new segment's coordinates.
							cur_segment_start_bps = cur_chrom_var_regs->at(i_reg)->start;
							cur_segment_start_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
							cur_segment_start_var_i = i_reg;

							TRACKED_HAPLO_STATE = cur_sampled_haplo_i;
						} // non-self sampling check.
#endif // __SELF_OTHER_SAMPLER__

#ifdef __TOTAL_LOOP__
						// Use this random value for sampling next haplotype state.
						double rand_cumul_val = rng->random_double_ran3();

						// Following is the full loop that samples the next haplotype.
						double cur_cumul_val = 0.0;
						for (int i_hap_state_i = 0; i_hap_state_i < n_original_haps; i_hap_state_i++)
						{
							int hap_state_i = haplotype_indices_2_sample->at(i_hap_state_i);

							// Update the cumulative.
							if (hap_state_i == TRACKED_HAPLO_STATE)
							{
								cur_cumul_val += self_prob;
							}
							else
							{
								cur_cumul_val += other_prob;
							}

							// Check the cumulative.
							if (cur_cumul_val > rand_cumul_val)
							{
								if (TRACKED_HAPLO_STATE != hap_state_i)
								{
									// Update the recombinations.
									cur_reg_n_recombs[i_s]++;
									n_recombs++;

									// Update the new segment's coordinates.
									cur_segment_start_bps = cur_chrom_var_regs->at(i_reg)->start;
									cur_segment_start_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
									cur_segment_start_var_i = i_reg;
								}

								TRACKED_HAPLO_STATE = hap_state_i;

								if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
								{
									fprintf(stderr, "SAMPLED %d[%d] var %d: state: %d\n", i_s, i_hap,
										cur_chrom_var_regs->at(i_reg)->start,
										TRACKED_HAPLO_STATE);
								}

								break;
							}
						} // i_hap_state_i loop.
#endif // __TOTAL_LOOP__
					} // region check for initing the haplo state.
					else
					{
						fprintf(stderr, "Sanity check failed: We are not supposed to be @ %s(%d)\n", __FILE__, __LINE__);
						exit(1);
					}

					// Save the resampled haplotype index, make sure this is stored here so that we do not miss it out.
					cur_reg_per_hap_resampled_hap_i[i_hap][i_s] = TRACKED_HAPLO_STATE;
				} // i_reg loop.

				//int thousand_fact = 1000 / n_threads;
				//if ((i_s - thread_i) % (n_threads * thousand_fact) == 0)
				if (__DUMP_OPTIM_HAPRESAMPLING_MSGS__)
				{
					if (i_s % 100 == 0 && thread_i == 0)
					{
						t_string::print_padded_string(stderr, '\r', 100, "Thread %d: Re-sampled sample %d (%d): %d recombinations (%d allelic errors).", thread_i, i_s, i_hap, n_recombs, n_erroneous_alleles);
					}
				}
			} // i_hap loop.
		} // i_s loop.
	} // i_chr loop.

	return NULL;
} // thread_callback_haplotype_state_only_resampling function.

static void* thread_callback_hapstat_transition_memalloc(void* __thread_params)
{
	void** thread_params = (void**)(__thread_params);
	int* int_vals = (int*)(thread_params[0]);
	int thread_i = int_vals[0];
	int n_threads = int_vals[1];
	int n_resampled_sample_size = int_vals[2];

	vector<t_annot_region*>* cur_chrom_recomb_var_regs = (vector<t_annot_region*>*)(thread_params[1]);

	int* n_recombs_mem_pool = (int*)(thread_params[2]);
	int* per_hap_sampled_hap_i_mem_pool = (int*)(thread_params[3]);

	size_t n_regs = vecsize_t(cur_chrom_recomb_var_regs);

	for (size_t i_reg = 0; i_reg < n_regs; i_reg++)
	{
		if (i_reg % n_threads != (size_t)thread_i)
		{
			continue;
		}

		void** cur_reg_info = (void**)(cur_chrom_recomb_var_regs->at(i_reg)->data);
		void** new_reg_info = new void* [10];
		new_reg_info[0] = cur_reg_info[0];

		int* n_recombs_per_var = n_recombs_mem_pool + i_reg * (size_t)(n_resampled_sample_size);
		new_reg_info[1] = n_recombs_per_var;

		int** per_hap_sampled_hap_i = new int* [2];
		int* base_per_hap_sampled_hap0_i = per_hap_sampled_hap_i_mem_pool + 2 * i_reg * (size_t)(n_resampled_sample_size);
		int* base_per_hap_sampled_hap1_i = per_hap_sampled_hap_i_mem_pool + (2 * i_reg + 1) * (size_t)(n_resampled_sample_size);

		per_hap_sampled_hap_i[0] = base_per_hap_sampled_hap0_i;
		per_hap_sampled_hap_i[1] = base_per_hap_sampled_hap1_i;

		new_reg_info[2] = per_hap_sampled_hap_i;
		// Replace the info ptr.
		cur_chrom_recomb_var_regs->at(i_reg)->data = new_reg_info;
	} // i_reg loop.

	return(NULL);
}

// THIS FUNCTION RESAMPLES TAGS AND SAVES THEM: We need to modify this to samples recombination 
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
	char* op_prefix)
{
	auto overall_resampling_start_chrono = std::chrono::high_resolution_clock::now();

	fprintf(stderr, "%d-thread Haplotype-state centric Re-sampling:\n\
genotype matrix: %s\n\
n_resampled_size: %d\n\
N_e/N_ref_hap: %.3f\n\
allelic error: %.6f\n\
var_selection_segment_length_min_cM: %.4f\n\
segment_length_in_bp: %.0f\n\
segment_length_cutoff_cM: %.3f\n\
segment_length_cutoff_in_var_number: %.0f\n\
genomic interval: [%d-%d]\n\
use_haplocoded_genotypes: %d\n\
save_recombination_patterns: %d\n",
		n_threads,
		haplocoded_genotype_matrix_fp,
		n_resampled_sample_size,
		N_e_2_n_ref_haplotypes,
		allele_error_prob,
		var_selection_segment_length_min_cM,
		segment_length_cutoff_bp, segment_length_cutoff_cM, segment_length_cutoff_in_var_number,
		start_pos, end_pos,
		(int)use_haplocoded_genotypes, save_recombination_patterns);

	vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_genotype_matrix_fp, sample_ids_list_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n", (int)haplocoded_geno_regs->size(), (int)sample_ids->size());

	vector<int>* haplotype_indices_2_sample = NULL;
	int n_original_haps = 0;
	if (all_haplotype_indices_2_sample != NULL &&
		all_haplotype_indices_2_sample->size() > 0)
	{
		haplotype_indices_2_sample = all_haplotype_indices_2_sample;
		n_original_haps = (int)haplotype_indices_2_sample->size();
	}
	else
	{
		haplotype_indices_2_sample = new vector<int>();
		for (int i_hap = 0; i_hap < 2 * (int)sample_ids->size(); i_hap++)
		{
			haplotype_indices_2_sample->push_back(i_hap);
		} // i_hap loop.
		n_original_haps = 2 * (int)sample_ids->size();
	}

	int max_genotype = get_max_genotype_value(haplocoded_geno_regs, sample_ids);
	if (max_genotype != 3)
	{
		fprintf(stderr, "***Could not detect phased genotypes, ProxyTyper will continue but this is highly unusual. Make sure that the panel is phased.***\n");
		fprintf(stderr, "***Could not detect phased genotypes, ProxyTyper will continue but this is highly unusual. Make sure that the panel is phased.***\n");
		fprintf(stderr, "***Could not detect phased genotypes, ProxyTyper will continue but this is highly unusual. Make sure that the panel is phased.***\n");
	}


	double N_e = N_e_2_n_ref_haplotypes * n_original_haps;

	// This is used to speed-up resampling by setting a minimum distance between recombination "hotspots".
	double min_cM_delta_per_recomb = var_selection_segment_length_min_cM;

	fprintf(stderr, "Haplotype-state Sampling from %d haplotypes:\n\
N_e=%.3f\n\
allele_eps=%.3f\n\
min_cM_per_recombable var=%.4f.\n", (int)(n_original_haps), N_e, allele_error_prob, min_cM_delta_per_recomb);

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	//vector<char**>** per_chr_per_subj_haplotypes = new vector<char**>*[restr_geno_regs->chr_ids->size()];
	vector<t_annot_region*>** per_chrom_recombable_regs = new vector<t_annot_region*>*[restr_geno_regs->chr_ids->size()];;
	for (int i_chr = 0; i_chr < vecsize(restr_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates.\n");
		for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Assign variant clusterings based on cM distance.
		double cur_cM = -100; // Initialize to a very far position so we make sure to add 1st variant as recombable.
		int n_recombable_vars = 0;
		vector<t_annot_region*>* cur_chrom_recombable_vars = new vector<t_annot_region*>();
		for (int i_reg = 0; i_reg < vecsize(cur_chrom_var_regs); i_reg++)
		{
			// IF the current variant is further away than the minimum cM cutoff, set it as a recomb. hotspot.
			if (fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_cM) > min_cM_delta_per_recomb)
			{
				cur_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
				cur_chrom_var_regs->at(i_reg)->score = i_reg; // This enables us to know the region's index on the original sorted array.

				cur_chrom_recombable_vars->push_back(cur_chrom_var_regs->at(i_reg));
				n_recombable_vars++;
			}
			else
			{
				cur_chrom_var_regs->at(i_reg)->score = -1;
			}
		} // i_reg loop.

		fprintf(stderr, "%d (%d) recombable variants..\n", n_recombable_vars, vecsize(cur_chrom_recombable_vars));

		per_chrom_recombable_regs[i_chr] = cur_chrom_recombable_vars;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set the resampled genotypes for the current variant.
		auto hapstate_alloc_start_chrono = std::chrono::high_resolution_clock::now();

		// These two are large allocations.
		fprintf(stderr, "Allocating haplotype state information for %d recombable variants (%d all vars.)..\n", vecsize(cur_chrom_recombable_vars), vecsize(cur_chrom_var_regs));
		if (n_recombable_vars != vecsize(cur_chrom_recombable_vars))
		{
			fprintf(stderr, "Sanity check failed @ %s(%d)\n", __FILE__, __LINE__);
			exit(1);
		}
		size_t n_recombs_mem_pool_size = (size_t)(n_resampled_sample_size) * (size_t(n_recombable_vars)) + (size_t)(2);
		int* n_recombs_mem_pool = new int[n_recombs_mem_pool_size];
		memset(n_recombs_mem_pool, 0, (n_recombs_mem_pool_size) * sizeof(int));

		size_t per_hap_sampled_hap_i_mem_pool_size = (size_t)(n_recombs_mem_pool_size + n_recombs_mem_pool_size);
		int* per_hap_sampled_hap_i_mem_pool = new int[per_hap_sampled_hap_i_mem_pool_size];
		memset(per_hap_sampled_hap_i_mem_pool, 0xff, (per_hap_sampled_hap_i_mem_pool_size) * sizeof(int));
		size_t total_mem = sizeof(int) * (n_recombs_mem_pool_size + n_recombs_mem_pool_size + n_recombs_mem_pool_size);
		double total_mem_in_GB = (double)total_mem / (1024.0 * 1024.0 * 1024.0);
		fprintf(stderr, "Allocated %.3f gigabytes of haplotype state information..\n", total_mem_in_GB);

		vector<t_ansi_thread*>* allocing_threads = new vector<t_ansi_thread*>();
		for (int thread_i = 0; thread_i < n_threads; thread_i++)
		{
			void** thread_info_ptr = new void* [10];
			int* int_vals = new int[10];
			int_vals[0] = thread_i;
			int_vals[1] = n_threads;
			int_vals[2] = n_resampled_sample_size;
			thread_info_ptr[0] = int_vals;
			thread_info_ptr[1] = cur_chrom_recombable_vars;
			thread_info_ptr[2] = n_recombs_mem_pool;
			thread_info_ptr[3] = per_hap_sampled_hap_i_mem_pool;

			t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_hapstat_transition_memalloc, thread_info_ptr);
			allocing_threads->push_back(cur_thread);
			cur_thread->run_thread();
		} // thread_i loop.

		fprintf(stderr, "Waiting for %d allocing threads..\n", vecsize(allocing_threads));
		for (int thread_i = 0; thread_i < vecsize(allocing_threads); thread_i++)
		{
			allocing_threads->at(thread_i)->wait_thread();
		} // thread_i loop.

		auto hapstate_alloc_end_chrono = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> hapstate_alloc_duration = hapstate_alloc_end_chrono - hapstate_alloc_start_chrono;
		fprintf(stderr, "Finished hap-state allocations in %.4f seconds.\n", hapstate_alloc_duration.count());
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	} // chromosome index.

	t_restr_annot_region_list* restr_recombable_regs = new t_restr_annot_region_list();
	restr_recombable_regs->chr_ids = restr_geno_regs->chr_ids;
	restr_recombable_regs->regions_per_chrom = per_chrom_recombable_regs;

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		t_string::print_padded_string(stderr, '\r', 80, "Starting %d. re-sampling thread..", thread_i);

		void** thread_ptrs_list = new void* [20];

		int* int_vals = new int[20];
		int_vals[0] = thread_i;
		int_vals[1] = n_threads;
		int_vals[2] = start_pos;
		int_vals[3] = end_pos;
		int_vals[4] = n_original_haps;
		int_vals[5] = n_resampled_sample_size;
		thread_ptrs_list[0] = int_vals;

		double* dbl_vals = new double[20];
		dbl_vals[0] = N_e;
		dbl_vals[1] = allele_error_prob;
		dbl_vals[2] = segment_length_cutoff_bp;
		dbl_vals[3] = segment_length_cutoff_cM;
		dbl_vals[4] = segment_length_cutoff_in_var_number;
		thread_ptrs_list[1] = dbl_vals;

		thread_ptrs_list[2] = restr_recombable_regs;

		thread_ptrs_list[3] = recombination_rate_dir;

		thread_ptrs_list[4] = haplotype_indices_2_sample;

		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_haplotype_state_only_resampling, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Started %d/%d threads; waiting..\n", (int)resampling_threads->size(), n_threads);
	auto hapstate_sampling_start_chrono = std::chrono::high_resolution_clock::now();
	for (int thread_i = 0; thread_i < (int)resampling_threads->size(); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		t_string::print_padded_string(stderr, '\r', 100, "%d. thread finished.", thread_i);
	} // thread_i waiting loop.
	auto hapstate_sampling_end_chrono = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> hapstate_resampling_duration = hapstate_sampling_end_chrono - hapstate_sampling_start_chrono;
	fprintf(stderr, "Finished hap-state re-sampling in %.4f seconds.\n", hapstate_resampling_duration.count());

	///////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Saving re-sampled genotypes and information.\n");

	char resampled_sample_ids_fp[1000];
	sprintf(resampled_sample_ids_fp, "%s_resampled_sample_ids.list", op_prefix);
	FILE* f_samples = open_f(resampled_sample_ids_fp, "w");
	vector<char*>* resampled_sample_ids = new vector<char*>();
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		fprintf(f_samples, "sample_%d\n", i_s);

		char cur_sample_id[100];
		sprintf(cur_sample_id, "sample_%d\n", i_s);
		resampled_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
	} // i_s loop.
	close_f(f_samples, NULL);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Save the recombination pattern signals.
	if (save_recombination_patterns > 0 &&
		vecsize(sample_ids) < 5000)
	{
		for (int i_chr = 0; i_chr < vecsize(restr_recombable_regs->chr_ids); i_chr++)
		{
			if (vecsize(restr_recombable_regs->regions_per_chrom[i_chr]) < 500000)
			{
				vector<t_annot_region*>* resamp_info_regs = new vector<t_annot_region*>();
				char resampling_patterns_sigbed_fp[1000];
				sprintf(resampling_patterns_sigbed_fp, "%s_resampling_pattern_signal_var_regs_%s.sigbed.gz", op_prefix, restr_recombable_regs->chr_ids->at(i_chr));

				vector<t_annot_region*>* cur_chr_recombable_regs = restr_recombable_regs->regions_per_chrom[i_chr];
				for (int i_var = 0; i_var < vecsize(cur_chr_recombable_regs); i_var++)
				{
					void** cur_var_info = (void**)(cur_chr_recombable_regs->at(i_var)->data);
					t_annot_region* dup_reg = duplicate_region(cur_chr_recombable_regs->at(i_var));
					void** dup_reg_info = new void* [5];
					dup_reg_info[0] = cur_var_info[2];

					if (dup_reg->name == NULL)
					{
						char name_str[100];
						sprintf(name_str, "recombable_%s_%d", restr_recombable_regs->chr_ids->at(i_chr), cur_chr_recombable_regs->at(i_var)->start);
						dup_reg->name = t_string::copy_me_str(name_str);
					}

					dup_reg->strand = '+';

					dup_reg->data = dup_reg_info;
					resamp_info_regs->push_back(dup_reg);
				} // i_var loop.
				save_resampling_pattern_signal_var_regs(resamp_info_regs, vecsize(sample_ids), n_resampled_sample_size, N_e_2_n_ref_haplotypes, resampling_patterns_sigbed_fp);
			}
		} // i_chr loop.
	} // saving check flag.
	///////////////////////////////////////////////////////////////////////////////////////////
	// Save the # of recombinations at each variant.
	char per_var_n_recombs_fp[1000];
	sprintf(per_var_n_recombs_fp, "%s_per_var_n_recombs.txt", op_prefix);
	FILE* f_per_pos_n_recombs = open_f(per_var_n_recombs_fp, "w");
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		// Write n_recombs only for recombable variants.
		if (haplocoded_geno_regs->at(i_reg)->score !=  (unsigned int)(-1))
		{
			void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);
			int* cur_var_n_recombs_per_sample = (int*)(cur_reg_data[1]);

			double cur_var_n_recombs = 0;
			for (int i_s = 0; i_s < (int)resampled_sample_ids->size(); i_s++)
			{
				cur_var_n_recombs += cur_var_n_recombs_per_sample[i_s];
			} // i_s loop.

			fprintf(f_per_pos_n_recombs, "%s\t%d\t%d\t%s\t%.1f\t+\n",
				haplocoded_geno_regs->at(i_reg)->chrom,
				translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				haplocoded_geno_regs->at(i_reg)->name,
				cur_var_n_recombs);
		}
		else
		{
			double cur_var_n_recombs = 0;
			fprintf(f_per_pos_n_recombs, "%s\t%d\t%d\t%s\t%.1f\t+\n",
				haplocoded_geno_regs->at(i_reg)->chrom,
				translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				haplocoded_geno_regs->at(i_reg)->name,
				cur_var_n_recombs);
		}
	} // i_reg loop.
	fclose(f_per_pos_n_recombs);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Start saving the genotypes.
	for (int i_chr = 0; i_chr < vecsize(restr_geno_regs->chr_ids); i_chr++)
	{
		save_resampled_genotypes_per_recombination_patterns_multithreaded(restr_geno_regs->regions_per_chrom[i_chr], restr_recombable_regs->regions_per_chrom[i_chr], sample_ids, n_resampled_sample_size, use_haplocoded_genotypes, n_threads, op_prefix);
		break;
	} // i_chr loop.

	auto overall_resampling_end_chrono = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> overall_resampling_duration = overall_resampling_end_chrono - overall_resampling_start_chrono;

	fprintf(stderr, "All resampling finished in %.4f seconds.\n", overall_resampling_duration.count());
}

