#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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

#include <vector>
#include <algorithm>

using namespace std;

const bool __DUMP_HAPRESAMPLING_MSGS__ = false;

vector<char*>* extract_geno_kmers_per_genome(vector<t_annot_region*>* genocoded_var_regs, vector<char*>* genocoded_panel_sample_ids, int win_start, int win_end)
{
	vector<char*>* per_subj_genocoded_kmers = new vector<char*>();

	for (int i_s = 0; i_s < vecsize(genocoded_panel_sample_ids); i_s++)
	{
		// For this subject, do comparison to all k-mers; but we need to compare only hets???
		char* cur_genocoded_subject_geno_kmer = new char[win_end - win_start + 3];
		for (int j_var = win_start; j_var <= win_end; j_var++)
		{
			void** genocoded_var_reg_info = (void**)(genocoded_var_regs->at(j_var)->data);
			char* geno_sig = (char*)(genocoded_var_reg_info[0]);

			cur_genocoded_subject_geno_kmer[j_var - win_start] = geno_sig[i_s];
		} // j_var loop.

		// Finish the kmer.
		cur_genocoded_subject_geno_kmer[win_end + 1 - win_start] = -1;

		per_subj_genocoded_kmers->push_back(cur_genocoded_subject_geno_kmer);
	} // i_s loop.

	return(per_subj_genocoded_kmers);
} // extract_geno_kmers_per_genome option.


static void* thread_callback_save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns(void* __thread_info_ptr)
{
	void** cur_thread_params = (void**)(__thread_info_ptr);

	//int* int_val = new int[10];
	int* int_val = (int*)(cur_thread_params[0]);
	int i_t = int_val[0];
	int n_threads = int_val[1];
	int n_resampled_sample_size = int_val[2];
	int n_vic_tag_vars = int_val[3];
	int start_target_region = int_val[4];
	int end_target_region = int_val[5];

	double* dbl_val = (double*)(cur_thread_params[1]);
	double n_original_haps = dbl_val[0];
	double rel_N_e = dbl_val[1];

	t_restr_annot_region_list* restr_per_var_resampling_patterns = (t_restr_annot_region_list*)(cur_thread_params[2]);

	//vector<t_annot_region*>* all_per_var_resampling_patterns = (vector<t_annot_region*>*)(cur_thread_params[3]);

	t_restr_annot_region_list* restr_tag_target_var_regs = (t_restr_annot_region_list*)(cur_thread_params[4]);

	char* op_dir = (char*)(cur_thread_params[5]);

	for (int i_chr = 0; i_chr < vecsize(restr_per_var_resampling_patterns->chr_ids); i_chr++)
	{
		// This is the list of sorted regions.
		vector<t_annot_region*>* per_var_resampling_patterns = restr_per_var_resampling_patterns->regions_per_chrom[i_chr];

		// Get the tag/target variants.
		int tag_target_chr_i = t_string::get_i_str(restr_tag_target_var_regs->chr_ids, restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		vector<t_annot_region*>* tag_target_regs = restr_tag_target_var_regs->regions_per_chrom[tag_target_chr_i];

		for (int target_var_i = 0; target_var_i < vecsize(tag_target_regs); target_var_i++)
		{
			if (target_var_i % n_threads != i_t)
			{
				continue;
			}

			t_string::print_padded_string(stderr, '\r', 100, "@ %d/%d. target", target_var_i, vecsize(tag_target_regs));

			// Skip this if we are to the left of start region.
			if (tag_target_regs->at(target_var_i)->end < start_target_region)
			{
				continue;
			}

			// Stop processing if we are past the end.
			if (tag_target_regs->at(target_var_i)->start > end_target_region)
			{
				break;
			}

			// Find the next target variant.
			if (tag_target_regs->at(target_var_i)->dbl_score == 0)
			{
				continue;
			}
			else
			{
				if (target_var_i % 1000 == 0)
				{
					fprintf(stderr, "Saving scores for %s:%d\n",
						tag_target_regs->at(target_var_i)->chrom, tag_target_regs->at(target_var_i)->start);
				}

				// Find the tags around this target:
				vector<t_annot_region*>* cur_target_predicting_tag_regions = new vector<t_annot_region*>();
				int n_left_tags = 0;
				for (int j_reg = target_var_i;
					j_reg >= 0 && (n_left_tags < n_vic_tag_vars);
					j_reg--)
				{
					if (tag_target_regs->at(j_reg)->dbl_score == 0)
					{
						//cur_target_predicting_tag_regions->push_back(tag_target_regs->at(j_reg));
						cur_target_predicting_tag_regions->push_back(duplicate_region(tag_target_regs->at(j_reg)));
						n_left_tags++;
					}
				} // j_reg loop.

				int n_right_tags = 0;
				for (int j_reg = target_var_i;
					j_reg < vecsize(tag_target_regs) && (n_right_tags < n_vic_tag_vars);
					j_reg++)
				{
					if (tag_target_regs->at(j_reg)->dbl_score == 0)
					{
						//cur_target_predicting_tag_regions->push_back(tag_target_regs->at(j_reg));
						cur_target_predicting_tag_regions->push_back(duplicate_region(tag_target_regs->at(j_reg)));

						n_right_tags++;
					}
				} // j_reg loop.

				if (__DUMP_HAPRESAMPLING_MSGS__)
				{
					fprintf(stderr, "Found %d predicting tag regions (L/R)=(%d/%d)\n", vecsize(cur_target_predicting_tag_regions), n_left_tags, n_right_tags);
				}

				double*** cur_tag_hap_recomb_scores = score_resampled_haplotypes_on_target_regions_per_ProxyTyper_sampling_patterns(cur_target_predicting_tag_regions,
					per_var_resampling_patterns,
					rel_N_e,
					n_resampled_sample_size,
					n_original_haps);

				double** cur_tag_hap_scores = cur_tag_hap_recomb_scores[0];
				double** cur_tag_hap_n_recombs = cur_tag_hap_recomb_scores[1];

				// Save the scores.
				char scores_fp[1000];
				sprintf(scores_fp, "%s/%s_%d_%d.hap_scores", op_dir, tag_target_regs->at(target_var_i)->chrom, tag_target_regs->at(target_var_i)->start,
					tag_target_regs->at(target_var_i)->end);
				FILE* f_cur_target_scores = open_f(scores_fp, "w");
				for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
				{
					for (int i_hap = 0; i_hap < 2; i_hap++)
					{
						fprintf(f_cur_target_scores, "%d\t%d\t%.5f\t%.0f\n", i_s, i_hap, cur_tag_hap_scores[i_s][i_hap], cur_tag_hap_n_recombs[i_s][i_hap]);
					} // i_hap loop.
				} // i_s loop.
				close_f(f_cur_target_scores, scores_fp);
			}
		} // target_var_i loop.
	} // i_chr loop.

	return(NULL);
} // thread_callback_save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns function.

void save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns_multithreaded(char* haplocoded_tag_genotype_matrix_fp, char* haplocoded_target_genotype_matrix_fp, char* sample_ids_list_fp,
	char* resampling_hap_info_sigbed_fp,
	char* recombination_rate_dir,
	int n_vic_tag_vars,
	int start_target_region,
	int end_target_region,
	int n_threads,
	char* op_dir)
{
	fprintf(stderr, "Calculating haplotypes scores using %d threads..\n", n_threads);

	vector<t_annot_region*>* haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(haplocoded_tag_genotype_matrix_fp, sample_ids_list_fp);
	sort(haplocoded_tag_geno_regs->begin(), haplocoded_tag_geno_regs->end(), sort_regions);

	vector<t_annot_region*>* haplocoded_target_geno_regs = load_variant_signal_regions_wrapper(haplocoded_target_genotype_matrix_fp, sample_ids_list_fp);
	sort(haplocoded_target_geno_regs->begin(), haplocoded_target_geno_regs->end(), sort_regions);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n",
		vecsize(haplocoded_tag_geno_regs),
		vecsize(sample_ids));

	double n_original_haps = 2.0 * sample_ids->size();

	//t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);
	vector<char*>* chr_ids = get_chr_ids(haplocoded_tag_geno_regs);
	if (vecsize(chr_ids) != 1)
	{
		fprintf(stderr, "Found multiple chromosomes on tag genotypes matrix.\n");
		exit(1);
	}

	for (int i_reg = 0; i_reg < vecsize(haplocoded_tag_geno_regs); i_reg++)
	{
		haplocoded_tag_geno_regs->at(i_reg)->dbl_score = 0;
	} // i_reg loop.

	for (int i_reg = 0; i_reg < vecsize(haplocoded_target_geno_regs); i_reg++)
	{
		haplocoded_target_geno_regs->at(i_reg)->dbl_score = 1;
	} // i_reg loop.

	vector<t_annot_region*>* all_tag_target_regs = new vector<t_annot_region*>();
	all_tag_target_regs->insert(all_tag_target_regs->end(), haplocoded_tag_geno_regs->begin(), haplocoded_tag_geno_regs->end());
	all_tag_target_regs->insert(all_tag_target_regs->end(), haplocoded_target_geno_regs->begin(), haplocoded_target_geno_regs->end());

	t_restr_annot_region_list* restr_tag_target_var_regs = restructure_annot_regions(all_tag_target_regs);

	////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loading resampling patterns.\n");
	double rel_N_e = 0;
	int n_resampled_sample_size = 0;
	int n_original_sample_size = 0;

	vector<t_annot_region*>* all_per_var_resampling_patterns = load_resampling_pattern_signal_var_regs(resampling_hap_info_sigbed_fp, rel_N_e, n_original_sample_size, n_resampled_sample_size);
	fprintf(stderr, "Loaded resampling patterns for %d variants with rel_N_e=%.3f; n_original_sample_size=%d; n_resampled_sample_size=%d\n", vecsize(all_per_var_resampling_patterns),
		rel_N_e, n_original_sample_size, n_resampled_sample_size);

	if (vecsize(sample_ids) != n_original_sample_size)
	{
		fprintf(stderr, "Sanity check failed: %d vs %d on the original sample sizes.\n",
			vecsize(sample_ids), n_original_sample_size);

		exit(1);
	}

	t_restr_annot_region_list* restr_per_var_resampling_patterns = restructure_annot_regions(all_per_var_resampling_patterns);

	fprintf(stderr, "Assigning genetic distances/recombination rates..\n");
	for (int i_chr = 0; i_chr < vecsize(restr_per_var_resampling_patterns->chr_ids); i_chr++)
	{
		// This is the list of sorted regions.
		vector<t_annot_region*>* per_var_resampling_patterns = restr_per_var_resampling_patterns->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates on chromosome %s.\n", restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		for (int i_reg = 0; i_reg < (int)per_var_resampling_patterns->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(per_var_resampling_patterns->at(i_reg), cur_chrom_recomb_regs);
			per_var_resampling_patterns->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.
	} // i_chr loop.


	vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
	for (int i_t = 0; i_t < n_threads; i_t++)
	{
		void** cur_thread_params = new void* [10];

		// Integer:
		int* int_val = new int[10];
		int_val[0] = i_t;
		int_val[1] = n_threads;
		int_val[2] = n_resampled_sample_size;
		int_val[3] = n_vic_tag_vars;
		int_val[4] = start_target_region;
		int_val[5] = end_target_region;
		cur_thread_params[0] = int_val;

		// double:
		double* dbl_val = new double[10];
		dbl_val[0] = n_original_haps;
		dbl_val[1] = rel_N_e;
		cur_thread_params[1] = dbl_val;

		cur_thread_params[2] = restr_per_var_resampling_patterns;

		cur_thread_params[3] = all_per_var_resampling_patterns;

		cur_thread_params[4] = restr_tag_target_var_regs;

		cur_thread_params[5] = op_dir;

		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns, 
														cur_thread_params);
		threads->push_back(cur_thread);
		cur_thread->run_thread();
	} // i_t loop.

	for (int i_t = 0; i_t < n_threads; i_t++)
	{
		threads->at(i_t)->wait_thread();
		t_string::print_padded_string(stderr, '\r', 100, "Thread %d finished.", i_t);

	} // i_t loop.

	fprintf(stderr, "\n");
} // save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns_multithreaded function.

// This is aimed at using the scores for training purposes.
void save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns(char* haplocoded_tag_genotype_matrix_fp, char* haplocoded_target_genotype_matrix_fp, char* sample_ids_list_fp,
	char* resampling_hap_info_sigbed_fp, 
	char* recombination_rate_dir,
	int n_vic_tag_vars,
	int start_target_region, 
	int end_target_region,
	char* op_dir)
{
	vector<t_annot_region*>* haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(haplocoded_tag_genotype_matrix_fp, sample_ids_list_fp);
	sort(haplocoded_tag_geno_regs->begin(), haplocoded_tag_geno_regs->end(), sort_regions);

	vector<t_annot_region*>* haplocoded_target_geno_regs = load_variant_signal_regions_wrapper(haplocoded_target_genotype_matrix_fp, sample_ids_list_fp);
	sort(haplocoded_target_geno_regs->begin(), haplocoded_target_geno_regs->end(), sort_regions);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n",
		vecsize(haplocoded_tag_geno_regs),
		vecsize(sample_ids));

	double n_original_haps = 2.0 * sample_ids->size();

	//t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);
	vector<char*>* chr_ids = get_chr_ids(haplocoded_tag_geno_regs);
	if (vecsize(chr_ids) != 1)
	{
		fprintf(stderr, "Found multiple chromosomes on tag genotypes matrix.\n");
		exit(1);
	}

	for (int i_reg = 0; i_reg < vecsize(haplocoded_tag_geno_regs); i_reg++)
	{
		haplocoded_tag_geno_regs->at(i_reg)->dbl_score = 0;
	} // i_reg loop.

	for (int i_reg = 0; i_reg < vecsize(haplocoded_target_geno_regs); i_reg++)
	{
		haplocoded_target_geno_regs->at(i_reg)->dbl_score = 1;
	} // i_reg loop.

	vector<t_annot_region*>* all_tag_target_regs = new vector<t_annot_region*>();
	all_tag_target_regs->insert(all_tag_target_regs->end(), haplocoded_tag_geno_regs->begin(), haplocoded_tag_geno_regs->end());
	all_tag_target_regs->insert(all_tag_target_regs->end(), haplocoded_target_geno_regs->begin(), haplocoded_target_geno_regs->end());

	t_restr_annot_region_list* restr_tag_target_var_regs = restructure_annot_regions(all_tag_target_regs);

	////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loading resampling patterns.\n");
	double rel_N_e = 0;
	int n_resampled_sample_size = 0;
	int n_original_sample_size = 0;

	vector<t_annot_region*>* all_per_var_resampling_patterns = load_resampling_pattern_signal_var_regs(resampling_hap_info_sigbed_fp, rel_N_e, n_original_sample_size, n_resampled_sample_size);
	fprintf(stderr, "Loaded resampling patterns for %d variants with rel_N_e=%.3f; n_original_sample_size=%d; n_resampled_sample_size=%d\n", vecsize(all_per_var_resampling_patterns),
			rel_N_e, n_original_sample_size, n_resampled_sample_size);

	if (vecsize(sample_ids) != n_original_sample_size)
	{
		fprintf(stderr, "Sanity check failed: %d vs %d on the original sample sizes.\n", 
				vecsize(sample_ids), n_original_sample_size);

		exit(1);
	}

	t_restr_annot_region_list* restr_per_var_resampling_patterns = restructure_annot_regions(all_per_var_resampling_patterns);

	fprintf(stderr, "Assigning genetic distances/recombination rates..\n");
	for (int i_chr = 0; i_chr < vecsize(restr_per_var_resampling_patterns->chr_ids); i_chr++)
	{
		// This is the list of sorted regions.
		vector<t_annot_region*>* per_var_resampling_patterns = restr_per_var_resampling_patterns->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates on chromosome %s.\n", restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		for (int i_reg = 0; i_reg < (int)per_var_resampling_patterns->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(per_var_resampling_patterns->at(i_reg), cur_chrom_recomb_regs);
			per_var_resampling_patterns->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Get the tag/target variants.
		int tag_target_chr_i = t_string::get_i_str(restr_tag_target_var_regs->chr_ids, restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		vector<t_annot_region*>* tag_target_regs = restr_tag_target_var_regs->regions_per_chrom[tag_target_chr_i];

		for (int target_var_i = 0; target_var_i < vecsize(tag_target_regs); target_var_i++)
		{
			t_string::print_padded_string(stderr, '\r', 100, "@ %d/%d. target", target_var_i, vecsize(tag_target_regs));

			// Skip this if we are to the left of start region.
			if (tag_target_regs->at(target_var_i)->end < start_target_region)
			{
				continue;
			}

			// Stop processing if we are past the end.
			if (tag_target_regs->at(target_var_i)->start > end_target_region)
			{
				break;
			}

			// Find the next target variant.
			if (tag_target_regs->at(target_var_i)->dbl_score == 0)
			{
				continue;
			}
			else
			{
				if (target_var_i % 1000 == 0)
				{
					fprintf(stderr, "Saving scores for %s:%d\n",
						tag_target_regs->at(target_var_i)->chrom, tag_target_regs->at(target_var_i)->start);
				}

				// Find the tags around this target:
				vector<t_annot_region*>* cur_target_predicting_tag_regions = new vector<t_annot_region*>();
				int n_left_tags = 0;
				for (int j_reg = target_var_i;
					j_reg >= 0 && (n_left_tags < n_vic_tag_vars);
					j_reg--)
				{
					if (tag_target_regs->at(j_reg)->dbl_score == 0)
					{
						//cur_target_predicting_tag_regions->push_back(tag_target_regs->at(j_reg));
						cur_target_predicting_tag_regions->push_back(duplicate_region(tag_target_regs->at(j_reg)));
						n_left_tags++;
					}
				} // j_reg loop.

				int n_right_tags = 0;
				for (int j_reg = target_var_i;
					j_reg < vecsize(tag_target_regs) && (n_right_tags < n_vic_tag_vars);
					j_reg++)
				{
					if (tag_target_regs->at(j_reg)->dbl_score == 0)
					{
						//cur_target_predicting_tag_regions->push_back(tag_target_regs->at(j_reg));
						cur_target_predicting_tag_regions->push_back(duplicate_region(tag_target_regs->at(j_reg)));

						n_right_tags++;
					}
				} // j_reg loop.

				if (__DUMP_HAPRESAMPLING_MSGS__)
				{
					fprintf(stderr, "Found %d predicting tag regions (L/R)=(%d/%d)\n", vecsize(cur_target_predicting_tag_regions), n_left_tags, n_right_tags);
				}

				double*** cur_tag_hap_recomb_scores = score_resampled_haplotypes_on_target_regions_per_ProxyTyper_sampling_patterns(cur_target_predicting_tag_regions,
					per_var_resampling_patterns,
					rel_N_e,
					n_resampled_sample_size,
					n_original_haps);

				double** cur_tag_hap_scores = cur_tag_hap_recomb_scores[0];
				double** cur_tag_hap_n_recombs = cur_tag_hap_recomb_scores[1];

				// Save the scores.
				char scores_fp[1000];
				sprintf(scores_fp, "%s/%s_%d_%d.hap_scores", op_dir, tag_target_regs->at(target_var_i)->chrom, tag_target_regs->at(target_var_i)->start,
					tag_target_regs->at(target_var_i)->end);
				FILE* f_cur_target_scores = open_f(scores_fp, "w");
				for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
				{
					for (int i_hap = 0; i_hap < 2; i_hap++)
					{
						fprintf(f_cur_target_scores, "%d\t%d\t%.5f\t%.0f\n", i_s, i_hap, cur_tag_hap_scores[i_s][i_hap], cur_tag_hap_n_recombs[i_s][i_hap]);
					} // i_hap loop.
				} // i_s loop.
				close_f(f_cur_target_scores, scores_fp);
			}
		} // target_var_i loop.
	} // i_chr loop.
} // save_per_imputation_target_resampled_haplotype_scores_per_sampling_patterns_files function.

// Note that this function does not rely on ProxyTyper's sampling patterns, this functionality exists in variation_tools now as a forked copy.
double*** score_resampled_haplotypes_on_target_regions_per_ProxyTyper_sampling_patterns(vector<t_annot_region*>* target_regions,
	vector<t_annot_region*>* all_per_var_resampling_patterns,
	double rel_N_e,
	int n_resampled_sample_size,
	int n_original_haps)
{
	vector<t_annot_region*>* intersects = intersect_annot_regions(all_per_var_resampling_patterns, target_regions, false);
	vector<t_annot_region*>* target_region_vars_resampling_patterns = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		target_region_vars_resampling_patterns->push_back(int_info->src_reg);
	} // i_int loop.

	//fprintf(stderr, "Found %d target variants with resampling patterns.\n", vecsize(target_region_vars_resampling_patterns));

	double** per_subj_per_hap_LOD_scores = new double* [n_resampled_sample_size + 2];
	double** per_subj_n_recombs = new double* [n_resampled_sample_size + 2];
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		per_subj_n_recombs[i_s] = new double[2];
		per_subj_n_recombs[i_s][0] = 0;
		per_subj_n_recombs[i_s][1] = 0;

		per_subj_per_hap_LOD_scores[i_s] = new double[2];
		per_subj_per_hap_LOD_scores[i_s][0] = 0;
		per_subj_per_hap_LOD_scores[i_s][1] = 0;
	} // i_s loop.

	// This is necessary to correctly track the transitions and identify segments.
	t_restr_annot_region_list* restr_per_var_resampling_patterns = restructure_annot_regions(target_region_vars_resampling_patterns);

	for (int i_chr = 0; i_chr < vecsize(restr_per_var_resampling_patterns->chr_ids); i_chr++)
	{
		// This is the list of sorted regions.
		vector<t_annot_region*>* per_var_resampling_patterns = restr_per_var_resampling_patterns->regions_per_chrom[i_chr];

		//char cur_chr_recombination_rate_fp[1000];
		//sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		//vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		//if (cur_chrom_recomb_regs == NULL)
		//{
		//	fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		//	exit(1);
		//}

		//// Assign the recomb rates.
		//fprintf(stderr, "Setting recombination rates on chromosome %s.\n", restr_per_var_resampling_patterns->chr_ids->at(i_chr));
		//for (int i_reg = 0; i_reg < (int)per_var_resampling_patterns->size(); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(per_var_resampling_patterns->at(i_reg), cur_chrom_recomb_regs);
		//	per_var_resampling_patterns->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		// Buffer the self and other transition probabilities.
		double* per_var_self_probs = new double[vecsize(per_var_resampling_patterns) + 2];
		double* per_var_other_probs = new double[vecsize(per_var_resampling_patterns) + 2];

		per_var_self_probs[0] = 1.0;
		per_var_other_probs[0] = 0.0;

		//fprintf(stderr, "Setting up transition probabilities..\n");

		double N_e = n_original_haps * rel_N_e;

		// Get the recombination rate between these positions.
		for (int i_reg = 1; i_reg < vecsize(per_var_resampling_patterns); i_reg++)
		{
			double r_m = fabs(per_var_resampling_patterns->at(i_reg)->dbl_score - per_var_resampling_patterns->at(i_reg - 1)->dbl_score);
			double rho_m = 4 * N_e * r_m;

			double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

			double other_prob = tau_m / n_original_haps;
			double self_prob = (1 - tau_m) + (tau_m / n_original_haps);

			per_var_self_probs[i_reg] = self_prob;
			per_var_other_probs[i_reg] = other_prob;
		} // i_reg loop.

		// For each resampled subject, go over ever
		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				double cur_hap_baseline_LOD_score = 0;
				double n_total_recombs = 0;

				int cur_resampled_sample_i = -1;
				int cur_resampled_hap_i = -1;
				//int cur_resampled_seg_start = -1;
				for (int i_var = 1; i_var < vecsize(per_var_resampling_patterns); i_var++)
				{
					// Calculate the LOD between self vs transition from var-{i-1} to var-i.
					double cur_var_other_LOD = xlog((n_original_haps-1) * per_var_other_probs[i_var] / per_var_self_probs[i_var]);
					double cur_var_self_LOD = 0.0;

					void** cur_var_info = (void**)(per_var_resampling_patterns->at(i_var)->data);
					int** per_hap_resampled_hap_i = (int**)(cur_var_info[0]);

					int cur_haplo_state = per_hap_resampled_hap_i[i_hap][i_s];
					int sample_i = (cur_haplo_state - (cur_haplo_state % 2)) / 2;
					int sampled_hap_i = (cur_haplo_state % 2);

					if (cur_resampled_sample_i == sample_i &&
						cur_resampled_hap_i == sampled_hap_i)
					{
						// We don't do anything in this case.
						cur_hap_baseline_LOD_score += cur_var_self_LOD;
					}
					else
					{
						if (cur_resampled_sample_i != -1)
						{
							cur_hap_baseline_LOD_score += cur_var_other_LOD;
							n_total_recombs += 1.0;

							/*	
							fprintf(f_per_hap_resampled_hap_i, "%d\t%d\t%s[%d]\t%s\t%d\t%d\t%d\t%d\n",
								i_s, i_hap,
								generating_sample_ids->at(cur_resampled_sample_i), cur_resampled_hap_i,
								per_var_resampling_patterns->at(i_var - 1)->chrom, cur_resampled_seg_start, i_var - 1,
								per_var_resampling_patterns->at(cur_resampled_seg_start)->start,
								per_var_resampling_patterns->at(i_var - 1)->end);
							*/
						}

						cur_resampled_sample_i = sample_i;
						cur_resampled_hap_i = sampled_hap_i;
						//cur_resampled_seg_start = i_var;
					}
				} // i_var loop.

				//fprintf(f_per_hap_resampled_hap_i, "\n");
				//fprintf(f_hap_score_op, "%d\t%d\t%.4f\n", i_s, i_hap, cur_hap_baseline_LOD_score);
				per_subj_per_hap_LOD_scores[i_s][i_hap] = cur_hap_baseline_LOD_score;
				per_subj_n_recombs[i_s][i_hap] = n_total_recombs;
			} // i_hap loop.
		} // i_s loop.
	} // i_chr loop.
	//close_f(f_hap_score_op, hap_score_op_fp);

	double*** per_subj_recomb_info = new double** [5];
	per_subj_recomb_info[0] = per_subj_per_hap_LOD_scores;
	per_subj_recomb_info[1] = per_subj_n_recombs;


	return(per_subj_recomb_info);
}


// Following write the sampling patterns into a parseable text file from all variants.
void summarize_sampled_segments_per_resampled_haplotype_info(char* resampling_hap_info_sigbed_fp, char* resampling_sample_list_fp, char* generating_sample_list_fp, char* op_fp)
{
	//vector<t_annot_region*>* all_per_var_resampling_patterns = load_resampling_pattern_signal_var_regs(resampling_hap_info_sigbed_fp);
	double read_rel_N_e = 0;
	int read_n_original_sample_size = 0;
	int read_n_resampled_sample_size = 0;
	vector<t_annot_region*>* all_per_var_resampling_patterns = load_resampling_pattern_signal_var_regs(resampling_hap_info_sigbed_fp, read_rel_N_e, read_n_original_sample_size, read_n_resampled_sample_size);
	fprintf(stderr, "Loaded resampling patterns for %d variants.\n", vecsize(all_per_var_resampling_patterns));

	// This is necessary to correctly track the transitions and identify segments.
	t_restr_annot_region_list* restr_per_var_resampling_patterns = restructure_annot_regions(all_per_var_resampling_patterns);

	vector<char*>* resampled_sample_ids = buffer_file(resampling_sample_list_fp);
	fprintf(stderr, "Loaded %d resampled subject ids.\n", (int)resampled_sample_ids->size());

	vector<char*>* generating_sample_ids = buffer_file(generating_sample_list_fp);
	fprintf(stderr, "Loaded %d generating panel subject ids.\n", (int)generating_sample_ids->size());

	int n_resampled_sample_size = (int)resampled_sample_ids->size();

	FILE* f_per_hap_resampled_hap_i = open_f(op_fp, "w");

	for (int i_chr = 0; i_chr < vecsize(restr_per_var_resampling_patterns->chr_ids); i_chr++)
	{
		vector<t_annot_region*>* per_var_resampling_patterns = restr_per_var_resampling_patterns->regions_per_chrom[i_chr];

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				int cur_resampled_sample_i = -1;
				int cur_resampled_hap_i = -1;
				int cur_resampled_seg_start = -1;
				for (int i_var = 0; i_var < vecsize(per_var_resampling_patterns); i_var++)
				{
					void** cur_var_info = (void**)(per_var_resampling_patterns->at(i_var)->data);
					int** per_hap_resampled_hap_i = (int**)(cur_var_info[0]);

					int cur_haplo_state = per_hap_resampled_hap_i[i_hap][i_s];
					int sample_i = (cur_haplo_state - (cur_haplo_state % 2)) / 2;
					int sampled_hap_i = (cur_haplo_state % 2);

					if (cur_resampled_sample_i == sample_i &&
						cur_resampled_hap_i == sampled_hap_i)
					{

					}
					else
					{
						if (cur_resampled_sample_i != -1)
						{
							fprintf(f_per_hap_resampled_hap_i, "%d\t%d\t%s[%d]\t%s\t%d\t%d\t%d\t%d\n",
								i_s, i_hap,
								generating_sample_ids->at(cur_resampled_sample_i), cur_resampled_hap_i,
								per_var_resampling_patterns->at(i_var - 1)->chrom, cur_resampled_seg_start, i_var - 1,
								per_var_resampling_patterns->at(cur_resampled_seg_start)->start,
								per_var_resampling_patterns->at(i_var - 1)->end);
						}

						cur_resampled_sample_i = sample_i;
						cur_resampled_hap_i = sampled_hap_i;
						cur_resampled_seg_start = i_var;
					}
				} // i_var loop.

				//fprintf(f_per_hap_resampled_hap_i, "\n");
			} // i_hap loop.
		} // i_s loop.
	} // i_chr loop.
	fclose(f_per_hap_resampled_hap_i);
}


vector<t_annot_region*>* load_resampling_pattern_signal_var_regs(char* resampling_pattern_sig_BED_fp, double& rel_resampling_N_e, int& n_sampling_sample_size, int& n_resampled_sample_size)
{
	vector<t_annot_region*>* resamp_info_regs = new vector<t_annot_region*>();

	FILE* f_resamp_regs = open_f(resampling_pattern_sig_BED_fp, "rb");
	if (f_resamp_regs == NULL)
	{
		return(NULL);
	}

	//fwrite(&(n_vars), sizeof(int), 1, f_resamp_regs);
	//fwrite(&(n_sampling_sample_size), sizeof(int), 1, f_resamp_regs);
	//fwrite(&(n_resampled_sample_size), sizeof(int), 1, f_resamp_regs);
	//fwrite(&(rel_resampling_N_e), sizeof(double), 1, f_resamp_regs);

	int n_regs = read_bin_int(f_resamp_regs);
	n_sampling_sample_size = read_bin_int(f_resamp_regs);
	n_resampled_sample_size = read_bin_int(f_resamp_regs);

	rel_resampling_N_e = -1;
	read_bin_double_array(f_resamp_regs, &rel_resampling_N_e, 1);

	fprintf(stderr, "Loading %d regions for %d samples (%d haplotypes) from %s\n", n_regs, n_resampled_sample_size, n_resampled_sample_size * 2, resampling_pattern_sig_BED_fp);
	for (int i_reg = 0; i_reg < n_regs; i_reg++)
	{
		// 1st col: Chrom.
		int l_chrom = -1;
		//fread(&l_chrom, sizeof(int), 1, f_resamp_regs);
		l_chrom = read_bin_int(f_resamp_regs);

		char* chrom = new char[l_chrom + 2];
		memset(chrom, 0, sizeof(char) * (l_chrom + 2));
		//fread(chrom, sizeof(char), l_chrom, f_resamp_regs);
		read_bin_char_array(f_resamp_regs, chrom, l_chrom);

		// 2-3 cols: Start-end
		int reg_start = -1;
		int reg_end = -1;
		//fread(&(reg_start), sizeof(int), 1, f_resamp_regs);
		reg_start = read_bin_int(f_resamp_regs);
		//fread(&(reg_end), sizeof(int), 1, f_resamp_regs);
		reg_end = read_bin_int(f_resamp_regs);

		// 4th col: name
		int l_name_str = -1;
		//fread(&(l_name_str), sizeof(int), 1, f_resamp_regs);
		l_name_str = read_bin_int(f_resamp_regs);

		char* name_str = new char[l_name_str + 2];
		memset(name_str, 0, sizeof(char) * (l_name_str + 2));
		//fread(name_str, sizeof(char), l_name_str, f_resamp_regs);
		read_bin_char_array(f_resamp_regs, name_str, l_name_str);

		// 5th column
		unsigned int reg_score = -1;
		//fread(&(reg_score), sizeof(unsigned int), 1, f_resamp_regs);
		read_bin_uint_array(f_resamp_regs, &reg_score, 1);

		// 6th column
		char strand = -1;
		//fread(&(strand), sizeof(char), 1, f_resamp_regs);
		strand = read_bin_char(f_resamp_regs);

		// Save the resampled haplotype indices for this variant: Save haps 0 and 1 in order.
		int** per_hap_resampled_hap_i = new int* [2];
		per_hap_resampled_hap_i[0] = new int[n_resampled_sample_size + 2];
		per_hap_resampled_hap_i[1] = new int[n_resampled_sample_size + 2];

		// Write the resampled haplotype indices.
		//fread(per_hap_resampled_hap_i[0], sizeof(int), n_resampled_sample_size, f_resamp_regs);
		read_bin_int_array(f_resamp_regs, per_hap_resampled_hap_i[0], n_resampled_sample_size);
		//fread(per_hap_resampled_hap_i[1], sizeof(int), n_resampled_sample_size, f_resamp_regs);
		read_bin_int_array(f_resamp_regs, per_hap_resampled_hap_i[1], n_resampled_sample_size);

		// Allocate and set the regions, there is no translation to these as they are written and read in BED coordinate system.
		t_annot_region* new_reg = get_empty_region();
		new_reg->chrom = chrom;
		new_reg->start = reg_start;
		new_reg->end = reg_end;
		new_reg->name = name_str;
		new_reg->score = reg_score;
		new_reg->strand = strand;

		void** reg_info = new void* [3];
		reg_info[0] = per_hap_resampled_hap_i;

		new_reg->data = reg_info;

		resamp_info_regs->push_back(new_reg);
	} // i_reg loop.
	close_f(f_resamp_regs, resampling_pattern_sig_BED_fp);

	fprintf(stderr, "Loaded resampling pattern information for %d regions.\n", (int)resamp_info_regs->size());
	return(resamp_info_regs);
} // load_resampling_pattern_signal_var_regs function.

void save_resampling_pattern_signal_var_regs(vector<t_annot_region*>* haplocoded_geno_regs, int n_sampling_sample_size, int n_resampled_sample_size, double rel_resampling_N_e, const char* resampling_pattern_sig_BED_fp)
{
	fprintf(stderr, "Saving the resampled haplotype indices for %d (%d subjects) regions to %s\n", (int)haplocoded_geno_regs->size(), n_resampled_sample_size, resampling_pattern_sig_BED_fp);
	FILE* f_resamp_regs = open_f(resampling_pattern_sig_BED_fp, "wb");

	int n_vars = (int)haplocoded_geno_regs->size();
	fwrite(&(n_vars), sizeof(int), 1, f_resamp_regs);
	fwrite(&(n_sampling_sample_size), sizeof(int), 1, f_resamp_regs);
	fwrite(&(n_resampled_sample_size), sizeof(int), 1, f_resamp_regs);
	fwrite(&(rel_resampling_N_e), sizeof(double), 1, f_resamp_regs);
	
	for (int i_reg = 0; i_reg < n_vars; i_reg++)
	{
		// 1st col: Chrom.
		char* chrom = haplocoded_geno_regs->at(i_reg)->chrom;
		int l_chrom = t_string::string_length(chrom);
		fwrite(&l_chrom, sizeof(int), 1, f_resamp_regs);
		fwrite(chrom, sizeof(char), l_chrom, f_resamp_regs);

		// 2-3 cols: Start-end
		fwrite(&(haplocoded_geno_regs->at(i_reg)->start), sizeof(int), 1, f_resamp_regs);
		fwrite(&(haplocoded_geno_regs->at(i_reg)->end), sizeof(int), 1, f_resamp_regs);

		// 4th col: name
		char* name_str = haplocoded_geno_regs->at(i_reg)->name;
		int l_name_str = t_string::string_length(name_str);
		fwrite(&(l_name_str), sizeof(int), 1, f_resamp_regs);
		fwrite(name_str, sizeof(char), l_name_str, f_resamp_regs);

		// 5th column
		fwrite(&(haplocoded_geno_regs->at(i_reg)->score), sizeof(unsigned int), 1, f_resamp_regs);

		// 6th column
		fwrite(&(haplocoded_geno_regs->at(i_reg)->strand), sizeof(char), 1, f_resamp_regs);

		// Save the resampled haplotype indices for this variant: Save haps 0 and 1 in order.
		void** cur_var_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		int** per_hap_resampled_hap_i = (int**)(cur_var_info[0]);

		// Write the resampled haplotype indices.
		fwrite(per_hap_resampled_hap_i[0], sizeof(int), n_resampled_sample_size, f_resamp_regs);
		fwrite(per_hap_resampled_hap_i[1], sizeof(int), n_resampled_sample_size, f_resamp_regs);
	} // i_reg loop.
	close_f(f_resamp_regs, resampling_pattern_sig_BED_fp);

	fprintf(stderr, "Finished saving the resampled haplotype indices..\n");
} // save_resampling_pattern_signal_var_regs function.


static void* resampling_thread_callback_length_cutoff_save_recomb(void* thread_info_ptr)
{
	void** thread_ptrs_list = (void**)(thread_info_ptr);

	t_restr_annot_region_list* restr_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[0]);

	//char* recombination_rate_dir = (char*)(thread_ptrs_list[1]);

	int* thread_i_ptr = (int*)(thread_ptrs_list[2]);
	int thread_i = thread_i_ptr[0];

	int* n_threads_ptr = (int*)(thread_ptrs_list[3]);
	int n_threads = n_threads_ptr[0];

	double* N_e_ptr = (double*)(thread_ptrs_list[4]);
	double N_e = N_e_ptr[0];

	double* allele_error_ptr = (double*)(thread_ptrs_list[5]);
	double allele_error = allele_error_ptr[0];

	double* length_cutoff_in_bps_ptr = (double*)(thread_ptrs_list[6]);
	double length_cutoff_in_bps = length_cutoff_in_bps_ptr[0];

	double* length_cutoff_in_cM_ptr = (double*)(thread_ptrs_list[7]);
	double length_cutoff_in_cM = length_cutoff_in_cM_ptr[0];

	double* length_cutoff_in_var_number_ptr = (double*)(thread_ptrs_list[8]);
	double length_cutoff_in_var_number = length_cutoff_in_var_number_ptr[0];

	int* n_original_haps_ptr = (int*)(thread_ptrs_list[9]);
	int n_original_haps = n_original_haps_ptr[0];

	int* n_resampled_sample_size_ptr = (int*)(thread_ptrs_list[10]);
	int n_resampled_sample_size = n_resampled_sample_size_ptr[0];

	int* start_end_coords = (int*)(thread_ptrs_list[11]);
	int start_coord = start_end_coords[0];
	int end_coord = start_end_coords[1];

	vector<int>* haplotype_indices_2_sample = (vector<int>*)(thread_ptrs_list[12]);

	vector<char**>** per_chr_per_subj_haplotypes = (vector<char**>**)(thread_ptrs_list[13]);

	if ((int)haplotype_indices_2_sample->size() != n_original_haps)
	{
		fprintf(stderr, "Could not match the number of haplotypes: %d, %d\n", (int)haplotype_indices_2_sample->size(), (int)n_original_haps);
		exit(1);
	}

	fprintf(stderr, "Length cutoffs:\n\
Max segment in bps: %.1f\n\
Max segment in cMs: %.1f\n\
Max segment in n vars.: %.1f\n", length_cutoff_in_bps, length_cutoff_in_cM, length_cutoff_in_var_number);

	// Load the recombination rates.
	int cur_thread_seed = time(NULL) + thread_i;
	//t_rng* rng = new t_rng(cur_thread_seed);
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	fprintf(stderr, "Setting thread %d/%d: N_e: %d, error_prob: %.5f, n_ref_haps: %d, n_resampled_size: %d; Interval: [%d-%d]; random_seed: %d\n",
		thread_i,
		n_threads,
		(int)N_e,
		allele_error,
		n_original_haps,
		n_resampled_sample_size,
		start_coord, end_coord,
		cur_thread_seed);

	fprintf(stderr, "First 5 random numbers: ");
	for (int r_i = 0; r_i < 5; r_i++)
	{
		fprintf(stderr, "%.4f, ", rng->random_double_ran3());
	} // r_i loop.
	fprintf(stderr, "\n");

	for (int i_chr = 0; i_chr < (int)restr_geno_regs->chr_ids->size(); i_chr++)
	{
		//fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];

		vector<char**>* cur_chr_per_subj_haplotypes = per_chr_per_subj_haplotypes[i_chr];

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

			if (thread_i == 0)
			{
				fprintf(stderr, "Sampling Subject %d..\n", i_s);
			}

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
					double* cur_reg_n_recombs = (double*)(cur_reg_info[2]);
					//double* cur_reg_n_errors = (double*)(cur_reg_info[3]);
					int** cur_reg_per_hap_resampled_hap_i = (int**)(cur_reg_info[4]);

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
					else if (cur_chrom_var_regs->at(i_reg)->score == 1) // Make sure that tihs is a recombable variant.
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

								if (__DUMP_HAPRESAMPLING_MSGS__)
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

					// Save the resampled haplotype index, make sure this is stored here so that we do not miss it out.
					cur_reg_per_hap_resampled_hap_i[i_hap][i_s] = TRACKED_HAPLO_STATE;

					// Copy the sampled allele.
					int sample_i = (TRACKED_HAPLO_STATE - (TRACKED_HAPLO_STATE % 2)) / 2;
					int sampled_hap_i = (TRACKED_HAPLO_STATE % 2);
					char* cur_reg_sampled_geno = (char*)(cur_reg_info[1]);

					int cur_allele = cur_chr_per_subj_haplotypes->at(sample_i)[sampled_hap_i][i_reg];

					////////////////////////////////////////////
					//// Following is a sanity check to make sure we are reading the correct allele: Comment out to cut down computes.
					//char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
					//int cur_allele_check = get_allele_per_haplotype(cur_reg_geno_sig[sample_i], sampled_hap_i);
					//if (cur_allele_check != cur_allele)
					//{
					//	fprintf(stderr, "Sanity check failed: Buffered allele does not match original allele.\n");
					//	exit(1);
					//}
					////////////////////////////////////////////

					// Allelic errors are skippped.
					//if (rng->random_double_ran3() < allele_error)
					//{
					//	n_erroneous_alleles++;
					//	cur_allele = 1 - cur_allele; // Flip the allele.

					//	cur_reg_n_errors[i_s]++;
					//}

					// Copy the allele to the haplotype 
					char cur_val = (char)(cur_reg_sampled_geno[i_s]);
					cur_reg_sampled_geno[i_s] = cur_val | (cur_allele << i_hap);
				} // i_reg loop.

				int thousand_fact = 1000 / n_threads;
				if ((i_s - thread_i) % (n_threads * thousand_fact) == 0)
				{
					fprintf(stderr, "Thread %d: Re-sampled sample %d (%d): %d recombinations (%d allelic errors).\n", thread_i, i_s, i_hap, n_recombs, n_erroneous_alleles);
				}
			} // i_hap loop.
		} // i_s loop.
	} // i_chr loop.

	return NULL;
} // resampling_thread_callback function.

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
	char* op_prefix)
{
	fprintf(stderr, "%d-thread AF-prior-based Re-sampling %d recombination-aware genotypes using N_e/N_ref_hap=%.3f and allelic error=%.6f in genomic interval [%d-%d]\n",
		n_threads,
		n_resampled_sample_size,
		N_e_2_n_ref_haplotypes,
		allele_error_prob,
		start_pos, end_pos);

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

	double N_e = N_e_2_n_ref_haplotypes * n_original_haps;

	// This is used to speed-up resampling by setting a minimum distance between recombination "hotspots".
	double min_cM_delta_per_recomb = 0.1;

	fprintf(stderr, "Sampling from %d haplotypes using N_e=%.3f, allele_eps=%.3f with min cM per recombable var=%.4f.\n", (int)(n_original_haps), N_e, allele_error_prob, min_cM_delta_per_recomb);

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	//int n_resampled_sample_size = (int)(upsampling_rate * sample_ids->size());
	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	vector<char**>** per_chr_per_subj_haplotypes = new vector<char**>*[restr_geno_regs->chr_ids->size()];

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

		// This is not used for now.
		//// Assign variant clusterings based on cM distance.
		double cur_cM = -1;
		int n_recombable_vars = 0;
		for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
		{
			// IF the current variant is further away than the minimum cM cutoff, set it as a recomb. hotspot.
			if (fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_cM) > min_cM_delta_per_recomb)
			{
				cur_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
				cur_chrom_var_regs->at(i_reg)->score = 1;
				n_recombable_vars++;
			}
			else
			{
				cur_chrom_var_regs->at(i_reg)->score = 0;
			}
		} // i_reg loop.

		fprintf(stderr, "%d recombable variants\n", n_recombable_vars);

		fprintf(stderr, "Extracting per chromosome per subject haplotypes.\n");
		vector<char**>* per_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(cur_chrom_var_regs, sample_ids);
		per_chr_per_subj_haplotypes[i_chr] = per_subj_haplotypes;

		// Set the resampled genotypes for the current variant.
		for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
			void** new_reg_info = new void* [10];
			new_reg_info[0] = cur_reg_info[0];

			char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
			memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
			new_reg_info[1] = cur_reg_resampled_geno;

			double* n_recombs_per_var = new double[n_resampled_sample_size + 2];
			memset(n_recombs_per_var, 0, sizeof(double) * (n_resampled_sample_size + 2));
			new_reg_info[2] = n_recombs_per_var;

			double* n_errors_per_var = new double[n_resampled_sample_size + 2];
			memset(n_errors_per_var, 0, sizeof(double) * (n_resampled_sample_size + 2));
			new_reg_info[3] = n_errors_per_var;

			int** per_hap_sampled_hap_i = new int* [2];
			per_hap_sampled_hap_i[0] = new int[n_resampled_sample_size];
			per_hap_sampled_hap_i[1] = new int[n_resampled_sample_size];

			// Initialize all to -1, this is for sanity checking later.
			for (int i_re_s = 0; i_re_s < n_resampled_sample_size; i_re_s++)
			{
				per_hap_sampled_hap_i[0][i_re_s] = -1;
				per_hap_sampled_hap_i[1][i_re_s] = -1;
			} // i_re_s loop.
			new_reg_info[4] = per_hap_sampled_hap_i;

			// Replace the info ptr.
			cur_chrom_var_regs->at(i_reg)->data = new_reg_info;
		} // i_reg loop.
	} // chromosome index.

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting %d. re-sampling thread..\n", thread_i);

		void** thread_ptrs_list = new void* [20];

		thread_ptrs_list[0] = restr_geno_regs;

		thread_ptrs_list[1] = recombination_rate_dir;

		int* thread_i_ptr = new int[2];
		thread_i_ptr[0] = thread_i;
		thread_ptrs_list[2] = thread_i_ptr;

		int* n_threads_ptr = new int[2];
		n_threads_ptr[0] = n_threads;
		thread_ptrs_list[3] = n_threads_ptr;

		double* N_e_ptr = new double[2];
		N_e_ptr[0] = N_e;
		thread_ptrs_list[4] = N_e_ptr;

		double* allele_error_ptr = new double[2];
		allele_error_ptr[0] = allele_error_prob;
		thread_ptrs_list[5] = allele_error_ptr;

		double* length_cutoff_in_bps_ptr = new double[2];
		length_cutoff_in_bps_ptr[0] = segment_length_cutoff_bp;
		thread_ptrs_list[6] = length_cutoff_in_bps_ptr;

		double* length_cutoff_in_cM_ptr = new double[2];
		length_cutoff_in_cM_ptr[0] = segment_length_cutoff_cM;
		thread_ptrs_list[7] = length_cutoff_in_cM_ptr;

		double* length_cutoff_in_var_number_ptr = new double[2];
		length_cutoff_in_var_number_ptr[0] = segment_length_cutoff_in_var_number;
		thread_ptrs_list[8] = length_cutoff_in_var_number_ptr;

		int* n_original_haps_ptr = new int[2];
		n_original_haps_ptr[0] = n_original_haps;
		thread_ptrs_list[9] = n_original_haps_ptr;

		int* n_resampled_sample_size_ptr = new int[2];
		n_resampled_sample_size_ptr[0] = n_resampled_sample_size;
		thread_ptrs_list[10] = n_resampled_sample_size_ptr;

		int* start_end_coords = new int[2];
		start_end_coords[0] = start_pos;
		start_end_coords[1] = end_pos;
		thread_ptrs_list[11] = start_end_coords;

		thread_ptrs_list[12] = haplotype_indices_2_sample;

		thread_ptrs_list[13] = per_chr_per_subj_haplotypes;

		t_ansi_thread* cur_thread = new t_ansi_thread(resampling_thread_callback_length_cutoff_save_recomb, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Started %d/%d threads; waiting.\n", (int)resampling_threads->size(), n_threads);

	for (int thread_i = 0; thread_i < (int)resampling_threads->size(); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		fprintf(stderr, "%d. thread finished.\n", thread_i);
	} // thread_i waiting loop.

	fprintf(stderr, "Saving re-sampled genotypes and information.\n");

	///////////////////////////////////////////////////////////////////////////////////////////
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
	char resampling_patterns_sigbed_fp[1000];
	sprintf(resampling_patterns_sigbed_fp, "%s_resampling_pattern_signal_var_regs.sigbed.gz", op_prefix);
	vector<t_annot_region*>* resamp_info_regs = new vector<t_annot_region*>();
	for (int i_var = 0; i_var < (int)haplocoded_geno_regs->size(); i_var++)
	{
		void** cur_var_info = (void**)(haplocoded_geno_regs->at(i_var)->data);
		t_annot_region* dup_reg = duplicate_region(haplocoded_geno_regs->at(i_var));
		void** dup_reg_info = new void* [5];
		dup_reg_info[0] = cur_var_info[4];
		dup_reg->data = dup_reg_info;
		resamp_info_regs->push_back(dup_reg);
	} // i_var loop.
	save_resampling_pattern_signal_var_regs(resamp_info_regs, vecsize(sample_ids), n_resampled_sample_size, N_e_2_n_ref_haplotypes, resampling_patterns_sigbed_fp);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Save the # of recombinations at each variant.
	char per_var_n_recombs_fp[1000];
	sprintf(per_var_n_recombs_fp, "%s_per_var_n_recombs.txt", op_prefix);
	FILE* f_per_pos_n_recombs = open_f(per_var_n_recombs_fp, "w");
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		double* cur_var_n_recombs_per_sample = (double*)(cur_reg_data[2]);

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
	} // i_reg loop.
	fclose(f_per_pos_n_recombs);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Write the number of errors for each variant.
	char per_var_n_errors_fp[1000];
	sprintf(per_var_n_errors_fp, "%s_per_var_n_errors.txt", op_prefix);
	FILE* f_per_pos_n_errors = open_f(per_var_n_errors_fp, "w");
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		double* cur_var_n_errors_per_sample = (double*)(cur_reg_data[3]);
		char* cur_var_geno = (char*)(cur_reg_data[0]);

		double cur_var_n_errors = 0;
		for (int i_s = 0; i_s < (int)resampled_sample_ids->size(); i_s++)
		{
			cur_var_n_errors += cur_var_n_errors_per_sample[i_s];
		} // i_s loop.

		double total_geno = 0;
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			total_geno += get_genotype_per_haplocoded_genotype(cur_var_geno[i_s]);
		}

		fprintf(f_per_pos_n_errors, "%s\t%d\t%d\t%s\t%.1f\t%.1f\t+\n",
			haplocoded_geno_regs->at(i_reg)->chrom,
			translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			haplocoded_geno_regs->at(i_reg)->name,
			cur_var_n_errors,
			total_geno);
	} // i_reg loop.
	fclose(f_per_pos_n_errors);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Fixed the data index of sampled genotypes.
	vector<t_annot_region*>* resampled_tag_geno_var_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		t_annot_region* dup_reg = duplicate_region(haplocoded_geno_regs->at(i_reg));
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);

		void** dup_reg_info = new void* [2];
		dup_reg_info[0] = cur_reg_data[1];
		dup_reg->data = dup_reg_info;
		resampled_tag_geno_var_regs->push_back(dup_reg);
	} // i_reg loop.

	// Save in binarized format.
	char tag_resampled_geno_matbed_fp[1000];
	sprintf(tag_resampled_geno_matbed_fp, "%s_resampled_tags.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(resampled_tag_geno_var_regs, NULL, resampled_sample_ids, tag_resampled_geno_matbed_fp);
}

static void* thread_callback_tag_anchored_target_sampling(void* __thread_info_ptr_ptr)
{
	void** __thread_info_ptr = (void**)(__thread_info_ptr_ptr);

	////////////
	int* int_vals = (int*)(__thread_info_ptr[0]);
	int thread_i = int_vals[0];
	int n_threads = int_vals[1];
	int n_blocks = int_vals[2];
	int n_resampled_sample_size = int_vals[3];
	int n_original_haps = int_vals[4];

	////////////
	//double* dbl_vals = (double*)(__thread_info_ptr[1]);

	////////////
	vector<t_annot_region*>** per_tag2tag_target_blocks = (vector<t_annot_region*>**)(__thread_info_ptr[2]);
	vector<t_annot_region*>* per_tag2tag_target_block_starting_tag_reg = (vector<t_annot_region*>*)(__thread_info_ptr[3]);
	vector<t_annot_region*>* per_tag2tag_target_block_ending_tag_reg = (vector<t_annot_region*>*)(__thread_info_ptr[4]);

	////////////

	fprintf(stderr, "Thread %d: Setting target genotypes..           \r", thread_i);
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		if (i_s % n_threads != thread_i)
		{
			continue;
		}

		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			if (i_s % 100 == 0 && i_hap == 0)
			{
				fprintf(stderr, "Setting %d/%d [%d]        \r", i_s, n_resampled_sample_size, i_hap);
			}

			// Allow at most one recombination in the target region.
			// Process each block, this can be done in parallel.
			for (int i_bl = 0; i_bl < n_blocks; i_bl++)
			{
				vector<t_annot_region*>* cur_block_target_regs = per_tag2tag_target_blocks[i_bl];
				t_annot_region* cur_block_starting_tag_reg = per_tag2tag_target_block_starting_tag_reg->at(i_bl);
				t_annot_region* cur_block_ending_tag_reg = per_tag2tag_target_block_ending_tag_reg->at(i_bl);

				// Handle the start and end blocks at the same haplotype.
				if (cur_block_starting_tag_reg == NULL &&
					cur_block_ending_tag_reg == NULL)
				{
					fprintf(stderr, "Block %d: Both start and end tags are NULL.\n", i_bl);
					exit(1);
				}

				if (cur_block_starting_tag_reg == NULL)
				{
					cur_block_starting_tag_reg = cur_block_ending_tag_reg;
				}

				if (cur_block_ending_tag_reg == NULL)
				{
					cur_block_ending_tag_reg = cur_block_starting_tag_reg;
				}

				// Check the states at the start/end blocks.
				void** starting_tag_info = (void**)(cur_block_starting_tag_reg->data);
				void** ending_tag_info = (void**)(cur_block_ending_tag_reg->data);

				int** start_tag_per_hap_sampled_hap_i = (int**)(starting_tag_info[4]);
				int** end_tag_per_hap_sampled_hap_i = (int**)(ending_tag_info[4]);

				int* start_tag_cur_hap_sampled_hap_i = start_tag_per_hap_sampled_hap_i[i_hap];
				int* end_tag_cur_hap_sampled_hap_i = end_tag_per_hap_sampled_hap_i[i_hap];

				int cur_subj_start_tag_i_hap = start_tag_cur_hap_sampled_hap_i[i_s];
				int cur_subj_end_tag_i_hap = end_tag_cur_hap_sampled_hap_i[i_s];

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Check the tag genotypes for the starting and ending tags.
				// Get the resampled alleles for these tags.
				char* starting_tag_res_geno = (char*)(starting_tag_info[1]);
				char starting_var_hap_res_allele = get_allele_per_haplotype(starting_tag_res_geno[i_s], i_hap);

				char* ending_tag_res_geno = (char*)(ending_tag_info[1]);
				char ending_var_hap_res_allele = get_allele_per_haplotype(ending_tag_res_geno[i_s], i_hap);

				// Now, extract the original alleles for these tags.
				int starting_reg_i_s = start_tag_cur_hap_sampled_hap_i[i_s] / 2;
				int starting_reg_i_hap = start_tag_cur_hap_sampled_hap_i[i_s] % 2;
				char* starting_tag_orig_geno = (char*)(starting_tag_info[0]);
				char starting_var_hap_orig_allele = get_allele_per_haplotype(starting_tag_orig_geno[starting_reg_i_s], starting_reg_i_hap);

				int ending_reg_i_s = end_tag_cur_hap_sampled_hap_i[i_s] / 2;
				int ending_reg_i_hap = end_tag_cur_hap_sampled_hap_i[i_s] % 2;
				char* ending_tag_orig_geno = (char*)(ending_tag_info[0]);
				char ending_var_hap_orig_allele = get_allele_per_haplotype(ending_tag_orig_geno[ending_reg_i_s], ending_reg_i_hap);

				// Compare whether the res/orig alleles match.
				if (starting_var_hap_res_allele != starting_var_hap_orig_allele ||
					ending_var_hap_res_allele != ending_var_hap_orig_allele)
				{
					fprintf(stderr, "Sanity check failed: Start/End Tag resampled/original alleles could not be matched.\n");
					exit(1);
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				if (cur_subj_start_tag_i_hap == -1 ||
					cur_subj_end_tag_i_hap == -1)
				{
					fprintf(stderr, "Sanity check failed: Resampled index is uninitialized.\n");
					exit(1);
				}

				if (cur_subj_start_tag_i_hap >= n_original_haps ||
					cur_subj_end_tag_i_hap >= n_original_haps)
				{
					fprintf(stderr, "found illegal haplotype: %d, %d out of %d\n",
						cur_subj_start_tag_i_hap, cur_subj_end_tag_i_hap, n_original_haps);
					exit(1);
				}

				if (__DUMP_HAPRESAMPLING_MSGS__)
				{
					fprintf(stderr, "%d[%d]::Block %d: %d-%d [%d targets]: %d--%d ;; %.4f--%.4f\n",
						i_s, i_hap,
						i_bl, cur_block_starting_tag_reg->start, cur_block_ending_tag_reg->end,
						vecsize(cur_block_target_regs), cur_subj_start_tag_i_hap, cur_subj_end_tag_i_hap,
						cur_block_starting_tag_reg->dbl_score, cur_block_ending_tag_reg->dbl_score);
				}

				// If the haplotypes are the same, just copy alleles.
				if (cur_subj_start_tag_i_hap == cur_subj_end_tag_i_hap)
				{
					for (int i_target = 0; i_target < vecsize(cur_block_target_regs); i_target++)
					{
						int cur_target_resampled_i_hap = cur_subj_start_tag_i_hap;

						void** cur_target_info = (void**)(cur_block_target_regs->at(i_target)->data);
						char* cur_target_geno = (char*)(cur_target_info[0]);
						char* cur_target_resampled_geno = (char*)(cur_target_info[1]);

						int resampled_i_s = (cur_target_resampled_i_hap - (cur_target_resampled_i_hap % 2)) / 2;
						int resampled_i_hap = cur_target_resampled_i_hap % 2;
						int resampled_allele = get_allele_per_haplotype(cur_target_geno[resampled_i_s], resampled_i_hap);

						// Copy this allele.
						int prev_geno = cur_target_resampled_geno[i_s];
						int updated_geno = prev_geno | (resampled_allele << i_hap);
						cur_target_resampled_geno[i_s] = updated_geno;
					} // i_target loop.
				} // same haplotype check.
				else
				{
					// Tags are at different haplotypes.
					// Look for the half genetic distance.
					double cur_block_targets_half_genetic_cM = (cur_block_starting_tag_reg->dbl_score + cur_block_ending_tag_reg->dbl_score) / 2;

					for (int i_target = 0; i_target < vecsize(cur_block_target_regs); i_target++)
					{
						// Select the i_hap based on genetic distance.
						int cur_target_resampled_i_hap = cur_subj_start_tag_i_hap;
						if (cur_block_target_regs->at(i_target)->dbl_score > cur_block_targets_half_genetic_cM)
						{
							cur_target_resampled_i_hap = cur_subj_end_tag_i_hap;
						}

						if (cur_block_target_regs->at(i_target)->dbl_score > cur_block_ending_tag_reg->dbl_score ||
							cur_block_target_regs->at(i_target)->dbl_score < cur_block_starting_tag_reg->dbl_score)
						{
							fprintf(stderr, "Target is out of block's genetic range: %.4f\n", cur_block_target_regs->at(i_target)->dbl_score);
							exit(1);
						}

						if (cur_block_target_regs->at(i_target)->start > cur_block_ending_tag_reg->end ||
							cur_block_target_regs->at(i_target)->start < cur_block_starting_tag_reg->start)
						{
							fprintf(stderr, "Target is out of block's position range: %d\n", cur_block_target_regs->at(i_target)->start);
							exit(1);
						}

						if (__DUMP_HAPRESAMPLING_MSGS__)
						{
							fprintf(stderr, "%d. Target: %d [%.4f]: State: %d [Half genetic distance: %.4f]\n",
								i_target,
								cur_block_target_regs->at(i_target)->start,
								cur_block_target_regs->at(i_target)->dbl_score,
								cur_target_resampled_i_hap,
								cur_block_targets_half_genetic_cM);
						}

						// Extract the target variant information.
						void** cur_target_info = (void**)(cur_block_target_regs->at(i_target)->data);
						char* cur_target_geno = (char*)(cur_target_info[0]); // This is the actual genotypes for this variant.
						char* cur_target_resampled_geno = (char*)(cur_target_info[1]); // This is the resampled genotype.

						// Select the resampled individual.
						int resampled_i_s = (cur_target_resampled_i_hap - (cur_target_resampled_i_hap % 2)) / 2;
						int resampled_i_hap = cur_target_resampled_i_hap % 2;
						int resampled_allele = get_allele_per_haplotype(cur_target_geno[resampled_i_s], resampled_i_hap);

						// Copy this allele.
						int prev_geno = cur_target_resampled_geno[i_s];
						int updated_geno = prev_geno | (resampled_allele << i_hap);
						cur_target_resampled_geno[i_s] = updated_geno;
					} // i_target loop.
				}
			} // i_bl loop.
		} // i_hap loop.
	} // i_s loop.

	return(NULL);
} // thread_callback_tag_anchored_target_sampling function.

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
	char* op_prefix)
{
	fprintf(stderr, "%d-thread AF-prior-based Re-sampling %d recombination-aware genotypes using N_e/N_ref_hap=%.3f and allelic error=%.6f in genomic interval [%d-%d]\n",
		n_threads,
		n_resampled_sample_size,
		N_e_2_n_ref_haplotypes,
		allele_error_prob,
		start_pos, end_pos);

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

	double N_e = N_e_2_n_ref_haplotypes * n_original_haps;

	// This is used to speed-up resampling by setting a minimum distance between recombination "hotspots".
	//double min_cM_delta_per_recomb = 0.01;

	fprintf(stderr, "Sampling from %d haplotypes using N_e=%.3f, allele_eps=%.3f with min cM per recombable var=%.4f.\n", (int)(n_original_haps), N_e, allele_error_prob, min_cM_delta_per_recomb);

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	//int n_resampled_sample_size = (int)(upsampling_rate * sample_ids->size());
	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	vector<char**>** per_chr_per_subj_haplotypes = new vector<char**>*[restr_geno_regs->chr_ids->size()];

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
			//fprintf(stderr, "%d: %.4f\n", cur_chrom_var_regs->at(i_reg)->start, cur_chrom_var_regs->at(i_reg)->dbl_score);
		} // i_reg loop.

		// Assign variant clusterings based on cM distance, this speeds up-sampling quite a lot.
		double cur_cM = -1;
		int n_recombable_vars = 0;
		for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
		{
			// IF the current variant is further away than the minimum cM cutoff, set it as a recomb. hotspot.
			if (fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_cM) > min_cM_delta_per_recomb)
			{
				cur_cM = cur_chrom_var_regs->at(i_reg)->dbl_score;
				cur_chrom_var_regs->at(i_reg)->score = 1;
				n_recombable_vars++;
			}
			else
			{
				cur_chrom_var_regs->at(i_reg)->score = 0;
			}
		} // i_reg loop.

		fprintf(stderr, "%d recombable variants\n", n_recombable_vars);

		fprintf(stderr, "Extracting per chromosome per subject haplotypes.\n");
		vector<char**>* per_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(cur_chrom_var_regs, sample_ids);
		per_chr_per_subj_haplotypes[i_chr] = per_subj_haplotypes;

		// Set the resampled genotypes for the current variant.
		for (int i_reg = 0; i_reg < (int)cur_chrom_var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
			void** new_reg_info = new void* [10];
			new_reg_info[0] = cur_reg_info[0];

			char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
			memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
			new_reg_info[1] = cur_reg_resampled_geno;

			double* n_recombs_per_var = new double[n_resampled_sample_size + 2];
			memset(n_recombs_per_var, 0, sizeof(double) * (n_resampled_sample_size + 2));
			new_reg_info[2] = n_recombs_per_var;

			double* n_errors_per_var = new double[n_resampled_sample_size + 2];
			memset(n_errors_per_var, 0, sizeof(double) * (n_resampled_sample_size + 2));
			new_reg_info[3] = n_errors_per_var;

			int** per_hap_sampled_hap_i = new int* [2];
			per_hap_sampled_hap_i[0] = new int[n_resampled_sample_size];
			per_hap_sampled_hap_i[1] = new int[n_resampled_sample_size];

			// Initialize all to -1, this is for sanity checking later.
			for (int i_re_s = 0; i_re_s < n_resampled_sample_size; i_re_s++)
			{
				per_hap_sampled_hap_i[0][i_re_s] = -1;
				per_hap_sampled_hap_i[1][i_re_s] = -1;
			} // i_re_s loop.
			new_reg_info[4] = per_hap_sampled_hap_i;

			// Replace the info ptr.
			cur_chrom_var_regs->at(i_reg)->data = new_reg_info;
		} // i_reg loop.
	} // chromosome index.

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting %d. re-sampling thread..\n", thread_i);

		void** thread_ptrs_list = new void* [20];

		thread_ptrs_list[0] = restr_geno_regs;

		thread_ptrs_list[1] = recombination_rate_dir;

		int* thread_i_ptr = new int[2];
		thread_i_ptr[0] = thread_i;
		thread_ptrs_list[2] = thread_i_ptr;

		int* n_threads_ptr = new int[2];
		n_threads_ptr[0] = n_threads;
		thread_ptrs_list[3] = n_threads_ptr;

		double* N_e_ptr = new double[2];
		N_e_ptr[0] = N_e;
		thread_ptrs_list[4] = N_e_ptr;

		double* allele_error_ptr = new double[2];
		allele_error_ptr[0] = allele_error_prob;
		thread_ptrs_list[5] = allele_error_ptr;

		double* length_cutoff_in_bps_ptr = new double[2];
		length_cutoff_in_bps_ptr[0] = segment_length_cutoff_bp;
		thread_ptrs_list[6] = length_cutoff_in_bps_ptr;

		double* length_cutoff_in_cM_ptr = new double[2];
		length_cutoff_in_cM_ptr[0] = segment_length_cutoff_cM;
		thread_ptrs_list[7] = length_cutoff_in_cM_ptr;

		double* length_cutoff_in_var_number_ptr = new double[2];
		length_cutoff_in_var_number_ptr[0] = segment_length_cutoff_in_var_number;
		thread_ptrs_list[8] = length_cutoff_in_var_number_ptr;

		int* n_original_haps_ptr = new int[2];
		n_original_haps_ptr[0] = n_original_haps;
		thread_ptrs_list[9] = n_original_haps_ptr;

		int* n_resampled_sample_size_ptr = new int[2];
		n_resampled_sample_size_ptr[0] = n_resampled_sample_size;
		thread_ptrs_list[10] = n_resampled_sample_size_ptr;

		int* start_end_coords = new int[2];
		start_end_coords[0] = start_pos;
		start_end_coords[1] = end_pos;
		thread_ptrs_list[11] = start_end_coords;

		thread_ptrs_list[12] = haplotype_indices_2_sample;

		thread_ptrs_list[13] = per_chr_per_subj_haplotypes;

		t_ansi_thread* cur_thread = new t_ansi_thread(resampling_thread_callback_length_cutoff_save_recomb, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Started %d/%d threads; waiting.\n", (int)resampling_threads->size(), n_threads);

	for (int thread_i = 0; thread_i < (int)resampling_threads->size(); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		fprintf(stderr, "%d. thread finished.\n", thread_i);
	} // thread_i waiting loop.

	fprintf(stderr, "Saving re-sampled genotypes and information.\n");

	///////////////////////////////////////////////////////////////////////////////////////////
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
	char resampling_patterns_sigbed_fp[1000];
	sprintf(resampling_patterns_sigbed_fp, "%s_resampling_pattern_signal_var_regs.sigbed.gz", op_prefix);
	vector<t_annot_region*>* resamp_info_regs = new vector<t_annot_region*>();
	for (int i_var = 0; i_var < (int)haplocoded_geno_regs->size(); i_var++)
	{
		void** cur_var_info = (void**)(haplocoded_geno_regs->at(i_var)->data);
		t_annot_region* dup_reg = duplicate_region(haplocoded_geno_regs->at(i_var));
		void** dup_reg_info = new void* [5];
		dup_reg_info[0] = cur_var_info[4];
		dup_reg->data = dup_reg_info;
		resamp_info_regs->push_back(dup_reg);
	} // i_var loop.
	//save_resampling_pattern_signal_var_regs(resamp_info_regs, n_resampled_sample_size, resampling_patterns_sigbed_fp);
	save_resampling_pattern_signal_var_regs(resamp_info_regs, vecsize(sample_ids), n_resampled_sample_size, N_e_2_n_ref_haplotypes, resampling_patterns_sigbed_fp);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Save the # of recombinations at each variant.
	char per_var_n_recombs_fp[1000];
	sprintf(per_var_n_recombs_fp, "%s_per_var_n_recombs.txt", op_prefix);
	FILE* f_per_pos_n_recombs = open_f(per_var_n_recombs_fp, "w");
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		double* cur_var_n_recombs_per_sample = (double*)(cur_reg_data[2]);

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
	} // i_reg loop.
	fclose(f_per_pos_n_recombs);
	///////////////////////////////////////////////////////////////////////////////////////////
	// Write the number of errors for each variant.
	char per_var_n_errors_fp[1000];
	sprintf(per_var_n_errors_fp, "%s_per_var_n_errors.txt", op_prefix);
	FILE* f_per_pos_n_errors = open_f(per_var_n_errors_fp, "w");
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		double* cur_var_n_errors_per_sample = (double*)(cur_reg_data[3]);
		char* cur_var_geno = (char*)(cur_reg_data[0]);

		double cur_var_n_errors = 0;
		for (int i_s = 0; i_s < (int)resampled_sample_ids->size(); i_s++)
		{
			cur_var_n_errors += cur_var_n_errors_per_sample[i_s];
		} // i_s loop.

		double total_geno = 0;
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			total_geno += get_genotype_per_haplocoded_genotype(cur_var_geno[i_s]);
		}

		fprintf(f_per_pos_n_errors, "%s\t%d\t%d\t%s\t%.1f\t%.1f\t+\n",
			haplocoded_geno_regs->at(i_reg)->chrom,
			translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			haplocoded_geno_regs->at(i_reg)->name,
			cur_var_n_errors,
			total_geno);
	} // i_reg loop.
	fclose(f_per_pos_n_errors);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save resampled tag genotypes.
	vector<t_annot_region*>* resampled_tag_geno_var_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		t_annot_region* dup_reg = duplicate_region(haplocoded_geno_regs->at(i_reg));
		void** cur_reg_data = (void**)(haplocoded_geno_regs->at(i_reg)->data);

		void** dup_reg_info = new void* [2];
		dup_reg_info[0] = cur_reg_data[1];
		dup_reg->data = dup_reg_info;
		resampled_tag_geno_var_regs->push_back(dup_reg);
	} // i_reg loop.

	// Save in binarized format.
	char tag_resampled_geno_matbed_fp[1000];
	sprintf(tag_resampled_geno_matbed_fp, "%s_resampled_tags.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(resampled_tag_geno_var_regs, NULL, resampled_sample_ids, tag_resampled_geno_matbed_fp);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Sampling targets for each haplotype using anchored tags.\n");

	// Following code assume there is only 1 chromosome; this can be changed in the future since tag resampling can handle multiple chromosomes.
	if (restr_geno_regs->chr_ids->size() > 1)
	{
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");

		exit(1);
	}

	vector<t_annot_region*>* target_geno_regs = load_variant_signal_regions_wrapper(haplocoded_target_genotype_matrix_fp, sample_ids_list_fp);
	sort(target_geno_regs->begin(), target_geno_regs->end(), sort_regions);
	vector<char*>* target_chr_ids = get_chr_ids(target_geno_regs);
	if (vecsize(target_chr_ids) != 1)
	{
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");
		fprintf(stderr, "** ERROR: Currently we cannot sample multi-chromosomal genotype regions. Separate the data into different chromosomes and do sampling again. **\n");

		exit(1);
	}

	///////////////////////////
	// Validate tag and target chromosomes match.
	char* target_chrom = target_chr_ids->at(0);
	if (!t_string::compare_strings(restr_geno_regs->chr_ids->at(0), target_chrom))
	{
		fprintf(stderr, "Tags and targets are not on the same chromosome: %s vs %s\n",
			restr_geno_regs->chr_ids->at(0), target_geno_regs->at(0)->chrom);

		exit(1);
	}
	///////////////////////////

	// Assign resampled target genotypes.
	for (int i_reg = 0; i_reg < vecsize(target_geno_regs); i_reg++)
	{
		// Allocate the resampled genotypes for the target variants.
		char* resampled_geno_sig = new char[n_resampled_sample_size + 2];
		memset(resampled_geno_sig, 0, n_resampled_sample_size * sizeof(char));

		void** target_reg_info = (void**)(target_geno_regs->at(i_reg)->data);
		target_reg_info[1] = resampled_geno_sig;
	} // i_reg loop.
	///////////////////////////

	fprintf(stderr, "Assigning genetic distances to target variants.\n");
	char cur_chr_recombination_rate_fp[1000];
	sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, target_chrom);
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		exit(1);
	}

	// Assign the recomb rates.
	fprintf(stderr, "Setting recombination rates.\n");
	for (int i_reg = 0; i_reg < (int)target_geno_regs->size(); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(target_geno_regs->at(i_reg), cur_chrom_recomb_regs);
		target_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.
	///////////////////////////

	// Pool tag/target regions and build the target blocks
	vector<t_annot_region*>* tag_target_geno_regs = new vector<t_annot_region*>();

	// Mark tags and targets; strand distinguishes the tag/targets. We dont use score since it is used before for marking recombable tags.
	for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
	{
		haplocoded_geno_regs->at(i_reg)->strand = 0;
	} // i_Reg loop.
	for (int i_reg = 0; i_reg < (int)target_geno_regs->size(); i_reg++)
	{
		target_geno_regs->at(i_reg)->strand = 1;
	} // i_Reg loop.

	tag_target_geno_regs->insert(tag_target_geno_regs->end(), haplocoded_geno_regs->begin(), haplocoded_geno_regs->end());
	tag_target_geno_regs->insert(tag_target_geno_regs->end(), target_geno_regs->begin(), target_geno_regs->end());

	// Sort tag+target list before finding blocks.
	sort(tag_target_geno_regs->begin(), tag_target_geno_regs->end(), sort_regions);
	///////////////////////////

	// Find blocks of tag-tag blocks of target regions: We always have 1+n_tags blocks.
	fprintf(stderr, "Extracting tag-tag blocks of target variants.\n");
	vector<t_annot_region*>** per_tag2tag_target_blocks = new vector<t_annot_region*>*[vecsize(haplocoded_geno_regs) + 3];

	// Initialize tag2tag block.
	vector<t_annot_region*>* per_tag2tag_target_block_starting_tag_reg = new vector<t_annot_region*>();
	per_tag2tag_target_block_starting_tag_reg->push_back(NULL);
	vector<t_annot_region*>* per_tag2tag_target_block_ending_tag_reg = new vector<t_annot_region*>();

	int i_block = 0;
	vector<t_annot_region*>* cur_block_target_regs = new vector<t_annot_region*>();

	int i_reg = 0;
	while (i_reg < vecsize(tag_target_geno_regs))
	{
		// Empty the current block target regions.
		cur_block_target_regs->clear();

		// Find the next tag region for the current block.
		while (i_reg < vecsize(tag_target_geno_regs))
		{
			// Is this a tag; regardless of whether we have target regions or not, finding a new tag closes the block.
			if (tag_target_geno_regs->at(i_reg)->strand == 0)
			{
				// Ending region.
				per_tag2tag_target_block_ending_tag_reg->push_back(tag_target_geno_regs->at(i_reg));

				// Initiate the next block: this block starts with the closing variant.
				per_tag2tag_target_block_starting_tag_reg->push_back(tag_target_geno_regs->at(i_reg));

				if (__DUMP_HAPRESAMPLING_MSGS__)
				{
					fprintf(stderr, ">>>>Closing tag (%d) in %d. block: %s:%d (%s)\n",
						i_reg,
						i_block,
						tag_target_geno_regs->at(i_reg)->chrom,
						tag_target_geno_regs->at(i_reg)->start, tag_target_geno_regs->at(i_reg)->name);
				}

				// Move to the next variant.
				i_reg++;
				break;
			} // tag check.
			else if (tag_target_geno_regs->at(i_reg)->strand == 1)
			{
				// This is not a tag, add it to the targets.
				cur_block_target_regs->push_back(tag_target_geno_regs->at(i_reg));

				if (__DUMP_HAPRESAMPLING_MSGS__)
				{
					fprintf(stderr, "Target (%d) in %d. block: %s:%d (%s)\n",
						i_reg,
						i_block,
						tag_target_geno_regs->at(i_reg)->chrom,
						tag_target_geno_regs->at(i_reg)->start, tag_target_geno_regs->at(i_reg)->name);
				}
			} // target check.

			// Move to the next variant.
			i_reg++;
		} // i_reg loop.

		// For the currently finished block, set the target regions.
		vector<t_annot_region*>* block_target_regs = new vector<t_annot_region*>();
		block_target_regs->insert(block_target_regs->end(), cur_block_target_regs->begin(), cur_block_target_regs->end());
		per_tag2tag_target_blocks[i_block] = block_target_regs;
		i_block++;
	} // i_reg loop.

	// Add the final block.
	per_tag2tag_target_blocks[i_block] = cur_block_target_regs;
	per_tag2tag_target_block_ending_tag_reg->push_back(NULL); // This is the non-existing tag at the end of the chromosome.

	if (i_block != vecsize(per_tag2tag_target_block_ending_tag_reg) ||
		i_block != vecsize(per_tag2tag_target_block_ending_tag_reg))
	{
		fprintf(stderr, "Could not match the number of blocks to the number of start/end tag regions per block.\n");
		exit(1);
	}

	// Count the number of target variants in all blocks.
	int n_total_targets_in_blocks = 0;
	for (int i_bl = 0; i_bl < i_block; i_bl++)
	{
		n_total_targets_in_blocks += vecsize(per_tag2tag_target_blocks[i_bl]);
	} // i_bl loop.

	if (n_total_targets_in_blocks != vecsize(target_geno_regs))
	{
		fprintf(stderr, "# of targets do not match the number of variants in blocks: %d vs %d\n",
			n_total_targets_in_blocks, vecsize(target_geno_regs));

		exit(1);
	}

	int n_blocks = i_block;
	fprintf(stderr, "Setup %d blocks for %d target variants.\n", n_blocks, vecsize(target_geno_regs));

	///////////////////////////

	fprintf(stderr, "Starting target resampling threads..\n");
	vector<t_ansi_thread*>* target_sampling_threads = new vector<t_ansi_thread*>();
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		void** thread_info_ptr = new void* [20];

		////////////
		int* int_vals = new int[10];
		thread_info_ptr[0] = int_vals;
		int_vals[0] = i_thread;
		int_vals[1] = n_threads;
		int_vals[2] = n_blocks;
		int_vals[3] = n_resampled_sample_size;
		int_vals[4] = n_original_haps;

		////////////
		double* dbl_vals = new double[10];
		thread_info_ptr[1] = dbl_vals;

		////////////
		thread_info_ptr[2] = per_tag2tag_target_blocks;
		thread_info_ptr[3] = per_tag2tag_target_block_starting_tag_reg;
		thread_info_ptr[4] = per_tag2tag_target_block_ending_tag_reg;

		t_ansi_thread* new_thread = new t_ansi_thread(thread_callback_tag_anchored_target_sampling, thread_info_ptr);
		new_thread->run_thread();
		target_sampling_threads->push_back(new_thread);
		t_string::print_padded_string(stderr, '\r', 100, "Started %d. thread..", i_thread);
	} // i_thread loop.

	t_string::print_padded_string(stderr, '\n', 100, "Waiting for threads..");
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		target_sampling_threads->at(i_thread)->wait_thread();
	} // i_thread loop.

	///////////////////////////
	// Fixe the data index of resampled target genotypes.
	for (int i_reg = 0; i_reg < (int)target_geno_regs->size(); i_reg++)
	{
		void** cur_reg_data = (void**)(target_geno_regs->at(i_reg)->data);
		cur_reg_data[0] = cur_reg_data[1];
	} // i_reg loop.

	// Save targets in binarized format.
	char target_resampled_geno_matbed_fp[1000];
	sprintf(target_resampled_geno_matbed_fp, "%s_resampled_targets.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(target_geno_regs, NULL, resampled_sample_ids, target_resampled_geno_matbed_fp);
} // resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling function.

