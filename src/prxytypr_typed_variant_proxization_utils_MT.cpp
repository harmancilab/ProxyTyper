#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "prxytypr_file_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_variation_tools.h"
//#include "../../lib/genomics_utils/variation/imputation_utils.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_nomenclature.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_histogram.h"
//#include "../../lib/genomics_utils/signal_track/signal_track_tools.h"
//#include "../../lib/utils/xmath/matrix/matrix_linalg_utils.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_vector_macros.h"
#include "prxytypr_proxytyper.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

const bool __DUMP_PROXIZATION_MSGS__ = false;

void random_per_sample_hap_switch(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp, double switch_prob, char* op_matbed_fp)
{
	fprintf(stderr, "Random haplotype switching with %.4f probability.\n", switch_prob);

	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	sort(haplocoded_var_regs->begin(), haplocoded_var_regs->end(), sort_regions);

	fprintf(stderr, "Loaded %d variants for %d subjects.\n", vecsize(haplocoded_var_regs), vecsize(sample_ids));

	for (int i_var = 0; i_var < vecsize(haplocoded_var_regs); i_var++)
	{
		void** var_reg_info = (void**)(haplocoded_var_regs->at(i_var)->data);

		void** new_var_reg_info = new void* [4];

		new_var_reg_info[0] = var_reg_info[0];

		char* new_geno_sig = new char[sample_ids->size() + 2];
		new_var_reg_info[1] = new_geno_sig;

		haplocoded_var_regs->at(i_var)->data = new_var_reg_info;
	} // i_var loop.

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Variants do not look like they are phased..\n");
		exit(1);
	}

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	for (int i_s = 0; i_s < vecsize(sample_ids); i_s++)
	{
		int n_switches = 0;
		int cur_allele_haps[2];
		cur_allele_haps[0] = 0;
		cur_allele_haps[1] = 1;
		for (int i_var = 0; i_var < vecsize(haplocoded_var_regs); i_var++)
		{
			if (rng->random_double_ran3() < switch_prob)
			{
				// Make a switch.
				cur_allele_haps[0] = 1 - cur_allele_haps[0];
				cur_allele_haps[1] = 1 - cur_allele_haps[1];

				n_switches++;
			}

			void** var_info = (void**)(haplocoded_var_regs->at(i_var)->data);
			char* orig_geno_sig = (char*)(var_info[0]);
			char* new_geno_sig = (char*)(var_info[1]);

			// Copy alleles.
			int alleles[2];
			alleles[0] = get_allele_per_haplotype(orig_geno_sig[i_s], 0);
			alleles[1] = get_allele_per_haplotype(orig_geno_sig[i_s], 1);
			
			int new_geno = alleles[0] * (1 << cur_allele_haps[0]) + alleles[1] * (1 << cur_allele_haps[1]);
			new_geno_sig[i_s] = (char)(new_geno);
		} // i_var loop.

		fprintf(stderr, "%s: %d/%d switches.\n", sample_ids->at(i_s), n_switches, vecsize(haplocoded_var_regs));
	} // i_s loop.

	// Replace the genotype signals.
	for (int i_var = 0; i_var < vecsize(haplocoded_var_regs); i_var++)
	{
		void** var_reg_info = (void**)(haplocoded_var_regs->at(i_var)->data);
		var_reg_info[0] = var_reg_info[1];
	} // i_var loop.

	// Save.
	binarize_variant_genotype_signal_regions(haplocoded_var_regs, NULL, sample_ids, op_matbed_fp);
}

static void* callback_filter_proxize_variants_per_vicinity(void* _thread_info_ptr)
{
	void** thread_info = (void**)(_thread_info_ptr);

	int* which_outof = (int*)(thread_info[0]);
	int which = which_outof[0];
	int outof = which_outof[1];
	
	vector<char*>* sample_ids = (vector<char*>*)(thread_info[1]);
	
	vector<t_annot_region*>* cur_chr_haplocoded_var_regs = (vector<t_annot_region*>*)(thread_info[2]);
	
	t_rng* rng = (t_rng*)(thread_info[3]);

	int* int_pars_list = (int*)(thread_info[4]);
	int avg_geno_mod = int_pars_list[0];
	int n_vicinity = int_pars_list[1];

	double* dbl_pars_list = (double*)(thread_info[5]);
	double allele_err_eps = dbl_pars_list[0];
	
	// Process each sample.
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		// Process each haplotype
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			for (int i_var = 0; i_var < (int)cur_chr_haplocoded_var_regs->size(); i_var++)
			{
				if (i_var % outof != which)
				{
					continue;
				}

				// Store the var_i info.
				void** var_i_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(i_var)->data);
				//char* var_i_geno_sig = var_i_reg_info[0];
				char* var_i_proxized_geno_sig = (char*)(var_i_reg_info[1]);

				// Select the proxization parameters for this variant.
				int* per_var_weight = (int*)(var_i_reg_info[2]);
				int** per_var2var_interaction_weight = (int**)(var_i_reg_info[3]);
				int*** per_var2var2var_interaction_weight = (int***)(var_i_reg_info[4]);

				if (per_var_weight == NULL ||
					per_var2var_interaction_weight == NULL)
				{
					fprintf(stderr, "Sanity check failed: %s:%d-%d does not have a vicinity model.\n",
						cur_chr_haplocoded_var_regs->at(i_var)->chrom,
						cur_chr_haplocoded_var_regs->at(i_var)->start,
						cur_chr_haplocoded_var_regs->at(i_var)->end);

					exit(1);
				}

				// Set window bounds.
				int cur_win_start = MAX(0, (i_var - n_vicinity));
				int cur_win_end = MIN((int)cur_chr_haplocoded_var_regs->size() - 1, (i_var + n_vicinity));

				if (__DUMP_PROXIZATION_MSGS__)
				{
					// Write only for 1st subject.
					if (i_s == 0)
					{
						fprintf(stderr, "Variant @%d (%s:%d-%d::%s): Using window %d-%d\n",
							i_var,
							cur_chr_haplocoded_var_regs->at(i_var)->chrom,
							cur_chr_haplocoded_var_regs->at(i_var)->start,
							cur_chr_haplocoded_var_regs->at(i_var)->end,
							cur_chr_haplocoded_var_regs->at(i_var)->name,
							cur_win_start, cur_win_end);
					} // first subject check.
				} // msg check.

				// Calculate average.
				int cur_var_i_vicinity_allele_sum = 0;
				for (int j_var = cur_win_start; j_var <= cur_win_end; j_var++)
				{
					void** var_j_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(j_var)->data);
					char* var_j_geno_sig = (char*)(var_j_reg_info[0]);

					// Extract the allele on the current haplotype.
					int cur_j_allele = get_allele_per_haplotype(var_j_geno_sig[i_s], i_hap);

					int cur_int_term = per_var_weight[j_var - cur_win_start];
					if (cur_int_term == 1)
					{
						cur_var_i_vicinity_allele_sum += (cur_j_allele);
					}
					else if (cur_int_term == -1)
					{
						cur_var_i_vicinity_allele_sum += (1 - cur_j_allele);
					}

					// Check allele validity.
					if (cur_j_allele > 1 || cur_j_allele < 0)
					{
						fprintf(stderr, "Sanity check failed @ %s allele is not valid: %d\n",
							cur_chr_haplocoded_var_regs->at(j_var)->name,
							cur_j_allele);
						exit(1);
					}

					// Add the interaction terms of this allele: Make sure they cover only one half.
					for (int k_var = cur_win_start; k_var < j_var; k_var++)
					{
						void** var_k_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(k_var)->data);
						char* var_k_geno_sig = (char*)(var_k_reg_info[0]);

						int cur_k_allele = get_allele_per_haplotype(var_k_geno_sig[i_s], i_hap);

						int cur_int_term = per_var2var_interaction_weight[j_var - cur_win_start][k_var - cur_win_start];

						if (cur_int_term == 1)
						{
							cur_var_i_vicinity_allele_sum += (cur_j_allele * cur_k_allele);
						}
						else if (cur_int_term == -1)
						{
							cur_var_i_vicinity_allele_sum += (1 - (cur_j_allele * cur_k_allele));
						}

						// Add the final interaction.
						for (int l_var = cur_win_start; l_var < k_var; l_var++)
						{
							void** var_l_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(l_var)->data);
							char* var_l_geno_sig = (char*)(var_l_reg_info[0]);

							int cur_l_allele = get_allele_per_haplotype(var_l_geno_sig[i_s], i_hap);

							int cur_int_term = per_var2var2var_interaction_weight[j_var - cur_win_start][k_var - cur_win_start][l_var - cur_win_start];

							if (cur_int_term == 1)
							{
								cur_var_i_vicinity_allele_sum += (cur_j_allele * cur_k_allele * cur_l_allele);
							}
							else if (cur_int_term == -1)
							{
								cur_var_i_vicinity_allele_sum += (1 - (cur_j_allele * cur_k_allele * cur_l_allele));
							}

						} // l_var loop.
					} // k_var loop.
				} // j_var loop.

				// Add error:
				int allele_err = 0;
				if (rng->random_double_ran3() < allele_err_eps)
				{
					cur_chr_haplocoded_var_regs->at(i_var)->score++;
					allele_err = 1;
				}

				// Add error.
				cur_var_i_vicinity_allele_sum += allele_err;

				// Take the modulus.
				int cur_vic_mod_allele = cur_var_i_vicinity_allele_sum % avg_geno_mod;

				// Set the allele in the genotype.
				var_i_proxized_geno_sig[i_s] = var_i_proxized_geno_sig[i_s] | (cur_vic_mod_allele << i_hap);
			} // i_var loop.
		} // i_hap loop.
	} // i_s loop.

	return NULL;
} // callback_filter_proxize_variants_per_vicinity function.

// This function tests the proxization of variants (genotypes) using modular average of surrounding variants with interacting terms, i.e., AND operations.
void MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
	int n_threads,
	char* proxized_haplocoded_var_regs_op_fp)
{
	fprintf(stderr, "Loading per-variant vicinity parameters.\n");
	int loaded_n_vicinity = -1;
	int loaded_avg_geno_mod = -1;
	vector<t_annot_region*>* per_variant_vicinity_params = load_per_variant_site_mixing_parameters(per_var_vicinity_weight_params_fp, loaded_n_vicinity, loaded_avg_geno_mod);

	fprintf(stderr, "Loaded parameters for %d-vicinity and mod-%d for %d variants.\n", loaded_n_vicinity, loaded_avg_geno_mod, (int)per_variant_vicinity_params->size());

	int n_vicinity = loaded_n_vicinity;
	int avg_geno_mod = loaded_avg_geno_mod;

	fprintf(stderr, "Proxizing %s (%s) per variant-specific non-linear modular combinations using:\n\\sum_{vic @ [-%d, +%d]) (w_i x g_i)+\\sum_{vic @ [-%d, +%d] w_ij x g_i x g_j mod %d using %d variant parameters.\n",
		haplocoded_geno_signal_fp, sample_ids_list_fp,
		n_vicinity, n_vicinity, n_vicinity, n_vicinity, avg_geno_mod, (int)per_variant_vicinity_params->size());

	// Instantiate the RNG.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(0);
	}

	// Divide the variants into chromosomes.
	t_restr_annot_region_list* restr_haplocoded_var_regs = restructure_annot_regions(haplocoded_var_regs);

	fprintf(stderr, "Allocating and setting up the proxized variant genotypes.\n");
	for (int i_var = 0; i_var < (int)haplocoded_var_regs->size(); i_var++)
	{
		void** cur_var_reg_info = (void**)(haplocoded_var_regs->at(i_var)->data);
		char* geno_sig = (char*)(cur_var_reg_info[0]);
		char* proxized_geno_sig = new char[(int)sample_ids->size() + 2];

		// Reset all proxized genotypes to 0.
		memset(proxized_geno_sig, 0, sizeof(char) * (int)sample_ids->size());

		void** new_reg_info = new void* [10];
		new_reg_info[0] = geno_sig;
		new_reg_info[1] = proxized_geno_sig;
		new_reg_info[2] = NULL; // First degree vicinity weights.
		new_reg_info[3] = NULL; // Second degree vicinity weights.
		new_reg_info[4] = NULL; // Third degree vicinity weights.

		// Replace the variant information.
		haplocoded_var_regs->at(i_var)->data = new_reg_info;
		haplocoded_var_regs->at(i_var)->score = 0;
	} // i_var loop.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Assigning per-variant vicinity parameters.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(haplocoded_var_regs, per_variant_vicinity_params, true);
	fprintf(stderr, "Processing %d intersects.\n", (int)intersects->size());
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* var_reg = int_info->src_reg;
		t_annot_region* vic_params_reg = int_info->dest_reg;

		void** var_reg_info = (void**)(var_reg->data);
		void** vic_params_reg_info = (void**)(vic_params_reg->data);

		// Assign the vicinity parameters.
		var_reg_info[2] = vic_params_reg_info[2];
		var_reg_info[3] = vic_params_reg_info[3];
		var_reg_info[4] = vic_params_reg_info[4];
	} // i_int loop.

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start processing the proxized genotypes.
	for (int i_chr = 0; i_chr < (int)restr_haplocoded_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_haplocoded_var_regs = restr_haplocoded_var_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Sorting the variants per position..\n");

		sort(cur_chr_haplocoded_var_regs->begin(), cur_chr_haplocoded_var_regs->end(), sort_regions);

		// ADD THREADING HERE.
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int i_th = 0; i_th < n_threads; i_th++)
		{
			void** cur_thread_info_ptr = new void* [10];
			int* which_outof = new int[2];
			which_outof[0] = i_th;
			which_outof[1] = n_threads;

			cur_thread_info_ptr[0] = which_outof;
			cur_thread_info_ptr[1] = sample_ids;
			cur_thread_info_ptr[2] = cur_chr_haplocoded_var_regs;
			cur_thread_info_ptr[3] = new t_rng(t_seed_manager::seed_me_getrandom());

			int* int_pars_list = new int[10];
			int_pars_list[0] = avg_geno_mod;
			int_pars_list[1] = n_vicinity;
			cur_thread_info_ptr[4] = int_pars_list;

			double* dbl_pars_list = new double[10];
			dbl_pars_list[0] = allele_err_eps;
			cur_thread_info_ptr[5] = dbl_pars_list;

			//vector<char*>* sample_ids = (vector<char*>*)(thread_info[0]);
			//vector<t_annot_region*>* cur_chr_haplocoded_var_regs = (vector<char*>*)(thread_info[1]);
			//t_rng* rng = (vector<char*>*)(thread_info[2]);

			t_ansi_thread* cur_thread = new t_ansi_thread(callback_filter_proxize_variants_per_vicinity, cur_thread_info_ptr);

			fprintf(stderr, "Started thread %d..          \r", i_th);
			cur_thread->run_thread();

			threads->push_back(cur_thread);
		} // i_th loop.

		fprintf(stderr, "\nWaiting for threads..\n");
		for (int i_th = 0; i_th < n_threads; i_th++)
		{
			threads->at(i_th)->wait_thread();
			fprintf(stderr, "Thread %d finished..          \r", i_th);
		} // i_th loop.

		fprintf(stderr, "\nAll Threads finished..\n");
	} // i_chr loop.

	// Save.
	for (int i_var = 0; i_var < (int)haplocoded_var_regs->size(); i_var++)
	{
		void** cur_var_reg_info = (void**)(haplocoded_var_regs->at(i_var)->data);
		//char* geno_sig = (char*)(cur_var_reg_info[0]);
		char* proxized_geno_sig = (char*)(cur_var_reg_info[1]);

		cur_var_reg_info[0] = proxized_geno_sig;
	} // i_var loop.

	// Write the error statistics.
	char var_stats_fp[1000];
	sprintf(var_stats_fp, "%s_var_stats.txt", proxized_haplocoded_var_regs_op_fp);
	FILE* f_var_stats = open_f(var_stats_fp, "w");
	for (int i_var = 0; i_var < (int)haplocoded_var_regs->size(); i_var++)
	{
		fprintf(f_var_stats, "%s\t%d\t%d\t%s\t%d\n", haplocoded_var_regs->at(i_var)->chrom,
			haplocoded_var_regs->at(i_var)->start, haplocoded_var_regs->at(i_var)->end,
			haplocoded_var_regs->at(i_var)->name,
			haplocoded_var_regs->at(i_var)->score);
	} // i_var loop.
	close_f(f_var_stats, NULL);

	// Save the genotype signals.
	binarize_variant_genotype_signal_regions(haplocoded_var_regs, NULL, sample_ids, proxized_haplocoded_var_regs_op_fp);
} // proxize_variants_per_vicinity_modular_average function.

