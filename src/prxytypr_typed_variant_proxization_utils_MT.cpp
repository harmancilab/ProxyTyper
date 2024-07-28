#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <zlib.h>
#include <iostream>
#include "prxytypr_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_variation_tools.h"
#include "prxytypr_imputation_utils.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_nomenclature.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_histogram.h"
#include "prxytypr_matrix_linalg_utils.h"
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

//				if (__DUMP_PROXIZATION_MSGS__)
				{
					// Write only for 1st subject.
					if (i_s == 0)
					{
						//fprintf(stderr, "Variant @%d (%s:%d-%d::%s): Using window %d-%d\n",
						//	i_var,
						//	cur_chr_haplocoded_var_regs->at(i_var)->chrom,
						//	cur_chr_haplocoded_var_regs->at(i_var)->start,
						//	cur_chr_haplocoded_var_regs->at(i_var)->end,
						//	cur_chr_haplocoded_var_regs->at(i_var)->name,
						//	cur_win_start, cur_win_end);

						if (__DUMP_PROXIZATION_MSGS__)
						{
							fprintf(stderr, "==================================\nVar-%d:\n", i_var);
							for (int z_var = cur_win_start; z_var <= cur_win_end; z_var++)
							{
								fprintf(stderr, "%d\t%d\t%d-%d\t%s\n",
									cur_chr_haplocoded_var_regs->at(i_var)->start,
									cur_chr_haplocoded_var_regs->at(i_var)->end,
									cur_chr_haplocoded_var_regs->at(z_var)->start,
									cur_chr_haplocoded_var_regs->at(z_var)->end,
									cur_chr_haplocoded_var_regs->at(z_var)->name);
							}
						}
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

					//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(i_var)->start == 63244)
					if (__DUMP_PROXIZATION_MSGS__)
					{
						fprintf(stderr, "@ %d: j_var: %d [weight: %d; j_all: %d]; cumulative: %d\n",
							cur_chr_haplocoded_var_regs->at(i_var)->start,
							cur_chr_haplocoded_var_regs->at(j_var)->start,
							cur_int_term, cur_j_allele, 
							cur_var_i_vicinity_allele_sum);
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

						//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(i_var)->start == 63244)
						if (__DUMP_PROXIZATION_MSGS__)
						{
							fprintf(stderr, "@ %d: {j-k}: {%d-%d} [j-k weight: %d]; cumulative: %d\n",
								cur_chr_haplocoded_var_regs->at(i_var)->start,
								cur_chr_haplocoded_var_regs->at(j_var)->start,
								cur_chr_haplocoded_var_regs->at(k_var)->start,
								cur_int_term, cur_var_i_vicinity_allele_sum);
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

							//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(i_var)->start == 63244)
							if (__DUMP_PROXIZATION_MSGS__)
							{
								fprintf(stderr, "@ %d: {j-k-l}: {%d-%d-%d} [j-k-l weight: %d]; cumulative: %d\n",
									cur_chr_haplocoded_var_regs->at(i_var)->start,
									cur_chr_haplocoded_var_regs->at(j_var)->start,
									cur_chr_haplocoded_var_regs->at(k_var)->start,
									cur_chr_haplocoded_var_regs->at(l_var)->start,
									cur_int_term, cur_var_i_vicinity_allele_sum);
							}

						} // l_var loop.
					} // k_var loop.
				} // j_var loop.

				// Add error:
				int allele_err = 0;
				if (rng->random_double_ran3() < allele_err_eps)
				{
					cur_chr_haplocoded_var_regs->at(i_var)->dbl_score++;
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
// UPDATE 6/24: We added a new parameter to focus on regions of interest to be proxized; this is necessary when we are doing encoding only on, for example, typed variants. Rather than separating from the panel, this function can perform this operation in place.
// UPDATE 6/24: This function does not anonymize the coordinates/genetic maps, this should be handled in other mechanisms.
// UPDATE 6/24: We do testing on matrix representation to optimize performance. Use unstructured data if possible to improve performance.
// UPDATE 6/24: Output can be a prefix to save as matrix representation, which should make it faster to process.
void MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
	int n_threads,
	char* proxized_haplocoded_var_regs_op_fp)
{
	if (!check_file(per_var_vicinity_weight_params_fp))
	{
		fprintf(stderr, "Model weights are not found @ \"%s\"\n", per_var_vicinity_weight_params_fp);
		exit(1);
	}

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

	//// Instantiate the RNG.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* all_haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	// Validate the phasing of variants.
	int max_geno = get_max_genotype_value(all_haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector<t_annot_region*>* focus_haplocoded_regions = new vector<t_annot_region*>();

	size_t proxized_mem_pool_size = vecsize_t(per_variant_vicinity_params) * vecsize_t(sample_ids) * sizeof(unsigned char);
	fprintf(stderr, "Allocating the proxized genotype memory pool of %.3f Gb..\n", (double)proxized_mem_pool_size / (1024.0 * 1024.0 * 1024.0));
	unsigned char* proxized_geno_sig_mem_pool = new unsigned char[proxized_mem_pool_size];
	memset(proxized_geno_sig_mem_pool, 0, proxized_mem_pool_size);

	size_t n_ptrs_per_var_info = 10;
	void** focus_reg_info_mem_pool = new void* [n_ptrs_per_var_info * vecsize_t(per_variant_vicinity_params)];
	memset(focus_reg_info_mem_pool, 0, n_ptrs_per_var_info * vecsize_t(per_variant_vicinity_params) * sizeof(void*));

	fprintf(stderr, "Assigning per-variant vicinity parameters.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(per_variant_vicinity_params, all_haplocoded_var_regs, true);
	fprintf(stderr, "Processing %d intersects.\n", (int)intersects->size());
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* vic_params_reg = int_info->src_reg;
		t_annot_region* var_reg = int_info->dest_reg;

		if (vic_params_reg->start == var_reg->start &&
			vic_params_reg->end == var_reg->end)
		{
			void** vic_params_reg_info = (void**)(vic_params_reg->data);

			void** var_reg_info = (void**)(var_reg->data);

			// We are doing pointer arithmetic here, be careful with void** element.
			size_t var_i_pool = focus_haplocoded_regions->size();
			char* proxized_geno_sig = (char*)(proxized_geno_sig_mem_pool + var_i_pool * vecsize_t(sample_ids));
			void** new_reg_info = (void**)(focus_reg_info_mem_pool + var_i_pool * n_ptrs_per_var_info);

			new_reg_info[0] = var_reg_info[0];
			new_reg_info[1] = proxized_geno_sig;
			new_reg_info[2] = vic_params_reg_info[0]; // First degree vicinity weights.
			new_reg_info[3] = vic_params_reg_info[1]; // Second degree vicinity weights.
			new_reg_info[4] = vic_params_reg_info[2]; // Third degree vicinity weights.

			var_reg->data = new_reg_info;
			var_reg->dbl_score = 0; // # of errors.

			focus_haplocoded_regions->push_back(var_reg);
		}
	} // i_int loop.
	fprintf(stderr, "Proxizing %d/%d variant regions.\n", vecsize(focus_haplocoded_regions), vecsize(all_haplocoded_var_regs));

	// Divide the variants into chromosomes.
	t_restr_annot_region_list* restr_focus_haplocoded_var_regs = restructure_annot_regions(focus_haplocoded_regions);

	// Start processing the proxized genotypes.
	for (int i_chr = 0; i_chr < vecsize(restr_focus_haplocoded_var_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_focus_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_focus_haplocoded_var_regs = restr_focus_haplocoded_var_regs->regions_per_chrom[i_chr];

		// This is very important -- if variants are not sorted, algorithm does not perform proxizing on neighboring variants.
		fprintf(stderr, "Sorting the variants per position..\n");
		sort(cur_chr_focus_haplocoded_var_regs->begin(), cur_chr_focus_haplocoded_var_regs->end(), sort_regions);

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
			cur_thread_info_ptr[2] = cur_chr_focus_haplocoded_var_regs;
			//cur_thread_info_ptr[3] = new t_rng(t_seed_manager::seed_me_getrandom());
			cur_thread_info_ptr[3] = new t_rng(t_seed_manager::seed_me() + i_th);

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

			t_string::print_padded_string(stderr, '\r', 100, "Started thread %d..", i_th);

			cur_thread->run_thread();

			threads->push_back(cur_thread);
		} // i_th loop.

		fprintf(stderr, "\nWaiting for threads..\n");
		for (int i_th = 0; i_th < n_threads; i_th++)
		{
			threads->at(i_th)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Thread %d finished..", i_th);
		} // i_th loop.

		fprintf(stderr, "\nAll Threads finished..\n");
	} // i_chr loop.

	// Save: Save all variants.
	for (int i_var = 0; i_var < vecsize(focus_haplocoded_regions); i_var++)
	{
		void** cur_var_reg_info = (void**)(focus_haplocoded_regions->at(i_var)->data);

		// Replace the genotype signal with proxized genotypes; this is done only on the proxized variants.
		//char* geno_sig = (char*)(cur_var_reg_info[0]);
		char* proxized_geno_sig = (char*)(cur_var_reg_info[1]);
		cur_var_reg_info[0] = proxized_geno_sig;
	} // i_var loop.

	// Write the error statistics.
	char var_stats_fp[1000];
	sprintf(var_stats_fp, "%s_var_stats.txt", proxized_haplocoded_var_regs_op_fp);
	FILE* f_var_stats = open_f(var_stats_fp, "w");
	for (int i_var = 0; i_var < vecsize(focus_haplocoded_regions); i_var++)
	{
		fprintf(f_var_stats, "%s\t%d\t%d\t%s\t%d\n", focus_haplocoded_regions->at(i_var)->chrom,
			focus_haplocoded_regions->at(i_var)->start, focus_haplocoded_regions->at(i_var)->end,
			focus_haplocoded_regions->at(i_var)->name,
			(int)(focus_haplocoded_regions->at(i_var)->dbl_score));
	} // i_var loop.
	close_f(f_var_stats, NULL);

	// Save the genotype signals: Make sure to save all of the haplocoded variant regions.
	//binarize_variant_genotype_signal_regions(all_haplocoded_var_regs, NULL, sample_ids, proxized_haplocoded_var_regs_op_fp);
	//binarize_variant_signal_regions_wrapper(all_haplocoded_var_regs, sample_ids, proxized_haplocoded_var_regs_op_fp);
	binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(all_haplocoded_var_regs, sample_ids, true, n_threads, proxized_haplocoded_var_regs_op_fp);
} // proxize_variants_per_vicinity_modular_average function.

#define ERRMSG(msg) \
    do { \
        std::cout << "File: " << __FILE__ << ", Line: " << __LINE__ << " - " << msg << std::endl; \
    } while (0)

static void* callback_filter_proxize_variants_per_vicinity_no_buffer(void* _thread_info_ptr)
{
	void** thread_info = (void**)(_thread_info_ptr);

	int* int_vals = (int*)(thread_info[0]);
	int thread_i = int_vals[0];
	int n_threads = int_vals[1];
	int cur_thread_var_start_i = int_vals[2];
	int cur_thread_var_end_i = int_vals[3];
	int avg_geno_mod = int_vals[4];
	int n_vicinity = int_vals[5];

	double* dbl_vals = (double*)(thread_info[1]);
	double allele_err_eps = dbl_vals[0];

	vector<t_annot_region*>* cur_chr_haplocoded_var_regs = (vector<t_annot_region*>*)(thread_info[2]);

	// This is the list of typed variants that are used for hashing, it is necessary to keep this separate since we want to use only tags for hashing.
	vector<t_annot_region*>* cur_chr_vic_params_var_regs = (vector<t_annot_region*>*)(thread_info[3]);

	for (int i_reg = 1; i_reg < vecsize(cur_chr_vic_params_var_regs); i_reg++)
	{
		if (cur_chr_vic_params_var_regs->at(i_reg)->start <= cur_chr_vic_params_var_regs->at(i_reg - 1)->start)
		{
			fprintf(stderr, "SANITY CHECK FAILED!!!\n");
			exit(1);
		}
	}

	vector<char*>* sample_ids = (vector<char*>*)(thread_info[4]);

	t_rng* rng = (t_rng*)(thread_info[5]);

	// Instantiate an rng: Dont instantiate rng like this.
	//t_rng* rng = new t_rng(t_seed_manager::seed_me() + thread_i);

	//if (__DUMP_PROXIZATION_MSGS__)
	{
		t_string::print_padded_string(stderr, '\r', 100, "Thread-%d/%d: [%d-%d] / %d (%d subjects); %d-modular; %d-vicinity",
			thread_i, n_threads,
			cur_thread_var_start_i, cur_thread_var_end_i, vecsize(cur_chr_haplocoded_var_regs),
			vecsize(sample_ids),
			avg_geno_mod, n_vicinity);
	}

	// Open the file for writing.
	char cur_thread_proxized_geno_mat_fp[1000];
	sprintf(cur_thread_proxized_geno_mat_fp, "proxized_geno_sig_%d.matrix.gz", thread_i);
	gzFile f_cur_thread_proxized_geno_mat = gzopen(cur_thread_proxized_geno_mat_fp, "wb");
	if (f_cur_thread_proxized_geno_mat == NULL)
	{
		fprintf(stderr, "Unable to open %s for writing\n", cur_thread_proxized_geno_mat_fp);
		exit(1);
	}

	// Process each variant.
	char* cur_subj_var_i_proxized_geno_sig = new char[vecsize_t(sample_ids)];
	for (int all_i_var = cur_thread_var_start_i; all_i_var < cur_thread_var_end_i; all_i_var++)
	{
		//fprintf(stderr, "Thread-%d: var_i: %d\n", thread_i, i_var);
		/////////////////////
		// Store the var_i info.
		void** var_i_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(all_i_var)->data);
		char* var_geno_sig = (char* )(var_i_reg_info[0]);

		// Select the proxization parameters for this variant.
		int* per_var_weight = (int*)(var_i_reg_info[2]);
		int** per_var2var_interaction_weight = (int**)(var_i_reg_info[3]);
		int*** per_var2var2var_interaction_weight = (int***)(var_i_reg_info[4]);
		int vicin_param_reg_i = cur_chr_haplocoded_var_regs->at(all_i_var)->score;

		bool has_var_weight = (per_var_weight != NULL);
		bool has_var2var_weight = (per_var2var_interaction_weight != NULL);
		bool has_var2var2var_weight = (per_var2var2var_interaction_weight != NULL);

		if ((has_var_weight && has_var2var_weight && has_var2var2var_weight) ||
			(!has_var_weight && !has_var2var_weight && !has_var2var2var_weight)) 
		{
			if (__DUMP_PROXIZATION_MSGS__)
			{
				fprintf(stderr, "We are at a tag variant..");
			}

			if (has_var_weight && 
				vicin_param_reg_i == -1)
			{
				fprintf(stderr, "Sanity check failed @ %s (%d) : Variant has parameters but no associated vicinity model region..\n", __FILE__, __LINE__);
				exit(1);
			}
		}
		else 
		{
			fprintf(stderr, "Sanity check failed @ %s(%d): One of the weights exist, memory corruption??\n", __FILE__, __LINE__);
			exit(1);
		}		

		// If weights don't exist (untyped variant), copy the genotype signal and move on.
		if (!has_var_weight)
		{
			memcpy(cur_subj_var_i_proxized_geno_sig, var_geno_sig, vecsize(sample_ids) * sizeof(char));
		}
		else
		{
			// This is the index we will loop around on the tag variants.
			if (cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->start != cur_chr_haplocoded_var_regs->at(all_i_var)->start ||
				((void**)cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->data)[2] != per_var_weight ||
				((void**)cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->data)[3] != per_var2var_interaction_weight ||
				((void**)cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->data)[4] != per_var2var2var_interaction_weight)
			{
				ERRMSG("Tag variant weights do not match the matching untyped variant..");
				exit(1);
			}

			// Set window bounds.
			int cur_win_start = MAX(0, (vicin_param_reg_i - n_vicinity));
			int cur_win_end = MIN((int)cur_chr_vic_params_var_regs->size() - 1, (vicin_param_reg_i + n_vicinity));

			if (__DUMP_PROXIZATION_MSGS__)
			{
				fprintf(stderr, "Var-%d: Vicinity: [%d-%d]\n", all_i_var, cur_win_start, cur_win_end);
			}

			if (__DUMP_PROXIZATION_MSGS__)
			{
				fprintf(stderr, "==================================\nVar-%d:\n", all_i_var);
				for (int z_var = cur_win_start; z_var <= cur_win_end; z_var++)
				{
					fprintf(stderr, "%d\t%d\t%d-%d\t%s\n",
						cur_chr_haplocoded_var_regs->at(all_i_var)->start,
						cur_chr_haplocoded_var_regs->at(all_i_var)->end,
						cur_chr_vic_params_var_regs->at(z_var)->start,
						cur_chr_vic_params_var_regs->at(z_var)->end,
						cur_chr_vic_params_var_regs->at(z_var)->name);
				}
			}
			/////////////////////

			// Reset the genotypes for this variant.
			memset(cur_subj_var_i_proxized_geno_sig, 0, sizeof(char) * (vecsize_t(sample_ids)));

			// Process each subject.
			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				// Process each haplotype.
				for (int i_hap = 0; i_hap < 2; i_hap++)
				{
					if (__DUMP_PROXIZATION_MSGS__)
					{
						// Write only for 1st subject.
						if (i_s == 0)
						{
							fprintf(stderr, "Variant @%d (%s:%d-%d::%s): Using window %d-%d\n",
								vicin_param_reg_i,
								cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->chrom,
								cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->start,
								cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->end,
								cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->name,
								cur_win_start, cur_win_end);
						} // first subject check.
					} // msg check.

					// Calculate average.
					int cur_var_i_vicinity_allele_sum = 0;
					for (int j_var = cur_win_start; j_var <= cur_win_end; j_var++)
					{
						void** var_j_reg_info = (void**)(cur_chr_vic_params_var_regs->at(j_var)->data);
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
								cur_chr_vic_params_var_regs->at(j_var)->name,
								cur_j_allele);
							exit(1);
						}

						//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(all_i_var)->start == 63244)
						if (__DUMP_PROXIZATION_MSGS__)
						{
							fprintf(stderr, "@ %d: j_var: %d [weight: %d; j_all: %d]; cumulative: %d\n",
								cur_chr_haplocoded_var_regs->at(all_i_var)->start,
								cur_chr_vic_params_var_regs->at(j_var)->start,
								cur_int_term, cur_j_allele,
								cur_var_i_vicinity_allele_sum);
						}

						// Add the interaction terms of this allele: Make sure they cover only one half.
						for (int k_var = cur_win_start; k_var < j_var; k_var++)
						{
							void** var_k_reg_info = (void**)(cur_chr_vic_params_var_regs->at(k_var)->data);
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

							//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(all_i_var)->start == 63244)
							if (__DUMP_PROXIZATION_MSGS__)
							{
								fprintf(stderr, "@ %d: {j-k}: {%d-%d} [j-k weight: %d]; cumulative: %d\n",
									cur_chr_haplocoded_var_regs->at(all_i_var)->start,
									cur_chr_vic_params_var_regs->at(j_var)->start,
									cur_chr_vic_params_var_regs->at(k_var)->start,
									cur_int_term, cur_var_i_vicinity_allele_sum);
							}

							// Add the final interaction.
							for (int l_var = cur_win_start; l_var < k_var; l_var++)
							{
								void** var_l_reg_info = (void**)(cur_chr_vic_params_var_regs->at(l_var)->data);
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

								//if (i_s == 0 && cur_chr_haplocoded_var_regs->at(all_i_var)->start == 63244)
								if(__DUMP_PROXIZATION_MSGS__)
								{
									fprintf(stderr, "@ %d: {j-k-l}: {%d-%d-%d} [j-k-l weight: %d]; cumulative: %d\n",
										cur_chr_haplocoded_var_regs->at(all_i_var)->start,
										cur_chr_vic_params_var_regs->at(j_var)->start,
										cur_chr_vic_params_var_regs->at(k_var)->start,
										cur_chr_vic_params_var_regs->at(l_var)->start,
										cur_int_term, cur_var_i_vicinity_allele_sum);
								}
							} // l_var loop.
						} // k_var loop.
					} // j_var loop.

					// Add error:
					int allele_err = 0;
					if (rng->random_double_ran3() < allele_err_eps)
					{
						cur_chr_vic_params_var_regs->at(vicin_param_reg_i)->dbl_score++;
						allele_err = 1;
					}

					// Add error.
					cur_var_i_vicinity_allele_sum += allele_err;

					// Take the modulus.
					int cur_vic_mod_allele = cur_var_i_vicinity_allele_sum % avg_geno_mod;

					// Set the allele in the genotype.
					cur_subj_var_i_proxized_geno_sig[i_s] = cur_subj_var_i_proxized_geno_sig[i_s] | (cur_vic_mod_allele << i_hap);
				} // i_hap loop.
			} // i_s loop.
		} // weight existence check; i.e., untyped variants.
		// Save this variant.
		int bytes_written = gzwrite(f_cur_thread_proxized_geno_mat, cur_subj_var_i_proxized_geno_sig, vecsize(sample_ids));
		if (bytes_written != vecsize(sample_ids))
		{
			int err_no = 0;
			fprintf(stderr, "Error during compression: %s", gzerror(f_cur_thread_proxized_geno_mat, &err_no));
			gzclose(f_cur_thread_proxized_geno_mat);
			exit(1);
		}
	} // i_var loop.

	t_string::print_padded_string(stderr, '\r', 100, "Thread-%d finished, closing..", thread_i);
	gzclose(f_cur_thread_proxized_geno_mat);

	return NULL;
} // callback_filter_proxize_variants_per_vicinity_no_buffer function.





// This function tests the proxization of variants (genotypes) using modular average of surrounding variants with interacting terms, i.e., AND operations.
// UPDATE 6/24: We added a new parameter to focus on regions of interest to be proxized; this is necessary when we are doing encoding only on, for example, typed variants. Rather than separating from the panel, this function can perform this operation in place.
// UPDATE 6/24: This function does not anonymize the coordinates/genetic maps, this should be handled in other mechanisms.
// UPDATE 6/24: We do testing on matrix representation to optimize performance. Use unstructured data if possible to improve performance.
// UPDATE 6/24: Output can be a prefix to save as matrix representation, which should make it faster to process.
void MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
	int n_threads,
	char* proxized_haplocoded_var_regs_op_fp)
{
	if (!check_file(per_var_vicinity_weight_params_fp))
	{
		fprintf(stderr, "Model weights are not found @ \"%s\"\n", per_var_vicinity_weight_params_fp);
		exit(1);
	}

	fprintf(stderr, "Loading per-variant vicinity parameters.\n");
	int loaded_n_vicinity = -1;
	int loaded_avg_geno_mod = -1;
	vector<t_annot_region*>* per_variant_vicinity_params = load_per_variant_site_mixing_parameters(per_var_vicinity_weight_params_fp, loaded_n_vicinity, loaded_avg_geno_mod);

	fprintf(stderr, "Loaded parameters for %d-vicinity and mod-%d for %d variants.\n", loaded_n_vicinity, loaded_avg_geno_mod, (int)per_variant_vicinity_params->size());

	int n_vicinity = loaded_n_vicinity;
	int avg_geno_mod = loaded_avg_geno_mod;

	fprintf(stderr, "Proxizing %s (%s) per variant-specific non-linear modular combinations using:\n\\sum_{vic @ [-%d, +%d]) (w_i x g_i)+\\sum_{vic @ [-%d, +%d] w_ij x g_i x g_j+\\sum_{vic @ [-%d, +%d] w_ijk x g_i x g_j x g_k mod %d using %d variant parameters.\n",
		haplocoded_geno_signal_fp, sample_ids_list_fp,
		n_vicinity, n_vicinity, n_vicinity, n_vicinity, n_vicinity, n_vicinity, avg_geno_mod, (int)per_variant_vicinity_params->size());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* all_haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	for (int i_var = 0; i_var < vecsize(all_haplocoded_var_regs); i_var++)
	{
		t_annot_region* var_reg = all_haplocoded_var_regs->at(i_var);

		void** var_reg_info = (void**)(var_reg->data);

		// We have to expand memory.
		void** new_reg_info = new void* [10];

		new_reg_info[0] = var_reg_info[0];
		new_reg_info[1] = NULL;
		new_reg_info[2] = NULL; // First degree vicinity weights.
		new_reg_info[3] = NULL; // Second degree vicinity weights.
		new_reg_info[4] = NULL; // Third degree vicinity weights.
		new_reg_info[5] = NULL; // Set the vic params region to null.

		var_reg->data = new_reg_info;

		var_reg->dbl_score = 0; // # of errors.
	} // i_var loop.

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	// Validate the phasing of variants.
	int max_geno = get_max_genotype_value(all_haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Divide the variants into chromosomes: We restructure all variants here.
	t_restr_annot_region_list* restr_all_haplocoded_var_regs = restructure_annot_regions(all_haplocoded_var_regs);

	// Start processing the proxized genotypes.
	for (int i_chr = 0; i_chr < vecsize(restr_all_haplocoded_var_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_all_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_haplocoded_var_regs = restr_all_haplocoded_var_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chr_vic_params_var_regs = get_regions_per_chromosome(per_variant_vicinity_params, restr_all_haplocoded_var_regs->chr_ids->at(i_chr));

		// This is very important -- if variants are not sorted, algorithm does not perform proxizing on neighboring variants.
		fprintf(stderr, "Sorting the variants per position..\n");
		sort(cur_chr_haplocoded_var_regs->begin(), cur_chr_haplocoded_var_regs->end(), sort_regions);
		sort(cur_chr_vic_params_var_regs->begin(), cur_chr_vic_params_var_regs->end(), sort_regions);

		// Set all genotype region model region indices to -1.
		for (int i_var = 0; i_var < vecsize(cur_chr_haplocoded_var_regs); i_var++)
		{
			cur_chr_haplocoded_var_regs->at(i_var)->score = -1;
		} // i_var loop.

		for (int i_var = 0; i_var < vecsize(cur_chr_vic_params_var_regs); i_var++)
		{
			cur_chr_vic_params_var_regs->at(i_var)->score = i_var;
		} // i_var loop.

		fprintf(stderr, "Assigning per-variant vicinity parameters.\n");
		vector<t_annot_region*>* intersects = intersect_annot_regions(cur_chr_vic_params_var_regs, cur_chr_haplocoded_var_regs, true);
		if (vecsize(cur_chr_vic_params_var_regs) != vecsize(intersects))
		{
			fprintf(stderr, "Sanity check failed @ %s(%d): Some of the vicinity hashing tag variants are missing or duplicated: %d/%d\n",
				__FILE__, __LINE__, vecsize(cur_chr_vic_params_var_regs), vecsize(intersects));
			exit(1);
		}

		fprintf(stderr, "Processing %d intersects.\n", (int)intersects->size());
		vector<t_annot_region*>* cur_chr_geno_var_regs_w_vic_params = new vector<t_annot_region*>();
		for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* vic_params_reg = int_info->src_reg;
			t_annot_region* var_reg = int_info->dest_reg;

			if (vic_params_reg->start == var_reg->start &&
				vic_params_reg->end == var_reg->end)
			{
				void** vic_params_reg_info = (void**)(vic_params_reg->data);

				void** var_reg_info = (void**)(var_reg->data);

				void** new_reg_info = new void* [10];

				new_reg_info[0] = var_reg_info[0];
				new_reg_info[1] = NULL;
				new_reg_info[2] = vic_params_reg_info[0]; // First degree vicinity weights.
				new_reg_info[3] = vic_params_reg_info[1]; // Second degree vicinity weights.
				new_reg_info[4] = vic_params_reg_info[2]; // Third degree vicinity weights.

				// This region is typed, has vicinity parameters. Below we set the index for each of these regions, which is used while saving genotype signals.
				cur_chr_geno_var_regs_w_vic_params->push_back(var_reg);
				
				var_reg->data = new_reg_info;
				var_reg->dbl_score = 0; // # of errors.

				//focus_haplocoded_regions->push_back(var_reg);
			}
		} // i_int loop.
		fprintf(stderr, "Proxizing %d/%d variant regions.\n", vecsize(cur_chr_geno_var_regs_w_vic_params), vecsize(all_haplocoded_var_regs));

		// Assign the vicinity parameter index.
		for (int i_var = 0; i_var < vecsize(cur_chr_geno_var_regs_w_vic_params); i_var++)
		{
			cur_chr_geno_var_regs_w_vic_params->at(i_var)->score = i_var;
		} // i_var loop.

		int n_vars_per_thread = MAX(1, ceil(vecsize(cur_chr_haplocoded_var_regs) / ((double)n_threads)));

		// ADD THREADING HERE.
		int cur_var_i_start = 0;
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int thread_i = 0; thread_i < n_threads; thread_i++)
		{
			if (cur_var_i_start > vecsize(cur_chr_haplocoded_var_regs))
			{
				break;
			}

			int cur_var_i_end = MIN((cur_var_i_start + n_vars_per_thread), vecsize(cur_chr_haplocoded_var_regs));

			void** cur_thread_info_ptr = new void* [10];
			int* int_vals = new int[10];
			int_vals[0] = thread_i;
			int_vals[1] = n_threads;
			int_vals[2] = cur_var_i_start;
			int_vals[3] = cur_var_i_end;
			int_vals[4] = avg_geno_mod;
			int_vals[5] = n_vicinity;
			cur_thread_info_ptr[0] = int_vals;

			double* dbl_vals = new double[10];
			cur_thread_info_ptr[1] = dbl_vals;
			dbl_vals[0] = allele_err_eps;

			cur_thread_info_ptr[2] = cur_chr_haplocoded_var_regs;

			cur_thread_info_ptr[3] = cur_chr_geno_var_regs_w_vic_params;

			cur_thread_info_ptr[4] = sample_ids;

			cur_thread_info_ptr[5] = new t_rng(t_seed_manager::seed_me() + thread_i);

			t_ansi_thread* cur_thread = new t_ansi_thread(callback_filter_proxize_variants_per_vicinity_no_buffer, cur_thread_info_ptr);

			t_string::print_padded_string(stderr, '\r', 100, "Starting thread %d: [%d-%d]", thread_i, cur_var_i_start, cur_var_i_end);

			cur_thread->run_thread();
			threads->push_back(cur_thread);

			cur_var_i_start += n_vars_per_thread;
		} // i_th loop.

		t_string::print_padded_string(stderr, '\n', 100, "Waiting for threads..");
		vector<char*>* per_thread_matrix_file = new vector<char*>();
		for (int thread_i = 0; thread_i < vecsize(threads); thread_i++)
		{
			threads->at(thread_i)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Thread %d finished..", thread_i);

			char cur_thread_proxized_geno_mat_fp[1000];
			sprintf(cur_thread_proxized_geno_mat_fp, "proxized_geno_sig_%d.matrix.gz", thread_i);

			if (!check_file(cur_thread_proxized_geno_mat_fp))
			{
				fprintf(stderr, "Sanity check failed @ %s(%d): Matrix file was not found @ %s", __FILE__, __LINE__, cur_thread_proxized_geno_mat_fp);
				exit(1);
			}

			per_thread_matrix_file->push_back(t_string::copy_me_str(cur_thread_proxized_geno_mat_fp));
			//fprintf(stderr, "Added %d. file name: %s          \r", thread_i, per_thread_matrix_file->back());
		} // thread_i loop.

		// Start saving genotypes.
		t_string::print_padded_string(stderr, '\n', 100, "Concatenating %d matrix files..", vecsize(per_thread_matrix_file));
		char op_matrix_fp[1000];
		sprintf(op_matrix_fp, "%s_genotypes.matrix.gz", proxized_haplocoded_var_regs_op_fp);
		concatenateGzipFiles(op_matrix_fp, per_thread_matrix_file);

		// Delete matrix files.
		fprintf(stderr, "Deleting temp matrix files..\n");
		for (int thread_i = 0; thread_i < vecsize(threads); thread_i++)
		{
			char cur_thread_proxized_geno_mat_fp[1000];
			sprintf(cur_thread_proxized_geno_mat_fp, "proxized_geno_sig_%d.matrix.gz", thread_i);
			//fprintf(stderr, "Deleting %s\n", cur_thread_proxized_geno_mat_fp);
			delete_file(cur_thread_proxized_geno_mat_fp);
		} // thread_i loop.

		// Save the variants.
		char cur_matrix_variants_BED_fp[1000];
		sprintf(cur_matrix_variants_BED_fp, "%s_variants.bed", proxized_haplocoded_var_regs_op_fp);
		dump_BED(cur_matrix_variants_BED_fp, cur_chr_haplocoded_var_regs);

		// Save the subject identifiers.
		char cur_matrix_subject_list_fp[1000];
		sprintf(cur_matrix_subject_list_fp, "%s_subjects.list", proxized_haplocoded_var_regs_op_fp);
		save_lines(sample_ids, cur_matrix_subject_list_fp);

		fprintf(stderr, "Processing only 1 chromosome..\n");
		break;
	} // i_chr loop.

	//// Write the error statistics.
	//char var_stats_fp[1000];
	//sprintf(var_stats_fp, "%s_var_stats.txt", proxized_haplocoded_var_regs_op_fp);
	//FILE* f_var_stats = open_f(var_stats_fp, "w");
	//for (int i_var = 0; i_var < vecsize(focus_haplocoded_regions); i_var++)
	//{
	//	fprintf(f_var_stats, "%s\t%d\t%d\t%s\t%d\n", focus_haplocoded_regions->at(i_var)->chrom,
	//		focus_haplocoded_regions->at(i_var)->start, focus_haplocoded_regions->at(i_var)->end,
	//		focus_haplocoded_regions->at(i_var)->name,
	//		(int)(focus_haplocoded_regions->at(i_var)->dbl_score));
	//} // i_var loop.
	//close_f(f_var_stats, NULL);
} // proxize_variants_per_vicinity_modular_average function.