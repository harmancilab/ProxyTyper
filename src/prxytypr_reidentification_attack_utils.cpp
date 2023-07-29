#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "prxytypr_file_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_vector_macros.h"
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
#include "prxytypr_proxytyper.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

const bool __DUMP_REIDENTIFICATION_ATTACK_MSGS__ = false;

void linking_attack_per_haplotype_frequency_signatures_per_ref_panel(char* target_query_panel_matbed_fp, char* target_query_sample_list_fp,
	char* ref_panel_matbed_fp, char* ref_panel_sample_list_fp,
	char* target_proxy_panel_matbed_fp, char* target_proxy_panel_sample_list_fp,
	int n_vicinity,
	int n_step,
	char* op_fp)
{
	vector<char*>* query_sample_ids = buffer_file(target_query_sample_list_fp);
	vector<t_annot_region*>* query_var_regs = load_variant_signal_regions_wrapper(target_query_panel_matbed_fp, target_query_sample_list_fp);

	vector<char*>* ref_sample_ids = buffer_file(ref_panel_sample_list_fp);
	vector<t_annot_region*>* ref_var_regs = load_variant_signal_regions_wrapper(ref_panel_matbed_fp, ref_panel_sample_list_fp);

	vector<char*>* proxy_sample_ids = buffer_file(target_proxy_panel_sample_list_fp);
	vector<t_annot_region*>* proxy_var_regs = load_variant_signal_regions_wrapper(target_proxy_panel_matbed_fp, target_proxy_panel_sample_list_fp);

	// This is necessary to make sure we are following the haplotypes correctly.
	sort(query_var_regs->begin(), query_var_regs->end(), sort_regions);
	sort(ref_var_regs->begin(), ref_var_regs->end(), sort_regions);
	sort(proxy_var_regs->begin(), proxy_var_regs->end(), sort_regions);

	int max_geno_proxy = get_max_genotype_value(proxy_var_regs, proxy_sample_ids);
	int max_geno_ref = get_max_genotype_value(ref_var_regs, ref_sample_ids);
	int max_geno_query = get_max_genotype_value(query_var_regs, query_sample_ids);

	if (max_geno_proxy == 3 &&
		max_geno_query == 3 &&
		max_geno_ref == 3)
	{
		fprintf(stderr, "All panels are haplocoded\n");
	}
	else
	{
		fprintf(stderr, "One of the panels is likely genocoded (%d, %d), make sure they are both haplocoded (i.e., phased)\n", max_geno_ref, max_geno_query);
		exit(0);
	}

	if (proxy_var_regs->size() != ref_var_regs->size() ||
		proxy_var_regs->size() != query_var_regs->size())
	{
		fprintf(stderr, "Panel variant sizes are not matching (%d vs %d)\n", vecsize(ref_var_regs), vecsize(query_var_regs));

		exit(1);
	}

	// Reference and query must match in coordinates.
	for (int i_var = 0; i_var < (int)ref_var_regs->size(); i_var++)
	{
		if (ref_var_regs->at(i_var)->start != query_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Query/Reference variants are not matching @ %d. variant: %d vs %d\n",
				i_var,
				ref_var_regs->at(i_var)->start, query_var_regs->at(i_var)->start);

			exit(1);
		}
	} // i_var loop.

	vector<char**>* per_ref_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_var_regs, query_sample_ids);
	vector<char**>* per_proxy_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(proxy_var_regs, proxy_sample_ids);

	double** query_per_hap_per_var_ref_hap_freqs = new double* [2 * vecsize(query_sample_ids) + 2];
	for (int i_hap = 0; i_hap < 2 * vecsize(query_sample_ids); i_hap++)
	{
		query_per_hap_per_var_ref_hap_freqs[i_hap] = new double[vecsize(query_var_regs) + 2];
		memset(query_per_hap_per_var_ref_hap_freqs[i_hap], 0, sizeof(double) * vecsize(query_var_regs));
	} // i_hap loop.

	double** proxy_per_hap_per_var_hap_freqs = new double* [2 * vecsize(proxy_sample_ids) + 2];
	for (int i_hap = 0; i_hap < 2 * vecsize(proxy_sample_ids); i_hap++)
	{
		proxy_per_hap_per_var_hap_freqs[i_hap] = new double[vecsize(proxy_var_regs) + 2];
		memset(proxy_per_hap_per_var_hap_freqs[i_hap], 0, sizeof(double) * vecsize(proxy_var_regs));
	} // i_hap loop.

	for (int var_i = 0; var_i < vecsize(ref_var_regs); var_i += n_step)
	{
		if (var_i % 100 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vicinity, 0);
		int win_end_i = MIN(var_i + n_vicinity, vecsize(ref_var_regs) - 1);

		int n_vars_per_win = win_end_i - win_start_i + 1;

		// Get the reference window haplotypes: These are used for assigning frequencies to query haplotypes.
		vector<char*>* cur_win_ref_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < vecsize(ref_sample_ids); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref_haps->push_back(per_ref_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref_haps->begin(), cur_win_ref_haps->end(), sort_haplotypes);

		////////////////////////////////////////////////////////////////////////////////////////////////////

		// Get query window haplotypes.
		vector<char*>* cur_win_query_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < vecsize(query_sample_ids); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_query_unsorted_haps->push_back(per_query_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the proxy window haplotypes; we need both the sorted and unsorted haplotypes.
		vector<char*>* cur_win_proxy_haps = new vector<char*>();
		vector<char*>* cur_win_proxy_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < vecsize(proxy_sample_ids); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_proxy_haps->push_back(per_proxy_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_proxy_unsorted_haps->push_back(per_proxy_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_proxy_haps->begin(), cur_win_proxy_haps->end(), sort_haplotypes);

		////////////////////////////////////////////////////////////////////////////////

		// Count the sorted haplotypes to get frequencies.
		vector<int>* per_unique_ref_haplotype_cnts = new vector<int>();
		vector<char*>* ref_unique_haplotypes = new vector<char*>();
		int cur_hap_cnt = 1;
		char* cur_uniq_hap = cur_win_ref_haps->at(0);
		for (int i_hap = 1; i_hap < vecsize(cur_win_ref_haps); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref_haps->at(i_hap - 1), cur_win_ref_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref_haplotype_cnts->push_back(cur_hap_cnt);
				ref_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref_haplotype_cnts->push_back(cur_hap_cnt);
			ref_unique_haplotypes->push_back(cur_uniq_hap);
		}

		// Following does a check on the total of unique haplotype counts from previous loop.
		double hap_count_check1 = 0;
		for (int i_hap = 0; i_hap < (int)ref_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check1 += per_unique_ref_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check1 != vecsize(ref_sample_ids) * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the unique proxy haplotypes.
		vector<int>* per_unique_proxy_haplotype_cnts = new vector<int>();
		vector<char*>* proxy_unique_haplotypes = new vector<char*>();
		cur_hap_cnt = 1;
		cur_uniq_hap = cur_win_proxy_haps->at(0);
		for (int i_hap = 1; i_hap < vecsize(cur_win_proxy_haps); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_proxy_haps->at(i_hap - 1), cur_win_proxy_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_proxy_haplotype_cnts->push_back(cur_hap_cnt);
				proxy_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_proxy_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_proxy_haplotype_cnts->push_back(cur_hap_cnt);
			proxy_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check2 = 0;
		for (int i_hap = 0; i_hap < vecsize(proxy_unique_haplotypes); i_hap++)
		{
			hap_count_check2 += per_unique_proxy_haplotype_cnts->at(i_hap);
		} // i_hap loop.

		if (hap_count_check2 != (int)proxy_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		// We need to map the proxy haplotypes to unique haplotypes.
		vector<vector<int>*>* per_unique_proxy_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_proxy_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_proxy_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_proxy_unsorted_haps->at(i_hap), proxy_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.

			per_unique_proxy_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_proxy_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref2: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_proxy_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update the hapfreq signatures at this variant:
		// For the query, use the reference frequencies, for the proxy, use its own frequencies.
		for (int i_hap = 0; i_hap < vecsize(cur_win_query_unsorted_haps); i_hap++)
		{
			// Go over all reference haplotypes and compare.
			for (int i_uniq_hap = 0; i_uniq_hap < vecsize(per_unique_ref_haplotype_cnts); i_uniq_hap++)
			{
				if (compare_haplotypes(ref_unique_haplotypes->at(i_uniq_hap), cur_win_query_unsorted_haps->at(i_hap), n_vars_per_win))
				{
					query_per_hap_per_var_ref_hap_freqs[i_hap][var_i] = per_unique_ref_haplotype_cnts->at(i_uniq_hap);
				}
			} // i_uniq_hap loop.
		} // i_hap loop.

		// For proxy, use its own frequencies.
		for (int i_uniq_hap = 0; i_uniq_hap < vecsize(per_unique_proxy_haplotype_cnts); i_uniq_hap++)
		{
			vector<int>* cur_orig_hap_indices = per_unique_proxy_haplotype_subject_i->at(i_uniq_hap);

			for (int cur_i_hap = 0; cur_i_hap < vecsize(cur_orig_hap_indices); cur_i_hap++)
			{
				proxy_per_hap_per_var_hap_freqs[cur_orig_hap_indices->at(cur_i_hap)][var_i] = per_unique_proxy_haplotype_cnts->at(i_uniq_hap);
			}
		} // i_uniq_hap loop.

		// Free memory.
		delete cur_win_ref_haps;
		delete cur_win_query_unsorted_haps;

		delete per_unique_ref_haplotype_cnts;
		delete ref_unique_haplotypes;

		delete per_unique_proxy_haplotype_cnts;
		delete proxy_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_proxy_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_proxy_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_proxy_haplotype_subject_i;
	} // var_i loop.

	// Now do linking based on the frequency signatures.
	fprintf(stderr, "Comparing haplotype frequency signatures.\n");
	FILE* f_op = open_f(op_fp, "w");
	for (int i_hap = 0; i_hap < 2 * vecsize(query_sample_ids); i_hap++)
	{
		double cur_hap_max_corr = -1;
		int cur_hap_max_corr_j_hap = -1;
		double cur_hap_max_corr_2nd = -1;
		for (int j_hap = 0; j_hap < 2 * vecsize(proxy_sample_ids); j_hap++)
		{
			double cur_corr;
			get_correlation(query_per_hap_per_var_ref_hap_freqs[i_hap], proxy_per_hap_per_var_hap_freqs[j_hap],
							vecsize(ref_var_regs),
							cur_corr);

			if (cur_corr > cur_hap_max_corr)
			{
				cur_hap_max_corr_2nd = cur_hap_max_corr;

				cur_hap_max_corr = cur_corr;
				cur_hap_max_corr_j_hap = j_hap;
			}
		} // j_hap

		fprintf(f_op, "%d\t%d\t%.5f\t%.5f\n", i_hap, cur_hap_max_corr_j_hap, cur_hap_max_corr, cur_hap_max_corr_2nd);
		fprintf(stderr, "%d\t%d\t%.5f\t%.5f\n", i_hap, cur_hap_max_corr_j_hap, cur_hap_max_corr, cur_hap_max_corr_2nd);
	} // j_hap

	close_f(f_op, op_fp);
} // linking_attack_per_haplotype_frequency_signatures_per_ref_panel function.



void linking_attack_per_haplotype_frequency_signatures(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vicinity,
	int n_step,
	char* op_fp)
{
	vector<char*>* ref1_sample_ids = buffer_file(ref1_sample_list_fp);
	vector<t_annot_region*>* ref1_var_regs = load_variant_signal_regions_wrapper(ref1_panel_matbed_fp, ref1_sample_list_fp);

	vector<char*>* ref2_sample_ids = buffer_file(ref2_sample_list_fp);
	vector<t_annot_region*>* ref2_var_regs = load_variant_signal_regions_wrapper(ref2_panel_matbed_fp, ref2_sample_list_fp);

	// This is necessary to make sure we are following the haplotypes correctly.
	sort(ref1_var_regs->begin(), ref1_var_regs->end(), sort_regions);
	sort(ref2_var_regs->begin(), ref2_var_regs->end(), sort_regions);

	int max_geno_ref1 = get_max_genotype_value(ref1_var_regs, ref1_sample_ids);
	int max_geno_ref2 = get_max_genotype_value(ref2_var_regs, ref2_sample_ids);

	if (max_geno_ref1 == 3 &&
		max_geno_ref2 == 3)
	{
		fprintf(stderr, "Both panels are haplocoded\n");
	}
	else
	{
		fprintf(stderr, "One of the panels is likely genocoded (%d, %d), make sure they are both haplocoded (i.e., phased)\n", max_geno_ref1, max_geno_ref2);
		exit(0);
	}

	if (ref1_var_regs->size() != ref2_var_regs->size())
	{
		fprintf(stderr, "Panel variant sizes are not matching (%d vs %d)\n", (int)ref1_var_regs->size(), (int)ref2_var_regs->size());

		exit(1);
	}

	// We are relaxing this condition.
	//for (int i_var = 0; i_var < (int)ref1_var_regs->size(); i_var++)
	//{
	//	if (ref1_var_regs->at(i_var)->start != ref2_var_regs->at(i_var)->start)
	//	{
	//		fprintf(stderr, "Variants are not matching @ %d. variant: %d vs %d\n",
	//			i_var,
	//			ref1_var_regs->at(i_var)->start, ref2_var_regs->at(i_var)->start);

	//		exit(1);
	//	}
	//} // i_var loop.

	vector<char**>* per_ref1_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref1_var_regs, ref1_sample_ids);
	vector<char**>* per_ref2_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref2_var_regs, ref2_sample_ids);

	//t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	double** ref1_per_hap_per_var_hap_freqs = new double* [2*vecsize(ref1_sample_ids) + 2];
	for (int i_hap = 0; i_hap < 2 * vecsize(ref1_sample_ids); i_hap++)
	{
		ref1_per_hap_per_var_hap_freqs[i_hap] = new double [vecsize(ref1_var_regs) + 2];
		memset(ref1_per_hap_per_var_hap_freqs[i_hap], 0, sizeof(double) * vecsize(ref1_var_regs));
	} // i_hap loop.

	double** ref2_per_hap_per_var_hap_freqs = new double* [2 * vecsize(ref2_sample_ids) + 2];
	for (int i_hap = 0; i_hap < 2 * vecsize(ref2_sample_ids); i_hap++)
	{
		ref2_per_hap_per_var_hap_freqs[i_hap] = new double[vecsize(ref2_var_regs) + 2];
		memset(ref2_per_hap_per_var_hap_freqs[i_hap], 0, sizeof(double) * vecsize(ref2_var_regs));
	} // i_hap loop.


	FILE* f_op = open_f(op_fp, "w");
	for (int var_i = 0; var_i < (int)ref1_var_regs->size(); var_i += n_step)
	{
		if (var_i % 100 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vicinity, 0);
		int win_end_i = MIN(var_i + n_vicinity, (int)ref1_var_regs->size() - 1);

		int n_vars_per_win = win_end_i - win_start_i + 1;

		vector<char*>* cur_win_ref1_haps = new vector<char*>();
		vector<char*>* cur_win_ref1_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref1_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref1_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref1_unsorted_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref1_haps->begin(), cur_win_ref1_haps->end(), sort_haplotypes);

		// Get ref2 haplotypes.
		vector<char*>* cur_win_ref2_haps = new vector<char*>();
		vector<char*>* cur_win_ref2_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref2_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref2_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref2_unsorted_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref2_haps->begin(), cur_win_ref2_haps->end(), sort_haplotypes);

		////////////////////////////////////////////////////////////////////////////////

		// Count the sorted haplotypes to get frequencies.
		vector<int>* per_unique_ref1_haplotype_cnts = new vector<int>();
		vector<char*>* ref1_unique_haplotypes = new vector<char*>();
		int cur_hap_cnt = 1;
		char* cur_uniq_hap = cur_win_ref1_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref1_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref1_haps->at(i_hap - 1), cur_win_ref1_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
				ref1_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref1_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
			ref1_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check1 = 0;
		for (int i_hap = 0; i_hap < (int)ref1_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check1 += per_unique_ref1_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check1 != (int)ref1_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		// Assign the subject indices to each of the unique haplotypes.
		vector<vector<int>*>* per_unique_ref1_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref1_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref1_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref1: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref1_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		vector<int>* per_unsorted_ref1_hap_freq = new vector<int>();
		for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
		{
			int cur_hap_cnt = -1;
			for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					if (cur_hap_cnt != -1)
					{
						fprintf(stderr, "Sanity check failed: Ref1 haplotype %d matches to multiple unique haplotypes.\n", i_hap);
						exit(1);
					}
					cur_hap_cnt = per_unique_ref1_haplotype_cnts->at(i_uniq_hap);
				}
			} // i_uniq_hap loop.

			if (cur_hap_cnt == -1)
			{
				fprintf(stderr, "Sanity check failed: Ref1 haplotype %d does not match to any unique haplotypes.\n", i_hap);
				exit(1);
			}
			per_unsorted_ref1_hap_freq->push_back(cur_hap_cnt);
		} // i_hap loop.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the unique ref2 haplotypes.
		vector<int>* per_unique_ref2_haplotype_cnts = new vector<int>();
		vector<char*>* ref2_unique_haplotypes = new vector<char*>();
		cur_hap_cnt = 1;
		cur_uniq_hap = cur_win_ref2_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref2_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref2_haps->at(i_hap - 1), cur_win_ref2_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
				ref2_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref2_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
			ref2_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check2 = 0;
		for (int i_hap = 0; i_hap < (int)ref2_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check2 += per_unique_ref2_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check2 != (int)ref2_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		vector<vector<int>*>* per_unique_ref2_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref2_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref2_unsorted_haps->at(i_hap), ref2_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref2_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref2_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref2: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref2_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Update the hapfreq signatures at this variant:
		for (int i_uniq_hap = 0; i_uniq_hap < vecsize(per_unique_ref1_haplotype_cnts); i_uniq_hap++)
		{
			vector<int>* cur_orig_hap_indices = per_unique_ref1_haplotype_subject_i->at(i_uniq_hap);

			for (int cur_i_hap = 0; cur_i_hap < vecsize(cur_orig_hap_indices); cur_i_hap++)
			{
				ref1_per_hap_per_var_hap_freqs[cur_orig_hap_indices->at(cur_i_hap)][var_i] = per_unique_ref1_haplotype_cnts->at(i_uniq_hap);
			}
		} // i_uniq_hap loop.


		for (int i_uniq_hap = 0; i_uniq_hap < vecsize(per_unique_ref2_haplotype_cnts); i_uniq_hap++)
		{
			vector<int>* cur_orig_hap_indices = per_unique_ref2_haplotype_subject_i->at(i_uniq_hap);

			for (int cur_i_hap = 0; cur_i_hap < vecsize(cur_orig_hap_indices); cur_i_hap++)
			{
				ref2_per_hap_per_var_hap_freqs[cur_orig_hap_indices->at(cur_i_hap)][var_i] = per_unique_ref2_haplotype_cnts->at(i_uniq_hap);
			}
		} // i_uniq_hap loop.

		// Free memory.
		delete cur_win_ref1_haps;
		delete cur_win_ref1_unsorted_haps;
		delete cur_win_ref2_haps;
		delete cur_win_ref2_unsorted_haps;

		delete per_unique_ref1_haplotype_cnts;
		delete ref1_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref1_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref1_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref1_haplotype_subject_i;

		delete per_unsorted_ref1_hap_freq;

		delete per_unique_ref2_haplotype_cnts;
		delete ref2_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref2_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref2_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref2_haplotype_subject_i;
	} // var_i loop.

	// Now do linking based on the frequency signatures.
	for (int i_hap = 0; i_hap < 2 * vecsize(ref1_sample_ids); i_hap++)
	{
		double cur_hap_max_corr = -1;
		int cur_hap_max_corr_j_hap = -1;
		double cur_hap_max_corr_2nd = -1;
		for (int j_hap = 0; j_hap < 2 * vecsize(ref1_sample_ids); j_hap++)
		{
			double cur_corr;
			get_correlation(ref1_per_hap_per_var_hap_freqs[i_hap], ref2_per_hap_per_var_hap_freqs[j_hap], 
							vecsize(ref1_var_regs),
							cur_corr);

			if (cur_corr > cur_hap_max_corr)
			{
				cur_hap_max_corr_2nd = cur_hap_max_corr;

				cur_hap_max_corr = cur_corr;
				cur_hap_max_corr_j_hap = j_hap;
			}
		} // j_hap

		fprintf(f_op, "%d\t%d\t%.5f\t%.5f\n", i_hap, cur_hap_max_corr_j_hap, cur_hap_max_corr, cur_hap_max_corr_2nd);
		fprintf(stderr, "%d\t%d\t%.5f\t%.5f\n", i_hap, cur_hap_max_corr_j_hap, cur_hap_max_corr, cur_hap_max_corr_2nd);
	} // j_hap

	close_f(f_op, op_fp);
}

void get_consecutive_block_variant_correlations(char* panel_target_geno_sig_regs_fp, char* panel_sample_list_fp, int l_block, int l_step, char* op_fp)
{
	vector<t_annot_region*>* panel_var_geno_sig_regs = load_variant_signal_regions_wrapper(panel_target_geno_sig_regs_fp, panel_sample_list_fp);
	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);

	sort(panel_var_geno_sig_regs->begin(), panel_var_geno_sig_regs->end(), sort_regions);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_reg = 0; i_reg < (int)panel_var_geno_sig_regs->size(); i_reg += l_step)
	{
		if (i_reg % 100 == 0)
		{
			fprintf(stderr, "i_reg = %d        \r", i_reg);
		}

		//size_t min_j_reg = MAX(0, i_reg - l_block);
		size_t min_j_reg = i_reg + 1;
		size_t max_j_reg = MIN((int)panel_var_geno_sig_regs->size()-1, i_reg + l_block);
		for (size_t j_reg = min_j_reg; j_reg < max_j_reg; j_reg += l_step)
		{
			//if (j_reg <= i_reg)
			//{
			//	continue;
			//}
			// Correlate the variant genotypes.
			t_annot_region* cur_reg = panel_var_geno_sig_regs->at(i_reg);
			void** cur_reg_info = (void**)(cur_reg->data);
			char* cur_reg_geno = (char*)(cur_reg_info[0]);
			t_annot_region* prev_reg = panel_var_geno_sig_regs->at(j_reg);
			void** prev_reg_info = (void**)(prev_reg->data);
			char* prev_reg_geno = (char*)(prev_reg_info[0]);

			double* cur_reg_geno_sig_dbl = new double[panel_sample_ids->size() + 2];
			double* prev_reg_geno_sig_dbl = new double[panel_sample_ids->size() + 2];
			for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
			{
				cur_reg_geno_sig_dbl[i_s] = (double)(get_genotype_per_haplocoded_genotype(cur_reg_geno[i_s]));
				prev_reg_geno_sig_dbl[i_s] = (double)(get_genotype_per_haplocoded_genotype(prev_reg_geno[i_s]));
			} // i_s loop.

			double var2var_corr = 0;
			get_correlation(cur_reg_geno_sig_dbl, prev_reg_geno_sig_dbl, (int)panel_sample_ids->size(), var2var_corr);
			fprintf(f_op, "%d\t%d\t%.4f\n",
					(int)i_reg, (int)j_reg, var2var_corr);

			fprintf(f_op, "%d\t%d\t%.4f\n",
				(int)j_reg, (int)i_reg, var2var_corr);

			delete[] cur_reg_geno_sig_dbl;
			delete[] prev_reg_geno_sig_dbl;
		} // j_reg loop.
	} // i_reg loop.

	close_f(f_op, op_fp);
}

void get_consecutive_variant_correlations(char* panel_target_geno_sig_regs_fp, char* panel_sample_list_fp)
{
	vector<t_annot_region*>* panel_var_geno_sig_regs = load_variant_signal_regions_wrapper(panel_target_geno_sig_regs_fp, panel_sample_list_fp);
	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);

	sort(panel_var_geno_sig_regs->begin(), panel_var_geno_sig_regs->end(), sort_regions);

	FILE* f_op = open_f("var2var_corrs.bed", "w");
	for (int i_reg = 1; i_reg < (int)panel_var_geno_sig_regs->size(); i_reg++)
	{
		// Correlate the variant genotypes.
		t_annot_region* cur_reg = panel_var_geno_sig_regs->at(i_reg);
		void** cur_reg_info = (void**)(cur_reg->data);
		char* cur_reg_geno = (char*)(cur_reg_info[0]);
		t_annot_region* prev_reg = panel_var_geno_sig_regs->at(i_reg - 1);
		void** prev_reg_info = (void**)(prev_reg->data);
		char* prev_reg_geno = (char*)(prev_reg_info[0]);

		double* cur_reg_geno_sig_dbl = new double[panel_sample_ids->size() + 2];
		double* prev_reg_geno_sig_dbl = new double[panel_sample_ids->size() + 2];
		for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
		{
			cur_reg_geno_sig_dbl[i_s] = (double)(cur_reg_geno[i_s]);
			prev_reg_geno_sig_dbl[i_s] = (double)(prev_reg_geno[i_s]);
		} // i_s loop.

		double var2var_corr = 0;
		get_correlation(cur_reg_geno_sig_dbl, prev_reg_geno_sig_dbl, (int)panel_sample_ids->size(), var2var_corr);
		fprintf(f_op, "%s\t%d\t%d\t%.4f\n",
			panel_var_geno_sig_regs->at(i_reg)->chrom,
			panel_var_geno_sig_regs->at(i_reg)->start,
			panel_var_geno_sig_regs->at(i_reg)->end,
			var2var_corr);

		delete[] cur_reg_geno_sig_dbl;
		delete[] prev_reg_geno_sig_dbl;
	} // i_reg loop.

	close_f(f_op, NULL);
}

// Correlaton of haplotype freq's in different proxizations.
void get_proxization1_2_proxization_mapped_haplotype_frequency_correlation(char* haplocoded_geno_signal_fp,
	char* proxy1_haplocoded_geno_signal_fp,
	char* proxy2_haplocoded_geno_signal_fp,
	char* sample_ids_list_fp,
	int n_vicinity,
	char* op_prefix)
{
	fprintf(stderr, "Calculating the correlation between the haplotype frequencies of same encodings of kmers between proxizations of %s (%s):\n\
	Proxization 1: %s\n\
	Proxization 2: %s\n", haplocoded_geno_signal_fp, sample_ids_list_fp,
proxy1_haplocoded_geno_signal_fp, proxy2_haplocoded_geno_signal_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	vector<t_annot_region*>* cleartext_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	
	vector<t_annot_region*>* proxy1_var_regs = load_variant_signal_regions_wrapper(proxy1_haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<t_annot_region*>* proxy2_var_regs = load_variant_signal_regions_wrapper(proxy2_haplocoded_geno_signal_fp, sample_ids_list_fp);

	sort(proxy1_var_regs->begin(), proxy1_var_regs->end(), sort_regions);
	sort(proxy2_var_regs->begin(), proxy2_var_regs->end(), sort_regions);
	sort(cleartext_var_regs->begin(), cleartext_var_regs->end(), sort_regions);

	vector<char**>* per_cleartext_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(cleartext_var_regs, sample_ids);
	vector<char**>* per_proxy1_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(proxy1_var_regs, sample_ids);
	vector<char**>* per_proxy2_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(proxy2_var_regs, sample_ids);

	// At each variant's vicinity, get the frequencies of proxy1 and proxy2 kmers aligned to cleartext kmers and correlate them.
	char model1_proxy_2_model2_proxy_hapfreq_corr_stats_fp[1000];
	sprintf(model1_proxy_2_model2_proxy_hapfreq_corr_stats_fp, "%s_model1_proxy_2_model2_proxy_hapfreq_corrs.bed", op_prefix);
	FILE* f_model1proxy_2_model2proxy_hapfreq_corr_stats = open_f(model1_proxy_2_model2_proxy_hapfreq_corr_stats_fp, "w");
	for (int var_i = 0; var_i < vecsize(cleartext_var_regs); var_i++)
	{
		if (var_i % 1000 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vicinity, 0);
		int win_end_i = MIN(var_i + n_vicinity, vecsize(cleartext_var_regs) - 1);

		vector<t_kmer*>* cur_win_cleartext_kmers = extract_kmers_per_haplotype(per_cleartext_subj_haplotypes, win_start_i, win_end_i);
		vector<t_kmer*>* cur_win_proxy1_kmers = extract_kmers_per_haplotype(per_proxy1_subj_haplotypes, win_start_i, win_end_i);
		vector<t_kmer*>* cur_win_proxy2_kmers = extract_kmers_per_haplotype(per_proxy2_subj_haplotypes, win_start_i, win_end_i);

		vector<t_kmer*>* cur_win_unique_cleartext_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_cleartext_kmer_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_cleartext_kmers, cur_win_unique_cleartext_kmers, cur_win_unique_cleartext_kmer_cnts);

		vector<t_kmer*>* cur_win_unique_proxy1_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_proxy1_kmer_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_proxy1_kmers, cur_win_unique_proxy1_kmers, cur_win_unique_proxy1_kmer_cnts);

		vector<t_kmer*>* cur_win_unique_proxy2_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_proxy2_kmer_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_proxy2_kmers, cur_win_unique_proxy2_kmers, cur_win_unique_proxy2_kmer_cnts);

		if (vecsize(cur_win_cleartext_kmers) != 2 * vecsize(sample_ids) ||
			vecsize(cur_win_proxy2_kmers) != 2 * vecsize(sample_ids) ||
			vecsize(cur_win_proxy2_kmers) != 2 * vecsize(sample_ids))
		{
			fprintf(stderr, "Sanity check failed: The kmer size is not as expected.\n");
			exit(0);
		}

		double proxy1_cnts[1000];
		double proxy2_cnts[1000];
		for (int i_unique_clear_kmer = 0; i_unique_clear_kmer < vecsize(cur_win_unique_cleartext_kmers); i_unique_clear_kmer++)
		{
			// Which unsorted kmer is this unique cleartext kmer residing at?
			int cur_clear_kmer_sorted_hap_i = -1;
			for (int i_unsorted_hap = 0; i_unsorted_hap < vecsize(cur_win_cleartext_kmers); i_unsorted_hap++)
			{
				if (compare_kmers(cur_win_cleartext_kmers->at(i_unsorted_hap), cur_win_unique_cleartext_kmers->at(i_unique_clear_kmer)))
				{
					cur_clear_kmer_sorted_hap_i = i_unsorted_hap;
					break;
				}
			} // i_p1_kmer loop.

			if (cur_clear_kmer_sorted_hap_i == -1)
			{
				fprintf(stderr, "Could not find cleartext kmer in unsorted.\n");
				exit(1);
			}

			// For the current unique cleartext kmer, use the corresponding model1 proxized panel and find which unique proxy kmer it maps to
			// Then store the frequency.
			proxy1_cnts[i_unique_clear_kmer] = 0;
			for (int i_p1_kmer = 0; i_p1_kmer < vecsize(cur_win_unique_proxy1_kmers); i_p1_kmer++)
			{
				if (compare_kmers(cur_win_proxy1_kmers->at(cur_clear_kmer_sorted_hap_i), cur_win_unique_proxy1_kmers->at(i_p1_kmer)))
				{
					proxy1_cnts[i_unique_clear_kmer] = cur_win_unique_proxy1_kmer_cnts->at(i_p1_kmer);

					if (__DUMP_REIDENTIFICATION_ATTACK_MSGS__)
					{
						fprintf(stderr, "Found %d. cleartext kmer @ %d. kmer in model 1 proxy\n", i_unique_clear_kmer, i_p1_kmer);
					}

					break;
				}
			} // i_p1_kmer loop.

			if (proxy1_cnts[i_unique_clear_kmer] == 0)
			{
				fprintf(stderr, "Could not find cleartext kmer in proxy1.\n");
				exit(1);
			}

			// For the current unique cleartext kmer, use the corresponding model2 proxized panel and find which unique proxy kmer it maps to
			// Then store the frequency.
			proxy2_cnts[i_unique_clear_kmer] = 0;
			for (int i_p2_kmer = 0; i_p2_kmer < vecsize(cur_win_unique_proxy2_kmers); i_p2_kmer++)
			{
				if (compare_kmers(cur_win_proxy2_kmers->at(cur_clear_kmer_sorted_hap_i), cur_win_unique_proxy2_kmers->at(i_p2_kmer)))
				{
					proxy2_cnts[i_unique_clear_kmer] = cur_win_unique_proxy2_kmer_cnts->at(i_p2_kmer);

					if (__DUMP_REIDENTIFICATION_ATTACK_MSGS__)
					{
						fprintf(stderr, "Found %d. cleartext kmer @ %d. kmer in model 2 proxy\n", i_unique_clear_kmer, i_p2_kmer);
					}
					
					break;
				}
			} // i_p2_kmer loop.

			if (proxy2_cnts[i_unique_clear_kmer] == 0)
			{
				fprintf(stderr, "Could not find cleartext kmer in proxy2.\n");
				exit(1);
			}
		} // i_clear_kmer loop.

		double cur_corr = 0;
		get_correlation(proxy1_cnts, proxy2_cnts, vecsize(cur_win_unique_cleartext_kmers), cur_corr);
		//fprintf(stderr, "%d\t%lf\n", var_i, cur_corr);
		fprintf(f_model1proxy_2_model2proxy_hapfreq_corr_stats, "%d\t%lf\n", var_i, cur_corr);

		////////////////////////////////////////////////////////////////////////////////
	} // i_var loop.

	close_f(f_model1proxy_2_model2proxy_hapfreq_corr_stats, model1_proxy_2_model2_proxy_hapfreq_corr_stats_fp);
} // get_proxy2proxy_mapped_haplotype_frequency_correlation function.

void get_hapfreq_confusion_matrix(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vicinity,
	double min_HF_ref1_freq_2_quantify,
	double max_LF_ref2_freq_2_quantify,
	char* op_fp)
{
	vector<char*>* ref1_sample_ids = buffer_file(ref1_sample_list_fp);
	vector<t_annot_region*>* ref1_var_regs = load_variant_signal_regions_wrapper(ref1_panel_matbed_fp, ref1_sample_list_fp);

	vector<char*>* ref2_sample_ids = buffer_file(ref2_sample_list_fp);
	vector<t_annot_region*>* ref2_var_regs = load_variant_signal_regions_wrapper(ref2_panel_matbed_fp, ref2_sample_list_fp);

	fprintf(stderr, "Loaded %d/%d variants for %d/%d subjects.\n", vecsize(ref1_var_regs), vecsize(ref2_var_regs), vecsize(ref1_sample_ids), vecsize(ref2_sample_ids));

	fprintf(stderr, "Calculating haplotype confusion for min(ref1 freq)=%.4f ; max(ref2 ref)=%.4f over +/-%d var windows.\n", 
		min_HF_ref1_freq_2_quantify, max_LF_ref2_freq_2_quantify, n_vicinity);

	// This is necessary to make sure we are following the haplotypes correctly.
	sort(ref1_var_regs->begin(), ref1_var_regs->end(), sort_regions);
	sort(ref2_var_regs->begin(), ref2_var_regs->end(), sort_regions);

	//double min_HF_ref1_freq_2_quantify = 0.2; // Do not try to match haplotypes lower than this frequency.
	//double max_LF_ref2_freq_2_quantify = 0.05; // Do not try to match haplotypes lower than this frequency.

	int max_geno_ref1 = get_max_genotype_value(ref1_var_regs, ref1_sample_ids);
	int max_geno_ref2 = get_max_genotype_value(ref2_var_regs, ref2_sample_ids);

	if (max_geno_ref1 == 3 &&
		max_geno_ref2 == 3)
	{
		fprintf(stderr, "Both panels are haplocoded\n");
	}
	else
	{
		fprintf(stderr, "One of the panels is likely genocoded (%d, %d), make sure they are both haplocoded (i.e., phased)\n", max_geno_ref1, max_geno_ref2);
		exit(0);
	}

	if (ref1_var_regs->size() != ref2_var_regs->size())
	{
		fprintf(stderr, "Panel variant sizes are not matching (%d vs %d)\n", vecsize(ref1_var_regs), vecsize(ref2_var_regs));

		exit(1);
	}

	for (int i_var = 0; i_var < (int)ref1_var_regs->size(); i_var++)
	{
		if (ref1_var_regs->at(i_var)->start != ref2_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Variants are not matching @ %d. variant: %d vs %d\n",
				i_var,
				ref1_var_regs->at(i_var)->start, ref2_var_regs->at(i_var)->start);

			exit(1);
		}
	} // i_var loop.

	vector<char**>* per_ref1_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref1_var_regs, ref1_sample_ids);
	vector<char**>* per_ref2_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref2_var_regs, ref2_sample_ids);

	//t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	FILE* f_op = open_f(op_fp, "w");

	int* per_ref1_hap_n_HF_2_LF_mapping_haps = new int[ref1_sample_ids->size() * 2];
	memset(per_ref1_hap_n_HF_2_LF_mapping_haps, 0, sizeof(int) * (ref1_sample_ids->size() * 2));

	#define GET_HAP_FREQ_BIN_INDEX(freq,n_bins) (MIN(((int)(freq * n_bins)), n_bins-1))

	//int N_FREQ_BINS = (2 * vecsize(ref1_sample_ids)) / 5;
	int N_FREQ_BINS = 20;
	double*** per_ref1_hap_hap2hap_confusion_matrix = new double**[2 * vecsize(ref1_sample_ids)];

	for (int i_hap = 0; i_hap < vecsize(ref1_sample_ids) * 2; i_hap++)
	{
		per_ref1_hap_hap2hap_confusion_matrix[i_hap] = new double* [N_FREQ_BINS + 2];
		for (int i_freq_bin = 0; i_freq_bin <= N_FREQ_BINS; i_freq_bin++)
		{
			per_ref1_hap_hap2hap_confusion_matrix[i_hap][i_freq_bin] = new double[N_FREQ_BINS + 2];
		} // i_freq_bin loop.
	} // i_hap loop.

	double* ref1_sorted_unique_hap_freq_histogram = new double[2 * vecsize(ref1_sample_ids) + 2];
	memset(ref1_sorted_unique_hap_freq_histogram, 0, sizeof(double) * (2 * vecsize(ref1_sample_ids)));
	double* ref2_sorted_unique_hap_freq_histogram = new double[2 * vecsize(ref2_sample_ids) + 2];
	memset(ref2_sorted_unique_hap_freq_histogram, 0, sizeof(double) * (2 * vecsize(ref2_sample_ids)));

	//// For each subject haplotype, count the number of uniq haplotypes in the freq. bin.
	//double** per_hap_n_uniq_haps_per_freq_hist_bin = new double*[2 * vecsize(ref2_sample_ids) + 2];
	//for (int i_hap = 0; i_hap < N_FREQ_BINS; i_hap++)
	//{
	//	per_hap_n_uniq_haps_per_freq_hist_bin[i_hap] = new double[N_FREQ_BINS];
	//	memset(per_hap_n_uniq_haps_per_freq_hist_bin[i_hap], 0, sizeof(double) * N_FREQ_BINS);
	//} // i_hap loop.

	for (int var_i = 0; var_i < (int)ref1_var_regs->size(); var_i++)
	{
		if (var_i % 100 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vicinity, 0);
		int win_end_i = MIN(var_i + n_vicinity, (int)ref1_var_regs->size() - 1);

		int n_vars_per_win = win_end_i - win_start_i + 1;

		vector<char*>* cur_win_ref1_haps = new vector<char*>();
		vector<char*>* cur_win_ref1_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < vecsize(ref1_sample_ids); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref1_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref1_unsorted_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref1_haps->begin(), cur_win_ref1_haps->end(), sort_haplotypes);

		// Get ref2 haplotypes.
		vector<char*>* cur_win_ref2_haps = new vector<char*>();
		vector<char*>* cur_win_ref2_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < vecsize(ref2_sample_ids); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref2_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref2_unsorted_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref2_haps->begin(), cur_win_ref2_haps->end(), sort_haplotypes);

		////////////////////////////////////////////////////////////////////////////////

		// Count the sorted haplotypes to get frequencies.
		vector<int>* per_unique_ref1_haplotype_cnts = new vector<int>();
		vector<char*>* ref1_unique_haplotypes = new vector<char*>();
		int cur_hap_cnt = 1;
		char* cur_uniq_hap = cur_win_ref1_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref1_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref1_haps->at(i_hap - 1), cur_win_ref1_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
				ref1_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref1_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
			ref1_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check1 = 0;
		for (int i_hap = 0; i_hap < (int)ref1_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check1 += per_unique_ref1_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check1 != (int)ref1_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		vector<vector<int>*>* per_unique_ref1_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref1_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref1_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref1: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref1_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		vector<int>* per_unsorted_ref1_hap_freq = new vector<int>();
		for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
		{
			int cur_hap_cnt = -1;
			for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					if (cur_hap_cnt != -1)
					{
						fprintf(stderr, "Sanity check failed: Ref1 haplotype %d matches to multiple unique haplotypes.\n", i_hap);
						exit(1);
					}
					cur_hap_cnt = per_unique_ref1_haplotype_cnts->at(i_uniq_hap);
				}
			} // i_uniq_hap loop.

			if (cur_hap_cnt == -1)
			{
				fprintf(stderr, "Sanity check failed: Ref1 haplotype %d does not match to any unique haplotypes.\n", i_hap);
				exit(1);
			}
			per_unsorted_ref1_hap_freq->push_back(cur_hap_cnt);
		} // i_hap loop.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the unique ref2 haplotypes.
		vector<int>* per_unique_ref2_haplotype_cnts = new vector<int>();
		vector<char*>* ref2_unique_haplotypes = new vector<char*>();
		cur_hap_cnt = 1;
		cur_uniq_hap = cur_win_ref2_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref2_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref2_haps->at(i_hap - 1), cur_win_ref2_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
				ref2_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref2_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
			ref2_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check2 = 0;
		for (int i_hap = 0; i_hap < (int)ref2_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check2 += per_unique_ref2_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check2 != (int)ref2_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		vector<vector<int>*>* per_unique_ref2_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref2_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref2_unsorted_haps->at(i_hap), ref2_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref2_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref2_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref2: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref2_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		// Following aims at quantifying the probability discordance between haplotypes matching in frequency.
		// Calculate and save the total probability difference:
		int n_HF_ref1_2_LF_ref2_matching_haps = 0;
		int n_unique_ref1_haps_evaluated = 0;
		int n_ref1_subjects_evaluated = 0;
		for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
		{
			// Following is used for directing the discordance estimation to certain positions.
			double cur_uniq_ref1_freq = per_unique_ref1_haplotype_cnts->at(i_uniq1) / ((double)(cur_win_ref1_haps->size()));
			if (cur_uniq_ref1_freq < min_HF_ref1_freq_2_quantify)
			{
				continue;
			}

			n_unique_ref1_haps_evaluated++;
			n_ref1_subjects_evaluated += per_unique_ref1_haplotype_cnts->at(i_uniq1);

			vector<int>* cur_uniq_ref1_hap_subj_i = per_unique_ref1_haplotype_subject_i->at(i_uniq1);

			//double cur_uniq_ref1_hap_discordance_score = 0;
			//double cur_uniq_ref1_hap_rand_discordance_score = 0;

			//vector<int>* cur_rand_i_uniq2_indices = rng->fast_permute_indices(0, (int)per_unique_ref2_haplotype_cnts->size());

			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				// Count the number of individuals that are matching between these haplotypes.
				vector<int>* cur_uniq_ref2_hap_subj_i = per_unique_ref2_haplotype_subject_i->at(i_uniq2);

				double cur_uniq_ref2_freq = per_unique_ref2_haplotype_cnts->at(i_uniq2) / ((double)(cur_win_ref2_haps->size()));
				if (cur_uniq_ref2_freq > max_LF_ref2_freq_2_quantify)
				{
					continue;
				}

				// For the current uniq ref1 and ref2 haplotypes, we check below if there are matching individuals.
				// If there are matching individuals, we need to keep 
				int n_matching = 0;
				for (int i_hap = 0; i_hap < (int)cur_uniq_ref1_hap_subj_i->size(); i_hap++)
				{
					for (int j_hap = 0; j_hap < (int)cur_uniq_ref2_hap_subj_i->size(); j_hap++)
					{
						if (cur_uniq_ref1_hap_subj_i->at(i_hap) == cur_uniq_ref2_hap_subj_i->at(j_hap))
						{
							// Update the current haplotype's count.
							per_ref1_hap_n_HF_2_LF_mapping_haps[cur_uniq_ref1_hap_subj_i->at(i_hap)]++;
							n_matching++;

							int ref1_hapfreq_bin_i = GET_HAP_FREQ_BIN_INDEX(cur_uniq_ref1_freq, N_FREQ_BINS);
							int ref2_hapfreq_bin_i = GET_HAP_FREQ_BIN_INDEX(cur_uniq_ref2_freq, N_FREQ_BINS);
							per_ref1_hap_hap2hap_confusion_matrix[cur_uniq_ref1_hap_subj_i->at(i_hap)][ref1_hapfreq_bin_i][ref2_hapfreq_bin_i]++;

							break;
						}
					} // j_s loop.
				} // i_s loop.

				// Here, the number of matches signifies number of individuals that are mapped between the same haplotypes.
				if (n_matching > 0)
				{
					n_HF_ref1_2_LF_ref2_matching_haps++;
				} // match check.
			} // i_uniq2 loop.

			//delete cur_rand_i_uniq2_indices;

			//getc(stdin);
		} // i_uniq1 loop.

		fprintf(f_op, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
				ref1_var_regs->at(var_i)->chrom, ref1_var_regs->at(var_i)->start,
				n_HF_ref1_2_LF_ref2_matching_haps,
				n_unique_ref1_haps_evaluated,
				n_ref1_subjects_evaluated,
				(int)per_unique_ref1_haplotype_cnts->size(), (int)per_unique_ref2_haplotype_cnts->size());

		vector<int>* sorted_ref1_uniq_counts = new vector<int>();
		sorted_ref1_uniq_counts->insert(sorted_ref1_uniq_counts->end(), per_unique_ref1_haplotype_cnts->begin(), per_unique_ref1_haplotype_cnts->end());
		sort(sorted_ref1_uniq_counts->begin(), sorted_ref1_uniq_counts->end(), [](int a, int b) {return a > b; });

		vector<int>* sorted_ref2_uniq_counts = new vector<int>();
		sorted_ref2_uniq_counts->insert(sorted_ref2_uniq_counts->end(), per_unique_ref2_haplotype_cnts->begin(), per_unique_ref2_haplotype_cnts->end());
		sort(sorted_ref2_uniq_counts->begin(), sorted_ref2_uniq_counts->end(), [](int a, int b) {return a > b; });

		for (int i_hap = 0; i_hap < vecsize(sorted_ref1_uniq_counts); i_hap++)
		{
			ref1_sorted_unique_hap_freq_histogram[i_hap] += (double)sorted_ref1_uniq_counts->at(i_hap) / (2*vecsize(ref1_sample_ids));
		} // i_hap loop.
		
		for (int i_hap = 0; i_hap < vecsize(sorted_ref2_uniq_counts); i_hap++)
		{
			ref2_sorted_unique_hap_freq_histogram[i_hap] += (double)sorted_ref2_uniq_counts->at(i_hap) / (2 * vecsize(ref2_sample_ids));
		} // i_hap loop.


		fprintf(stderr, "Win: %d (%d variants)\n", var_i, n_vars_per_win);
		for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
		{
			fprintf(stderr, "%d, ", sorted_ref1_uniq_counts->at(i_uniq1));
		} // i_uniq1 loop.
		fprintf(stderr, "\n");

		for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
		{
			fprintf(stderr, "%d, ", sorted_ref2_uniq_counts->at(i_uniq2));
		}
		fprintf(stderr, "\n");

		delete sorted_ref1_uniq_counts;
		delete sorted_ref2_uniq_counts;

		if (__DUMP_REIDENTIFICATION_ATTACK_MSGS__)
		{
			// Report the haplotype confusion matrix:
			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "Win: %d (%d variants)\n", var_i, n_vars_per_win);

			fprintf(stderr, "Sorted Ref1 Histogram:\n");
			vector<int>* sorted_ref1_uniq_counts = new vector<int>();
			sorted_ref1_uniq_counts->insert(sorted_ref1_uniq_counts->end(), per_unique_ref1_haplotype_cnts->begin(), per_unique_ref1_haplotype_cnts->end());
			sort(sorted_ref1_uniq_counts->begin(), sorted_ref1_uniq_counts->end());
			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "%d, ", sorted_ref1_uniq_counts->at(i_uniq1));
			} // i_uniq1 loop.
			fprintf(stderr, "\n");

			fprintf(stderr, "Sorted Ref2 Histogram:\n");
			vector<int>* sorted_ref2_uniq_counts = new vector<int>();
			sorted_ref2_uniq_counts->insert(sorted_ref2_uniq_counts->end(), per_unique_ref2_haplotype_cnts->begin(), per_unique_ref2_haplotype_cnts->end());
			sort(sorted_ref2_uniq_counts->begin(), sorted_ref2_uniq_counts->end());
			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				fprintf(stderr, "%d, ", sorted_ref2_uniq_counts->at(i_uniq2));
			}
			fprintf(stderr, "\n");

			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "Ref1: Distribution of %d uniq:\n", per_unique_ref1_haplotype_cnts->at(i_uniq1));

				vector<int>* cur_uniq_ref1_hap_subj_i = per_unique_ref1_haplotype_subject_i->at(i_uniq1);
				for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
				{
					// Count the number of individuals that are matching between these haplotypes.
					vector<int>* cur_uniq_ref2_hap_subj_i = per_unique_ref2_haplotype_subject_i->at(i_uniq2);

					// Count the number of subjects matching.
					int n_matching = 0;
					for (int i_s = 0; i_s < (int)cur_uniq_ref1_hap_subj_i->size(); i_s++)
					{
						for (int j_s = 0; j_s < (int)cur_uniq_ref2_hap_subj_i->size(); j_s++)
						{
							if (cur_uniq_ref1_hap_subj_i->at(i_s) == cur_uniq_ref2_hap_subj_i->at(j_s))
							{
								n_matching++;
								break;
							}
						} // j_s loop.
					} // i_s loop.

					if (n_matching > 0)
					{
						fprintf(stderr, "%d-vs-%d\t%d\t%d\t%d\n", i_uniq1, i_uniq2,
							per_unique_ref1_haplotype_cnts->at(i_uniq1), per_unique_ref2_haplotype_cnts->at(i_uniq2),
							n_matching);
					}
				} // i_uniq2 loop.
			} // i_uniq1 loop.

			fprintf(stderr, "Win: %d (%d variants)\n", var_i, n_vars_per_win);

			fprintf(stderr, "Sorted Ref1 Histogram:\n");
			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "%d, ", sorted_ref1_uniq_counts->at(i_uniq1));
			} // i_uniq1 loop.
			fprintf(stderr, "\n");

			fprintf(stderr, "Sorted Ref2 Histogram:\n");
			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				fprintf(stderr, "%d, ", sorted_ref2_uniq_counts->at(i_uniq2));
			}
			fprintf(stderr, "\n");

			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "------------------------------------\n");

			getc(stdin);
		} // Message dump check.

		// Free memory.
		delete cur_win_ref1_haps;
		delete cur_win_ref1_unsorted_haps;
		delete cur_win_ref2_haps;
		delete cur_win_ref2_unsorted_haps;

		delete per_unique_ref1_haplotype_cnts;
		delete ref1_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref1_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref1_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref1_haplotype_subject_i;

		delete per_unsorted_ref1_hap_freq;

		delete per_unique_ref2_haplotype_cnts;
		delete ref2_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref2_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref2_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref2_haplotype_subject_i;
	} // var_i loop.

	close_f(f_op, op_fp);

	// Write the per haplotype HF-LF mapping statistics.
	FILE* f_n_HF_2_LF_mapping_haps_per_hap = open_f("n_HF_2_LF_mapping_haps_per_hap.txt", "w");
	for (int i_hap = 0; i_hap < (2 * vecsize(ref1_sample_ids)); i_hap++)
	{
		fprintf(f_n_HF_2_LF_mapping_haps_per_hap, "%d\t%d\n", i_hap, per_ref1_hap_n_HF_2_LF_mapping_haps[i_hap]);
	} // i_hap loop.
	close_f(f_n_HF_2_LF_mapping_haps_per_hap, NULL);

	FILE* f_hap2hap_matrix = open_f("hap2hap_confusion_matrix.txt", "w");
	for (int i_hap = 0; i_hap < 2 * vecsize(ref1_sample_ids); i_hap++)
	{
		for (int i_bin = 0; i_bin < N_FREQ_BINS; i_bin++)
		{
			for (int j_bin = 0; j_bin < N_FREQ_BINS; j_bin++)
			{
				fprintf(f_hap2hap_matrix, "%d\t%d\t%d\t%.3f\n", i_hap, i_bin, j_bin, per_ref1_hap_hap2hap_confusion_matrix[i_hap][i_bin][j_bin]);
			} // j_bin loop.
		} // i_bin loop.
	} // i_hap loop.
	close_f(f_hap2hap_matrix, NULL);

	FILE* f_ref1_sorted_hap_hist = open_f("ref1_sorted_hap_hist.txt", "w");
	FILE* f_ref2_sorted_hap_hist = open_f("ref2_sorted_hap_hist.txt", "w");
	for (int i_hap = 0; i_hap < 2 * vecsize(ref1_sample_ids); i_hap++)
	{
		fprintf(f_ref1_sorted_hap_hist, "%d\t%.4f\n", i_hap, ref1_sorted_unique_hap_freq_histogram[i_hap]);
		fprintf(f_ref2_sorted_hap_hist, "%d\t%.4f\n", i_hap, ref2_sorted_unique_hap_freq_histogram[i_hap]);
	} // i_hap loop.

	close_f(f_ref1_sorted_hap_hist, NULL);
	close_f(f_ref2_sorted_hap_hist, NULL);
}

// This should be by far the most effective attack:
// Adversary matches the haplotype frequencies between proxy and a reference panel
// For the matching bins of frequency, generates a mapping
void calculate_haplotype_frequency_discordance_score(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
	char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
	int n_vicinity,
	char* op_fp)
{
	vector<char*>* ref1_sample_ids = buffer_file(ref1_sample_list_fp);
	vector<t_annot_region*>* ref1_var_regs = load_variant_signal_regions_wrapper(ref1_panel_matbed_fp, ref1_sample_list_fp);

	vector<char*>* ref2_sample_ids = buffer_file(ref2_sample_list_fp);
	vector<t_annot_region*>* ref2_var_regs = load_variant_signal_regions_wrapper(ref2_panel_matbed_fp, ref2_sample_list_fp);

	// This is necessary to make sure we are following the haplotypes correctly.
	sort(ref1_var_regs->begin(), ref1_var_regs->end(), sort_regions);
	sort(ref2_var_regs->begin(), ref2_var_regs->end(), sort_regions);

	int max_geno_ref1 = get_max_genotype_value(ref1_var_regs, ref1_sample_ids);
	int max_geno_ref2 = get_max_genotype_value(ref2_var_regs, ref2_sample_ids);

	if (max_geno_ref1 == 3 &&
		max_geno_ref2 == 3)
	{
		fprintf(stderr, "Both panels are haplocoded\n");
	}
	else
	{
		fprintf(stderr, "One of the panels is likely genocoded (%d, %d), make sure they are both haplocoded (i.e., phased)\n", max_geno_ref1, max_geno_ref2);
		exit(0);
	}

	if (ref1_var_regs->size() != ref2_var_regs->size())
	{
		fprintf(stderr, "Panel variant sizes are not matching (%d vs %d)\n", (int)ref1_var_regs->size(), (int)ref2_var_regs->size());

		exit(1);
	}

	for (int i_var = 0; i_var < (int)ref1_var_regs->size(); i_var++)
	{
		if (ref1_var_regs->at(i_var)->start != ref2_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Variants are not matching @ %d. variant: %d vs %d\n",
				i_var,
				ref1_var_regs->at(i_var)->start, ref2_var_regs->at(i_var)->start);

			exit(1);
		}
	} // i_var loop.

	vector<char**>* per_ref1_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref1_var_regs, ref1_sample_ids);
	vector<char**>* per_ref2_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref2_var_regs, ref2_sample_ids);

	double* per_hap_low_freq_ref1_haplotypes = new double[(int)per_ref1_subj_haplotypes->size() + 2];
	memset(per_hap_low_freq_ref1_haplotypes, 0, sizeof(double) * ((int)per_ref1_subj_haplotypes->size() + 2));

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	FILE* f_op = open_f(op_fp, "w");
	for (int var_i = 0; var_i < (int)ref1_var_regs->size(); var_i++)
	{
		if (var_i % 100 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vicinity, 0);
		int win_end_i = MIN(var_i + n_vicinity, (int)ref1_var_regs->size() - 1);

		int n_vars_per_win = win_end_i - win_start_i + 1;

		vector<char*>* cur_win_ref1_haps = new vector<char*>();
		vector<char*>* cur_win_ref1_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref1_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref1_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref1_unsorted_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref1_haps->begin(), cur_win_ref1_haps->end(), sort_haplotypes);

		// Get ref2 haplotypes.
		vector<char*>* cur_win_ref2_haps = new vector<char*>();
		vector<char*>* cur_win_ref2_unsorted_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref2_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref2_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
				cur_win_ref2_unsorted_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref2_haps->begin(), cur_win_ref2_haps->end(), sort_haplotypes);

		////////////////////////////////////////////////////////////////////////////////

		// Count the sorted haplotypes to get frequencies.
		vector<int>* per_unique_ref1_haplotype_cnts = new vector<int>();		
		vector<char*>* ref1_unique_haplotypes = new vector<char*>();
		int cur_hap_cnt = 1;
		char* cur_uniq_hap = cur_win_ref1_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref1_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref1_haps->at(i_hap - 1), cur_win_ref1_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
				ref1_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref1_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref1_haplotype_cnts->push_back(cur_hap_cnt);
			ref1_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check1 = 0;
		for (int i_hap = 0; i_hap < (int)ref1_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check1 += per_unique_ref1_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check1 != (int)ref1_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		vector<vector<int>*>* per_unique_ref1_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref1_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref1_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref1: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref1_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		vector<int>* per_unsorted_ref1_hap_freq = new vector<int>();
		for (int i_hap = 0; i_hap < (int)cur_win_ref1_unsorted_haps->size(); i_hap++)
		{
			int cur_hap_cnt = -1;
			for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq_hap++)
			{
				if (compare_haplotypes(cur_win_ref1_unsorted_haps->at(i_hap), ref1_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					if (cur_hap_cnt != -1)
					{
						fprintf(stderr, "Sanity check failed: Ref1 haplotype %d matches to multiple unique haplotypes.\n", i_hap);
						exit(1);
					}
					cur_hap_cnt = per_unique_ref1_haplotype_cnts->at(i_uniq_hap);
				}
			} // i_uniq_hap loop.

			if (cur_hap_cnt == -1)
			{
				fprintf(stderr, "Sanity check failed: Ref1 haplotype %d does not match to any unique haplotypes.\n", i_hap);
				exit(1);
			}
			per_unsorted_ref1_hap_freq->push_back(cur_hap_cnt);
		} // i_hap loop.

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the unique ref2 haplotypes.
		vector<int>* per_unique_ref2_haplotype_cnts = new vector<int>();
		vector<char*>* ref2_unique_haplotypes = new vector<char*>();
		cur_hap_cnt = 1;
		cur_uniq_hap = cur_win_ref2_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_ref2_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_ref2_haps->at(i_hap - 1), cur_win_ref2_haps->at(i_hap), n_vars_per_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
				ref2_unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_ref2_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_ref2_haplotype_cnts->push_back(cur_hap_cnt);
			ref2_unique_haplotypes->push_back(cur_uniq_hap);
		}

		double hap_count_check2 = 0;
		for (int i_hap = 0; i_hap < (int)ref2_unique_haplotypes->size(); i_hap++)
		{
			hap_count_check2 += per_unique_ref2_haplotype_cnts->at(i_hap);
		}

		if (hap_count_check2 != (int)ref2_sample_ids->size() * 2)
		{
			fprintf(stderr, "# of unique hap counts do not match sample size.\n");
			exit(0);
		}

		vector<vector<int>*>* per_unique_ref2_haplotype_subject_i = new vector<vector<int>*>();
		for (int i_uniq_hap = 0; i_uniq_hap < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq_hap++)
		{
			vector<int>* cur_uniq_hap_subject_i = new vector<int>();
			for (int i_hap = 0; i_hap < (int)cur_win_ref2_unsorted_haps->size(); i_hap++)
			{
				if (compare_haplotypes(cur_win_ref2_unsorted_haps->at(i_hap), ref2_unique_haplotypes->at(i_uniq_hap), n_vars_per_win))
				{
					cur_uniq_hap_subject_i->push_back(i_hap);
				}
			} // i_s loop.
			per_unique_ref2_haplotype_subject_i->push_back(cur_uniq_hap_subject_i);

			if ((int)cur_uniq_hap_subject_i->size() != per_unique_ref2_haplotype_cnts->at(i_uniq_hap))
			{
				fprintf(stderr, "Ref2: Sanity check failed for uniq haplotype %d: %d vs %d.\n", i_uniq_hap,
					(int)cur_uniq_hap_subject_i->size(), per_unique_ref2_haplotype_cnts->at(i_uniq_hap));
				exit(1);
			}
		} // i_uniq_hap loop.

		// Following aims at quantifying the probability discordance between haplotypes matching in frequency.
		// Calculate and save the total probability difference:
		double cur_var_total_discordance_score = 0;
		double cur_var_total_rand_discordance_score = 0;
		int n_unique_ref1_haps_evaluated = 0;
		int n_ref1_subjects_evaluated = 0;
		double max_ref1_freq_2_quantify = 0.1; // Do not try to match haplotypes lower than this frequency.
		for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
		{
			// Following is used for directing the discordance estimation to certain positions.
			double cur_uniq_ref1_freq = per_unique_ref1_haplotype_cnts->at(i_uniq1) / ((double)(cur_win_ref1_haps->size()));
			if (cur_uniq_ref1_freq < max_ref1_freq_2_quantify)
			{
				continue;
			}

			n_unique_ref1_haps_evaluated++;
			n_ref1_subjects_evaluated += per_unique_ref1_haplotype_cnts->at(i_uniq1);

			vector<int>* cur_uniq_ref1_hap_subj_i = per_unique_ref1_haplotype_subject_i->at(i_uniq1);

			double cur_uniq_ref1_hap_discordance_score = 0;
			double cur_uniq_ref1_hap_rand_discordance_score = 0;

			vector<int>* cur_rand_i_uniq2_indices = rng->fast_permute_indices(0, (int)per_unique_ref2_haplotype_cnts->size());

			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				// Count the number of individuals that are matching between these haplotypes.
				vector<int>* cur_uniq_ref2_hap_subj_i = per_unique_ref2_haplotype_subject_i->at(i_uniq2);

				// For the current uniq ref1 and ref2 haplotypes, we check below if there are matching individuals.
				// If there are matching individuals, we need to keep 
				int n_matching = 0;
				for (int i_hap = 0; i_hap < (int)cur_uniq_ref1_hap_subj_i->size(); i_hap++)
				{
					for (int j_hap = 0; j_hap < (int)cur_uniq_ref2_hap_subj_i->size(); j_hap++)
					{
						if (cur_uniq_ref1_hap_subj_i->at(i_hap) == cur_uniq_ref2_hap_subj_i->at(j_hap))
						{
							n_matching++;
							break;
						} 
					} // j_s loop.
				} // i_s loop.

				// Here, the number of matches signifies number of individuals that are mapped between the same haplotypes.
				if(n_matching > 0)
				{
					//double match_freq = (double)n_matching / (ref1_sample_ids->size() * 2.0);
					double match_freq = (double)n_matching / (double)(per_unique_ref1_haplotype_cnts->at(i_uniq1));
					if (match_freq > 1.0)
					{
						fprintf(stderr, "Sanity check failed: The matching freq is above total freq.\n");
						exit(1);
					}

					double weighting_factor = match_freq;
					//double weighting_factor = cur_ref1_freq;

					double cur_ref1_freq = (double)(per_unique_ref1_haplotype_cnts->at(i_uniq1)) / ((int)ref1_sample_ids->size() * 2);
					double cur_ref2_freq = (double)(per_unique_ref2_haplotype_cnts->at(i_uniq2)) / ((int)ref2_sample_ids->size() * 2);

					double cur_ref12_freq_diff = fabs(cur_ref1_freq - cur_ref2_freq);

					// Weight the difference with ref1 hap's frequency.
					cur_uniq_ref1_hap_discordance_score += (weighting_factor * cur_ref12_freq_diff);

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// Now, select a random unique haplotype 2 and assign the random matching score.
					// Explicitly simulate matching of the subjects to ref2 randomly.
					// This way, we are explicitly matching the site1's individuals to site2.
					for (int i_match = 0; i_match < n_matching; i_match++)
					{
						// This score matches what we have above.
						double rand_match_factor = 1.0 / (double)(per_unique_ref1_haplotype_cnts->at(i_uniq1));
						double rand_weighting_factor = rand_match_factor;

						// Get the current randomly selected unique haplotype from ref2 and extract its frequency.
						int rand_i_uniq2 = floor(rng->random_double_ran3() * (int)per_unique_ref2_haplotype_cnts->size());
						rand_i_uniq2 = MIN(rand_i_uniq2, (int)per_unique_ref2_haplotype_cnts->size() - 1);
						double cur_ref2_rand_freq = (double)(per_unique_ref2_haplotype_cnts->at(rand_i_uniq2)) / ((int)ref2_sample_ids->size() * 2);

						double cur_ref12_rand_freq_diff = fabs(cur_ref1_freq - cur_ref2_rand_freq);
						cur_uniq_ref1_hap_rand_discordance_score += (rand_weighting_factor * cur_ref12_rand_freq_diff);
					} // i_match simulation.

					//fprintf(stderr, "%d-vs-%d\t%d\t%d\t%d\t%.4f\t%.4f\n", i_uniq1, i_uniq2,
					//	per_unique_ref1_haplotype_cnts->at(i_uniq1), per_unique_ref2_haplotype_cnts->at(i_uniq2),
					//	n_matching,
					//	cur_uniq_ref1_hap_discordance_score,
					//	cur_uniq_ref1_hap_rand_discordance_score);
				} // match check.
			} // i_uniq2 loop.

			delete cur_rand_i_uniq2_indices;

			// How do we pool the discordance scores?
			cur_var_total_discordance_score += cur_uniq_ref1_hap_discordance_score;
			cur_var_total_rand_discordance_score += cur_uniq_ref1_hap_rand_discordance_score;

			//getc(stdin);
		} // i_uniq1 loop.

		fprintf(f_op, "%s\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\n", 
				ref1_var_regs->at(var_i)->chrom, ref1_var_regs->at(var_i)->start,
				cur_var_total_discordance_score, cur_var_total_rand_discordance_score, 
				n_unique_ref1_haps_evaluated,
				n_ref1_subjects_evaluated,
				(int)per_unique_ref1_haplotype_cnts->size(), (int)per_unique_ref2_haplotype_cnts->size());

		if (__DUMP_REIDENTIFICATION_ATTACK_MSGS__)
		{
			// Report the haplotype confusion matrix:
			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "Win: %d (%d variants)\n", var_i, n_vars_per_win);

			fprintf(stderr, "Sorted Ref1 Histogram:\n");
			vector<int>* sorted_ref1_uniq_counts = new vector<int>();
			sorted_ref1_uniq_counts->insert(sorted_ref1_uniq_counts->end(), per_unique_ref1_haplotype_cnts->begin(), per_unique_ref1_haplotype_cnts->end());
			sort(sorted_ref1_uniq_counts->begin(), sorted_ref1_uniq_counts->end());
			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "%d, ", sorted_ref1_uniq_counts->at(i_uniq1));
			} // i_uniq1 loop.
			fprintf(stderr, "\n");
			
			fprintf(stderr, "Sorted Ref2 Histogram:\n");
			vector<int>* sorted_ref2_uniq_counts = new vector<int>();
			sorted_ref2_uniq_counts->insert(sorted_ref2_uniq_counts->end(), per_unique_ref2_haplotype_cnts->begin(), per_unique_ref2_haplotype_cnts->end());
			sort(sorted_ref2_uniq_counts->begin(), sorted_ref2_uniq_counts->end());
			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				fprintf(stderr, "%d, ", sorted_ref2_uniq_counts->at(i_uniq2));
			}
			fprintf(stderr, "\n");

			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "Ref1: Distribution of %d uniq:\n", per_unique_ref1_haplotype_cnts->at(i_uniq1));

				vector<int>* cur_uniq_ref1_hap_subj_i = per_unique_ref1_haplotype_subject_i->at(i_uniq1);
				for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
				{
					// Count the number of individuals that are matching between these haplotypes.
					vector<int>* cur_uniq_ref2_hap_subj_i = per_unique_ref2_haplotype_subject_i->at(i_uniq2);

					// Count the number of subjects matching.
					int n_matching = 0;
					for (int i_s = 0; i_s < (int)cur_uniq_ref1_hap_subj_i->size(); i_s++)
					{
						for (int j_s = 0; j_s < (int)cur_uniq_ref2_hap_subj_i->size(); j_s++)
						{
							if (cur_uniq_ref1_hap_subj_i->at(i_s) == cur_uniq_ref2_hap_subj_i->at(j_s))
							{
								n_matching++;
								break;
							} 
						} // j_s loop.
					} // i_s loop.

					if(n_matching > 0)
					{ 
						fprintf(stderr, "%d-vs-%d\t%d\t%d\t%d\n", i_uniq1, i_uniq2,
							per_unique_ref1_haplotype_cnts->at(i_uniq1), per_unique_ref2_haplotype_cnts->at(i_uniq2),
							n_matching);
					}
				} // i_uniq2 loop.
			} // i_uniq1 loop.

			fprintf(stderr, "Win: %d (%d variants)\n", var_i, n_vars_per_win);

			fprintf(stderr, "Sorted Ref1 Histogram:\n");
			for (int i_uniq1 = 0; i_uniq1 < (int)per_unique_ref1_haplotype_cnts->size(); i_uniq1++)
			{
				fprintf(stderr, "%d, ", sorted_ref1_uniq_counts->at(i_uniq1));
			} // i_uniq1 loop.
			fprintf(stderr, "\n");

			fprintf(stderr, "Sorted Ref2 Histogram:\n");
			for (int i_uniq2 = 0; i_uniq2 < (int)per_unique_ref2_haplotype_cnts->size(); i_uniq2++)
			{
				fprintf(stderr, "%d, ", sorted_ref2_uniq_counts->at(i_uniq2));
			}
			fprintf(stderr, "\n");

			fprintf(stderr, "------------------------------------\n");
			fprintf(stderr, "------------------------------------\n");

			 getc(stdin);
		} // Message dump check.

		//// Count, for each individual, the number of rare haplotypes they have.
		//for (int i_hap = 0; i_hap < per_ref1_subj_haplotypes->size(); i_hap++)
		//{
		//	if (per_unsorted_ref1_hap_freq->at(i_hap) < 20)
		//	{
		//		per_hap_low_freq_ref1_haplotypes[i_hap]++;
		//	}
		//} // i_s loop.

		// Free memory.
		delete cur_win_ref1_haps;
		delete cur_win_ref1_unsorted_haps;
		delete cur_win_ref2_haps;
		delete cur_win_ref2_unsorted_haps;

		delete per_unique_ref1_haplotype_cnts;
		delete ref1_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref1_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref1_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref1_haplotype_subject_i;

		delete per_unsorted_ref1_hap_freq;

		delete per_unique_ref2_haplotype_cnts;
		delete ref2_unique_haplotypes;
		for (int i_hap = 0; i_hap < (int)per_unique_ref2_haplotype_subject_i->size(); i_hap++)
		{
			delete per_unique_ref2_haplotype_subject_i->at(i_hap);
		} // i_hap loop.
		delete per_unique_ref2_haplotype_subject_i;
	} // var_i loop.

	//for (int i_hap = 0; i_hap < per_ref1_subj_haplotypes->size(); i_hap++)
	//{
	//	fprintf(f_op, "%d\t%.1f\n", i_hap, per_hap_low_freq_ref1_haplotypes[i_hap]);
	//} // i_hap loop.

	close_f(f_op, op_fp);
}

void calculate_proxy2clear_var2var_correlation_stats(char* proxized_panel_matbed, char* proxized_panel_sample_list_fp,
	char* cleartext_panel_matbed, char* cleartext_panel_sample_list_fp)
{
	fprintf(stderr, "Calculating var2var correlations between proxy and clear panels..\n");

	vector<t_annot_region*>* proxy_var_regs = load_variant_signal_regions_wrapper(proxized_panel_matbed, proxized_panel_sample_list_fp);
	sort(proxy_var_regs->begin(), proxy_var_regs->end(), sort_regions);
	vector<t_annot_region*>* orig_var_regs = load_variant_signal_regions_wrapper(cleartext_panel_matbed, cleartext_panel_sample_list_fp);
	sort(orig_var_regs->begin(), orig_var_regs->end(), sort_regions);

	vector<char*>* proxy_sample_ids = buffer_file(proxized_panel_sample_list_fp);
	vector<char*>* orig_sample_ids = buffer_file(cleartext_panel_sample_list_fp);

	if (proxy_sample_ids->size() != orig_sample_ids->size())
	{
		fprintf(stderr, "Sample sizes are not matching: %d, %d\n",
				(int)proxy_sample_ids->size(), (int)orig_sample_ids->size());
		exit(0);
	}

	FILE* f_per_var2var_corrs = open_f("var2var_correlations.bed", "w");
	double* proxy_geno_sig = new double[(int)proxy_sample_ids->size() + 2];
	double* orig_geno_sig = new double[(int)proxy_sample_ids->size() + 2];
	for (int i_var = 0; i_var < (int)proxy_var_regs->size(); i_var++)
	{
		if (proxy_var_regs->at(i_var)->start != orig_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: Coordinates are not matching between original and proxy..\n");
			exit(0);
		}

		void** proxy_var_info = (void**)(proxy_var_regs->at(i_var)->data);
		char* cur_var_proxy_geno_sig = (char*)(proxy_var_info[0]);
		void** orig_var_info = (void**)(orig_var_regs->at(i_var)->data);
		char* cur_var_orig_geno_sig = (char*)(orig_var_info[0]);

		for (int i_s = 0; i_s < (int)proxy_sample_ids->size(); i_s++)
		{
			proxy_geno_sig[i_s] = get_genotype_per_haplocoded_genotype(cur_var_proxy_geno_sig[i_s]);
			orig_geno_sig[i_s] = get_genotype_per_haplocoded_genotype(cur_var_orig_geno_sig[i_s]);
		}// i_s loop.

		double cur_var_corr = 0;
		get_correlation(proxy_geno_sig, orig_geno_sig, (int)proxy_sample_ids->size(), cur_var_corr);
		fprintf(f_per_var2var_corrs, "%s\t%d\t%.4f\n", proxy_var_regs->at(i_var)->chrom, proxy_var_regs->at(i_var)->start, cur_var_corr);
	} // i_var loop.
	close_f(f_per_var2var_corrs, NULL);
}

void calculate_proxy2clear_pairwise_distance_stats(char* proxized_panel_matbed, char* proxized_panel_sample_list_fp,
	char* cleartext_panel_matbed, char* cleartext_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	int l_cleartext_windowizing_win,
	int l_proxy_windowizing_win,
	int l_corr_win,
	char* op_prefix)
{
	vector<t_annot_region*>* proxized_panel_geno_var_regs = load_variant_signal_regions_wrapper(proxized_panel_matbed, proxized_panel_sample_list_fp);
	vector<char*>* proxized_proxy_panel_sample_ids = buffer_file(proxized_panel_sample_list_fp);
	sort(proxized_panel_geno_var_regs->begin(), proxized_panel_geno_var_regs->end(), sort_regions);

	int max_proxy_geno = get_max_genotype_value(proxized_panel_geno_var_regs, proxized_proxy_panel_sample_ids);
	if (max_proxy_geno == 3)
	{
		fprintf(stderr, "Proxized panel is haplocoded..\n");
	}
	else if (max_proxy_geno == 2)
	{
		fprintf(stderr, "Proxized panel is genocoded..\n");
	}
	else
	{
		fprintf(stderr, "Sanity check failed: Could not map coding in proxy panel.\n");
		exit(0);
	}

	vector<t_annot_region*>* cleartext_panel_geno_var_regs = load_variant_signal_regions_wrapper(cleartext_panel_matbed, cleartext_panel_sample_list_fp);
	vector<char*>* cleartext_panel_sample_ids = buffer_file(cleartext_panel_sample_list_fp);
	sort(cleartext_panel_geno_var_regs->begin(), cleartext_panel_geno_var_regs->end(), sort_regions);

	int max_cleartext_geno = get_max_genotype_value(proxized_panel_geno_var_regs, proxized_proxy_panel_sample_ids);
	if (max_cleartext_geno == 3)
	{
		fprintf(stderr, "Cleartext panel is haplocoded..\n");
	}
	else if (max_cleartext_geno == 2)
	{
		fprintf(stderr, "Cleartext panel is genocoded..\n");
	}
	else
	{
		fprintf(stderr, "Sanity check failed: Could not map coding in cleartext panel.\n");
		exit(0);
	}



	double** per_cleartext_windowized_signal = new double* [(int)cleartext_panel_sample_ids->size() + 2];
	for (int j_s = 0; j_s < (int)cleartext_panel_sample_ids->size(); j_s++)
	{
		per_cleartext_windowized_signal[j_s] = new double[(int)cleartext_panel_geno_var_regs->size() + 2];
		memset(per_cleartext_windowized_signal[j_s], 0, sizeof(double) * ((int)cleartext_panel_geno_var_regs->size() + 2));

		// Calculate the vicinity statistic for the current sample.
		for (int mid_i_var = 0; mid_i_var < (int)cleartext_panel_geno_var_regs->size(); mid_i_var++)
		{
			int win_start_i = MAX(0, mid_i_var - l_cleartext_windowizing_win);
			int win_end_i = MIN((int)cleartext_panel_geno_var_regs->size() - 1, mid_i_var + l_cleartext_windowizing_win);

			double cur_var_sig = 0;
			for (int var_i = win_start_i; var_i <= win_end_i; var_i++)
			{
				void** cur_var_info = (void**)(cleartext_panel_geno_var_regs->at(var_i)->data);
				char* cur_var_geno_sig = (char*)(cur_var_info[0]);

				int cur_geno = (int)(cur_var_geno_sig[j_s]);
				if (max_cleartext_geno == 3)
				{
					cur_geno = get_genotype_per_haplocoded_genotype(cur_geno);
				}

				// This can be changed to other summarization options to change attack statistic.
				cur_var_sig += cur_geno;
			} // var_i loop.

			// Set the windowized value.
			per_cleartext_windowized_signal[j_s][mid_i_var] = cur_var_sig;
		} // i_var loop.
	} // j_s loop.

	double*** per_hap_cleartext_windowized_signal = new double** [2];
	for (int j_hap = 0; j_hap < 2; j_hap++)
	{
		per_hap_cleartext_windowized_signal[j_hap] = new double* [(int)cleartext_panel_sample_ids->size() + 2];
		for (int j_s = 0; j_s < (int)cleartext_panel_sample_ids->size(); j_s++)
		{
			per_hap_cleartext_windowized_signal[j_hap][j_s] = new double[(int)cleartext_panel_geno_var_regs->size() + 2];
			memset(per_hap_cleartext_windowized_signal[j_hap][j_s], 0, sizeof(double) * ((int)cleartext_panel_geno_var_regs->size() + 2));

			// Calculate the vicinity statistic for the current sample.
			for (int mid_i_var = 0; mid_i_var < (int)cleartext_panel_geno_var_regs->size(); mid_i_var++)
			{
				int win_start_i = MAX(0, mid_i_var - l_cleartext_windowizing_win);
				int win_end_i = MIN((int)cleartext_panel_geno_var_regs->size() - 1, mid_i_var + l_cleartext_windowizing_win);

				double cur_var_sig = 0;
				for (int var_i = win_start_i; var_i <= win_end_i; var_i++)
				{
					void** cur_var_info = (void**)(cleartext_panel_geno_var_regs->at(var_i)->data);
					char* cur_var_geno_sig = (char*)(cur_var_info[0]);

					int cur_allele = get_allele_per_haplotype(cur_var_geno_sig[j_s], j_hap);

					// This can be changed to other summarization options to change attack statistic.
					cur_var_sig += cur_allele;
				} // var_i loop.

				// Set the windowized value.
				per_hap_cleartext_windowized_signal[j_hap][j_s][mid_i_var] = cur_var_sig;
			} // i_var loop.
		} // j_hap loop.
	} // j_s loop.


	///////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Windowizing the proxy panel.\n");
	// Proxize the windowized signal for proxy panel.
	double** per_proxy_windowized_signal = new double* [(int)proxized_proxy_panel_sample_ids->size() + 2];
	for (int i_s = 0; i_s < (int)proxized_proxy_panel_sample_ids->size(); i_s++)
	{
		per_proxy_windowized_signal[i_s] = new double[(int)proxized_panel_geno_var_regs->size() + 2];
		memset(per_proxy_windowized_signal[i_s], 0, sizeof(double) * ((int)proxized_panel_geno_var_regs->size() + 2));

		// Calculate the vicinity statistic for the current sample.
		for (int mid_i_var = 0; mid_i_var < (int)proxized_panel_geno_var_regs->size(); mid_i_var++)
		{
			int win_start_i = MAX(0, mid_i_var - l_proxy_windowizing_win);
			int win_end_i = MIN((int)proxized_panel_geno_var_regs->size() - 1, mid_i_var + l_proxy_windowizing_win);

			double cur_var_sig = 0;
			for (int var_i = win_start_i; var_i <= win_end_i; var_i++)
			{
				void** cur_var_info = (void**)(proxized_panel_geno_var_regs->at(var_i)->data);
				char* cur_var_geno_sig = (char*)(cur_var_info[0]);

				int cur_geno = (int)(cur_var_geno_sig[i_s]);
				if (max_cleartext_geno == 3)
				{
					cur_geno = get_genotype_per_haplocoded_genotype(cur_geno);
				}

				// This can be changed to other summarization options to change attack statistic.
				cur_var_sig += cur_geno;
			} // var_i loop.

			// Set the windowized value.
			per_proxy_windowized_signal[i_s][mid_i_var] = cur_var_sig;
		} // i_var loop.
	} // i_s loop.

	double*** per_hap_proxy_windowized_signal = new double** [2];
	for (int i_hap = 0; i_hap < 2; i_hap++)
	{
		per_hap_proxy_windowized_signal[i_hap] = new double* [(int)proxized_proxy_panel_sample_ids->size() + 2];
		for (int i_s = 0; i_s < (int)proxized_proxy_panel_sample_ids->size(); i_s++)
		{
			per_hap_proxy_windowized_signal[i_hap][i_s] = new double[(int)proxized_panel_geno_var_regs->size() + 2];
			memset(per_hap_proxy_windowized_signal[i_hap][i_s], 0, sizeof(double) * ((int)proxized_panel_geno_var_regs->size() + 2));

			// Calculate the vicinity statistic for the current sample.
			for (int mid_i_var = 0; mid_i_var < (int)proxized_panel_geno_var_regs->size(); mid_i_var++)
			{
				int win_start_i = MAX(0, mid_i_var - l_proxy_windowizing_win);
				int win_end_i = MIN((int)proxized_panel_geno_var_regs->size() - 1, mid_i_var + l_proxy_windowizing_win);

				double cur_var_sig = 0;
				for (int var_i = win_start_i; var_i <= win_end_i; var_i++)
				{
					void** cur_var_info = (void**)(proxized_panel_geno_var_regs->at(var_i)->data);
					char* cur_var_geno_sig = (char*)(cur_var_info[0]);

					int cur_allele = get_allele_per_haplotype(cur_var_geno_sig[i_s], i_hap);

					// This can be changed to other summarization options to change attack statistic.
					cur_var_sig += cur_allele;
				} // var_i loop.

				// Set the windowized value.
				per_hap_proxy_windowized_signal[i_hap][i_s][mid_i_var] = cur_var_sig;
			} // i_var loop.
		} // i_hap loop.
	} // i_s loop.

	// End of windowizing the proxy panel.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (l_corr_win <= 0)
	{
		l_corr_win = vecsize(proxized_panel_geno_var_regs);
	}

	fprintf(stderr, "Correlating genotype and haplotype level..\n");
	// Correlate all the statistics.
	char corr_stats_fp[1000];
	sprintf(corr_stats_fp, "%s_corr_stats.txt", op_prefix);
	FILE* f_corr_stats = open_f(corr_stats_fp, "w");

	char per_hap_corr_stats_fp[1000];
	sprintf(per_hap_corr_stats_fp, "%s_per_hap_corr_stats.txt", op_prefix);
	FILE* f_per_hap_corr_stats = open_f(per_hap_corr_stats_fp, "w");

	for (int i_s = 0; i_s < (int)proxized_proxy_panel_sample_ids->size(); i_s++)
	{
		for (int cur_win_start_i = 0;
			cur_win_start_i <= ((int)proxized_panel_geno_var_regs->size() - l_corr_win);
			cur_win_start_i += l_corr_win)
		{
			// Genotype level::
			{
				double max_corr_j_s = 0;
				double max_corr = -1000;
				double second_max_corr = 0;
				for (int j_s = 0; j_s < (int)cleartext_panel_sample_ids->size(); j_s++)
				{
					double cur_ij_corr = 0;
					//get_correlation(cur_proxy_subj_geno_sig, per_cleartext_windowized_signal[j_s], proxized_panel_geno_var_regs->size(), cur_ij_corr);
					double* cur_proxy_subj_geno_sig = per_proxy_windowized_signal[i_s];
					double* cur_win_proxy_geno_sig = &(cur_proxy_subj_geno_sig[cur_win_start_i]);

					double* cur_cleartext_subj_geno_sig = per_cleartext_windowized_signal[j_s];
					double* cur_win_cleartext_geno_sig = &(cur_cleartext_subj_geno_sig[cur_win_start_i]);

					get_correlation(cur_win_proxy_geno_sig, cur_win_cleartext_geno_sig, l_corr_win, cur_ij_corr);

					if (cur_ij_corr > max_corr)
					{
						second_max_corr = max_corr;
						max_corr = cur_ij_corr;
						max_corr_j_s = j_s;
					}
					else if (cur_ij_corr > second_max_corr)
					{
						second_max_corr = cur_ij_corr;
					}
				} // j_s loop.

				fprintf(f_corr_stats, "%s\t%d\t%s\t%.4f\t%.4f\n",
					proxized_proxy_panel_sample_ids->at(i_s),
					cur_win_start_i,
					cleartext_panel_sample_ids->at(max_corr_j_s), max_corr, second_max_corr);
			} // geno-level corrs block.


			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				int max_corr_j_s = 0;
				int max_corr_j_hap = 0;
				double max_corr = -1000;

				for (int j_hap = 0; j_hap < 2; j_hap++)
				{
					for (int j_s = 0; j_s < (int)cleartext_panel_sample_ids->size(); j_s++)
					{
						double cur_ij_corr = 0;
						double* cur_proxy_subj_geno_sig = per_hap_proxy_windowized_signal[i_hap][i_s];
						double* cur_win_proxy_geno_sig = &(cur_proxy_subj_geno_sig[cur_win_start_i]);

						double* cur_cleartext_subj_geno_sig = per_hap_cleartext_windowized_signal[j_hap][j_s];
						double* cur_win_cleartext_geno_sig = &(cur_cleartext_subj_geno_sig[cur_win_start_i]);

						get_correlation(cur_win_proxy_geno_sig, cur_win_cleartext_geno_sig, l_corr_win, cur_ij_corr);

						if (cur_ij_corr > max_corr)
						{
							max_corr_j_hap = j_hap;
							max_corr = cur_ij_corr;
							max_corr_j_s = j_s;
						}
					} // j_s loop.
				} // j_hap loop.

				fprintf(f_per_hap_corr_stats, "%s\t%d\t%d\t%s\t%d\t%.4f\n",
					proxized_proxy_panel_sample_ids->at(i_s),
					i_hap,
					cur_win_start_i,
					cleartext_panel_sample_ids->at(max_corr_j_s), max_corr_j_hap, max_corr);
			} // i_hap loop.
		} // cur_win_start_i loop.
	} // i_s loop.
	close_f(f_per_hap_corr_stats, NULL);
	close_f(f_corr_stats, NULL);
} // calculate_proxy2clear_pairwise_distance_stats

// This function tests different aspects of Homer's t on panels with proxized genotypes.
void calculate_windowed_Homer_t_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	int n_vicinity,
	char* op_prefix)
{
	// Load the panels.
	vector<t_annot_region*>* target_query_geno_var_regs = load_variant_signal_regions_wrapper(target_query_panel_matbed, target_query_panel_sample_list_fp);
	vector<char*>* target_query_sample_ids = buffer_file(target_query_panel_sample_list_fp);

	// Following are the AF databases that are used for the attack.
	vector<t_annot_region*>* reference_geno_var_regs = load_variant_signal_regions_wrapper(ref_AF_panel_matbed, ref_AF_panel_sample_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_AF_panel_sample_list_fp);

	vector<t_annot_region*>* database_geno_var_regs = load_variant_signal_regions_wrapper(target_AF_database_panel_matbed, target_AF_database_panel_sample_list_fp);
	vector<char*>* database_sample_ids = buffer_file(target_AF_database_panel_sample_list_fp);

	// Do a global coordinate check on all variants.
	for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
	{
		if (target_query_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start ||
			target_query_geno_var_regs->at(i_var)->start != database_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: Coordinates are not matching among query, ref, mixture database.\n");
			exit(0);
		}
	} // i_var loop.

	// Assign the allele frequencies to all variant regions; these are assumed to be proxy values.
	assign_variant_AFs_to_var_regs(reference_geno_var_regs, ref_sample_ids);
	assign_variant_AFs_to_var_regs(database_geno_var_regs, database_sample_ids);

	double** per_query_sample_windowized_signal = new double* [(int)target_query_sample_ids->size() + 2];
	for (int query_i_s = 0; query_i_s < (int)target_query_sample_ids->size(); query_i_s++)
	{
		per_query_sample_windowized_signal[query_i_s] = new double[(int)target_query_geno_var_regs->size() + 2];
		memset(per_query_sample_windowized_signal[query_i_s], 0, sizeof(double) * ((int)target_query_geno_var_regs->size() + 2));

		// Calculate the vicinity statistic for the current sample.
		for (int mid_i_var = 0; mid_i_var < (int)target_query_geno_var_regs->size(); mid_i_var++)
		{
			int win_start_i = MAX(0, mid_i_var - n_vicinity);
			int win_end_i = MIN((int)target_query_geno_var_regs->size() - 1, mid_i_var + n_vicinity);

			double cur_var_sig = 0;
			for (int var_i = win_start_i; var_i <= win_end_i; var_i++)
			{
				void** cur_var_info = (void**)(target_query_geno_var_regs->at(var_i)->data);
				char* cur_var_geno_sig = (char*)(cur_var_info[0]);

				double cur_geno = (double)(get_genotype_per_haplocoded_genotype(cur_var_geno_sig[query_i_s]));

				// This can be changed to other summarization options to change attack statistic.
				cur_var_sig += cur_geno;
			} // var_i loop.

			// Set the windowized value as the average over the current window.
			per_query_sample_windowized_signal[query_i_s][mid_i_var] = cur_var_sig / (win_end_i - win_start_i + 1);

			if (per_query_sample_windowized_signal[query_i_s][mid_i_var] > 2 ||
				per_query_sample_windowized_signal[query_i_s][mid_i_var] < 0)
			{
				fprintf(stderr, "Sanity check failed: Invalid genotype.\n");
				exit(0);
			}
		} // i_var loop.
	} // j_s loop.

	// Calculate the correlation of each variant to the t statistic: |query-db|-|query-ref| over the individuals?
	// The main question here is that the adversary has access to the original genotypes of the target individuals (query)
	double* per_query_t_stats = new double[(int)target_query_sample_ids->size() + 2];
	char per_query_t_stats_fp[1000];
	sprintf(per_query_t_stats_fp, "%s_per_query_t_stats.txt", op_prefix);
	FILE* f_per_query_t_stats = open_f(per_query_t_stats_fp, "w");
	for (int query_i_s = 0; query_i_s < (int)target_query_sample_ids->size(); query_i_s++)
	{
		// Calculate the t statistic for this sample: Calculate the AFs for all variants.

		// This attack assumes that adversary sequenced the genotypes of the query individuals in original form.
		double* per_var_Dyij = new double[(int)target_query_geno_var_regs->size() + 1];
		for (int mid_i_var = 0; mid_i_var < (int)target_query_geno_var_regs->size(); mid_i_var++)
		{
			//void** cur_var_info = (void**)(target_query_geno_var_regs->at(mid_i_var)->data);
			//char* cur_var_original_geno_sig = (char*)(cur_var_info[0]);
			double cur_query_sample_geno_dbl = per_query_sample_windowized_signal[query_i_s][mid_i_var];

			////double cur_query_sample_geno_dbl = (double)(get_genotype_per_haplocoded_genotype(cur_var_original_geno_sig[query_i_s]));

			//if (cur_query_sample_geno_dbl != 0 &&
			//	cur_query_sample_geno_dbl != 1 &&
			//	cur_query_sample_geno_dbl != 2)
			//{
			//	fprintf(stderr, "Sanity check failed: Invalid genotype!\n");
			//	exit(0);
			//}

			double cur_Dyij = fabs(cur_query_sample_geno_dbl / 2 - reference_geno_var_regs->at(mid_i_var)->dbl_score) -
				fabs(cur_query_sample_geno_dbl / 2 - database_geno_var_regs->at(mid_i_var)->dbl_score);

			per_var_Dyij[mid_i_var] = cur_Dyij;
			//cur_sample_t_stat += cur_var_t_cont;
		} // i_var loop.

		double mean_Dyi = 0;
		double var_Dyi = 0;
		get_stats(per_var_Dyij, (int)target_query_geno_var_regs->size(), mean_Dyi, var_Dyi);
		//double SD_Dyi = pow(var_Dyi, .5);
		//double sqrt_var_num = pow((double)(target_query_geno_var_regs->size()), .5);
		//double cur_sample_t_stat = mean_Dyi / (SD_Dyi / sqrt_var_num);

		//per_query_t_stats[query_i_s] = cur_sample_t_stat;
		per_query_t_stats[query_i_s] = mean_Dyi;

		fprintf(f_per_query_t_stats, "%s\t%.5f\n", target_query_sample_ids->at(query_i_s), mean_Dyi);
	} // query_i_s loop.
	close_f(f_per_query_t_stats, NULL);
} // calculate_Homer_t_statistics_on_proxized_panels function.

// This function tests different aspects of Homer's t on panels with proxized genotypes.
void calculate_Homer_t_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	char* op_prefix)
{
	// Load the panels.
	vector<t_annot_region*>* target_query_geno_var_regs = load_variant_signal_regions_wrapper(target_query_panel_matbed, target_query_panel_sample_list_fp);
	vector<char*>* target_query_sample_ids = buffer_file(target_query_panel_sample_list_fp);
	sort(target_query_geno_var_regs->begin(), target_query_geno_var_regs->end(), sort_regions);

	// Following are the AF databases that are used for the attack.
	vector<t_annot_region*>* reference_geno_var_regs = load_variant_signal_regions_wrapper(ref_AF_panel_matbed, ref_AF_panel_sample_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_AF_panel_sample_list_fp);
	sort(reference_geno_var_regs->begin(), reference_geno_var_regs->end(), sort_regions);

	vector<t_annot_region*>* database_geno_var_regs = load_variant_signal_regions_wrapper(target_AF_database_panel_matbed, target_AF_database_panel_sample_list_fp);
	vector<char*>* database_sample_ids = buffer_file(target_AF_database_panel_sample_list_fp);
	sort(database_geno_var_regs->begin(), database_geno_var_regs->end(), sort_regions);

	// Do a global coordinate check on all variants.
	for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
	{
		if (target_query_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start ||
			target_query_geno_var_regs->at(i_var)->start != database_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: Coordinates are not matching among query, ref, mixture database.\n");
			exit(0);
		}
	} // i_var loop.

	// Assign the allele frequencies to all variant regions.
	assign_variant_AFs_to_var_regs(reference_geno_var_regs, ref_sample_ids);
	assign_variant_AFs_to_var_regs(database_geno_var_regs, database_sample_ids);

	// Calculate the correlation of each variant to the t statistic: |query-db|-|query-ref| over the individuals?
	// The main question here is that the adversary has access to the original genotypes of the target individuals (query)
	double* per_query_t_stats = new double[(int)target_query_sample_ids->size() + 2];
	char per_query_t_stats_fp[1000];
	sprintf(per_query_t_stats_fp, "%s_per_query_t_stats.txt", op_prefix);
	FILE* f_per_query_t_stats = open_f(per_query_t_stats_fp, "w");

	char per_var_per_sample_t_stat_Dyij_fp[1000];
	sprintf(per_var_per_sample_t_stat_Dyij_fp, "%s_per_query_per_var_t_stat_Dyij.txt.gz", op_prefix);
	FILE* f_per_query_per_var_t_stat_Dyij = open_f(per_var_per_sample_t_stat_Dyij_fp, "w");

	for (int query_i_s = 0; query_i_s < (int)target_query_sample_ids->size(); query_i_s++)
	{
		// Calculate the t statistic for this sample: Calculate the AFs for all variants.

		// This attack assumes that adversary sequenced the genotypes of the query individuals in original form.
		double* per_var_Dyij = new double[(int)target_query_geno_var_regs->size() + 1];
		for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
		{
			void** cur_var_info = (void**)(target_query_geno_var_regs->at(i_var)->data);
			char* cur_var_original_geno_sig = (char*)(cur_var_info[0]);

			double cur_query_sample_geno_dbl = (double)(get_genotype_per_haplocoded_genotype(cur_var_original_geno_sig[query_i_s]));

			if (cur_query_sample_geno_dbl != 0 &&
				cur_query_sample_geno_dbl != 1 &&
				cur_query_sample_geno_dbl != 2)
			{
				fprintf(stderr, "Sanity check failed: Invalid genotype!\n");
				exit(0);
			}

			double cur_Dyij = fabs(cur_query_sample_geno_dbl / 2 - reference_geno_var_regs->at(i_var)->dbl_score) -
								fabs(cur_query_sample_geno_dbl / 2 - database_geno_var_regs->at(i_var)->dbl_score);

			per_var_Dyij[i_var] = cur_Dyij;

			// This is saved for pooling and summarizing later.
			fprintf(f_per_query_per_var_t_stat_Dyij, "%d\t%d\t%.5f\n", query_i_s, i_var, cur_Dyij);
		} // i_var loop.

		// t-statistic described in https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000167
		double mean_Dyi = 0;
		double var_Dyi = 0;
		get_stats(per_var_Dyij, (int)target_query_geno_var_regs->size(), mean_Dyi, var_Dyi);
		double SD_Dyi = pow(var_Dyi, .5);
		double sqrt_var_num = pow((double)((int)target_query_geno_var_regs->size()), .5);
		double cur_sample_t_stat = mean_Dyi / (SD_Dyi / sqrt_var_num);

		per_query_t_stats[query_i_s] = cur_sample_t_stat;

		fprintf(f_per_query_t_stats, "%s\t%.5f\n", target_query_sample_ids->at(query_i_s), cur_sample_t_stat);
	} // query_i_s loop.

	// Close summary files.
	close_f(f_per_query_t_stats, NULL);
	close_f(f_per_query_per_var_t_stat_Dyij, per_var_per_sample_t_stat_Dyij_fp);

} // calculate_Homer_t_statistics_on_proxized_panels function.

void pool_summarize_Homer_t_statistics_per_query(char* per_query_per_variant_files_list_fp, char* sample_ids_list_fp, char* summarized_results_op_fp)
{
	vector<char*>* per_query_variant_files = buffer_file(per_query_per_variant_files_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	fprintf(stderr, "Pooling %d query Dyij statistics files for %d subjects.\n", 
		(int)per_query_variant_files->size(), (int)sample_ids->size());

	vector<double>** per_sample_Dyij_stats = new vector<double>*[(int)sample_ids->size() + 2];
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		per_sample_Dyij_stats[i_s] = new vector<double>();
	} // i_s loop.

	// Start reading.
	for (int i_f = 0; i_f < (int)per_query_variant_files->size(); i_f++)
	{
		fprintf(stderr, "Processing Dyij statistics on %s\n", per_query_variant_files->at(i_f));
		FILE* f_Dyij_stats = open_f(per_query_variant_files->at(i_f), "r");
		while (1)
		{
			char* cur_line = getline(f_Dyij_stats);
			if (cur_line == NULL)
			{
				break;
			}

			int cur_i_s = 0;
			int cur_i_var = 0;
			double cur_Dyij = 0.0;
			if (sscanf(cur_line, "%d %d %lf", &cur_i_s, &cur_i_var, &cur_Dyij) != 3)
			{
				fprintf(stderr, "Could not parse the line: %s\n", cur_line);
				exit(1);
			}

			per_sample_Dyij_stats[cur_i_s]->push_back(cur_Dyij);
		} // file reading loop.
		close_f(f_Dyij_stats, per_query_variant_files->at(i_f));
	} // i_f loop.

	fprintf(stderr, "Loaded %d variants for %d subjects.\n", (int)per_sample_Dyij_stats[0]->size(), (int)sample_ids->size());

	FILE* f_per_query_t_stats = open_f(summarized_results_op_fp, "w");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		int n_vars = (int)per_sample_Dyij_stats[i_s]->size();
		double mean_Dyi = 0;
		double SD_Dyi = 0.0;
		get_stats(per_sample_Dyij_stats[i_s], mean_Dyi, SD_Dyi);
		double sqrt_var_num = pow((double)(n_vars), .5);
		double cur_sample_t_stat = mean_Dyi / (SD_Dyi / sqrt_var_num);

		fprintf(f_per_query_t_stats, "%s\t%.5f\t%d\n", sample_ids->at(i_s), cur_sample_t_stat, n_vars);
	} // i_s loop.
	close_f(f_per_query_t_stats, summarized_results_op_fp);
}

void write_securegenome_input_files(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* op_prefix)
{
	vector<t_annot_region*>* target_query_geno_var_regs = load_variant_signal_regions_wrapper(target_query_panel_matbed, target_query_panel_sample_list_fp);
	vector<char*>* target_query_sample_ids = buffer_file(target_query_panel_sample_list_fp);
	sort(target_query_geno_var_regs->begin(), target_query_geno_var_regs->end(), sort_regions);

	// Following are the AF databases that are used for the attack.
	vector<t_annot_region*>* reference_geno_var_regs = load_variant_signal_regions_wrapper(ref_AF_panel_matbed, ref_AF_panel_sample_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_AF_panel_sample_list_fp);
	sort(reference_geno_var_regs->begin(), reference_geno_var_regs->end(), sort_regions);

	//vector<t_annot_region*>* database_geno_var_regs = load_variant_signal_regions_wrapper(target_AF_database_panel_matbed, target_AF_database_panel_sample_list_fp);
	//vector<char*>* database_sample_ids = buffer_file(target_AF_database_panel_sample_list_fp);
	//sort(database_geno_var_regs->begin(), database_geno_var_regs->end(), sort_regions);

	// Do a global coordinate check on all variants.
	for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
	{
		//if (target_query_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start ||
		//	target_query_geno_var_regs->at(i_var)->start != database_geno_var_regs->at(i_var)->start)
		if (target_query_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: Coordinates are not matching among query, ref, mixture database.\n");
			exit(0);
		}
	} // i_var loop.

	fprintf(stderr, "Writing pool genotypes.\n");
	char pool_genotypes_fp[1000];
	sprintf(pool_genotypes_fp, "%s_pool_genotypes.txt", op_prefix);
	FILE* f_pool_genotypes = open_f(pool_genotypes_fp, "w");
	for (int i_s = 0; i_s < (int)target_query_sample_ids->size(); i_s++)
	{
		for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
		{
			void** cur_var_info = (void**)(target_query_geno_var_regs->at(i_var)->data);
			char* cur_var_geno_sig = (char*)(cur_var_info[0]);

			if (i_var > 0)
			{
				fprintf(f_pool_genotypes, "\t");
			}

			int cur_sample_geno = get_genotype_per_haplocoded_genotype(cur_var_geno_sig[i_s]);

			fprintf(f_pool_genotypes, "%d", cur_sample_geno);
		} // i_var loop.
		fprintf(f_pool_genotypes, "\n");
	} // i_s loop.
	close_f(f_pool_genotypes, NULL);

	fprintf(stderr, "Writing reference genotypes.\n");
	char ref_genotypes_fp[1000];
	sprintf(ref_genotypes_fp, "%s_ref_genotypes.txt", op_prefix);
	FILE* f_ref_genotypes = open_f(ref_genotypes_fp, "w");
	for (int i_s = 0; i_s < (int)ref_sample_ids->size(); i_s++)
	{
		for (int i_var = 0; i_var < (int)reference_geno_var_regs->size(); i_var++)
		{
			void** cur_var_info = (void**)(reference_geno_var_regs->at(i_var)->data);
			char* cur_var_geno_sig = (char*)(cur_var_info[0]);

			if (i_var > 0)
			{
				fprintf(f_ref_genotypes, "\t");
			}

			int cur_sample_geno = get_genotype_per_haplocoded_genotype(cur_var_geno_sig[i_s]);

			fprintf(f_ref_genotypes, "%d", cur_sample_geno);
		} // i_var loop.
		fprintf(f_ref_genotypes, "\n");
	} // i_s loop.
	close_f(f_ref_genotypes, NULL);

	fprintf(stderr, "Writing snp ranks.\n");
	char ranks_fp[1000];
	sprintf(ranks_fp, "%s_ranks.txt", op_prefix);
	FILE* f_ranks = open_f(ranks_fp, "w");
	for (int i_var = 0; i_var < (int)reference_geno_var_regs->size(); i_var++)
	{
		fprintf(f_ranks, "0.000001\n");
	} // i_var loop.
	close_f(f_ranks, NULL);
} // write_securegenome_input_files option.

// This function tests different aspects of Homer's t on panels with proxized genotypes.
void calculate_Sankararaman_LRT_statistics_on_proxized_panels(char* target_query_panel_matbed, char* target_query_panel_sample_list_fp,
	char* ref_AF_panel_matbed, char* ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
	char* target_AF_database_panel_matbed, char* target_AF_database_panel_sample_list_fp,
	double maf_cutoff,
	int min_var2var_dist,
	char* op_prefix)
{
	// Load the panels.
	vector<t_annot_region*>* target_query_geno_var_regs = load_variant_signal_regions_wrapper(target_query_panel_matbed, target_query_panel_sample_list_fp);
	vector<char*>* target_query_sample_ids = buffer_file(target_query_panel_sample_list_fp);
	sort(target_query_geno_var_regs->begin(), target_query_geno_var_regs->end(), sort_regions);

	// Following are the AF databases that are used for the attack.
	vector<t_annot_region*>* reference_geno_var_regs = load_variant_signal_regions_wrapper(ref_AF_panel_matbed, ref_AF_panel_sample_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_AF_panel_sample_list_fp);
	sort(reference_geno_var_regs->begin(), reference_geno_var_regs->end(), sort_regions);

	vector<t_annot_region*>* database_geno_var_regs = load_variant_signal_regions_wrapper(target_AF_database_panel_matbed, target_AF_database_panel_sample_list_fp);
	vector<char*>* database_sample_ids = buffer_file(target_AF_database_panel_sample_list_fp);
	sort(database_geno_var_regs->begin(), database_geno_var_regs->end(), sort_regions);

	// Do a global coordinate check on all variants.
	for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
	{
		if (target_query_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start ||
			target_query_geno_var_regs->at(i_var)->start != database_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: Coordinates are not matching among query, ref, mixture database.\n");
			exit(0);
		}
	} // i_var loop.

	// Assign the allele frequencies to all variant regions.
	assign_variant_AFs_to_var_regs(reference_geno_var_regs, ref_sample_ids);
	assign_variant_AFs_to_var_regs(database_geno_var_regs, database_sample_ids);

	// Select the variants: MAF cutoff & distance.
	for (int i_var = 0; i_var < (int)reference_geno_var_regs->size(); i_var++)
	{
		reference_geno_var_regs->at(i_var)->score = 0;
	} // i_var loop.

	t_annot_region* last_selected_var = NULL;
	int n_selected_vars = 0;
	for (int i_var = 0; i_var < (int)reference_geno_var_regs->size(); i_var++)
	{
		if (last_selected_var == NULL)
		{
			if(reference_geno_var_regs->at(i_var)->dbl_score > maf_cutoff &&
				reference_geno_var_regs->at(i_var)->dbl_score < (1.0 - maf_cutoff))
			{
				last_selected_var = reference_geno_var_regs->at(i_var);
				reference_geno_var_regs->at(i_var)->score = 1;
			}
		}
		else
		{
			if (reference_geno_var_regs->at(i_var)->dbl_score > maf_cutoff &&
				reference_geno_var_regs->at(i_var)->dbl_score < (1.0 - maf_cutoff) &&
				get_reg2reg_distance(reference_geno_var_regs->at(i_var), last_selected_var) > min_var2var_dist)
			{
				last_selected_var = reference_geno_var_regs->at(i_var);
				reference_geno_var_regs->at(i_var)->score = 1;
				n_selected_vars++;
			}
		}
	} // i_var loop.
	fprintf(stderr, "Selected %d variants per MAF/distance cutoff.\n", n_selected_vars);

	// Calculate the correlation of each variant to the LRT statistic.
	// The main question here is that the adversary has access to the original genotypes of the target individuals (query)
	//double* per_query_LRT_stats = new double[(int)target_query_sample_ids->size() + 2];
	char per_query_LRT_stats_fp[1000];
	sprintf(per_query_LRT_stats_fp, "%s_per_query_Sankararaman_LRT_stats.txt", op_prefix);
	FILE* f_per_query_LRT_stats = open_f(per_query_LRT_stats_fp, "w");
	for (int query_i_s = 0; query_i_s < (int)target_query_sample_ids->size(); query_i_s++)
	{
		// Calculate the t statistic for this sample: Calculate the AFs for all variants.

		// This attack assumes that adversary sequenced the genotypes of the query individuals in original form.
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			double cur_query_hap_LRT = 0;
			for (int i_var = 0; i_var < (int)target_query_geno_var_regs->size(); i_var++)
			{
				if (reference_geno_var_regs->at(i_var)->score != 1)
				{
					continue;
				}

				if (database_geno_var_regs->at(i_var)->start != reference_geno_var_regs->at(i_var)->start ||
					database_geno_var_regs->at(i_var)->start != target_query_geno_var_regs->at(i_var)->start)
				{
					fprintf(stderr, "Sanity check failed: Coordinate mismatch.\n");
					exit(1);
				}

				void** cur_var_info = (void**)(target_query_geno_var_regs->at(i_var)->data);
				char* cur_var_original_geno_sig = (char*)(cur_var_info[0]);

				int cur_allele = get_allele_per_haplotype(cur_var_original_geno_sig[query_i_s], i_hap);

				double p_hat = database_geno_var_regs->at(i_var)->dbl_score;
				double p = reference_geno_var_regs->at(i_var)->dbl_score;

				// Make sure these are not non-existing variants.
				if (p_hat > 0 && p > maf_cutoff &&
					p_hat < (1.0) && p < (1.0 - maf_cutoff))
				{
					cur_query_hap_LRT += cur_allele * log(p_hat / p) + (1 - cur_allele) * log((1 - p_hat) / (1 - p));
				}
			} // i_var loop.

			//fprintf(f_per_query_LRT_stats, "%s\t%d\t%.5f\n", target_query_sample_ids->at(query_i_s), cur_query_hap_LRT);
			fprintf(f_per_query_LRT_stats, "%d\t%d\t%.5f\t%d\n", query_i_s, i_hap, cur_query_hap_LRT, n_selected_vars);
		} // i_hap loop.
	} // query_i_s loop.

	close_f(f_per_query_LRT_stats, NULL);
} // calculate_Sankararaman_LRT_statistics_on_proxized_panels function.

void pool_summarize_Sankararaman_LRT_statistics_per_query(char* per_query_per_variant_files_list_fp, char* sample_ids_list_fp, char* summarized_results_op_fp)
{
	vector<char*>* per_query_variant_files = buffer_file(per_query_per_variant_files_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	fprintf(stderr, "Pooling %d query LRT statistics files for %d subjects.\n",
		(int)per_query_variant_files->size(), (int)sample_ids->size());

	vector<double>*** per_sample_per_hap_per_file_LRT_stats = new vector<double>**[(int)sample_ids->size() + 2];
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		per_sample_per_hap_per_file_LRT_stats[i_s] = new vector<double>*();
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			per_sample_per_hap_per_file_LRT_stats[i_s][i_hap] = new vector<double>();
		} // i_hap loop.
	} // i_s loop.

	vector<int>* per_file_n_vars = new vector<int>();

	// Start reading.
	for (int i_f = 0; i_f < (int)per_query_variant_files->size(); i_f++)
	{
		fprintf(stderr, "Processing LRT statistics on %s\n", per_query_variant_files->at(i_f));
		FILE* f_Dyij_stats = open_f(per_query_variant_files->at(i_f), "r");
		while (1)
		{
			char* cur_line = getline(f_Dyij_stats);
			if (cur_line == NULL)
			{
				break;
			}

			int cur_i_s = 0;
			int cur_i_hap = 0;
			double cur_LRT_stat = 0.0;
			int cur_n_vars = 0;
			if (sscanf(cur_line, "%d %d %lf %d", &cur_i_s, &cur_i_hap, &cur_LRT_stat, &cur_n_vars) != 4)
			{
				fprintf(stderr, "Could not parse the line: %s\n", cur_line);
				exit(1);
			}

			per_sample_per_hap_per_file_LRT_stats[cur_i_s][cur_i_hap]->push_back(cur_LRT_stat);
			per_file_n_vars->push_back(cur_n_vars);
		} // file reading loop.
		close_f(f_Dyij_stats, per_query_variant_files->at(i_f));
	} // i_f loop.

	fprintf(stderr, "Loaded %d LRT stat summary files for %d subjects.\n", (int)per_sample_per_hap_per_file_LRT_stats[0][0]->size(), (int)sample_ids->size());

	FILE* f_per_query_per_hap_LRT_stats = open_f(summarized_results_op_fp, "w");
	int n_files = (int)per_sample_per_hap_per_file_LRT_stats[0][0]->size();
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			if (n_files != (int)per_sample_per_hap_per_file_LRT_stats[i_s][i_hap]->size())
			{
				fprintf(stderr, "Inconsistent # of files: %d/%d\n", 
					n_files, (int)per_sample_per_hap_per_file_LRT_stats[i_s][i_hap]->size());

				exit(0);
			}
			double total_LRT_stat = 0;
			int n_total_vars = 0;
			for(int i_file = 0; i_file < n_files; i_file++)
			{ 
				total_LRT_stat += per_sample_per_hap_per_file_LRT_stats[i_s][i_hap]->at(i_file);
				n_total_vars += per_file_n_vars->at(i_file);
			} // i_var loop.

			fprintf(f_per_query_per_hap_LRT_stats, "%s\t%.5f\t%d\t%d\n", sample_ids->at(i_s), total_LRT_stat, n_total_vars, n_files);
		} // i_hap loop.
	} // i_s loop.
	close_f(f_per_query_per_hap_LRT_stats, summarized_results_op_fp);
}