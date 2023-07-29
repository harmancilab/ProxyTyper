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

using namespace std;

const bool __DUMP_PROXYTYPER_MSGS__ = false;

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp)
{
	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", vecsize(known_genotype_regs), vecsize(known_sample_ids));

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", vecsize(imputed_sample_ids));

	vector<t_annot_region*>* imputed_genotype_regs = load_variant_signal_regions_wrapper(imputed_genotypes_fp, imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", vecsize(imputed_genotype_regs));
	for (int i_reg = 0; i_reg < vecsize(imputed_genotype_regs); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	  // Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < vecsize(imputed_sample_ids); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < vecsize(known_sample_ids))
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	fprintf(stderr, "Intersecting known regions with imputed regions.\n");
	FILE* f_op = open_f("R2_stats.txt", "w");
	// Intersect and process.
	double** known_imp_sample_geno = new double* [2];
	known_imp_sample_geno[0] = new double[n_matched_samples + 2];
	known_imp_sample_geno[1] = new double[n_matched_samples + 2];
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);
	fprintf(stderr, "Found %d intersections\n", vecsize(intersects));
	int n_processed_imputed_targets = 0;
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "@ %d. intersect         \r", i_int);
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imp_reg = int_info->src_reg;
		t_annot_region* known_reg = int_info->dest_reg;

		void** known_reg_info = (void**)(known_reg->data);
		char* known_reg_geno = (char*)(known_reg_info[0]);

		void** imp_reg_info = (void**)(imp_reg->data);
		char* imp_reg_geno = (char*)(imp_reg_info[0]);

		if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			imp_reg->score == 0)
		{
			t_string_tokens* toks = t_string::tokenize_by_chars(imp_reg->name, "_");

			double AF = atof(toks->at(3)->str());

			imp_reg->score = 1;
			n_processed_imputed_targets++;

			double total_n_non_refs = 0;
			double n_matching_non_refs = 0;
			double n_matching_all = 0;
			double n_all = 0;

			int geno_i = 0;
			for (int imp_i = 0; imp_i < vecsize(imputed_sample_ids); imp_i++)
			{
				if (imp_2_known_sample_i->at(imp_i) < vecsize(known_sample_ids))
				{
					double known_geno = (double)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);
					double imp_geno = (double)(imp_reg_geno[imp_i]);

					known_geno = MAX(0, known_geno);
					imp_geno = MAX(0, imp_geno);

					known_imp_sample_geno[0][geno_i] = known_geno;
					known_imp_sample_geno[1][geno_i] = imp_geno;
					geno_i++;

					// If this known genotype contains a minor allele, use it to quantify MAF allele accuracy.
					bool known_geno_has_minor_allele = false;
					if (AF > 0.5 &&
						known_geno != 2) // If the alternate is the major allele, we will not use the homozygous alternates.
					{
						known_geno_has_minor_allele = true;
					}
					else if (AF < 0.5 &&
						known_geno != 0) // If the alternate is the minor allele, we will not use the homozygous references.
					{
						known_geno_has_minor_allele = true;
					}

					// Update non-ref concordance.
					if (known_geno_has_minor_allele)
					{
						if (known_geno == imp_geno)
						{
							n_matching_non_refs++;
						}

						total_n_non_refs++;
					}

					// Update all concordance.
					if (known_geno == imp_geno)
					{
						n_matching_all++;
					}
					n_all++;
				} // matching check.
			} // imp_i loop.			

			double cur_geno_corr = 0;
			get_correlation(known_imp_sample_geno[0], known_imp_sample_geno[1], geno_i, cur_geno_corr);

			fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t%.0f\t%.0f\t%.0f\t%.0f\n", imp_reg->chrom,
				translate_coord(imp_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(imp_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				imp_reg->name, cur_geno_corr * cur_geno_corr,
				n_matching_all, n_all,
				n_matching_non_refs, total_n_non_refs);
		} // overlapping region name comparison.
	} // i_int loop.
	fclose(f_op);

	fprintf(stderr, "\nDone.\n");
}

bool sort_decreasing(int a, int b)
{
	return(a > b);
}

void assign_variant_AFs_to_var_regs(vector<t_annot_region*>* panel_var_regs, vector<char*>* sample_ids)
{
	int max_geno = get_max_genotype_value(panel_var_regs, sample_ids);

	fprintf(stderr, "Assigning AFs to %d variants (Max geno: %d)\n", (int)panel_var_regs->size(), max_geno);

	// For each resampled subject, calculate the AF mismatch compared to the generating panel.
	// Calculate the allele frequencies of the generating panel.
	//for (int i_var = 0; i_var < (int)panel_var_regs->size(); i_var++)
	for (int i_var = 0; i_var < vecsize(panel_var_regs); i_var++)
	{
		void** var_reg_info = (void**)(panel_var_regs->at(i_var)->data);
		char* geno_sig = (char*)(var_reg_info[0]);

		double total_geno = 0;
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			if (max_geno == 3)
			{
				total_geno += get_genotype_per_haplocoded_genotype(geno_sig[i_s]);
			}
			else if (max_geno == 2)
			{
				total_geno += geno_sig[i_s];
			}
		} // i_s loop.

		panel_var_regs->at(i_var)->dbl_score = (total_geno / ((int)sample_ids->size() * 2));

		if (__DUMP_PROXYTYPER_MSGS__)
		{
			fprintf(stderr, "%s: %.4f\n", panel_var_regs->at(i_var)->name, panel_var_regs->at(i_var)->dbl_score);
		}
	} // i_var loop.
}


bool sort_per_win_hap_dist_nodes(double* node1, double* node2)
{
	return(node1[0] < node2[0]);
}

bool sort_haplotypes(char* hap1, char* hap2)
{
	int i = 0;
	while (hap1[i] != -1)
	{
		if (hap1[i] == hap2[i])
		{
			i++;
		}
		else
		{
			return(hap1[i] < hap2[i]);
		}
	} // i loop.

	return(false);
}

bool compare_haplotypes(char* hap1, char* hap2, int n_vars_per_win)
{
	for(int i = 0; i < n_vars_per_win; i++)
	{
		if (hap1[i] != hap2[i])
		{
			return(false);
		}
	} // i loop.

	return(true);
}

vector<char**>* get_per_subject_haplotypes_per_haplocoded_var_regs(vector<t_annot_region*>* haplocoded_panel_var_regs, vector<char*>* panel_sample_ids)
{
	fprintf(stderr, "Generating per subject haplotypes from %d variants and %d subjects.\n", (int)haplocoded_panel_var_regs->size(), (int)panel_sample_ids->size());

	vector<char**>* per_panel_subj_haplotypes = new vector<char**>[(int)panel_sample_ids->size() + 2];

	// Allocate the per subject haplotypes.
	for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
	{
		char** cur_subj_per_hap_alleles = new char*[2];
		cur_subj_per_hap_alleles[0] = new char[(int)haplocoded_panel_var_regs->size() + 2];
		cur_subj_per_hap_alleles[1] = new char[(int)haplocoded_panel_var_regs->size() + 2];

		per_panel_subj_haplotypes->push_back(cur_subj_per_hap_alleles);
	} // i_s loop.

	// Copy the per subject haplotypes.
	for (int i_var = 0; i_var < (int)haplocoded_panel_var_regs->size(); i_var++)
	{
		void** var_info = (void**)(haplocoded_panel_var_regs->at(i_var)->data);
		char* cur_var_geno = (char*)(var_info[0]);

		for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
		{
			per_panel_subj_haplotypes->at(i_s)[0][i_var] = get_allele_per_haplotype(cur_var_geno[i_s], 0);
			per_panel_subj_haplotypes->at(i_s)[1][i_var] = get_allele_per_haplotype(cur_var_geno[i_s], 1);
		} // i_s loop.
	} // i_var loop.

	for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
	{
		per_panel_subj_haplotypes->at(i_s)[0][(int)haplocoded_panel_var_regs->size()] = -1;
		per_panel_subj_haplotypes->at(i_s)[1][(int)haplocoded_panel_var_regs->size()] = -1;
	}

	return(per_panel_subj_haplotypes);
}

double get_hap_2_hap_distance(char* subj1_win_hap, char* subj2_win_hap, int n_vars_per_win)
{
	double tot_dist = 0;
	for (int var_i = 0; var_i < n_vars_per_win; var_i++)
	{
		double cur_var_diff = (double)(subj1_win_hap[var_i] - subj2_win_hap[var_i]) * (subj1_win_hap[var_i] - subj2_win_hap[var_i]);
		tot_dist += cur_var_diff;
	} // var_i loop.

	return(tot_dist);
}

void get_query_haplotype_frequency_per_reference(char* query_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* op_fp)
{
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* query_var_regs = load_variant_signal_regions_wrapper(query_haplocoded_geno_fp, query_sample_ids_fp);

	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);
	vector<t_annot_region*>* all_ref_var_regs = load_variant_signal_regions_wrapper(ref_haplocoded_geno_fp, ref_sample_ids_fp);

	// This is necessary to make sure we are following the haplotypes correctly.
	sort(query_var_regs->begin(), query_var_regs->end(), sort_regions);
	sort(all_ref_var_regs->begin(), all_ref_var_regs->end(), sort_regions);

	int max_geno_query = get_max_genotype_value(query_var_regs, query_sample_ids);
	int max_geno_ref = get_max_genotype_value(all_ref_var_regs, ref_sample_ids);

	if (max_geno_query == 3 &&
		max_geno_ref == 3)
	{
		fprintf(stderr, "Both panels are haplocoded\n");
	}
	else
	{
		fprintf(stderr, "One of the panels is likely genocoded (%d, %d), make sure they are both haplocoded (i.e., phased)\n", max_geno_query, max_geno_ref);
		exit(0);
	}

	vector<t_annot_region*>* intersects = intersect_annot_regions(query_var_regs, all_ref_var_regs, false, false);
	vector<t_annot_region*>* ref_var_regs = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* ref_var_reg = int_info->dest_reg;
		ref_var_regs->push_back(ref_var_reg);
	} // i_int loop.

	// Sort the overlapping reference regions.
	sort(ref_var_regs->begin(), ref_var_regs->end(), sort_regions);

	if (ref_var_regs->size() != query_var_regs->size())
	{
		fprintf(stderr, "Panel variant sizes are not matching (%d vs %d)\n", (int)ref_var_regs->size(), (int)query_var_regs->size());

		exit(1);
	}

	for (int i_var = 0; i_var < (int)query_var_regs->size(); i_var++)
	{
		if (query_var_regs->at(i_var)->start != ref_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Variants are not matching @ %d. variant: %d vs %d\n",
				i_var,
				query_var_regs->at(i_var)->start, ref_var_regs->at(i_var)->start);

			exit(1);
		}
	} // i_var loop.

	vector<char**>* per_query_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_var_regs, query_sample_ids);
	vector<char**>* per_ref_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_var_regs, ref_sample_ids);

	FILE* f_op = open_f(op_fp, "w");

	//if (var_i % 100 == 0)
	//{
	//	fprintf(stderr, "@ variant %d           \r", var_i);
	//}

	int win_start_i = 0;
	int n_vars_per_win = vecsize(ref_var_regs);
	//int win_end_i = MIN(var_i + n_vars_per_win, (int)ref1_var_regs->size() - 1);

	vector<char*>* cur_win_ref_haps = new vector<char*>();
	for (int i_r_s = 0; i_r_s < (int)ref_sample_ids->size(); i_r_s++)
	{
		for (int i_r_h = 0; i_r_h < 2; i_r_h++)
		{
			cur_win_ref_haps->push_back(per_ref_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
		} // i_r_h loop
	} // i_r_s loop.

	sort(cur_win_ref_haps->begin(), cur_win_ref_haps->end(), sort_haplotypes);

	////////////////////////////////////////////////////////////////////////////////
	// Get query haplotypes.
	vector<char*>* cur_win_query_haps = new vector<char*>();
	vector<char*>* cur_win_unsorted_query_haps = new vector<char*>();
	for (int i_r_s = 0; i_r_s < (int)query_sample_ids->size(); i_r_s++)
	{
		for (int i_r_h = 0; i_r_h < 2; i_r_h++)
		{
			cur_win_query_haps->push_back(per_query_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			cur_win_unsorted_query_haps->push_back(per_query_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
		} // i_r_h loop
	} // i_r_s loop.

	// Don't sort the query haplotypes.
	sort(cur_win_query_haps->begin(), cur_win_query_haps->end(), sort_haplotypes);
	////////////////////////////////////////////////////////////////////////////////

	// Count the sorted haplotypes to get frequencies.
	vector<int>* per_unique_ref_haplotype_cnts = new vector<int>();
	vector<char*>* ref_unique_haplotypes = new vector<char*>();
	int cur_hap_cnt = 1;
	char* cur_uniq_hap = cur_win_ref_haps->at(0);
	for (int i_hap = 1; i_hap < (int)cur_win_ref_haps->size(); i_hap++)
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

			// This is a new unique haplotype. This only gets added when we reach to a new haplotype.
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

	/////////////////////
	// Get the unique query haplotypes; we report these in the final output.
	vector<int>* per_unique_query_haplotype_cnts = new vector<int>();
	vector<char*>* query_unique_haplotypes = new vector<char*>();
	cur_hap_cnt = 1;
	cur_uniq_hap = cur_win_query_haps->at(0);
	for (int i_hap = 1; i_hap < (int)cur_win_query_haps->size(); i_hap++)
	{
		// Is this haplotype the same as the previous one?
		if (compare_haplotypes(cur_win_query_haps->at(i_hap - 1), cur_win_query_haps->at(i_hap), n_vars_per_win))
		{
			cur_hap_cnt++;
		}
		else
		{
			per_unique_query_haplotype_cnts->push_back(cur_hap_cnt);
			query_unique_haplotypes->push_back(cur_uniq_hap);

			cur_hap_cnt = 1;
			cur_uniq_hap = cur_win_query_haps->at(i_hap);
		}
	} // i_hap loop.

	// Add the current haplotype:
	if (cur_hap_cnt > 0)
	{
		per_unique_query_haplotype_cnts->push_back(cur_hap_cnt);
		query_unique_haplotypes->push_back(cur_uniq_hap);
	}

	////////////////////////////////////////////////////////////////////////////////
	//// Start comparing the unique haplotypes.
	//fprintf(stderr, "Var %d: Ref1: %d/%d unique haplotypes; Ref2: %d/%d unique haplotypes.\n", 
	//		var_i, 
	//		ref1_unique_haplotypes->size(), ref1_sample_ids->size(),
	//		ref2_unique_haplotypes->size(), ref2_sample_ids->size());

	////////////////////////////////////////////////////////////////////////////////
	for (int i_hap_q = 0; i_hap_q < (int)cur_win_unsorted_query_haps->size(); i_hap_q++)
	{
		int cur_query_hap_freq = 0;
		for (int i_uhap_r = 0; i_uhap_r < (int)ref_unique_haplotypes->size(); i_uhap_r++)
		{
			if (compare_haplotypes(cur_win_unsorted_query_haps->at(i_hap_q), ref_unique_haplotypes->at(i_uhap_r), n_vars_per_win))
			{
				cur_query_hap_freq = per_unique_ref_haplotype_cnts->at(i_uhap_r);
				break;
			}
		} // i_uhap_r loop.

		int q_i_s = i_hap_q / 2;
		int q_i_hap = i_hap_q % 2;

		fprintf(f_op, "%d\t%d\t%d\t%d\t%d\n", q_i_s, q_i_hap,
				cur_query_hap_freq, vecsize(per_unique_ref_haplotype_cnts), vecsize(query_unique_haplotypes));
	} // i_hap_q loop.

	close_f(f_op, op_fp);
}

void compare_per_win_haplotypes(char* ref1_panel_matbed_fp, char* ref1_sample_list_fp,
								char* ref2_panel_matbed_fp, char* ref2_sample_list_fp,
								int n_vars_per_win,
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

	FILE* f_op = open_f(op_fp, "w");
	for (int var_i = 0; var_i < (int)ref1_var_regs->size(); var_i++)
	{
		if (var_i % 100 == 0)
		{
			fprintf(stderr, "@ variant %d           \r", var_i);
		}

		int win_start_i = MAX(var_i - n_vars_per_win, 0);
		//int win_end_i = MIN(var_i + n_vars_per_win, (int)ref1_var_regs->size() - 1);

		vector<char*>* cur_win_ref1_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref1_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref1_haps->push_back(per_ref1_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_ref1_haps->begin(), cur_win_ref1_haps->end(), sort_haplotypes);

		// Get ref2 haplotypes.
		vector<char*>* cur_win_ref2_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref2_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_ref2_haps->push_back(per_ref2_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
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

		/////////////////////
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

		////////////////////////////////////////////////////////////////////////////////
		//// Start comparing the unique haplotypes.
		//fprintf(stderr, "Var %d: Ref1: %d/%d unique haplotypes; Ref2: %d/%d unique haplotypes.\n", 
		//		var_i, 
		//		ref1_unique_haplotypes->size(), ref1_sample_ids->size(),
		//		ref2_unique_haplotypes->size(), ref2_sample_ids->size());

		////////////////////////////////////////////////////////////////////////////////
		int n_ref1_spec_haps = 0;
		int n_ref1_spec_samples = 0;
		int n_ref2_spec_haps = 0;
		int n_ref2_spec_samples = 0;
		for (int i_hap1 = 0; i_hap1 < (int)ref1_unique_haplotypes->size(); i_hap1++)
		{
			bool found_hap1 = false;
			for (int i_hap2 = 0; i_hap2 < (int)ref2_unique_haplotypes->size(); i_hap2++)
			{
				if (compare_haplotypes(ref1_unique_haplotypes->at(i_hap1), ref2_unique_haplotypes->at(i_hap2), n_vars_per_win))
				{
					found_hap1 = true;
					break;
				}
			} // i_hap2 loop.

			if (!found_hap1)
			{
				n_ref1_spec_haps++;
				n_ref1_spec_samples += per_unique_ref1_haplotype_cnts->at(i_hap1);
			}
		} // i_hap1 loop.

		for (int i_hap2 = 0; i_hap2 < (int)ref2_unique_haplotypes->size(); i_hap2++)
		{
			bool found_hap2 = false;
			for (int i_hap1 = 0; i_hap1 < (int)ref1_unique_haplotypes->size(); i_hap1++)
			{
				if (compare_haplotypes(ref1_unique_haplotypes->at(i_hap1), ref2_unique_haplotypes->at(i_hap2), n_vars_per_win))
				{
					found_hap2 = true;
					break;
				}
			} // i_hap1 loop.

			if (!found_hap2)
			{
				n_ref2_spec_haps++;
				n_ref2_spec_samples += per_unique_ref2_haplotype_cnts->at(i_hap2);
			}
		} // i_hap2 loop.

		fprintf(f_op, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", var_i, 
				n_ref1_spec_haps, n_ref2_spec_haps,
				n_ref1_spec_samples, n_ref2_spec_samples,
			(int)ref1_unique_haplotypes->size(), (int)ref1_sample_ids->size(),
			(int)ref2_unique_haplotypes->size(), (int)ref2_sample_ids->size());

		sort(per_unique_ref1_haplotype_cnts->begin(), per_unique_ref1_haplotype_cnts->end());
		sort(per_unique_ref2_haplotype_cnts->begin(), per_unique_ref2_haplotype_cnts->end());

		fprintf(stderr, "Ref1 sorted histogram:\n");
		for (int i_hap = 0; i_hap < (int)ref1_unique_haplotypes->size(); i_hap++)
		{
			fprintf(stderr, "%d, ", per_unique_ref1_haplotype_cnts->at(i_hap));
		}

		fprintf(stderr, "\nRef2 sorted histogram:\n");
		for (int i_hap = 0; i_hap < (int)ref2_unique_haplotypes->size(); i_hap++)
		{
			fprintf(stderr, "%d, ", per_unique_ref2_haplotype_cnts->at(i_hap));
		}

		getc(stdin);
	} // var_i loop.

	close_f(f_op, op_fp);
}

void get_per_win_haplotype_frequencies_per_reference(char* ref_panel_matbed_fp, char* ref_sample_list_fp, int n_vars_per_win, char* op_fp)
{
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_list_fp);
	vector<t_annot_region*>* ref_var_regs = load_variant_signal_regions_wrapper(ref_panel_matbed_fp, ref_sample_list_fp);

	int max_geno = get_max_genotype_value(ref_var_regs, ref_sample_ids);
	if (max_geno == 3)
	{
		fprintf(stderr, "Reference is likely haplocoded\n");
	}
	else if (max_geno == 2)
	{
		fprintf(stderr, "Reference is likely genocoded\n");
		exit(0);
	}
	else
	{
		fprintf(stderr, "Max geno is %d, not expected.\n", max_geno);
		exit(0);
	}

	vector<char**>* per_ref_subj_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_var_regs, ref_sample_ids);

	FILE* f_per_win_hap_freqs = open_f(op_fp, "w");
	for (int var_i = 0; var_i < (int)ref_var_regs->size(); var_i++)
	{
		//if (var_i % 1000 == 0)
		//{
		//	fprintf(stderr, "@ variant %d           \r", var_i);
		//}

		int win_start_i = MAX(var_i - n_vars_per_win, 0);
		int win_end_i = MIN(var_i + n_vars_per_win, (int)ref_var_regs->size() - 1);

		int n_vars_per_cur_win = win_end_i - win_start_i + 1;

		vector<char*>* cur_win_haps = new vector<char*>();
		for (int i_r_s = 0; i_r_s < (int)ref_sample_ids->size(); i_r_s++)
		{
			for (int i_r_h = 0; i_r_h < 2; i_r_h++)
			{
				cur_win_haps->push_back(per_ref_subj_haplotypes->at(i_r_s)[i_r_h] + win_start_i);
			} // i_r_h loop
		} // i_r_s loop.

		sort(cur_win_haps->begin(), cur_win_haps->end(), sort_haplotypes);

		// Count the sorted haplotypes to get frequencies.
		vector<int>* per_unique_haplotype_cnts = new vector<int>();
		vector<char*>* unique_haplotypes = new vector<char*>();
		int cur_hap_cnt = 1;
		char* cur_uniq_hap = cur_win_haps->at(0);
		for (int i_hap = 1; i_hap < (int)cur_win_haps->size(); i_hap++)
		{
			// Is this haplotype the same as the previous one?
			if (compare_haplotypes(cur_win_haps->at(i_hap - 1), cur_win_haps->at(i_hap), n_vars_per_cur_win))
			{
				cur_hap_cnt++;
			}
			else
			{
				per_unique_haplotype_cnts->push_back(cur_hap_cnt);
				unique_haplotypes->push_back(cur_uniq_hap);

				cur_hap_cnt = 1;
				cur_uniq_hap = cur_win_haps->at(i_hap);
			}
		} // i_hap loop.

		// Add the current haplotype:
		if (cur_hap_cnt > 0)
		{
			per_unique_haplotype_cnts->push_back(cur_hap_cnt);
			unique_haplotypes->push_back(cur_uniq_hap);
		}

		if(var_i % 1000 == 0)
		{
			fprintf(stderr, "Var %d: %d/%d unique haplotypes.\n", var_i, (int)unique_haplotypes->size(), (int)ref_sample_ids->size());
		}
		
		fprintf(f_per_win_hap_freqs, "%d\t%d\t%d\n", MAX(0, var_i - n_vars_per_win), (int)unique_haplotypes->size(), (int)ref_sample_ids->size());
	} // var_i loop.

	close_f(f_per_win_hap_freqs, op_fp);
}

void convert_BEAGLE_imputed_meta_panel_2_matbed(char* vcf_fp, char* sample_ids_list_fp, double imp_geno_rand_weight, char* op_fp)
{
	fprintf(stderr, "Converting BEAGLE imputed genotypes to noisy matbed.\n");

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d subject ids.\n", (int)sample_ids->size());

	char tok_buff[1000];

	char chrom[1000];
	char var_id[1000];

	char ref_all[100];
	char alt_all[100];

	char info_str[1000];

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	vector<t_annot_region*>* geno_sig_regs = new vector<t_annot_region*>();

	FILE* f_geno_rand_stats = open_f("imputed_geno_randomization_stats.txt", "w");

	FILE* f_vcf = open_f(vcf_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete [] cur_line;
			continue;
		}

		char* cur_var_geno_sig = new char[(int)sample_ids->size() + 2];

		// Start parsing:
		// 22      18227077        rs8190313       A       G       .       PASS    DR2=1.00;AF=0.0469;IMP  0|0:0:0:0:1,0,0
		int cur_char_i = 0;
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		strcpy(chrom, tok_buff);
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		int posn = atoi(tok_buff);
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		strcpy(var_id, tok_buff);

		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		strcpy(ref_all, tok_buff);
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		strcpy(alt_all, tok_buff);

		// This is the empty column.
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);

		// PASS string.
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);

		// Get format string.
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);

		// Info string.
		t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i);
		strcpy(info_str, tok_buff);

		double abs_geno_dist = 0;
		double avg_geno_entropy = 0;
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			if (!t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i))
			{
				fprintf(stderr, "Could not read %d. subject's genotypes in variant:\n%s\n", i_s, cur_line);
				exit(0);
			}

			//0|0:0:0:0:1,0,0
			char geno_tok_buff[1000];
			int geno_char_i = 0;
			t_string::get_next_token(tok_buff, geno_tok_buff, 1000, ":", geno_char_i);
			int orig_M_allele = geno_tok_buff[0] - '0';
			int orig_F_allele = geno_tok_buff[2] - '0';
			double orig_imp_geno = orig_M_allele * 2 + orig_F_allele;

			t_string::get_next_token(tok_buff, geno_tok_buff, 1000, ":", geno_char_i);

			t_string::get_next_token(tok_buff, geno_tok_buff, 1000, ":", geno_char_i);
			double allele_M_AAF_prob = atof(geno_tok_buff);
			t_string::get_next_token(tok_buff, geno_tok_buff, 1000, ":", geno_char_i);
			double allele_F_AAF_prob = atof(geno_tok_buff);

			if (allele_F_AAF_prob > 0 &&
				allele_M_AAF_prob > 0)
			{
				avg_geno_entropy += (log(allele_F_AAF_prob)*allele_F_AAF_prob+ log(allele_M_AAF_prob)*allele_M_AAF_prob) / log(2);
			}

			// p(allele)=0.99
			// High randomization weight means we will make randomized distribution more like original calls.
			double scaled_M_AAF_prob = pow(allele_M_AAF_prob, imp_geno_rand_weight);
			double scaled_M_RAF_prob = pow(1 - allele_M_AAF_prob, imp_geno_rand_weight);
			double norm_M_AAF_prob = scaled_M_AAF_prob / (scaled_M_AAF_prob + scaled_M_RAF_prob);

			double scaled_F_AAF_prob = pow(allele_F_AAF_prob, imp_geno_rand_weight);
			double scaled_F_RAF_prob = pow(1 - allele_F_AAF_prob, imp_geno_rand_weight);
			double norm_F_AAF_prob = scaled_F_AAF_prob / (scaled_F_AAF_prob + scaled_F_RAF_prob);

			double M_rand = rng->random_double_ran3();
			char M_allele = (M_rand > norm_M_AAF_prob) ? (0) : (1);
			double F_rand = rng->random_double_ran3();
			char F_allele = (F_rand > norm_F_AAF_prob) ? (0) : (1);
			//char M_allele = (rng->random_double_ran3() > allele_M_AAF_prob) ? (0) : (1);
			//char F_allele = (rng->random_double_ran3() > allele_F_AAF_prob) ? (0) : (1);

			if (__DUMP_PROXYTYPER_MSGS__)
			{
				fprintf(stderr, "%s: %s: orig_geno: %.0f; AFs: %.3f / %.3f;; Rands=(%.3f / %.3f); (Res: %d|%d)\n", var_id, tok_buff, orig_imp_geno, norm_M_AAF_prob, norm_F_AAF_prob, M_rand, F_rand, M_allele, F_allele);
			}

			// Sample allele1 and allele2.
			char cur_geno = M_allele*2 + F_allele;

			abs_geno_dist += fabs(get_genotype_per_haplocoded_genotype(cur_geno) - get_genotype_per_haplocoded_genotype(orig_imp_geno));

			//cur_var_geno_sig[i_s] = cur_geno;
			cur_var_geno_sig[i_s] = orig_imp_geno;
		} // i_s loop.

		fprintf(f_geno_rand_stats, "%s\t%.3f\t%.4f\t%d\n", var_id, abs_geno_dist, avg_geno_entropy, (int)sample_ids->size());

		// We must have read all of the genotypes at this point.
		if (t_string::get_next_token(cur_line, tok_buff, 1000, "\t", cur_char_i))
		{
			fprintf(stderr, "Still has genotypes after reading %d subjects: %s: %s\n", (int)sample_ids->size(), var_id, tok_buff);
			exit(0);
		}

		int l_var = t_string::string_length(ref_all);

		t_annot_region* cur_geno_reg = get_empty_region();
		cur_geno_reg->chrom = t_string::copy_me_str(chrom);
		cur_geno_reg->start = posn;
		cur_geno_reg->end = posn + l_var - 1;
		cur_geno_reg->strand = '+';
		cur_geno_reg->name = t_string::copy_me_str(var_id);

		void** reg_info = new void*[4];
		reg_info[0] = cur_var_geno_sig;
		reg_info[1] = NULL;

		cur_geno_reg->data = reg_info;

		geno_sig_regs->push_back(cur_geno_reg);

		if ((int)geno_sig_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "@ %d. variant.     \r", (int)geno_sig_regs->size());
		}

		delete[] cur_line;
	} // vcf file reading loop.
	
	fclose(f_geno_rand_stats);
	close_f(f_vcf, vcf_fp);

	// Save the genotype regions.
	binarize_variant_genotype_signal_regions(geno_sig_regs, NULL, sample_ids, op_fp);
} // convert_BEAGLE_imputed_meta_panel_2_matbed

inline double get_self(double val)
{
	return(val);
}

