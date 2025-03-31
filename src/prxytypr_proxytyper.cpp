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

	int ref_max_geno = get_max_genotype_value(known_genotype_regs, known_sample_ids);

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", vecsize(imputed_sample_ids));

	vector<t_annot_region*>* imputed_genotype_regs = load_variant_signal_regions_wrapper(imputed_genotypes_fp, imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", vecsize(imputed_genotype_regs));

	int imputed_max_geno = get_max_genotype_value(imputed_genotype_regs, imputed_sample_ids);

	if (ref_max_geno != imputed_max_geno)
	{
		fprintf(stderr, "Genotype encodings are not matching between ref/imputed panels: %d/%d (%s/%s)\n", 
			ref_max_geno,  imputed_max_geno, imputed_genotypes_fp, known_genotypes_fp);
		exit(1);
	}

	int max_geno = ref_max_geno;

	if (max_geno != 2 && max_geno != 3)
	{
		fprintf(stderr, "Genocoding is not recognized: %d\n", max_geno);
		exit(1);
	}

	for (int i_reg = 0; i_reg < vecsize(imputed_genotype_regs); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	for (int i_reg = 0; i_reg < vecsize(known_genotype_regs); i_reg++)
	{
		double AAC = 0;
		void** reg_info = (void**)(known_genotype_regs->at(i_reg)->data);
		char* geno_sig = (char*)(reg_info[0]);
		for (int i_s = 0; i_s < vecsize(known_sample_ids); i_s++)
		{
			if (max_geno == 2)
			{
				AAC += geno_sig[i_s];
			}
			else
			{
				int all1 = get_allele_per_haplotype(geno_sig[i_s], 0);
				int all2 = get_allele_per_haplotype(geno_sig[i_s], 1);
				AAC += (all1 + all2);
			}
		} // i_s loop.

		double n_total_haps = 2.0 * vecsize(known_sample_ids);

		known_genotype_regs->at(i_reg)->dbl_score = AAC / n_total_haps;
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
			//double AF = known_reg->dbl_score;
			t_string_tokens* toks = t_string::tokenize_by_chars(known_reg->name, "_");
			if (vecsize(toks) != 4)
			{
				fprintf(stderr, "Could not parse the reference AF from the variant identifier: %s\n", known_reg->name);
				exit(1);
			}

			double AF = atof(toks->at(3)->str());
			t_string::clean_tokens(toks);

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
				imp_reg->name, 
				cur_geno_corr * cur_geno_corr,
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

	fprintf(stderr, "Assigning AFs to %d variants (Max geno: %d)\n", vecsize(panel_var_regs), max_geno);

	// For each resampled subject, calculate the AF mismatch compared to the generating panel.
	// Calculate the allele frequencies of the generating panel.
	//for (int i_var = 0; i_var < vecsize(panel_var_regs); i_var++)
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
		exit(1);
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
		exit(1);
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
		exit(1);
	}
	else
	{
		fprintf(stderr, "Max geno is %d, not expected.\n", max_geno);
		exit(1);
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
				exit(1);
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
			exit(1);
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

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(char* reference_haplocoded_tag_target_geno_regs_fp,
	char* ref_sample_ids_list_fp,
	int genotype_sequence_index_in_info,
	int n_vars_per_block)
{
	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_haplocoded_tag_target_geno_regs_fp,
		ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);

	vector<t_var_block*>* blocks = generate_reduced_state_blocks_constant_size_blocks(reference_haplocoded_tag_target_geno_regs,
		ref_sample_ids,
		genotype_sequence_index_in_info,
		n_vars_per_block);

	return(blocks);
}

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs,
	vector<char*>* ref_sample_ids,
	int genotype_sequence_index_in_info,
	int n_vars_per_block)
{
	fprintf(stderr, "Generating reduced state haplotype blocks for %d variant regions with %d samples using blocks of %d variants.\n",
		vecsize(reference_haplocoded_tag_target_geno_regs),
		vecsize(ref_sample_ids),
		n_vars_per_block);

	// Setup the tag matrix for the current vicinity.
	int n_ref_haplotypes = 2 * vecsize(ref_sample_ids);

	vector<t_var_block*>* all_haplo_var_blocks = new vector<t_var_block*>();

	double** per_haplo_per_var_alleles = new double* [n_ref_haplotypes + 2];
	for (int hap_i = 0; hap_i < n_ref_haplotypes; hap_i++)
	{
		per_haplo_per_var_alleles[hap_i] = new double[2 * vecsize_t(reference_haplocoded_tag_target_geno_regs) + 5];

		// Reset all entries to -1.
		for (int i_var = 0;
			i_var < vecsize(reference_haplocoded_tag_target_geno_regs);
			i_var++)
		{
			per_haplo_per_var_alleles[hap_i][i_var] = -1;
		} // i_var loop.				
	} // hap_i loop.

	// Set the haplotypes.
	for (int j_var = 0;
		j_var < vecsize(reference_haplocoded_tag_target_geno_regs);
		j_var++)
	{
		void** cur_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(j_var)->data);
		char* cur_train_sample_haplocoded_geno = (char*)(cur_var_info[genotype_sequence_index_in_info]);

		for (int i_s = 0; i_s < vecsize(ref_sample_ids); i_s++)
		{
			int rel_var_i = j_var;
			per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 0);
			per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 1);
		} // i_s loop.
	} // j_var loop.

	vector<double*>* all_haplotypes = new vector<double*>();
	for (int i_hap = 0; i_hap < n_ref_haplotypes; i_hap++)
	{
		all_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
	} // i_hap loop.	

	fprintf(stderr, "Setup %d haplotypes on %d regions.\n", vecsize(all_haplotypes), vecsize(reference_haplocoded_tag_target_geno_regs));

	// Copy the current haplotypes.
	vector<double*>* cur_block_haplotypes = new vector<double*>();
	for (int i_hap = 0; i_hap < vecsize(all_haplotypes); i_hap++)
	{
		double* cur_block = new double[n_vars_per_block + 2];
		cur_block_haplotypes->push_back(cur_block);
	} // i_hap loop.

	int i_reg = 0;
	while (i_reg < vecsize(reference_haplocoded_tag_target_geno_regs))
	{
		// Copy the current haplotypes.
		if (__DUMP_PROXYTYPER_MSGS__)
		{
			fprintf(stderr, "Copying the current block haplotypes @ i_reg=%d\n", i_reg);
		}

		int n_var_per_cur_block = 0;
		for (int i_hap = 0; i_hap < vecsize(all_haplotypes); i_hap++)
		{
			n_var_per_cur_block = 0;
			for (int i_var = i_reg;
				i_var < MIN(vecsize(reference_haplocoded_tag_target_geno_regs), i_reg + n_vars_per_block);
				i_var++)
			{
				cur_block_haplotypes->at(i_hap)[i_var - i_reg] = all_haplotypes->at(i_hap)[i_var];
				n_var_per_cur_block++;
			} // i_var loop.

			cur_block_haplotypes->at(i_hap)[n_var_per_cur_block] = -1;
		} // i_hap loop.

		// Generate the unique haplotype blocks; assign the 
		if (__DUMP_PROXYTYPER_MSGS__)
		{
			fprintf(stderr, "Computing the unique haplotypes @ i_reg=%d\n", i_reg);
		}

		vector<vector<int>*>* per_uniq_haplo_indices = new vector<vector<int>*>();
		get_unique_haplotype_indices(cur_block_haplotypes, per_uniq_haplo_indices);

		t_var_block* cur_var_block = new t_var_block();
		cur_var_block->start_var_i = i_reg + 1; // These are 1-based.
		cur_var_block->end_var_i = i_reg + n_var_per_cur_block; // These are 1-based.
		cur_var_block->haplogroup_haplo_indices = per_uniq_haplo_indices;

		all_haplo_var_blocks->push_back(cur_var_block);

		if (__DUMP_PROXYTYPER_MSGS__)
		{
			fprintf(stderr, "Block %d-%d: %d unique haplotypes.\n", i_reg, i_reg + n_var_per_cur_block - 1, vecsize(per_uniq_haplo_indices));
		}

		// Move to the next block.
		i_reg = i_reg + n_var_per_cur_block;
	} // i_reg loop.

	for (int i_hap = 0; i_hap < vecsize(all_haplotypes); i_hap++)
	{
		delete[] cur_block_haplotypes->at(i_hap);
	} // i_hap loop.
	delete cur_block_haplotypes;

	for (int hap_i = 0; hap_i < n_ref_haplotypes; hap_i++)
	{
		delete[] per_haplo_per_var_alleles[hap_i];
	} // hap_i loop.
	delete[] per_haplo_per_var_alleles;

	delete all_haplotypes;

	return(all_haplo_var_blocks);
}

inline double get_self(double val)
{
	return(val);
}

void calculate_haplotype_emission_probabilities_per_reference_State_Reduced(char* reference_tag_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* recombination_rate_dir,
	double N_e_2_n_ref_haplotypes,
	double allele_err,
	int l_blocks,
	char* geno_probs_op_fp)
{
	//double N_e = pow(10, 4);
	//double allele_err = pow(10, -4);

	double(*XSUM)(double num1, double num2) = NULL;
	double(*XMUL)(double num1, double num2) = NULL;
	//double(*XDIV)(double num1, double num2) = NULL;
	//double(*XSUB)(double num1, double num2) = NULL;
	//bool(*XCOMPARE)(double num1, double num2) = NULL;
	double(*XCONVERT_LIN)(double num1) = NULL;
	//double(*XCONVERT_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LOG)(double num1) = NULL;
	//double(*XCONVERT_2_LIN)(double num1) = NULL;
	double ZERO = 0.0;

	fprintf(stderr, "Math mode log.\n");
	XSUM = xlog_sum;
	XMUL = xlog_mul;
	//XDIV = xlog_div;
	//XSUB = xlog_sub;
	//XCOMPARE = xlog_comp;
	XCONVERT_LIN = xlog;
	//XCONVERT_LOG = get_self;
	XCONVERT_2_LOG = get_self;
	//XCONVERT_2_LIN = exp;
	ZERO = xlog(0);

	// This is the tag+target region genotypes of the reference panel.
	vector<t_annot_region*>* reference_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		vecsize(reference_haplocoded_tag_geno_regs),
		vecsize(ref_sample_ids),
		n_ref_haplotypes);

	double N_e = N_e_2_n_ref_haplotypes * n_ref_haplotypes;

	fprintf(stderr, "Calculating LL using %d reference haplotypes using N_e=%.3f, allele_eps=%.3f.\n", (int)(n_ref_haplotypes), N_e, allele_err);

	// These are the typed genotypes of the testing panel.
	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", vecsize(testing_haplocoded_tag_geno_regs), vecsize(testing_sample_ids));

	// Assign the typed variants to the variants of the reference panel.
	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", vecsize(intersects));
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	int n_states = n_ref_haplotypes;
	//int n_symbols = 2;

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Running forward-backward on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Generate the blocks and the reduced states: These must be generated before ends are added.
		vector<t_annot_region*>* cur_win_var_regs = cur_chr_testing_haplocoded_tag_geno_regs;
		vector<t_var_block*>* haplo_var_blocks = generate_reduced_state_blocks_constant_size_blocks(cur_win_var_regs,
			ref_sample_ids,
			1,
			l_blocks);

		// Set the start and end state regions.
		t_annot_region* start_state_reg = NULL;
		cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		int n_vars = vecsize(cur_win_var_regs) - 1;

		// Add the end state region.
		t_annot_region* end_state_reg = NULL;
		cur_win_var_regs->push_back(end_state_reg);

		//vector<t_annot_region*>* prev_win_var_regs = NULL;
		for (int test_sample_i = 0; test_sample_i < vecsize(testing_sample_ids); test_sample_i++)
		{
			fprintf(stderr, "Calculating probabilities for %d. subject.\n", test_sample_i);

			//fprintf(stderr, "Window (%d/%d/%d tags)::%d targets: %s:%d-%d (%d-%d) spanning %.4f cM between %.4f-%.4f cMs [%d haplo-blocks of %d variants]\n",
			//	vecsize(cur_win_var_regs), max_n_tag_vars_per_window, win_end_var_i - win_start_var_i + 1,
			//	n_target_vars_in_cur_win,
			//	cur_win_var_regs->at(0)->chrom, cur_win_var_regs->at(0)->start, cur_win_var_regs->back()->start,
			//	win_start_var_i, win_end_var_i,
			//	fabs(start_cM - end_cM),
			//	start_cM, end_cM,
			//	vecsize(haplo_var_blocks), l_blocks);

			if (__DUMP_PROXYTYPER_MSGS__)
			{
				fprintf(stderr, "%d haplotype blocks:\n", vecsize(haplo_var_blocks));
				for (int i_block = 0; i_block < vecsize(haplo_var_blocks); i_block++)
				{
					fprintf(stderr, "%d-%d: %d haplotypes\n",
						haplo_var_blocks->at(i_block)->start_var_i, haplo_var_blocks->at(i_block)->end_var_i,
						vecsize(haplo_var_blocks->at(i_block)->haplogroup_haplo_indices));
				} // i_block loop.
			}

			// Allocate and initialize the forward/backward arrays.
			double*** fore_scores_per_hap = new double** [n_vars + 2];
			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				fore_scores_per_hap[var_i] = new double* [2];

				for (int test_hap = 0; test_hap < 2; test_hap++)
				{
					fore_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
					memset(fore_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 1));
				}
			} // i loop.

			// Initialize the state probabilities for both haplotypes.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				fore_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				fore_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
			} // state_i loop.

			for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
			{
				// Process each block sequentially.
				for (int block_i = 0; block_i < vecsize(haplo_var_blocks); block_i++)
				{
					// Process the current block.
					int cur_block_start_var_i = haplo_var_blocks->at(block_i)->start_var_i;
					int cur_block_end_var_i = haplo_var_blocks->at(block_i)->end_var_i;

					int n_cur_block_haplogrps = vecsize(haplo_var_blocks->at(block_i)->haplogroup_haplo_indices);
					vector<vector<int>*>* per_hplgrp_hap_i = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices;

					if (__DUMP_PROXYTYPER_MSGS__)
					{
						fprintf(stderr, "Processing test_hap_i:%d;; Block:%d-%d with %d haplogroups.\n",
							test_hap_i,
							cur_block_start_var_i, cur_block_end_var_i,
							n_cur_block_haplogrps);
					}

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// The first haplotype in each haplogroup is the representative. No specific reason.
					vector<int>* per_hplgrp_representative_hap_i = new vector<int>();
					for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
					{
						int cur_grp_representative_hap_i = per_hplgrp_hap_i->at(i_grp)->at(0);

						per_hplgrp_representative_hap_i->push_back(cur_grp_representative_hap_i);

						if (__DUMP_PROXYTYPER_MSGS__)
						{
							fprintf(stderr, "Block %d::Haplogroup %d; representative_hap_i: %d\n",
								block_i,
								i_grp,
								cur_grp_representative_hap_i);
						}
					} // i_grp loop.

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// 2. Compute the total other transition for all groups in this block.
					for (int var_i = cur_block_start_var_i;
						var_i <= cur_block_end_var_i;
						var_i++)
					{
						if (__DUMP_PROXYTYPER_MSGS__)
						{
							if (var_i % 10 == 0)
							{
								fprintf(stderr, "Forward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, vecsize(testing_sample_ids), var_i);
							}
						}

						//void** cur_tag_var_info = NULL;
						char* test_sample_geno = NULL;
						char* ref_sample_geno = NULL;

						if (var_i <= n_vars)
						{
							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							test_sample_geno = (char*)(cur_tag_var_info[0]);
							ref_sample_geno = (char*)(cur_tag_var_info[1]);
						}

						/////////////////////////////////////////////////////
						// Pre-compute the transition probabilities:
						double prev_var_cM = cur_win_var_regs->at(1)->dbl_score;
						if (var_i > 1)
						{
							prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
						}

						double cur_var_cM = cur_win_var_regs->at(n_vars)->dbl_score;
						if (var_i <= n_vars)
						{
							cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
						}

						double r_m = fabs(cur_var_cM - prev_var_cM);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

						double other_prob = tau_m / n_ref_haplotypes;
						double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

						//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
						/////////////////////////////////////////////////////
						// Compute the total other-trans probabilities for each haplogroup.
						vector<double>* per_hpl_grp_other_trans_total = new vector<double>();
						for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
						{
							// The transition and emission probabilities are constant for each haplogroup.
							double cur_grp_total_other_trans_prob = ZERO;

							// Transition prob: Other trans.
							double trans_prob = other_prob;

							for (int i_hap_i = 0;
								i_hap_i < vecsize(per_hplgrp_hap_i->at(i_grp));
								i_hap_i++)
							{
								int cur_grp_state_i = per_hplgrp_hap_i->at(i_grp)->at(i_hap_i);

								cur_grp_total_other_trans_prob = XSUM(cur_grp_total_other_trans_prob,
									XMUL(XCONVERT_LIN(trans_prob), fore_scores_per_hap[var_i - 1][test_hap_i][cur_grp_state_i]));
							} // i_sample_i loop.

							per_hpl_grp_other_trans_total->push_back(cur_grp_total_other_trans_prob);

							if (__DUMP_PROXYTYPER_MSGS__)
							{
								fprintf(stderr, "Block %d::Haplogroup %d::Total Other Transition: %.4f\n",
									block_i,
									i_grp,
									cur_grp_total_other_trans_prob);
							}
						} // i_grp loop.

						// Loop over all the states.
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							fore_scores_per_hap[var_i][test_hap_i][state_i] = ZERO;

							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// Do a self transition with a delta probability.
							double self_min_trans_prob = self_prob - other_prob;

							int ref_sample_i = state_i / 2;
							int ref_hap_i = state_i % 2;
							int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
							int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

							double emit_prob_per_ref_allele[2];
							emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
							emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
							double emit_prob = emit_prob_per_ref_allele[cur_test_allele];

							double trans_emit_prob = XMUL(XCONVERT_LIN(self_min_trans_prob), XCONVERT_LIN(emit_prob));
							//trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

							fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i],
								XMUL(trans_emit_prob, fore_scores_per_hap[var_i - 1][test_hap_i][state_i]));

							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// Now do a transition to all the top scoring haplotypes.
							for (int hpl_grp_i = 0; hpl_grp_i < n_cur_block_haplogrps; hpl_grp_i++)
							{
								double cur_haplogrp_total_other_trans = per_hpl_grp_other_trans_total->at(hpl_grp_i);

								// Apply the scaler.
								double emit_prob = emit_prob_per_ref_allele[cur_test_allele];
								double trans_emit_prob = XMUL(XCONVERT_LIN(1.0), XCONVERT_LIN(emit_prob));
								//trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

								fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i],
									XMUL(trans_emit_prob, cur_haplogrp_total_other_trans));
							} // hpl_grp_i loop.
						} // state_i loop.

						delete per_hpl_grp_other_trans_total;
					} // var_i loop.

					delete per_hplgrp_representative_hap_i;
				} // block_i loop.
			} // test_hap_i loop.

			// Compute the total log forward and backward probabilities.
			double per_hap_total_log_fore_prob[2];
			per_hap_total_log_fore_prob[0] = XCONVERT_LIN(0.0);
			per_hap_total_log_fore_prob[1] = XCONVERT_LIN(0.0);

			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				per_hap_total_log_fore_prob[0] = XSUM(per_hap_total_log_fore_prob[0], fore_scores_per_hap[n_vars][0][state_i]);
				per_hap_total_log_fore_prob[1] = XSUM(per_hap_total_log_fore_prob[1], fore_scores_per_hap[n_vars][1][state_i]);
			} // state_i loop.

			fprintf(stderr, "Sample %d: Haplotype 0 total probabilities: fore=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[0]));
			fprintf(stderr, "Sample %d: Haplotype 1 total probabilities: fore=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[1]));

			//// Free haplogroup and block memories.
			//for (int i_block = 0; i_block < vecsize(haplo_var_blocks); i_block++)
			//{
			//	for (int i_grp = 0; i_grp < vecsize(haplo_var_blocks->at(i_block)->haplogroup_haplo_indices); i_grp++)
			//	{
			//		delete haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->at(i_grp);
			//	}
			//} // i_block loop.

			// Allocate and initialize the forward/backward arrays.
			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				for (int test_hap = 0; test_hap < 2; test_hap++)
				{
					delete[] fore_scores_per_hap[var_i][test_hap];
				} // test_hap loop.

				delete[] fore_scores_per_hap[var_i];
			} // i loop.

			delete[] fore_scores_per_hap;
		} // test_sample_i loop.
	} // i_chr loop.

	fprintf(stderr, "Done!\n");
}


static void* emission_probability_reduced_state_computing_thread_callback(void* thread_info_ptr)
{
	double(*XSUM)(double num1, double num2) = NULL;
	double(*XMUL)(double num1, double num2) = NULL;
	//double(*XDIV)(double num1, double num2) = NULL;
	//double(*XSUB)(double num1, double num2) = NULL;
	//bool(*XCOMPARE)(double num1, double num2) = NULL;
	double(*XCONVERT_LIN)(double num1) = NULL;
	//double(*XCONVERT_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LOG)(double num1) = NULL;
	//double(*XCONVERT_2_LIN)(double num1) = NULL;
	double ZERO = 0.0;

	fprintf(stderr, "Math mode log.\n");
	XSUM = xlog_sum;
	XMUL = xlog_mul;
	//XDIV = xlog_div;
	//XSUB = xlog_sub;
	//XCOMPARE = xlog_comp;
	XCONVERT_LIN = xlog;
	//XCONVERT_LOG = get_self;
	XCONVERT_2_LOG = get_self;
	//XCONVERT_2_LIN = exp;
	ZERO = xlog(0);

	void** thread_ptrs_list = (void**)(thread_info_ptr);

	// Copy the thread information.
	int* which_out_of_ptr = (int*)(thread_ptrs_list[0]);
	int thread_i = which_out_of_ptr[0];
	int n_threads = which_out_of_ptr[1];

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[1]);
	vector<char*>* testing_sample_ids = (vector<char*>*)(thread_ptrs_list[2]);

	//t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[3]);
	vector<char*>* ref_sample_ids = (vector<char*>*)(thread_ptrs_list[4]);

	//char* recombination_rate_dir = (char*)(thread_ptrs_list[5]);

	double* N_e_allele_err_ptr = (double*)(thread_ptrs_list[6]);
	double N_e = N_e_allele_err_ptr[0];
	double allele_err = N_e_allele_err_ptr[1];
	vector<double*>* per_sample_per_hap_emission_probs = (vector<double*>*)(thread_ptrs_list[7]);

	vector<t_var_block*>** per_chrom_haplo_blocks = (vector<t_var_block*>**)(thread_ptrs_list[8]);

	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	int n_states = n_ref_haplotypes;

	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Running forward-backward on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		//char cur_chr_recombination_rate_fp[1000];
		//sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		//vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		//if (cur_chrom_recomb_regs == NULL)
		//{
		//	fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		//	exit(1);
		//}

		//sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		//vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		//// Assign the recomb rates.
		//fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		//for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		//for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		// Generate the blocks and the reduced states: These must be generated before ends are added.
		vector<t_annot_region*>* cur_win_var_regs = cur_chr_testing_haplocoded_tag_geno_regs;
		vector<t_var_block*>* haplo_var_blocks = per_chrom_haplo_blocks[i_chr];

		//// Set the start and end state regions.
		//t_annot_region* start_state_reg = NULL;
		//cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		int n_vars = vecsize(cur_win_var_regs) - 2;

		//// Add the end state region.
		//t_annot_region* end_state_reg = NULL;
		//cur_win_var_regs->push_back(end_state_reg);

		fprintf(stderr, "Thread-%d allocating scores array..\n", thread_i);
		double*** fore_scores_per_hap = new double** [2];
		for (int test_hap = 0; test_hap < 2; test_hap++)
		{
			double* cur_mem_pool = new double[n_states * (n_vars + 2)];
			memset(cur_mem_pool, 0, sizeof(double) * n_states * (n_vars + 2));

			fore_scores_per_hap[test_hap] = new double* [n_vars + 2];

			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				fore_scores_per_hap[test_hap][var_i] = cur_mem_pool + var_i * n_states;
			} // var_i loop.
		} // test_hap loop.
		fprintf(stderr, "Thread-%d completed allocating scores array..\n", thread_i);

		for (int test_sample_i = 0; test_sample_i < vecsize(testing_sample_ids); test_sample_i++)
		{
			if (test_sample_i % n_threads != thread_i)
			{
				continue;
			}

			if (__DUMP_PROXYTYPER_MSGS__)
			{
				fprintf(stderr, "Calculating probabilities for %d. subject.\n", test_sample_i);
			}

			//fprintf(stderr, "Window (%d/%d/%d tags)::%d targets: %s:%d-%d (%d-%d) spanning %.4f cM between %.4f-%.4f cMs [%d haplo-blocks of %d variants]\n",
			//	vecsize(cur_win_var_regs), max_n_tag_vars_per_window, win_end_var_i - win_start_var_i + 1,
			//	n_target_vars_in_cur_win,
			//	cur_win_var_regs->at(0)->chrom, cur_win_var_regs->at(0)->start, cur_win_var_regs->back()->start,
			//	win_start_var_i, win_end_var_i,
			//	fabs(start_cM - end_cM),
			//	start_cM, end_cM,
			//	vecsize(haplo_var_blocks), l_blocks);

			if (__DUMP_PROXYTYPER_MSGS__)
			{
				fprintf(stderr, "%d haplotype blocks:\n", vecsize(haplo_var_blocks));
				for (int i_block = 0; i_block < vecsize(haplo_var_blocks); i_block++)
				{
					fprintf(stderr, "%d-%d: %d haplotypes\n",
						haplo_var_blocks->at(i_block)->start_var_i, haplo_var_blocks->at(i_block)->end_var_i,
						vecsize(haplo_var_blocks->at(i_block)->haplogroup_haplo_indices));
				} // i_block loop.
			}

			// Allocate and initialize the forward/backward arrays: This section is blocks over all threads while allocating memory.

			//double*** fore_scores_per_hap = new double** [n_vars + 2];
			//for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			//{
			//	double* per_hap_mem_pool = new double [2 * (n_states+2)];
			//	//fore_scores_per_hap[var_i] = &per_hap_mem_pool;
			//	memset(per_hap_mem_pool, 0, sizeof(double) * (2 * (n_states + 2)));

			//	for (int test_hap = 0; test_hap < 2; test_hap++)
			//	{
			//		//fore_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
			//		//memset(fore_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 1));
			//		fore_scores_per_hap[var_i][test_hap] = per_hap_mem_pool + test_hap * n_states;
			//	}
			//} // i loop.


			// Initialize the state probabilities for both haplotypes.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				/*fore_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				fore_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);*/
				fore_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				fore_scores_per_hap[1][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes); 
			} // state_i loop.

			for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
			{
				if ((test_sample_i - thread_i) % (5 * n_threads) == 0)
				{
					fprintf(stderr, "Thread %d: Computing fore-ward scores for sample %d/%d; hap: %d\n", thread_i, test_sample_i, vecsize(testing_sample_ids), test_hap_i);
				}

				// Process each block sequentially.
				for (int block_i = 0; block_i < vecsize(haplo_var_blocks); block_i++)
				{
					// Process the current block.
					int cur_block_start_var_i = haplo_var_blocks->at(block_i)->start_var_i;
					int cur_block_end_var_i = haplo_var_blocks->at(block_i)->end_var_i;

					int n_cur_block_haplogrps = vecsize(haplo_var_blocks->at(block_i)->haplogroup_haplo_indices);
					vector<vector<int>*>* per_hplgrp_hap_i = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices;

					if (__DUMP_PROXYTYPER_MSGS__)
					{
						fprintf(stderr, "Processing test_hap_i:%d;; Block:%d-%d with %d haplogroups.\n",
							test_hap_i,
							cur_block_start_var_i, cur_block_end_var_i,
							n_cur_block_haplogrps);
					}

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// The first haplotype in each haplogroup is the representative. No specific reason.
					vector<int>* per_hplgrp_representative_hap_i = new vector<int>();
					for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
					{
						int cur_grp_representative_hap_i = per_hplgrp_hap_i->at(i_grp)->at(0);

						per_hplgrp_representative_hap_i->push_back(cur_grp_representative_hap_i);

						if (__DUMP_PROXYTYPER_MSGS__)
						{
							fprintf(stderr, "Block %d::Haplogroup %d; representative_hap_i: %d\n",
								block_i,
								i_grp,
								cur_grp_representative_hap_i);
						}
					} // i_grp loop.

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// 2. Compute the total other transition for all groups in this block.
					for (int var_i = cur_block_start_var_i;
						var_i <= cur_block_end_var_i;
						var_i++)
					{
						if (__DUMP_PROXYTYPER_MSGS__)
						{
							if (var_i % 10 == 0)
							{
								fprintf(stderr, "Forward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, vecsize(testing_sample_ids), var_i);
							}
						}

						//void** cur_tag_var_info = NULL;
						char* test_sample_geno = NULL;
						char* ref_sample_geno = NULL;

						if (var_i <= n_vars)
						{
							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							test_sample_geno = (char*)(cur_tag_var_info[0]);
							ref_sample_geno = (char*)(cur_tag_var_info[1]);
						}

						/////////////////////////////////////////////////////
						// Pre-compute the transition probabilities:
						double prev_var_cM = cur_win_var_regs->at(1)->dbl_score;
						if (var_i > 1)
						{
							prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
						}

						double cur_var_cM = cur_win_var_regs->at(n_vars)->dbl_score;
						if (var_i <= n_vars)
						{
							cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
						}

						double r_m = fabs(cur_var_cM - prev_var_cM);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

						double other_prob = tau_m / n_ref_haplotypes;
						double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

						//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
						/////////////////////////////////////////////////////
						// Compute the total other-trans probabilities for each haplogroup.
						vector<double>* per_hpl_grp_other_trans_total = new vector<double>();
						for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
						{
							// The transition and emission probabilities are constant for each haplogroup.
							double cur_grp_total_other_trans_prob = ZERO;

							// Transition prob: Other trans.
							double trans_prob = other_prob;

							for (int i_hap_i = 0;
								i_hap_i < vecsize(per_hplgrp_hap_i->at(i_grp));
								i_hap_i++)
							{
								int cur_grp_state_i = per_hplgrp_hap_i->at(i_grp)->at(i_hap_i);

								cur_grp_total_other_trans_prob = XSUM(cur_grp_total_other_trans_prob,
									XMUL(XCONVERT_LIN(trans_prob), fore_scores_per_hap[test_hap_i][var_i - 1][cur_grp_state_i]));
									//XMUL(XCONVERT_LIN(trans_prob), fore_scores_per_hap[var_i - 1][test_hap_i][cur_grp_state_i]));
							} // i_sample_i loop.

							per_hpl_grp_other_trans_total->push_back(cur_grp_total_other_trans_prob);

							if (__DUMP_PROXYTYPER_MSGS__)
							{
								fprintf(stderr, "Block %d::Haplogroup %d::Total Other Transition: %.4f\n",
									block_i,
									i_grp,
									cur_grp_total_other_trans_prob);
							}
						} // i_grp loop.

						// Loop over all the states.
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							//fore_scores_per_hap[var_i][test_hap_i][state_i] = ZERO;
							fore_scores_per_hap[test_hap_i][var_i][state_i] = ZERO;

							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// Do a self transition with a delta probability.
							double self_min_trans_prob = self_prob - other_prob;

							int ref_sample_i = state_i / 2;
							int ref_hap_i = state_i % 2;
							int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
							int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

							double emit_prob_per_ref_allele[2];
							emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
							emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
							double emit_prob = emit_prob_per_ref_allele[cur_test_allele];

							double trans_emit_prob = XMUL(XCONVERT_LIN(self_min_trans_prob), XCONVERT_LIN(emit_prob));
							//trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

							fore_scores_per_hap[test_hap_i][var_i][state_i] = XSUM(fore_scores_per_hap[test_hap_i][var_i][state_i],
								XMUL(trans_emit_prob, fore_scores_per_hap[test_hap_i][var_i - 1][state_i]));
							//fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i],
							//	XMUL(trans_emit_prob, fore_scores_per_hap[var_i - 1][test_hap_i][state_i]));

							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// Now do a transition to all the top scoring haplotypes.
							for (int hpl_grp_i = 0; hpl_grp_i < n_cur_block_haplogrps; hpl_grp_i++)
							{
								double cur_haplogrp_total_other_trans = per_hpl_grp_other_trans_total->at(hpl_grp_i);

								// Apply the scaler.
								double emit_prob = emit_prob_per_ref_allele[cur_test_allele];
								double trans_emit_prob = XMUL(XCONVERT_LIN(1.0), XCONVERT_LIN(emit_prob));
								//trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

								fore_scores_per_hap[test_hap_i][var_i][state_i] = XSUM(fore_scores_per_hap[test_hap_i][var_i][state_i],
									XMUL(trans_emit_prob, cur_haplogrp_total_other_trans));
								/*fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i],
									XMUL(trans_emit_prob, cur_haplogrp_total_other_trans));*/
							} // hpl_grp_i loop.
						} // state_i loop.

						delete per_hpl_grp_other_trans_total;
					} // var_i loop.

					delete per_hplgrp_representative_hap_i;
				} // block_i loop.
			} // test_hap_i loop.

			// Compute the total log forward and backward probabilities.
			double per_hap_total_log_fore_prob[2];
			per_hap_total_log_fore_prob[0] = XCONVERT_LIN(0.0);
			per_hap_total_log_fore_prob[1] = XCONVERT_LIN(0.0);

			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				//per_hap_total_log_fore_prob[0] = XSUM(per_hap_total_log_fore_prob[0], fore_scores_per_hap[n_vars][0][state_i]);
				//per_hap_total_log_fore_prob[1] = XSUM(per_hap_total_log_fore_prob[1], fore_scores_per_hap[n_vars][1][state_i]);
				per_hap_total_log_fore_prob[0] = XSUM(per_hap_total_log_fore_prob[0], fore_scores_per_hap[0][n_vars][state_i]);
				per_hap_total_log_fore_prob[1] = XSUM(per_hap_total_log_fore_prob[1], fore_scores_per_hap[1][n_vars][state_i]);
			} // state_i loop.

			if (__DUMP_PROXYTYPER_MSGS__)
			{
				fprintf(stderr, "Sample %d: Haplotype 0 total probabilities: fore=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[0]));
				fprintf(stderr, "Sample %d: Haplotype 1 total probabilities: fore=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[1]));
			}

			per_sample_per_hap_emission_probs->at(test_sample_i)[0] = per_hap_total_log_fore_prob[0];
			per_sample_per_hap_emission_probs->at(test_sample_i)[1] = per_hap_total_log_fore_prob[1];

			//// Free haplogroup and block memories.
			//for (int i_block = 0; i_block < vecsize(haplo_var_blocks); i_block++)
			//{
			//	for (int i_grp = 0; i_grp < vecsize(haplo_var_blocks->at(i_block)->haplogroup_haplo_indices); i_grp++)
			//	{
			//		delete haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->at(i_grp);
			//	}
			//} // i_block loop.

			//// Allocate and initialize the forward/backward arrays.
			//for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			//{
			//	for (int test_hap = 0; test_hap < 2; test_hap++)
			//	{
			//		delete[] fore_scores_per_hap[var_i][test_hap];
			//	} // test_hap loop.

			//	delete[] fore_scores_per_hap[var_i];
			//} // i loop.

			//delete[] fore_scores_per_hap;
		} // test_sample_i loop.
	} // i_chr loop.

	return(NULL);
}

void calculate_haplotype_emission_probabilities_per_reference_State_Reduced_multithreaded(char* reference_tag_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* recombination_rate_dir,
	double N_e_2_n_ref_haplotypes,
	double allele_err,
	int l_blocks,
	int n_threads,
	char* geno_probs_op_fp)
{
	//double N_e = pow(10, 4);
	//double allele_err = pow(10, -4);

	//double(*XSUM)(double num1, double num2) = xlog_sum;
	//double(*XMUL)(double num1, double num2) = xlog_mul;
	//double(*XDIV)(double num1, double num2) = xlog_div;
	//double(*XSUB)(double num1, double num2) = xlog_sub;
	//bool(*XCOMPARE)(double num1, double num2) = xlog_comp;
	//double(*XCONVERT_LIN)(double num1) = xlog;
	//double(*XCONVERT_LOG)(double num1) = get_self;
	//double(*XCONVERT_2_LOG)(double num1) = get_self;
	//double(*XCONVERT_2_LIN)(double num1) = exp;
	//double ZERO = xlog(0);

	//fprintf(stderr, "Math mode log.\n");
	//XSUM = xlog_sum;
	//XMUL = xlog_mul;
	//XDIV = xlog_div;
	//XSUB = xlog_sub;
	//XCOMPARE = xlog_comp;
	//XCONVERT_LIN = xlog;
	//XCONVERT_LOG = get_self;
	//XCONVERT_2_LOG = get_self;
	//XCONVERT_2_LIN = exp;
	//ZERO = xlog(0);

	// This is the tag+target region genotypes of the reference panel.
	vector<t_annot_region*>* reference_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		vecsize(reference_haplocoded_tag_geno_regs),
		vecsize(ref_sample_ids),
		n_ref_haplotypes);

	double N_e = N_e_2_n_ref_haplotypes * n_ref_haplotypes;

	fprintf(stderr, "Calculating LL using %d reference haplotypes using N_e=%.3f, allele_eps=%.3f.\n", (int)(n_ref_haplotypes), N_e, allele_err);

	// These are the typed genotypes of the testing panel.
	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", vecsize(testing_haplocoded_tag_geno_regs), vecsize(testing_sample_ids));

	// Assign the typed variants to the variants of the reference panel.
	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", vecsize(intersects));
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	//int n_states = n_ref_haplotypes;
	//int n_symbols = 2;

	// Set the emission probabilities.
	vector<double*>* per_sample_per_hap_emission_probs = new vector<double*>();
	for (int i_s = 0; i_s < vecsize(testing_sample_ids); i_s++)
	{
		double* cur_subj_per_allele_probs = new double[2];
		cur_subj_per_allele_probs[0] = xlog(0);
		cur_subj_per_allele_probs[1] = xlog(0);
		per_sample_per_hap_emission_probs->push_back(cur_subj_per_allele_probs);
	} // i_s loop.

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	vector<t_var_block*>** per_chrom_haplo_blocks = new vector<t_var_block*>*[vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids)];

	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Setting up variant info on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Generate the blocks and the reduced states: These must be generated before ends are added.
		vector<t_annot_region*>* cur_win_var_regs = cur_chr_testing_haplocoded_tag_geno_regs;
		vector<t_var_block*>* haplo_var_blocks = generate_reduced_state_blocks_constant_size_blocks(cur_win_var_regs,
			ref_sample_ids,
			1,
			l_blocks);

		per_chrom_haplo_blocks[i_chr] = haplo_var_blocks;

		// Set the start and end state regions.
		t_annot_region* start_state_reg = NULL;
		cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		//int n_vars = vecsize(cur_win_var_regs) - 1;

		// Add the end state region.
		t_annot_region* end_state_reg = NULL;
		cur_win_var_regs->push_back(end_state_reg);
	} // i_chr loop.

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting %d. forward variable thread..\n", thread_i);

		void** thread_ptrs_list = new void* [12];

		int* which_out_of_ptr = new int[10];
		which_out_of_ptr[0] = thread_i;
		which_out_of_ptr[1] = n_threads;
		thread_ptrs_list[0] = which_out_of_ptr;

		//t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = NULL;
		thread_ptrs_list[1] = restr_testing_haplocoded_tag_geno_regs;

		//vector<char*>* testing_sample_ids = NULL;
		thread_ptrs_list[2] = testing_sample_ids;

		//t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = NULL;
		thread_ptrs_list[3] = restr_reference_haplocoded_tag_target_geno_regs;

		//vector<char*>* ref_sample_ids = NULL;
		thread_ptrs_list[4] = ref_sample_ids;

		//char* recombination_rate_dir = NULL;
		thread_ptrs_list[5] = recombination_rate_dir;

		//double* N_e_allele_err_ptr = NULL;
		double* N_e_allele_err_ptr = new double[5];
		N_e_allele_err_ptr[0] = N_e;
		N_e_allele_err_ptr[1] = allele_err;
		thread_ptrs_list[6] = N_e_allele_err_ptr;

		thread_ptrs_list[7] = per_sample_per_hap_emission_probs;

		thread_ptrs_list[8] = per_chrom_haplo_blocks;

		t_ansi_thread* cur_thread = new t_ansi_thread(emission_probability_reduced_state_computing_thread_callback, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Started %d/%d threads; waiting.\n", vecsize(resampling_threads), n_threads);

	for (int thread_i = 0; thread_i < vecsize(resampling_threads); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		fprintf(stderr, "%d. thread finished.\n", thread_i);
	} // thread_i waiting loop. 

	fprintf(stderr, "Saving results to %s\n", geno_probs_op_fp);
	FILE* f_per_test_subj_emission_prob = open_f(geno_probs_op_fp, "w");
	for (int i_s = 0; i_s < vecsize(testing_sample_ids); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			fprintf(f_per_test_subj_emission_prob, "%s\t%d\t%d\t%.6f\n", testing_sample_ids->at(i_s), i_s, i_hap, per_sample_per_hap_emission_probs->at(i_s)[i_hap]);
		} // i_hap loop.
	} // i_s loop.
	close_f(f_per_test_subj_emission_prob, geno_probs_op_fp);

	fprintf(stderr, "Done!\n");
}

static void* emission_probability_computing_thread_callback(void* thread_info_ptr)
{
	void** thread_ptrs_list = (void**)(thread_info_ptr);

	// Copy the thread information.
	int* which_out_of_ptr = (int*)(thread_ptrs_list[0]);
	int thread_i = which_out_of_ptr[0];
	int n_threads = which_out_of_ptr[1];

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[1]);
	vector<char*>* testing_sample_ids = (vector<char*>*)(thread_ptrs_list[2]);

	//t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[3]);
	vector<char*>* ref_sample_ids = (vector<char*>*)(thread_ptrs_list[4]);

	//char* recombination_rate_dir = (char*)(thread_ptrs_list[5]);

	double* N_e_allele_err_ptr = (double*)(thread_ptrs_list[6]);
	double N_e = N_e_allele_err_ptr[0];
	double allele_err = N_e_allele_err_ptr[1];
	vector<double*>* per_sample_per_hap_emission_probs = (vector<double*>*)(thread_ptrs_list[7]);

	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	int n_states = n_ref_haplotypes;

	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Running forward on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		//char cur_chr_recombination_rate_fp[1000];
		//sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		//vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		//if (cur_chrom_recomb_regs == NULL)
		//{
		//	fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		//	exit(1);
		//}

		//// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		//vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		//// Assign the recomb rates.
		//fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		//for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		//for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.		

		//// Set the start and end state regions.
		//t_annot_region* start_state_reg = get_empty_region();
		//start_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		//start_state_reg->start = 1;
		//start_state_reg->end = 1;
		//start_state_reg->strand = '+';
		//start_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->at(0)->dbl_score;
		//cur_chr_testing_haplocoded_tag_geno_regs->insert(cur_chr_testing_haplocoded_tag_geno_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		int n_vars = vecsize(cur_chr_testing_haplocoded_tag_geno_regs) - 2;

		//// Add the end state region.
		//t_annot_region* end_state_reg = get_empty_region();
		//end_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		//end_state_reg->start = 1;
		//end_state_reg->end = 1;
		//end_state_reg->strand = '+';
		//end_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->back()->dbl_score;
		//cur_chr_testing_haplocoded_tag_geno_regs->push_back(end_state_reg);

		for (int test_sample_i = 0; test_sample_i < vecsize(testing_sample_ids); test_sample_i++)
		{
			//fprintf(stderr, "Computing forward probabilities for sample %d\n", test_sample_i);

			if (test_sample_i % n_threads != thread_i)
			{
				continue;
			}

			for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
			{
				//if ((test_sample_i - thread_i) % 100 == 0)
				{
					fprintf(stderr, "Thread %d: Computing fore-ward scores for sample %d; hap: %d\n", thread_i, test_sample_i, test_hap_i);
				}

				// Allocate and initialize the forward array.
				double** fore_scores = new double* [n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					fore_scores[var_i] = new double[n_states + 2];
					memset(fore_scores[var_i], 0, sizeof(double) * (n_states + 1));
				} // i loop.

				// Initialize the state probabilities.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					fore_scores[0][state_i] = xlog((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				// Start recursing over the variants.
				for (int var_i = 1; var_i <= n_vars + 1; var_i++)
				{
					//if (var_i % 100 == 0)
					//{
					//	fprintf(stderr, "Forward: sample_i: %d/%d [%d]: var_i: %d         \r", test_sample_i, vecsize(testing_sample_ids), test_hap_i, var_i);
					//}

					void** cur_tag_var_info = (void**)(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->data);

					char* test_sample_geno = NULL;
					char* ref_sample_geno = NULL;
					if (var_i <= n_vars)
					{
						test_sample_geno = (char*)(cur_tag_var_info[0]);
						ref_sample_geno = (char*)(cur_tag_var_info[1]);
					}

					/////////////////////////////////////////////////////
					// Pre-compute the transition probabilities:
					double prev_var_cM = 0;
					if (var_i > 1)
					{
						prev_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i - 1)->dbl_score;
					}

					double cur_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score;

					double r_m = fabs(cur_var_cM - prev_var_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

					//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
					/////////////////////////////////////////////////////

					// Loop over all the states.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						// Recurse over the previous states.
						fore_scores[var_i][state_i] = xlog(0);
						for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
						{
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the transition probabilities.
							double trans_prob = 0;
							if (state_i == prev_state_i)
							{
								trans_prob = self_prob;
							}
							else
							{
								trans_prob = other_prob;
							}

							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
							double emit_prob = 0;

							if (var_i == n_vars + 1)
							{
								emit_prob = 1.0;
							}
							else
							{
								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);
								if (cur_ref_allele == cur_test_allele)
								{
									emit_prob = 1 - allele_err;
								}
								else
								{
									emit_prob = allele_err;
								}
							}
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

							double trans_emit_prob = xlog_mul(xlog(trans_prob), xlog(emit_prob));

							fore_scores[var_i][state_i] = xlog_sum(fore_scores[var_i][state_i], xlog_mul(trans_emit_prob, fore_scores[var_i - 1][prev_state_i]));

							if (__DUMP_PROXYTYPER_MSGS__)
							{
								//fprintf(stderr, "fore[%d][%d]: %.5f: P_emit(%d, %d)=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, fore_scores[var_i][state_i], cur_ref_allele, cur_test_allele, emit_prob, prev_state_i, state_i, trans_prob);
							}
						} // prev_state_i loop.
					} // state_i loop.
				} // var_i loop.

				// Compute the total log forward and backward probabilities.
				double total_log_fore_prob = xlog(0.0);
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					total_log_fore_prob = xlog_sum(total_log_fore_prob, fore_scores[n_vars + 1][state_i]);
				} // state_i loop.

				fprintf(stderr, "Total probabilities: fore=%.5f\n", total_log_fore_prob);

				// Save the probability for this sample.
				per_sample_per_hap_emission_probs->at(test_sample_i)[test_hap_i] = total_log_fore_prob;
			} // test_hap_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	return(NULL);
}

void calculate_haplotype_emission_probabilities_per_reference_Full_States_multithreaded(char* reference_tag_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* recombination_rate_dir,
	double min_tar2tag_cM,
	double N_e_2_n_ref_haplotypes,
	double allele_err,
	int n_threads,
	char* fore_back_output_fp)
{
	//fprintf(stderr, "Running fore-back with minimum tar2tag distance of %.3f cMs\n", min_tar2tag_cM);
	fprintf(stderr, "Calculating the haplotype emission probabilities.\n");

	vector<t_annot_region*>* reference_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		vecsize(reference_haplocoded_tag_geno_regs),
		vecsize(ref_sample_ids),
		n_ref_haplotypes);

	double N_e = N_e_2_n_ref_haplotypes * n_ref_haplotypes;

	fprintf(stderr, "Calculating LL using %d reference haplotypes using N_e=%.3f, allele_eps=%.3f.\n", (int)(n_ref_haplotypes), N_e, allele_err);

	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", vecsize(testing_haplocoded_tag_geno_regs), vecsize(testing_sample_ids));

	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", vecsize(intersects));
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	//int n_states = n_ref_haplotypes;
	//int n_symbols = 2;

	//double N_e = pow(10, 6);
	//double allele_err = pow(10, -4);

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	// Assign the recombination rates.
	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Setting the recombination rates on %s.\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Set the start and end state regions.
		t_annot_region* start_state_reg = get_empty_region();
		start_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		start_state_reg->start = 1;
		start_state_reg->end = 1;
		start_state_reg->strand = '+';
		start_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->at(0)->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->insert(cur_chr_testing_haplocoded_tag_geno_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		//int n_vars = vecsize(cur_chr_testing_haplocoded_tag_geno_regs) - 1;

		// Add the end state region.
		t_annot_region* end_state_reg = get_empty_region();
		end_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		end_state_reg->start = 1;
		end_state_reg->end = 1;
		end_state_reg->strand = '+';
		end_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->back()->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->push_back(end_state_reg);
	} // i_Chr loop.

	// Set the emission probabilities.
	vector<double*>* per_sample_per_hap_emission_probs = new vector<double*>();
	for (int i_s = 0; i_s < vecsize(testing_sample_ids); i_s++)
	{
		double* cur_subj_per_allele_probs = new double[2];
		cur_subj_per_allele_probs[0] = xlog(0);
		cur_subj_per_allele_probs[1] = xlog(0);
		per_sample_per_hap_emission_probs->push_back(cur_subj_per_allele_probs);
	} // i_s loop.

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting %d. forward variable thread..\n", thread_i);

		void** thread_ptrs_list = new void* [12];

		int* which_out_of_ptr = new int[10];
		which_out_of_ptr[0] = thread_i;
		which_out_of_ptr[1] = n_threads;
		thread_ptrs_list[0] = which_out_of_ptr;

		//t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = NULL;
		thread_ptrs_list[1] = restr_testing_haplocoded_tag_geno_regs;

		//vector<char*>* testing_sample_ids = NULL;
		thread_ptrs_list[2] = testing_sample_ids;

		//t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = NULL;
		thread_ptrs_list[3] = restr_reference_haplocoded_tag_target_geno_regs;

		//vector<char*>* ref_sample_ids = NULL;
		thread_ptrs_list[4] = ref_sample_ids;

		//char* recombination_rate_dir = NULL;
		thread_ptrs_list[5] = recombination_rate_dir;

		//double* N_e_allele_err_ptr = NULL;
		double* N_e_allele_err_ptr = new double[5];
		N_e_allele_err_ptr[0] = N_e;
		N_e_allele_err_ptr[1] = allele_err;
		thread_ptrs_list[6] = N_e_allele_err_ptr;

		thread_ptrs_list[7] = per_sample_per_hap_emission_probs;

		t_ansi_thread* cur_thread = new t_ansi_thread(emission_probability_computing_thread_callback, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Started %d/%d threads; waiting.\n", vecsize(resampling_threads), n_threads);

	for (int thread_i = 0; thread_i < vecsize(resampling_threads); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		fprintf(stderr, "%d. thread finished.\n", thread_i);
	} // thread_i waiting loop. 

	fprintf(stderr, "Saving results to %s\n", fore_back_output_fp);
	FILE* f_per_test_subj_emission_prob = open_f(fore_back_output_fp, "w");
	for (int i_s = 0; i_s < vecsize(testing_sample_ids); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			fprintf(f_per_test_subj_emission_prob, "%s\t%d\t%d\t%.6f\n", testing_sample_ids->at(i_s), i_s, i_hap, per_sample_per_hap_emission_probs->at(i_s)[i_hap]);
		} // i_hap loop.
	} // i_s loop.
	close_f(f_per_test_subj_emission_prob, fore_back_output_fp);

	fprintf(stderr, "Done!\n");
}

void calculate_haplotype_emission_probabilities_per_reference_Full_States(char* reference_tag_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* recombination_rate_dir,
	double min_tar2tag_cM,
	double N_e,
	double allele_err,
	char* fore_back_output_fp)
{
	//fprintf(stderr, "Running fore-back with minimum tar2tag distance of %.3f cMs\n", min_tar2tag_cM);
	fprintf(stderr, "Calculating the haplotype emission probabilities.\n");

	vector<t_annot_region*>* reference_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = vecsize(ref_sample_ids) * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		vecsize(reference_haplocoded_tag_geno_regs),
		vecsize(ref_sample_ids),
		n_ref_haplotypes);

	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", vecsize(testing_haplocoded_tag_geno_regs), vecsize(testing_sample_ids));

	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", vecsize(intersects));
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	int n_states = n_ref_haplotypes;
	//int n_symbols = 2;

	//double N_e = pow(10, 6);
	//double allele_err = pow(10, -4);

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	for (int i_chr = 0; i_chr < vecsize(restr_testing_haplocoded_tag_geno_regs->chr_ids); i_chr++)
	{
		fprintf(stderr, "Running forward on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
		for (int i_reg = 0; i_reg < vecsize(cur_chr_testing_haplocoded_tag_geno_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < vecsize(cur_chrom_ref_tag_target_var_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.		

		// Set the start and end state regions.
		t_annot_region* start_state_reg = get_empty_region();
		start_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		start_state_reg->start = 1;
		start_state_reg->end = 1;
		start_state_reg->strand = '+';
		start_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->at(0)->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->insert(cur_chr_testing_haplocoded_tag_geno_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		int n_vars = vecsize(cur_chr_testing_haplocoded_tag_geno_regs) - 1;

		// Add the end state region.
		t_annot_region* end_state_reg = get_empty_region();
		end_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		end_state_reg->start = 1;
		end_state_reg->end = 1;
		end_state_reg->strand = '+';
		end_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->back()->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->push_back(end_state_reg);

		for (int test_sample_i = 0; test_sample_i < vecsize(testing_sample_ids); test_sample_i++)
		{
			fprintf(stderr, "Computing forward probabilities for sample %d\n", test_sample_i);

			for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
			{
				fprintf(stderr, "Computing fore-ward scores for sample %d; hap: %d\n", test_sample_i, test_hap_i);

				// Allocate and initialize the forward array.
				double** fore_scores = new double* [n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					fore_scores[var_i] = new double[n_states + 2];
					memset(fore_scores[var_i], 0, sizeof(double) * (n_states + 1));
				} // i loop.

				// Initialize the state probabilities.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					fore_scores[0][state_i] = xlog((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				// Start recursing over the variants.
				for (int var_i = 1; var_i <= n_vars + 1; var_i++)
				{
					if (var_i % 100 == 0)
					{
						fprintf(stderr, "Forward: sample_i: %d/%d [%d]: var_i: %d         \r", test_sample_i, vecsize(testing_sample_ids), test_hap_i, var_i);
					}

					void** cur_tag_var_info = (void**)(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->data);

					char* test_sample_geno = NULL;
					char* ref_sample_geno = NULL;
					if (var_i <= n_vars)
					{
						test_sample_geno = (char*)(cur_tag_var_info[0]);
						ref_sample_geno = (char*)(cur_tag_var_info[1]);
					}

					/////////////////////////////////////////////////////
					// Pre-compute the transition probabilities:
					double prev_var_cM = 0;
					if (var_i > 1)
					{
						prev_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i - 1)->dbl_score;
					}

					double cur_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score;

					double r_m = fabs(cur_var_cM - prev_var_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

					//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
					/////////////////////////////////////////////////////

					// Loop over all the states.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						// Recurse over the previous states.
						fore_scores[var_i][state_i] = xlog(0);
						for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
						{
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the transition probabilities.
							double trans_prob = 0;
							if (state_i == prev_state_i)
							{
								trans_prob = self_prob;
							}
							else
							{
								trans_prob = other_prob;
							}

							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
							double emit_prob = 0;

							if (var_i == n_vars + 1)
							{
								emit_prob = 1.0;
							}
							else
							{
								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);
								if (cur_ref_allele == cur_test_allele)
								{
									emit_prob = 1 - allele_err;
								}
								else
								{
									emit_prob = allele_err;
								}
							}
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

							double trans_emit_prob = xlog_mul(xlog(trans_prob), xlog(emit_prob));

							fore_scores[var_i][state_i] = xlog_sum(fore_scores[var_i][state_i], xlog_mul(trans_emit_prob, fore_scores[var_i - 1][prev_state_i]));
							
							if (__DUMP_PROXYTYPER_MSGS__)
							{
								//fprintf(stderr, "fore[%d][%d]: %.5f: P_emit(%d, %d)=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, fore_scores[var_i][state_i], cur_ref_allele, cur_test_allele, emit_prob, prev_state_i, state_i, trans_prob);
							}
						} // prev_state_i loop.
					} // state_i loop.
				} // var_i loop.

				// Compute the total log forward and backward probabilities.
				double total_log_fore_prob = xlog(0.0);
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					total_log_fore_prob = xlog_sum(total_log_fore_prob, fore_scores[n_vars + 1][state_i]);
				} // state_i loop.
				fprintf(stderr, "Total probabilities: fore=%.5f\n", total_log_fore_prob);
			} // test_hap_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	fprintf(stderr, "Done!\n");
}
