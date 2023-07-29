#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "prxytypr_file_utils.h"
#include "prxytypr_vector_macros.h"
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
#include "prxytypr_proxytyper.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

const bool __DUMP_UNTYPED_PROXIZATION_MSGS__ = false;

void combine_BEAGLE_imputed_decomposed_genotype_probabilities(char* beagle_imputed_vcf_fp, 
												char* untyped_variant_decomposed_intervals_fp,
												char* panel_sample_list_fp,
												char* op_prefix)
{
	vector<t_annot_region*>* untyped_variant_prob_geno_regs = load_decomposed_variant_Intervals(untyped_variant_decomposed_intervals_fp);

	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);

	fprintf(stderr, "Combining genotype probabilities for %d variants on %d subjects.\n", 
			vecsize(untyped_variant_prob_geno_regs),
			vecsize(panel_sample_ids));

	// Setup the per haplotype alternate allele probability.
	for (int i_reg = 0; i_reg < vecsize(untyped_variant_prob_geno_regs); i_reg++)
	{
		// For the current variant, set the probabilities for each allele.
		double** per_hap_probs = new double*[2];
		per_hap_probs[0] = new double[vecsize(panel_sample_ids)];
		memset(per_hap_probs[0], 0, sizeof(double) * vecsize(panel_sample_ids));

		per_hap_probs[1] = new double[vecsize(panel_sample_ids)];
		memset(per_hap_probs[1], 0, sizeof(double) * vecsize(panel_sample_ids));

		untyped_variant_prob_geno_regs->at(i_reg)->data = per_hap_probs;
	} // i_reg loop.

	// Assign the decomposed regions back to original regions.
	vector<t_annot_region*>* decomposed_var_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < vecsize(untyped_variant_prob_geno_regs); i_reg++)
	{
		for (int j_reg = 0; j_reg < vecsize(untyped_variant_prob_geno_regs->at(i_reg)->intervals); j_reg++)
		{
			t_annot_region* cur_decomp_var_reg = untyped_variant_prob_geno_regs->at(i_reg)->intervals->at(j_reg);

			// Link the original variant to the decomp variant.
			cur_decomp_var_reg->data = untyped_variant_prob_geno_regs->at(i_reg);

			// Add the decomp variant to list of decomp variants to intersect with vcf regions.
			decomposed_var_regs->push_back(cur_decomp_var_reg);
		} // j_reg loop.
	} // i_reg loop.

	// Set sorting info for decomposed regions.
	sort_set_sorting_info(decomposed_var_regs, sort_regions);

	// Start loading and adding the probabilities for the original regions.
	fprintf(stderr, "Reading BEAGLE VCF @ %s\n", beagle_imputed_vcf_fp);
	FILE* f_vcf = open_f(beagle_imputed_vcf_fp, "r");
	char tok_buffer[1000];
	int n_processed_vars = 0;
	int n_unmatched_vcf_regs = 0;
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			char first_token[1000];
			sscanf(cur_line, "%s", first_token);
			if (t_string::compare_strings(first_token, "#CHROM"))
			{
				// Parse the subject identifiers.
				t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");

				if (vecsize(toks) != (vecsize(panel_sample_ids) + 9))
				{
					fprintf(stderr, "Subject ids are not matching.\n");
					exit(1);
				}

				for (int i_t = 9; i_t < vecsize(toks); i_t++)
				{
					int i_s = (i_t - 9);
					if (!t_string::compare_strings(panel_sample_ids->at(i_s), toks->at(i_t)->str()))
					{
						fprintf(stderr, "Subject ids are not matching.\n");
						exit(1);
					}
				} // i_t loop.

				fprintf(stderr, "Loaded matching subjects ids from VCF.\n");
			} // subject ids line check.
		} // comment check.
		else
		{
			n_processed_vars++;

			if (n_processed_vars % 1000 == 0)
			{
				fprintf(stderr, "@ %d. variant (%d unmatched)..           \r", n_processed_vars, n_unmatched_vcf_regs);
			}

			// This is a variant region.
			//22      16050275        rs376238049     C       T.PASS    DR2 = 0.00; AF = 0.0091; IMP  GT : DS:AP1:AP2:GP        0 | 0 : 0.02 : 0.01 : 0.01 : 0.98, 0.02, 0
			int var_posn = 0;
			char var_id[1000];
			char ref_all;
			char alt_all;

			int i_next_char = 0;
			for (int i_t = 0; i_t < 9; i_t++)
			{
				if (!t_string::get_next_token(cur_line, tok_buffer, 1000, "\t", i_next_char))
				{
					fprintf(stderr, "Could not read the variant info tokens.\n");
					exit(1);
				}

				if (i_t == 1)
				{
					var_posn = translate_coord(atoi(tok_buffer), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
				}
				else if (i_t == 2)
				{
					strcpy(var_id, tok_buffer);
				}
				else if (i_t == 3)
				{
					ref_all = tok_buffer[0];
				}
				else if (i_t == 4)
				{
					alt_all = tok_buffer[0];
				}
			} // i_t loop.

			ref_all += 1;
			alt_all += 1;

			////////////////////////////////
			// Locate the position of this variant.
			int decomp_var_i = locate_posn_region_per_region_starts(var_posn, decomposed_var_regs, 0, vecsize(decomposed_var_regs));
			while (decomp_var_i > 0 &&
				var_posn <= decomposed_var_regs->at(decomp_var_i)->start)
			{
				decomp_var_i--;
			} // var_i loop.

			bool found_decomp_reg = false;
			while (decomp_var_i < vecsize(decomposed_var_regs) &&
				var_posn >= decomposed_var_regs->at(decomp_var_i)->end)
			{
				if (var_posn == decomposed_var_regs->at(decomp_var_i)->start)
				{
					found_decomp_reg = true;
					break;
				}
				decomp_var_i++;
			} // decomp_var_i loop.
			////////////////////////////////

			if (found_decomp_reg)
			{
				// First, get the untyped variant region for this decomposed untyped variant.
				t_annot_region* untyped_var_reg = (t_annot_region*)(decomposed_var_regs->at(decomp_var_i)->data);
				int cur_decomp_var_bias = decomposed_var_regs->at(decomp_var_i)->score;

				// Get the per-haplotype probabilities.
				double** per_hap_probs = (double**)(untyped_var_reg->data);

				// Start loading the remaining tokens as genotypes.
				for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
				{
					if (!t_string::get_next_token(cur_line, tok_buffer, 1000, "\t", i_next_char))
					{
						fprintf(stderr, "Could not read the variant genotype tokens for %s @ %d. subject.\n", var_id, i_s);
						exit(1);
					}

					// 0|0:0.02:0.01:0.01:0.98,0.02,0
					int cur_geno_char_i = 0;
					char geno_tok_buffer[1000];
					char assigned_alleles[1000];
					double AP1 = 0;
					double AP2 = 0;
					t_string::get_next_token(tok_buffer, geno_tok_buffer, 1000, ":", cur_geno_char_i);
					strcpy(assigned_alleles, geno_tok_buffer);
					t_string::get_next_token(tok_buffer, geno_tok_buffer, 1000, ":", cur_geno_char_i);
					t_string::get_next_token(tok_buffer, geno_tok_buffer, 1000, ":", cur_geno_char_i);
					AP1 = atof(geno_tok_buffer);
					t_string::get_next_token(tok_buffer, geno_tok_buffer, 1000, ":", cur_geno_char_i);
					AP2 = atof(geno_tok_buffer);

					// Resolve 50%'s (0.5) by looking at the assigned alleles. This becomes necessary as BEAGLE does not report more than 2 significant digits of allele probability.
					if (AP1 == 0.5)
					{
						if (geno_tok_buffer[0] == '0')
						{
							AP1 = 0.49;
						}
						else
						{
							AP1 = 0.51;
						}
					}

					if (AP2 == 0.5)
					{
						if (geno_tok_buffer[2] == '0')
						{
							AP2 = 0.49;
						}
						else
						{
							AP2 = 0.51;
						}
					}

					// Switch with respect to bias.
					AP1 = (cur_decomp_var_bias == 0) ? (AP1) : (1 - AP1);
					AP2 = (cur_decomp_var_bias == 0) ? (AP2) : (1 - AP2);

					// Update the probabilities for this sample.
					// It is very important to code these correctly to be consistent with VCF loading function.
					// First allele in VCF is shifted 1 bit so below we assign AP1 to second position, which is 
					// coded with 2 bit shifting.
					per_hap_probs[0][i_s] += AP2;
					per_hap_probs[1][i_s] += AP1;
				} // i_s loop.
			} // check for finding the variant.
			else
			{
				n_unmatched_vcf_regs++;
			}
		} // comment check.

		delete[] cur_line;
	} // vcf reading loop.
	close_f(f_vcf, beagle_imputed_vcf_fp);

	vector<t_annot_region*>* untyped_var_hard_geno_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < vecsize(untyped_variant_prob_geno_regs); i_reg++)
	{
		t_annot_region* cur_var_reg = duplicate_region(untyped_variant_prob_geno_regs->at(i_reg));
		void** var_reg_info = new void* [5];
		char* var_geno_sig = new char[vecsize(panel_sample_ids) + 2];
		memset(var_geno_sig, 0, sizeof(char)* (vecsize(panel_sample_ids) + 2));
		var_reg_info[0] = var_geno_sig;
		cur_var_reg->data = var_reg_info;

		untyped_var_hard_geno_regs->push_back(cur_var_reg);

		double** per_hap_allele_probs = (double**)(untyped_variant_prob_geno_regs->at(i_reg)->data);

		// Hard decisions on the allele probabilities.
		for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
		{
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				if (per_hap_allele_probs[i_hap][i_s] >= 0.5)
				{
					var_geno_sig[i_s] = var_geno_sig[i_s] | (1 << i_hap);
				}
			} // i_hap loop.
		} // i_s loop.
	} // i_reg loop.

	char geno_sig_op_fp[1000];
	sprintf(geno_sig_op_fp, "%s_untyped_undecomp.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(untyped_var_hard_geno_regs, NULL, panel_sample_ids, geno_sig_op_fp);
} // combine_BEAGLE_imputed_decomposed_genotype_probabilities function.

vector<t_annot_region*>* load_decomposed_variant_Intervals(char* interval_fp)
{
	vector<t_annot_region*>* composite_interval_regions = new vector<t_annot_region*>();

	FILE* f_interval = open_f(interval_fp, "r");
	char* cur_line = getline(f_interval);
	while (cur_line != NULL)
	{
		bool skip_line = check_line_skip(cur_line);

		if (!skip_line)
		{
			// ENST00000469418.1|ENST00000480075.1     chr7    -       19757   35479   4       19757,20834,31060,35335 19895,21029,31606,35479
			char region_name[1000];
			char cur_chrom[1000];
			char cur_strand;
			int cur_start;
			int cur_end;
			int n_intervals;
			char int_starts[100000];
			char int_ends[100000];
			char scores[100000];
			char dbl_scores[100000];

			if (sscanf(cur_line, "%s %s %c %d %d %d %s %s %s %s", region_name, cur_chrom,
				&cur_strand, &cur_start, &cur_end, &n_intervals, int_starts, int_ends,
				scores, dbl_scores) == 10)
			{
				t_annot_region* new_int_comp_region = new t_annot_region();
				new_int_comp_region->chrom = t_string::copy_me_str(cur_chrom);
				normalize_chr_id(new_int_comp_region->chrom);
				new_int_comp_region->start = translate_coord(cur_start, INTERVAL_COORDS::start_base, CODEBASE_COORDS::start_base);
				new_int_comp_region->end = translate_coord(cur_end, INTERVAL_COORDS::end_base, CODEBASE_COORDS::end_base);
				new_int_comp_region->strand = cur_strand;
				new_int_comp_region->name = t_string::copy_me_str(region_name);
				new_int_comp_region->intervals = new vector<t_annot_region*>();

				// Allocate the interval.
				vector<t_annot_region*>* intervals = new vector<t_annot_region*>();
				t_string_tokens* starts_toks = t_string::tokenize_by_chars(int_starts, ",");
				t_string_tokens* ends_toks = t_string::tokenize_by_chars(int_ends, ",");
				t_string_tokens* scores_toks = t_string::tokenize_by_chars(scores, ",");
				t_string_tokens* dbl_scores_toks = t_string::tokenize_by_chars(dbl_scores, ",");

				if (starts_toks->size() != ends_toks->size() ||
					starts_toks->size() != scores_toks->size() ||
					starts_toks->size() != dbl_scores_toks->size())
				{
					fprintf(stderr, "The number of starts is not the same as ends in:\n%s\n", cur_line);
					exit(0);
				}

				// Parse all the intervals.
				for (int i_int = 0; i_int < (int)starts_toks->size(); i_int++)
				{
					t_annot_region* new_interval = new t_annot_region();
					new_interval->chrom = t_string::copy_me_str(cur_chrom);
					normalize_chr_id(new_interval->chrom);
					new_interval->start = translate_coord(atoi(starts_toks->at(i_int)->str()), INTERVAL_COORDS::start_base, CODEBASE_COORDS::start_base);
					new_interval->end = translate_coord(atoi(ends_toks->at(i_int)->str()), INTERVAL_COORDS::end_base, CODEBASE_COORDS::end_base);
					new_interval->score = atoi(scores_toks->at(i_int)->str());
					new_interval->dbl_score = atof(dbl_scores_toks->at(i_int)->str());

					if (new_interval->end < new_interval->start)
					{
						fprintf(stderr, "One of the starts is larger than the end for:\n%s\n", cur_line);
						exit(1);
					}

					if (new_interval->score > 1 ||
						new_interval->score < 0)
					{
						fprintf(stderr, "One of the scores is invalid:\n%s\n", cur_line);
						exit(1);
					}


					if (new_interval->dbl_score > 1 ||
						new_interval->dbl_score < 0)
					{
						fprintf(stderr, "One of the dbl_scores is invalid:\n%s\n", cur_line);
						exit(1);
					}

					new_interval->strand = cur_strand;

					intervals->push_back(new_interval);
				} // i_int loop.

				new_int_comp_region->intervals = intervals;

				composite_interval_regions->push_back(new_int_comp_region);
			}
			else
			{
				fprintf(stderr, "Could not parse:\n%s\n", cur_line);
				exit(0);
			}
		} // skip_line check.

		delete[] cur_line;

		// Read the next line.
		cur_line = getline(f_interval);
	}
	fclose(f_interval);

	return(composite_interval_regions);
} // load_decomposed_variant_Intervals function.

void dump_decomposed_variant_Intervals(char* interval_fp, vector<t_annot_region*>* annot_regions)
{
	FILE* f_Interval = open_f(interval_fp, "w");

	//00021   tars = readTarsFromBedFile ("-");
	//00022   for (i = 0; i < arrayMax (tars); i++) {
	//00023     currTar = arrp (tars,i,Tar);
	//00024     printf ("BED_%d\t%s\t.\t%d\t%d\t1\t%d\t%d\n",
	//00025             i + 1,
						//currTar->targetName,
						//currTar->start,
						//currTar->end,
						//currTar->start,
						//currTar->end);
	//00026   }
	//00027   return 0;
	//00028 }

	for (int i_reg = 0; i_reg < (int)annot_regions->size(); i_reg++)
	{
		/*
		1.   Name of the interval
		2.   Chromosome
		3.   Strand
		4.   Interval start (with respect to the "+")
		5.   Interval end (with respect to the "+")
		6.   Number of sub-intervals
		7.   Sub-interval starts (with respect to the "+", comma-delimited)
		8.   Sub-interval end (with respect to the "+", comma-delimited)
		*/
		char default_name_buffer[100];
		char* cur_int_name = NULL;
		if (annot_regions->at(i_reg)->name == NULL ||
			strcmp(annot_regions->at(i_reg)->name, ".") == 0)
		{
			sprintf(default_name_buffer, "Interval_%d", i_reg);
			cur_int_name = default_name_buffer;
		}
		else
		{
			cur_int_name = annot_regions->at(i_reg)->name;
		}

		// There needs to be interval information, if there is not, add it.
		if (annot_regions->at(i_reg)->intervals == NULL)
		{
			annot_regions->at(i_reg)->intervals = new vector<t_annot_region*>();
			annot_regions->at(i_reg)->intervals->push_back(duplicate_region(annot_regions->at(i_reg)));
		}

		//fprintf(f_Interval, "Interval_%d\t%s\t%c\t%d\t%d\t%d\t", 
		fprintf(f_Interval, "%s\t%s\t%c\t%d\t%d\t%d\t",
			cur_int_name,
			annot_regions->at(i_reg)->chrom,
			annot_regions->at(i_reg)->strand,
			translate_coord(annot_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, INTERVAL_COORDS::start_base),
			translate_coord(annot_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, INTERVAL_COORDS::end_base),
			(int)annot_regions->at(i_reg)->intervals->size());

		int i_e = 0;
		fprintf(f_Interval, "%d",
			translate_coord(annot_regions->at(i_reg)->intervals->at(i_e)->start, CODEBASE_COORDS::start_base, INTERVAL_COORDS::start_base));
		for (i_e = 1; i_e < (int)annot_regions->at(i_reg)->intervals->size(); i_e++)
		{
			fprintf(f_Interval, ",%d",
				translate_coord(annot_regions->at(i_reg)->intervals->at(i_e)->start, CODEBASE_COORDS::start_base, INTERVAL_COORDS::start_base));
		} // i_e loop.

		fprintf(f_Interval, "\t");

		i_e = 0;
		fprintf(f_Interval, "%d",
			translate_coord(annot_regions->at(i_reg)->intervals->at(i_e)->end, CODEBASE_COORDS::end_base, INTERVAL_COORDS::end_base));
		for (i_e = 1; i_e < (int)annot_regions->at(i_reg)->intervals->size(); i_e++)
		{
			//fprintf(f_Interval, ",%d", annot_regions->at(i_reg)->intervals->at(i_e)->end);
			fprintf(f_Interval, ",%d",
				translate_coord(annot_regions->at(i_reg)->intervals->at(i_e)->end, CODEBASE_COORDS::end_base, INTERVAL_COORDS::end_base));
		} // i_e loop.

		fprintf(f_Interval, "\t");

		// Save the scores.
		i_e = 0;
		fprintf(f_Interval, "%d", annot_regions->at(i_reg)->intervals->at(i_e)->score);
		for (i_e = 1; i_e < (int)annot_regions->at(i_reg)->intervals->size(); i_e++)
		{
			fprintf(f_Interval, ",%d", annot_regions->at(i_reg)->intervals->at(i_e)->score);
		} // i_e loop.

		fprintf(f_Interval, "\t");

		// Save the dbl_scores.
		i_e = 0;
		fprintf(f_Interval, "%.3f", annot_regions->at(i_reg)->intervals->at(i_e)->dbl_score);
		for (i_e = 1; i_e < (int)annot_regions->at(i_reg)->intervals->size(); i_e++)
		{
			fprintf(f_Interval, ",%.3f", annot_regions->at(i_reg)->intervals->at(i_e)->dbl_score);
		} // i_e loop.

		fprintf(f_Interval, "\n");
	} // i_reg loop.

	fclose(f_Interval);
} // dump_decomposed_variant_Intervals function.


// TODO::Add a bias (i.e., allele flip) to each decomposed variant.
// TODO::Add multiple variant decomposition for each variant.
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4716681/pdf/main.pdf
void simple_decompose_untyped_variants(char* typed_var_geno_sig_regs_fp, char* untyped_var_geno_sig_regs_fp, char* panel_sample_list_fp, double min_AAF_per_decomp_var, 
	bool shuffle_decomp_variants, char* op_prefix)
{
	fprintf(stderr, "Decomposing untyped variants with AAF > %.4f.\n", min_AAF_per_decomp_var);

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	vector<t_annot_region*>* typed_geno_sig_regs = load_variant_signal_regions_wrapper(typed_var_geno_sig_regs_fp, panel_sample_list_fp);
	vector<t_annot_region*>* untyped_geno_sig_regs = load_variant_signal_regions_wrapper(untyped_var_geno_sig_regs_fp, panel_sample_list_fp);
	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);

	vector<t_annot_region*>* all_var_geno_regs = new vector<t_annot_region*>();
	all_var_geno_regs->insert(all_var_geno_regs->end(), typed_geno_sig_regs->begin(), typed_geno_sig_regs->end());
	all_var_geno_regs->insert(all_var_geno_regs->end(), untyped_geno_sig_regs->begin(), untyped_geno_sig_regs->end());

	sort(all_var_geno_regs->begin(), all_var_geno_regs->end(), sort_regions);
	fprintf(stderr, "Processing %d typed+untyped variants per typed-typed blocks.\n", vecsize(all_var_geno_regs));

	for (int i_reg = 0; i_reg < vecsize(typed_geno_sig_regs); i_reg++)
	{
		typed_geno_sig_regs->at(i_reg)->score = 1;
	} // i_reg loop.

	int n_decomposed_vars = 0;
	for (int i_reg = 0; i_reg < vecsize(untyped_geno_sig_regs); i_reg++)
	{
		untyped_geno_sig_regs->at(i_reg)->score = 0;

		void** cur_reg_info = (void**)(untyped_geno_sig_regs->at(i_reg)->data);
		char* cur_geno_sig = (char*)(cur_reg_info[0]);

		double total_AF = 0;
		for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
		{
			total_AF += cur_geno_sig[i_s];
		}

		total_AF = total_AF / (2 * vecsize(panel_sample_ids));

		untyped_geno_sig_regs->at(i_reg)->dbl_score = total_AF;
	} // i_reg loop.

	dump_BED("typed_untyped_variants.bed", all_var_geno_regs);

	int i_reg = 0;
	vector<t_annot_region*>* cur_untyped_block = new vector<t_annot_region*>();
	int i_block = 0;
	t_annot_region* start_typed_var_reg = NULL;
	t_annot_region* end_typed_var_reg = NULL;

	int n_interval_too_short = 0;

	// Each block starts here; we may have reached the end of all variants, check it here.
	// Otherwise, we have a typed followed by untyped variants.
	while (i_reg < vecsize(all_var_geno_regs))
	{
		// Skip over the typed variants.
		while (i_reg < vecsize(all_var_geno_regs) &&
			all_var_geno_regs->at(i_reg)->score == 1)
		{
			start_typed_var_reg = all_var_geno_regs->at(i_reg);
			i_reg++;
		}

		// Fill the next block of untyped variants.
		while (i_reg < vecsize(all_var_geno_regs) &&
			all_var_geno_regs->at(i_reg)->score == 0)
		{
			cur_untyped_block->push_back(all_var_geno_regs->at(i_reg));
			i_reg++;
		}

		end_typed_var_reg = NULL;
		if (i_reg < vecsize(all_var_geno_regs))
		{
			end_typed_var_reg = all_var_geno_regs->at(i_reg);
		}

		if (i_reg < vecsize(all_var_geno_regs) &&
			all_var_geno_regs->at(i_reg)->score != 1)
		{
			fprintf(stderr, "The end typed var does not have correct id: %s:%d\n",
				all_var_geno_regs->at(i_reg)->chrom, all_var_geno_regs->at(i_reg)->start);

			exit(1);
		}

		if (cur_untyped_block->size() == 0)
		{
			continue;
		}

		sort(cur_untyped_block->begin(), cur_untyped_block->end(), sort_regions);

		int block_start_coord = cur_untyped_block->at(0)->start + 1;
		if (start_typed_var_reg != NULL)
		{
			block_start_coord = start_typed_var_reg->start + 1;
		}

		int block_end_coord = cur_untyped_block->back()->end - 1;
		if (end_typed_var_reg != NULL)
		{
			block_end_coord = end_typed_var_reg->end - 1;
		}

		// We divide the current block into intervals that we will add the decomposed variants into each interval (1 or 2).
		int per_var_interval = (block_end_coord - block_start_coord + 1) / (vecsize(cur_untyped_block) + 1);
		if (per_var_interval < 3)
		{
			fprintf(stderr, "Interval is very short @ Variant %d::%s:%d-%d (%d-%d) for %d variants.\n", 
				i_reg, 
				cur_untyped_block->at(0)->chrom, cur_untyped_block->at(0)->start, cur_untyped_block->back()->end,
				block_start_coord, block_end_coord,
				vecsize(cur_untyped_block));
			n_interval_too_short++;
//			exit(1);
		}

		// Find the current block of variants.
		fprintf(stderr, "@ block %d: %d variants [%s: %d-%d] (%d decomposed so far..)                    \r", i_block,
			vecsize(cur_untyped_block),
			untyped_geno_sig_regs->at(0)->chrom,
			block_start_coord,
			block_end_coord,
			n_decomposed_vars);

		vector<t_annot_region*>* cur_block_decomp_vars = new vector<t_annot_region*>();

		for (int i_reg = 0; i_reg < vecsize(cur_untyped_block); i_reg++)
		{
			// This is the normalized coordinate for this untyped variant in the block. We will generate decomposed variants around this.
			// Note that these coordinates are later shuffled.
			int cur_var_norm_coord = block_start_coord + (i_reg + 1) * per_var_interval;

			if (cur_untyped_block->at(i_reg)->dbl_score > min_AAF_per_decomp_var &&
				per_var_interval > 3)
			{
				n_decomposed_vars++;

				void** reg_info = (void**)(cur_untyped_block->at(i_reg)->data);
				char* geno_sig = (char*)(reg_info[0]);

				// Divide the alleles.
				char* left_geno_sig = new char[vecsize(panel_sample_ids) + 1];
				memset(left_geno_sig, 0, sizeof(char) * (vecsize(panel_sample_ids) + 1));
				char* right_geno_sig = new char[vecsize(panel_sample_ids) + 1];
				memset(right_geno_sig, 0, sizeof(char) * (vecsize(panel_sample_ids) + 1));

				// Select a uniform probability.
				double cur_var_left_decomp_prob = rng->random_double_ran3();

				// For each sample, randomly assign the allele to one of the decomposed variants.
				// We need to store, for each subject's haplotype, the allele that we assigned.
				for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
				{
					for (int i_hap = 0; i_hap < 2; i_hap++)
					{
						int cur_allele = get_allele_per_haplotype(geno_sig[i_s], i_hap);

						if (cur_allele == 1)
						{
							if (rng->random_double_ran3() < cur_var_left_decomp_prob)
							{
								// Assign to left
								left_geno_sig[i_s] = left_geno_sig[i_s] | (cur_allele << i_hap);
							}
							else
							{
								// Assign to right.
								right_geno_sig[i_s] = right_geno_sig[i_s] | (cur_allele << i_hap);
							}
						}
					} // i_hap loop.
				} // i_s loop.

				///////////////////////////////////////////////////////////////////////////////////////////////////
				// Add biases: This is done after partitioning the alleles to left and right variants.
				// Select the biases for left and right regions.
				int left_bias = (rng->random_double_ran3() > 0.5) ? (1) : (0);
				int right_bias = (rng->random_double_ran3() > 0.5) ? (1) : (0);

				char* bias_left_geno_sig = new char[vecsize(panel_sample_ids) + 1];
				memset(bias_left_geno_sig, 0, sizeof(char) * (vecsize(panel_sample_ids) + 1));
				char* bias_right_geno_sig = new char[vecsize(panel_sample_ids) + 1];
				memset(bias_right_geno_sig, 0, sizeof(char) * (vecsize(panel_sample_ids) + 1));

				// Add biases to both alleles.
				for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
				{
					for (int i_hap = 0; i_hap < 2; i_hap++)
					{
						int left_allele = get_allele_per_haplotype(left_geno_sig[i_s], i_hap);
						int bias_left_allele = (left_allele + left_bias) % 2;
						bias_left_geno_sig[i_s] = bias_left_geno_sig[i_s] | (bias_left_allele << i_hap);

						int right_allele = get_allele_per_haplotype(right_geno_sig[i_s], i_hap);
						int bias_right_allele = (right_allele + right_bias) % 2;
						bias_right_geno_sig[i_s] = bias_right_geno_sig[i_s] | (bias_right_allele << i_hap);
					} // i_hap loop.
				} // i_s loop.
				
				// Clear memory for non-biased geno signals.
				delete[] left_geno_sig;
				delete[] right_geno_sig;

				////////////////////////////////////////////////////////////////////////////////////////

				// Add two variants left and right that are 1/3rd of the per variant interval length of this variant.
				t_annot_region* left_reg = duplicate_region(cur_untyped_block->at(i_reg));
				//left_reg->start = MAX(block_start_coord, cur_untyped_block->at(i_reg)->start - 100);
				//left_reg->end = MAX(block_start_coord, cur_untyped_block->at(i_reg)->end - 100);
				left_reg->start = MAX(block_start_coord, cur_var_norm_coord - per_var_interval / 3);
				left_reg->end = MAX(block_start_coord, cur_var_norm_coord - per_var_interval / 3);
				left_reg->score = left_bias;
				left_reg->dbl_score = cur_var_left_decomp_prob;
				void** left_reg_info = new void* [5];
				left_reg_info[0] = bias_left_geno_sig;
				left_reg->data = left_reg_info;

				if (left_reg->name == NULL)
				{
					fprintf(stderr, "Sanity check failed.\n");
					exit(1);
				}

				char new_name[1000];
				sprintf(new_name, "left-%s", left_reg->name);
				delete[] left_reg->name;
				left_reg->name = t_string::copy_me_str(new_name);

				t_annot_region* right_reg = duplicate_region(cur_untyped_block->at(i_reg));
				//right_reg->start = MIN(block_end_coord, cur_untyped_block->at(i_reg)->start + 100);
				//right_reg->end = MIN(block_end_coord, cur_untyped_block->at(i_reg)->end + 100);
				right_reg->start = MAX(block_start_coord, cur_var_norm_coord + per_var_interval / 3);
				right_reg->end = MAX(block_start_coord, cur_var_norm_coord + per_var_interval / 3);
				right_reg->score = right_bias;
				right_reg->dbl_score = 1-cur_var_left_decomp_prob;
				void** right_reg_info = new void* [5];
				right_reg_info[0] = bias_right_geno_sig;
				right_reg->data = right_reg_info;

				if (left_reg->name == NULL)
				{
					fprintf(stderr, "Sanity check failed.\n");
					exit(1);
				}

				sprintf(new_name, "right-%s", right_reg->name);
				delete[] right_reg->name;
				right_reg->name = t_string::copy_me_str(new_name);

				vector<t_annot_region*>* cur_reg_intervals = new vector<t_annot_region*>();
				cur_reg_intervals->push_back(left_reg);
				cur_reg_intervals->push_back(right_reg);

				cur_untyped_block->at(i_reg)->intervals = cur_reg_intervals;

				cur_block_decomp_vars->push_back(left_reg);
				cur_block_decomp_vars->push_back(right_reg);
			} // AF check.
			else
			{
				void** reg_info = (void**)(cur_untyped_block->at(i_reg)->data);
				char* geno_sig = (char*)(reg_info[0]);

				int self_reg_bias = (rng->random_double_ran3() > 0.5) ? (1) : (0);

				char* self_reg_geno_sig = new char[vecsize(panel_sample_ids) + 1];
				memset(self_reg_geno_sig, 0, sizeof(char) * vecsize(panel_sample_ids));

				// Following adds the bias only for this variant.
				for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
				{
					for (int i_hap = 0; i_hap < 2; i_hap++)
					{
						int cur_allele = get_allele_per_haplotype(geno_sig[i_s], i_hap);

						int bias_allele = (cur_allele + self_reg_bias) % 2;

						self_reg_geno_sig[i_s] = self_reg_geno_sig[i_s] | (bias_allele << i_hap);
					} // i_hap loop.
				} // i_s loop.

				// This variant stays as is.
				t_annot_region* self_reg = duplicate_region(cur_untyped_block->at(i_reg));
				self_reg->start = cur_var_norm_coord;
				self_reg->end = cur_var_norm_coord;
				self_reg->score = self_reg_bias;
				self_reg->dbl_score = 1.0;
				//self_reg->data = cur_untyped_block->at(i_reg)->data;
				void** self_reg_info = new void* [5];
				self_reg_info[0] = self_reg_geno_sig;
				self_reg->data = self_reg_info;
				vector<t_annot_region*>* cur_reg_intervals = new vector<t_annot_region*>();
				//cur_reg_intervals->push_back(cur_untyped_block->at(i_reg));
				cur_reg_intervals->push_back(self_reg);
				cur_untyped_block->at(i_reg)->intervals = cur_reg_intervals;

				cur_block_decomp_vars->push_back(self_reg);
			}
		} // i_reg loop.

		// This should be a user specified parameter.
		if (shuffle_decomp_variants)
		{
			//////////////////////////////////////////////////////////////////
			// Shuffle all decomposed variants to remove their positional association with their parent variants.
			vector<int>* shuff_ind = rng->fast_permute_indices(0, cur_block_decomp_vars->size());
			vector<int>* shuff_var_start_pos = new vector<int>();
			vector<int>* shuff_var_end_pos = new vector<int>();
			for (int i_reg = 0; i_reg < vecsize(cur_block_decomp_vars); i_reg++)
			{
				int shuff_i_reg = shuff_ind->at(i_reg);
				shuff_var_start_pos->push_back(cur_block_decomp_vars->at(shuff_i_reg)->start);
				shuff_var_end_pos->push_back(cur_block_decomp_vars->at(shuff_i_reg)->end);
			} // i_reg loop.

			// Assign shuffled positions.
			for (int i_reg = 0; i_reg < vecsize(cur_block_decomp_vars); i_reg++)
			{
				cur_block_decomp_vars->at(i_reg)->start = shuff_var_start_pos->at(i_reg);
				cur_block_decomp_vars->at(i_reg)->end = shuff_var_end_pos->at(i_reg);
			} // i_reg loop.

			delete shuff_var_start_pos;
			delete shuff_var_end_pos;
			delete shuff_ind;
			delete cur_block_decomp_vars;
			//////////////////////////////////////////////////////////////////
		}

		i_block++;

		// Clear the current block for the next block.
		cur_untyped_block->clear();
	} // i_reg loop.
	fprintf(stderr, "Finished reading %d blocks (%d short intervals)\n", i_block, n_interval_too_short);

	// Save the decomposed genotypes and the interval file.
	fprintf(stderr, "Saving decomposing intervals list..\n");
	char decomposing_map_fp[1000];
	sprintf(decomposing_map_fp, "%s_untyped_decomposing.interval", op_prefix);

	// TODO::Replace this to include more information on the decomposed variants for each original variant.
	//dump_Interval(decomposing_map_fp, untyped_geno_sig_regs);
	dump_decomposed_variant_Intervals(decomposing_map_fp, untyped_geno_sig_regs);

	// Save the genotypes.
	vector<t_annot_region*>* decomposed_geno_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < vecsize(untyped_geno_sig_regs); i_reg++)
	{
		if (untyped_geno_sig_regs->at(i_reg)->intervals->size() == 0 ||
			untyped_geno_sig_regs->at(i_reg)->intervals->size() > 2)
		{
			fprintf(stderr, "Sanity check failed: %d decomposed variants @ %d. untyped variant with AF=%.4f.\n",
				vecsize(untyped_geno_sig_regs->at(i_reg)->intervals),
				i_reg,
				untyped_geno_sig_regs->at(i_reg)->dbl_score);
			exit(0);
		}

		decomposed_geno_regs->insert(decomposed_geno_regs->end(),
			untyped_geno_sig_regs->at(i_reg)->intervals->begin(),
			untyped_geno_sig_regs->at(i_reg)->intervals->end());
	} // i_reg loop.

	fprintf(stderr, "Checking uniqueness of decomposed untyped regions and the typed variants.\n");
	vector<t_annot_region*>* check_regs = new vector<t_annot_region*>();
	check_regs->insert(check_regs->end(), decomposed_geno_regs->begin(), decomposed_geno_regs->end());
	check_regs->insert(check_regs->end(), typed_geno_sig_regs->begin(), typed_geno_sig_regs->end());
	sort(check_regs->begin(), check_regs->end(), sort_regions);
	for (int i_reg = 1; i_reg < vecsize(check_regs); i_reg++)
	{
		if (check_regs->at(i_reg)->start == check_regs->at(i_reg - 1)->start)
		{
			fprintf(stderr, "Regions are not unique!\n");
			exit(1);
		}
	} // i_reg loop.
	fprintf(stderr, "Regions are unique!\n");

	fprintf(stderr, "Saving %d decomposed variant genotype regions out of %d original untyped variants.\n",
			vecsize(decomposed_geno_regs), vecsize(untyped_geno_sig_regs));

	char decomposed_geno_regs_fp[1000];
	sprintf(decomposed_geno_regs_fp, "%s_decomposed.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(decomposed_geno_regs, NULL, panel_sample_ids, decomposed_geno_regs_fp);

	fprintf(stderr, "Done!\n");
} // simple_decompose_untyped_variants function.

void get_untyped_variant_LD_statistics(char* typed_var_geno_sig_regs_fp, char* untyped_var_geno_sig_regs_fp, char* panel_sample_list_fp, int n_blocks_2_process, char* op_prefix)
{
	fprintf(stderr, "Generating LD statistics for %d blocks.\n", n_blocks_2_process);

	vector<t_annot_region*>* typed_geno_sig_regs = load_variant_signal_regions_wrapper(typed_var_geno_sig_regs_fp, panel_sample_list_fp);
	vector<t_annot_region*>* untyped_geno_sig_regs = load_variant_signal_regions_wrapper(untyped_var_geno_sig_regs_fp, panel_sample_list_fp);
	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);

	vector<t_annot_region*>* all_var_geno_regs = new vector<t_annot_region*>();
	all_var_geno_regs->insert(all_var_geno_regs->end(), typed_geno_sig_regs->begin(), typed_geno_sig_regs->end());
	all_var_geno_regs->insert(all_var_geno_regs->end(), untyped_geno_sig_regs->begin(), untyped_geno_sig_regs->end());

	sort(all_var_geno_regs->begin(), all_var_geno_regs->end(), sort_regions);
	fprintf(stderr, "Processing %d typed+untyped variants per typed-typed blocks.\n", vecsize(all_var_geno_regs));

	for (int i_reg = 0; i_reg < vecsize(typed_geno_sig_regs); i_reg++)
	{
		typed_geno_sig_regs->at(i_reg)->score = 1;
	} // i_reg loop.

	for (int i_reg = 0; i_reg < vecsize(untyped_geno_sig_regs); i_reg++)
	{
		untyped_geno_sig_regs->at(i_reg)->score = 0;

		void** cur_reg_info = (void**)(untyped_geno_sig_regs->at(i_reg)->data);
		char* cur_geno_sig = (char*)(cur_reg_info[0]);

		double total_AF = 0;
		for (int i_s = 0; i_s < vecsize(panel_sample_ids); i_s++)
		{
			total_AF += cur_geno_sig[i_s];
		}

		total_AF = total_AF / (2 * vecsize(panel_sample_ids));

		untyped_geno_sig_regs->at(i_reg)->dbl_score = total_AF;
	} // i_reg loop.

	int i_reg = 0;
	vector<t_annot_region*>* cur_untyped_block = new vector<t_annot_region*>();
	int i_block = 0;
	while (1)
	{
		// Find the current block of variants.
		fprintf(stderr, "Filling block %d\n", i_block);

		// Skip over the typed variants.
		while (i_reg < vecsize(all_var_geno_regs) &&
			all_var_geno_regs->at(i_reg)->score == 1)
		{
			i_reg++;
		}

		// Fill the next block of untyped variants.
		while (i_reg < vecsize(all_var_geno_regs) &&
			all_var_geno_regs->at(i_reg)->score == 0)
		{
			cur_untyped_block->push_back(all_var_geno_regs->at(i_reg));
			i_reg++;
		}

		fprintf(stderr, "Found a block of size %d\n", vecsize(cur_untyped_block));

		if (cur_untyped_block->size() == 0)
		{
			continue;
		}

		if (i_reg >= vecsize(all_var_geno_regs))
		{
			break;
		}

		sort(cur_untyped_block->begin(), cur_untyped_block->end(), sort_regions);

		fprintf(stderr, "block %d: %s:%d-%d [%d untyped vars]\n",
			i_block,
			cur_untyped_block->at(0)->chrom,
			cur_untyped_block->at(0)->start,
			cur_untyped_block->back()->end,
			vecsize(cur_untyped_block));

		// Calculate the pairwise LD in this block.
		char op_fp[1000];
		sprintf(op_fp, "%s_block_%d.txt", op_prefix, i_block);
		FILE* f_op = open_f(op_fp, "w");
		for (int i_reg = 0; i_reg < vecsize(cur_untyped_block); i_reg++)
		{
			if (i_reg % 100 == 0)
			{
				fprintf(stderr, "i_reg = %d        \r", i_reg);
			}

			//size_t min_j_reg = MAX(0, i_reg - l_block);
			for (int j_reg = 0; j_reg < vecsize(cur_untyped_block); j_reg++)
			{
				t_annot_region* cur_reg = cur_untyped_block->at(i_reg);
				void** cur_reg_info = (void**)(cur_reg->data);
				char* cur_reg_geno = (char*)(cur_reg_info[0]);
				t_annot_region* prev_reg = cur_untyped_block->at(j_reg);
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
				//fprintf(stderr, "%d[%.2f]\t%d[%.2f]\t%.4f\n",
				//	cur_untyped_block->at(i_reg)->start, 
				//	cur_untyped_block->at(i_reg)->dbl_score,
				//	cur_untyped_block->at(j_reg)->start, 
				//	cur_untyped_block->at(j_reg)->dbl_score,
				//	var2var_corr * var2var_corr);

				fprintf(f_op, "%d\t%d\t%.2f\t%d\t%d\t%.2f\t%.4f\n",
					i_reg,
					cur_untyped_block->at(i_reg)->start,
					cur_untyped_block->at(i_reg)->dbl_score,
					j_reg,
					cur_untyped_block->at(j_reg)->start,
					cur_untyped_block->at(j_reg)->dbl_score,
					var2var_corr * var2var_corr);

				//fprintf(stderr, "%d\t%d\t%.4f\n",
				//	(int)j_reg, (int)i_reg, var2var_corr * var2var_corr);

				delete[] cur_reg_geno_sig_dbl;
				delete[] prev_reg_geno_sig_dbl;
			} // j_reg loop.
		} // i_reg loop.
		close_f(f_op, NULL);

		//getc(stdin);
		i_block++;

		if (i_block >= n_blocks_2_process)
		{
			break;
		}

		// Clear the current block for the next block.
		cur_untyped_block->clear();
	} // i_reg loop.
} // get_untyped_variant_LD_statistics

vector<int>* locally_permute_indices(int n_elements, int n_vicinity, t_rng* rng)
{
	vector<int>* permuted_indices = new vector<int>();
	for (int i_el = 0; i_el < n_elements; i_el++)
	{
		permuted_indices->push_back(i_el);
	} // i_el loop.

	// Start permuting.
	for (int mid_i = 0; mid_i < (int)permuted_indices->size(); mid_i++)
	{
		int win_start_i = MAX(0, mid_i - n_vicinity);
		int win_end_i = MIN(n_elements - 1, mid_i + n_vicinity);

		int n_els = (win_end_i - win_start_i + 1);

		vector<int>* shuffled_local_rel_i = rng->permute_indices(n_els, n_els);
		for (int i_el = 0; i_el < n_els; i_el++)
		{
			int temp_i = permuted_indices->at(shuffled_local_rel_i->at(i_el) + win_start_i);

			permuted_indices->at(shuffled_local_rel_i->at(i_el) + win_start_i) = permuted_indices->at(i_el + win_start_i);
			permuted_indices->at(i_el + win_start_i) = temp_i;
		} // i_el loop.

		delete shuffled_local_rel_i;
	} // mid_i loop.

	return(permuted_indices);
}

void save_target_proxy_mapping_bed_per_anchoring_regions(vector<t_annot_region*>* panel_target_geno_sig_regs, char* target_proxy_mapping_bed_fp)
{
	FILE* f_bed = open_f(target_proxy_mapping_bed_fp, "w");
	if (f_bed == NULL)
	{
		fprintf(stderr, "Could not open %s for writing..\n", target_proxy_mapping_bed_fp);
		exit(1);
	}

	for (int i_reg = 0; i_reg < (int)panel_target_geno_sig_regs->size(); i_reg++)
	{
		if ((int)panel_target_geno_sig_regs->at(i_reg)->intervals->size() != 1)
		{
			fprintf(stderr, "%s:%d-%d::Sanity check failed: %d intervals..\n",
				panel_target_geno_sig_regs->at(i_reg)->chrom,
				panel_target_geno_sig_regs->at(i_reg)->start,
				panel_target_geno_sig_regs->at(i_reg)->end,
				(int)panel_target_geno_sig_regs->at(i_reg)->intervals->size());
			exit(1);
		}
		t_annot_region* cur_proxy_reg = panel_target_geno_sig_regs->at(i_reg)->intervals->at(0);

		void** proxy_reg_info = (void**)(cur_proxy_reg->data);
		int* allele_coding = (int*)(proxy_reg_info[1]);
		t_annot_region* permuted_variant_reg = (t_annot_region*)(proxy_reg_info[2]); // This is the actual variant whose alleles are stored.

		if (allele_coding[0] > 1 || allele_coding[0] < 0 ||
			allele_coding[1] > 1 || allele_coding[1] < 0)
		{
			fprintf(stderr, "Sanity check failed: Illegal allele coding in proxy target.\n");
			exit(0);
		}

		fprintf(f_bed, "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n",
			cur_proxy_reg->chrom,
			translate_coord(cur_proxy_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(cur_proxy_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			translate_coord(permuted_variant_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(permuted_variant_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			permuted_variant_reg->name,
			allele_coding[0], allele_coding[1]);
	} // i_reg loop.

	close_f(f_bed, target_proxy_mapping_bed_fp);
}

vector<t_annot_region*>* load_proxy_2_target_mapping_BED(char* proxy_2_target_mapping_BED_fp)
{
	vector<t_annot_region*>* proxy_2_target_mapping_regs = new vector<t_annot_region*>();

	FILE* f_proxy_2_target_map = open_f(proxy_2_target_mapping_BED_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_proxy_2_target_map);
		if (cur_line == NULL)
		{
			break;
		}

		/*fprintf(f_bed, "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n",
			cur_proxy_reg->chrom,
			translate_coord(cur_proxy_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(cur_proxy_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			translate_coord(permuted_variant_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(permuted_variant_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			permuted_variant_reg->name,
			allele_coding[0], allele_coding[1]);*/

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");
		if ((int)toks->size() != 8)
		{
			fprintf(stderr, "Could not correctly parse proxy-2-target mapping: %s (%d tokens)\n", cur_line, (int)toks->size());
			exit(0);
		}

		char* chrom = toks->at(0)->str();
		int proxy_start = atoi(toks->at(1)->str());
		int proxy_end = atoi(toks->at(2)->str());
		int perm_start = atoi(toks->at(3)->str());
		int perm_end = atoi(toks->at(4)->str());
		char* perm_name = toks->at(5)->str();
		int coding_0 = atoi(toks->at(6)->str());
		int coding_1 = atoi(toks->at(7)->str());

		t_annot_region* cur_proxy_reg = get_empty_region();
		cur_proxy_reg->chrom = t_string::copy_me_str(chrom);
		cur_proxy_reg->start = translate_coord(proxy_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_proxy_reg->end = translate_coord(proxy_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_proxy_reg->strand = '+';

		t_annot_region* cur_perm_reg = get_empty_region();
		cur_perm_reg->chrom = t_string::copy_me_str(chrom);
		cur_perm_reg->start = translate_coord(perm_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_perm_reg->end = translate_coord(perm_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_perm_reg->name = t_string::copy_me_str(perm_name);
		cur_perm_reg->strand = '+';

		int* coding_scheme = new int[2];
		coding_scheme[0] = coding_0;
		coding_scheme[1] = coding_1;

		void** proxy_reg_info = new void* [5];
		proxy_reg_info[0] = cur_perm_reg;
		proxy_reg_info[1] = coding_scheme;

		cur_proxy_reg->data = proxy_reg_info;

		proxy_2_target_mapping_regs->push_back(cur_proxy_reg);

		t_string::clean_tokens(toks);
	} // file reading loop.
	close_f(f_proxy_2_target_map, proxy_2_target_mapping_BED_fp);

	return(proxy_2_target_mapping_regs);
} // load_proxy_2_target_mapping_BED option.


void decode_untyped_variant_reference_panel_per_target_permutation(char* recoded_ref_panel_target_haplocoded_matbed_fp,
	char* proxy_2_target_mapping_BED_fp,
	char* panel_sample_list_fp,
	char* op_matbed_fp)
{
	vector<t_annot_region*>* proxy_var_geno_sig_regs = load_variant_signal_regions_wrapper(recoded_ref_panel_target_haplocoded_matbed_fp, panel_sample_list_fp);
	for (int i_var = 0; i_var < (int)proxy_var_geno_sig_regs->size(); i_var++)
	{
		void** var_reg_info = (void**)(proxy_var_geno_sig_regs->at(i_var)->data);
		void** new_var_reg_info = new void* [5];
		new_var_reg_info[0] = var_reg_info[0];
		new_var_reg_info[1] = NULL;
		new_var_reg_info[2] = NULL;
		new_var_reg_info[3] = NULL;
		new_var_reg_info[4] = NULL;

		proxy_var_geno_sig_regs->at(i_var)->data = new_var_reg_info;
	} // i_var loop.

	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);
	fprintf(stderr, "Loaded %d encoded proxy variant regions for %d subjects.\n", (int)proxy_var_geno_sig_regs->size(), (int)panel_sample_ids->size());

	vector<t_annot_region*>* proxy_2_target_mapping_regs = load_proxy_2_target_mapping_BED(proxy_2_target_mapping_BED_fp);
	fprintf(stderr, "Loaded %d proxy-target mapping regions..\n", (int)proxy_2_target_mapping_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(proxy_var_geno_sig_regs, proxy_2_target_mapping_regs, true);
	fprintf(stderr, "Processing %d intersects..\n", (int)intersects->size());
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* proxy_var_sig_reg = int_info->src_reg;
		t_annot_region* proxy_mapping_reg = int_info->dest_reg;

		void** var_sig_reg_info = (void**)(proxy_var_sig_reg->data);
		void** mapping_reg_info = (void**)(proxy_mapping_reg->data);

		// Copy the region information back.
		var_sig_reg_info[1] = mapping_reg_info[0]; // This is the permuted region.
		var_sig_reg_info[2] = mapping_reg_info[1]; // This is the allele re-coding.

		// We are mapping the permuted variant back to the signal region.
	} // i_int loop.

	// Go over the variant geno sig regions, replace the coordinates with permuted info, recode alleles.
	for (int i_var = 0; i_var < (int)proxy_var_geno_sig_regs->size(); i_var++)
	{
		void** cur_reg_info = (void**)(proxy_var_geno_sig_regs->at(i_var)->data);

		char* perm_var_geno_sig = (char*)(cur_reg_info[0]);
		t_annot_region* perm_var_reg = (t_annot_region*)(cur_reg_info[1]);
		int* allele_coding = (int*)(cur_reg_info[2]);

		if (perm_var_reg == NULL || allele_coding == NULL)
		{
			fprintf(stderr, "Sanity check failed: %s:%d-%d does not have proxy mapping information.\n",
				proxy_var_geno_sig_regs->at(i_var)->chrom, proxy_var_geno_sig_regs->at(i_var)->start, proxy_var_geno_sig_regs->at(i_var)->end);

			exit(0);
		}

		int* allele_decoder = new int[2];
		allele_decoder[allele_coding[0]] = 0;
		allele_decoder[allele_coding[1]] = 1;

		// Replace the coordinates.
		proxy_var_geno_sig_regs->at(i_var)->start = perm_var_reg->start;
		proxy_var_geno_sig_regs->at(i_var)->end = perm_var_reg->end;
		proxy_var_geno_sig_regs->at(i_var)->name = perm_var_reg->name;

		// Re-code the alleles:
		for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
		{
			int cur_coded_allele_0 = get_allele_per_haplotype(perm_var_geno_sig[i_s], 0);
			int cur_coded_allele_1 = get_allele_per_haplotype(perm_var_geno_sig[i_s], 1);

			// Decode the alleles.
			int decoded_allele_0 = allele_decoder[cur_coded_allele_0];
			int decoded_allele_1 = allele_decoder[cur_coded_allele_1];

			// Update the genotype:
			int decoded_geno = decoded_allele_0 + decoded_allele_1 * 2;
			perm_var_geno_sig[i_s] = decoded_geno;
		} // i_s loop.
	} // i_var loop.

	// Save.
	fprintf(stderr, "Saving %d decoded genotypes to %s\n", (int)proxy_var_geno_sig_regs->size(), op_matbed_fp);
	binarize_variant_genotype_signal_regions(proxy_var_geno_sig_regs, NULL, panel_sample_ids, op_matbed_fp);
}

void recode_untyped_variant_reference_panel_per_target_permutation(char* ref_panel_tag_haplocoded_matbed_fp,
	char* ref_panel_target_haplocoded_matbed_fp,
	char* panel_sample_list_fp,
	int untyped_var_perm_n_vicinity,
	double allele_switch_prob,
	char* op_prefix)
{
	// Load each target per the tag-tag block
	vector<t_annot_region*>* panel_tag_geno_sig_regs = load_variant_signal_regions_wrapper(ref_panel_tag_haplocoded_matbed_fp, panel_sample_list_fp);
	vector<t_annot_region*>* panel_target_geno_sig_regs = load_variant_signal_regions_wrapper(ref_panel_target_haplocoded_matbed_fp, panel_sample_list_fp);
	vector<char*>* panel_sample_ids = buffer_file(panel_sample_list_fp);
	fprintf(stderr, "Loaded %d tag and %d target regions for %d subjects..\n",
		(int)panel_tag_geno_sig_regs->size(),
		(int)panel_target_geno_sig_regs->size(),
		(int)panel_sample_ids->size());

	int max_geno = get_max_genotype_value(panel_tag_geno_sig_regs, panel_sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Seems like the tag panel is not phased.\n");
		exit(1);
	}

	max_geno = get_max_genotype_value(panel_target_geno_sig_regs, panel_sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Seems like the target panel is not phased.\n");
		exit(1);
	}

	fprintf(stderr, "Both panels are haplocoded..\n");

	int n_panel_haplotypes = 2 * (int)panel_sample_ids->size();

	fprintf(stderr, "%d panel haplotypes..\n", n_panel_haplotypes);

	for (int i_tag = 0; i_tag < (int)panel_tag_geno_sig_regs->size(); i_tag++)
	{
		panel_tag_geno_sig_regs->at(i_tag)->score = 0;
	} // i_tag loop.

	for (int i_target = 0; i_target < (int)panel_target_geno_sig_regs->size(); i_target++)
	{
		panel_target_geno_sig_regs->at(i_target)->score = 1;

		// Setup the intervals for the current target.
		panel_target_geno_sig_regs->at(i_target)->intervals = new vector<t_annot_region*>();
	} // i_target loop.

	//process_regions_per_callback(panel_tag_geno_sig_regs, [](t_annot_region* reg) {reg->score = 0; });
	//process_regions_per_callback(panel_target_geno_sig_regs, [](t_annot_region* reg) {reg->score = 1; });

	// Add tag and target regions to the list of all regions.
	vector<t_annot_region*>* tag_target_geno_sig_regs = new vector<t_annot_region*>();
	tag_target_geno_sig_regs->insert(tag_target_geno_sig_regs->end(), panel_tag_geno_sig_regs->begin(), panel_tag_geno_sig_regs->end());
	tag_target_geno_sig_regs->insert(tag_target_geno_sig_regs->end(), panel_target_geno_sig_regs->begin(), panel_target_geno_sig_regs->end());

	// Sort and restructure regions.
	t_restr_annot_region_list* restr_tag_target_geno_sig_regs = restructure_annot_regions(tag_target_geno_sig_regs);

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	//double allele_switch_prob = 0.5;

	// For each variant, we randomly pick a variant in the same block, use random variant's genotypes after recoding.
	// The coordinate of the variant is randomly shifted within the block.
	for (int i_chr = 0; i_chr < (int)restr_tag_target_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_tag_target_regs = restr_tag_target_geno_sig_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Processing %d tag/targets on %s\n", (int)cur_chr_tag_target_regs->size(),
			restr_tag_target_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_block_target_regs = new vector<t_annot_region*>();
		t_annot_region* cur_block_starting_tag_reg = NULL;
		int untyped_block_i = 0;
		int i_reg = 0;
		//for (int i_reg = 0; i_reg < cur_chr_tag_target_regs->size(); i_reg++)
		while (i_reg <= (int)cur_chr_tag_target_regs->size())
		{
			// Find the next tag region or the end of regions.
			if ((i_reg == (int)cur_chr_tag_target_regs->size() || cur_chr_tag_target_regs->at(i_reg)->score == 0))
			{
				t_annot_region* cur_block_ending_tag_reg = NULL;
				if (i_reg < (int)cur_chr_tag_target_regs->size())
				{
					cur_block_ending_tag_reg = cur_chr_tag_target_regs->at(i_reg);
				}

				// cur_block_starting_tag_reg==NULL means we are at the first block.
				// cur_block_ending_tag_reg==NULL means we are at the last block.

				// The blocks don't mean anything. The important point here are the start and end coordinates for the blocks, which are corrected at the ends of the chromosome.
				int cur_tag_block_start = -1;
				int cur_tag_block_end = -1;
				if (cur_block_starting_tag_reg == NULL)
				{
					// If we are here, we did not process any blocks yet, i.e., there were no blocks before this newly finished one.
					// This also means that there were no tag variants to the left of this.
					// So we set the starting position to 1kb left of first variant's position.
					// This ensures that we are to the left of any target variants. 
					// But note that there may not be any untyped variants in the first block
					// If that is the case, we should set the starting position correctly.
					//cur_tag_block_start = MAX(1, cur_block_target_regs->at(0)->start - 1000);
					cur_tag_block_start = MAX(1, cur_chr_tag_target_regs->at(0)->start - 1000);
				}
				else
				{
					cur_tag_block_start = cur_block_starting_tag_reg->start;
				}

				if (cur_block_ending_tag_reg != NULL)
				{
					cur_tag_block_end = cur_block_ending_tag_reg->end;
				}
				else
				{
					cur_tag_block_end = cur_chr_tag_target_regs->back()->end + 1000;
				}

				if (cur_block_starting_tag_reg != NULL &&
					cur_block_ending_tag_reg != NULL)
				{
					fprintf(stderr, "Found block: %s:%d-%d [%s] -- %s:%d-%d [%s]             \r",
						cur_block_starting_tag_reg->chrom,
						cur_block_starting_tag_reg->start, cur_block_starting_tag_reg->end,
						cur_block_starting_tag_reg->name,
						cur_block_ending_tag_reg->chrom,
						cur_block_ending_tag_reg->start, cur_block_ending_tag_reg->end,
						cur_block_ending_tag_reg->name);
				}

				//// Process the targets in the current block.
				//fprintf(stderr, "Processing %d target regions (%s:%d-%d).\n", cur_block_target_regs->size(), 
				//		cur_chr_tag_target_regs->at(i_reg)->chrom, cur_tag_block_start, cur_tag_block_end);

				// This check is done here after the block end is updated for the current block.
				if ((int)cur_block_target_regs->size() > 0)
				{
					// Permute its location to a new location within the same block.
					int n_target_cur_block = (int)cur_block_target_regs->size();
					//vector<int>* permuted_var_i = rng->fast_permute_indices(0, n_target_cur_block);
					vector<int>* permuted_var_i = locally_permute_indices(n_target_cur_block, untyped_var_perm_n_vicinity, rng);


					if ((int)permuted_var_i->size() != n_target_cur_block)
					{
						fprintf(stderr, "Sanity check failed: permuted variant size does not match.\n");
						exit(0);
					}

					// Store the target-target mappings (intervals)
					for (int i_bl_var = 0; i_bl_var < n_target_cur_block; i_bl_var++)
					{
						// Coordinates come from the original variant.
						t_annot_region* cur_var_perm_var_reg = duplicate_region(cur_block_target_regs->at(i_bl_var));

						//// Recode the coordinates: Change the coordinates for the permuted variant using the coordinates of the original variant.
						// TODO: Make sure all of the target variants are close to tag variants.
						// TODO: Make sure all of the target variants are close to tag variants.
						// TODO: Make sure all of the target variants are close to tag variants.
						double l_cur_block = (double)(cur_tag_block_end - cur_tag_block_start);
						//double coord_perm_win = coord_perm_win_frac * l_cur_block;
						//double l_cur_rand_coord_translation = (rng->random_double_ran3() - 0.5) * coord_perm_win;
						//cur_var_perm_var_reg->start = MAX(cur_tag_block_start, cur_block_target_regs->at(i_bl_var)->start + l_cur_rand_coord_translation);
						//cur_var_perm_var_reg->start = MIN(cur_tag_block_end, cur_var_perm_var_reg->start);
						//cur_var_perm_var_reg->end = cur_var_perm_var_reg->start;
						double delta_per_target = (int)(l_cur_block / (n_target_cur_block + 1));
						cur_var_perm_var_reg->start = cur_tag_block_start + delta_per_target * (i_bl_var + 1);
						cur_var_perm_var_reg->end = cur_tag_block_start + delta_per_target * (i_bl_var + 1);
						//cur_var_perm_var_reg->start = cur_block_target_regs->at(i_bl_var)->start + 10;
						//cur_var_perm_var_reg->end = cur_block_target_regs->at(i_bl_var)->start + 10;

						if (__DUMP_UNTYPED_PROXIZATION_MSGS__)
						{
							fprintf(stderr, "Target: %s:%d-%d::Permuted coords: %s:%d-%d within block [%d-%d]\n",
								cur_block_target_regs->at(i_bl_var)->chrom, cur_block_target_regs->at(i_bl_var)->start, cur_block_target_regs->at(i_bl_var)->end,
								cur_var_perm_var_reg->chrom, cur_var_perm_var_reg->start, cur_var_perm_var_reg->end,
								cur_tag_block_start, cur_tag_block_end);
						}

						// Recode the alleles for the permuted variant.
						void** perm_reg_info = new void* [5];
						int* recoded_alleles = new int[2];
						recoded_alleles[0] = 0; // Use original mapping of the alleles.
						if (rng->random_double_ran3() < allele_switch_prob)
						{
							recoded_alleles[0] = 1;
						}
						recoded_alleles[1] = 1 - recoded_alleles[0];

						// 1st entry is the recoded alleles.
						perm_reg_info[1] = recoded_alleles;

						// Recode and store the genotypes of the permuted variant to be set as this variant's genotypes.
						int perm_var_i = permuted_var_i->at(i_bl_var);
						void** cur_target_reg_info = (void**)(cur_block_target_regs->at(perm_var_i)->data);
						char* cur_target_var_geno_sig = (char*)(cur_target_reg_info[0]);
						char* cur_target_var_recoded_geno_sig = new char[panel_sample_ids->size()];

						for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
						{
							int cur_allele0 = get_allele_per_haplotype(cur_target_var_geno_sig[i_s], 0);
							int cur_recoded_allele0 = recoded_alleles[cur_allele0];

							int cur_allele1 = get_allele_per_haplotype(cur_target_var_geno_sig[i_s], 1);
							int cur_recoded_allele1 = recoded_alleles[cur_allele1];

							// Replace the recoded alleles on both haplotypes.
							cur_target_var_recoded_geno_sig[i_s] = cur_recoded_allele0 + (cur_recoded_allele1 * 2);
						} // i_s loop.

						// The 2nd entry is the permuted variant index.
						perm_reg_info[2] = cur_block_target_regs->at(perm_var_i);

						// 0th entry is the recoded genotypes.
						perm_reg_info[0] = cur_target_var_recoded_geno_sig;

						cur_var_perm_var_reg->data = perm_reg_info;

						// Add the permuted variant to the interval list for the original target variant.
						cur_block_target_regs->at(i_bl_var)->intervals->push_back(cur_var_perm_var_reg);
					} // i_bl_var loop.

					// Clear the target regions.
					cur_block_target_regs->clear();

					untyped_block_i++;
				} // target variant count check.

				// Update the block starting tag variant for next block.
				// Make sure we did not hit the end of the chromosome.
				if (i_reg < (int)cur_chr_tag_target_regs->size())
				{
					cur_block_starting_tag_reg = cur_chr_tag_target_regs->at(i_reg);
				}
			} // current block tag check.
			else if (cur_chr_tag_target_regs->at(i_reg)->score == 1)
			{
				// This is a target variant.
				cur_block_target_regs->push_back(cur_chr_tag_target_regs->at(i_reg));
			}

			// Check for the end of the blocks. This is handled at the beginning.
			if (i_reg == (int)cur_chr_tag_target_regs->size())
			{
				break;
			}

			// Move to the next region.
			i_reg++;
		} // i_tag/target reg loop
	} // i_chr loop.

	fprintf(stderr, "Saving proxized untyped variant mapping..\n");

	// Save the proxy information.
	char target_proxy_mapping_bed_fp[1000];
	sprintf(target_proxy_mapping_bed_fp, "%s_target_proxy_mapping.bed", op_prefix);
	save_target_proxy_mapping_bed_per_anchoring_regions(panel_target_geno_sig_regs, target_proxy_mapping_bed_fp);

	fprintf(stderr, "Saving proxized untyped genotypes..\n");
	// Save the genotypes using the proxized variants.
	vector<t_annot_region*>* proxized_target_geno_sig_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)panel_target_geno_sig_regs->size(); i_reg++)
	{
		vector<t_annot_region*>* cur_target_intervals = panel_target_geno_sig_regs->at(i_reg)->intervals;
		if ((int)cur_target_intervals->size() != 1)
		{
			fprintf(stderr, "%s:%d-%d::Sanity check failed: %d intervals..\n",
				panel_target_geno_sig_regs->at(i_reg)->chrom,
				panel_target_geno_sig_regs->at(i_reg)->start,
				panel_target_geno_sig_regs->at(i_reg)->end,
				(int)cur_target_intervals->size());
			exit(1);
		}

		proxized_target_geno_sig_regs->push_back(cur_target_intervals->at(0));
	} // i_reg loop.

	// Save the target genotypes.
	char proxized_target_geno_sig_regs_fp[1000];
	sprintf(proxized_target_geno_sig_regs_fp, "%s_proxized_targets.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(proxized_target_geno_sig_regs, NULL, panel_sample_ids, proxized_target_geno_sig_regs_fp);
}

// Dissect each target variant into 2,3,.. variants, depending on allele frequency.
void generate_untyped_variant_recoded_reference_panel_per_target_variant_dissection(char* ref_panel_tag_haplocoded_matbed_fp, char* ref_panel_target_haplocoded_matbed_fp, char* panel_sample_list_fp,
	char* op_dir, char* op_prefix)
{
	// Load each target per the tag-tag block.

	// Divide the target variant if it is above a certain AF into random targets.

	// Re-assign the genotypes.

	// Store target-target mappings (interval file)

	// Store the new genotypes for the dissected targets.
}