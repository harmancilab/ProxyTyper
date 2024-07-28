#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prxytypr_variation_tools.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_genome_sequence_tools.h"
#include "prxytypr_utils.h"
#include "prxytypr_histogram.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_nucleotide.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_vector_macros.h"
#include <string.h>
#include <ctype.h>

#include "prxytypr_nomenclature.h"
#include <zlib.h>
#include "prxytypr_ansi_thread.h"
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "prxytypr_imputation_utils.h"
#include "prxytypr_proxytyper.h"
#include "prxytypr_variation_tools.h"

#ifdef __unix__
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#endif 

#include <vector>
#include <algorithm>
using namespace std;

bool __DUMP_INPUTATION_UTILS_MSGS__ = false;
//
//struct t_LMSE_imputation_thread_info
//{
//	int which;
//	int outof;
//
//	double target_normalizer;
//	double tag_normalizer;
//
//	//vector<t_annot_region*>* testing_tag_genotype_signal_regs;
//	vector<t_annot_region*>* tag_genotype_signal_regs;
//	//vector<t_annot_region*>* testing_target_genotype_signal_regs;
//	vector<t_annot_region*>* target_genotype_signal_regs;
//
//	t_restr_annot_region_list* restr_target_genotype_signal_regs;
//	t_restr_annot_region_list* restr_tag_genotype_signal_regs;
//
//	int testing_flag;
//	int start_pos;
//	int end_pos;
//	char* op_dir;
//
//	vector<char*>* tag_sample_ids;
//	vector<char*>* target_sample_ids;
//
//	vector<char*>* testing_tag_sample_ids;
//	vector<char*>* testing_target_sample_ids;
//
//	int n_tags_vars_per_side;
//};

bool compare_haplotypes(double* haplo1, double* haplo2)
{
	int i_var = 0;
	while (haplo1[i_var] != -1)
	{
		if (haplo1[i_var] != haplo2[i_var])
		{
			return(false);
		}
		i_var++;
	}

	return(true);
}

bool sort_haplotypes_dbl(double* haplo1, double* haplo2)
{
	int i_var = 0;
	while (haplo1[i_var] >= 0)
	{
		if (haplo1[i_var] < haplo2[i_var])
		{
			return(true);
		}
		if (haplo1[i_var] > haplo2[i_var])
		{
			return(false);
		}

		i_var++;
	} // i_var loop.

	return(false);
}

void random_phase_het_genotypes(char* genotype_regs_BED_fp, char* sample_ids_list_fp, char* haplocoded_genotype_matrix_op_fp)
{
	fprintf(stderr, "Randomly phasing het variants that are in %s\n", genotype_regs_BED_fp);
	vector<t_annot_region*>* geno_signal_regs = load_variant_signal_regions_wrapper(genotype_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d variant region genotypes for %d samples.\n", (int)geno_signal_regs->size(), (int)sample_ids->size());
	for (int i_reg = 0; i_reg < (int)geno_signal_regs->size(); i_reg++)
	{
		void** reg_info = (void**)(geno_signal_regs->at(i_reg)->data);
		reg_info[1] = new char[(int)sample_ids->size() + 2];
	} // i_reg loop.

	int max_genotype = get_max_genotype_value(geno_signal_regs, sample_ids);
	if (max_genotype == 3)
	{
		fprintf(stderr, "Panel is phased.\n");
	}
	else if (max_genotype == 2)
	{
		fprintf(stderr, "Panel is unphased.\n");
	}
	else
	{
		fprintf(stderr, "%s(%d): Could not determine coding of %s\n", __FILE__, __LINE__, genotype_regs_BED_fp);
		exit(1);
	}

	fprintf(stderr, "Loaded %d variants for %d individuals.\n", (int)geno_signal_regs->size(), (int)sample_ids->size());

	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for (int i_reg = 0; i_reg < (int)geno_signal_regs->size(); i_reg++)
	{
		void** reg_info = (void**)(geno_signal_regs->at(i_reg)->data);
		char* geno_sig = (char*)(reg_info[0]);
		char* phased_geno_sig = (char*)(reg_info[1]);

		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			phased_geno_sig[i_s] = 0;
			char cur_geno = geno_sig[i_s];

			int allele0 = 0;
			int allele1 = 0;
			if (max_genotype == 3)
			{
				allele0 = get_allele_per_haplotype(cur_geno, 0);
				allele1 = get_allele_per_haplotype(cur_geno, 1);
			} // haplocoded check.
			else if (max_genotype == 2)
			{
				// Genocoded panel: 
				allele0 = 0;
				allele1 = 0;
				if (cur_geno == 2)
				{
					allele0 = 1;
					allele1 = 1;
				}
				else if (cur_geno == 1)
				{
					allele0 = 0;
					allele1 = 1;
				}
			} // genocoded check

			// Randomize the phase.
			if (rng->random_double_ran3() < 0.5)
			{
				phased_geno_sig[i_s] = allele0 + allele1 * 2;
			}
			else
			{
				phased_geno_sig[i_s] = allele1 + allele0 * 2;
			}
		} // i_s loop.
	} // i_reg loop.

	// Replace the genotype signal with random phased genotypes.
	for (int i_reg = 0; i_reg < (int)geno_signal_regs->size(); i_reg++)
	{
		void** reg_info = (void**)(geno_signal_regs->at(i_reg)->data);
		reg_info[0] = reg_info[1]; // replace the phased genotypes with the original unphased genotypes.
	} // i_reg loop.

	//binarize_variant_genotype_signal_regions(geno_signal_regs, NULL, sample_ids, haplocoded_genotype_matrix_op_fp);
	binarize_variant_signal_regions_wrapper(geno_signal_regs, sample_ids, haplocoded_genotype_matrix_op_fp);
} // random_phase_het_genotypes function.

bool sort_haplotype_info(void** haplo1_info, void** haplo2_info)
{
	double* haplo1 = (double*)(haplo1_info[0]);
	double* haplo2 = (double*)(haplo2_info[0]);

	int i_var = 0;
	while (haplo1[i_var] >= 0)
	{
		if (haplo1[i_var] < haplo2[i_var])
		{
			return(true);
		}
		if (haplo1[i_var] > haplo2[i_var])
		{
			return(false);
		}

		i_var++;
	} // i_var loop.

	return(false);
}

bool compare_haplotype_info(void** haplo1_info, void** haplo2_info)
{
	double* haplo1 = (double*)(haplo1_info[0]);
	double* haplo2 = (double*)(haplo2_info[0]);

	int i_var = 0;
	while (haplo1[i_var] != -1)
	{
		if (haplo1[i_var] != haplo2[i_var])
		{
			return(false);
		}
		i_var++;
	}

	return(true);
}

void get_unique_haplotype_indices(vector<double*>* haplotype_alleles, vector<vector<int>*>* per_uniq_haplo_indices)
{
	vector<void**>* haplo_info = new vector<void**>();
	for (int haplo_i = 0; haplo_i < (int)haplotype_alleles->size(); haplo_i++)
	{
		void** cur_haplo_info = new void*[2];
		cur_haplo_info[0] = haplotype_alleles->at(haplo_i);
		int* haplo_i_ptr = new int[1];
		haplo_i_ptr[0] = haplo_i;
		cur_haplo_info[1] = haplo_i_ptr;

		haplo_info->push_back(cur_haplo_info);
	} // haplo_i loop.

	sort(haplo_info->begin(), haplo_info->end(), sort_haplotype_info);

	if (__DUMP_INPUTATION_UTILS_MSGS__)
	{
		for (int hap_i = 0; hap_i < (int)haplo_info->size(); hap_i++)
		{
			double* cur_hap = (double*)(((void**)(haplo_info->at(hap_i)))[0]);

			fprintf(stderr, "Hap %d: ", hap_i);
			int var_i = 0;
			while (1)
			{
				if (cur_hap[var_i] == -1)
				{
					break;
				}
				else
				{
					fprintf(stderr, "%d", (int)(cur_hap[var_i]));
				}
				var_i++;
			}
			fprintf(stderr, "\n");
		} // hap_i loop.
	}

	int i_hap = 0;
	while (i_hap < (int)haplo_info->size())
	{
		int j_hap = i_hap;

		//int cur_hap_cnt = 0;
		vector<int>* cur_unique_haplo_indices = new vector<int>();
		while (j_hap < (int)haplo_info->size() &&
				compare_haplotype_info(haplo_info->at(i_hap), haplo_info->at(j_hap)))
		{
			void** cur_haplo_info = (void**)(haplo_info->at(j_hap));
			int* haplo_i_ptr = (int*)(cur_haplo_info[1]);
			cur_unique_haplo_indices->push_back(haplo_i_ptr[0]);
			j_hap++;
		} // j_hap loop.

		per_uniq_haplo_indices->push_back(cur_unique_haplo_indices);

		i_hap = j_hap;
	} // i_hap loop.

	//vector<void**>* haplo_info = new vector<void**>();
	for (int haplo_i = 0; haplo_i < (int)haplo_info->size(); haplo_i++)
	{
		void** cur_hap_info = haplo_info->at(haplo_i);
		delete[] ((int*)(cur_hap_info[1]));

		delete[] cur_hap_info;
	} // haplo_i loop.
	delete haplo_info;
}


void count_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes)
{
	sort(lowMAF_containing_haplotypes->begin(), lowMAF_containing_haplotypes->end(), sort_haplotypes_dbl);

	if (__DUMP_INPUTATION_UTILS_MSGS__)
	{
		for (int hap_i = 0; hap_i < (int)lowMAF_containing_haplotypes->size(); hap_i++)
		{
			fprintf(stderr, "Hap %d: ", hap_i);
			int var_i = 0;
			while (1)
			{
				if (lowMAF_containing_haplotypes->at(hap_i)[var_i] == -1)
				{
					break;
				}
				else
				{
					fprintf(stderr, "%d", (int)(lowMAF_containing_haplotypes->at(hap_i)[var_i]));
				}
				var_i++;
			}
			fprintf(stderr, "\n");
		} // hap_i loop.
	}

	int i_hap = 0;
	while (i_hap < (int)lowMAF_containing_haplotypes->size())
	{
		int j_hap = i_hap;

		int cur_hap_cnt = 0;
		while (j_hap < (int)lowMAF_containing_haplotypes->size() &&
			compare_haplotypes(lowMAF_containing_haplotypes->at(i_hap), lowMAF_containing_haplotypes->at(j_hap)))
		{
			cur_hap_cnt++;
			j_hap++;
		} // j_hap loop.

		n_cnt_per_uniq_haplotypes->push_back(cur_hap_cnt);

		i_hap = j_hap;
	} // i_hap loop.
}

vector<double*>* get_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes)
{
	vector<double*>* unique_lowMAF_containing_haplotypes = new vector<double*>();
	sort(lowMAF_containing_haplotypes->begin(), lowMAF_containing_haplotypes->end(), sort_haplotypes_dbl);

	int i_hap = 0;
	while(i_hap < (int)lowMAF_containing_haplotypes->size())
	{
		int j_hap = i_hap;		

		int cur_hap_cnt = 0;
		while (j_hap < (int)lowMAF_containing_haplotypes->size() &&
				compare_haplotypes(lowMAF_containing_haplotypes->at(i_hap), lowMAF_containing_haplotypes->at(j_hap)))
		{
			cur_hap_cnt++;
			j_hap++;
		} // j_hap loop.

		unique_lowMAF_containing_haplotypes->push_back(lowMAF_containing_haplotypes->at(i_hap));
		n_cnt_per_uniq_haplotypes->push_back(cur_hap_cnt);

		i_hap = j_hap;
	} // i_hap loop.

	return(unique_lowMAF_containing_haplotypes);
}

void is_haplo_coded(vector<t_annot_region*>* geno_regs, int n_samples, bool& found_haplocoded_genotype)
{
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(geno_regs->at(i_reg)->data);
		char* genos = (char*)(old_info[0]);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			if (genos[i_s] == 3)
			{
				found_haplocoded_genotype = true;
				return;
			}
		}
	} // i_reg loop.

	found_haplocoded_genotype = false;
	return;
}

vector<t_annot_region*>* load_recombination_rates(char* recombination_rate_fp)
{
	if (!check_file(recombination_rate_fp))
	{
		return NULL;
	}

	vector<t_annot_region*>* recomb_regs = new vector<t_annot_region*>();

	// plink formatted, no more headers.
	FILE* f_recombination_rate = open_f(recombination_rate_fp, "r");
	//char* header_line = getline(f_recombination_rate);
	//fprintf(stderr, "Read header: %s\n", header_line);
	while (1)
	{
		char* cur_line = getline(f_recombination_rate);

		if (cur_line == NULL)
		{
			break;
		}

		// 16052618 9.7078062447 0.0123386217370137
		char chrom[1000];
		char snp_id[1000];
		int cur_posn = 0;
		double cur_cM = 0;
		//if (sscanf(cur_line, "%d %*s %lf", &cur_posn, &cur_cM) != 2)
		if (sscanf(cur_line, "%s %s %lf %d", chrom, snp_id, &cur_cM, &cur_posn) != 4)
		{
			fprintf(stderr, "Could not parse line: %s\n", cur_line);
			exit(1);
		}

		t_annot_region* cur_recomb_reg = get_empty_region();
		cur_recomb_reg->chrom = t_string::copy_me_str(chrom);
		cur_recomb_reg->start = cur_posn;
		cur_recomb_reg->end = cur_posn;
		cur_recomb_reg->dbl_score = cur_cM;
		cur_recomb_reg->strand = '+';

		recomb_regs->push_back(cur_recomb_reg);

		if (recomb_regs->size() < 10)
		{
			fprintf(stderr, "%s => %s:%d @ cM: %.6f\n", cur_line, chrom, cur_posn, cur_cM);
		}

		delete[] cur_line;
	} // file reading loop.

	fprintf(stderr, "Loaded %d recombination regions.\n", (int)recomb_regs->size());
	close_f(f_recombination_rate, recombination_rate_fp);

	return(recomb_regs);
}

double get_avg_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs)
{
	// If the variant is to the left of the first recomb variant, return the first recomb rate.
	if (var_reg->end < sorted_recomb_regs->at(0)->start)
	{
		return(sorted_recomb_regs->at(0)->dbl_score);
	}

	// If the region is beyond the end of all the regions, return the last region.
	if (sorted_recomb_regs->back()->end < var_reg->start)
	{
		return(sorted_recomb_regs->back()->dbl_score);
	}

	int closest_reg_i = locate_posn_region_per_region_starts(var_reg->start, sorted_recomb_regs, 0, (int)sorted_recomb_regs->size());

	int i_reg = closest_reg_i;
	while (i_reg > 0 &&
		i_reg < (int)sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->end < var_reg->start)
		{
			break;
		}

		i_reg--;
	} // i_reg loop.

	//double interpolated_cM = -1;
	double averaged_recomb_rate = -1;
	while (i_reg >= 0 &&
		(i_reg + 1) < (int)sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->start > var_reg->end)
		{
			break;
		}

		if (sorted_recomb_regs->at(i_reg)->start <= var_reg->start &&
			sorted_recomb_regs->at(i_reg + 1)->start >= var_reg->start)
		{
			int l_recomb_reg = sorted_recomb_regs->at(i_reg + 1)->start - sorted_recomb_regs->at(i_reg)->start;
			double delta_cM = sorted_recomb_regs->at(i_reg + 1)->dbl_score - sorted_recomb_regs->at(i_reg)->dbl_score;

			/*double cM_slope = delta_cM / l_recomb_reg;

			double linear_interpolated_cM = sorted_recomb_regs->at(i_reg)->dbl_score + (var_reg->start - sorted_recomb_regs->at(i_reg)->start) * cM_slope;*/

			//// Get the distances and scale.
			//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
			//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

			//bool cur_reg_is_not_side = (sorted_recomb_regs->at(i_reg)->start != var_reg->start &&
			//	sorted_recomb_regs->at(i_reg + 1)->start != var_reg->start);

			//if (cur_reg_is_not_side &&
			//	(linear_interpolated_cM < MIN(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score) ||
			//		linear_interpolated_cM > MAX(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score)))
			//{
			//	fprintf(stderr, "Sanity check failed for interpolation.\n");
			//	exit(1);
			//}

			averaged_recomb_rate = delta_cM;
			averaged_recomb_rate /= l_recomb_reg;
			averaged_recomb_rate *= (1000 * 1000);

			break;
		} // overlap check.

		i_reg++;
	} // i_reg loop.

	return(averaged_recomb_rate);
}


double get_cumulative_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs)
{
	// If the variant is to the left of the first recomb variant, return the first recomb rate.
	if (var_reg->end < sorted_recomb_regs->at(0)->start)
	{
		return(sorted_recomb_regs->at(0)->dbl_score);
	}

	// If the region is beyond the end of all the regions, return the last region.
	if (sorted_recomb_regs->back()->end < var_reg->start)
	{
		return(sorted_recomb_regs->back()->dbl_score);
	}

	int closest_reg_i = locate_posn_region_per_region_starts(var_reg->start, sorted_recomb_regs, 0, (int)sorted_recomb_regs->size());
	
	int i_reg = closest_reg_i;
	while (i_reg > 0 && 
		i_reg < (int)sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->end < var_reg->start)
		{
			break;
		}

		i_reg--;
	} // i_reg loop.

	double interpolated_cM = -1;
	while (i_reg >= 0 &&
			(i_reg+1) < (int)sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->start > var_reg->end)
		{
			break;
		}

		if (sorted_recomb_regs->at(i_reg)->start <= var_reg->start &&
			sorted_recomb_regs->at(i_reg + 1)->start >= var_reg->start)
		{
			int l_recomb_reg = sorted_recomb_regs->at(i_reg + 1)->start - sorted_recomb_regs->at(i_reg)->start;
			double delta_cM = sorted_recomb_regs->at(i_reg + 1)->dbl_score - sorted_recomb_regs->at(i_reg)->dbl_score;

			double cM_slope = delta_cM / l_recomb_reg;

			double linear_interpolated_cM = sorted_recomb_regs->at(i_reg)->dbl_score + (var_reg->start - sorted_recomb_regs->at(i_reg)->start) * cM_slope;

			//// Get the distances and scale.
			//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
			//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

			bool cur_reg_is_not_side = (sorted_recomb_regs->at(i_reg)->start != var_reg->start &&
										sorted_recomb_regs->at(i_reg + 1)->start != var_reg->start);

			if (cur_reg_is_not_side &&
				(linear_interpolated_cM < MIN(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score) ||
					linear_interpolated_cM > MAX(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score)))
			{
				fprintf(stderr, "Sanity check failed for interpolation.\n");
				exit(1);
			}

			interpolated_cM = linear_interpolated_cM;
			break;
		} // overlap check.

		i_reg++;
	} // i_reg loop.

	// Following is brute-force search.
	bool CHECK_RECOMB_VAL = false;
	if (CHECK_RECOMB_VAL)
	{
		for (int i_reg = 0; i_reg < (int)sorted_recomb_regs->size() - 1; i_reg++)
		{
			if (sorted_recomb_regs->at(i_reg)->start <= var_reg->start &&
				sorted_recomb_regs->at(i_reg + 1)->start >= var_reg->start)
			{
				int l_recomb_reg = sorted_recomb_regs->at(i_reg + 1)->start - sorted_recomb_regs->at(i_reg)->start;
				double delta_cM = sorted_recomb_regs->at(i_reg + 1)->dbl_score - sorted_recomb_regs->at(i_reg)->dbl_score;

				double cM_slope = delta_cM / l_recomb_reg;

				double linear_interpolated_cM = sorted_recomb_regs->at(i_reg)->dbl_score + (var_reg->start - sorted_recomb_regs->at(i_reg)->start) * cM_slope;

				//// Get the distances and scale.
				//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
				//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

				if (linear_interpolated_cM < MIN(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score) ||
					linear_interpolated_cM > MAX(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score))
				{
					fprintf(stderr, "Sanity check failed for interpolation.\n");
					exit(1);
				}

				if (interpolated_cM != linear_interpolated_cM)
				{
					fprintf(stderr, "Recomb rate is assigned incorrectly for %s:%d: %.5f, %.5f\n",
						var_reg->chrom, var_reg->start,
						interpolated_cM, linear_interpolated_cM);
					exit(1);
				}
				break;
			} // overlap check
		} // i_reg loop.
	} // check for assignment.

	return(interpolated_cM);
}

double get_cumulative_recomb_rate_per_variant(t_annot_region* var_reg, vector<t_annot_region*>* recomb_regs)
{
	for (int i_reg = 0; i_reg < (int)recomb_regs->size()-1; i_reg++)
	{
		if (recomb_regs->at(i_reg)->start <= var_reg->start &&
			recomb_regs->at(i_reg+1)->start >= var_reg->start)
		{
			int l_recomb_reg = recomb_regs->at(i_reg + 1)->start - recomb_regs->at(i_reg)->start;
			double delta_cM = recomb_regs->at(i_reg+1)->dbl_score - recomb_regs->at(i_reg)->dbl_score;

			double cM_slope = delta_cM / l_recomb_reg;

			double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score + (var_reg->start - recomb_regs->at(i_reg)->start) * cM_slope;

			//// Get the distances and scale.
			//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
			//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

			bool cur_reg_is_not_side = (recomb_regs->at(i_reg)->start != var_reg->start &&
										recomb_regs->at(i_reg + 1)->start != var_reg->start);

			if (cur_reg_is_not_side &&
				(linear_interpolated_cM < MIN(recomb_regs->at(i_reg + 1)->dbl_score, recomb_regs->at(i_reg)->dbl_score) ||
				linear_interpolated_cM > MAX(recomb_regs->at(i_reg + 1)->dbl_score, recomb_regs->at(i_reg)->dbl_score)))
			{
				fprintf(stderr, "Sanity check failed for interpolation: %lf, %lf, %lf\n",
					recomb_regs->at(i_reg + 1)->dbl_score, linear_interpolated_cM, recomb_regs->at(i_reg)->dbl_score);
				exit(1);
			}

			return(linear_interpolated_cM);
		}
	}

	return(recomb_regs->back()->dbl_score);
}

vector<t_annot_region*>* load_binarized_binarized_per_allele_recomb_counts(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids)
{
	FILE* f_bin_geno_sig_regs = open_f(bin_geno_sig_bed_fp, "rb");

	// Load the chromosomes.
	int n_chrs = 0;
	//fread(&n_chrs, sizeof(int), 1, f_bin_geno_sig_regs);
	n_chrs = read_bin_int(f_bin_geno_sig_regs);
	vector<char*>* chr_ids = new vector<char*>();
	for (int i_chr = 0; i_chr < n_chrs; i_chr++)
	{
		char cur_chr[1000];
		//fread(cur_chr, sizeof(char), 1000, f_bin_geno_sig_regs);
		read_bin_char_array(f_bin_geno_sig_regs, cur_chr, 1000);
		chr_ids->push_back(t_string::copy_me_str(cur_chr));
	} // i_chr loop.

	fprintf(stderr, "Loaded %d chromosomes.\n", (int)chr_ids->size());
	int sample_size = 0;
	//fread(&sample_size, sizeof(int), 1, f_bin_geno_sig_regs);
	sample_size = read_bin_int(f_bin_geno_sig_regs);
	fprintf(stderr, "Reading sample size of %d.\n", sample_size);

	if (sample_size != (int)geno_sample_ids->size())
	{
		fprintf(stderr, "Sanity check failed: Sample sizes do not match: %d, %d\n", sample_size, (int)geno_sample_ids->size());
		exit(1);
	}

	int n_regs = 0;
	//fread(&n_regs, sizeof(int), 1, f_bin_geno_sig_regs);
	n_regs = read_bin_int(f_bin_geno_sig_regs);
	fprintf(stderr, "Reading %d regions.\n", n_regs);

	vector<t_annot_region*>* geno_sig_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < n_regs; i_reg++)
	{
		int i_chr = 0;
		int reg_BED_start = 0;
		int reg_BED_end = 0;
		//fread(&i_chr, sizeof(int), 1, f_bin_geno_sig_regs);
		i_chr = read_bin_int(f_bin_geno_sig_regs);
		//fread(&(reg_BED_start), sizeof(int), 1, f_bin_geno_sig_regs);
		reg_BED_start = read_bin_int(f_bin_geno_sig_regs);
		//fread(&(reg_BED_end), sizeof(int), 1, f_bin_geno_sig_regs);
		reg_BED_end = read_bin_int(f_bin_geno_sig_regs);

		// Read the region's name.
		int l_reg_name_str = 0;
		//fread(&l_reg_name_str, sizeof(int), 1, f_bin_geno_sig_regs);
		l_reg_name_str = read_bin_int(f_bin_geno_sig_regs);
		char* cur_reg_name = new char[l_reg_name_str + 2];
		memset(cur_reg_name, 0, sizeof(char) * (l_reg_name_str + 2));
		//fread(cur_reg_name, sizeof(char), l_reg_name_str, f_bin_geno_sig_regs);
		read_bin_char_array(f_bin_geno_sig_regs, cur_reg_name, l_reg_name_str);

		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
		reg->start = translate_coord(reg_BED_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		reg->end = translate_coord(reg_BED_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		reg->strand = '+';
		reg->name = cur_reg_name;

		void** cur_reg_info = new void*[2];
		int** cur_reg_geno_sig = new int*[2];
		cur_reg_geno_sig[0] = new int[sample_size + 2];
		cur_reg_geno_sig[1] = new int[sample_size + 2];
		if (fread(cur_reg_geno_sig[0], sizeof(int), sample_size, f_bin_geno_sig_regs) != (size_t)sample_size)
		{
			fprintf(stderr, "Could not read counts @ %d. reg\n", i_reg);
			exit(1);
		}

		if (fread(cur_reg_geno_sig[1], sizeof(int), sample_size, f_bin_geno_sig_regs) != (size_t)sample_size)
		{
			fprintf(stderr, "Could not read counts @ %d. reg\n", i_reg);
			exit(1);
		}

		cur_reg_info[0] = cur_reg_geno_sig;
		cur_reg_info[1] = NULL;

		reg->data = cur_reg_info;

		geno_sig_regs->push_back(reg);
	} // i_reg loop.

	// Close the file.
	close_f(f_bin_geno_sig_regs, bin_geno_sig_bed_fp);

	return(geno_sig_regs);
}

t_annot_region* copy_geno_signal_region(t_annot_region* src_geno_signal_reg, int sample_size)
{
	t_annot_region* copy_signal_reg = duplicate_region(src_geno_signal_reg);

	void** src_info = (void**)(src_geno_signal_reg->data);
	char* src_geno_sig = (char*)(src_info[0]);

	void** copy_info = new void* [2];
	char* copy_geno_sig = new char[sample_size];
	memcpy(copy_geno_sig, src_geno_sig, sizeof(char) * sample_size);

	copy_info[0] = copy_geno_sig;
	copy_signal_reg->data = copy_info;

	return(copy_signal_reg);
} // copy_geno_signal_region function.

void delete_geno_signal_regions(vector<t_annot_region*>* geno_sig_regs)
{
	for (int i_reg = 0; i_reg < vecsize(geno_sig_regs); i_reg++)
	{
		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		delete[] cur_reg_geno_sig;
		delete[] cur_reg_info;
	} // i_reg loop.

	delete_annot_regions(geno_sig_regs);
} // delete_geno_signal_regions function.

void resample_phased_haplotypes_per_recombination_rates(char* haplocoded_genotype_matrix_fp, 
														char* sample_ids_list_fp, 
														char* recombination_rate_dir, 
														double upsampling_rate,
														double N_e,
														double allele_error_prob,
														char* op_fp)
{
	//double N_e = 10 ^ 6;
	//double allele_error_prob = pow(10, -4);

	fprintf(stderr, "Re-sampling recombination-aware genotypes using N_e=%d and allelic error=%.6f\n", (int)N_e, allele_error_prob);
	vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_genotype_matrix_fp, sample_ids_list_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n", vecsize(haplocoded_geno_regs), vecsize(sample_ids));
	double n_original_haps = 2.0 * vecsize(sample_ids);

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	int n_resampled_sample_size = (int)(upsampling_rate * vecsize(sample_ids));
	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	// Load the recombination rates.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for(int i_chr = 0; i_chr < vecsize(restr_geno_regs->chr_ids); i_chr++)
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
		for (int i_reg = 0; i_reg < vecsize(cur_chrom_var_regs); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < vecsize(cur_chrom_var_regs); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
			char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
			memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
			cur_reg_info[1] = cur_reg_resampled_geno;
		} // i_reg loop.

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			// Sample the two haplotypes.
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				fprintf(stderr, "Re-sampling sample %d (%d)\n", i_s, i_hap);				

				int cur_haplo_state = -1;
				int n_recombs = 0;
				for (int i_reg = 0; i_reg < vecsize(cur_chrom_var_regs); i_reg++)
				{
					if (i_reg == 0)
					{
						cur_haplo_state = MIN((n_original_haps - 1), floor(rng->random_double_ran3() * n_original_haps));
					}
					else
					{
						// Get the recombination rate between these positions.
						double r_m = fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_chrom_var_regs->at(i_reg -1 )->dbl_score);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

						double other_prob = tau_m / n_original_haps;
						double self_prob = (1 - tau_m) + (tau_m / n_original_haps);

						double rand_cumul_val = rng->random_double_ran3();

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Sample %d: Var %d (%s:%d); Haplo_state: %d; dcM:%.4f (%.5f - %.5f), other_prob=%.3f, self_prob=%.3f\n",
								i_s, i_reg,
								cur_chrom_var_regs->at(i_reg)->chrom, cur_chrom_var_regs->at(i_reg)->start,
								cur_haplo_state, r_m,
								cur_chrom_var_regs->at(i_reg)->dbl_score, cur_chrom_var_regs->at(i_reg - 1)->dbl_score,
								other_prob, self_prob);
						}

						double cur_cumul_val = 0;
						for (int hap_state_i = 0; hap_state_i < n_original_haps; hap_state_i++)
						{
							// Update the cumulative.
							if (hap_state_i == cur_haplo_state)
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
								if (cur_haplo_state != hap_state_i)
								{
									n_recombs++;
								}

								cur_haplo_state = hap_state_i;
								break;
							}
						} // hap_state_i loop.
					} // region check for initing the haplo state.

					// Copy the allele.
					int sample_i = (cur_haplo_state - (cur_haplo_state % 2)) / 2;
					int sampled_hap_i = (cur_haplo_state % 2);
					void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
					char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
					char* cur_reg_sampled_geno = (char*)(cur_reg_info[1]);
					int cur_allele = get_allele_per_haplotype(cur_reg_geno_sig[sample_i], sampled_hap_i);

					// Copy the allele to the haplotype 
					char cur_val = (char)(cur_reg_sampled_geno[i_s]);
					cur_reg_sampled_geno[i_s] = cur_val | (cur_allele << i_hap);
				} // i_reg loop.

				fprintf(stderr, "Re-sampling sample %d (%d): %d recombinations.\n", i_s, i_hap, n_recombs);
			} // i_hap loop.
 		} // i_s loop.
	} // i_chr loop.

	fprintf(stderr, "Saving re-sampled genotypes.\n");

	// Save the genotypes.
	FILE* f_op = open_f(op_fp, "w");
	for (int i_reg = 0; i_reg < vecsize(haplocoded_geno_regs); i_reg++)
	{
		void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		char* cur_resampled_geno_sig = (char*)(cur_reg_info[1]);

		fprintf(f_op, "%s\t%d\t%d\t%s", haplocoded_geno_regs->at(i_reg)->chrom, 
				translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				haplocoded_geno_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			fprintf(f_op, "\t%d", (int)(cur_resampled_geno_sig[i_s]));
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);

	FILE* f_samples = open_f("sample_ids.list", "w");
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		fprintf(f_samples, "sample_%d\n", i_s);
	} // i_s loop.
	fclose(f_samples);
}



char** get_match_by_identifiers()
{
	char** match_by_identifiers = new char*[N_MATCH_BY_IDS + 2];
	match_by_identifiers[MATCH_BY_NAME] = t_string::copy_me_str("byName");
	match_by_identifiers[MATCH_BY_START_POSN] = t_string::copy_me_str("byStartPosn");
	return(match_by_identifiers);
}

static void* thread_callback_extract_BEAGLE_ref_option(void* __thread_info_ptr)
{
	void** cur_thread_info = (void**)(__thread_info_ptr);

	int* int_vals = (int*)(cur_thread_info[0]);
	int thr_i = int_vals[0];
	//int n_thrs = int_vals[1];
	int var_start_i = int_vals[2];
	int var_end_i = int_vals[3];
	
	//double* dbl_vals = (double*)(cur_thread_info[1]);

	vector<t_annot_region*>* cur_chr_ref_panel_var_regs = (vector<t_annot_region*>*)(cur_thread_info[2]);
	vector<char*>* reference_haplo_sample_ids = (vector<char*>*)(cur_thread_info[3]);

	//fprintf(stderr, "Thread-%d: %d-%d // %d\n", thr_i, var_start_i, var_end_i, vecsize(cur_chr_ref_panel_var_regs));
	t_string::print_padded_string(stderr, '\r', 100, "Ref Thread-%d: %d-%d // %d", thr_i, var_start_i, var_end_i, vecsize(cur_chr_ref_panel_var_regs));

	char cur_thread_file_name[1000];
	sprintf(cur_thread_file_name, "input_option_%d.vcf.gz", thr_i);
	//FILE* f_ref_option = open_f(cur_thread_file_name, "w");
	gzFile f_ref_option = gzopen(cur_thread_file_name, "wb");
	if (!f_ref_option)
	{
		fprintf(stderr, "%s(%d): Could not open %s for writing..\n", __FILE__, __LINE__, cur_thread_file_name);
		exit(1);
	}

	size_t L_LINE_BUFFER = 100 * 1000;
	char* line_buffer = new char[L_LINE_BUFFER];
	t_string* MAIN_LINE_STR_BUFFER = new t_string();

	for (int i_reg = var_start_i; i_reg < MIN(var_end_i, vecsize(cur_chr_ref_panel_var_regs)); i_reg++)
	{
		// Empty the main line buffer string.
		MAIN_LINE_STR_BUFFER->empty();

		// Write the legend entry.			
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
		if (toks->size() < 3)
		{
			fprintf(stderr, "Could not parse region name: %s\n", cur_chr_ref_panel_var_regs->at(i_reg)->name);
			exit(1);
		}

		char* cur_var_name = toks->at(0)->str();
		char* cur_var_ref_str = toks->at(1)->str();
		char* cur_var_alt_str = toks->at(2)->str();

		// Save the reference haplotype file (-h).
		void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
		char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

		// Write the alleles for each sample.
		// 22      20000086        rs138720731     T       C       100     PASS    . GT
		//fprintf(f_ref_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
		//	cur_chr_ref_panel_var_regs->at(i_reg)->chrom,
		//	cur_chr_ref_panel_var_regs->at(i_reg)->start,
		//	cur_var_name,
		//	cur_var_ref_str, cur_var_alt_str);
		sprintf(line_buffer, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
			cur_chr_ref_panel_var_regs->at(i_reg)->chrom,
			cur_chr_ref_panel_var_regs->at(i_reg)->start,
			cur_var_name,
			cur_var_ref_str, cur_var_alt_str);

		MAIN_LINE_STR_BUFFER->concat_string(line_buffer);

		// Write the alleles for each sample.
		for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
		{
			int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
			int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

			// Setup the phased genotype string.
			// Following seems to match getallele(0) | getallele(1)
			char geno_str[10];
			sprintf(geno_str, "\t%d|%d", cur_hap0, cur_hap1);

			MAIN_LINE_STR_BUFFER->concat_string(geno_str);
		} // i_s loop.

		//fprintf(f_ref_option, "\n");
		MAIN_LINE_STR_BUFFER->concat_char('\n');

		if (gzwrite(f_ref_option, MAIN_LINE_STR_BUFFER->str(), MAIN_LINE_STR_BUFFER->length()) != MAIN_LINE_STR_BUFFER->length())
		{
			fprintf(stderr, "%s(%d): Could not write line buffer in thread-%d..\n", __FILE__, __LINE__, thr_i);
			exit(1);
		}

		t_string::clean_tokens(toks);
	} // i_reg loop.

	gzclose(f_ref_option);

	return(NULL);
}

static void* thread_callback_extract_BEAGLE_gt_option(void* __thread_info_ptr)
{
	void** cur_thread_info = (void**)(__thread_info_ptr);

	int* int_vals = (int*)(cur_thread_info[0]);
	int thr_i = int_vals[0];
	//int n_thrs = int_vals[1];
	int var_start_i = int_vals[2];
	int var_end_i = int_vals[3];
	int save_phased_gt_option = int_vals[4];
	int max_input_geno = int_vals[5];

	vector<t_annot_region*>* cur_chr_input_geno_sig_regs = (vector<t_annot_region*>*)(cur_thread_info[2]);
	vector<char*>* input_geno_sample_ids = (vector<char*>*)(cur_thread_info[3]);

	//fprintf(stderr, "GT Thread-%d: %d-%d // %d [max_geno=%d]\n", thr_i, var_start_i, var_end_i, vecsize(cur_chr_input_geno_sig_regs), max_input_geno);
	t_string::print_padded_string(stderr, '\r', 100, "GT Thread-%d: %d-%d // %d [max_geno=%d]", thr_i, var_start_i, var_end_i, vecsize(cur_chr_input_geno_sig_regs), max_input_geno);

	char cur_thread_file_name[1000];
	sprintf(cur_thread_file_name, "gt_option_%d.vcf.gz", thr_i);
	//FILE* f_gt_option = open_f(cur_thread_file_name, "w");
	gzFile f_gt_option = gzopen(cur_thread_file_name, "wb");
	if (!f_gt_option)
	{
		fprintf(stderr, "%s(%d): Could not open %s for writing..\n", __FILE__, __LINE__, cur_thread_file_name);
		exit(1);
	}

	size_t L_LINE_BUFFER = 100 * 1000;
	char* line_buffer = new char[L_LINE_BUFFER];
	t_string* MAIN_LINE_STR_BUFFER = new t_string();

	//vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
	for (int i_reg = var_start_i; i_reg < MIN(var_end_i, (int)cur_chr_input_geno_sig_regs->size()); i_reg++)
	{
		// Empty the main line buffer string.
		MAIN_LINE_STR_BUFFER->empty();

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
		char* cur_var_name = toks->at(0)->str();
		char* cur_var_ref_str = toks->at(1)->str();
		char* cur_var_alt_str = toks->at(2)->str();

		// Apparently, BEAGLE does not want this.
		// Write the alleles for each sample.
		// 22      20000086        rs138720731     T       C       100     PASS    . GT
		//fprintf(f_gt_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
		//	cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
		//	cur_chr_input_geno_sig_regs->at(i_reg)->start,
		//	cur_var_name,
		//	cur_var_ref_str, cur_var_alt_str);
		sprintf(line_buffer, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
			cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
			cur_chr_input_geno_sig_regs->at(i_reg)->start,
			cur_var_name,
			cur_var_ref_str, cur_var_alt_str);

		MAIN_LINE_STR_BUFFER->concat_string(line_buffer);

		// Save the reference haplotype file (-h).
		void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
		char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

		for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
		{
			if (max_input_geno == 3)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				// Following seems to match to ref option when haplocoded option is chosen: getallele(0) | getallele(1)
				if (save_phased_gt_option)
				{
					//fprintf(f_gt_option, "\t%d|%d", cur_hap0, cur_hap1);
					char geno_str[10];
					sprintf(geno_str, "\t%d|%d", cur_hap0, cur_hap1);
					MAIN_LINE_STR_BUFFER->concat_string(geno_str);
				}
				else
				{
					char geno_str[10];

					int geno = cur_hap0 + cur_hap1;

					if (geno == 0)
					{
						//fprintf(f_gt_option, "\t0/0");
						sprintf(geno_str, "\t0/0");
					}
					else if (geno == 1)
					{
						//fprintf(f_gt_option, "\t0/1");
						sprintf(geno_str, "\t0/1");
					}
					else if (geno == 2)
					{
						//fprintf(f_gt_option, "\t1/1");
						sprintf(geno_str, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(1);
					}

					MAIN_LINE_STR_BUFFER->concat_string(geno_str);
				} // phased gt option check.
			} // max_geno=3 check.
			else
			{
				if (save_phased_gt_option)
				{
					fprintf(stderr, "%s(%d): We are not supposed to be here; saving GT as phased for an unphased panel..\n", __FILE__, __LINE__);
					exit(1);
				}

				int geno = per_sample_haplocoded_geno[i_s];

				char geno_str[10];

				if (geno == 0)
				{
					//fprintf(f_gt_option, "\t0/0");
					sprintf(geno_str, "\t0/0");
				}
				else if (geno == 1)
				{
					//fprintf(f_gt_option, "\t0/1");
					sprintf(geno_str, "\t0/1");
				}
				else if (geno == 2)
				{
					//fprintf(f_gt_option, "\t1/1");
					sprintf(geno_str, "\t1/1");
				}
				else
				{
					fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d)\n",
						cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
						(int)(per_sample_haplocoded_geno[i_s]), geno);

					exit(1);
				}

				MAIN_LINE_STR_BUFFER->concat_string(geno_str);
			} // max_geno!=3 check.

		} // i_s loop.

		//fprintf(f_gt_option, "\n");
		MAIN_LINE_STR_BUFFER->concat_char('\n');

		if (gzwrite(f_gt_option, MAIN_LINE_STR_BUFFER->str(), MAIN_LINE_STR_BUFFER->length()) != MAIN_LINE_STR_BUFFER->length())
		{
			fprintf(stderr, "%s(%d): Could not write line buffer in thread-%d..\n", __FILE__, __LINE__, thr_i);
			exit(1);
		}

		t_string::clean_tokens(toks);
	} // i_reg loop.

	//close_f(f_gt_option, NULL);
	gzclose(f_gt_option);

	return(NULL);
}

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_w_header_multithreaded(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	char* header_fp,
	bool save_phased_gt_option,
	int n_threads)
{
	fprintf(stderr, "%d-thread Saving the BEAGLE files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-ref %s\n\
	-gt %s\n\
	header_file %s",
		n_threads,
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		ref_option_fp, gt_option_fp,
		header_fp);

	if (save_phased_gt_option)
	{
		fprintf(stderr, "**Saving phased gt option.**\n");
	}

	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	vector<char*>* header_lines = buffer_file(header_fp);
	fprintf(stderr, "Loaded %d header lines.\n", vecsize(header_lines));

	int max_ref_geno = get_max_genotype_value(reference_haplo_geno_sig_regs, reference_haplo_sample_ids);
	if (max_ref_geno != 3)
	{
		fprintf(stderr, "%s(%d): Reference panel does not seem to be phased, cannot save unphased reference option for BEAGLE.\n", __FILE__, __LINE__);
		exit(1);
	}

	// Saving the legend file.
	t_string* ref_header_line_str = new t_string();
	for (int i_h = 0; i_h < vecsize(header_lines); i_h++)
	{
		ref_header_line_str->concat_string(header_lines->at(i_h));
		ref_header_line_str->concat_string("\n");
	} // i_h loop.
	ref_header_line_str->concat_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
	{
		ref_header_line_str->concat_string("\t");
		ref_header_line_str->concat_string(reference_haplo_sample_ids->at(i_s));
	} // i_s loop.
	//fprintf(f_ref_header, "\n");
	ref_header_line_str->concat_string("\n");

	gzFile f_ref_header = gzopen("REF_HEADER.vcf.gz", "wb");
	gzwrite(f_ref_header, ref_header_line_str->str(), ref_header_line_str->length());
	delete(ref_header_line_str);
	gzclose(f_ref_header);

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Process all chromosomes.
	for (int i_chr = 0; i_chr < (int)restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];

		int n_vars_per_thread = MAX(1, ceil(vecsize(cur_chr_ref_panel_var_regs) / ((double)n_threads)));

		int cur_var_start_i = 0;
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int i_thr = 0; i_thr < n_threads; i_thr++)
		{
			if (cur_var_start_i > vecsize(cur_chr_ref_panel_var_regs))
			{
				fprintf(stderr, "Early stopping @ thread_i=%d\n", i_thr);
				break;
			}

			void** cur_thread_info = new void* [10];
			int* int_vals = new int[10];
			cur_thread_info[0] = int_vals;
			int_vals[0] = i_thr;
			int_vals[1] = n_threads;
			int_vals[2] = cur_var_start_i;
			int_vals[3] = MIN(cur_var_start_i + n_vars_per_thread, vecsize(cur_chr_ref_panel_var_regs));
			double* dbl_vals = new double[10];
			cur_thread_info[1] = dbl_vals;
			cur_thread_info[2] = cur_chr_ref_panel_var_regs;
			cur_thread_info[3] = reference_haplo_sample_ids;

			t_ansi_thread* thread = new t_ansi_thread(thread_callback_extract_BEAGLE_ref_option, cur_thread_info);
			threads->push_back(thread);
			thread->run_thread();

			cur_var_start_i += n_vars_per_thread;

		} // i_thr loop.

		vector<char*>* per_thread_ref_files = new vector<char*>();
		per_thread_ref_files->push_back(t_string::copy_me_str("REF_HEADER.vcf.gz"));
		for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
		{
			threads->at(i_thr)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Ref-Thread-%d finished..", i_thr);

			char cur_thread_file_name[1000];
			sprintf(cur_thread_file_name, "input_option_%d.vcf.gz", i_thr);

			per_thread_ref_files->push_back(t_string::copy_me_str(cur_thread_file_name));
		} // i_thr loop.

		t_string::print_padded_string(stderr, '\n', 100, "Concatenating %d Ref matrices..", vecsize(per_thread_ref_files));
		concatenateGzipFiles(ref_option_fp, per_thread_ref_files);

		for (int i_thr = 0; i_thr < vecsize(per_thread_ref_files); i_thr++)
		{
			// Delete this file to make sure it does not interfere??
			t_string::print_padded_string(stderr, '\n', 100, "Deleting %s", per_thread_ref_files->at(i_thr));
			delete_file(per_thread_ref_files->at(i_thr));
		} // i_thr loop.

		// Do not process multiple chromosome here.
		break;
	} // i_chr loop.

	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", (int)input_geno_sample_ids->size());
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_signal_regions_wrapper(input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	// Make sure input genotypes are phased if they are going to be saved as phased.
	int max_input_geno = get_max_genotype_value(input_geno_sig_regs, input_geno_sample_ids);
	if (max_input_geno != 3 &&
		save_phased_gt_option)
	{
		fprintf(stderr, "%s(%d): Study panel does not seem to be phased, cannot save unphased gt with phased option for BEAGLE.\n", __FILE__, __LINE__);
		exit(1);
	}

	// Reset the scores to 0 to indicate these are not assigned a ref region, yet.
	for (int i_reg = 0; i_reg < (int)input_geno_sig_regs->size(); i_reg++)
	{
		input_geno_sig_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	// Write the header.
	t_string* gt_header_line_str = new t_string();
	for (int i_h = 0; i_h < vecsize(header_lines); i_h++)
	{
		gt_header_line_str->concat_string(header_lines->at(i_h));
		gt_header_line_str->concat_string("\n");
	} // i_h loop.
	gt_header_line_str->concat_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
	{
		gt_header_line_str->concat_string("\t");
		gt_header_line_str->concat_string(input_geno_sample_ids->at(i_s));
	} // i_s loop.
	//fprintf(f_ref_header, "\n");
	gt_header_line_str->concat_string("\n");

	gzFile f_gt_header = gzopen("GT_HEADER.vcf.gz", "wb");
	gzwrite(f_gt_header, gt_header_line_str->str(), gt_header_line_str->length());
	delete(gt_header_line_str);
	gzclose(f_gt_header);

	// Restructure the genotype signal regions.
	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		int n_vars_per_thread = MAX(1, ceil(vecsize(cur_chr_input_geno_sig_regs) / ((double)n_threads)));

		int cur_var_start_i = 0;
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int i_thr = 0; i_thr < n_threads; i_thr++)
		{
			if (cur_var_start_i > vecsize(cur_chr_input_geno_sig_regs))
			{
				fprintf(stderr, "Early stopping @ thread_i=%d\n", i_thr);
				break;
			}

			void** cur_thread_info = new void* [10];
			int* int_vals = new int[10];
			cur_thread_info[0] = int_vals;
			int_vals[0] = i_thr;
			int_vals[1] = n_threads;
			int_vals[2] = cur_var_start_i;
			int_vals[3] = MIN(cur_var_start_i + n_vars_per_thread, vecsize(cur_chr_input_geno_sig_regs));
			int_vals[4] = save_phased_gt_option;
			int_vals[5] = max_input_geno;
			double* dbl_vals = new double[10];
			cur_thread_info[1] = dbl_vals;
			cur_thread_info[2] = cur_chr_input_geno_sig_regs;
			cur_thread_info[3] = input_geno_sample_ids;

			t_ansi_thread* thread = new t_ansi_thread(thread_callback_extract_BEAGLE_gt_option, cur_thread_info);
			threads->push_back(thread);
			thread->run_thread();

			cur_var_start_i += n_vars_per_thread;
		} // i_thr loop.

		vector<char*>* per_thread_gt_files = new vector<char*>();
		per_thread_gt_files->push_back(t_string::copy_me_str("GT_HEADER.vcf.gz"));
		for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
		{
			threads->at(i_thr)->wait_thread();

			t_string::print_padded_string(stderr, '\r', 100, "Ref-Thread-%d finished..", i_thr);

			char cur_thread_file_name[1000];
			sprintf(cur_thread_file_name, "gt_option_%d.vcf.gz", i_thr);

			per_thread_gt_files->push_back(t_string::copy_me_str(cur_thread_file_name));
		} // i_thr loop.

		fprintf(stderr, "Concatenating %d GT files..                       \n", vecsize(per_thread_gt_files));
		concatenateGzipFiles(gt_option_fp, per_thread_gt_files);

		for (int i_thr = 0; i_thr < vecsize(per_thread_gt_files); i_thr++)
		{
			// Delete this file to make sure it does not interfere??
			delete_file(per_thread_gt_files->at(i_thr));
		} // i_thr loop.

		break;
	} // i_chr loop. 
}

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_multithreaded(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option,
	int n_threads)
{
	fprintf(stderr, "%d-thread Saving the BEAGLE files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-ref %s\n\
	-gt %s\n",
		n_threads,
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		ref_option_fp, gt_option_fp);

	if (save_phased_gt_option)
	{
		fprintf(stderr, "**Saving phased gt option.**\n");
	}

	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	int max_ref_geno = get_max_genotype_value(reference_haplo_geno_sig_regs, reference_haplo_sample_ids);
	if (max_ref_geno != 3)
	{
		fprintf(stderr, "%s(%d): Reference panel does not seem to be phased, cannot save unphased reference option for BEAGLE.\n", __FILE__, __LINE__);
		exit(1);
	}

	//// Saving the legend file.
	//FILE* f_ref_header = open_f("REF_HEADER.vcf", "w");
	//fprintf(f_ref_header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	//for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
	//{
	//	fprintf(f_ref_header, "\t%s", reference_haplo_sample_ids->at(i_s));
	//} // i_s loop.
	//fprintf(f_ref_header, "\n");
	//close_f(f_ref_header, NULL);
	t_string* ref_header_line_str = new t_string();
	ref_header_line_str->concat_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
	{
		ref_header_line_str->concat_string("\t");
		ref_header_line_str->concat_string(reference_haplo_sample_ids->at(i_s));
	} // i_s loop.
	//fprintf(f_ref_header, "\n");
	ref_header_line_str->concat_string("\n");

	gzFile f_ref_header = gzopen("REF_HEADER.vcf.gz", "wb");
	gzwrite(f_ref_header, ref_header_line_str->str(), ref_header_line_str->length());
	delete(ref_header_line_str);
	gzclose(f_ref_header);

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Process all chromosomes.
	for (int i_chr = 0; i_chr < (int)restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];

		int n_vars_per_thread = MAX(1, ceil(vecsize(cur_chr_ref_panel_var_regs) / ((double)n_threads)));

		int cur_var_start_i = 0;
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int i_thr = 0; i_thr < n_threads; i_thr++)
		{
			if (cur_var_start_i > vecsize(cur_chr_ref_panel_var_regs))
			{
				fprintf(stderr, "Early stopping @ thread_i=%d\n", i_thr);
				break;
			}

			void** cur_thread_info = new void* [10];
			int* int_vals = new int[10];
			cur_thread_info[0] = int_vals;
			int_vals[0] = i_thr;
			int_vals[1] = n_threads;
			int_vals[2] = cur_var_start_i;
			int_vals[3] = MIN(cur_var_start_i + n_vars_per_thread, vecsize(cur_chr_ref_panel_var_regs));
			double* dbl_vals = new double[10];
			cur_thread_info[1] = dbl_vals;
			cur_thread_info[2] = cur_chr_ref_panel_var_regs;
			cur_thread_info[3] = reference_haplo_sample_ids;

			t_ansi_thread* thread = new t_ansi_thread(thread_callback_extract_BEAGLE_ref_option, cur_thread_info);
			threads->push_back(thread);
			thread->run_thread();

			cur_var_start_i += n_vars_per_thread;

		} // i_thr loop.
		
		vector<char*>* per_thread_ref_files = new vector<char*>();
		per_thread_ref_files->push_back(t_string::copy_me_str("REF_HEADER.vcf.gz"));
		for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
		{
			threads->at(i_thr)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Ref-Thread-%d finished..", i_thr);

			char cur_thread_file_name[1000];
			sprintf(cur_thread_file_name, "input_option_%d.vcf.gz", i_thr);

			per_thread_ref_files->push_back(t_string::copy_me_str(cur_thread_file_name));
		} // i_thr loop.

		fprintf(stderr, "Concatenating %d Ref matrices..                       \n", vecsize(per_thread_ref_files));
		concatenateGzipFiles(ref_option_fp, per_thread_ref_files);

		for (int i_thr = 0; i_thr < vecsize(per_thread_ref_files); i_thr++)
		{
			// Delete this file to make sure it does not interfere??
			t_string::print_padded_string(stderr, '\r', 80, "Deleting %s", per_thread_ref_files->at(i_thr));

			delete_file(per_thread_ref_files->at(i_thr));
		} // i_thr loop.

		// Do not process multiple chromosome here.
		break;
	} // i_chr loop.
	
	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", (int)input_geno_sample_ids->size());
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_signal_regions_wrapper(input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	// Make sure input genotypes are phased if they are going to be saved as phased.
	int max_input_geno = get_max_genotype_value(input_geno_sig_regs, input_geno_sample_ids);
	if (max_input_geno != 3 &&
		save_phased_gt_option)
	{
		fprintf(stderr, "%s(%d): Study panel does not seem to be phased, cannot save unphased gt with phased option for BEAGLE.\n", __FILE__, __LINE__);
		exit(1);
	}

	// Reset the scores to 0 to indicate these are not assigned a ref region, yet.
	for (int i_reg = 0; i_reg < (int)input_geno_sig_regs->size(); i_reg++)
	{
		input_geno_sig_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	//// Write the header.
	//FILE* f_gt_header = open_f("GT_HEADER.vcf", "w");
	//fprintf(f_gt_header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	//for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
	//{
	//	fprintf(f_gt_header, "\t%s", input_geno_sample_ids->at(i_s));
	//} // i_s loop.
	//fprintf(f_gt_header, "\n");
	//close_f(f_gt_header, NULL);
	t_string* gt_header_line_str = new t_string();
	gt_header_line_str->concat_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
	{
		gt_header_line_str->concat_string("\t");
		gt_header_line_str->concat_string(input_geno_sample_ids->at(i_s));
	} // i_s loop.
	//fprintf(f_ref_header, "\n");
	gt_header_line_str->concat_string("\n");

	gzFile f_gt_header = gzopen("GT_HEADER.vcf.gz", "wb");
	gzwrite(f_gt_header, gt_header_line_str->str(), gt_header_line_str->length());
	delete(gt_header_line_str);
	gzclose(f_gt_header);

	// Restructure the genotype signal regions.
	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		int n_vars_per_thread = MAX(1, ceil(vecsize(cur_chr_input_geno_sig_regs) / ((double)n_threads)));

		int cur_var_start_i = 0;
		vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
		for (int i_thr = 0; i_thr < n_threads; i_thr++)
		{
			if (cur_var_start_i > vecsize(cur_chr_input_geno_sig_regs))
			{
				fprintf(stderr, "Early stopping @ thread_i=%d\n", i_thr);
				break;
			}

			void** cur_thread_info = new void* [10];
			int* int_vals = new int[10];
			cur_thread_info[0] = int_vals;
			int_vals[0] = i_thr;
			int_vals[1] = n_threads;
			int_vals[2] = cur_var_start_i;
			int_vals[3] = MIN(cur_var_start_i + n_vars_per_thread, vecsize(cur_chr_input_geno_sig_regs));
			int_vals[4] = save_phased_gt_option;
			int_vals[5] = max_input_geno;
			double* dbl_vals = new double[10];
			cur_thread_info[1] = dbl_vals;
			cur_thread_info[2] = cur_chr_input_geno_sig_regs;
			cur_thread_info[3] = input_geno_sample_ids;

			t_ansi_thread* thread = new t_ansi_thread(thread_callback_extract_BEAGLE_gt_option, cur_thread_info);
			threads->push_back(thread);
			thread->run_thread();

			cur_var_start_i += n_vars_per_thread;
		} // i_thr loop.

		vector<char*>* per_thread_gt_files = new vector<char*>();
		per_thread_gt_files->push_back(t_string::copy_me_str("GT_HEADER.vcf.gz"));
		for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
		{
			threads->at(i_thr)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Ref-Thread-%d finished..", i_thr);

			char cur_thread_file_name[1000];
			sprintf(cur_thread_file_name, "gt_option_%d.vcf.gz", i_thr);

			per_thread_gt_files->push_back(t_string::copy_me_str(cur_thread_file_name));
		} // i_thr loop.

		fprintf(stderr, "Concatenating %d GT files..\n", vecsize(per_thread_gt_files));
		concatenateGzipFiles(gt_option_fp, per_thread_gt_files);

		for (int i_thr = 0; i_thr < vecsize(per_thread_gt_files); i_thr++)
		{
			// Delete this file to make sure it does not interfere??
			delete_file(per_thread_gt_files->at(i_thr));
		} // i_thr loop.

		break;
	} // i_chr loop. 
}

// ref and gt options must have the same number of variants.
// Make sure to code the gt and ref options with the same haplotype representation, i.e., 1 vs 2.
void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option)
{
	fprintf(stderr, "Saving the BEAGLE files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-ref %s\n\
	-gt %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		ref_option_fp, gt_option_fp);

	if (save_phased_gt_option)
	{
		fprintf(stderr, "**Saving phased gt option.**\n");
	}

	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	//vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	FILE* f_ref_option = open_f(ref_option_fp, "w");
	fprintf(f_ref_option, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
	{
		fprintf(f_ref_option, "\t%s", reference_haplo_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_ref_option, "\n");

	for (int i_chr = 0; i_chr < (int)restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_ref_panel_var_regs->size(); i_reg++)
		{
			// Write the legend entry.			
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
			if (toks->size() < 3)
			{
				fprintf(stderr, "Could not parse region name: %s\n", cur_chr_ref_panel_var_regs->at(i_reg)->name);
				exit(1);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_ref_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_ref_panel_var_regs->at(i_reg)->chrom,
				cur_chr_ref_panel_var_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Write the alleles for each sample.
			for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				// Setup the phased genotype string.
				// Following seems to match getallele(0) | getallele(1)
				char geno_str[10];
				sprintf(geno_str, "%d|%d", cur_hap0, cur_hap1);

				fprintf(f_ref_option, "\t%s", geno_str);
			} // i_s loop.

			fprintf(f_ref_option, "\n");
			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop.
	close_f(f_ref_option, ref_option_fp);

	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", (int)input_geno_sample_ids->size());
	//vector<t_annot_region*>* input_geno_sig_regs = load_variant_genotype_signal_regions(input_genotype_matrix_matbed_fp, input_geno_sample_ids);
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_signal_regions_wrapper(input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	// Reset the scores to 0 to indicate these are not assigned a ref region, yet.
	for (int i_reg = 0; i_reg < (int)input_geno_sig_regs->size(); i_reg++)
	{
		input_geno_sig_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	FILE* f_gt_option = open_f(gt_option_fp, "w");

	// Write the header.
	fprintf(f_gt_option, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
	{
		fprintf(f_gt_option, "\t%s", input_geno_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_gt_option, "\n");

	// Restructure the genotype signal regions.
	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Apparently, BEAGLE does not want this.
			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_gt_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				// Following seems to match to ref option when haplocoded option is chosen: getallele(0) | getallele(1)
				if(save_phased_gt_option)
				{ 
					fprintf(f_gt_option, "\t%d|%d", cur_hap0, cur_hap1);
				}
				else
				{
					int geno = cur_hap0 + cur_hap1;

					if (geno == 0)
					{
						fprintf(f_gt_option, "\t0/0");
					}
					else if (geno == 1)
					{
						fprintf(f_gt_option, "\t0/1");
					}
					else if (geno == 2)
					{
						fprintf(f_gt_option, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(1);
					}
				}

			} // i_s loop.

			fprintf(f_gt_option, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_gt_option, gt_option_fp);
}
