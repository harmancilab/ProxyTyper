#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <ctype.h>
#include <zlib.h>
#include <algorithm>
#include "prxytypr_variation_tools.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_genome_sequence_tools.h"
#include "prxytypr_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_nucleotide.h"

#include "prxytypr_annot_region_tools.h"
//#include "prxytypr_structure.h"
#include "prxytypr_nomenclature.h"
//#include "../../../lib/genomics_utils/alignment/alignment_tools.h"
//#include "prxytypr_multi_sequence_feature_utils.h"
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_vector_macros.h"

using namespace std;

bool __DUMP_VARIATION_TOOLS_MSGS__ = false;

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

void add_noise_to_geno_signal_regions(char* geno_sig_regs_fp, char* sample_ids_list_fp, double allele_error_prob, char* op_fp)
{
	if (!check_file(geno_sig_regs_fp) ||
		!check_file(sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find genotype matbed or sample ids @ %s (%s)\n", geno_sig_regs_fp, sample_ids_list_fp);
		exit(1);
	}

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	fprintf(stderr, "Loaded %d variant genotypes for %d subjects, adding allele noise with %.4f probability.\n", vecsize(geno_sig_regs), vecsize(sample_ids), allele_error_prob);

	int max_geno = get_max_genotype_value(geno_sig_regs, sample_ids);

	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, note that this may introduce biases for genotypes.\n");
	}
	else
	{
		fprintf(stderr, "Genotypes seems to be haplocoded.\n");
	}

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	for (int i_var = 0; i_var < vecsize(geno_sig_regs); i_var++)
	{
		void** var_info = (void**)(geno_sig_regs->at(i_var)->data);
		char* geno_sig = (char*)(var_info[0]);

		for (int i_s = 0; i_s < vecsize(sample_ids); i_s++)
		{
			if (max_geno == 3)
			{
				int hap0_err = (rng->random_double_ran3() < allele_error_prob) ? (1) : (0);
				int hap1_err = (rng->random_double_ran3() < allele_error_prob) ? (1) : (0);

				int all0 = get_allele_per_haplotype(geno_sig[i_s], 0);
				int all1 = get_allele_per_haplotype(geno_sig[i_s], 1);

				// Add the errors.
				all0 = (all0 + hap0_err) % 2;
				all1 = (all1 + hap1_err) % 2;

				// Replace the genotype.
				geno_sig[i_s] = all0 + all1 * 2;
			}
		} // i_s loop.
	} // i_var loop.

	// Save the noisy genotypes.
	binarize_variant_genotype_signal_regions(geno_sig_regs, NULL, sample_ids, op_fp);
} // add_noise_to_geno_signal_regions function.

void uniquefy_genotype_signal_regions(char* geno_reg_sig_file, char* sample_list_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(geno_reg_sig_file, sample_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_list_fp);

	fprintf(stderr, "Unique-fying %d variant regions on %d subjects..\n", vecsize(geno_regs), vecsize(sample_ids));

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(geno_regs);

	// Get unique regions on each chromosome separately.
	vector<t_annot_region*>* unique_coord_regs = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < vecsize(restr_geno_regs->chr_ids); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_cat_regs = restr_geno_regs->regions_per_chrom[i_chr];

		// Add the first region.
		unique_coord_regs->push_back(cur_chr_cat_regs->at(0));
		for (int i_reg = 1; i_reg < vecsize(cur_chr_cat_regs); i_reg++)
		{
			// If this region overlaps exactly with previous, remove it.
			if (cur_chr_cat_regs->at(i_reg)->start != cur_chr_cat_regs->at(i_reg - 1)->start ||
				cur_chr_cat_regs->at(i_reg)->end != cur_chr_cat_regs->at(i_reg - 1)->end)
			{
				unique_coord_regs->push_back(cur_chr_cat_regs->at(i_reg));
			}
		} // i_reg loop.
	} // i_chr loop.

	fprintf(stderr, "Saving %d/%d unique genotype regions...\n", vecsize(geno_regs), vecsize(unique_coord_regs));
	binarize_variant_signal_regions_wrapper(unique_coord_regs, sample_ids, op_fp);
}

void concat_variant_wide_genotype_signal_regions_per_geno_sig_list(char* matbed_file_list_fp, char* sample_list_fp, bool remove_unique_coord_vars_flag, char* op_fp)
{
	if (!check_file(matbed_file_list_fp) ||
		!check_file(sample_list_fp))
	{
		fprintf(stderr, "Could not find the list files @ %s, %s\n", matbed_file_list_fp, sample_list_fp);
		exit(1);
	}

	vector<char*>* matbed_files = buffer_file(matbed_file_list_fp);
	vector<char*>* sample_id_list_files = buffer_file(sample_list_fp);

	if (vecsize(matbed_files) != vecsize(sample_id_list_files))
	{
		fprintf(stderr, "Number of genotype matrices do not match number of sample lists.\n");
		exit(1);
	}

	fprintf(stderr, "Concatting %d lists.\n", vecsize(matbed_files));

	vector<char*>* cat_var_sample_ids = new vector<char*>();
	vector<t_annot_region*>* cat_var_geno_sig_regs = new vector<t_annot_region*>();

	for (int i_l = 0; i_l < vecsize(matbed_files); i_l++)
	{
		fprintf(stderr, "Concatting %d. genotype regions: %s (%s)\n", i_l, matbed_files->at(i_l), sample_id_list_files->at(i_l));

		vector<t_annot_region*>* cur_geno_regs = load_variant_signal_regions_wrapper(matbed_files->at(i_l), sample_id_list_files->at(i_l));
		vector<char*>* cur_geno_sample_ids = buffer_file(sample_id_list_files->at(i_l));
			
		if (cat_var_sample_ids->size() == 0)
		{
			cat_var_sample_ids->insert(cat_var_sample_ids->end(), cur_geno_sample_ids->begin(), cur_geno_sample_ids->end());
		}
		else
		{
			if (cat_var_sample_ids->size() != cur_geno_sample_ids->size())
			{
				fprintf(stderr, "Sample id list sizes are not the same: %d, %d\n", vecsize(cat_var_sample_ids), vecsize(cur_geno_sample_ids));

				exit(1);
			}

			for (int i_s = 0; i_s < vecsize(cat_var_sample_ids); i_s++)
			{
				if (!t_string::compare_strings(cat_var_sample_ids->at(i_s), cur_geno_sample_ids->at(i_s)))
				{
					fprintf(stderr, "Sample id's do not match; make sure they are same size and same order.");
					exit(1);
				}
			} // i_s loop.
		}

		fprintf(stderr, "Sample ids are good..\n");

		// Concat the regions.
		cat_var_geno_sig_regs->insert(cat_var_geno_sig_regs->end(), cur_geno_regs->begin(), cur_geno_regs->end());
		fprintf(stderr, "Added %d genotype regions: %d regions for %d subjects..\n", vecsize(cur_geno_regs), vecsize(cat_var_geno_sig_regs), vecsize(cat_var_sample_ids));
	} // i_l loop.

	t_restr_annot_region_list* restr_concat_regs = restructure_annot_regions(cat_var_geno_sig_regs);

	fprintf(stderr, "%d unique chromosomes..\n", vecsize(restr_concat_regs->chr_ids));

	vector<t_annot_region*>* unique_coord_regs = NULL;
	if (remove_unique_coord_vars_flag)
	{
		unique_coord_regs = new vector<t_annot_region*>();

		// Get unique regions on each chromosome separately.
		for (int i_chr = 0; i_chr < vecsize(restr_concat_regs->chr_ids); i_chr++)
		{
			vector<t_annot_region*>* cur_chr_cat_regs = restr_concat_regs->regions_per_chrom[i_chr];

			// Add the first region.
			unique_coord_regs->push_back(cur_chr_cat_regs->at(0));
			for (int i_reg = 1; i_reg < vecsize(cur_chr_cat_regs); i_reg++)
			{
				if (cur_chr_cat_regs->at(i_reg)->start != cur_chr_cat_regs->at(i_reg - 1)->start ||
					cur_chr_cat_regs->at(i_reg)->end != cur_chr_cat_regs->at(i_reg - 1)->end)
				{
					unique_coord_regs->push_back(cur_chr_cat_regs->at(i_reg));
				}
			} // i_reg loop.
		} // i_chr loop.
	}
	else
	{
		unique_coord_regs = cat_var_geno_sig_regs;
	}

	fprintf(stderr, "%d/%d unique variants.\n", (int)unique_coord_regs->size(), (int)cat_var_geno_sig_regs->size());

	//binarize_variant_genotype_signal_regions(unique_coord_regs, NULL, cat_var_sample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(unique_coord_regs, cat_var_sample_ids, op_fp);
} // concat_variant_wide_genotype_signal_regions_per_geno_sig_list function.

void concat_variant_wide_genotype_signal_regions(char* matbed1_fp, char* sample1_list_fp, char* matbed2_fp, char* sample2_list_fp, bool remove_unique_coord_vars_flag, char* op_fp)
{
	vector<t_annot_region*>* geno1_regs = load_variant_signal_regions_wrapper(matbed1_fp, sample1_list_fp);
	vector<char*>* panel1_sample_ids = buffer_file(sample1_list_fp);

	vector<t_annot_region*>* geno2_regs = load_variant_signal_regions_wrapper(matbed2_fp, sample1_list_fp);
	vector<char*>* panel2_sample_ids = buffer_file(sample2_list_fp);

	if (panel1_sample_ids->size() != panel2_sample_ids->size())
	{
		fprintf(stderr, "Sample id's do not match; make sure they are same size and same order.");
		exit(1);
	}

	for (int i_s = 0; i_s < (int)panel1_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(panel1_sample_ids->at(i_s), panel2_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Sample id's do not match; make sure they are same size and same order.");
			exit(1);
		}
	} // i_s loop.

	vector<t_annot_region*>* concat_regs = new vector<t_annot_region*>();
	concat_regs->insert(concat_regs->end(), geno1_regs->begin(), geno1_regs->end());
	concat_regs->insert(concat_regs->end(), geno2_regs->begin(), geno2_regs->end());

	t_restr_annot_region_list* restr_concat_regs = restructure_annot_regions(concat_regs);

	vector<t_annot_region*>* unique_coord_regs = NULL;
	if (remove_unique_coord_vars_flag)
	{
		unique_coord_regs = new vector<t_annot_region*>();

		// Get unique regions on each chromosome separately.
		for (int i_chr = 0; i_chr < vecsize(restr_concat_regs->chr_ids); i_chr++)
		{
			vector<t_annot_region*>* cur_chr_cat_regs = restr_concat_regs->regions_per_chrom[i_chr];

			// Add the first region.
			unique_coord_regs->push_back(cur_chr_cat_regs->at(0));
			for (int i_reg = 1; i_reg < vecsize(cur_chr_cat_regs); i_reg++)
			{
				if (cur_chr_cat_regs->at(i_reg)->start != cur_chr_cat_regs->at(i_reg - 1)->start ||
					cur_chr_cat_regs->at(i_reg)->end != cur_chr_cat_regs->at(i_reg - 1)->end)
				{
					unique_coord_regs->push_back(cur_chr_cat_regs->at(i_reg));
				}
			} // i_reg loop.
		} // i_chr loop.
	}
	else
	{
		unique_coord_regs = concat_regs;
	}

	fprintf(stderr, "%d/%d unique variants.\n", (int)unique_coord_regs->size(), (int)concat_regs->size());

	//binarize_variant_genotype_signal_regions(unique_coord_regs, NULL, panel1_sample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(unique_coord_regs, panel1_sample_ids, op_fp);
} // concat_variant_wide_genotype_signal_regions function.

void extract_hapgen2_files_per_single_reference(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp)
{
	fprintf(stderr, "Saving the hapgen2 files using the reference haplotype data @ %s (%s). Outputs:\n\
	-h %s\n\
	-l %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		h_option_fp,
		l_option_fp);

	/*
	./impute2 \
	-m ./Example/example.chr22.map :: Recombination map.
	-h ./Example/example.chr22.1kG.haps :: List of haplotypes: one haplotype per row: Generate this for the individuals only: Matches the legend (.legend) file.
	-l ./Example/example.chr22.1kG.legend :: rsID position a0 a1 :: rs4821114 20300810 G C
	-Ne 20000 :: Do not change
	-o ./Example/example.chr22.one.phased.impute2 :: Output file.
	*/
	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	//vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_binarized_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	FILE* f_l_option = open_f(l_option_fp, "w");
	FILE* f_h_option = open_f(h_option_fp, "w");

	fprintf(f_l_option, "rsID position a0 a1\n");
	for (int i_chr = 0; i_chr < (int)restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_ref_panel_var_regs->size(); i_reg++)
		{
			// Write the legend entry.			
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
			//char* cur_var_name = toks->at(0)->str();
			char* cur_var_name = cur_chr_ref_panel_var_regs->at(i_reg)->name;
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			fprintf(f_l_option, "%s %d %s %s\n",
				cur_var_name,
				cur_chr_ref_panel_var_regs->at(i_reg)->start,
				cur_var_ref_str,
				cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				if (i_s > 0)
				{
					fprintf(f_h_option, " ");
				}

				fprintf(f_h_option, "%d %d", cur_hap0, cur_hap1);
			} // i_s loop.

			fprintf(f_h_option, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop.
	close_f(f_l_option, l_option_fp);
	close_f(f_h_option, h_option_fp);
}

void convert_HAP_LEGEND_2_haplocoded_matbed(char* hap_fp, char* legend_fp, char* chrom_str, bool add_AF_info_2_id, char* output_matbed_fp)
{
	fprintf(stderr, "Converting %s/%s file with chromosome %s and saving to %s [add_AF_info_2_id=%d]\n", hap_fp, legend_fp, chrom_str, output_matbed_fp, add_AF_info_2_id);

	FILE* f_hap = open_f(hap_fp, "r");
	FILE* f_legend = open_f(legend_fp, "r");

	// Read legend.
	char* legend_line = getline(f_legend); 
	legend_line[0] = 0;
	delete[] legend_line;

	int sample_size = -1;
	char cur_tok[100];
	vector<double>* cur_var_geno = new vector<double>();
	vector<t_annot_region*>* geno_var_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* legend_line = getline(f_legend);
		if (legend_line == NULL)
		{
			break;
		}

		int pos;
		char id[100];
		char ref_all[100];
		char alt_all[100];
		if (sscanf(legend_line, "%s %d %s %s", id, &pos, ref_all, alt_all) != 4)
		{
			fprintf(stderr, "Could not parse legend line: %s\n", legend_line);
			exit(1);
		}

		// Clear the current variant's genotypes list.
		cur_var_geno->clear();

		// Start parsing the genotypes, each 3 entry corresponds to a subject's genotypes.
		char* hap_line = getline(f_hap);
		if (hap_line == NULL)
		{
			fprintf(stderr, "hap file ended before legend file..\n");
			exit(1);
		}

		int cur_char_i = 0;
		while (1)
		{
			if (!t_string::get_next_token(hap_line, cur_tok, 100, " ", cur_char_i))
			{
				// We are finished.
				break;
			}
			int all0 = atoi(cur_tok);

			if (!t_string::get_next_token(hap_line, cur_tok, 100, " ", cur_char_i))
			{
				fprintf(stderr, "Could not parse geno line: %s\n", hap_line);
				exit(1);
			}
			int all1 = atoi(cur_tok);

			if (all0 != 0 && all0 != 1)
			{
				fprintf(stderr, "Illegal allele.\n");
				exit(1);
			}

			if (all1 != 0 && all1 != 1)
			{
				fprintf(stderr, "Illegal allele.\n");
				exit(1);
			}

			cur_var_geno->push_back(all1 * 2 + all0);
		} // line parsing.

		// First, check the assignment of sample size.
		if (sample_size == -1)
		{
			fprintf(stderr, "Setting sample size to %d\n", vecsize(cur_var_geno));
			sample_size = vecsize(cur_var_geno);
		}
		else if (sample_size != vecsize(cur_var_geno))
		{
			fprintf(stderr, "Could not parse consistent # of subjects: %d/%d\n", sample_size, vecsize(cur_var_geno));
			exit(1);
		}

		// Copy the genotype.
		double total_geno = 0;
		char* geno_sig = new char[sample_size + 2];
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			int genocoded_geno = (get_allele_per_haplotype(cur_var_geno->at(i_s), 0) + get_allele_per_haplotype(cur_var_geno->at(i_s), 1));
			geno_sig[i_s] = cur_var_geno->at(i_s);
			total_geno += genocoded_geno;
		} // i_s loop.

		double AAF = (total_geno / (2 * sample_size));

		// We make no assumptions on the identifier. It is required to have 4 entries and it may already be as needed.
		char full_id[1000];
		if (add_AF_info_2_id)
		{
			sprintf(full_id, "%s_%s_%s_%.4f", id, ref_all, alt_all, AAF);
		}
		else
		{
			strcpy(full_id, id);
		}		

		// Copy region.
		t_annot_region* cur_var_reg = get_empty_region();
		cur_var_reg->chrom = t_string::copy_me_str(chrom_str);
		cur_var_reg->start = pos; // gen file positions are 1 based.
		cur_var_reg->end = pos; // gen file positions are 1 based.
		cur_var_reg->name = t_string::copy_me_str(full_id);
		cur_var_reg->strand = '+';

		void** reg_info = new void* [10];
		reg_info[0] = geno_sig;
		cur_var_reg->data = reg_info;

		geno_var_regs->push_back(cur_var_reg);

		delete[] hap_line;
	} // file reading loop.

	close_f(f_legend, legend_fp);
	close_f(f_hap, hap_fp);

	// save the matbed.
	vector<char*>* sample_ids = new vector<char*>();
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		char cur_sample_id[1000];
		sprintf(cur_sample_id, "SUBJ_%d", i_s);
		sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
	} // i_s loop.

	fprintf(stderr, "Sorting and saving the %d variant genotype matrix on %d subjects to %s..\n", vecsize(geno_var_regs), vecsize(sample_ids), output_matbed_fp);
	sort(geno_var_regs->begin(), geno_var_regs->end(), sort_regions);
	binarize_variant_signal_regions_wrapper(geno_var_regs, sample_ids, output_matbed_fp);
} // convert_HAP_LEGEND_2_haplocoded_matbed function.


void convert_GEN_2_matbed(char* gen_fp, char* chrom_str, char* output_matbed_fp)
{
	fprintf(stderr, "Converting GEN file to %s with chromosome %s and saving to %s\n", gen_fp, chrom_str, output_matbed_fp);

	FILE* f_gen = open_f(gen_fp, "r");

	int sample_size = -1;
	char cur_tok[100];
	vector<double>* cur_var_geno = new vector<double>();
	vector<t_annot_region*>* geno_var_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* gen_line = getline(f_gen);
		if (gen_line == NULL)
		{
			break;
		}

		// Clear the current variant's genotypes list.
		cur_var_geno->clear();

		// Read the legend entries: ID, pos, alleles.
		int pos;
		char id[100];
		char ref_all[100];
		char alt_all[100];
		int cur_char_i = 0;
		t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i);
		t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i);
		strcpy(id, cur_tok);
		t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i);
		pos = atoi(cur_tok);
		t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i);
		strcpy(ref_all, cur_tok);
		t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i);
		strcpy(alt_all, cur_tok);

		// Start parsing the genotypes, each 3 entry corresponds to a subject's genotypes.
		double geno_probs[3];
		while (1)
		{			
			if (!t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i))
			{
				// We are finished.
				break;
			}
			geno_probs[0] = atof(cur_tok);

			if (!t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i))
			{
				fprintf(stderr, "Could not parse geno line: %s\n", gen_line);
				exit(1);
			}
			geno_probs[1] = atof(cur_tok);


			if (!t_string::get_next_token(gen_line, cur_tok, 100, " ", cur_char_i))
			{
				fprintf(stderr, "Could not parse geno line: %s\n", gen_line);
				exit(1);
			}
			geno_probs[2] = atof(cur_tok);

			if (MAX(geno_probs[0], MAX(geno_probs[1], geno_probs[2])) != 1)
			{
				fprintf(stderr, "Could not find max entry of 1 @ snp_pos=%d in %s\n", pos, gen_fp);
				exit(1);
			}

			if (geno_probs[0] == 1)
			{
				cur_var_geno->push_back(0);
			}
			else if (geno_probs[1] == 1)
			{
				cur_var_geno->push_back(1);
			}
			else if (geno_probs[2] == 1)
			{
				cur_var_geno->push_back(2);
			}
			else
			{
				fprintf(stderr, "Could not parse geno: snp_pos=%d in %s\n", pos, gen_fp);
				exit(1);
			}
		} // line parsing.

		// First, check the assignment of sample size.
		if (sample_size == -1)
		{
			fprintf(stderr, "Setting sample size to %d\n", vecsize(cur_var_geno));
			sample_size = vecsize(cur_var_geno);
		}
		else if (sample_size != vecsize(cur_var_geno))
		{
			fprintf(stderr, "Could not parse consistent # of subjects: %d/%d\n", sample_size, vecsize(cur_var_geno));
			exit(1);
		}

		// Copy the genotype.
		double total_geno = 0;
		char* geno_sig = new char[sample_size + 2];
		for (int i_s = 0; i_s < sample_size; i_s++)
		{
			geno_sig[i_s] = cur_var_geno->at(i_s);
			total_geno += cur_var_geno->at(i_s);
		} // i_s loop.

		double AAF = (total_geno / (2 * sample_size));

		char full_id[1000];
		sprintf(full_id, "%s_%s_%s_%.4f", id, ref_all, alt_all, AAF);


		// Copy region.
		t_annot_region* cur_var_reg = get_empty_region();
		cur_var_reg->chrom = t_string::copy_me_str(chrom_str);
		cur_var_reg->start = pos; // gen file positions are 1 based.
		cur_var_reg->end = pos; // gen file positions are 1 based.
		cur_var_reg->name = t_string::copy_me_str(full_id);
		cur_var_reg->strand = '+';

		void** reg_info = new void* [10];
		reg_info[0] = geno_sig;
		cur_var_reg->data = reg_info;

		geno_var_regs->push_back(cur_var_reg);

		delete[] gen_line;
	} // file reading loop.

	close_f(f_gen, gen_fp);

	// save the matbed.
	vector<char*>* sample_ids = new vector<char*>();
	for (int i_s = 0; i_s < sample_size; i_s++)
	{
		char cur_sample_id[1000];
		sprintf(cur_sample_id, "SUBJ_%d", i_s);
		sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
	} // i_s loop.

	fprintf(stderr, "Sorting and saving the %d variant genotype matrix on %d subjects to %s..\n", vecsize(geno_var_regs), vecsize(sample_ids), output_matbed_fp);
	sort(geno_var_regs->begin(), geno_var_regs->end(), sort_regions);
	binarize_variant_signal_regions_wrapper(geno_var_regs, sample_ids, output_matbed_fp);
} // convert_GEN_2_matbed function.

void extract_hapgen2_files(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp)
{
	fprintf(stderr, "Saving the hapgen2 files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-h %s\n\
	-l %s\n\
	-g %s\n\
	-strand_g %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		h_option_fp,
		l_option_fp,
		g_option_fp,
		strand_g_option_fp);

	/*
	./impute2 \
	-m ./Example/example.chr22.map :: Recombination map.
	-h ./Example/example.chr22.1kG.haps :: List of haplotypes: one haplotype per row: Generate this for the individuals only: Matches the legend (.legend) file.
	-l ./Example/example.chr22.1kG.legend :: rsID position a0 a1 :: rs4821114 20300810 G C
	-g ./Example/example.chr22.study.gens :: Known genotypes with positions: Must be sorted with respect to position.
	-strand_g ./Example/example.chr22.study.strand :: Strand of each variant in '-g' option.
	-int 20.4e6 20.5e6 :: The start and end positions.
	-Ne 20000 :: Do not change
	-o ./Example/example.chr22.one.phased.impute2 :: Output file.
	*/
	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	//vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_binarized_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	FILE* f_l_option = open_f(l_option_fp, "w");
	FILE* f_h_option = open_f(h_option_fp, "w");

	fprintf(f_l_option, "rsID position a0 a1\n");
	for (int i_chr = 0; i_chr < (int)restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_ref_panel_var_regs->size(); i_reg++)
		{
			// Write the legend entry.			
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
			//char* cur_var_name = toks->at(0)->str();
			char* cur_var_name = cur_chr_ref_panel_var_regs->at(i_reg)->name;
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			fprintf(f_l_option, "%s %d %s %s\n",
				cur_var_name,
				cur_chr_ref_panel_var_regs->at(i_reg)->start,
				cur_var_ref_str,
				cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			for (int i_s = 0; i_s < (int)reference_haplo_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				if (i_s > 0)
				{
					fprintf(f_h_option, " ");
				}

				fprintf(f_h_option, "%d %d", cur_hap0, cur_hap1);
			} // i_s loop.

			fprintf(f_h_option, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop.
	close_f(f_l_option, l_option_fp);
	close_f(f_h_option, h_option_fp);

	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", (int)input_geno_sample_ids->size());
	//vector<t_annot_region*>* input_geno_sig_regs = load_binarized_variant_genotype_signal_regions(input_genotype_matrix_matbed_fp, input_geno_sample_ids);
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_signal_regions_wrapper(input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	FILE* f_g_option = open_f(g_option_fp, "w");
	FILE* f_strand_g_option = open_f(strand_g_option_fp, "w");
	for (int i_chr = 0; i_chr < (int)restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			// Write the legend entry.
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the strand.
			fprintf(f_strand_g_option, "%d +\n", cur_chr_input_geno_sig_regs->at(i_reg)->start);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			fprintf(f_g_option, "%s %s %d %s %s",
				cur_var_name, cur_var_name,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_ref_str, cur_var_alt_str);

			for (int i_s = 0; i_s < (int)input_geno_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				int geno = cur_hap0 + cur_hap1;

				if (geno == 0)
				{
					fprintf(f_g_option, " 1 0 0");
				}
				else if (geno == 1)
				{
					fprintf(f_g_option, " 0 1 0");
				}
				else if (geno == 2)
				{
					fprintf(f_g_option, " 0 0 1");
				}
				else
				{
					fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
						cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
						(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

					exit(1);
				}

			} // i_s loop.

			fprintf(f_g_option, "\n");
		} // i_reg loop.
	} // i_chr loop. 
	close_f(f_g_option, g_option_fp);
	close_f(f_strand_g_option, strand_g_option_fp);
}


void replace_variant_alleles_per_genotypes_signals(char* geno_regs_BED_fp, char* sample_ids_list_fp, char* new_alleles_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(geno_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d genotype regions with %d samples.\n", (int)geno_regs->size(), (int)sample_ids->size());

	// Set the scores.
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		geno_regs->at(i_reg)->score = 0;
	}

	vector<t_annot_region*>* new_allele_regs = load_BED(new_alleles_BED_fp);
	fprintf(stderr, "Loaded %d new allele regions.\n", (int)new_allele_regs->size());

	vector<t_annot_region*>* intersects = intersect_annot_regions(geno_regs, new_allele_regs, false);
	fprintf(stderr, "Processing %d intersections.\n", (int)intersects->size());

	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* geno_reg = int_info->src_reg;
		t_annot_region* new_allele_reg = int_info->dest_reg;

		vector<char*>* geno_reg_name_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(geno_reg->name, "_"));
		vector<char*>* new_allele_reg_name_toks = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(new_allele_reg->name, "_"));

		int n_geno_toks = (int)geno_reg_name_toks->size();
		int n_new_allele_toks = (int)new_allele_reg_name_toks->size();

		if (n_geno_toks != n_new_allele_toks)
		{
			fprintf(stderr, "Sanity check failed, name formatting does not match: %s, %s\n", geno_reg->name, new_allele_reg->name);
			exit(1);
		}

		char allele_codings[10];
		allele_codings[0] = new_allele_reg_name_toks->at(n_new_allele_toks - 2)[0];
		allele_codings[1] = new_allele_reg_name_toks->at(n_new_allele_toks - 1)[0];
		char allele_values[10];
		allele_values[0] = new_allele_reg_name_toks->at(n_new_allele_toks - 2)[2];
		allele_values[1] = new_allele_reg_name_toks->at(n_new_allele_toks - 1)[2];

		if (new_allele_reg->strand == '-')
		{
			allele_values[0] = get_complementary_dna_nuc_per_dna_nuc(allele_values[0]);
			allele_values[1] = get_complementary_dna_nuc_per_dna_nuc(allele_values[1]);
		}

		char new_ref_allele = 0;
		char new_alt_allele = 0;
		if (geno_reg_name_toks->at(n_geno_toks - 2)[0] == allele_codings[0])
		{
			new_ref_allele = allele_values[0];
		}
		else if (geno_reg_name_toks->at(n_geno_toks - 2)[0] == allele_codings[1])
		{
			new_ref_allele = allele_values[1];
		}
		else
		{
			fprintf(stderr, "Could not map reference allele (%c): %s ;; %s\n", geno_reg_name_toks->at(n_geno_toks - 2)[0], geno_reg->name, new_allele_reg->name);
			exit(1);
		}

		// Map alternate allele.
		if (geno_reg_name_toks->at(n_geno_toks - 1)[0] == allele_codings[0])
		{
			new_alt_allele = allele_values[0];
		}
		else if (geno_reg_name_toks->at(n_geno_toks - 1)[0] == allele_codings[1])
		{
			new_alt_allele = allele_values[1];
		}
		else
		{
			fprintf(stderr, "Could not map alternate allele (%c): %s ;; %s\n", geno_reg_name_toks->at(n_geno_toks - 1)[0], geno_reg->name, new_allele_reg->name);
			exit(1);
		}

		if (new_ref_allele == new_alt_allele)
		{
			fprintf(stderr, "Sanity check failed, mapped ref/alt to same allele: %c, %c:: %s ;; %s\n", new_ref_allele, new_alt_allele, 
					geno_reg->name, new_allele_reg->name);

			exit(1);
		}

		// Everything looks good, replace new ref/alt alleles, then move on.
		geno_reg_name_toks->at(n_geno_toks - 2)[0] = new_ref_allele;
		geno_reg_name_toks->at(n_geno_toks - 1)[0] = new_alt_allele;

		// Concatenate to get the final identifier.
		char concat_new_geno_id[1000];
		strcpy(concat_new_geno_id, geno_reg_name_toks->at(0));

		for (int i_tok = 1; i_tok < n_geno_toks; i_tok++)
		{
			strcat(concat_new_geno_id, "_");
			strcat(concat_new_geno_id, geno_reg_name_toks->at(i_tok));
		} // i_tok loop.

		geno_reg->score = 1;
		fprintf(stderr, "New id: %s:%d::%s / %s\n", geno_reg->chrom, geno_reg->start, concat_new_geno_id, geno_reg->name);
		geno_reg->name = t_string::copy_me_str(concat_new_geno_id);
	} // i_int loop.

	vector<t_annot_region*>* reset_geno_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		if (geno_regs->at(i_reg)->score == 1)
		{
			reset_geno_regs->push_back(geno_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Saving allele recoded regions for %d/%d regions.\n", (int)reset_geno_regs->size(), (int)geno_regs->size());

	binarize_variant_genotype_signal_regions(reset_geno_regs, NULL, sample_ids, op_fp);
} // replace_variant_alleles_per_genotypes_signals 

double get_recomb_rate_difference_per_regions(t_restr_annot_region_list* restr_recomb_rate_regs, 
											t_annot_region* reg1, t_annot_region* reg2)
{
	if (!t_string::compare_strings(reg1->chrom, reg2->chrom))
	{
		return(1);
	}

	int i_chr = t_string::get_i_str(restr_recomb_rate_regs->chr_ids, reg1->chrom);

	t_annot_region* left_reg = NULL;
	t_annot_region* right_reg = NULL;
	if (reg1->start < reg2->start)
	{
		left_reg = reg1;
		right_reg = reg2;
	}
	else
	{
		left_reg = reg2;
		right_reg = reg1;
	}
	
	// Overlap the left region:
	t_annot_region* left_oling_recomb_reg = NULL;
	if (left_reg->start <= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->start)
	{
		left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0);
	}
	else if(left_reg->start >= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->end)
	{
		left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->back();
	}
	else
	{
		int reg_i = locate_posn_region_per_region_starts(left_reg->start, restr_recomb_rate_regs->regions_per_chrom[i_chr], 0, (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size());
		while (reg_i > 0 &&
			reg_i < (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_end >= left_reg->start)
		{
			reg_i--;
		} // reg_i loop.

		while (reg_i >= 0 &&
			reg_i < (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_start <= left_reg->end)
		{
			int ol_start = MAX(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->start, left_reg->start);
			int ol_end = MIN(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->end, left_reg->end);
			if (ol_end >= ol_start)
			{
				// Found overlap, set it.
				left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Found right reg @ %s:%d:: %s:%d-%d\n",
						left_reg->chrom, left_reg->start,
						left_oling_recomb_reg->chrom, left_oling_recomb_reg->start, left_oling_recomb_reg->end);
				}

				break;
			}
			reg_i++;
		} // reg_i loop.
	}

	// Overlap the right region:
	t_annot_region* right_oling_recomb_reg = NULL;
	if (right_reg->start <= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->start)
	{
		right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0);
	}
	else if (right_reg->start >= restr_recomb_rate_regs->regions_per_chrom[i_chr]->back()->end)
	{
		right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->back();
	}
	else
	{
		int reg_i = locate_posn_region_per_region_starts(right_reg->start, restr_recomb_rate_regs->regions_per_chrom[i_chr], 0, (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size());
		while (reg_i > 0 &&
			reg_i < (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_end >= right_reg->start)
		{
			reg_i--;
		} // reg_i loop.

		while (reg_i >= 0 &&
			reg_i < (int)restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_start <= right_reg->end)
		{
			int ol_start = MAX(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->start, right_reg->start);
			int ol_end = MIN(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->end, right_reg->end);
			if (ol_end >= ol_start)
			{
				// Found overlap, set it.
				right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Found right reg @ %s:%d:: %s:%d-%d\n",
						right_reg->chrom, right_reg->start,
						right_oling_recomb_reg->chrom, right_oling_recomb_reg->start, right_oling_recomb_reg->end);
				}

				break;
			}
			reg_i++;
		} // reg_i loop.
	}

	double recomb_diff = -1;
	if (left_oling_recomb_reg != NULL &&
		right_oling_recomb_reg != NULL)
	{		
		double* first_reg_recomb_info = (double*)(left_oling_recomb_reg->data);
		double* last_reg_recomb_info = (double*)(right_oling_recomb_reg->data);
		recomb_diff = last_reg_recomb_info[0] - first_reg_recomb_info[1];

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "delta recomb(%s:%d, %s:%d)::[%s:%d-%d (%.4f) ---> %s:%d-%d (%.4f)]\n",
				reg1->chrom, reg1->start, reg2->chrom, reg2->start,
				left_oling_recomb_reg->chrom, left_oling_recomb_reg->start, left_oling_recomb_reg->end, first_reg_recomb_info[1],
				right_oling_recomb_reg->chrom, right_oling_recomb_reg->start, right_oling_recomb_reg->end, last_reg_recomb_info[0]);
		}
	}
	else
	{
		fprintf(stderr, "Could not match one of the regions for delta recomb(%s:%d, %s:%d)\n",
			reg1->chrom, reg1->start, reg2->chrom, reg2->start);
		exit(1);
	}

	return(recomb_diff);
}

double get_recomb_prob_per_delta_cM(double delta_cM)
{
	double recomb_prob = 0.5 * (1 - exp(-2 * delta_cM / 100));
	return(recomb_prob);
}

vector<t_annot_region*>* load_recombination_rate_regions_per_map_file(char* recombination_rate_dir)
{
	char recomb_rate_chr_ids_fp[1000];
	sprintf(recomb_rate_chr_ids_fp, "%s/chr_ids.txt", recombination_rate_dir);
	vector<char*>* recomb_chr_ids = buffer_file(recomb_rate_chr_ids_fp);
	if (recomb_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from recombination directory @ %s\n", recomb_rate_chr_ids_fp);
		exit(1);
	}

	fprintf(stderr, "Loaded %d chr ids for the recombination rate.\n", (int)recomb_chr_ids->size());
	vector<t_annot_region*>* recomb_rate_regs = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < (int)recomb_chr_ids->size(); i_chr++)
	{
		char cur_chr_recomb_rate_fp[1000];
		sprintf(cur_chr_recomb_rate_fp, "%s/%s.map", recombination_rate_dir, recomb_chr_ids->at(i_chr));

		fprintf(stderr, "Loading the recombination rates on %s from %s\n", recomb_chr_ids->at(i_chr), cur_chr_recomb_rate_fp);
		FILE* f_recomb_map = open_f(cur_chr_recomb_rate_fp, "r");
		vector<t_annot_region*>* cur_chr_recomb_regs = new vector<t_annot_region*>();
		while (1)
		{
			char* cur_line = getline(f_recomb_map);
			if (cur_line == NULL)
			{
				break;
			}

			// 22      .       0.000000        16051347
			char read_chr_id[100];
			double cM = 0;
			int cur_posn = 0;
			if (sscanf(cur_line, "%s %*s %lf %d", read_chr_id, &cM, &cur_posn) != 3)
			{
				fprintf(stderr, "Could not parse: %s\n", cur_line);
				exit(1);
			}

			t_annot_region* cur_recomb_reg = get_empty_region();
			cur_recomb_reg->chrom = t_string::copy_me_str(read_chr_id);
			cur_recomb_reg->start = cur_posn;
			cur_recomb_reg->end = cur_posn;
			cur_recomb_reg->dbl_score = cM;
			cur_recomb_reg->strand = '+';

			recomb_rate_regs->push_back(cur_recomb_reg);
			cur_chr_recomb_regs->push_back(cur_recomb_reg);
			delete[] cur_line;
		} // file reading loop.
		fclose(f_recomb_map);

		fprintf(stderr, "Resetting the coordinates and the recombination info for the %d recombination regions.\n", (int)cur_chr_recomb_regs->size());

		// Fix the coordinates by neighboring regions positions.
		sort(cur_chr_recomb_regs->begin(), cur_chr_recomb_regs->end(), sort_regions);

		for (int i_reg = 1; i_reg < (int)cur_chr_recomb_regs->size(); i_reg++)
		{
			cur_chr_recomb_regs->at(i_reg)->start = cur_chr_recomb_regs->at(i_reg - 1)->end+1;

			// Setup the recombination rates for this region.
			double* cur_reg_recomb_info = new double[2];
			cur_chr_recomb_regs->at(i_reg)->data = cur_reg_recomb_info;
			cur_reg_recomb_info[0] = cur_chr_recomb_regs->at(i_reg)->dbl_score;
			cur_reg_recomb_info[1] = cur_chr_recomb_regs->at(i_reg-1)->dbl_score;
		} // i_reg loop.
	} // i_chr loop.

	fprintf(stderr, "Loaded %d recombination regions in total.\n", (int)recomb_rate_regs->size());

	return(recomb_rate_regs);
}

// The "unstructured" functions are useful if the code is developed with MatRegSubj format (rather than usual signal-per-region construct). This should provide a faster interface.
void save_MatSubjReg_unstructured(vector<t_annot_region*>* var_regs, vector<char*>* subj_ids, const unsigned char* geno_mat, bool compress, char* op_prefix)
{
	char subj_id_fp[1000];
	sprintf(subj_id_fp, "%s_subjects.list", op_prefix);
	save_lines(subj_ids, subj_id_fp);

	char regs_BED_fp[1000];
	sprintf(regs_BED_fp, "%s_variants.bed", op_prefix);
	dump_BED(regs_BED_fp, var_regs);

	size_t n_bytes_2_write = vecsize_t(var_regs) * vecsize_t(subj_ids);
	if (compress)
	{
		char geno_mat_fp[1000];
		sprintf(geno_mat_fp, "%s_genotypes.matrix.gz", op_prefix);
		compressSaveBuffer(geno_mat, n_bytes_2_write, geno_mat_fp);
	}
	else
	{
		char geno_mat_fp[1000];
		sprintf(geno_mat_fp, "%s_genotypes.matrix", op_prefix);

		FILE* f_geno_mat = open_f(geno_mat_fp, "wb");
		fwrite(geno_mat, sizeof(char), n_bytes_2_write, f_geno_mat);
		close_f(f_geno_mat, geno_mat_fp);
	}
}


// This is a fast accession function for mid-large size genotype matrices. Skip assigning to regions.
unsigned char* load_MatSubjReg_unstructured(char* matsubjreg_prefix, vector<char*>* subj_ids, vector<t_annot_region*>* var_regs)
{
	auto zlib_load_start_chrono = std::chrono::high_resolution_clock::now();
	char sample_ids_list_fp[1000];
	sprintf(sample_ids_list_fp, "%s_subjects.list", matsubjreg_prefix);
	if (!check_file(sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find %s\n", sample_ids_list_fp);
		return(NULL);
	}
	vector<char*>* read_subj_ids = buffer_file(sample_ids_list_fp);
	subj_ids->insert(subj_ids->end(), read_subj_ids->begin(), read_subj_ids->end());

	char regs_bed_fp[1000];
	sprintf(regs_bed_fp, "%s_variants.bed", matsubjreg_prefix);
	if (!check_file(sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find %s\n", sample_ids_list_fp);
		return(NULL);
	}
	vector<t_annot_region*>* read_zlib_regs = load_BED(regs_bed_fp);
	var_regs->insert(var_regs->end(), read_zlib_regs->begin(), read_zlib_regs->end());

	char geno_mat_fp[1000];
	sprintf(geno_mat_fp, "%s_genotypes.matrix.gz", matsubjreg_prefix);
	if (!check_file(geno_mat_fp))
	{
		fprintf(stderr, "could not find geontype matrix..\n");
		exit(1);
	}

	size_t geno_mat_size = (size_t)(vecsize(subj_ids)) * (size_t)(vecsize(read_zlib_regs));
	unsigned char* geno_matrix = new unsigned char[geno_mat_size];
	size_t n_bytes_read = 0;
	decompressGzipFile(geno_mat_fp, geno_matrix, n_bytes_read);

	auto zlib_load_end_chrono = std::chrono::high_resolution_clock::now();

	if (n_bytes_read != geno_mat_size)
	{
		fprintf(stderr, "Sanity check failed: zlib did not read expected bytes: %lu vs %lu\n", geno_mat_size, n_bytes_read);
		exit(1);
	}

	std::chrono::duration<double> zlib_load_duration = zlib_load_end_chrono - zlib_load_start_chrono;
	fprintf(stderr, "zlib loading finished in %.3f seconds..\n", zlib_load_duration.count());

	return(geno_matrix);
}

void binarize_variant_signal_regions_wrapper(vector<t_annot_region*>* roi_regs_w_signals, vector<char*>* sample_ids, char* op_fp)
{
	if (t_string::compare_strings(op_fp, "stdin") ||
		t_string::ends_with(op_fp, ".bed") ||
		t_string::ends_with(op_fp, ".txt") ||
		t_string::ends_with(op_fp, ".txt.gz") ||
		t_string::ends_with(op_fp, ".bed.gz"))
	{
		dump_geno_sig_regs_plain(roi_regs_w_signals, sample_ids, false, op_fp);
	}
	else if (t_string::ends_with(op_fp, ".bedmat") ||
		t_string::ends_with(op_fp, ".matbed") ||
		t_string::ends_with(op_fp, ".matbed.gz") ||
		t_string::ends_with(op_fp, ".bedmat.gz") ||
		t_string::ends_with(op_fp, ".bedsig") ||
		t_string::ends_with(op_fp, ".bedsig.gz") ||
		t_string::ends_with(op_fp, ".sigbed") ||
		t_string::ends_with(op_fp, ".sigbed.gz"))
	{
		binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
	}
	else
	{
		// Treat the output as MatSubjReg formatted.
		//binarize_variant_genotype_signal_regions_per_matrix_subjects_regions(roi_regs_w_signals, sample_ids, true, op_fp);
#define USE_MT_SAVE_MSR
//#undef USE_MT_SAVE_MSR

#ifdef USE_MT_SAVE_MSR
		int n_saving_threads = 20;
		binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(roi_regs_w_signals, sample_ids, true, n_saving_threads, op_fp);
#endif

#ifndef USE_MT_SAVE_MSR
		binarize_variant_genotype_signal_regions_per_matrix_subjects_regions(roi_regs_w_signals, sample_ids, true, op_fp);
#endif
	}
}

vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp)
{
	vector<t_annot_region*>* var_sig_regs = NULL;

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	if (sample_ids == NULL || !check_file(sample_ids_list_fp))
	{
		fprintf(stderr, "Could not load the sample ids from %s", sample_ids_list_fp);
		exit(1);
	}

	if (t_string::compare_strings(geno_sig_regs_BED_fp, "stdin") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed") || 
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt.gz") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed.gz"))
	{
		var_sig_regs = load_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else if (t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed.gz"))
	{
		var_sig_regs = load_binarized_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else
	{
		// Check MatSubjReg formatted data.
		char subj_test_fp[1000];
		sprintf(subj_test_fp, "%s_subjects.list", geno_sig_regs_BED_fp);
		char var_regs_test_fp[1000];
		sprintf(var_regs_test_fp, "%s_variants.bed", geno_sig_regs_BED_fp);
		char geno_mat_test_fp[1000];

		if (check_file(var_regs_test_fp) &&
			check_file(subj_test_fp))
		{
			bool valid_MRS = false;
			sprintf(geno_mat_test_fp, "%s_genotypes.matrix", geno_sig_regs_BED_fp);
			if (check_file(geno_mat_test_fp))
			{
				vector<char*>* loaded_sample_ids = new vector<char*>();
				var_sig_regs = load_variant_genotype_signal_regions_per_matrix_subjects_regions(geno_sig_regs_BED_fp, loaded_sample_ids);
				t_string::clean_string_list(loaded_sample_ids);
				valid_MRS = true;
			}

			sprintf(geno_mat_test_fp, "%s_genotypes.matrix.gz", geno_sig_regs_BED_fp);
			if (check_file(geno_mat_test_fp))
			{
				vector<char*>* loaded_sample_ids = new vector<char*>();
				var_sig_regs = load_variant_genotype_signal_regions_per_matrix_subjects_regions(geno_sig_regs_BED_fp, loaded_sample_ids);
				t_string::clean_string_list(loaded_sample_ids);
				valid_MRS = true;
			}

			if (!valid_MRS)
			{
				fprintf(stderr, "%s(%d): Could not find the genotypes matrix for %s\n", __FILE__, __LINE__, geno_sig_regs_BED_fp);
				exit(1);
			}
		}
		else
		{
			fprintf(stderr, "%s(%d): Could not identify variant signal regions file type: %s\n", __FILE__, __LINE__, geno_sig_regs_BED_fp);
			exit(1);
		}
	}

	t_string::clean_string_list(sample_ids);

	return(var_sig_regs);
}

void generate_random_sampled_genotype_matrix(char* genotype_regs_BED_fp, char* sample_ids_list_fp, 
	int n_samples_2_generate, 
	char* op_fp)
{
	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(genotype_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_regs->size(), (int)sample_ids->size());
	for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
	{
		void** reg_info = (void**)(geno_regs->at(i_reg)->data);
		char* sampled_geno_sig = new char[n_samples_2_generate + 2];
		memset(sampled_geno_sig, 0, n_samples_2_generate * sizeof(char));
		char* ref_geno_sig = (char*)(reg_info[0]);
		reg_info[0] = sampled_geno_sig;
		reg_info[1] = ref_geno_sig;
	} // i_reg loop.

	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	vector<char*>* generated_sample_ids = new vector<char*>();
	for (int i_s = 0; i_s < n_samples_2_generate; i_s++)
	{
		fprintf(stderr, "Generating %d. samples     \r", n_samples_2_generate);
		for (int i_reg = 0; i_reg < (int)geno_regs->size(); i_reg++)
		{
			int rand_i_s = floor(rng->random_double_ran3() * (int)sample_ids->size());
			if (rand_i_s >= (int)sample_ids->size())
			{
				rand_i_s = (int)sample_ids->size() - 1;
			}

			void** cur_reg_info = (void**)(geno_regs->at(i_reg)->data);
			char* cur_geno = (char*)(cur_reg_info[1]);
			char* sampled_cur_geno = (char*)(cur_reg_info[0]);
			sampled_cur_geno[i_s] = cur_geno[rand_i_s];
		} // i_reg loop.

		char cur_sample_id[1000];
		sprintf(cur_sample_id, "sample_%d", i_s);
		generated_sample_ids->push_back(t_string::copy_me_str(cur_sample_id));
	} // i_s loop.

	fprintf(stderr, "Saving generated sample of %d individuals.\n", (int)generated_sample_ids->size());
	binarize_variant_genotype_signal_regions(geno_regs, NULL, generated_sample_ids, op_fp);
}

void merge_samples_per_regions_list(char* regions_list_fp,
									char* samples_list_fp,
									bool match_region_names,
									char* regions_op_fp,
									char* samples_op_fp)
{
	vector<char*>* regions_list = buffer_file(regions_list_fp);
	fprintf(stderr, "Pooling %d regions.\n", (int)regions_list->size());

	vector<char*>* samples_list = buffer_file(samples_list_fp);
	fprintf(stderr, "Pooling %d sample lists.\n", (int)samples_list->size());

	if ((int)samples_list->size() != (int)regions_list->size())
	{
		fprintf(stderr, "Samples list does not match to regions list.\n");
		exit(1);
	}

	vector<char*>* pooled_sample_ids = new vector<char*>();
	for (int i_f = 0; i_f < (int)regions_list->size(); i_f++)
	{
		vector<char*>* cur_sample_ids = buffer_file(samples_list->at(i_f));
		pooled_sample_ids->insert(pooled_sample_ids->end(), cur_sample_ids->begin(), cur_sample_ids->end());
	} // i_f loop.

	fprintf(stderr, "Loaded %d samples, pooling.\n", (int)pooled_sample_ids->size());

	fprintf(stderr, "Loading first sample as the reference panel...\n");
	vector<t_annot_region*>* geno_regs = load_variant_signal_regions_wrapper(regions_list->at(0), samples_list->at(0));
	vector<char*>* first_samples_list = buffer_file(samples_list->at(0));
	int max_coded_geno = get_max_genotype_value(geno_regs, first_samples_list);
	
	size_t pooled_signal_mempool_size = vecsize_t(geno_regs) * vecsize_t(pooled_sample_ids);
	fprintf(stderr, "Re-allocating pooled genotype signals [%.3f Gb]..\n", pooled_signal_mempool_size / (1024.0 * 1024.0 * 1024.0));
	char* pooled_signal_mempool = new char[pooled_signal_mempool_size + 2];
	memset(pooled_signal_mempool, 0, pooled_signal_mempool_size);
	fprintf(stderr, "Assigning the pooled signal buffers to regions..\n");
	for (size_t i_reg = 0; i_reg < vecsize_t(geno_regs); i_reg++)
	{
		void** cur_reg_info = (void**)(geno_regs->at(i_reg)->data);
		cur_reg_info[0] = pooled_signal_mempool + (i_reg * vecsize_t(pooled_sample_ids));
		//char* cur_geno = (char*)(cur_reg_info[0]);
		//char* pooled_geno = new char[(int)pooled_sample_ids->size() + 2];
		//memset(pooled_geno, 0, sizeof(char) * (int)pooled_sample_ids->size());

		// Don't delete this signal because it may be allocated as a chunk. This will create some amount of leak, we should be able to live with this.
		//delete[] cur_geno;

		// Replace the genotypes.
		//cur_reg_info[0] = pooled_geno;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d regions as the reference.\n", (int)geno_regs->size());

	vector<char*>* samples_so_far = new vector<char*>();
	for (int i_f = 0; i_f < (int)regions_list->size(); i_f++)
	{
		fprintf(stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		fprintf(stderr, "@ %d/%d sample;\n", i_f, vecsize(regions_list));
		vector<t_annot_region*>* cur_geno_regs = load_variant_signal_regions_wrapper(regions_list->at(i_f), samples_list->at(i_f));
		vector<char*>* cur_samples_list = buffer_file(samples_list->at(i_f));

		int cur_max_coded_geno = get_max_genotype_value(cur_geno_regs, cur_samples_list);

		if (cur_max_coded_geno != max_coded_geno)
		{
			fprintf(stderr, "Sanity check failed: Genocoding of variants are not the same @ %s (%d/%d).\n", regions_list->at(i_f), max_coded_geno, cur_max_coded_geno);
			exit(1);
		}

		vector<t_annot_region*>* intersects = intersect_annot_regions(cur_geno_regs, geno_regs, false);
		fprintf(stderr, "Pooling %s (%d samples) using %d/%d intersects\n", regions_list->at(i_f), (int)cur_samples_list->size(),
			(int)cur_geno_regs->size(), (int)geno_regs->size());

		for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
		{
			t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* cur_geno_reg = cur_int_info->src_reg;
			t_annot_region* cur_pooled_reg = cur_int_info->dest_reg;

			void** pooled_reg_info = (void**)(cur_pooled_reg->data);
			char* pooled_geno = (char*)(pooled_reg_info[0]);

			void** cur_reg_info = (void**)(cur_geno_reg->data);
			char* cur_geno = (char*)(cur_reg_info[0]);

			// Copy the genotypes.
			int cur_starting_i_s = (int)samples_so_far->size();
			for (int i_s = cur_starting_i_s; 
				i_s < cur_starting_i_s + (int)cur_samples_list->size();
				i_s++)
			{
				pooled_geno[i_s] = cur_geno[i_s - cur_starting_i_s];
			} // i_s loop.

			delete cur_int_info;
		} // i_int loop.
		delete_annot_regions(intersects);

		// Dont free memory.
		//// Delete memory for the current matrix.
		//for (int i_reg = 0; i_reg < (int)cur_geno_regs->size(); i_reg++)
		//{
		//	void** cur_reg_info = (void**)(cur_geno_regs->at(i_reg)->data);
		//	char* cur_reg_geno = (char*)(cur_reg_info[0]);
		//	delete[] cur_reg_geno;
		//	delete[] cur_reg_info;
		//} // i_reg loop.
		//delete_annot_regions(cur_geno_regs);

		samples_so_far->insert(samples_so_far->end(), cur_samples_list->begin(), cur_samples_list->end());
	} // i_f loop.

	fprintf(stderr, "Saving %d samples (%d samples).\n", (int)samples_so_far->size(), (int)pooled_sample_ids->size());
	//binarize_variant_genotype_signal_regions(geno_regs, NULL, pooled_sample_ids, regions_op_fp);
	binarize_variant_signal_regions_wrapper(geno_regs, pooled_sample_ids, regions_op_fp);

	// Save sample id's.
	FILE* f_samples_list = open_f(samples_op_fp, "w");
	for (int i_s = 0; i_s < (int)pooled_sample_ids->size(); i_s++)
	{
		fprintf(f_samples_list, "%s\n", pooled_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_samples_list);
}

void merge_samples_per_matching_genotype_signal_regions(char* sample1_regs_fp, char* sample1_list_fp, 
														char* sample2_regs_fp, char* sample2_list_fp, bool match_region_names, char* op_fp)
{
	vector<char*>* sample1_ids = buffer_file(sample1_list_fp);
	vector<t_annot_region*>* geno1_sig_regs = load_variant_signal_regions_wrapper(sample1_regs_fp, sample1_list_fp);

	int max_coded_geno1 = get_max_genotype_value(geno1_sig_regs, sample1_ids);

	vector<char*>* sample2_ids = buffer_file(sample2_list_fp);
	vector<t_annot_region*>* geno2_sig_regs = load_variant_signal_regions_wrapper(sample2_regs_fp, sample2_list_fp);

	int max_coded_geno2 = get_max_genotype_value(geno2_sig_regs, sample2_ids);

	if (max_coded_geno1 != max_coded_geno2)
	{
		fprintf(stderr, "Genocoding of variants seem to be different (%d/%d), make sure they are the same.\n", max_coded_geno1, max_coded_geno2);
		exit(1);
	}

	fprintf(stderr, "Merging %d regions of %d samples in %s and %d regions of %d samples and saving to %s; %s\n",
		(int)geno1_sig_regs->size(), (int)sample1_ids->size(), sample1_regs_fp,
		(int)geno2_sig_regs->size(), (int)sample2_ids->size(), sample2_regs_fp,
		op_fp);

	vector<t_annot_region*>* intersects = intersect_annot_regions(geno1_sig_regs, geno2_sig_regs, true);
	fprintf(stderr, "Processing %d matching regions.\n", (int)intersects->size());

	vector<t_annot_region*>* merged_sample_genotype_regions = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* geno1_reg = int_info->src_reg;
		t_annot_region* geno2_reg = int_info->dest_reg;

		if (geno1_reg->start == geno2_reg->start &&
			geno1_reg->end == geno2_reg->end)
		{
			if (!match_region_names ||
				t_string::compare_strings(geno1_reg->name, geno2_reg->name))
			{
				t_annot_region* new_reg = get_empty_region();
				new_reg->chrom = t_string::copy_me_str(geno1_reg->chrom);
				new_reg->start = geno1_reg->start;
				new_reg->end = geno1_reg->end;
				new_reg->strand = geno1_reg->strand;

				if (match_region_names)
				{
					new_reg->name = t_string::copy_me_str(geno1_reg->name);
				}

				void** reg1_info = (void**)(geno1_reg->data);
				char* reg1_geno_sig = (char*)(reg1_info[0]);
				void** reg2_info = (void**)(geno2_reg->data);
				char* reg2_geno_sig = (char*)(reg2_info[0]);

				char* pooled_geno_sig = new char[(int)sample1_ids->size() + (int)sample2_ids->size() + 2];
				for (int i_s = 0; i_s < (int)sample1_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s] = reg1_geno_sig[i_s];
				}

				for (int i_s = 0; i_s < (int)sample2_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s + (int)sample1_ids->size()] = reg2_geno_sig[i_s];
				}

				// Copy the region info.
				void** new_reg_info = new void*[2];
				new_reg_info[0] = pooled_geno_sig;
				new_reg->data = new_reg_info;

				merged_sample_genotype_regions->push_back(new_reg);
			} // region name check.
		} // coordinate check.
	} // i_int loop.

	vector<char*>* merged_sample_ids = new vector<char*>();
	merged_sample_ids->insert(merged_sample_ids->end(), sample1_ids->begin(), sample1_ids->end());
	merged_sample_ids->insert(merged_sample_ids->end(), sample2_ids->begin(), sample2_ids->end());

	fprintf(stderr, "Saving the %d regions with %d samples to %s\n", (int)merged_sample_genotype_regions->size(), (int)merged_sample_ids->size(), op_fp);

	// Write the sample ids.
	char merged_sample_ids_fp[1000];
	sprintf(merged_sample_ids_fp, "%s_samples.list", op_fp);
	FILE* f_merged_sample_ids = open_f(merged_sample_ids_fp, "w");
	for (int i_s = 0; i_s < (int)merged_sample_ids->size(); i_s++)
	{
		fprintf(f_merged_sample_ids, "%s\n", merged_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_merged_sample_ids);

	binarize_variant_genotype_signal_regions(merged_sample_genotype_regions, NULL, merged_sample_ids, op_fp);
}

void extract_genotype_signals_per_chr_ids(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* chr_ids_list_fp, char* op_fp)
{
	//vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	vector<char*>* all_chr_ids_2_use = buffer_file(chr_ids_list_fp);
	vector<char*>* chr_ids_2_use = new vector<char*>();
	for (int i_chr = 0; i_chr < (int)all_chr_ids_2_use->size(); i_chr++)
	{
		char* cur_chr_id = t_string::copy_me_str(all_chr_ids_2_use->at(i_chr));
		normalize_chr_id(cur_chr_id);
		chr_ids_2_use->push_back(cur_chr_id);
		fprintf(stderr, "%s (%s)\n", cur_chr_id, all_chr_ids_2_use->at(i_chr));
	} // i_reg loop.
	fprintf(stderr, "Extracting %d chromosomes.\n", (int)chr_ids_2_use->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	//int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		if(t_string::get_i_str(chr_ids_2_use, geno_sig_regs->at(i_reg)->chrom) < (int)chr_ids_2_use->size())
		{
			roi_regs_w_signals->push_back(geno_sig_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", (int)roi_regs_w_signals->size());

	// Save the regions with signals on them.
	binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
} // extract_genotype_signals_per_chr_ids

void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	for (int i_reg = 0; i_reg < (int)roi_regs->size(); i_reg++)
	{
		roi_regs->at(i_reg)->data = NULL;
	} // i_reg loop.
	fprintf(stderr, "Extracting genotype signals for %d ROIs.\n", (int)roi_regs->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	fprintf(stderr, "Intersecting ROI's with genotype signal regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(roi_regs, geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", (int)intersects->size());
	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_roi_reg = int_info->src_reg;
		t_annot_region* cur_geno_sig_reg = int_info->dest_reg;

		if (t_string::compare_strings(cur_roi_reg->chrom, cur_geno_sig_reg->chrom) &&
			cur_roi_reg->start == cur_geno_sig_reg->start &&
			cur_roi_reg->end == cur_geno_sig_reg->end)
		{
			if (cur_roi_reg->data != NULL)
			{
				fprintf(stderr, "***WARNING: Already assigned: %s:%d-%d ***\n", cur_geno_sig_reg->chrom, cur_geno_sig_reg->start, cur_geno_sig_reg->end);
				//exit(1);
			}

			cur_roi_reg->data = cur_geno_sig_reg->data;
		}
	} // i_int loop.

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	//int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < (int)roi_regs->size(); i_reg++)
	{
		if (roi_regs->at(i_reg)->data != NULL)
		{
			roi_regs_w_signals->push_back(roi_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", (int)roi_regs_w_signals->size());

	// Save the regions with signals on them.
	//binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(roi_regs_w_signals, sample_ids, op_fp);
}

// These are the defacto standards for mapping genotypes to vcf entries. These should be used in any function that converts from genotype to vcf string.
void geno2vcf_entry(const char geno, char* vcf_geno_str, bool haplocoded)
{
	if (haplocoded)
	{
		// Entries for 0,1,2,3.
		char per_geno_phased_vcf_geno[4][4] = { "0|0", "0|1", "1|0", "1|1" };
		strcpy(vcf_geno_str, per_geno_phased_vcf_geno[(unsigned char)geno]);
	}
	else
	{
		// Entries for 0,1,2.
		char per_geno_unphased_vcf_geno[3][4] = { "0/0", "0/1", "1/1" };
		strcpy(vcf_geno_str, per_geno_unphased_vcf_geno[(unsigned char)geno]);
	}
}

// Matches geno2vcf_entry(...).
// This function assumes that vcf_entry is a valid genotype string. It does not do checks.
char vcf_entry2geno(const char* vcf_entry, bool haplocoded)
{
	if (vcf_entry[0] == '.' || vcf_entry[2] == '.')
	{
		return(-1);
	}

	char geno0 = vcf_entry[0] - '0';
	char geno1 = vcf_entry[2] - '0';

	if (haplocoded)
	{
		return(2 * geno0 + geno1);
	}
	else
	{
		return(geno0 + geno1);
	}
}

char alleles2geno(int all0, int all1, bool haplocoded)
{
	char geno = all0 + all1;
	if (haplocoded)
	{
		geno = all0 * 2 + all1;
	}

	return(geno);
}

char alleles2geno(int* all01, bool haplocoded)
{
	return(alleles2geno(all01[0], all01[1], haplocoded));
}

void alleles_2_vcf_entry(int all0, int all1, char* vcf_entry, bool haplocoded)
{
	char geno = alleles2geno(all0, all1, haplocoded);

	geno2vcf_entry(geno, vcf_entry, haplocoded);
}


// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF_memoptimized(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* var_regions_BED_fp,			// Regions to focus on while extracting.
	char* chr_id_2_process,				// Chromosome to process.
	char* bin_seq_dir,					// Sequences directory.
	bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
	bool match_region_names_flag,		// This is a flag to enforce matching of region names.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(1);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", (int)vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == (int)chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", (int)var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[(int)chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			int l_chrom = 0;
			fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
			char cur_chr_seq_fp[1000];
			sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

			if (check_file(cur_chr_seq_fp))
			{
				per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
			}
			else
			{
				sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

				if (check_file(cur_chr_seq_fp))
				{
					per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
				}
				else
				{
					fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
					exit(1);
				}
			}
		} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	// Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < (int)var_regions->size(); i_v++)
	{
		//char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	// Set sorting info per variant regions on each chromosome.
	t_restr_annot_region_list* restr_var_regions = restructure_annot_regions(var_regions);
	for (int i_chr = 0; i_chr < (int)restr_var_regions->chr_ids->size(); i_chr++)
	{
		sort_set_sorting_info(restr_var_regions->regions_per_chrom[i_chr], sort_regions);
	} // i_chr loop.

	char* buff = new char[100 * 1000];
	int n_processed_regions = 0;
	int n_assigned_regions = 0;
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if (n_processed_regions % 10000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region (%d).           \r", n_processed_regions, n_assigned_regions);
		}

		n_processed_regions++;

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
		char chrom[100];
		char posn_str[100];
		char id[100];
		char ref[100];
		char alt[100];
		char qual[100];
		char filter[100];
		char info[10000];
		char format[100];

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(format, buff);

		// Check for intersect with variant regions.
		int var_i_chr = t_string::get_i_str(restr_var_regions->chr_ids, chrom);
		if (var_i_chr == (int)restr_var_regions->chr_ids->size())
		{
			delete[] cur_line;
			continue;
		}

		int cur_var_posn = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
		int cur_var_end = cur_var_posn + t_string::string_length(ref) - 1;
		vector<t_annot_region*>* var_regs_per_cur_chr = restr_var_regions->regions_per_chrom[var_i_chr];
		int i_leftmost_reg = locate_posn_region_per_region_starts(cur_var_posn, var_regs_per_cur_chr, 0, (int)var_regs_per_cur_chr->size() - 1);

		// Go back till the cumulative end for the reg2 is to the left of reg1.
		while (i_leftmost_reg > 0 &&
				var_regs_per_cur_chr->at(i_leftmost_reg)->sort_info->cumulative_sorted_end >= cur_var_posn)
		{
			i_leftmost_reg--;
		} // i_leftmost_reg loop.

		// Check if there is any overlap. Check while we are not strictly passed over the region.
		bool all_checks_pass = false;
		int ol_start = 0;
		int ol_end = 0;
		int matching_i_reg = -1;
		while (i_leftmost_reg < (int)var_regs_per_cur_chr->size() &&
				var_regs_per_cur_chr->at(i_leftmost_reg)->start <= cur_var_end)
		{
			ol_start = MAX(var_regs_per_cur_chr->at(i_leftmost_reg)->start, cur_var_posn);
			ol_end = MIN(var_regs_per_cur_chr->at(i_leftmost_reg)->end, cur_var_posn);

			if (ol_end >= ol_start)
			{
				bool coord_check_pass = true;
				if (cur_var_posn != var_regs_per_cur_chr->at(i_leftmost_reg)->start ||
					cur_var_posn != var_regs_per_cur_chr->at(i_leftmost_reg)->end)
				{
					coord_check_pass = false;
				}

				// Do we need to have name matching?
				bool name_check_pass = true;
				if (match_region_names_flag &&
					var_regs_per_cur_chr->at(i_leftmost_reg)->name != NULL &&
					id != NULL)
				{
					int i_char = 0;

					// If there is a dot for the name, skip.
					if (!t_string::compare_strings(var_regs_per_cur_chr->at(i_leftmost_reg)->name, ".") &&
						!t_string::compare_strings(id, ".") &&
						!t_string::compare_substrings_ci(var_regs_per_cur_chr->at(i_leftmost_reg)->name, id, i_char) &&
						!t_string::compare_substrings_ci(id, var_regs_per_cur_chr->at(i_leftmost_reg)->name, i_char))
					{
						name_check_pass = false;
					}
				}

				if (name_check_pass &&
					coord_check_pass)
				{
					all_checks_pass = true;
					matching_i_reg = i_leftmost_reg;
					break;
				}
			} // overlap check.

			i_leftmost_reg++;
		} // left most position tracing loop.

		if (!all_checks_pass)
		{
			delete[] cur_line;
			continue;
		}

		if (matching_i_reg == -1)
		{
			fprintf(stderr, "Matching region index is invalid.\n");
			exit(1);
		}

		// Set the vcf region.
		t_annot_region* cur_vcf_reg = var_regs_per_cur_chr->at(matching_i_reg);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];

		// It is important to provide the actual length of buffer not the string length.
		while (t_string::get_next_token(format, cur_tok, l_format + 2, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.
		delete[] cur_tok;

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}

		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(1);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the genotype correctly for entry %d in line: %s:%s: %s. sample: %s (Potentially multiallelic)\n", i_s, chrom, posn_str, id, buff);
				//exit(1);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
					//char cur_geno = vcf_entry2geno(buff, true);
					//cur_var_geno_sig[i_s] = cur_geno;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;
					//char cur_geno = vcf_entry2geno(buff, true);

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;
					//char cur_geno = vcf_entry2geno(buff, true);

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(1);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Update assigned regions count.
			n_assigned_regions++;

			// Intersect with the variant regions.
			int l_ref_allele = t_string::string_length(ref);

			// Don't change the name of the target region.
			//cur_vcf_reg->name = t_string::copy_me_str(id);
			
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// Set the region info.
			cur_vcf_reg->data = vcf_reg_info;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(1);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	fprintf(stderr, "Assigned to %d regions in total out of %d VCF regions.\n", n_assigned_regions, n_processed_regions);

	// Copy the current haplotype alleles: Make sure every region has some signal.
	for (int i_reg = 0; i_reg < (int)var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size()];
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	//binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(var_regions, vcf_sample_ids, op_fp);
}


// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
										char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
										char* var_regions_BED_fp,			// Regions to focus on while extracting.
										char* chr_id_2_process,				// Chromosome to process.
										char* bin_seq_dir,					// Sequences directory.
										bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
										bool match_region_names_flag,		// This is a flag to enforce matching of region names.
										bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
										char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(1);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", (int)vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == (int)chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", (int)var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[(int)chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			int l_chrom = 0;
			fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
			char cur_chr_seq_fp[1000];
			sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

			if (check_file(cur_chr_seq_fp))
			{
				per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
			}
			else
			{
				sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

				if (check_file(cur_chr_seq_fp))
				{
					per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
				}
				else
				{
					fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
					exit(1);
				}
			}
		} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	  // Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < (int)var_regions->size(); i_v++)
	{
		//char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	int L_BUFF = 100 * 1000;
	char* buff = new char[L_BUFF];
	char* chrom = new char[L_BUFF];
	char* posn_str = new char[L_BUFF];
	char* id = new char[L_BUFF];
	char* ref = new char[L_BUFF];
	char* alt = new char[L_BUFF];
	char* qual = new char[L_BUFF];
	char* filter = new char[L_BUFF];
	char* info = new char[L_BUFF];
	char* format = new char[L_BUFF];

	vector<t_annot_region*>* vcf_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if ((int)vcf_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region.           \r", (int)vcf_regs->size());
		}

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);
		strcpy(format, buff);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];

		// It is important to provide the actual length of buffer not the string length.
		while (t_string::get_next_token(format, cur_tok, l_format+2, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.
		delete[] cur_tok;

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}
		
		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, L_BUFF, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(1);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the %d. genotype correctly for entry %s:%s: %s. sample: %s (Potentially multiallelic)\n", i_s, chrom, posn_str, id, buff);
				//exit(1);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(1);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Intersect with the variant regions.
			t_annot_region* cur_vcf_reg = get_empty_region();
			cur_vcf_reg->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(cur_vcf_reg->chrom);
			cur_vcf_reg->start = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);

			int l_ref_allele = t_string::string_length(ref);
			cur_vcf_reg->end = cur_vcf_reg->start + l_ref_allele - 1;
			cur_vcf_reg->name = t_string::copy_me_str(id);
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(1);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.

			cur_vcf_reg->data = vcf_reg_info;

			vcf_regs->push_back(cur_vcf_reg);
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regions, vcf_regs, true);
	fprintf(stderr, "\nProcessing %d intersections using coordinate matching.               \n", (int)intersects->size());

	if (match_region_names_flag)
	{
		fprintf(stderr, "Enforcing substring based name matching.\n");
	}

	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* pooled_var_reg = cur_int_info->src_reg;
		t_annot_region* vcf_reg = cur_int_info->dest_reg;

		bool coord_check_pass = true;
		if (pooled_var_reg->start != vcf_reg->start ||
			pooled_var_reg->end != vcf_reg->end)
		{
			coord_check_pass = false;
		}

		// Do we need to have name matching?
		bool name_check_pass = true;
		if (match_region_names_flag &&
			pooled_var_reg->name != NULL &&
			vcf_reg->name != NULL)
		{
			int i_char = 0;

			// If there is a dot for the name, skip.
			if (!t_string::compare_strings(pooled_var_reg->name, ".") &&
				!t_string::compare_strings(vcf_reg->name, ".") &&
				!t_string::compare_substrings_ci(pooled_var_reg->name, vcf_reg->name, i_char) &&
				!t_string::compare_substrings_ci(vcf_reg->name, pooled_var_reg->name, i_char))
			{
				name_check_pass = false;
			}
		}

		if (name_check_pass &&
			coord_check_pass)
		{
			//int i_chr = t_string::get_i_str(chr_ids, pooled_var_reg->chrom);
			//char* chr_seq = per_chr_seq[i_chr];

			// Get the genotype signal for the vcf region.
			void** vcf_reg_info = (void**)(vcf_reg->data);
			char* cur_vcf_reg_geno_sig = (char*)(vcf_reg_info[0]);

			// Set the genotype signal.
			void** cur_var_reg_info = (void**)(pooled_var_reg->data);
			cur_var_reg_info[0] = cur_vcf_reg_geno_sig;
		}
		else
		{

		}

		delete cur_int_info;
	} // i_int loop.
	delete_annot_regions(vcf_regs);

	// Copy the current haplotype alleles.
	for (int i_reg = 0; i_reg < (int)var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size()];
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	//binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(var_regions, vcf_sample_ids, op_fp);

	//// Dump the plain matrix, too.
	//dump_geno_sig_regs_plain(var_regions, vcf_sample_ids, "op.bed");
	
	//// Now dump the variant signal profiles.
	//FILE* f_geno_sig = open_f("op.bed", "w");
	//for (int i_v = 0; i_v < var_regions->size(); i_v++)
	//{
	//	fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
	//		var_regions->at(i_v)->chrom, 
	//		translate_coord(var_regions->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//		translate_coord(var_regions->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		var_regions->at(i_v)->name);

	//	void** cur_var_reg_info = (void**)(var_regions->at(i_v)->data);
	//	int* cur_var_signal = (int*)(cur_var_reg_info[0]);
	//	for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_geno_sig, "\t%d", cur_var_signal[i_s]);
	//	} // i_s loop.

	//	fprintf(f_geno_sig, "\n");
	//} // i_v loop.
	//fclose(f_geno_sig);
}

void convert_genocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from genocoded %s genotypes.\n", geno_sigs_reg_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);

	int max_geno = get_max_genotype_value(geno_sig_regs, sample_ids);
	if (max_geno != 2)
	{
		fprintf(stderr, "##### WARNING : Could not detect genocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect genocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect genocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect genocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
	}

	fprintf(stderr, "Loaded %d regions.\n", (int)geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				//fprintf(stderr, "Converting %d/%d variant.             \r", i_reg, (int)cur_chr_input_geno_sig_regs->size());
				t_string::print_padded_string(stderr, '\r', 100, "Converting %d/%d variant.", i_reg, (int)cur_chr_input_geno_sig_regs->size());
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if ((int)toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(1);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_genocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				int geno = (int)(per_sample_genocoded_geno[i_s]);

				if (geno == 0)
				{
					fprintf(f_VCF, "\t0/0");
				}
				else if (geno == 1)
				{
					fprintf(f_VCF, "\t0/1");
				}
				else if (geno == 2)
				{
					fprintf(f_VCF, "\t1/1");
				}
				else
				{
					fprintf(f_VCF, "\t./.");
					//fprintf(stderr, "Could not parse the genotype: %s:%d::%d\n",
					//	cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
					//	(int)(per_sample_haplocoded_geno[i_s]));

					//exit(1);
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}

static void* thread_callback_convert_haplocoded_signal_regs_2_VCF_w_header(void* __thread_info_ptr)
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

	//fprintf(stderr, "VCF Export Thread-%d: %d-%d // %d [max_geno=%d]\n", thr_i, var_start_i, var_end_i, vecsize(cur_chr_input_geno_sig_regs), max_input_geno);
	t_string::print_padded_string(stderr, '\r', 100, "VCF Export Thread-%d: %d-%d // %d [max_geno=%d]", thr_i, var_start_i, var_end_i, vecsize(cur_chr_input_geno_sig_regs), max_input_geno);

	char cur_thread_file_name[1000];
	sprintf(cur_thread_file_name, "temp_vcf_export_%d.vcf.gz", thr_i);
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

void convert_haplocoded_signal_regs_2_VCF_w_header_multithreaded(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* header_file, bool phase_op, int n_threads, char* op_VCF_fp)
{
	fprintf(stderr, "%d-thread Exporting the VCF file using the panel data @ %s (%s) on %s using header in %s.\n",
		n_threads,
		geno_sigs_reg_fp, sample_ids_list_fp,
		op_VCF_fp,
		header_file);

	if (phase_op)
	{
		fprintf(stderr, "**Saving phased VCF.**\n");
	}

	vector<char*>* reference_haplo_sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", (int)reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	vector<char*>* header_lines = buffer_file(header_file);
	fprintf(stderr, "Loaded %d header lines.\n", vecsize(header_lines));

	int max_ref_geno = get_max_genotype_value(reference_haplo_geno_sig_regs, reference_haplo_sample_ids);
	if (max_ref_geno != 3 &&
		phase_op)
	{
		fprintf(stderr, "%s(%d): Reference panel does not seem to be phased, cannot save an unphased panel as a phased VCF.\n", __FILE__, __LINE__);
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

	gzFile f_ref_header = gzopen("TEMP_VCF_HEADER.vcf.gz", "wb");
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
			int_vals[4] = phase_op;
			int_vals[5] = max_ref_geno;
			double* dbl_vals = new double[10];
			cur_thread_info[1] = dbl_vals;
			cur_thread_info[2] = cur_chr_ref_panel_var_regs;
			cur_thread_info[3] = reference_haplo_sample_ids;

			t_ansi_thread* thread = new t_ansi_thread(thread_callback_convert_haplocoded_signal_regs_2_VCF_w_header, cur_thread_info);
			threads->push_back(thread);
			thread->run_thread();

			cur_var_start_i += n_vars_per_thread;

		} // i_thr loop.

		vector<char*>* per_thread_ref_files = new vector<char*>();
		per_thread_ref_files->push_back(t_string::copy_me_str("TEMP_VCF_HEADER.vcf.gz"));
		for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
		{
			threads->at(i_thr)->wait_thread();
			t_string::print_padded_string(stderr, '\r', 100, "Ref-Thread-%d finished..", i_thr);

			char cur_thread_file_name[1000];
			sprintf(cur_thread_file_name, "temp_vcf_export_%d.vcf.gz", i_thr);

			per_thread_ref_files->push_back(t_string::copy_me_str(cur_thread_file_name));
		} // i_thr loop.

		fprintf(stderr, "Concatenating %d Ref matrices..                       \n", vecsize(per_thread_ref_files));
		concatenateGzipFiles(op_VCF_fp, per_thread_ref_files);

		fprintf(stderr, "Cleaning up..\n");
		for (int i_thr = 0; i_thr < vecsize(per_thread_ref_files); i_thr++)
		{
			// Delete this file to make sure it does not interfere??
			t_string::print_padded_string(stderr, '\r', 100, "Deleting %s", per_thread_ref_files->at(i_thr));
			delete_file(per_thread_ref_files->at(i_thr));
		} // i_thr loop.

		// Do not process multiple chromosome here.
		break;
	} // i_chr loop.
}


void convert_haplocoded_signal_regs_2_VCF_w_header(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* header_file, bool phase_op, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from %s using phased=%d genotypes.\n", geno_sigs_reg_fp, phase_op);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	vector<char*>* header_lines = buffer_file(header_file);
	fprintf(stderr, "Loaded %d header lines..\n", vecsize(header_lines));

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);

	int max_geno = get_max_genotype_value(geno_sig_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
	}

	fprintf(stderr, "Loaded %d regions.\n", (int)geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	//fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	for (int i_h = 0; i_h < vecsize(header_lines); i_h++)
	{
		fprintf(f_VCF, "%s\n", header_lines->at(i_h));
	} // i_h loop.

	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				t_string::print_padded_string(stderr, '\r', 100, "Converting %d/%d variant.", i_reg, (int)cur_chr_input_geno_sig_regs->size());

			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if ((int)toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(1);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				if (!phase_op)
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int counted_geno = cur_hap0 + cur_hap1;
					int haplocoded_geno = per_sample_haplocoded_geno[i_s];

					if (haplocoded_geno > 3 ||
						haplocoded_geno < 0)
					{
						fprintf(f_VCF, "\t./.");
					}
					else if (counted_geno == 0)
					{
						fprintf(f_VCF, "\t0/0");
					}
					else if (counted_geno == 1)
					{
						fprintf(f_VCF, "\t0/1");
					}
					else if (counted_geno == 2)
					{
						fprintf(f_VCF, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(1);
					}
				}
				else
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int haplocoded_geno = per_sample_haplocoded_geno[i_s];

					if (haplocoded_geno > 3 ||
						haplocoded_geno < 0)
					{
						fprintf(f_VCF, "\t.|.");
					}
					else
					{
						char geno_str[10];
						strcpy(geno_str, "0|0");

						if (cur_hap0 == 1)
						{
							geno_str[0] = '1';
						}

						if (cur_hap1 == 1)
						{
							geno_str[2] = '1';
						}

						fprintf(f_VCF, "\t%s", geno_str);

						int geno = cur_hap0 + cur_hap1;
						if (geno > 2)
						{
							fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
								cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
								(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

							exit(1);
						}
					}
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}


void convert_haplocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, bool phase_op, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from %s using phased=%d genotypes.\n", geno_sigs_reg_fp, phase_op);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sigs_reg_fp, sample_ids_list_fp);

	int max_geno = get_max_genotype_value(geno_sig_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
		fprintf(stderr, "##### WARNING : Could not detect haplocoding of the genotype signal, will proceed with converting %s to VCF but make sure this is the intended behavior. ###########\n", geno_sigs_reg_fp);
	}

	fprintf(stderr, "Loaded %d regions.\n", (int)geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < (int)cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				t_string::print_padded_string(stderr, '\r', 100, "Converting %d/%d variant.", i_reg, (int)cur_chr_input_geno_sig_regs->size());
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if ((int)toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(1);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				if (!phase_op)
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int counted_geno = cur_hap0 + cur_hap1;
					int haplocoded_geno = per_sample_haplocoded_geno[i_s];

					if (haplocoded_geno > 3 ||
						haplocoded_geno < 0)
					{
						fprintf(f_VCF, "\t./.");
					}
					else if (counted_geno == 0)
					{
						fprintf(f_VCF, "\t0/0");
					}
					else if (counted_geno == 1)
					{
						fprintf(f_VCF, "\t0/1");
					}
					else if (counted_geno == 2)
					{
						fprintf(f_VCF, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(1);
					}
				}
				else
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int haplocoded_geno = per_sample_haplocoded_geno[i_s];

					if (haplocoded_geno > 3 ||
						haplocoded_geno < 0)
					{
						fprintf(f_VCF, "\t.|.");
					}
					else
					{
						char geno_str[10];
						strcpy(geno_str, "0|0");

						if (cur_hap0 == 1)
						{
							geno_str[0] = '1';
						}

						if (cur_hap1 == 1)
						{
							geno_str[2] = '1';
						}

						fprintf(f_VCF, "\t%s", geno_str);

						int geno = cur_hap0 + cur_hap1;
						if (geno > 2)
						{
							fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
								cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
								(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

							exit(1);
						}
					}
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}

int get_genotype_per_haplocoded_genotype(char haplocoded_geno)
{
	if (haplocoded_geno > 3 || haplocoded_geno < 0)
	{
		return(0);
	}

	int all1 = get_allele_per_haplotype(haplocoded_geno, 0);
	int all2 = get_allele_per_haplotype(haplocoded_geno, 1);
	return(all1 + all2);
}

int get_max_genotype_value(vector<t_annot_region*>* regs, vector<char*>* sample_ids)
{
	int max_geno = 0;
	for (int i_reg = 0; i_reg < (int)regs->size(); i_reg++)
	{
		void** reg_info = (void**)(regs->at(i_reg)->data);
		char* reg_geno = (char*)(reg_info[0]);
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			max_geno = (max_geno < reg_geno[i_s]) ? (reg_geno[i_s]) : (max_geno);

			if (max_geno == 3)
			{
				return(max_geno);
			}
		} // i_s loop.
	} // i_reg loop.

	return(max_geno);
}

bool is_haplocoded_allele_het(char haplocoded_geno)
{
	if (haplocoded_geno == 1 ||
		haplocoded_geno == 2)
	{
		return true;
	}

	return false;
}

int get_allele_per_haplotype(char geno, int hap_i)
{
	// If the genotype is non-existing or not valid, return non-existing.
	if (geno > 3 || geno < 0)
	{
		return 0;
	}

	int allele = ((geno & (1 << hap_i)) >> hap_i);
	return(allele);
}

int get_allele_haploswitch_per_haplotype(char geno, int hap_i, char& haplo_switch)
{
	char switch_allele_2bit = (geno & (3 << (2 * hap_i))) >> (2 * hap_i);

	haplo_switch = (switch_allele_2bit >> 1);
	int allele = (switch_allele_2bit & 1);

	return(allele);
}

void dump_haplo_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, int haplo_2_extract, const char* op_fp)
{
	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < (int)geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);


		void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
		char* cur_var_signal = (char*)(cur_var_reg_info[0]);
		for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
		{
			char cur_sample_haplo_str[10];
			memset(cur_sample_haplo_str, 0, 10);

			if (haplo_2_extract == 2)
			{
				cur_sample_haplo_str[1] = '|';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 3)
			{
				cur_sample_haplo_str[1] = '\t';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 0)
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
			}
			else
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}

			fprintf(f_geno_sig, "\t%s", cur_sample_haplo_str);
		} // i_s loop.

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_geno_sig_regs_plain, const char* op_fp)
{
	if (dump_geno_sig_regs_plain)
	{
		fprintf(stderr, "Saving regions only.\n");
	}

	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < (int)geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);

		if (!dump_geno_sig_regs_plain)
		{
			void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
			char* cur_var_signal = (char*)(cur_var_reg_info[0]);
			for (int i_s = 0; i_s < (int)vcf_sample_ids->size(); i_s++)
			{
				fprintf(f_geno_sig, "\t%d", (int)cur_var_signal[i_s]);
			} // i_s loop.
		}

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	close_f(f_geno_sig, op_fp);
}

void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", (int)geno_sig_regs->size(), (int)sample_ids->size());

	vector<char*>* geno_subsample_ids = buffer_file(subsample_ids_list_fp);
	if (geno_subsample_ids == NULL)
	{
		fprintf(stderr, "Could not load subsample id's list from %s\n", subsample_ids_list_fp);
		exit(1);
	}
	else
	{
		fprintf(stderr, "Extracting genotype signals for %d individuals.\n", (int)geno_subsample_ids->size());
	}

	vector<int>* per_subsample_geno_sample_i = new vector<int>();
	int n_unmatched_samples = 0;
	for (int i_ss = 0; i_ss < (int)geno_subsample_ids->size(); i_ss++)
	{
		int cur_ss_sample_i = t_string::get_i_str(sample_ids, geno_subsample_ids->at(i_ss));
		per_subsample_geno_sample_i->push_back(cur_ss_sample_i);

		if (cur_ss_sample_i == (int)sample_ids->size())
		{
			n_unmatched_samples++;
		}
	} // i_sub_s loop.

	fprintf(stderr, "Found %d unmatched samples.\n", n_unmatched_samples);

	vector<t_annot_region*>* geno_sig_regs_w_subsample_signals = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			t_string::print_padded_string(stderr, '\r', 100, "Extracting %d. region", i_reg);
		}

		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		char* cur_reg_geno_sigs = (char*)(cur_reg_info[0]);

		// Copy the signals.
		char* cur_copy_reg_geno_sigs = new char[(int)geno_subsample_ids->size() + 2];
		for (int i_ss = 0; i_ss < (int)geno_subsample_ids->size(); i_ss++)
		{			
			if (per_subsample_geno_sample_i->at(i_ss) < (int)sample_ids->size())
			{
				cur_copy_reg_geno_sigs[i_ss] = cur_reg_geno_sigs[per_subsample_geno_sample_i->at(i_ss)];
			}
			else
			{
				cur_copy_reg_geno_sigs[i_ss] = -1;
			}
		} // i_ss loop.

		void** cur_copy_reg_info = new void*[2];		
		cur_copy_reg_info[0] = cur_copy_reg_geno_sigs;

		t_annot_region* copy_reg = duplicate_region(geno_sig_regs->at(i_reg));
		copy_reg->data = cur_copy_reg_info;
		geno_sig_regs_w_subsample_signals->push_back(copy_reg);
	} // i_reg loop.

	// Save the regions.
	//binarize_variant_genotype_signal_regions(geno_sig_regs_w_subsample_signals, NULL, geno_subsample_ids, op_fp);
	binarize_variant_signal_regions_wrapper(geno_sig_regs_w_subsample_signals, geno_subsample_ids, op_fp);
} // extract_genotype_signals_per_subsample_list function.

void filter_redundant_variants_per_genotype_distance(char* genotype_matbed_fp, char* sample_ids_list_fp,
	double max_geno_distance, int max_index_difference,
	char* unique_regs_BED_op_fp)
{
	vector<t_annot_region*>* regs = load_variant_signal_regions_wrapper(genotype_matbed_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	t_restr_annot_region_list* restr_regs = restructure_annot_regions(regs);

	fprintf(stderr, "Calculating genotype unique variants out of %d variants for %d subjects using %d index long local windows with maximum genotype distance of %.0f..\n", 
			vecsize(regs), vecsize(sample_ids),
			max_index_difference, max_geno_distance);

	int max_geno = get_max_genotype_value(regs, sample_ids);
	if (max_geno == 3)
	{
		fprintf(stderr, "Detected haplocoded genotypes..\n");
	}
	else if (max_geno == 2)
	{
		fprintf(stderr, "Detected genocoded genotypes..\n");
	}
	else
	{
		fprintf(stderr, "Could not detect the genotype coding (%d)\n", max_geno);
		exit(1);
	}

	char non_unique_vars_BED_fp[1000];
	sprintf(non_unique_vars_BED_fp, "%s_nonuniqueness.bed", genotype_matbed_fp);
	FILE* f_non_unique_vars_BED = open_f(non_unique_vars_BED_fp, "w");

	vector<t_annot_region*>* unique_geno_signal_regs = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < vecsize(restr_regs->chr_ids); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_regs = restr_regs->regions_per_chrom[i_chr];

		// Assign AC's:
		for (int i_reg = 0; i_reg < vecsize(cur_chr_regs); i_reg++)
		{
			cur_chr_regs->at(i_reg)->dbl_score = 0;

			void** cur_reg_info = (void**)(cur_chr_regs->at(i_reg)->data);
			char* cur_reg_geno = (char*)(cur_reg_info[0]);
			for (int i_s = 0; i_s < vecsize(sample_ids); i_s++)
			{
				int cur_geno = (int)(cur_reg_geno[i_s]);
				if (max_geno == 3)
				{
					cur_geno = get_genotype_per_haplocoded_genotype(cur_geno);
				}

				cur_chr_regs->at(i_reg)->dbl_score += cur_geno;
			} // i_s loop.
		} // i_reg loop.

		// Start comparing variants.
		for (int i_reg = 0; i_reg < vecsize(cur_chr_regs); i_reg++)
		{
			void** cur_i_reg_info = (void**)(cur_chr_regs->at(i_reg)->data);
			char* cur_i_reg_geno = (char*)(cur_i_reg_info[0]);

			// Note that we are always comparing with previous regions; i.e., 
			int min_j_reg = MAX((i_reg - max_index_difference), 0);
			int max_j_reg = i_reg - 1;

			bool i_reg_unique = true;

			for (int j_reg = min_j_reg; j_reg < max_j_reg; j_reg++)
			{
				// AC count check first. If i and j do not have the same AC count, i cannot match j and we dont have to compare.
				if (max_geno_distance > 0 ||
					(cur_chr_regs->at(i_reg)->dbl_score == cur_chr_regs->at(j_reg)->dbl_score))
				{
					// Compare
					void** cur_j_reg_info = (void**)(cur_chr_regs->at(j_reg)->data);
					char* cur_j_reg_geno = (char*)(cur_j_reg_info[0]);

					int cur_ij_dist = 0;
					for (int i_s = 0; i_s < vecsize(sample_ids); i_s++)
					{
						if (cur_i_reg_geno[i_s] != cur_j_reg_geno[i_s])
						{
							cur_ij_dist++;

							// We don't have to run all subjects over the max difference.
							if (cur_ij_dist > max_geno_distance)
							{
								break;
							}

						}
					} // i_s loop.

					// If j matches i exactly, this variant is not unique.
					if (cur_ij_dist <= max_geno_distance)
					{
						fprintf(f_non_unique_vars_BED, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n",
								cur_chr_regs->at(i_reg)->chrom, 
								translate_coord(cur_chr_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
								translate_coord(cur_chr_regs->at(i_reg)->start, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
								cur_chr_regs->at(j_reg)->chrom,
								translate_coord(cur_chr_regs->at(j_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
								translate_coord(cur_chr_regs->at(j_reg)->start, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
								cur_ij_dist,
								i_reg, j_reg);

						i_reg_unique = false;
						break;
					}
				} // AC check.
			} // j_reg loop.

			if (i_reg_unique)
			{
				unique_geno_signal_regs->push_back(cur_chr_regs->at(i_reg));
			}
		} // i_reg loop.
	} // i_chr loop

	close_f(f_non_unique_vars_BED, NULL);

	dump_BED(unique_regs_BED_op_fp, unique_geno_signal_regs);
} // filter_redundant_variants_per_genotype_distance option.


vector<t_annot_region*>* load_variant_genotype_signal_regions_per_matrix_subjects_regions(const char* op_prefix, vector<char*>* geno_sample_ids)
{
	if (geno_sample_ids == NULL ||
		vecsize(geno_sample_ids) != 0)
	{
		fprintf(stderr, "Genotype sample id's is not as expected or it does not exist, yet.\n");
		exit(1);
	}

	char subjects_ids_list_fp[1000];
	sprintf(subjects_ids_list_fp, "%s_subjects.list", op_prefix);
	vector<char*>* read_geno_sample_ids = buffer_file(subjects_ids_list_fp);
	//fprintf(stderr, "Loading %d. subject genotypes..\n", vecsize(read_geno_sample_ids));
	geno_sample_ids->insert(geno_sample_ids->end(), read_geno_sample_ids->begin(), read_geno_sample_ids->end());


	char regs_BED_fp[1000];
	sprintf(regs_BED_fp, "%s_variants.bed", op_prefix);
	vector<t_annot_region*>* genotype_signal_regions = load_BED(regs_BED_fp);
	//fprintf(stderr, "Loading genotypes for %d variants..\n", vecsize(genotype_signal_regions));

	fprintf(stderr, "Loading %d subj-vs-%d variant genotype matrix..\n", vecsize(read_geno_sample_ids), vecsize(genotype_signal_regions));

	FILE* f_geno_mat = NULL;
	char geno_matrix_op_fp[1000];
	sprintf(geno_matrix_op_fp, "%s_genotypes.matrix.gz", op_prefix);

	// This the the main buffer that we fill below.
	size_t total_n_bytes = (size_t)(vecsize(geno_sample_ids)) * (size_t)(vecsize(genotype_signal_regions));
	char* genotype_matrix_buff = new char[total_n_bytes + 2];
	if(check_file(geno_matrix_op_fp))
	{
		// Compressed load via zlib.
		gzFile f_geno_mat = gzopen(geno_matrix_op_fp, "rb");
		if (!f_geno_mat)
		{
			perror("gzopen");
			exit(1);
		}

		size_t CHUNKSIZE = (size_t)1024 * (size_t)1024 * (size_t)1024;
		size_t ptr_pos = 0;
		while (1)
		{
			size_t n_bytes_2_read = MIN(CHUNKSIZE, total_n_bytes - ptr_pos);
			if (n_bytes_2_read == 0)
			{
				break;
			}

			t_string::print_padded_string(stderr, '\r', 100, "Loading %.2f Gb.. [%.2f Gb left]", n_bytes_2_read / (double)CHUNKSIZE, (total_n_bytes - ptr_pos) / (double)CHUNKSIZE);
			if (gzread(f_geno_mat, genotype_matrix_buff + ptr_pos, n_bytes_2_read) != (int)n_bytes_2_read)
			{
				fprintf(stderr, "%s(%d): gzread error\n", __FILE__, __LINE__);
				fprintf(stderr, "Could not read genotypes from %s\n", geno_matrix_op_fp);
				gzclose(f_geno_mat);
				exit(1);
			}

			ptr_pos += n_bytes_2_read;
		} // chunk writing loop.

		//// 1 gig chunk size.
		//size_t ptr_pos = 0;
		//size_t n_bytes_2_read = vecsize_t(geno_sample_ids);
		//for(int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
		//{			
		//	//fprintf(stderr, "Writing %.3f megabytes (%.3f megabytes left).             \r", (double)n_bytes_2_write / (1024.0 * 1024.0), (double)(geno_sig_mem_size - ptr_pos) / (1024.0 * 1024.0));
		//	if (gzread(f_geno_mat, genotype_matrix_buff + ptr_pos, n_bytes_2_read) != (int)n_bytes_2_read)
		//	{
		//		fprintf(stderr, "%s(%d): gzread error\n", __FILE__, __LINE__);
		//		fprintf(stderr, "Could not read genotypes from %s\n", geno_matrix_op_fp);
		//		gzclose(f_geno_mat);
		//		exit(1);
		//	}

		//	ptr_pos += n_bytes_2_read;
		//} // chunk writing loop.
		fprintf(stderr, "Completed loading..                       \n");
		gzclose(f_geno_mat);
	}
	else
	{
		sprintf(geno_matrix_op_fp, "%s_genotypes.matrix", op_prefix);

		if (check_file(geno_matrix_op_fp))
		{
			f_geno_mat = open_f(geno_matrix_op_fp, "rb");
		}
		else
		{
			fprintf(stderr, "Could not resolve the genotype matrix file name for %s_genotypes.matrix[.gz].\n", op_prefix);
			exit(1);
		}
		
		if (fread(genotype_matrix_buff, sizeof(char), total_n_bytes, f_geno_mat) != (size_t)(total_n_bytes))
		{
			fprintf(stderr, "Could not read the genotype signal for %dx%d (%lu bytes) from %s..\n",
				vecsize(geno_sample_ids),
				vecsize(genotype_signal_regions),
				total_n_bytes,
				geno_matrix_op_fp);
			exit(1);
		}

		close_f(f_geno_mat, geno_matrix_op_fp);
	}

	fprintf(stderr, "Assigning the genotype vectors to regions..\n");
	for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
	{
		char* geno_sig = genotype_matrix_buff + i_reg * (size_t)(vecsize(geno_sample_ids));
		void** reg_info = new void* [10];
		reg_info[0] = geno_sig;
		genotype_signal_regions->at(i_reg)->data = reg_info;
	} // i_reg loop.
	//for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
	//{
	//	char* geno_sig = new char[vecsize(geno_sample_ids) + 2];
	//	if (fread(geno_sig, sizeof(char), vecsize(geno_sample_ids), f_geno_mat) != (size_t)vecsize(geno_sample_ids))
	//	{
	//		fprintf(stderr, "Could not read the genotype signal for %d/%d. region..\n", i_reg, vecsize(genotype_signal_regions));
	//		exit(1);
	//	}

	//	// Copy the genotype signal.
	//	void** reg_info = new void* [10];
	//	reg_info[0] = geno_sig;
	//	genotype_signal_regions->at(i_reg)->data = reg_info;
	//} // i_reg loop.

	return(genotype_signal_regions);
} // load_variant_genotype_signal_regions_per_matrix_subjects_regions option.

void binarize_variant_genotype_signal_regions_per_matrix_subjects_regions(vector<t_annot_region*>* genotype_signal_regions, vector<char*>* geno_sample_ids, 
	const bool compress_geno_matrix,
	const char* op_prefix)
{
	// Save the subject identifiers.
	char subjects_ids_list_fp[1000];
	sprintf(subjects_ids_list_fp, "%s_subjects.list", op_prefix);

	FILE* f_subjects_ids_list = open_f(subjects_ids_list_fp, "w");
	for (int i_s = 0; i_s < vecsize(geno_sample_ids); i_s++)
	{
		fprintf(f_subjects_ids_list, "%s\n", geno_sample_ids->at(i_s));
	} // i_s loop.
	close_f(f_subjects_ids_list, subjects_ids_list_fp);

	// Save the genotype matrix.
	//vector<char*>* chr_ids = get_chr_ids(genotype_signal_regions);

	// This should depend on the type of file from extension.
	char geno_matrix_op_fp[1000];
	if(compress_geno_matrix)
	{
		sprintf(geno_matrix_op_fp, "%s_genotypes.matrix.gz", op_prefix);
	}
	else
	{
		sprintf(geno_matrix_op_fp, "%s_genotypes.matrix", op_prefix);
	}

	// Save the regions.
	char regs_BED_fp[1000];
	sprintf(regs_BED_fp, "%s_variants.bed", op_prefix);
	FILE* f_BED = open_f(regs_BED_fp, "w");

	fprintf(stderr, "Copying genotype pool..\n");
	size_t geno_sig_mem_size = vecsize_t(genotype_signal_regions) * vecsize_t(geno_sample_ids);
	char* total_geno_mem_array = new char[vecsize_t(genotype_signal_regions) * vecsize_t(geno_sample_ids)];
	for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
	{
		// First, write the signal region to the bed file.
		if (genotype_signal_regions->at(i_reg)->strand == '-' ||
			genotype_signal_regions->at(i_reg)->strand == '+')
		{
			if (genotype_signal_regions->at(i_reg)->name == NULL)
			{
				genotype_signal_regions->at(i_reg)->name = new char[5];
				strcpy(genotype_signal_regions->at(i_reg)->name, ".");
			}

			// Translate the start and end to CODEBASE's start and end.
			fprintf(f_BED, "%s\t%d\t%d\t%s\t.\t%c\n", genotype_signal_regions->at(i_reg)->chrom,
					translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					genotype_signal_regions->at(i_reg)->name,
					genotype_signal_regions->at(i_reg)->strand);
		}
		else
		{
			fprintf(f_BED, "%s\t%d\t%d\n", genotype_signal_regions->at(i_reg)->chrom,
					translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
		}

		// Second, write the signal region to the bed file.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		//fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_geno_matrix_op);
		memcpy((total_geno_mem_array + i_reg * vecsize_t(geno_sample_ids)), cur_reg_geno_sig, vecsize_t(geno_sample_ids));
	} // i_reg loop.

	//fprintf(stderr, "\nClosing the file.\n");
	close_f(f_BED, regs_BED_fp);	
	//fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	if (compress_geno_matrix)
	{
		fprintf(stderr, "Saving genotype pool size of %.3f megabytes\n", (double)geno_sig_mem_size / (1024.0*1024));
		// Save the pool in one go:
		gzFile f_geno_matrix_op = gzopen(geno_matrix_op_fp, "wb");
		if (!f_geno_matrix_op)
		{
			perror("gzopen");
			return;
		}

		// 1 gig chunk size.
		size_t CHUNKSIZE = 1024 * 1024 * 1024;
		size_t ptr_pos = 0;
		while (1)
		{
			size_t n_bytes_2_write = MIN(CHUNKSIZE, geno_sig_mem_size - ptr_pos);
			if (n_bytes_2_write == 0)
			{
				break;
			}

			t_string::print_padded_string(stderr, '\r', 100, "Writing %.3f megabytes (%.3f megabytes left).", (double)n_bytes_2_write / (1024.0 * 1024.0), (double)(geno_sig_mem_size - ptr_pos) / (1024.0 * 1024.0));

			if (gzwrite(f_geno_matrix_op, total_geno_mem_array + ptr_pos, n_bytes_2_write) != (int)n_bytes_2_write)
			{
				perror("gzwrite");
				fprintf(stderr, "Could not save genotypes to %s\n", geno_matrix_op_fp);
				gzclose(f_geno_matrix_op);
				exit(1);
			}

			ptr_pos += n_bytes_2_write;
		} // chunk writing loop.
		fprintf(stderr, "Completed saving..\n");
		gzclose(f_geno_matrix_op);

		//// Compress to gz.
		//gzFile f_geno_matrix_op = gzopen(geno_matrix_op_fp, "wb");
		//if (!f_geno_matrix_op)
		//{
		//	perror("gzopen");
		//	return;
		//}

		//for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
		//{
		//	// Second, write the signal region to the bed file.
		//	void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		//	char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		//	//fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_geno_matrix_op);

		//	if (gzwrite(f_geno_matrix_op, cur_reg_geno_sig, vecsize_t(geno_sample_ids)) != vecsize_t(geno_sample_ids))
		//	{
		//		perror("gzwrite");
		//		gzclose(f_geno_matrix_op);
		//		return;
		//	}
		//} // i_reg loop.

		//gzclose(f_geno_matrix_op);
	}
	else
	{
		FILE* f_geno_matrix_op = open_f(geno_matrix_op_fp, "wb");
		fwrite(total_geno_mem_array, sizeof(char), geno_sig_mem_size, f_geno_matrix_op);
		close_f(f_geno_matrix_op, geno_matrix_op_fp);
		//// Write plain.
		//FILE* f_geno_matrix_op = open_f(geno_matrix_op_fp, "wb");
		//for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
		//{
		//	// Second, write the signal region to the bed file.
		//	void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		//	char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		//	fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_geno_matrix_op);
		//} // i_reg loop.
		//close_f(f_geno_matrix_op, geno_matrix_op_fp);
	} // uncompressed output check.

	//bool perform_check = true;
	bool perform_check = false;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	//clock_t cur_time = clock();
	vector<char*>* loaded_geno_sample_ids = new vector<char*>();
	vector<t_annot_region*>* loaded_sig_regs = load_variant_genotype_signal_regions_per_matrix_subjects_regions(op_prefix, loaded_geno_sample_ids);

	//fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (vecsize(loaded_geno_sample_ids) != vecsize(geno_sample_ids))
	{
		fprintf(stderr, "Dumping/Loading failed at subject identifiers.\n");
		exit(1);
	}

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(1);
	}

	for (int i_reg = 0; i_reg < (int)loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n",
				loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end,
				genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(1);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n", i_s,
					cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(1);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
} // binarize_variant_genotype_signal_regions_per_matrix_subjects_regions function.

static void* thread_callback_annot_region_genotype_save(void* __thread_info_ptr)
{
	void** thread_info = (void**)(__thread_info_ptr);
	int* int_vals = (int*)(thread_info[0]);
	int thr_i = int_vals[0];
	//int n_threads = int_vals[1];
	int var_i_start = int_vals[2];
	int var_i_end = int_vals[3];
	int compress_geno_matrix = int_vals[4];

	//double* dbl_vals = (double*)(thread_info[1]);

	if (compress_geno_matrix)
	{
		vector<t_annot_region*>* genotype_signal_regions = (vector<t_annot_region*>*)(thread_info[2]);
		vector<char*>* geno_sample_ids = (vector<char*>*)(thread_info[3]);

		t_string::print_padded_string(stderr, '\r', 100, "Genotype Saving Thread-%d: [%d-%d]", thr_i, var_i_start, var_i_end);

		char geno_matrix_op_fp[1000];
		sprintf(geno_matrix_op_fp, "genotype_saving_thread_%d.matrix.gz", thr_i);
		gzFile f_geno_matrix_op = gzopen(geno_matrix_op_fp, "wb");
		if (!f_geno_matrix_op)
		{
			perror("gzopen error..");
			exit(1);
			//return NULL;
		}

		// Go over the list of variants for this thread.
		for (int i_reg = var_i_start; i_reg < MIN(var_i_end, vecsize(genotype_signal_regions)); i_reg++)
		{
			// Second, write the signal region to the bed file.
			void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
			char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

			if (gzwrite(f_geno_matrix_op, cur_reg_geno_sig, vecsize_t(geno_sample_ids)) != vecsize(geno_sample_ids))
			{
				perror("gzwrite error..");
				exit(1);
				//gzclose(f_geno_matrix_op);
				//return NULL;
			}
		} // i_reg loop.

		gzclose(f_geno_matrix_op);
	} // compress data check
	else
	{
		vector<t_annot_region*>* genotype_signal_regions = (vector<t_annot_region*>*)(thread_info[2]);
		vector<char*>* geno_sample_ids = (vector<char*>*)(thread_info[3]);

		fprintf(stderr, "Thread-%d: [%d-%d]\n", thr_i, var_i_start, var_i_end);

		char geno_matrix_op_fp[1000];
		sprintf(geno_matrix_op_fp, "genotype_saving_thread_%d.matrix.gz", thr_i);
		FILE* f_geno_matrix_op = open_f(geno_matrix_op_fp, "wb");
		if (!f_geno_matrix_op)
		{
			fprintf(stderr, "Could not open %s for writing..\n", geno_matrix_op_fp);
			exit(1);
		}

		// Go over the list of variants for this thread.
		for (int i_reg = var_i_start; i_reg < MIN(var_i_end, vecsize(genotype_signal_regions)); i_reg++)
		{
			// Second, write the signal region to the bed file.
			void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
			char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
			fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_geno_matrix_op);
		} // i_reg loop.

		close_f(f_geno_matrix_op, geno_matrix_op_fp);
	} // compress data check
	
	return(NULL);
} // thread_callback_annot_region_genotype_save function.

// This function saves MatSubjReg formatted data using multiple threads.
void binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(vector<t_annot_region*>* genotype_signal_regions, vector<char*>* geno_sample_ids,
	const bool compress_geno_matrix,
	int n_threads,
	const char* op_prefix)
{
	fprintf(stderr, "Binarizing %d regions on %d subjects using %d threads..\n", vecsize(genotype_signal_regions), vecsize(geno_sample_ids),  n_threads);

	// We cannot have less that 4 threads (see below).
	if (n_threads < 4)
	{
		fprintf(stderr, "**RESETTING THE NUMBER OF REQUESTED THREADS TO 4**\n");
		fprintf(stderr, "**RESETTING THE NUMBER OF REQUESTED THREADS TO 4**\n");
		fprintf(stderr, "**RESETTING THE NUMBER OF REQUESTED THREADS TO 4**\n");
		fprintf(stderr, "**RESETTING THE NUMBER OF REQUESTED THREADS TO 4**\n");
		n_threads = 10;
	}

	// Save the subject identifiers.
	char subjects_ids_list_fp[1000];
	sprintf(subjects_ids_list_fp, "%s_subjects.list", op_prefix);

	FILE* f_subjects_ids_list = open_f(subjects_ids_list_fp, "w");
	for (int i_s = 0; i_s < vecsize(geno_sample_ids); i_s++)
	{
		fprintf(f_subjects_ids_list, "%s\n", geno_sample_ids->at(i_s));
	} // i_s loop.
	close_f(f_subjects_ids_list, subjects_ids_list_fp);

	// This should depend on the type of file from extension.
	char geno_matrix_op_fp[1000];
	if (compress_geno_matrix)
	{
		sprintf(geno_matrix_op_fp, "%s_genotypes.matrix.gz", op_prefix);
	}
	else
	{
		sprintf(geno_matrix_op_fp, "%s_genotypes.matrix", op_prefix);
	}

	// Save the regions.
	char regs_BED_fp[1000];
	sprintf(regs_BED_fp, "%s_variants.bed", op_prefix);
	FILE* f_BED = open_f(regs_BED_fp, "w");

	for (int i_reg = 0; i_reg < vecsize(genotype_signal_regions); i_reg++)
	{
		// First, write the signal region to the bed file.
		if (genotype_signal_regions->at(i_reg)->strand == '-' ||
			genotype_signal_regions->at(i_reg)->strand == '+')
		{
			if (genotype_signal_regions->at(i_reg)->name == NULL)
			{
				genotype_signal_regions->at(i_reg)->name = new char[5];
				strcpy(genotype_signal_regions->at(i_reg)->name, ".");
			}

			// Translate the start and end to CODEBASE's start and end.
			fprintf(f_BED, "%s\t%d\t%d\t%s\t.\t%c\n", genotype_signal_regions->at(i_reg)->chrom,
				translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				genotype_signal_regions->at(i_reg)->name,
				genotype_signal_regions->at(i_reg)->strand);
		}
		else
		{
			fprintf(f_BED, "%s\t%d\t%d\n", genotype_signal_regions->at(i_reg)->chrom,
				translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
		}
	} // i_reg loop.

	//fprintf(stderr, "\nClosing the file.\n");
	close_f(f_BED, regs_BED_fp);
	//fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	int n_vars_per_thread = MAX(1, ceil(vecsize(genotype_signal_regions) / ((double)n_threads)));

	int cur_var_start_i = 0;
	vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();

	for (int i_thr = 0; i_thr < n_threads; i_thr++)
	{
		if (cur_var_start_i >= vecsize(genotype_signal_regions))
		{
			fprintf(stderr, "Early stopping @ thread_i=%d\n", i_thr);
			break;
		}

		void** thread_info = new void* [20];
		int* int_vals = new int[20];
		int_vals[0] = i_thr;
		int_vals[1] = n_threads;
		int_vals[2] = cur_var_start_i;
		int_vals[3] = MIN(cur_var_start_i + n_vars_per_thread, vecsize(genotype_signal_regions));
		int_vals[4] = compress_geno_matrix;
		thread_info[0] = int_vals;

		//double* dbl_vals = (double*)(thread_info[1]);
		double* dbl_vals = new double[20];
		thread_info[1] = dbl_vals;

		//vector<t_annot_region*>* genotype_signal_regions = (vector<t_annot_region*>*)(thread_info[2]);
		thread_info[2] = genotype_signal_regions;
		//vector<char*>* geno_sample_ids = (vector<char*>*)(thread_info[3]);
		thread_info[3] = geno_sample_ids;

		cur_var_start_i += n_vars_per_thread;

		t_ansi_thread* thread = new t_ansi_thread(thread_callback_annot_region_genotype_save, thread_info);
		threads->push_back(thread);
		thread->run_thread();
	} // i_thr loop.

	vector<char*>* per_thread_matrix_files = new vector<char*>();
	for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
	{
		threads->at(i_thr)->wait_thread();
		t_string::print_padded_string(stderr, '\r', 80, "Geno-Save-Thread-%d finished..", i_thr);
		char cur_mat_file[1000];
		sprintf(cur_mat_file, "genotype_saving_thread_%d.matrix.gz", i_thr);

		per_thread_matrix_files->push_back(t_string::copy_me_str(cur_mat_file));
	} // i_thr loop.
		
	fprintf(stderr, "Concatenating genotype matrices..                          \n");
	concatenateGzipFiles(geno_matrix_op_fp, per_thread_matrix_files);

	for (int i_thr = 0; i_thr < vecsize(threads); i_thr++)
	{
		// Delete this file to make sure it does not interfere??
		delete_file(per_thread_matrix_files->at(i_thr));
	} // i_thr loop.

	//bool perform_check = true;
	bool perform_check = false;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	fprintf(stderr, "***Sanity checking multithread saved matrix genotype matrix***\n");

	//clock_t cur_time = clock();
	vector<char*>* loaded_geno_sample_ids = new vector<char*>();
	vector<t_annot_region*>* loaded_sig_regs = load_variant_genotype_signal_regions_per_matrix_subjects_regions(op_prefix, loaded_geno_sample_ids);

	//fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (vecsize(loaded_geno_sample_ids) != vecsize(geno_sample_ids))
	{
		fprintf(stderr, "Dumping/Loading failed at subject identifiers.\n");
		exit(1);
	}

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(1);
	}

	for (int i_reg = 0; i_reg < (int)loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n",
				loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end,
				genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(1);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n", i_s,
					cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(1);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
} // binarize_variant_genotype_signal_regions_per_matrix_subjects_regions function.

// Following can dump the supplied regions or load them then dump.
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp)
{
	if (genotype_signal_regions == NULL)
	{
		if (t_string::compare_strings(variant_geno_sig_regs_BED_fp, "stdin") || check_file(variant_geno_sig_regs_BED_fp))
		{
			genotype_signal_regions = load_variant_genotype_signal_regions(variant_geno_sig_regs_BED_fp, geno_sample_ids);
		}		
		else
		{
			fprintf(stderr, "Signal regions are not supplied and could not load them using %s.\n", variant_geno_sig_regs_BED_fp);
			exit(1);
		}
	}

	vector<char*>* chr_ids = get_chr_ids(genotype_signal_regions);

	// This should depend on the type of file from extension.
	FILE* f_op = open_f(op_fp, "wb");
	int n_chrs = (int)chr_ids->size();
	fwrite(&n_chrs, sizeof(int), 1, f_op);
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		char cur_chr[1000];
		strcpy(cur_chr, chr_ids->at(i_chr));
		fwrite(cur_chr, sizeof(char), 1000, f_op);
	} // i_chr loop.

	int sample_size = (int)geno_sample_ids->size();
	int n_regs = (int)genotype_signal_regions->size();

	fwrite(&sample_size, sizeof(int), 1, f_op);
	fwrite(&n_regs, sizeof(int), 1, f_op);
	for (int i_reg = 0; i_reg < (int)genotype_signal_regions->size(); i_reg++)
	{
		// Write the chromosome index.
		int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		fwrite(&i_chr, sizeof(int), 1, f_op);
		fwrite(&(reg_BED_start), sizeof(int), 1, f_op);
		fwrite(&(reg_BED_end), sizeof(int), 1, f_op);

		// Write the regions's name.
		if (genotype_signal_regions->at(i_reg)->name == NULL)
		{
			genotype_signal_regions->at(i_reg)->name = t_string::copy_me_str(".");
		}

		int l_reg_name_str = t_string::string_length(genotype_signal_regions->at(i_reg)->name);
		fwrite(&l_reg_name_str, sizeof(int), 1, f_op);
		fwrite(genotype_signal_regions->at(i_reg)->name, sizeof(char), l_reg_name_str, f_op);

		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		fwrite(cur_reg_geno_sig, sizeof(char), (int)geno_sample_ids->size(), f_op);
	} // i_reg loop.

	fprintf(stderr, "Closing the file.\n");
	close_f(f_op, op_fp);
	fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	bool perform_check = false;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	clock_t cur_time = clock();
	vector<t_annot_region*>* loaded_sig_regs = load_binarized_variant_genotype_signal_regions(op_fp, geno_sample_ids);
	fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(1);
	}

	for (int i_reg = 0; i_reg < (int)loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n", 
					loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end, 
					genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(1);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n", i_s,
						cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(1);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
}

vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids)
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
		//int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		//int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		//int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

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
		char* cur_reg_geno_sig = new char[sample_size + 2];
		//fread(cur_reg_geno_sig, sizeof(char), sample_size, f_bin_geno_sig_regs);
		read_bin_char_array(f_bin_geno_sig_regs, cur_reg_geno_sig, sample_size);

		cur_reg_info[0] = cur_reg_geno_sig;
		cur_reg_info[1] = NULL;

		reg->data = cur_reg_info;

		geno_sig_regs->push_back(reg);
	} // i_reg loop.

	// Close the file.
	close_f(f_bin_geno_sig_regs, bin_geno_sig_bed_fp);

	return(geno_sig_regs);
}

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids)
{
	vector<t_annot_region*>* link_var_regs = new vector<t_annot_region*>();
	char tok_buff[1000];
	FILE* f_link_variant_genotype_signal = open_f(link_variant_genotype_signal_fp, "r");

	if (f_link_variant_genotype_signal == NULL)
	{
		fprintf(stderr, "Could not open %s\n", link_variant_genotype_signal_fp);
		exit(1);
	}

	int i_reg = 0;
	while (1)
	{
		char* cur_reg_line = getline(f_link_variant_genotype_signal);
		if (cur_reg_line == NULL)
		{
			break;
		}
		else
		{
			i_reg++;
		}

		if (i_reg % 10000 == 0)
		{
			t_string::print_padded_string(stderr, '\r', 100, "Parsing genotype signal for %d. region.", i_reg);
		}

		//char* cur_reg_line = (char*)(link_var_regs->at(i_reg)->data);
		//t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		char* cur_var_geno_signal = new char[(int)geno_sample_ids->size()];

		t_annot_region* cur_reg = get_empty_region();

		// Copy the name.
		int i_cur_char = 0;
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->chrom = t_string::copy_me_str(tok_buff);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->start = translate_coord(atoi(tok_buff), BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->end = translate_coord(atoi(tok_buff), BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_reg->strand = '+';

		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->name = t_string::copy_me_str(tok_buff);

		//fprintf(stderr, "%s\t%d\t%d: %s            \r", link_var_regs->at(i_reg)->chrom, link_var_regs->at(i_reg)->start, link_var_regs->at(i_reg)->end, link_var_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < (int)geno_sample_ids->size(); i_s++)
		{
			if (t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char) == false)
			{
				fprintf(stderr, "Could not parse the %d. sample's genotype signal entry.\n", i_s);
				exit(1);
			}

			cur_var_geno_signal[i_s] = (char)(atoi(tok_buff));
		} // i_tok loop.

		void** link_var_reg_info = new void*[3];
		link_var_reg_info[0] = cur_var_geno_signal;
		link_var_reg_info[1] = NULL;

		cur_reg->data = link_var_reg_info;

		link_var_regs->push_back(cur_reg);
	} // i_reg loop.

	close_f(f_link_variant_genotype_signal, link_variant_genotype_signal_fp);

	return(link_var_regs);
}

char* copy_nucs(char* seq)
{
	char* copy = new char[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(char));

	for (int i = 0; i < (int)strlen(seq); i++)
	{
		copy[i] = seq[i];
	} // i loop.

	return(copy);
}

int* copy_bps(char* seq, int* bps)
{
	int* copy = new int[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(int));

	for (int i = 0; i < (int)strlen(seq); i++)
	{
		copy[i] = bps[i];
	} // i loop.

	return(copy);
}

void parse_variants_per_info_string(char* info_str,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	t_string_tokens* info_str_tokens = t_string::tokenize_by_chars(info_str, "=;");

	int i_tok = 0;
	while (i_tok < (int)info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = (info_str_tokens->at(i_tok + 1)->str())[0];
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	t_string::clean_tokens(info_str_tokens);
}

char mutate_nuc(char nuc, t_rng* rng)
{
	char cur_rand_nuc = toupper(nuc);
	while (cur_rand_nuc == toupper(nuc))
	{
		int rand_i = (int)(floor(4 * rng->random_double_ran3()));
		char rand_nucs[] = "ACGT";
		cur_rand_nuc = rand_nucs[rand_i];
	}

	return(cur_rand_nuc);
}

bool mutate_GC_2_AU(char* seq)
{
	// Count GC content.
	int n_gc = 0;
	for (int i = 0; i < (int)strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			n_gc++;
		}
	} // i loop.

	if (n_gc == 0)
	{
		return(false);
	}

	int i_rand = rand() % n_gc;

	int i_gc = 0;
	for (int i = 0; i < (int)strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			if (i_gc == i_rand)
			{
				// Choose a nuc.
				int n_rand = rand() % 2;
				char rand_nuc = (n_rand == 0) ? ('A') : ('U');
				fprintf(stderr, "%c->%c @ %d\n", seq[i], rand_nuc, i);
				seq[i] = rand_nuc;
				break;
			}

			i_gc++;
		}
	} // i loop.

	return(true);
}

char mutate_nuc_2_GC(char nuc, t_rng* rng)
{
	if (toupper(nuc) != 'G' &&
		toupper(nuc) != 'C')
	{
		int val = (int)(floor(2 * rng->random_double_ran3()));
		if (val == 0)
		{
			return('G');
		}
		else
		{
			return('C');
		}
	}
	else
	{
		return(toupper(nuc));
	}
}

char* get_random_seq(int l)
{
	char* seq = new char[l + 2];
	memset(seq, 0, (l + 1) * sizeof(char));
	for (int i = 0; i < l; i++)
	{
		int random = rand() % 4;
		seq[i] = num_2_nuc(random);
	} // i loop

	return(seq);
}
