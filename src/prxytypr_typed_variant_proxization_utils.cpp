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
#include "prxytypr_seed_manager.h"
#include "prxytypr_ansi_thread.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_histogram.h"
#include "prxytypr_matrix_linalg_utils.h"
#include "prxytypr_vector_macros.h"
#include "prxytypr_proxytyper.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

const bool __DUMP_PROXIZATION_MSGS__ = false;

// This function simply maps the coordinates to a new coordinate system, it does not do any type of permutations.
void generic_map_regions(vector<t_annot_region*>* target_regs, char* source_coords_bed_fp, char* dest_coords_bed_fp)
{
	vector<t_annot_region*>* source_coord_regs = load_BED(source_coords_bed_fp);
	vector<t_annot_region*>* dest_coord_regs = load_BED(dest_coords_bed_fp);

	if (source_coord_regs->size() != dest_coord_regs->size())
	{
		fprintf(stderr, "Source and destination coordinate sizes are not the same: %d/%d\n", (int)source_coord_regs->size(), (int)dest_coord_regs->size());
		exit(1);
	}

	fprintf(stderr, "Generic mapping coordinates of %d regions using %d source regions to %d destination regions.\n", (int)target_regs->size(),
		(int)source_coord_regs->size(), (int)dest_coord_regs->size());

	for (int i_reg = 0; i_reg < (int)source_coord_regs->size(); i_reg++)
	{
		source_coord_regs->at(i_reg)->data = dest_coord_regs->at(i_reg);
	} // i_reg loop.

	for (int i_reg = 0; i_reg < (int)target_regs->size(); i_reg++)
	{
		void** new_info = new void* [3];
		new_info[0] = target_regs->at(i_reg)->data;
		new_info[1] = NULL;
		new_info[2] = NULL;
		target_regs->at(i_reg)->data = new_info;
	} // i_reg loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(target_regs, source_coord_regs, true);

	//// This ensures that every coordinate we have in the source is mapped to sth in the target regions
	//// This is reasonable because we expect that the source-destination coordiantes to be an strict subset of the target regions (???).
	//if (vecsize(intersects) != vecsize(source_coord_regs))
	//{
	//	fprintf(stderr, "Sanity check failed @ %s(%d): Some source coordinates are not matched %d/%d\n", __FILE__, __LINE__,
	//			vecsize(intersects), vecsize(source_coord_regs));

	//	exit(1);
	//}

	fprintf(stderr, "Mapping coordinates of %d/%d target regions.\n", (int)intersects->size(), vecsize(target_regs));

	for (int i_int = 0; i_int < (int)intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_trgt_reg = int_info->src_reg;
		t_annot_region* cur_src_reg = int_info->dest_reg;

		if (cur_trgt_reg->start == cur_src_reg->start &&
			cur_trgt_reg->end == cur_src_reg->end)
		{
			void** trgt_reg_info = (void**)(cur_trgt_reg->data);
			trgt_reg_info[1] = cur_src_reg->data;
		}
	} // i_int loop.

	for (int i_reg = 0; i_reg < (int)target_regs->size(); i_reg++)
	{
		void** trgt_reg_info = (void**)(target_regs->at(i_reg)->data);
		if (trgt_reg_info[1] == NULL)
		{
			fprintf(stderr, "Could not map %s:%d (This panel may be on another coordinate system?)\n", target_regs->at(i_reg)->chrom, target_regs->at(i_reg)->start);
			exit(1);
		}

		t_annot_region* remapped_reg = (t_annot_region*)(trgt_reg_info[1]);

		target_regs->at(i_reg)->start = remapped_reg->start;
		target_regs->at(i_reg)->end = remapped_reg->end;
		target_regs->at(i_reg)->name = t_string::copy_me_str(remapped_reg->name);

		// Restore the region information.
		void** orig_reg_info = (void**)(trgt_reg_info[0]);
		target_regs->at(i_reg)->data = orig_reg_info;
	} // i_reg loop.

	fprintf(stderr, "Mapped %d region coordinates.\n", (int)target_regs->size());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is the main function to anonymize coordinates and genetic maps.
// This function takes tags and targets separately since we need the genetic map only for tags.
void anonymize_tag_target_genetic_map_coordinates_per_ref_query_variants(char* tag_variants_BED_fp, char* tag_target_variants_BED_fp,
	char* genetic_map_file, double cM_noise_SD, unsigned int max_base_posn, char* op_dir)
{
	fprintf(stderr, "Generating coordinate/genetic map anonymization mappings:\n\
\\sigma_noise=%.4f cM\n\
Anonymized coordinates: [1-%u]\n", cM_noise_SD, max_base_posn);

	vector<t_annot_region*>* anon_cM_tag_var_regs = anonymize_genetic_distances_ret_regs(tag_variants_BED_fp, genetic_map_file, cM_noise_SD);

	for (int i_reg = 0; i_reg < (int)anon_cM_tag_var_regs->size(); i_reg++)
	{
		anon_cM_tag_var_regs->at(i_reg)->score = 1;

		// This is the original coords.
		anon_cM_tag_var_regs->at(i_reg)->data = duplicate_region(anon_cM_tag_var_regs->at(i_reg));
	} // i_reg loop.

	// Validate complete overlap between the tag variants and the tag+target variant set in reference panel.
	vector<t_annot_region*>* raw_tag_target_regs = load_BED(tag_target_variants_BED_fp);

	vector<t_annot_region*>* intersect_validator = intersect_annot_regions(raw_tag_target_regs, anon_cM_tag_var_regs, false);
	if (vecsize(intersect_validator) != vecsize(anon_cM_tag_var_regs))
	{
		fprintf(stderr, "%s(%d): The tags are not totally covered by reference panel's variant set.\n", __FILE__, __LINE__);
		exit(1);
	}

	vector<t_annot_region*>* target_regs = exclude_annot_regions(raw_tag_target_regs, anon_cM_tag_var_regs);

	//vector<t_annot_region*>* target_regs = load_BED(target_variants_BED_fp);
	fprintf(stderr, "Loaded %d target regions.\n", (int)target_regs->size());
	for (int i_reg = 0; i_reg < (int)target_regs->size(); i_reg++)
	{
		target_regs->at(i_reg)->score = 0;

		// This is the original coords.
		target_regs->at(i_reg)->data = duplicate_region(target_regs->at(i_reg));
	} // i_reg loop.

	vector<t_annot_region*>* tag_target_regs = new vector<t_annot_region*>();
	tag_target_regs->insert(tag_target_regs->end(), anon_cM_tag_var_regs->begin(), anon_cM_tag_var_regs->end());
	tag_target_regs->insert(tag_target_regs->end(), target_regs->begin(), target_regs->end());

	fprintf(stderr, "Processing %d tag/target regions.\n", (int)target_regs->size());

	sort(tag_target_regs->begin(), tag_target_regs->end(), sort_regions);

	int cur_var_coord = 1;
	int delta_pos = max_base_posn / ((int)tag_target_regs->size());
	for (int i_reg = 0; i_reg < (int)tag_target_regs->size(); i_reg++)
	{
		tag_target_regs->at(i_reg)->start = cur_var_coord;
		tag_target_regs->at(i_reg)->end = cur_var_coord;
		cur_var_coord += delta_pos;
	} // i_reg loop.

	char* chr_id = tag_target_regs->at(0)->chrom;

	// Write anonymized genetic map file in tag coordinates.
	char anon_genetic_map_fp[1000];
	sprintf(anon_genetic_map_fp, "%s/%s.map", op_dir, chr_id);
	FILE* f_anon_genetic_map = open_f(anon_genetic_map_fp, "w");

	//fprintf(f_anon_genetic_map, "position COMBINED_rate(cM / Mb) Genetic_Map(cM)\n");

	for (int i_reg = 0; i_reg < (int)anon_cM_tag_var_regs->size(); i_reg++)
	{
		//fprintf(f_anon_genetic_map, "%d 0.000 %.8f\n",
		//	anon_cM_tag_var_regs->at(i_reg)->end,
		//	anon_cM_tag_var_regs->at(i_reg)->dbl_score);

		fprintf(f_anon_genetic_map, "%s . %.4f %d\n",
			anon_cM_tag_var_regs->at(i_reg)->chrom,
			anon_cM_tag_var_regs->at(i_reg)->dbl_score,
			anon_cM_tag_var_regs->at(i_reg)->start);
	} // i_reg loop.

	close_f(f_anon_genetic_map, anon_genetic_map_fp);

	// Write tag-target mapping coordinates and anonymize the variant identifier.
	vector<t_annot_region*>* tag_target_orig_regs = new vector<t_annot_region*>();
	vector<t_annot_region*>* tag_target_mapping_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)tag_target_regs->size(); i_reg++)
	{
		// Anonymize the variant identifier.
		if (tag_target_regs->at(i_reg)->name != NULL)
		{
			delete[] tag_target_regs->at(i_reg)->name;
		}

		char anon_reg_name[1000];
		if (tag_target_regs->at(i_reg)->score == 0)
		{
			sprintf(anon_reg_name, "untyped%d_A_G_0.000", i_reg);
		}
		else
		{
			sprintf(anon_reg_name, "typed%d_A_G_0.000", i_reg);
		}

		tag_target_regs->at(i_reg)->name = t_string::copy_me_str(anon_reg_name);

		// Add the original coordinates.
		t_annot_region* orig_reg = (t_annot_region*)(tag_target_regs->at(i_reg)->data);
		tag_target_orig_regs->push_back(orig_reg);

		// Add the mapping coordinates.
		tag_target_mapping_regs->push_back(tag_target_regs->at(i_reg));
	} // i_reg loop.

	char tag_target_orig_coords_fp[1000];
	sprintf(tag_target_orig_coords_fp, "%s/original_coords.bed", op_dir);
	dump_BED(tag_target_orig_coords_fp, tag_target_orig_regs);

	char tag_target_mapping_coords_fp[1000];
	sprintf(tag_target_mapping_coords_fp, "%s/mapping_coords.bed", op_dir);
	dump_BED(tag_target_mapping_coords_fp, tag_target_mapping_regs);
} // anonymize_tag_target_genetic_map_coordinates function.

// This function is deprecated for anonymizing the coordinates/maps. Use anonymize_tag_target_genetic_map_coordinates_per_ref_query_variants that contains more sanity/error checks.
void anonymize_tag_target_genetic_map_coordinates(char* tag_variants_BED_fp, char* target_variants_BED_fp,
	char* genetic_map_file, double cM_noise_SD, unsigned int max_base_posn, char* op_dir)
{
	fprintf(stderr, "Generating coordinate/genetic map anonymization mappings:\n\
\\sigma_noise=%.4f cM\n\
Anonymized coordinates: [1-%u]\n", cM_noise_SD, max_base_posn);

	vector<t_annot_region*>* anon_cM_tag_var_regs = anonymize_genetic_distances_ret_regs(tag_variants_BED_fp, genetic_map_file, cM_noise_SD);

	for (int i_reg = 0; i_reg < (int)anon_cM_tag_var_regs->size(); i_reg++)
	{
		anon_cM_tag_var_regs->at(i_reg)->score = 1;

		// This is the original coords.
		anon_cM_tag_var_regs->at(i_reg)->data = duplicate_region(anon_cM_tag_var_regs->at(i_reg));
	} // i_reg loop.

	vector<t_annot_region*>* target_regs = load_BED(target_variants_BED_fp);
	fprintf(stderr, "Loaded %d target regions.\n", (int)target_regs->size());
	for (int i_reg = 0; i_reg < (int)target_regs->size(); i_reg++)
	{
		target_regs->at(i_reg)->score = 0;

		// This is the original coords.
		target_regs->at(i_reg)->data = duplicate_region(target_regs->at(i_reg));
	} // i_reg loop.

	vector<t_annot_region*>* tag_target_regs = new vector<t_annot_region*>();
	tag_target_regs->insert(tag_target_regs->end(), anon_cM_tag_var_regs->begin(), anon_cM_tag_var_regs->end());
	tag_target_regs->insert(tag_target_regs->end(), target_regs->begin(), target_regs->end());

	fprintf(stderr, "Processing %d tag/target regions.\n", (int)target_regs->size());

	sort(tag_target_regs->begin(), tag_target_regs->end(), sort_regions);

	int cur_var_coord = 1;
	int delta_pos = max_base_posn / ((int)tag_target_regs->size());
	for (int i_reg = 0; i_reg < (int)tag_target_regs->size(); i_reg++)
	{
		tag_target_regs->at(i_reg)->start = cur_var_coord;
		tag_target_regs->at(i_reg)->end = cur_var_coord;
		cur_var_coord += delta_pos;
	} // i_reg loop.

	char* chr_id = tag_target_regs->at(0)->chrom;

	// Write anonymized genetic map file in tag coordinates.
	char anon_genetic_map_fp[1000];
	sprintf(anon_genetic_map_fp, "%s/%s.map", op_dir, chr_id);
	FILE* f_anon_genetic_map = open_f(anon_genetic_map_fp, "w");

	//fprintf(f_anon_genetic_map, "position COMBINED_rate(cM / Mb) Genetic_Map(cM)\n");

	for (int i_reg = 0; i_reg < (int)anon_cM_tag_var_regs->size(); i_reg++)
	{
		//fprintf(f_anon_genetic_map, "%d 0.000 %.8f\n",
		//	anon_cM_tag_var_regs->at(i_reg)->end,
		//	anon_cM_tag_var_regs->at(i_reg)->dbl_score);

		fprintf(f_anon_genetic_map, "%s . %.4f %d\n",
				anon_cM_tag_var_regs->at(i_reg)->chrom,
				anon_cM_tag_var_regs->at(i_reg)->dbl_score,
				anon_cM_tag_var_regs->at(i_reg)->start);
	} // i_reg loop.

	close_f(f_anon_genetic_map, anon_genetic_map_fp);

	// Write tag-target mapping coordinates and anonymize the variant identifier.
	vector<t_annot_region*>* tag_target_orig_regs = new vector<t_annot_region*>();
	vector<t_annot_region*>* tag_target_mapping_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < (int)tag_target_regs->size(); i_reg++)
	{
		// Anonymize the variant identifier.
		if (tag_target_regs->at(i_reg)->name != NULL)
		{
			delete[] tag_target_regs->at(i_reg)->name;
		}

		char anon_reg_name[1000];
		if (tag_target_regs->at(i_reg)->score == 0)
		{
			sprintf(anon_reg_name, "untyped%d_A_G_0.000", i_reg);
		}
		else 
		{
			sprintf(anon_reg_name, "typed%d_A_G_0.000", i_reg);
		}

		tag_target_regs->at(i_reg)->name = t_string::copy_me_str(anon_reg_name);

		// Add the original coordinates.
		t_annot_region* orig_reg = (t_annot_region*)(tag_target_regs->at(i_reg)->data);
		tag_target_orig_regs->push_back(orig_reg);

		// Add the mapping coordinates.
		tag_target_mapping_regs->push_back(tag_target_regs->at(i_reg));
	} // i_reg loop.

	char tag_target_orig_coords_fp[1000];
	sprintf(tag_target_orig_coords_fp, "%s/original_coords.bed", op_dir);
	dump_BED(tag_target_orig_coords_fp, tag_target_orig_regs);

	char tag_target_mapping_coords_fp[1000];
	sprintf(tag_target_mapping_coords_fp, "%s/mapping_coords.bed", op_dir);
	dump_BED(tag_target_mapping_coords_fp, tag_target_mapping_regs);
} // anonymize_tag_target_genetic_map_coordinates function.

vector<t_annot_region*>* anonymize_genetic_distances_ret_regs(char* variants_BED_fp, char* genetic_map_file, double cM_noise_SD)
{
	vector<t_annot_region*>* anon_cM_var_regs = load_BED(variants_BED_fp);
	fprintf(stderr, "Loaded %d variant regions.\n", (int)anon_cM_var_regs->size());

	//char cur_chr_recombination_rate_fp[1000];
	//sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_sig_regs->chr_ids->at(i_chr));
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(genetic_map_file);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", genetic_map_file);
		exit(1);
	}

	fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", (int)cur_chrom_recomb_regs->size());
	for (int i_reg = 0; i_reg < (int)anon_cM_var_regs->size(); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(anon_cM_var_regs->at(i_reg), cur_chrom_recomb_regs);
		anon_cM_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.

	sort(anon_cM_var_regs->begin(), anon_cM_var_regs->end(), sort_regions);

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	for (int i_reg = 0; i_reg < (int)anon_cM_var_regs->size(); i_reg++)
	{
		// Add the noise.
		double cur_cM_noise = rng->random_gaussian_double_ran3() * cM_noise_SD;

		anon_cM_var_regs->at(i_reg)->dbl_score += cur_cM_noise;
	} // i_reg loop.

	sort(anon_cM_var_regs->begin(), anon_cM_var_regs->end(), sort_regions_per_dbl_score);

	vector<t_annot_region*>* coord_sorted_var_regs = load_BED(variants_BED_fp);
	sort(coord_sorted_var_regs->begin(), coord_sorted_var_regs->end(), sort_regions);

	double base_cM = anon_cM_var_regs->at(0)->dbl_score;

	// Save the final distances. Push everything back to 0;
	fprintf(stderr, "Assigning genetic distances.\n");
	for (int i_reg = 0; i_reg < (int)coord_sorted_var_regs->size(); i_reg++)
	{
		//anon_cM_var_regs->at(i_reg)->dbl_score = MAX(0, anon_cM_var_regs->at(i_reg)->dbl_score);

		//coord_sorted_var_regs->at(i_reg)->dbl_score = anon_cM_var_regs->at(i_reg)->dbl_score;
		coord_sorted_var_regs->at(i_reg)->dbl_score = anon_cM_var_regs->at(i_reg)->dbl_score - base_cM;
	} // i_reg loop.

	return(coord_sorted_var_regs);
} // anonymize_genetic_distances option.

void anonymize_genetic_distances(char* variants_BED_fp, char* genetic_map_file, double cM_noise_SD, char* op_fp)
{
	vector<t_annot_region*>* anon_dist_regs = anonymize_genetic_distances_ret_regs(variants_BED_fp, genetic_map_file, cM_noise_SD);

	// Write.
	FILE* f_op = open_f(op_fp, "w");

	//fprintf(f_op, "position COMBINED_rate(cM / Mb) Genetic_Map(cM)\n");

	for (int i_reg = 0; i_reg < (int)anon_dist_regs->size(); i_reg++)
	{
		//anon_cM_var_regs->at(i_reg)->dbl_score = MAX(0, anon_cM_var_regs->at(i_reg)->dbl_score);

		fprintf(f_op, "%s . %.4f %d\n",
				anon_dist_regs->at(i_reg)->chrom,
				anon_dist_regs->at(i_reg)->dbl_score,
				anon_dist_regs->at(i_reg)->start);
	} // i_reg loop.

	close_f(f_op, op_fp);
} // anonymize_genetic_distances option.

vector<t_annot_region*>* load_tag_permute_proxy_mapping_regs(char* permute_proxy_mappings_BED_fp)
{
	vector<t_annot_region*>* perm_var_orig_regs = new vector<t_annot_region*>();

	vector<t_annot_region*>* all_var_regs = load_BED_with_line_information(permute_proxy_mappings_BED_fp);
	fprintf(stderr, "Loaded %d permute-proxy mapping regions.\n", (int)all_var_regs->size());

	for (int i_var = 0; i_var < (int)all_var_regs->size(); i_var++)
	{
		char* cur_reg_line = (char*)(all_var_regs->at(i_var)->data);
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");

		/*
 		fprintf(f_op, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n", 
		geno_sig_regs->at(i_reg)->chrom,
		translate_coord(geno_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
		translate_coord(geno_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
		permuted_regs->at(i_reg)->chrom, 
		translate_coord(permuted_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
		translate_coord(permuted_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
		geno_sig_regs->at(i_reg)->score, permuted_regs->at(i_reg)->score,
		geno_sig_regs->at(i_reg)->dbl_score, permuted_regs->at(i_reg)->dbl_score,
		invert_geno);
		*/

		if ((int)toks->size() != 11)
		{
			fprintf(stderr, "Could not parse line: %s\n", cur_reg_line);
			exit(1);
		}

		char* orig_chrom = toks->at(0)->str();
		int orig_start = atoi(toks->at(1)->str());
		int orig_end = atoi(toks->at(2)->str());
		char* perm_chrom = toks->at(3)->str();
		int perm_start = atoi(toks->at(4)->str());
		int perm_end = atoi(toks->at(5)->str());
		//int orig_score = atoi(toks->at(6)->str());
		//int perm_score = atoi(toks->at(7)->str());
		double orig_dbl_score = atoi(toks->at(8)->str());
		double perm_dbl_score = atoi(toks->at(9)->str());
		int invert_geno = atoi(toks->at(10)->str());

		t_annot_region* cur_orig_reg = get_empty_region();
		cur_orig_reg->chrom = t_string::copy_me_str(orig_chrom);
		cur_orig_reg->start = translate_coord(orig_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_orig_reg->end = translate_coord(orig_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_orig_reg->score = 0;
		cur_orig_reg->dbl_score = orig_dbl_score;
		cur_orig_reg->strand = '+';
		t_annot_region* cur_perm_reg = get_empty_region();
		cur_perm_reg->chrom = t_string::copy_me_str(perm_chrom);
		cur_perm_reg->start = translate_coord(perm_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_perm_reg->end = translate_coord(perm_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_perm_reg->score = invert_geno;
		cur_perm_reg->dbl_score = perm_dbl_score;
		cur_perm_reg->strand = '+';
		cur_perm_reg->data = NULL; // This is important for sanity checks.

		cur_orig_reg->data = cur_perm_reg;

		perm_var_orig_regs->push_back(cur_orig_reg);
		t_string::clean_tokens(toks);
	} // i_var loop.

	fprintf(stderr, "Loaded %d permute-proxy mapping regions.\n", (int)perm_var_orig_regs->size());

	return(perm_var_orig_regs);
}

void build_augmenting_variants_mapper(char* augmenting_variants_BED_fp, char* panel_variants_BED_fp, double augment_prob, char* augmentation_mapper_op_prefix)
{
	fprintf(stderr, "Building tag augmenter model with %.3f probability..\n", augment_prob);

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	vector<t_annot_region*>* augmenting_var_regs = load_BED(augmenting_variants_BED_fp);
	sort(augmenting_var_regs->begin(), augmenting_var_regs->end(), sort_regions);
	fprintf(stderr, "Loaded %d augmenting variant regions.\n", vecsize(augmenting_var_regs));

	vector<t_annot_region*>* panel_var_regs = load_BED(panel_variants_BED_fp);
	fprintf(stderr, "Loaded %d whole panel variant regions.\n", vecsize(panel_var_regs));
	sort(panel_var_regs->begin(), panel_var_regs->end(), sort_regions);
	int l_indicator = (size_t)250 * (size_t)1000 * (size_t)1000;
	char* panel_var_posn_indicator = new char[l_indicator];
	memset(panel_var_posn_indicator, 0, sizeof(char) * l_indicator);
	for (int i_reg = 0; i_reg < vecsize(panel_var_regs); i_reg++)
	{
		panel_var_regs->at(i_reg)->score = i_reg;
		panel_var_posn_indicator[panel_var_regs->at(i_reg)->start] = 1;
	} // i_reg loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(panel_var_regs, augmenting_var_regs, false);
	if (intersects->size() != augmenting_var_regs->size())
	{
		fprintf(stderr, "The augmenting variants do not completely overlap: %d/%d\n", vecsize(intersects), vecsize(augmenting_var_regs));
		exit(1);
	}

	// Place the augmenting regions.
	char original_augmenting_regs_BED_fp[1000];
	sprintf(original_augmenting_regs_BED_fp, "%s_original.bed", augmentation_mapper_op_prefix);
	char mapped_augmenting_regs_BED_fp[1000];
	sprintf(mapped_augmenting_regs_BED_fp, "%s_mapped.bed", augmentation_mapper_op_prefix);

	FILE* f_original_augmenting_regs_BED = open_f(original_augmenting_regs_BED_fp, "w");
	FILE* f_mapped_augmenting_regs_BED = open_f(mapped_augmenting_regs_BED_fp, "w");
	int n_augmenting_vars = 0;
	for (int i_reg = 0; i_reg < vecsize(augmenting_var_regs); i_reg++)
	{
		//int panel_reg_i = augmenting_var_regs->at(i_reg)->score;
		if (rng->random_double_ran3() < augment_prob)
		{
			int cur_augment_reg_posn = augmenting_var_regs->at(i_reg)->start;

			int prev_augment_reg_posn = MIN(panel_var_regs->at(0)->start, cur_augment_reg_posn-1000);
			if (i_reg > 0)
			{
				prev_augment_reg_posn = augmenting_var_regs->at(i_reg - 1)->start;
			}

			if ((cur_augment_reg_posn - prev_augment_reg_posn) > 100)
			{
				// Find a position that does not have any panel variants.
				int posn_candidate = int(prev_augment_reg_posn + (cur_augment_reg_posn - prev_augment_reg_posn) * rng->random_double_ran3());
				while (panel_var_posn_indicator[posn_candidate] == 1)
				{
					posn_candidate = int(prev_augment_reg_posn + (cur_augment_reg_posn - prev_augment_reg_posn) * rng->random_double_ran3());
				} // position searching loop.

				if (__DUMP_PROXIZATION_MSGS__)
				{
					fprintf(stderr, "Augmenting var-%d: %s:%d-%d ;; selected posn_candidate=%d\n", i_reg,
						augmenting_var_regs->at(i_reg)->chrom, augmenting_var_regs->at(i_reg)->start, augmenting_var_regs->at(i_reg)->end,
						posn_candidate);
				}

				fprintf(f_original_augmenting_regs_BED, "%s\t%d\t%d\t%s\t.\t+\n", augmenting_var_regs->at(i_reg)->chrom,
					translate_coord(augmenting_var_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(augmenting_var_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					augmenting_var_regs->at(i_reg)->name);

				fprintf(f_mapped_augmenting_regs_BED, "%s\t%d\t%d\t%s_AUGMENTED_%d\t+\n", augmenting_var_regs->at(i_reg)->chrom,
					translate_coord(posn_candidate, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(posn_candidate, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					augmenting_var_regs->at(i_reg)->name, i_reg);

				n_augmenting_vars++;
			}
		} // probability check.
	} // i_reg loop.	

	close_f(f_mapped_augmenting_regs_BED, NULL);
	close_f(f_original_augmenting_regs_BED, NULL);

	fprintf(stderr, "Saved %d augmenting tag variants..\n", n_augmenting_vars);
} // build_augmenting_variants_mapper function.

void augment_tag_variants(char* augmentation_mapper_prefix, char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	int geno_augment_type, char* op_matbed)
{
	char original_augmenting_regs_BED_fp[1000];
	sprintf(original_augmenting_regs_BED_fp, "%s_original.bed", augmentation_mapper_prefix);
	char mapped_augmenting_regs_BED_fp[1000];
	sprintf(mapped_augmenting_regs_BED_fp, "%s_mapped.bed", augmentation_mapper_prefix);

	vector<t_annot_region*>* original_augmenting_panel_regs = load_BED(original_augmenting_regs_BED_fp);
	vector<t_annot_region*>* mapped_augmenting_panel_regs = load_BED(mapped_augmenting_regs_BED_fp);

	if (vecsize(original_augmenting_panel_regs) != vecsize(mapped_augmenting_panel_regs))
	{
		fprintf(stderr, "%s(%d): Sanity check failed Mapped and augmenting regions do not match in numbers: %d/%d\n",
			__FILE__, __LINE__, vecsize(original_augmenting_panel_regs), vecsize(mapped_augmenting_panel_regs));

		exit(1);
	}

	fprintf(stderr, "Augmenting %d regions with augmentation type: %d (%d: ZERO, %d: COPY)\n", vecsize(original_augmenting_panel_regs), 
			geno_augment_type, VAR_AUGMENT_ZERO, VAR_AUGMENT_COPY);

	for (int i_reg = 0; i_reg < vecsize(original_augmenting_panel_regs); i_reg++)
	{
		original_augmenting_panel_regs->at(i_reg)->score = i_reg;
	} // i_reg loop.

	vector<t_annot_region*>* panel_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Augmenting %d region on %d regions on %d subjects.\n", vecsize(original_augmenting_panel_regs), vecsize(panel_var_regs), vecsize(sample_ids));

	vector<t_annot_region*>* intersects = intersect_annot_regions(panel_var_regs, original_augmenting_panel_regs, false);
	if (intersects->size() != original_augmenting_panel_regs->size())
	{
		fprintf(stderr, "The augmenting variants do not completely overlap: %d/%d\n", vecsize(intersects), vecsize(original_augmenting_panel_regs));
		exit(1);
	}

	vector<t_annot_region*>* augmenting_panel_var_regs = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);

		t_annot_region* panel_var_reg = int_info->src_reg;
		t_annot_region* augmenting_reg = int_info->dest_reg;

		// Add this augmenting region.
		t_annot_region* cur_augmenting_panel_var_reg = duplicate_region(mapped_augmenting_panel_regs->at(augmenting_reg->score));
		void** panel_var_reg_info = (void**)(panel_var_reg->data);
		char* panel_var_reg_geno_sig = (char*)(panel_var_reg_info[0]);
		char* cur_augmenting_panel_var_reg_geno = new char[vecsize(sample_ids) + 2];

		if (geno_augment_type == VAR_AUGMENT_COPY)
		{
			memcpy(cur_augmenting_panel_var_reg_geno, panel_var_reg_geno_sig, vecsize(sample_ids) * sizeof(char));
		}
		else if (geno_augment_type == VAR_AUGMENT_ZERO)
		{
			memset(cur_augmenting_panel_var_reg_geno, 0, vecsize(sample_ids) * sizeof(char));
		}
		else
		{
			fprintf(stderr, "%s(%d): Genotype augmentation type %d is not understood, use %d=COPY or %d=ZERO\n", 
					__FILE__, __LINE__,
					geno_augment_type,
					VAR_AUGMENT_COPY, VAR_AUGMENT_ZERO);
			exit(1);
		}

		void** augment_var_reg_info = new void* [10];
		augment_var_reg_info[0] = cur_augmenting_panel_var_reg_geno;
		cur_augmenting_panel_var_reg->data = augment_var_reg_info;

		augmenting_panel_var_regs->push_back(cur_augmenting_panel_var_reg);
	} // i_int loop.

	vector<t_annot_region*>* pooled_geno_regs = new vector<t_annot_region*>();
	pooled_geno_regs->insert(pooled_geno_regs->end(), augmenting_panel_var_regs->begin(), augmenting_panel_var_regs->end());
	pooled_geno_regs->insert(pooled_geno_regs->end(), panel_var_regs->begin(), panel_var_regs->end());
	sort(pooled_geno_regs->begin(), pooled_geno_regs->end(), sort_regions);

	fprintf(stderr, "Saving %d augmented+%d original regions=%d total regions.\n", 
			vecsize(augmenting_panel_var_regs), vecsize(panel_var_regs), vecsize(pooled_geno_regs));

	binarize_variant_signal_regions_wrapper(pooled_geno_regs, sample_ids, op_matbed);
} // augment_tag_variants option.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 6/24: This is the main function for locally permuting/inverting (typed) variants. This function is similar to the recoding function but it does not use anchors.
// 6/24: This function should work when the panel contains more variants than the ones with permutations.
void proxize_variants_per_locality_permutation(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* permute_proxy_mappings_BED_fp,
	int n_threads,
	char* proxized_haplocoded_var_regs_op_fp)
{
	if (!check_file(permute_proxy_mappings_BED_fp))
	{
		fprintf(stderr, "Could not find the permutation mapping BED file @ %s\n", permute_proxy_mappings_BED_fp);
		exit(1);
	}

	vector<t_annot_region*>* perm_var_orig_regs = load_tag_permute_proxy_mapping_regs(permute_proxy_mappings_BED_fp);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and set the memory for the permuted variant signals.
	//size_t proxized_mem_pool_size = vecsize_t(perm_var_orig_regs) * vecsize_t(sample_ids) * sizeof(unsigned char);
	//fprintf(stderr, "Allocating the proxized genotype memory pool of %.3f Gb..\n", (double)proxized_mem_pool_size / (1024.0 * 1024.0 * 1024.0));
	//unsigned char* proxized_geno_sig_mem_pool = new unsigned char[proxized_mem_pool_size];
	//memset(proxized_geno_sig_mem_pool, 0, proxized_mem_pool_size);

	size_t n_ptrs_per_var_info = 10;
	void** focus_reg_info_mem_pool = new void* [n_ptrs_per_var_info * vecsize_t(perm_var_orig_regs)];
	memset(focus_reg_info_mem_pool, 0, n_ptrs_per_var_info * vecsize_t(perm_var_orig_regs) * sizeof(void*));

	vector<t_annot_region*>* focus_haplocoded_regions = new vector<t_annot_region*>();
	vector<t_annot_region*>* focus_haplocoded_unlinked_perm_regions = new vector<t_annot_region*>();

	// Intersect the panel's regions with the permuting regions. We have to have a complete match to the permuting regions and their connections.
	// Note that if this holds, we are good to go for the permuted regions.
	vector<t_annot_region*>* intersects = intersect_annot_regions(haplocoded_var_regs, perm_var_orig_regs, false, false);

	if (vecsize(intersects) != vecsize(perm_var_orig_regs))
	{
		fprintf(stderr, "Sanity check failed: The permuting variants list does not match the number of overlapping regions in the panel: %d-vs-%d\n",
			vecsize(intersects), vecsize(perm_var_orig_regs));
		exit(1);
	}

	for (int i_int = 0; i_int < vecsize(intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* var_reg = int_info->src_reg;
		t_annot_region* perm_desc_reg = int_info->dest_reg;

		if (perm_desc_reg->start == var_reg->start &&
			perm_desc_reg->end == var_reg->end)
		{
			// This is the flag for inverting permuted genotype.
			var_reg->score = perm_desc_reg->score;

			void** cur_var_reg_info = (void** )(var_reg->data);
			t_annot_region* unlinked_perm_reg = (t_annot_region*)(perm_desc_reg->data);

			size_t var_i_pool = focus_haplocoded_regions->size();
			void** new_var_reg_info = (void**)(focus_reg_info_mem_pool + var_i_pool * n_ptrs_per_var_info);
			new_var_reg_info[0] = cur_var_reg_info[0]; // Keep the genotypes.
			new_var_reg_info[1] = unlinked_perm_reg; // Find out where this region is on the haplocoded regions.
			new_var_reg_info[2] = NULL; // This is where the permuted region maps to in the actual list with genotype signals.
			new_var_reg_info[3] = NULL; // This is where the permuted region maps to in the actual list with genotype signals.
			var_reg->data = new_var_reg_info;

			// We do not know if this maps anywhere
			if (unlinked_perm_reg->data != NULL)
			{
				fprintf(stderr, "Sanity check failed @ %s(%d): Unlinked region was already assigned a region with genotype signal\n", __FILE__, __LINE__);
				exit(1);
			}
			unlinked_perm_reg->data = var_reg;

			focus_haplocoded_regions->push_back(var_reg);
			focus_haplocoded_unlinked_perm_regions->push_back(unlinked_perm_reg);
		}
	} // i_int loop.

	// Now map the pairs to original regions.
	vector<t_annot_region*>* pair_perm_intersects = intersect_annot_regions(focus_haplocoded_regions, focus_haplocoded_unlinked_perm_regions, false, false);

	for (int i_int = 0; i_int < vecsize(pair_perm_intersects); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(pair_perm_intersects->at(i_int)->data);

		t_annot_region* perm_reg_matching_geno_sig_reg = int_info->src_reg;
		t_annot_region* unlinked_perm_reg = int_info->dest_reg;

		// Does this region with genotype signal match to the permuted target region?
		if (perm_reg_matching_geno_sig_reg->start == unlinked_perm_reg->start &&
			perm_reg_matching_geno_sig_reg->end == unlinked_perm_reg->end)
		{
			//void** focus_reg_info = (void**)(focus_reg->data);

			// This is source region with genotype signal. This is what we are looking for. We now need to set the target of this variant as the matching geno. signal variant.
			t_annot_region* unlinked_var_haplocoded_reg = (t_annot_region*)(unlinked_perm_reg->data);
			void** unlinked_var_haplocoded_reg_info = (void**)(unlinked_var_haplocoded_reg->data);

			// This region must not have received a permuting region, yet.
			if (unlinked_var_haplocoded_reg_info[2] != NULL)
			{
				fprintf(stderr, "Sanity check failed @ %s(%d): We already assigned a permuting region with genotype signal\n", __FILE__, __LINE__);
				exit(1);
			}
			unlinked_var_haplocoded_reg_info[2] = perm_reg_matching_geno_sig_reg;
		}
	} // i_int loop.

	// Now, we assign the genotype signal regions; we could have done this above but this also checks to make sure all regions have an assigned region.
	char* permuted_geno_sig_mem_pool = new char[vecsize_t(focus_haplocoded_regions) * vecsize_t(sample_ids)];
	for (int i_reg = 0; i_reg < vecsize(focus_haplocoded_regions); i_reg++)
	{
		void** reg_info = (void**)(focus_haplocoded_regions->at(i_reg)->data);
		t_annot_region* perm_pair_geno_sig_reg = (t_annot_region*)(reg_info[2]);
		t_annot_region* unlinked_perm_pair_reg = (t_annot_region*)(reg_info[1]);

		// If there is no paired region, don't do anything; i.e., we do not permute this one.
		if (perm_pair_geno_sig_reg != NULL)
		{
			char* perm_geno_sig = (char*)(((void**)(perm_pair_geno_sig_reg->data))[0]);
			char* geno_sig = (permuted_geno_sig_mem_pool + i_reg * vecsize_t(sample_ids));
			memcpy(geno_sig, perm_geno_sig, vecsize_t(sample_ids));
			reg_info[3] = geno_sig;
		}
		else
		{
			// This is a problematic case, the paired region is not found among permuted regions; does that mean we will not correctly permute?
			fprintf(stderr, "Sanity check failed; there are non-matching connections: %s:%d-%d connects to %s:%d-%d, which does not exist in signal regions.\n",
					focus_haplocoded_regions->at(i_reg)->chrom, focus_haplocoded_regions->at(i_reg)->start, focus_haplocoded_regions->at(i_reg)->end,
					unlinked_perm_pair_reg->chrom, unlinked_perm_pair_reg->start, unlinked_perm_pair_reg->end);

			exit(1);
		}
	} // i_reg loop.

	// After everything is copied to permuted positions, replace the final genotype signals.
	fprintf(stderr, "Replacing, random inverting, genotype signals before saving...\n");
	for (int i_reg = 0; i_reg < vecsize(focus_haplocoded_regions); i_reg++)
	{
		// We need to copy this signal.
		void** reg_info = (void**)(focus_haplocoded_regions->at(i_reg)->data);
		char* geno_sig = (char*)(reg_info[3]);
		reg_info[0] = geno_sig;

		// Invert?
		if (focus_haplocoded_regions->at(i_reg)->score == 1)
		{
			//fprintf(stderr, "Inverting: %s:%d-%d\n", cur_chr_haplocoded_var_regs->at(i_var)->chrom,
			//	cur_chr_haplocoded_var_regs->at(i_var)->start, cur_chr_haplocoded_var_regs->at(i_var)->end);
			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				// This effectively reverse the signal:
				geno_sig[i_s] = 3 - geno_sig[i_s];
			}
		} // inversion check.		
	} // i_reg loop.

	fprintf(stderr, "Saving genotypes..\n");
	// Save the genotype signals.
	//binarize_variant_signal_regions_wrapper(haplocoded_var_regs, sample_ids, proxized_haplocoded_var_regs_op_fp);
	binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(haplocoded_var_regs, sample_ids,
		true,
		n_threads,
		proxized_haplocoded_var_regs_op_fp);
} // proxize_variants_per_locality_permutation function.

////////////////////////////////////////////////////////////////////////////////////////////////
// THIS IS THE MAIN FUNCTION FOR GENERATING LOCAL PERMUTATION MAPS FOR TAG VARIANTS.
// 6/24: Can we make permutation probabilistic? Even small permutations would protect variants very strongly; at places with low recomb rates.
void generate_permute_proxizing_parameters(char* geno_regions_BED_fp, int n_vicinity, double local_permute_prob, double geno_inversion_prob, char* op_BED_fp)
{
	fprintf(stderr, "%d-variant permuting variant positions:\n\
local_permute_prob: %.4f\n\
geno_inversion_prob: %.4f\n\
output_BED: %s\n", n_vicinity, local_permute_prob, geno_inversion_prob, op_BED_fp);

	//vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(haplocoded_geno_sig_regs_fp, sample_ids_list_fp);
	vector<t_annot_region*>* geno_sig_regs = load_BED(geno_regions_BED_fp);
	fprintf(stderr, "Loaded %d genotype regions.\n", (int)geno_sig_regs->size());

	vector<t_annot_region*>* permuted_regs = new vector<t_annot_region*>();
	permuted_regs->insert(permuted_regs->end(), geno_sig_regs->begin(), geno_sig_regs->end());

	fprintf(stderr, "%d-locality permuting %d regions.\n", n_vicinity, (int)geno_sig_regs->size());

	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	sort(geno_sig_regs->begin(), geno_sig_regs->end(), sort_regions);
	sort(permuted_regs->begin(), permuted_regs->end(), sort_regions);

	for (int i_reg = 0; i_reg < (int)permuted_regs->size(); i_reg++)
	{
		permuted_regs->at(i_reg)->score = i_reg;
	} // i_reg loop.

	for (int i_reg = 0; i_reg < (int)permuted_regs->size(); i_reg++)
	{
		int win_i_start = MAX(0, i_reg - n_vicinity);
		int win_i_end = MIN(i_reg + n_vicinity, (int)geno_sig_regs->size()-1);

		// Do probabilistic permutation only.
		if (rng->random_double_ran3() < local_permute_prob)
		{
			int n_vars_in_block = win_i_end - win_i_start + 1;
			vector<int>* perm_i = rng->fast_permute_indices(0, n_vars_in_block);

			// Copy the variants in permuted positions.
			for (int j_reg = 0; j_reg < n_vars_in_block; j_reg++)
			{
				if (i_reg % 1000 == 0)
				{
					fprintf(stderr, "Window[%d-%d]: %d\n", win_i_start, win_i_end, (int)perm_i->size());
					fprintf(stderr, "Permuting %d, %d (%d)\n",
						win_i_start + j_reg,
						win_i_start + perm_i->at(j_reg),
						(int)permuted_regs->size());
				}

				// Exchange the regions in permutation.
				t_annot_region* temp_reg = permuted_regs->at(win_i_start + j_reg);
				permuted_regs->at(win_i_start + j_reg) = permuted_regs->at(win_i_start + perm_i->at(j_reg));
				permuted_regs->at(win_i_start + perm_i->at(j_reg)) = temp_reg;
			} // j_reg loop.
		}
	} // i_reg loop.

	FILE* f_op = open_f(op_BED_fp, "w");
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		int invert_geno = 0;
		if (rng->random_double_ran3() < geno_inversion_prob)
		{
			invert_geno = 1;
		}

		fprintf(f_op, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%d\n", 
				geno_sig_regs->at(i_reg)->chrom,
				translate_coord(geno_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(geno_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				permuted_regs->at(i_reg)->chrom, 
				translate_coord(permuted_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(permuted_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				geno_sig_regs->at(i_reg)->score, permuted_regs->at(i_reg)->score,
				geno_sig_regs->at(i_reg)->dbl_score, permuted_regs->at(i_reg)->dbl_score,
				invert_geno);
	} // i_reg loop.
	close_f(f_op, op_BED_fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// THIS IS THE MAIN FUNCTION FOR GENERATING HASHING PARAMETERS.
// UPDATE 6/24::Input regions are simple variant regions, not genotype signal'ed regions. This function does (should) not rely on genotypes, otherwise it would leak information.
void generate_save_per_site_mixing_parameters_LD_aware(vector<t_annot_region*>* all_geno_sig_regs,
	t_rng* rng,
	int n_vicinity,
	double per_var_weight_prob,
	double weight_inversion_prob,
	double per_var2var_interaction_prob,
	double per_var2var2var_interaction_prob,
	int avg_geno_mod,
	double requested_normalized_N_e,
	int min_n_params_per_var,
	char* recombination_rate_dir,
	char* parameter_op_fp)
{
	fprintf(stderr, "Generating mixing parameters for %d signal regions over:\n\
%d bp vicinity\n\
%.4f first, %.4f second, %.4f third order weight probs\n\
%.4f weight inversion prob.\n\
Min # of weights per site: %d\n\
%d-modularity,\n\
Var-var interaction N_e: %.4f\n",
(int)all_geno_sig_regs->size(),
n_vicinity,
per_var_weight_prob,
per_var2var_interaction_prob,
per_var2var2var_interaction_prob,
weight_inversion_prob,
min_n_params_per_var,
avg_geno_mod,
requested_normalized_N_e);

	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(all_geno_sig_regs);

	int n_total_params = 0;

	FILE* f_param_vic_stats = open_f("per_var_vicinity_stats.bed", "w");
	for (int i_chr = 0; i_chr < (int)restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Setting parameters on chromosome %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_sig_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(1);
		}

		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", (int)cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(geno_sig_regs->at(i_reg), cur_chrom_recomb_regs);
			geno_sig_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		// Start assigning the weights: Within the vicinity, select a weight probabilistically if it is close in genotype distance.
		double init_normalized_N_e = requested_normalized_N_e;
		double delta_N_e_relaxation = requested_normalized_N_e / 1000.0;

		fprintf(stderr, "Requested N_e=%.4f\n\
Relaxation delta N_e=%.4f\n", requested_normalized_N_e, delta_N_e_relaxation);

		for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
		{
			void** new_reg_info = new void* [10];
			geno_sig_regs->at(i_reg)->data = new_reg_info;

			// Initialize the mixing weights.
			int* per_var_weight = new int[2 * n_vicinity + 1];

			int** per_var2var_interaction_weight = new int* [2 * n_vicinity + 1];
			memset(per_var2var_interaction_weight, 0, sizeof(int*) * (2 * n_vicinity + 1));

			int*** per_var2var2var_interaction_weight = new int** [2 * n_vicinity + 1];
			memset(per_var2var2var_interaction_weight, 0, sizeof(int**) * (2 * n_vicinity + 1));

			int max_vic_i = -1;
			int max_vic_j = -1;
			int max_vic_k = -1;

			int n_first_order_weights = 0;
			int n_second_order_weights = 0;
			int n_third_order_weights = 0;

			double normalized_N_e = init_normalized_N_e;

			// This check makes sure we selected at least one vicinity weight for this variant.			
			// This step checks if we have enough weights for this position.
			while ((n_first_order_weights + n_second_order_weights + n_third_order_weights) < min_n_params_per_var)
			{
				int n_cur_total_weights = (n_first_order_weights + n_second_order_weights + n_third_order_weights);

				if (__DUMP_PROXIZATION_MSGS__)
				{
					fprintf(stderr, "Generating weights @ %s:%d: %d weights [Normalized N_e=%.4f]               \r",
						geno_sig_regs->at(i_reg)->chrom, geno_sig_regs->at(i_reg)->start, n_cur_total_weights, normalized_N_e);
				}

				// Reset the number of parameters since we are reselecting them.
				n_total_params = 0;
				n_first_order_weights = 0;
				n_second_order_weights = 0;
				n_third_order_weights = 0;

				// First, reset the first order weights if we allocated it before.
				if (per_var_weight != NULL)
				{
					delete[] per_var_weight;
					per_var_weight = new int[2 * n_vicinity + 1];
				}
				memset(per_var_weight, 0, sizeof(int) * (2 * n_vicinity + 1));

				// Go over each first order vicinity and set weights.
				for (int vic_i = 0; vic_i < (2 * n_vicinity + 1); vic_i++)
				{
					// This is the index of the proximal variant.
					int vic_i_i_reg = MAX(0, i_reg - n_vicinity + vic_i);
					vic_i_i_reg = MIN(vic_i_i_reg, (int)geno_sig_regs->size() - 1);
					double i_vic_i_genetic_dist = fabs(geno_sig_regs->at(i_reg)->dbl_score - geno_sig_regs->at(vic_i_i_reg)->dbl_score);

					// Should we use this vicinity variant? We select these probabilistically to introduce randomization.
					double r_m = fabs(i_vic_i_genetic_dist);
					double rho_m = 4 * normalized_N_e * r_m;

					double n_ref_haplotypes = 1.0;
					double tau_m_vic_i = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					if (__DUMP_PROXIZATION_MSGS__)
					{
						fprintf(stderr, "Var @ %d::%s:%d ;; Vic_i Variant @ %d [%d]: %s:%d;; gen-dist: %.4f; recomb-prob: %.4f\n",
							i_reg, geno_sig_regs->at(i_reg)->chrom, geno_sig_regs->at(i_reg)->start,
							vic_i_i_reg, vic_i, geno_sig_regs->at(vic_i_i_reg)->chrom, geno_sig_regs->at(vic_i_i_reg)->start,
							i_vic_i_genetic_dist, tau_m_vic_i);
					}

					// Tau_m is the recombination probability to other haplotypes. If it is high, we don't want to use it.
					// Following randomly selects small recombination probabilities among vicinity.
					// normalized_N_e tunes this option for selecting variants; small N_e always selects large vicinity variants. High N_e only selects a small vicinity subset.
					bool use_this_vic_var = (rng->random_double_ran3() > tau_m_vic_i);
					if (use_this_vic_var)
					{
						if (rng->random_double_ran3() < per_var_weight_prob)
						{
							//fprintf(stderr, "Single %d\n", i_var);
							max_vic_i = (int)(MAX(max_vic_i, fabs((double)(vic_i - n_vicinity))));
							n_first_order_weights++;

							per_var_weight[vic_i] = 1;
							if (rng->random_double_ran3() < weight_inversion_prob)
							{
								per_var_weight[vic_i] = -1;
							}
						}
					}

					// Update the total number of parameters that are evaluated.
					n_total_params++;

					// Go over all interaction terms for vic_i variant.
					if (per_var2var_interaction_weight[vic_i] != NULL)
					{
						delete[] per_var2var_interaction_weight[vic_i];
					}
					per_var2var_interaction_weight[vic_i] = new int[2 * n_vicinity + 1];
					memset(per_var2var_interaction_weight[vic_i], 0, sizeof(int) * (2 * n_vicinity + 1));

					if (per_var2var2var_interaction_weight[vic_i] != NULL)
					{
						delete[] per_var2var2var_interaction_weight[vic_i];
					}
					per_var2var2var_interaction_weight[vic_i] = new int* [2 * n_vicinity + 1];
					memset(per_var2var2var_interaction_weight[vic_i], 0, sizeof(int*) * (2 * n_vicinity + 1));

					for (int vic_j = 0; vic_j < vic_i; vic_j++)
					{
						int vic_j_i_reg = MAX(0, i_reg - n_vicinity + vic_j);
						vic_j_i_reg = MIN(vic_j_i_reg, (int)geno_sig_regs->size() - 1);
						double i_vic_j_genetic_dist = fabs(geno_sig_regs->at(i_reg)->dbl_score - geno_sig_regs->at(vic_j_i_reg)->dbl_score);
						r_m = fabs(i_vic_j_genetic_dist);
						rho_m = 4 * normalized_N_e * r_m;
						double tau_m_vic_j = 1 - exp(-1 * rho_m / n_ref_haplotypes);

						bool use_this_vic2vic_int = (rng->random_double_ran3() > tau_m_vic_i) && (rng->random_double_ran3() > tau_m_vic_j);

						// Update the total number of parameters that are evaluated.
						n_total_params++;

						if (use_this_vic2vic_int)
						{
							if (rng->random_double_ran3() < per_var2var_interaction_prob)
							{
								max_vic_i = (int)(MAX(max_vic_i, fabs((double)(vic_i - n_vicinity))));
								max_vic_j = (int)(MAX(max_vic_j, fabs((double)(vic_j - n_vicinity))));
								n_second_order_weights++;

								//fprintf(stderr, "Pairwise %d-%d\n", i_var, j_var);
								per_var2var_interaction_weight[vic_i][vic_j] = 1;

								if (rng->random_double_ran3() < weight_inversion_prob)
								{
									per_var2var_interaction_weight[vic_i][vic_j] = -1;
								}
							}
						}

						// Allocate the 3rd dim. interaction terms.
						if (per_var2var2var_interaction_weight[vic_i][vic_j] != NULL)
						{
							delete[] per_var2var2var_interaction_weight[vic_i][vic_j];
						}
						per_var2var2var_interaction_weight[vic_i][vic_j] = new int[2 * n_vicinity + 1];
						memset(per_var2var2var_interaction_weight[vic_i][vic_j], 0, sizeof(int) * (2 * n_vicinity + 1));

						// This is the third variant, we only allocate and assign these for vic_i>vic_j>vic_k.
						for (int vic_k = 0; vic_k < vic_j; vic_k++)
						{
							int vic_k_i_reg = MAX(0, i_reg - n_vicinity + vic_k);
							vic_k_i_reg = MIN(vic_k_i_reg, (int)geno_sig_regs->size() - 1);
							double i_vic_k_genetic_dist = fabs(geno_sig_regs->at(i_reg)->dbl_score - geno_sig_regs->at(vic_k_i_reg)->dbl_score);
							r_m = fabs(i_vic_k_genetic_dist);
							rho_m = 4 * normalized_N_e * r_m;
							double tau_m_vic_k = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							bool use_this_vic2vic2vic_int = (rng->random_double_ran3() > tau_m_vic_i) &&
								(rng->random_double_ran3() > tau_m_vic_j) &&
								(rng->random_double_ran3() > tau_m_vic_k);

							// Update the total number of parameters that are evaluated.
							n_total_params++;

							if (use_this_vic2vic2vic_int)
							{
								if (rng->random_double_ran3() < per_var2var2var_interaction_prob)
								{
									max_vic_i = (int)(MAX(max_vic_i, fabs((double)(vic_i - n_vicinity))));
									max_vic_j = (int)(MAX(max_vic_j, fabs((double)(vic_j - n_vicinity))));
									max_vic_k = (int)(MAX(max_vic_k, fabs((double)(vic_k - n_vicinity))));
									n_third_order_weights++;

									//fprintf(stderr, "Pairwise %d-%d\n", i_var, j_var);
									per_var2var2var_interaction_weight[vic_i][vic_j][vic_k] = 1;

									if (rng->random_double_ran3() < weight_inversion_prob)
									{
										per_var2var2var_interaction_weight[vic_i][vic_j][vic_k] = -1;
									}
								}
							}
						}
					} // vic_j loop.
				} // vic_i loop.

				// Check if we have enough variants for proxizing this variant.
				if ((n_first_order_weights + n_second_order_weights + n_third_order_weights) < min_n_params_per_var)
				{
					normalized_N_e -= delta_N_e_relaxation;
					normalized_N_e = MAX(0, normalized_N_e);
				}
			} // parameter selection loop.

			fprintf(f_param_vic_stats, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n",
				geno_sig_regs->at(i_reg)->chrom, geno_sig_regs->at(i_reg)->start, geno_sig_regs->at(i_reg)->end,
				max_vic_i, max_vic_j, max_vic_k,
				n_first_order_weights, n_second_order_weights, n_third_order_weights,
				normalized_N_e); // Report the normalized_N_e that is used to select this variant's weights.

			if (i_reg == 0)
			{
				fprintf(stderr, "%d total proxy-weights per variant.\n", n_total_params);
			}

			if (__DUMP_PROXIZATION_MSGS__)
			{
				if (i_reg < 10)
				{
					fprintf(stderr, "%s:%d-%d generated parameters:\n",
						geno_sig_regs->at(i_reg)->chrom,
						geno_sig_regs->at(i_reg)->start,
						geno_sig_regs->at(i_reg)->end);

					fprintf(stderr, "Per variant weights:\n");
					for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
					{
						fprintf(stderr, "%d=>%d ;; ", vic_i, per_var_weight[vic_i]);
					} // vic_i loop.

					fprintf(stderr, "\n");

					fprintf(stderr, "Var-2-var interaction weights:\n");
					for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
					{
						for (int vic_j = 0; vic_j < vic_i; vic_j++)
						{
							fprintf(stderr, "%dx%d => %d ;; ", vic_i, vic_j, per_var2var_interaction_weight[vic_i][vic_j]);
						} // vic_j loop.
					} // vic_i loop.

					fprintf(stderr, "\n\n");
				} // i_reg check.
			} // proxy message dump check.

			// Save the parameters.
			void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
			cur_reg_info[0] = per_var_weight;
			cur_reg_info[1] = per_var2var_interaction_weight;
			cur_reg_info[2] = per_var2var2var_interaction_weight;
		} // i_reg loop.
	} // i_chr loop.

	close_f(f_param_vic_stats, NULL);

	fprintf(stderr, "%d total parameters per variant.\n", n_total_params);

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Save the parameters: This saves the final file, below is sanity check for the loading function.
	save_per_variant_site_mixing_parameters(all_geno_sig_regs, n_vicinity, avg_geno_mod, parameter_op_fp);

	// Load/check.
	// Now load and test.
	fprintf(stderr, "Loading mixing parameters per variant..\n");
	int loaded_n_vicinity = 0;
	int loaded_avg_geno_mod = 0;
	vector<t_annot_region*>* loaded_per_var_proxy_params = load_per_variant_site_mixing_parameters(parameter_op_fp, loaded_n_vicinity, loaded_avg_geno_mod);

	// Sanity check:
	fprintf(stderr, "Sanity checking loaded parameters.. ");
	for (int i_var = 0; i_var < (int)loaded_per_var_proxy_params->size(); i_var++)
	{
		void** cur_var_loaded_reg_info = (void**)(loaded_per_var_proxy_params->at(i_var)->data);
		int* loaded_per_var_weights = (int*)(cur_var_loaded_reg_info[0]);
		int** loaded_per_var2var_weights = (int**)(cur_var_loaded_reg_info[1]);
		int*** loaded_per_var2var2var_weights = (int***)(cur_var_loaded_reg_info[2]);

		void** cur_var_reg_info = (void**)(all_geno_sig_regs->at(i_var)->data);
		int* per_var_weights = (int*)(cur_var_reg_info[0]);
		int** per_var2var_weights = (int**)(cur_var_reg_info[1]);
		int*** per_var2var2var_weights = (int***)(cur_var_reg_info[2]);

		for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
		{
			if (per_var_weights[vic_i] != loaded_per_var_weights[vic_i])
			{
				fprintf(stderr, "Sanity check failed: %s:%d-%d: Var %d; %d/%d\n",
					loaded_per_var_proxy_params->at(i_var)->chrom,
					loaded_per_var_proxy_params->at(i_var)->start,
					loaded_per_var_proxy_params->at(i_var)->end,
					vic_i,
					per_var_weights[vic_i], loaded_per_var_weights[vic_i]);

				exit(1);
			}

			//for (int vic_j = 0; vic_j < 2 * n_vicinity + 1; vic_j++)
			for (int vic_j = 0; vic_j < vic_i; vic_j++)
			{
				if (per_var2var_weights[vic_i][vic_j] != loaded_per_var2var_weights[vic_i][vic_j])
				{
					fprintf(stderr, "Sanity check failed: %s:%d-%d: Var2var[%dx%d]; %d/%d\n",
						loaded_per_var_proxy_params->at(i_var)->chrom,
						loaded_per_var_proxy_params->at(i_var)->start,
						loaded_per_var_proxy_params->at(i_var)->end,
						vic_i, vic_j,
						per_var2var_weights[vic_i][vic_j], loaded_per_var2var_weights[vic_i][vic_j]);

					exit(1);
				}

				for (int vic_k = 0; vic_k < 2 * n_vicinity + 1; vic_k++)
				{
					if (per_var2var2var_weights[vic_i][vic_j][vic_k] != loaded_per_var2var2var_weights[vic_i][vic_j][vic_k])
					{
						fprintf(stderr, "Sanity check failed: %s:%d-%d: Var2Var2Var[%dx%dx%d]; %d/%d\n",
							loaded_per_var_proxy_params->at(i_var)->chrom,
							loaded_per_var_proxy_params->at(i_var)->start,
							loaded_per_var_proxy_params->at(i_var)->end,
							vic_i, vic_j, vic_k,
							per_var2var2var_weights[vic_i][vic_j][vic_k], loaded_per_var2var2var_weights[vic_i][vic_j][vic_k]);

						exit(1);
					}
				} // vic_k loop.
			} // vic_j loop.
		} // vic_i loop.
	} // i_var loop.

	fprintf(stderr, "Successfully generated the mixing (hashing/proxizing...) parameters for the variant set..\n");
} // generate_save_per_site_mixing_parameters function.

void generate_save_per_site_mixing_parameters(vector<t_annot_region*>* geno_sig_regs, 
	t_rng* rng, 
	int n_vicinity,
	double per_var_weight_prob,
	double per_var2var_interaction_prob,
	int avg_geno_mod,
	char* parameter_op_fp)
{
	fprintf(stderr, "Generating mixing parameters for %d signal regions over %d bp vicinity using %.4f first, %.4f second order weight probs, and %d-modularity.\n",
		(int)geno_sig_regs->size(), n_vicinity, per_var_weight_prob, per_var2var_interaction_prob, avg_geno_mod);

	int n_total_params = 0;
	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		n_total_params = 0;

		// Initialize the mixing weights.
		int* per_var_weight = new int[2 * n_vicinity + 1];
		memset(per_var_weight, 0, sizeof(int) * (2 * n_vicinity + 1));

		int** per_var2var_interaction_weight = new int* [2 * n_vicinity + 1];

		// Params file doesnt exist, initialize the weights.
		//fprintf(stderr, "Initializing parameters from scratch.\n");
		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			if (rng->random_double_ran3() < per_var_weight_prob)
			{
				//fprintf(stderr, "Single %d\n", i_var);
				per_var_weight[i_var] = 1;
			}

			n_total_params++;

			per_var2var_interaction_weight[i_var] = new int[2 * n_vicinity + 1];
			memset(per_var2var_interaction_weight[i_var], 0, sizeof(int) * (2 * n_vicinity + 1));

			for (int j_var = 0; j_var < i_var; j_var++)
			{
				n_total_params++;
				if (rng->random_double_ran3() < per_var2var_interaction_prob)
				{
					//fprintf(stderr, "Pairwise %d-%d\n", i_var, j_var);
					per_var2var_interaction_weight[i_var][j_var] = 1;
				}
			} // j_var loop.
		} // i_var loop.

		if (i_reg == 0)
		{
			fprintf(stderr, "%d total parameters.\n", n_total_params);
		}

		if (__DUMP_PROXIZATION_MSGS__)
		{
			if (i_reg < 10)
			{
				fprintf(stderr, "%s:%d-%d generated parameters:\n",
						geno_sig_regs->at(i_reg)->chrom,
						geno_sig_regs->at(i_reg)->start,
						geno_sig_regs->at(i_reg)->end);

				fprintf(stderr, "Per variant weights:\n");
				for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
				{
					fprintf(stderr, "%d=>%d ;; ", vic_i, per_var_weight[vic_i]);
				} // vic_i loop.

				fprintf(stderr, "\n");

				fprintf(stderr, "Var-2-var interaction weights:\n");
				for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
				{
					for (int vic_j = 0; vic_j < vic_i; vic_j++)
					{
						fprintf(stderr, "%dx%d => %d ;; ", vic_i, vic_j, per_var2var_interaction_weight[vic_i][vic_j]);
					} // vic_j loop.
				} // vic_i loop.

				fprintf(stderr, "\n\n");
			} // i_reg check.
		} // proxy message dump check.

		// Save the parameters.
		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		cur_reg_info[2] = per_var_weight;
		cur_reg_info[3] = per_var2var_interaction_weight;
	} // i_reg loop.

	fprintf(stderr, "%d total parameters per variant.\n", n_total_params);

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Save the parameters: This saves the final file, below is sanity check for the loading function.
	save_per_variant_site_mixing_parameters(geno_sig_regs, n_vicinity, avg_geno_mod, parameter_op_fp);

	// Load/check.
	// Now load and test.
	fprintf(stderr, "Loading mixing parameters per variant..\n");
	int loaded_n_vicinity = 0;
	int loaded_avg_geno_mod = 0;
	vector<t_annot_region*>* loaded_per_var_proxy_params = load_per_variant_site_mixing_parameters(parameter_op_fp, loaded_n_vicinity, loaded_avg_geno_mod);

	// Sanity check:
	fprintf(stderr, "Sanity checking loaded parameters.. ");
	for (int i_var = 0; i_var < (int)loaded_per_var_proxy_params->size(); i_var++)
	{
		void** cur_var_loaded_reg_info = (void**)(loaded_per_var_proxy_params->at(i_var)->data);
		int* loaded_per_var_weights = (int*)(cur_var_loaded_reg_info[2]);
		int** loaded_per_var2var_weights = (int**)(cur_var_loaded_reg_info[3]);

		void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_var)->data);
		int* per_var_weights = (int*)(cur_var_reg_info[2]);
		int** per_var2var_weights = (int**)(cur_var_reg_info[3]);

		for (int vic_i = 0; vic_i < 2 * n_vicinity + 1; vic_i++)
		{
			if (per_var_weights[vic_i] != loaded_per_var_weights[vic_i])
			{
				fprintf(stderr, "Sanity check failed: %s:%d-%d: Var %d; %d/%d\n",
					loaded_per_var_proxy_params->at(i_var)->chrom,
					loaded_per_var_proxy_params->at(i_var)->start,
					loaded_per_var_proxy_params->at(i_var)->end,
					vic_i,
					per_var_weights[vic_i], loaded_per_var_weights[vic_i]);

				exit(1);
			}

			for (int vic_j = 0; vic_j < 2 * n_vicinity + 1; vic_j++)
			{
				if (per_var2var_weights[vic_i][vic_j] != loaded_per_var2var_weights[vic_i][vic_j])
				{
					fprintf(stderr, "Sanity check failed: %s:%d-%d: Var2var[%dx%d]; %d/%d\n",
						loaded_per_var_proxy_params->at(i_var)->chrom,
						loaded_per_var_proxy_params->at(i_var)->start,
						loaded_per_var_proxy_params->at(i_var)->end,
						vic_i, vic_j,
						per_var2var_weights[vic_i][vic_j], loaded_per_var2var_weights[vic_i][vic_j]);

					exit(1);
				}
			} // vic_j loop.
		} // vic_i loop.
	} // i_var loop.

	fprintf(stderr, "Success..\n");
} // generate_save_per_site_mixing_parameters function.

void save_per_variant_site_mixing_parameters(vector<t_annot_region*>* geno_sig_regs, int n_vicinity, int avg_geno_mod, char* parameter_op_fp)
{
	////////////////////////////////////////////////////////////////////////////////////////////////
// Save the parameters.
	fprintf(stderr, "Saving per-variant parameters to %s\n", parameter_op_fp);
	FILE* f_weight_params = open_f(parameter_op_fp, "wb");

	fwrite(&(n_vicinity), sizeof(int), 1, f_weight_params);
	fwrite(&(avg_geno_mod), sizeof(int), 1, f_weight_params);

	// Write the number of variants.
	int n_vars = (int)((int)geno_sig_regs->size());
	fwrite(&(n_vars), sizeof(int), 1, f_weight_params);

	for (int i_reg = 0; i_reg < (int)geno_sig_regs->size(); i_reg++)
	{
		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		int* per_var_weight = (int*)(cur_reg_info[0]);
		int** per_var2var_interaction_weight = (int**)(cur_reg_info[1]);
		int*** per_var2var2var_interaction_weight = (int***)(cur_reg_info[2]);

		int l_chrom_str = t_string::string_length(geno_sig_regs->at(i_reg)->chrom);
		fwrite(&(l_chrom_str), sizeof(int), 1, f_weight_params);
		fwrite(geno_sig_regs->at(i_reg)->chrom, sizeof(char), l_chrom_str, f_weight_params);

		fwrite(&(geno_sig_regs->at(i_reg)->start), sizeof(int), 1, f_weight_params);
		fwrite(&(geno_sig_regs->at(i_reg)->end), sizeof(int), 1, f_weight_params);

		// Write params
		fwrite(per_var_weight, sizeof(int), 2 * n_vicinity + 1, f_weight_params);

		// Write the pw interaction terms.
		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			fwrite(per_var2var_interaction_weight[i_var], sizeof(int), (2 * n_vicinity + 1), f_weight_params);
		} // i_vic loop.

		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			for (int j_var = 0; j_var < i_var; j_var++)
			{
				fwrite(per_var2var2var_interaction_weight[i_var][j_var], sizeof(int), (2 * n_vicinity + 1), f_weight_params);
			}
		} // i_vic loop.
	} // i_reg loop.
	close_f(f_weight_params, parameter_op_fp);
}

vector<t_annot_region*>* load_per_variant_site_mixing_parameters(char* parameters_fp,
	int& loaded_n_vicinity, int& loaded_avg_geno_mod)
{
	if (!check_file(parameters_fp))
	{
		return(NULL);
	}

	// Save the parameters.
	fprintf(stderr, "Reading per-variant parameters to %s\n", parameters_fp);
	FILE* f_weight_params = open_f(parameters_fp, "rb");

	vector<t_annot_region*>* geno_reg_w_params = new vector<t_annot_region*>();

	// Write the number of variants.
	int n_vicinity = 0;
	int avg_geno_mod = 0;
	//fread(&(n_vicinity), sizeof(int), 1, f_weight_params);
	n_vicinity = read_bin_int(f_weight_params);
	//fread(&(avg_geno_mod), sizeof(int), 1, f_weight_params);
	avg_geno_mod = read_bin_int(f_weight_params);

	// Set the loaded parameters.
	loaded_n_vicinity = n_vicinity;
	loaded_avg_geno_mod = avg_geno_mod;

	fprintf(stderr, "# vicinity vars: %d; geno_mod: %d\n", n_vicinity, avg_geno_mod);

	int n_vars = 0;
	//fread(&n_vars, sizeof(int), 1, f_weight_params);
	n_vars = read_bin_int(f_weight_params);

	fprintf(stderr, "Loading parameters for %d variants.\n", n_vars);

	for (int i_reg = 0; i_reg < n_vars; i_reg++)
	{
		t_annot_region* geno_reg = get_empty_region();
		void** cur_reg_info = new void*[10];

		// Allocate the mixing parameters.
		int* per_var_weight = new int[2 * n_vicinity + 1];
		memset(per_var_weight, 0, sizeof(int) * (2 * n_vicinity + 1));
		int** per_var2var_interaction_weight = new int* [2 * n_vicinity + 1];
		int*** per_var2var2var_interaction_weight = new int** [2 * n_vicinity + 1];
		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			per_var2var_interaction_weight[i_var] = new int[2 * n_vicinity + 1];
			memset(per_var2var_interaction_weight[i_var], 0, sizeof(int) * (2 * n_vicinity + 1));

			per_var2var2var_interaction_weight[i_var] = new int* [2 * n_vicinity + 1];

			//for (int j_var = 0; j_var < (2 * n_vicinity + 1); i_var++)
			for (int j_var = 0; j_var < i_var; j_var++)
			{
				per_var2var2var_interaction_weight[i_var][j_var] = new int [2 * n_vicinity + 1];
				memset(per_var2var2var_interaction_weight[i_var][j_var], 0, sizeof(int) * (2 * n_vicinity + 1));
			}
		} // j_reg loop.

		// Read the chromosome id.
		int l_chrom_str = 0;
		//fread(&(l_chrom_str), sizeof(int), 1, f_weight_params);
		l_chrom_str = read_bin_int(f_weight_params);
		char cur_chrom[100];
		memset(cur_chrom, 0, sizeof(char) * 100);
		//fread(cur_chrom, sizeof(char), l_chrom_str, f_weight_params);
		read_bin_char_array(f_weight_params, cur_chrom, l_chrom_str);

		int start_pos = 0;
		int end_pos = 0;
		//fread(&(start_pos), sizeof(int), 1, f_weight_params);
		start_pos = read_bin_int(f_weight_params);
		//fread(&(end_pos), sizeof(int), 1, f_weight_params);
		end_pos = read_bin_int(f_weight_params);

		geno_reg->chrom = t_string::copy_me_str(cur_chrom);
		geno_reg->start = start_pos;
		geno_reg->end = end_pos;
		geno_reg->strand = '+';
		geno_reg->data = cur_reg_info;

		cur_reg_info[0] = per_var_weight;
		cur_reg_info[1] = per_var2var_interaction_weight;
		cur_reg_info[2] = per_var2var2var_interaction_weight;

		// Write params
		//fread(per_var_weight, sizeof(int), 2 * n_vicinity + 1, f_weight_params);
		read_bin_int_array(f_weight_params, per_var_weight, 2 * n_vicinity + 1);

		// Read the pw interaction terms.
		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			//fread(per_var2var_interaction_weight[i_var], sizeof(int), (2 * n_vicinity + 1), f_weight_params);
			read_bin_int_array(f_weight_params, per_var2var_interaction_weight[i_var], (2 * n_vicinity + 1));
		} // i_vic loop.

		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			//for (int j_var = 0; j_var < (2 * n_vicinity + 1); j_var++)
			for (int j_var = 0; j_var < i_var; j_var++)
			{
				//fread(per_var2var2var_interaction_weight[i_var][j_var], sizeof(int), (2 * n_vicinity + 1), f_weight_params);
				read_bin_int_array(f_weight_params, per_var2var2var_interaction_weight[i_var][j_var], (2 * n_vicinity + 1));
			} // j_var loop.
		} // i_var loop.

		geno_reg_w_params->push_back(geno_reg);
	} // i_reg loop.
	close_f(f_weight_params, parameters_fp);

	if (__DUMP_PROXIZATION_MSGS__)
	{
		for (int i_var = 0; i_var < 10; i_var++)
		{
			fprintf(stderr, "%s:%d-%d loaded parameters:\n",
				geno_reg_w_params->at(i_var)->chrom,
				geno_reg_w_params->at(i_var)->start,
				geno_reg_w_params->at(i_var)->end);

			void** reg_info = (void**)(geno_reg_w_params->at(i_var)->data);
			int* per_var_weight = (int*)(reg_info[2]);
			int** per_var2var_interaction_weight = (int**)(reg_info[3]);

			fprintf(stderr, "Per variant weights:\n");
			for (int vic_i = 0; vic_i < 2 * loaded_n_vicinity + 1; vic_i++)
			{
				fprintf(stderr, "%d=>%d ;; ", vic_i, per_var_weight[vic_i]);
			} // vic_i loop.

			fprintf(stderr, "\n");

			fprintf(stderr, "Var-2-var interaction weights:\n");
			for (int vic_i = 0; vic_i < 2 * loaded_n_vicinity + 1; vic_i++)
			{
				for (int vic_j = 0; vic_j < vic_i; vic_j++)
				{
					fprintf(stderr, "%dx%d => %d ;; ", vic_i, vic_j, per_var2var_interaction_weight[vic_i][vic_j]);
				} // vic_j loop.
			} // vic_i loop.

			// Add extra line.
			fprintf(stderr, "\n\n");
		} // i_var loop.
	} // dump check.

	// Return the list of regions.
	return(geno_reg_w_params);
} // load_per_variant_site_mixing_parameters function.

// This function tests the proxization of variants (genotypes) using modular average of surrounding variants with interacting terms, i.e., AND operations.
void proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	char* per_var_vicinity_weight_params_fp,
	double allele_err_eps,
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
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
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
		var_reg_info[0] = vic_params_reg_info[0];
		var_reg_info[1] = vic_params_reg_info[1];
		var_reg_info[2] = vic_params_reg_info[2];
	} // i_int loop.
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start processing the proxized genotypes.
	for (int i_chr = 0; i_chr < (int)restr_haplocoded_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_haplocoded_var_regs = restr_haplocoded_var_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Sorting the variants per position..\n");

		sort(cur_chr_haplocoded_var_regs->begin(), cur_chr_haplocoded_var_regs->end(), sort_regions);

		// Process each sample.
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			// Process each haplotype
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				for (int i_var = 0; i_var < (int)cur_chr_haplocoded_var_regs->size(); i_var++)
				{
					// Store the var_i info.
					void** var_i_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(i_var)->data);
					//char* var_i_geno_sig = var_i_reg_info[0];
					char* var_i_proxized_geno_sig = (char*)(var_i_reg_info[1]);

					// Select the proxization parameters for this variant.
					int* per_var_weight = (int*)(var_i_reg_info[0]);
					int** per_var2var_interaction_weight = (int**)(var_i_reg_info[1]);
					int*** per_var2var2var_interaction_weight = (int***)(var_i_reg_info[2]);

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

// This function tests the proxization of variants (genotypes) using modular average of surrounding variants with interacting terms, i.e., AND operations.
void proxize_variants_per_vicinity_non_linear_modular_average_uniform(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
	int n_vicinity,
	char* weight_params_fp,
	double per_var_weight_prob,
	double per_var2var_interaction_prob,	
	int avg_geno_mod,
	double allele_err_eps,
	char* proxized_haplocoded_var_regs_op_fp)
{
	fprintf(stderr, "Proxizing %s (%s) per non-linear modular combinations using:\n\\sum_{vic @ [-%d, +%d]) (w_i x g_i)+\\sum_{vic @ [-%d, +%d] w_ij x g_i x g_j mod %d\n\
			prob(w_i)=%.4f; prob(w_ij)=%.4f\n",
			haplocoded_geno_signal_fp, sample_ids_list_fp,
			n_vicinity, n_vicinity, n_vicinity, n_vicinity, avg_geno_mod, per_var_weight_prob, per_var2var_interaction_prob);

	// Instantiate the RNG.
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load and validate the genotypes.
	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
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

		// Replace the variant information.
		haplocoded_var_regs->at(i_var)->data = new_reg_info;
		haplocoded_var_regs->at(i_var)->score = 0;
	} // i_var loop.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Initialize the mixing weights.
	int* per_var_weight = new int[2 * n_vicinity + 1];
	memset(per_var_weight, 0, sizeof(int) * (2 * n_vicinity + 1));

	int** per_var2var_interaction_weight = new int* [2 * n_vicinity + 1];

	// Initialize the weights from a file?
	if (check_file(weight_params_fp))
	{
		fprintf(stderr, "Loading parameters from %s\n", weight_params_fp);

		// Load the parameters:
		FILE* f_weight_params = open_f(weight_params_fp, "rb");		
		//fread(per_var_weight, sizeof(int), 2 * n_vicinity + 1, f_weight_params);
		read_bin_int_array(f_weight_params, per_var_weight, (2 * n_vicinity + 1));

		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			per_var2var_interaction_weight[i_var] = new int[2 * n_vicinity + 1];
			//memset(per_var2var_interaction_weight[i_var], 0, sizeof(int) * (2 * n_vicinity + 1));
			//fread(per_var2var_interaction_weight[i_var], sizeof(int), 2 * n_vicinity + 1, f_weight_params);
			read_bin_int_array(f_weight_params, per_var2var_interaction_weight[i_var], (2 * n_vicinity + 1));
		} // i_vic loop.

		close_f(f_weight_params, weight_params_fp);
	} // Weight loading check.
	else
	{
		// Params file doesnt exist, initialize the weights.
		fprintf(stderr, "Initializing parameters from scratch.\n");
		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			if (rng->random_double_ran3() < per_var_weight_prob)
			{
				fprintf(stderr, "Single %d\n", i_var);
				per_var_weight[i_var] = 1;
			}

			per_var2var_interaction_weight[i_var] = new int[2 * n_vicinity + 1];
			memset(per_var2var_interaction_weight[i_var], 0, sizeof(int) * (2 * n_vicinity + 1));

			for (int j_var = 0; j_var < i_var; j_var++)
			{
				if (rng->random_double_ran3() < per_var2var_interaction_prob)
				{
					fprintf(stderr, "Pairwise %d-%d\n", i_var, j_var);
					per_var2var_interaction_weight[i_var][j_var] = 1;
				}
			} // j_var loop.
		} // i_var loop.

		// Save the parameters:
		fprintf(stderr, "Saving parameters to %s\n", weight_params_fp);
		FILE* f_weight_params = open_f(weight_params_fp, "wb");
		fwrite(per_var_weight, sizeof(int), 2 * n_vicinity + 1, f_weight_params);

		for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
		{
			fwrite(per_var2var_interaction_weight[i_var], sizeof(int), (2 * n_vicinity + 1), f_weight_params);
		} // i_vic loop.

		close_f(f_weight_params, weight_params_fp);
	} // New Weight Initialization check.

	///////////////////////////////////////////////////////////
	// Report the weight statistics:
	int n_single_terms = 0;
	int n_pw_terms = 0;
	int n_total_params = 0;
	for (int i_var = 0; i_var < (2 * n_vicinity + 1); i_var++)
	{
		if (per_var_weight[i_var] == 1)
		{
			fprintf(stderr, "Weight %d = 1\n", i_var);
			n_single_terms++;
		}

		n_total_params++;

		for (int j_var = 0; j_var < i_var; j_var++)
		{
			if (per_var2var_interaction_weight[i_var][j_var] == 1)
			{
				fprintf(stderr, "Pairwise @ %d-%d = 1\n", i_var, j_var);
				n_pw_terms++;
			}

			n_total_params++;
		} // j_var loop.
	} // i_var loop.
	fprintf(stderr, "Using %d single, %d pairwise terms out of %d parameters.\n", n_single_terms, n_pw_terms, n_total_params);
	///////////////////////////////////////////////////////////


	for (int i_chr = 0; i_chr < (int)restr_haplocoded_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_haplocoded_var_regs = restr_haplocoded_var_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Sorting the variants per position..\n");

		sort(cur_chr_haplocoded_var_regs->begin(), cur_chr_haplocoded_var_regs->end(), sort_regions);

		// Process each sample.
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			// Process each haplotype
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				for (int i_var = 0; i_var < (int)cur_chr_haplocoded_var_regs->size(); i_var++)
				{
					// Store the var_i info.
					void** var_i_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(i_var)->data);
					//char* var_i_geno_sig = var_i_reg_info[0];
					char* var_i_proxized_geno_sig = (char*)(var_i_reg_info[1]);

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
						cur_var_i_vicinity_allele_sum += (cur_j_allele * per_var_weight[j_var - cur_win_start]);

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
							cur_var_i_vicinity_allele_sum += (cur_j_allele * cur_k_allele * cur_int_term);
						} // k_var loop.
					} // j_var loop.

					// Add error:
					int allele_err = 0;
					if (rng->random_double_ran3() < allele_err_eps)
					{
						haplocoded_var_regs->at(i_var)->score++;
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

// This function tests the simples proxization of variants (genotypes) using modular average of surrounding variants.
void proxize_variants_per_vicinity_modular_average(char* haplocoded_geno_signal_fp, char* sample_ids_list_fp,
													int n_vicinity, int per_var_weight, 
													int avg_geno_mod,
													double allele_err_eps,
													char* proxized_haplocoded_var_regs_op_fp)
{
	fprintf(stderr, "Proxizing %s (%s) using:\n\\sum_{vic @ [-%d, +%d]) (%d x g_i) mod %d\n", 
			haplocoded_geno_signal_fp, sample_ids_list_fp,
			n_vicinity, n_vicinity, per_var_weight, avg_geno_mod);

	// Instantiate the RNG.
	t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

	vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	int max_geno = get_max_genotype_value(haplocoded_var_regs, sample_ids);
	if (max_geno != 3)
	{
		fprintf(stderr, "Genotypes are not haplocoded, make sure they are phased..\n");
		exit(1);
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

		// Replace the variant information.
		haplocoded_var_regs->at(i_var)->data = new_reg_info;
		haplocoded_var_regs->at(i_var)->score = 0;
	} // i_var loop.

	for (int i_chr = 0; i_chr < (int)restr_haplocoded_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Proxizing variants on %s\n", restr_haplocoded_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_haplocoded_var_regs = restr_haplocoded_var_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Sorting the variants per position..\n");

		sort(cur_chr_haplocoded_var_regs->begin(), cur_chr_haplocoded_var_regs->end(), sort_regions);

		// Process each sample.
		for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
		{
			// Process each haplotype
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				for (int i_var = 0; i_var < (int)cur_chr_haplocoded_var_regs->size(); i_var++)
				{
					// Store the var_i info.
					void** var_i_reg_info = (void**)(cur_chr_haplocoded_var_regs->at(i_var)->data);
					//char* var_i_geno_sig = var_i_reg_info[0];
					char* var_i_proxized_geno_sig = (char*)(var_i_reg_info[1]);

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
						int cur_allele = get_allele_per_haplotype(var_j_geno_sig[i_s], i_hap);
						cur_var_i_vicinity_allele_sum += (cur_allele * per_var_weight);

						// Check allele validity.
						if (cur_allele > 1 || cur_allele < 0)
						{
							fprintf(stderr, "Sanity check failed @ %s allele is not valid: %d\n", 
								cur_chr_haplocoded_var_regs->at(j_var)->name,
									cur_allele);
							exit(1);
						}
					} // j_var loop.

					// Add error:
					int allele_err = 0;
					if (rng->random_double_ran3() < allele_err_eps)
					{
						haplocoded_var_regs->at(i_var)->score++;
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


