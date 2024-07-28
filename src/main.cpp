#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <thread>
#include <ctype.h>
#include <string.h>
#include "prxytypr_proxytyper.h"
#include "prxytypr_haplotype_resampling_utils.h"
#include "prxytypr_haplotype_resampling_utils_optimized.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_utils.h"
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_variation_tools.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_matrix_linalg_utils.h"
#include "prxytypr_unphased_panel_resampling_utils.h"
#include "prxytypr_vcf_io_utils.h"
#include "prxytypr_genomics_coords.h"
#include "prxytypr_vector_macros.h"
#include "prxytypr_imputation_utils.h"
#include "prxytypr_nomenclature.h"

#include <vector>
#include <algorithm>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Leakage statistics:\n\
	-get_proxization1_2_proxization_mapped_haplotype_frequency_correlation\n\
	-get_hapfreq_confusion_matrix\n\
	-calculate_proxy2clear_var2var_correlation_stats\n\
	-generic_map_regions_per_src_dest_regs\n\
	-get_consecutive_variant_correlations\n\
	-get_untyped_variant_LD_statistics\n\
	-get_consecutive_block_variant_correlations\n\
	-calculate_Homer_t_statistics_on_proxized_panels\n\
	-pool_summarize_Homer_t_statistics_per_query\n\
	-pool_summarize_Sankararaman_LRT_statistics_per_query\n\
	-calculate_Sankararaman_LRT_statistics_on_proxized_panels\n\
	-write_securegenome_input_files\n\
	-calculate_windowed_stats_Homer_t_statistics_on_proxized_panels\n\
	-calculate_proxy2clear_pairwise_distance_stats\n\
	-haplotype_frequency_attack\n\
	-linking_attack_per_haplotype_frequency_signatures\n\
	-linking_attack_per_haplotype_frequency_signatures_per_ref_panel\n\
Haplotype emission probabilities:\n\
	-calculate_haplotype_emission_probabilities_per_reference_Full_States_multithreaded\n\
	-calculate_haplotype_emission_probabilities_per_reference_Full_States\n\
	-calculate_haplotype_emission_probabilities_per_reference_state_reduced\n\
	-calculate_haplotype_emission_probabilities_per_reference_state_reduced_multithreaded\n\
Haplotype Decoding Attacks:\n\
	-get_allele_error_per_decoded_panel_known_panel\n\
	-decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack\n\
	-decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi\n\
	-decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching\n\
	-decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM\n\
	-decode_site_alleles_per_proxized_reference\n\
Variant Proxization:\n\
	-build_augmenting_variants_mapper\n\
	-augment_tag_variants\n\
	-context_specific_decompose_variants\n\
	-anonymize_genetic_distances\n\
	-generic_map_genotype_regs_per_src_dest_regs\n\
	-anonymize_tag_target_genetic_map_coords: THIS IS THE MAIN OPTION TO ANONYMIZE COORDINATES AND GENETIC MAPS\n\
	-random_per_sample_hap_switch\n\
	-modify_per_site_mixing_parameters\n\
	-generate_save_per_site_mixing_parameters\n\
	-generate_permute_proxizing_parameters: THIS IS THE MAIN OPTION TO GENERATE LOCAL WINDOW PERMUTE PROXIZING OPTION FOR TAG VARIANTS..\n\
	-generate_save_per_site_mixing_parameters_LD_aware: THE MAIN OPTION TO GENERATE THE HASHING/PROXIZATION PARAMETERS..\n\
	-permute_proxize_genotype_signal_regions\n\
	-proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters\n\
	-MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters: THE MAIN OPTION TO HASH/PROXIZE GENOTYPES WITH GIVEN WEIGHTS..\n\
	-MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer\n\
	-proxize_variants_per_vicinity_modular_average\n\
	-proxize_variants_per_vicinity_non_linear_modular_average_uniform\n\
	-generate_untyped_variant_recoded_reference_panel\n\
	-recode_untyped_variant_reference_panel_per_target_permutation: THIS IS THE FUNCTION FOR TARGET VARIANT PERMUTATIONS; THIS IS NOT USED ANY MORE SINCE IT IS REDUNDANT WITH DECOMPOSITION STEP.\n\
	-decode_untyped_variant_reference_panel_per_target_permutation\n\
	-check_singular_allele1_2_haplotype_matrix_right_inversion\n\
Untyped Variant Decomposition:\n\
	-simple_decompose_untyped_variants\n\
	-combine_BEAGLE_imputed_decomposed_genotype_probabilities\n\
	-combine_BEAGLE_imputed_decomposed_genotype_probabilities_multithreaded\n\
Resampling:\n\
	-sample_unphased_genotype_phasings\n\
	-resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb\n\
	-resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling\n\
	-summarize_sampled_segments_per_resampled_haplotype_info\n\
	-compute_AF_per_geno_signal_regions\n\
Misc:\n\
	-intersect\n\
	-extract_genotype_signals_per_VCF_no_buffer_multithreaded\n\
	-validate_system\n\
	-get_R2_per_imputed_genotypes\n\
	-get_query_haplotype_frequency_per_reference\n\
	-locally_permute_indices\n\
	-get_per_window_sampling_probability_per_reference\n\
	-get_per_win_haplotype_frequencies_per_reference\n\
	-compare_per_win_haplotypes\n\
	-convert_BEAGLE_imputed_meta_panel_2_matbed\n", argv[0]);
		exit(1);
	}

	clock_t start_c = clock();
	time_t start_t = time(NULL);

	if (strcmp(argv[1], "-random_phase_genotypes") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -random_phase_genotypes [Genotype signals file path] [Sample ids list file path] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* genotype_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* haplocoded_genotype_matrix_op_fp = argv[4];

		random_phase_het_genotypes(genotype_regs_BED_fp, sample_ids_list_fp, haplocoded_genotype_matrix_op_fp);
	} // -random_phase_haplocoded_genotypes option.
	else if (strcmp(argv[1], "-compute_AF_per_geno_signal_regions") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Genotype signals file path] [Sample ids list file path] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* geno_sigs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		vector<t_annot_region*>* var_regs = load_variant_signal_regions_wrapper(geno_sigs_BED_fp, sample_ids_list_fp);
		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		int max_geno = get_max_genotype_value(var_regs, sample_ids);
		if (max_geno == 3)
		{
			fprintf(stderr, "Genotype is likely haplocoded\n");
		}
		else if (max_geno == 2)
		{
			fprintf(stderr, "Genotype is likely genocoded\n");
		}
		else if (max_geno == 1)
		{
			fprintf(stderr, "***Max geno is %d, not expected***.\n", max_geno);
			fprintf(stderr, "***Max geno is %d, not expected***.\n", max_geno);
			fprintf(stderr, "***Assuming genocoded***.\n");
		}
		else
		{
			fprintf(stderr, "***Max geno is %d, not expected; no signal???***.\n", max_geno);
			exit(1);
		}

		FILE* f_op = open_f(op_fp, "w");
		for (int i_reg = 0; i_reg < (int)var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(var_regs->at(i_reg)->data);
			char* geno_sig = (char*)(cur_reg_info[0]);

			double tot_geno = 0;
			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				if (max_geno == 3)
				{
					tot_geno += get_genotype_per_haplocoded_genotype(geno_sig[i_s]);
				}
				else
				{
					tot_geno += (double)(geno_sig[i_s]);
				}
			} // i_s loop.

			fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t+\n", var_regs->at(i_reg)->chrom,
				translate_coord(var_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(var_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				var_regs->at(i_reg)->name,
				tot_geno / (2 * (int)sample_ids->size()));
		} // i_reg loop.

		close_f(f_op, op_fp);
	} // -compute_AF_per_geno_signal_regions option.
	else if (t_string::compare_strings(argv[1], "-resample_phased_haplotypes_per_state_only_sampling"))
	{
		if (argc != 18)
		{
			fprintf(stderr, "USAGE: %s %s \n\
[Haplocoded genotype signals BED file path (HC)] \n\
[Sample ids list file path] \n\
[3-column Recombination Rate files dir] \n\
[Resampled size] \n\
[N_e_2_n_ref_haplotypes (Ne / N_haplotypes)] \n\
[Min. delta cM genetic distance between anchoring variants] \n\
[Allele error prob (10e-4)] \n\
[Max segment length cutoff in bps] \n\
[Max segment length cutoff in cM] \n\
[Max segment length cutoff in # vars] \n\
[Haplocoded output? (1:Haplocoded/0:Genocoded) \n\
[# threads] \n\
[start position] [end position]\n\
[Save recombination patterns? (0/1)]\n\
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_genotype_matrix_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* recombination_rate_dir = argv[4];
		int n_resampled_sample_size = atof(argv[5]);
		double N_e_2_n_ref_haplotypes = atof(argv[6]);
		double min_cM_delta_per_recomb = atof(argv[7]);
		double allele_error_prob = atof(argv[8]);
		double segment_length_cutoff_bp = atof(argv[9]);
		double segment_length_cutoff_cM = atof(argv[10]);
		double segment_length_cutoff_in_var_number = atof(argv[11]);
		bool use_haplocoded_genotypes = (argv[12][0] == '1') ? (true) : (false);
		int n_threads = atoi(argv[13]);
		int start_posn = atoi(argv[14]);
		int end_posn = atoi(argv[15]);
		bool save_recombination_patterns = (argv[16][0] == '1') ? (true) : (false);
		char* op_prefix = argv[17];

		resample_phased_haplotypes_per_state_only_sampling(haplocoded_genotype_matrix_fp,
			sample_ids_list_fp,
			recombination_rate_dir,
			n_resampled_sample_size,
			N_e_2_n_ref_haplotypes,
			min_cM_delta_per_recomb,
			allele_error_prob,
			segment_length_cutoff_bp,
			segment_length_cutoff_cM,
			segment_length_cutoff_in_var_number,
			NULL,
			use_haplocoded_genotypes,
			n_threads,
			start_posn, end_posn,
			save_recombination_patterns,
			op_prefix);
	} // -resample_phased_haplotypes_per_state_only_sampling option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_VCF") == 0)
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_VCF [VCF file path] [VCF sample ids list file path (Use EpiLeak to extract)] \
[Variant regions BED file path] [chromosome id to process] \
[Binary sequence directory (Necessary for ref matching)] [Match reference allele? (0/1)] \
[Match region names? (0/1)] \
[Haplotype specific encoding (0/1)] \
[Output file path]\n", argv[0]);
			exit(1);
		}

		char* vcf_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		char* var_regions_BED_fp = argv[4];
		char* chr_id_2_process = argv[5];
		char* bin_seq_dir = argv[6];
		bool match_ref_alleles_flag = (argv[7][0] == '1');
		bool match_region_names_flag = (argv[8][0] == '1');
		bool haplotype_specific_encoding = (argv[9][0] == '1');
		char* op_fp = argv[10];

		extract_genotype_signals_per_VCF(vcf_fp,
			vcf_sample_ids_list_fp,
			var_regions_BED_fp,
			chr_id_2_process,
			bin_seq_dir,
			match_ref_alleles_flag,
			match_region_names_flag,
			haplotype_specific_encoding,
			op_fp);

		exit(0);
	} // -extract_genotype_signals_per_1kg_VCF option.
	else if (strcmp(argv[1], "-get_me_VCF_header_per_chrom_lengths") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [chrom sizes file path] [Assembly header identifier] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* chrom_sizes_fp = argv[2];
		char* assembly_header_id = argv[3];
		char* op_fp = argv[4];

		vector<char*>* chrom_size_lines = buffer_file(chrom_sizes_fp);
		fprintf(stderr, "Loaded %d chrom sizes.\n", (int)chrom_size_lines->size());

		FILE* f_op = open_f(op_fp, "w");
		fprintf(f_op, "##fileformat=VCFv4.0\n");
		for (int i_chr = 0; i_chr < (int)chrom_size_lines->size(); i_chr++)
		{
			char cur_chrom[1000];
			int l_cur_chr;

			if (sscanf(chrom_size_lines->at(i_chr), "%s %d", cur_chrom, &l_cur_chr) != 2)
			{
				fprintf(stderr, "Could not parse: %s\n", chrom_size_lines->at(i_chr));
				exit(1);
			}

			// Normalize chromosome identifier.
			normalize_chr_id(cur_chrom);

			fprintf(f_op, "##contig=<ID=%s,length=%d,assembly=%s>\n", cur_chrom, l_cur_chr, assembly_header_id);
		} // i_chr loop.

		fprintf(f_op, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		close_f(f_op, op_fp);

		exit(0);
	}
	else if (strcmp(argv[1], "-get_me_VCF_header_example") == 0)
	{
		fprintf(stderr, "##fileformat=VCFv4.1\n\
##fileDate=20150121\n\
##source=TCGA\n\
##contig=<ID=1,length=249250621,assembly=b37>\n\
##contig=<ID=2,length=243199373,assembly=b37>\n\
##contig=<ID=3,length=198022430,assembly=b37>\n\
##contig=<ID=4,length=191154276,assembly=b37>\n\
##contig=<ID=5,length=180915260,assembly=b37>\n\
##contig=<ID=6,length=171115067,assembly=b37>\n\
##contig=<ID=7,length=159138663,assembly=b37>\n\
##contig=<ID=8,length=146364022,assembly=b37>\n\
##contig=<ID=9,length=141213431,assembly=b37>\n\
##contig=<ID=10,length=135534747,assembly=b37>\n\
##contig=<ID=11,length=135006516,assembly=b37>\n\
##contig=<ID=12,length=133851895,assembly=b37>\n\
##contig=<ID=13,length=115169878,assembly=b37>\n\
##contig=<ID=14,length=107349540,assembly=b37>\n\
##contig=<ID=15,length=102531392,assembly=b37>\n\
##contig=<ID=16,length=90354753,assembly=b37>\n\
##contig=<ID=17,length=81195210,assembly=b37>\n\
##contig=<ID=18,length=78077248,assembly=b37>\n\
##contig=<ID=19,length=59128983,assembly=b37>\n\
##contig=<ID=20,length=63025520,assembly=b37>\n\
##contig=<ID=21,length=48129895,assembly=b37>\n\
##contig=<ID=22,length=51304566,assembly=b37>\n\
##contig=<ID=X,length=155270560,assembly=b37>\n\
##contig=<ID=Y,length=59373566,assembly=b37>\n\
##contig=<ID=MT,length=16569,assembly=b37>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

		exit(0);
	}
	else if (strcmp(argv[1], "-convert_genocoded_signal_regs_2_VCF") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -convert_genocoded_signal_regs_2_VCF [Genotype signal regions file path] [Sample ids list file path] [Output VCF file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_vcf_fp = argv[4];

		convert_genocoded_signal_regs_2_VCF(geno_sig_regs_fp, sample_ids_list_fp, op_vcf_fp);

		exit(0);
	} // -convert_genocoded_signal_regs_2_VCF option.
	else if (strcmp(argv[1], "-convert_genotype_signal_regions_2_VCF") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -convert_genotype_signal_regions_2_VCF [Genotype signal regions file path] [Sample ids list file path] [Phased output? (0/1)] [Output VCF file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		bool phase_op = (argv[4][0] == '1');
		char* op_vcf_fp = argv[5];

		convert_haplocoded_signal_regs_2_VCF(geno_sig_regs_fp, sample_ids_list_fp, phase_op, op_vcf_fp);

		exit(0);
	} // -convert_genotype_signal_regions_2_VCF option.
	else if (strcmp(argv[1], "-convert_genotype_signal_regions_2_VCF_w_header") == 0)
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -convert_genotype_signal_regions_2_VCF_w_header [Genotype signal regions file path] [Sample ids list file path] [Header file] [Phased output? (0/1)] [Output VCF file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* header_file = argv[4];
		bool phase_op = (argv[5][0] == '1');
		char* op_vcf_fp = argv[6];

		convert_haplocoded_signal_regs_2_VCF_w_header(geno_sig_regs_fp, sample_ids_list_fp, header_file, phase_op, op_vcf_fp);

		exit(0);
	} // -convert_genotype_signal_regions_2_VCF_w_header option.
	else if (strcmp(argv[1], "-get_me_VCF_header_per_chrom_lengths") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [chrom sizes file path] [Assembly header identifier] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* chrom_sizes_fp = argv[2];
		char* assembly_header_id = argv[3];
		char* op_fp = argv[4];

		vector<char*>* chrom_size_lines = buffer_file(chrom_sizes_fp);
		fprintf(stderr, "Loaded %d chrom sizes.\n", (int)chrom_size_lines->size());

		FILE* f_op = open_f(op_fp, "w");
		fprintf(f_op, "##fileformat=VCFv4.0\n");
		for (int i_chr = 0; i_chr < (int)chrom_size_lines->size(); i_chr++)
		{
			char cur_chrom[1000];
			int l_cur_chr;

			if (sscanf(chrom_size_lines->at(i_chr), "%s %d", cur_chrom, &l_cur_chr) != 2)
			{
				fprintf(stderr, "Could not parse: %s\n", chrom_size_lines->at(i_chr));
				exit(1);
			}

			// Normalize chromosome identifier.
			normalize_chr_id(cur_chrom);

			fprintf(f_op, "##contig=<ID=%s,length=%d,assembly=%s>\n", cur_chrom, l_cur_chr, assembly_header_id);
		} // i_chr loop.

		fprintf(f_op, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		close_f(f_op, op_fp);
	}
	else if (strcmp(argv[1], "-convert_genotype_signal_regions_2_VCF_w_header_multithreaded") == 0)
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s -convert_genotype_signal_regions_2_VCF_w_header_multithreaded [Genotype signal regions file path] [Sample ids list file path] [# threads] [Header file] [Phased output? (0/1)] [Output VCF file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* header_file = argv[4];
		bool phase_op = (argv[5][0] == '1');
		int n_threads = atoi(argv[6]);
		char* op_vcf_fp = argv[7];

		convert_haplocoded_signal_regs_2_VCF_w_header_multithreaded(geno_sig_regs_fp, sample_ids_list_fp, header_file, phase_op, n_threads, op_vcf_fp);

		exit(0);
	} // -convert_genotype_signal_regions_2_VCF_w_header_multithreaded option.
	else if (strcmp(argv[1], "-extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_multithreaded") == 0)
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s -extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix [Reference Panel matbed file path] [reference panel sample ids list file path] \
[Input Study Panel matbed file path] [Input Study panel sample ids list file path] \
[ref_option file path] [gt_option file path] [Save phased gt option (0/1)] [# threads]\n", argv[0]);
			exit(1);
		}

		char* reference_haplotype_genotype_matrix_matbed_fp = argv[2];
		char* reference_haplotype_genotype_matrix_sample_ids_list_fp = argv[3];
		char* input_genotype_matrix_matbed_fp = argv[4];
		char* input_genotype_matrix_sample_ids_list_fp = argv[5];
		char* ref_option_fp = argv[6];
		char* gt_option_fp = argv[7];
		bool save_phased_gt_option = (argv[8][0] == '1');
		int n_threads = atoi(argv[9]);

		extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix_multithreaded(reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
			input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
			ref_option_fp,
			gt_option_fp,
			save_phased_gt_option,
			n_threads);
	} // -extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix option.
	else if (strcmp(argv[1], "-sort_genosignal_matrix") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Genotype matrix path or prefix] [Geno. matrix sample ids list file] [Output file/prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* matbed_file = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(matbed_file, sample_ids_list_fp);
		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		fprintf(stderr, "Loaded %d regions on %d subjects, sorting and saving..\n", vecsize(geno_sig_regs), vecsize(sample_ids));

		sort(geno_sig_regs->begin(), geno_sig_regs->end(), sort_regions);

		binarize_variant_signal_regions_wrapper(geno_sig_regs, sample_ids, op_fp);
	} // -sort_genosignal_matrix option.
	else if (strcmp(argv[1], "-dump_plain_haplo_signal_regions") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -dump_plain_haplo_signal_regions [Genotype signals matBED file path] [sample ids list file path] [haplotype to extract: 0/1/2(phased together)/3(phased separate)] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		int haplo_2_extract = atoi(argv[4]);
		char* op_fp = argv[5];

		vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		dump_haplo_sig_regs_plain(geno_sig_regs, sample_ids, haplo_2_extract, op_fp);
		} // -dump_plain_haplo_signal_regions option.
	else if (strcmp(argv[1], "-dump_plain_geno_signal_regions") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -dump_plain_geno_signal_regions [Genotype signals matBED file path] [sample ids list file path] [Save regions only? (0/1)] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		bool regs_only_flag = (argv[4][0] == '1');
		char* op_fp = argv[5];

		vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		dump_geno_sig_regs_plain(geno_sig_regs, sample_ids, regs_only_flag, op_fp);
	} // -dump_plain_geno_signal_regions option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_chr_ids") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_chr_ids [Genotype signals BED file path] [sample ids list file path] [chr ids list path] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* chr_ids_list_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_chr_ids(geno_sig_regs_BED_fp, sample_ids_list_fp, chr_ids_list_fp, op_fp);
	}
	else if (strcmp(argv[1], "-extract_genotype_signals_per_region_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_region_list [Genotype signals BED file path] [sample ids list file path] [BED file with regions] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* regions_BED_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_region_list(geno_sig_regs_BED_fp, sample_ids_list_fp, regions_BED_fp, op_fp);
	} // -extract_genotype_signals_per_region_list option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_subsample_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_subsample_list [Genotype signals BED file path] [sample ids list file path] [Subsample ids list file path] [Output file path]\n", argv[0]);
			exit(1);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* subsample_ids_list_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_subsample_list(geno_sig_regs_BED_fp, sample_ids_list_fp, subsample_ids_list_fp, op_fp);
	} // -extract_genotype_signals_per_subsample_list option.
	else if (strcmp(argv[1], "-compute_AF_per_geno_signal_regions") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Genotype signals file path] [Sample ids list file path] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* geno_sigs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		vector<t_annot_region*>* var_regs = load_variant_signal_regions_wrapper(geno_sigs_BED_fp, sample_ids_list_fp);
		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		int max_geno = get_max_genotype_value(var_regs, sample_ids);
		if (max_geno == 3)
		{
			fprintf(stderr, "Genotype is likely haplocoded\n");
		}
		else if (max_geno == 2)
		{
			fprintf(stderr, "Genotype is likely genocoded\n");
		}
		else if (max_geno == 1)
		{
			fprintf(stderr, "***Max geno is %d, not expected***.\n", max_geno);
			fprintf(stderr, "***Max geno is %d, not expected***.\n", max_geno);
			fprintf(stderr, "***Assuming genocoded***.\n");
		}
		else
		{
			fprintf(stderr, "***Max geno is %d, not expected; no signal???***.\n", max_geno);
			exit(1);
		}

		FILE* f_op = open_f(op_fp, "w");
		for (int i_reg = 0; i_reg < (int)var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(var_regs->at(i_reg)->data);
			char* geno_sig = (char*)(cur_reg_info[0]);

			double tot_geno = 0;
			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				if (max_geno == 3)
				{
					tot_geno += get_genotype_per_haplocoded_genotype(geno_sig[i_s]);
				}
				else
				{
					tot_geno += (double)(geno_sig[i_s]);
				}
			} // i_s loop.

			fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t+\n", var_regs->at(i_reg)->chrom,
				translate_coord(var_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(var_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				var_regs->at(i_reg)->name,
				tot_geno / (2 * (int)sample_ids->size()));
		} // i_reg loop.

		close_f(f_op, op_fp);
	} // -compute_AF_per_geno_signal_regions option.
	else if (strcmp(argv[1], "-convert_haplocoded_2_genocoded") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -convert_haplocoded_2_genocoded [Haplocoded genotype matbed file path] \
[Haplocoded genotype matrix sample ids list path] \
[Output genotype signal regions bed file path]\n", argv[0]);
			exit(1);
		}

		char* haplocoded_geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_matbed_fp = argv[4];

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
		vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_geno_sig_regs_fp, sample_ids_list_fp);

		fprintf(stderr, "Loaded %d haplocoded genotype regions for %d samples.\n", (int)haplocoded_geno_regs->size(), (int)sample_ids->size());

		for (int i_reg = 0; i_reg < (int)haplocoded_geno_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
			char* haplocoded_geno_sigs = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < (int)sample_ids->size(); i_s++)
			{
				int genocoded_geno = get_genotype_per_haplocoded_genotype(haplocoded_geno_sigs[i_s]);
				haplocoded_geno_sigs[i_s] = genocoded_geno;
			} // i_s loop.
		} // i_reg loop.

		// Save.
		fprintf(stderr, "Saving to %s.\n", op_matbed_fp);
		//binarize_variant_genotype_signal_regions(haplocoded_geno_regs, NULL, sample_ids, op_matbed_fp);
		binarize_variant_signal_regions_wrapper(haplocoded_geno_regs, sample_ids, op_matbed_fp);

		exit(0);
	} // -convert_haplocoded_2_genocoded option.
	else if (strcmp(argv[1], "-exclude") == 0)
	{
		// Exclude the regions in second region file from the first region file.
		if (argc != 5)
		{
			printf("USAGE: %s -exclude [region1 file path] [region2 file path] [Strand specific (y/n)]\n", argv[0]);
			exit(1);
		}

		char* reg1_fp = argv[2];
		char* reg2_fp = argv[3];
		bool strand_specific = t_string::compare_strings_ci(argv[4], "y") ? (true) : (false);

		// Load the regions: Depends on the file format.
		vector<t_annot_region*>* region_list1 = load_BED_with_line_information(reg1_fp);
		vector<t_annot_region*>* region_list2 = load_BED_with_line_information(reg2_fp);

		vector<t_annot_region*>* excluded_region_list = exclude_annot_regions(region_list1,
			region_list2,
			strand_specific);

		fprintf(stderr, "Found %d excluded regions.\n", vecsize(excluded_region_list));

		printf("Dumping excluded regions.\n");
		FILE* f_exc = fopen("excluded.bed", "w");
		for (int i = 0; i < vecsize(excluded_region_list); i++)
		{
			t_annot_region* src_region = (t_annot_region*)(excluded_region_list->at(i)->data);
			fprintf(f_exc, "%s\n", (char*)(src_region->data));
		}
		fclose(f_exc);

		exit(0);
	}
	else if (strcmp(argv[1], "-intersect") == 0)
	{
		if (argc != 7)
		{
			printf("USAGE: %s -intersect [region1 file path] [region2 file path] [Strand specific (Yes/No)] [Find all overlaps? (Yes/No)] [Type of report: \"Reg1\"/\"Reg2\"/\"Reg12\"/\"Overlap\"]\n", argv[0]);
			exit(1);
		}

		char* reg1_fp = argv[2];
		char* reg2_fp = argv[3];
		bool strand_specific = t_string::compare_strings_ci(argv[4], "yes") ? (true) : (false);
		bool find_all_overlaps = t_string::compare_strings_ci(argv[5], "yes") ? (true) : (false);

		// Select report type.
		int report_type = 0;
		if (t_string::compare_strings_ci(argv[6], "reg1"))
		{
		}
		else if (t_string::compare_strings_ci(argv[6], "reg2"))
		{
			report_type = 1;
		}
		else if (t_string::compare_strings_ci(argv[6], "reg12"))
		{
			report_type = 2;
		}
		else if (t_string::compare_strings_ci(argv[6], "overlap"))
		{
			report_type = 3;
		}
		else
		{
			fprintf(stderr, "Unknown report type, use \"Reg1\", \"Reg2\", \"Reg12\", or \"Overlap\"");
			exit(1);
		}

		// Load the regions: Depends on the file format.
		//vector<t_annot_region*>* region_list1 = load_regions(reg1_format, reg1_fp);
		//vector<t_annot_region*>* region_list2 = load_regions(reg2_format, reg2_fp);
		vector<t_annot_region*>* region_list1 = load_BED_with_line_information(reg1_fp);
		vector<t_annot_region*>* region_list2 = load_BED_with_line_information(reg2_fp);

		vector<t_annot_region*>* intersected_region_list = intersect_annot_regions(region_list1,
			region_list2,
			strand_specific,
			find_all_overlaps);

		FILE* f_intersect = fopen("intersected.bed", "w");
		for (int i = 0; i < vecsize(intersected_region_list); i++)
		{
			t_annot_region* src_region = ((t_intersect_info*)(intersected_region_list->at(i)->data))->src_reg;
			t_annot_region* dest_region = ((t_intersect_info*)(intersected_region_list->at(i)->data))->dest_reg;

			//int src_start = translate_coord(src_region->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			//int src_end = translate_coord(src_region->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			//int dest_start = translate_coord(dest_region->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			//int dest_end = translate_coord(dest_region->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			int overlap_start = translate_coord(intersected_region_list->at(i)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
			int overlap_end = translate_coord(intersected_region_list->at(i)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

			// Depending on the report type, dump the region.
			if (report_type == 0)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, src_start, src_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\n", (char*)(src_region->data));
			}
			else if (report_type == 1)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, dest_start, dest_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\n", (char*)(dest_region->data));
			}
			else if (report_type == 2)
			{
				//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, dest_start, dest_end, intersected_region_list->at(i)->strand);
				fprintf(f_intersect, "%s\t%s\n", (char*)(src_region->data), (char*)(dest_region->data));
			}
			else
			{
				fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, overlap_start, overlap_end, intersected_region_list->at(i)->strand);
			}

			// This is a narrowPeak info.
			//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", src_region->chrom, src_region->start, src_region->end, src_region->strand);
			//fprintf(f_intersect, "%s\t%d\t%d\t.\t.\t%c\n", intersected_region_list->at(i)->chrom, intersected_region_list->at(i)->start, intersected_region_list->at(i)->end, intersected_region_list->at(i)->strand);
		}
		fclose(f_intersect);

		exit(0);
	} // -intersect option.
	else if (t_string::compare_strings(argv[1], "-build_augmenting_variants_mapper"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Augmenting variants BED] [Panel variants (all variants) BED] [Augmentation probability (0-1)] [Augmentation mapper prefix]\n",
				argv[0], argv[1]);
			exit(1);
		}
		
		char* augmenting_variants_BED_fp = argv[2];
		char* panel_variants_BED_fp = argv[3];
		double augment_prob = atof(argv[4]);
		char* augmentation_mapper_op_prefix = argv[5];

		build_augmenting_variants_mapper(augmenting_variants_BED_fp, panel_variants_BED_fp, augment_prob, augmentation_mapper_op_prefix);
	} // -build_augmenting_variants_mapper option.
	else if (t_string::compare_strings(argv[1], "-augment_tag_variants"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Augmentation mapper prefix] [Panel matbed] [Panel sample list] [Augmentation Type (0: Zero, 1: Copy] [Output augmented panel prefix]\n",
				argv[0], argv[1]);
			exit(1);
		}

		char* augmentation_mapper_prefix = argv[2];
		char* haplocoded_geno_signal_fp = argv[3];
		char* sample_ids_list_fp = argv[4];
		int geno_augment_type = atoi(argv[5]);
		char* op_matbed = argv[6];

		augment_tag_variants(augmentation_mapper_prefix, 
							haplocoded_geno_signal_fp, sample_ids_list_fp,
							geno_augment_type, op_matbed);
	} // -augment_tag_variants option.
	else if (t_string::compare_strings(argv[1], "-extract_genotype_signals_per_VCF_no_buffer_multithreaded"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s \n\
[VCF file] \n\
[Sample ids list file path] \n\
[Haplotype specific encoding? (0/1)] \n\
[Add Ref-Alt allele/AF info to name? (0/1)] \n\
[# threads] \n\
[Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* vcf_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		bool haplotype_specific_encoding = (argv[4][0] == '1');
		bool add_AF_info_2_id = (argv[5][0] == '1');
		int n_threads = atoi(argv[6]);
		char* op_prefix = argv[7];

		extract_genotype_signals_per_VCF_no_buffer_multithreaded(vcf_fp,							// This is the VCF file from which the genotypes are read.
			vcf_sample_ids_list_fp,		// This is the sample ids list file path.
			haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
			add_AF_info_2_id,
			n_threads,
			op_prefix);
		exit(0);
	}
	else if (strcmp(argv[1], "-uniquefy_genotype_signal_regions") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -uniquefy_genotype_signal_regions [Matbed files list file] [Sample IDs list file] [Output file]\n", argv[0]);
			exit(1);
		}

		char* geno_reg_sig_file = argv[2];
		char* sample_list_fp = argv[3];
		char* op_fp = argv[4];
		uniquefy_genotype_signal_regions(geno_reg_sig_file, sample_list_fp, op_fp);
	} // -uniquefy_genotype_signal_regions option.
	else if (t_string::compare_strings(argv[1], "-resample_unphased_genotypes_per_state_only_sampling"))
	{
		if (argc != 17)
		{
			fprintf(stderr, "USAGE: %s %s \n\
[Genocoded (unphased) genotype signals BED file path (HC)] \n\
[Sample ids list file path] \n\
[3-column Recombination Rate files dir] \n\
[Resampled size] \n\
[N_e_2_n_ref_haplotypes (Ne / N_haplotypes)] \n\
[Min. delta cM genetic distance between anchoring variants] \n\
[Allele error prob (10e-4)] \n\
[Max segment length cutoff in bps] \n\
[Max segment length cutoff in cM] \n\
[Max segment length cutoff in # vars] \n\
[# threads] \n\
[start position] [end position]\n\
[Save recombination patterns? (0/1)]\n\
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* genocoded_genotype_matrix_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* recombination_rate_dir = argv[4];
		int n_resampled_sample_size = atof(argv[5]);
		double N_e_2_n_ref_haplotypes = atof(argv[6]);
		double min_cM_delta_per_recomb = atof(argv[7]);
		double allele_error_prob = atof(argv[8]);
		double segment_length_cutoff_bp = atof(argv[9]);
		double segment_length_cutoff_cM = atof(argv[10]);
		double segment_length_cutoff_in_var_number = atof(argv[11]);
		int n_threads = atoi(argv[12]);
		int start_posn = atoi(argv[13]);
		int end_posn = atoi(argv[14]);
		bool save_recombination_patterns = (argv[15][0] == '1') ? (true) : (false);
		char* op_prefix = argv[16];

		 resample_unphased_genotypes_per_state_only_sampling(genocoded_genotype_matrix_fp,
			sample_ids_list_fp,
			recombination_rate_dir,
			n_resampled_sample_size,
			N_e_2_n_ref_haplotypes,
			 min_cM_delta_per_recomb, // This is the minimum distance between the variants to switch states.
			allele_error_prob,
			segment_length_cutoff_bp,
			segment_length_cutoff_cM,
			segment_length_cutoff_in_var_number,
			n_threads,
			start_posn, end_posn,
			save_recombination_patterns,
			op_prefix);
	} // -resample_unphased_genotypes_per_state_only_sampling option.
	else if (t_string::compare_strings(argv[1], "-sample_unphased_genotype_phasings"))
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s %s [Unphased geno matrix] [Unphased sample list] [Reference phased matrix] [Reference sample list] [# resamples per subject] [window length] [Distance weighter] [Output geno matrix/id]\n", argv[0], argv[1]);
			exit(1);
		}
		char* unphased_geno_matrix_matbed = argv[2];
		char* unphased_geno_matrix_sample_list = argv[3];
		char* ref_haplocoded_geno_matrix_matbed = argv[4];
		char* ref_haplocoded_geno_matrix_sample_list = argv[5];
		int n_resample_per_subject = atoi(argv[6]);
		int l_win = atoi(argv[7]);
		double dist_weighter = atof(argv[8]);
		char* op_fp = argv[9];

		sample_unphased_genotype_phasings(unphased_geno_matrix_matbed, unphased_geno_matrix_sample_list,
			ref_haplocoded_geno_matrix_matbed, ref_haplocoded_geno_matrix_sample_list,
			n_resample_per_subject,
			l_win, dist_weighter, op_fp);

		exit(1);
	} // -sample_unphased_genotype_phasings option.
	else if (t_string::compare_strings(argv[1], "-validate_system"))
	{
		if (argc != 3)
		{
			fprintf(stderr, "USAGE: %s %s [# requested threads]\n", argv[0], argv[1]);
			exit(1);
		}

		int n_requested_threads = atoi(argv[2]);

		////////////////////////////////////////////////////
		// 2's complement check:
		int test_val = 0;
		memset(&test_val, 0xff, sizeof(int));
		if (test_val != -1)
		{
			fprintf(stderr, "Sanity check failed for 2s complement.\n");
			exit(1);
		}
		else
		{
			fprintf(stderr, "2s complement seems to hold..\n");
		}

		fprintf(stderr, "sizeof(int)=%u bytes; sizeof(size_t)=%u bytes.\n", (unsigned int)sizeof(int), (unsigned int)sizeof(size_t));

		if (sizeof(size_t) != 8 ||
			sizeof(int) != 4)
		{
			fprintf(stderr, "sizeof(size_t) != 8 || sizeof(int) != 4\n");
			exit(1);
		}

		bool test_bool = true;
		fprintf(stderr, "int(bool(true))=%d; \n", (int)(test_bool));
		test_bool = false;
		fprintf(stderr, "int(bool(false))=%d; \n", (int)(test_bool));

		if ((int)(true) != 1 ||
			(int)(false) != 0)
		{
			fprintf(stderr, "Failed: int(bool(true))=%d; int(bool(false))=%d; \n", (int)(true), (int)(false));
			exit(1);
		}

		////////////////////////////////////////////////////
		unsigned int n_cpus = std::thread::hardware_concurrency();
		std::cout << "Number of hardware threads: " << n_cpus << std::endl;

		if (n_cpus < (unsigned int)(2 * n_requested_threads))
		{
			fprintf(stderr, "Too many threads are requested, decrease number to optimize performance: %u cpus vs %u requested.\n", 
					n_cpus, (unsigned int)(n_requested_threads));
			exit(1);
		}

		// Add a check for the number of variants processed by each thread too low or too large maybe problematic.
		// Add a check for the number of variants processed by each thread too low or too large maybe problematic.
		// Add a check for the number of variants processed by each thread too low or too large maybe problematic.
		// Add a check for the number of variants processed by each thread too low or too large maybe problematic.
		// Add a check for the number of variants processed by each thread too low or too large maybe problematic.


		fprintf(stderr, "All looks good (hopefully)..\n");
		exit(0);
	}
	else if (t_string::compare_strings(argv[1], "-get_allele_error_per_decoded_panel_known_panel"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Decoded genotypes matbed file] [Known genotypes matbed file] \
[Sample ids] [Output file]\n", argv[0], argv[1]);
			exit(1);
		}

		char* decoded_panel_matbed_fp = argv[2];
		char* known_panel_matbed_fp = argv[3];
		char* sample_ids_list_fp = argv[4];
		char* f_op = argv[5];
		get_allele_error_per_decoded_panel_known_panel(decoded_panel_matbed_fp, known_panel_matbed_fp, sample_ids_list_fp, f_op);
	}

	else if (t_string::compare_strings(argv[1], "-calculate_haplotype_emission_probabilities_per_reference_Full_States_multithreaded"))
	{
		if (argc != 12)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded reference matbed] [Reference sample list] \
[Haplocoded Query matbed] [Query sample list] \
[Recomb. rates directory] [Min tag-tag distance (cM)] \
[N_e] [allele eps.] [# threads] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* reference_tag_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* recombination_rate_dir = argv[6];
		double min_tar2tag_cM = atof(argv[7]);
		double N_e = atof(argv[8]);
		double allele_eps = atof(argv[9]);
		int n_threads = atoi(argv[10]);
		char* fore_back_output_fp = argv[11];

		calculate_haplotype_emission_probabilities_per_reference_Full_States_multithreaded(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			recombination_rate_dir,
			min_tar2tag_cM,
			N_e, allele_eps,
			n_threads,
			fore_back_output_fp);
	} // -calculate_haplotype_emission_probabilities_per_reference_Full_States option.
	else if (t_string::compare_strings(argv[1], "-calculate_haplotype_emission_probabilities_per_reference_state_reduced"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded reference matbed] [Reference sample list] \
	[Haplocoded Query matbed] [Query sample list] \
	[Recomb. rates directory] [# of variants per block] \
	[N_e] [allele eps.] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* reference_tag_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* recombination_rate_dir = argv[6];
		double l_blocks = atoi(argv[7]);
		double N_e_2_n_ref_haplotypes = atof(argv[8]);
		double allele_eps = atof(argv[9]);
		char* fore_back_output_fp = argv[10];

		calculate_haplotype_emission_probabilities_per_reference_State_Reduced(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			recombination_rate_dir,
			N_e_2_n_ref_haplotypes, allele_eps, l_blocks,
			fore_back_output_fp);
	} // -calculate_haplotype_emission_probabilities_per_reference_Full_States option.
	else if (t_string::compare_strings(argv[1], "-calculate_haplotype_emission_probabilities_per_reference_state_reduced_multithreaded"))
	{
		if (argc != 12)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded reference matbed] [Reference sample list] \
[Haplocoded Query matbed] [Query sample list] \
[Recomb. rates directory] [# variants per block] \
[N_e / # ref haplotypes] [allele eps.] [# threads] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* reference_tag_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* recombination_rate_dir = argv[6];
		int l_blocks = atoi(argv[7]);
		double N_e_2_n_ref_haplotypes = atof(argv[8]);
		double allele_eps = atof(argv[9]);
		int n_threads = atoi(argv[10]);
		char* fore_back_output_fp = argv[11];

		calculate_haplotype_emission_probabilities_per_reference_State_Reduced_multithreaded(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			recombination_rate_dir,
			N_e_2_n_ref_haplotypes, allele_eps,
			l_blocks,
			n_threads,
			fore_back_output_fp);
	} // -calculate_haplotype_emission_probabilities_per_reference_Full_States option.
	else if (t_string::compare_strings(argv[1], "-calculate_haplotype_emission_probabilities_per_reference_Full_States"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded reference matbed] [Reference sample list] \
[Haplocoded Query matbed] [Query sample list] \
[Recomb. rates directory] [Min tag-tag distance (cM)] \
[N_e] [allele eps.] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* reference_tag_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* recombination_rate_dir = argv[6];
		double min_tar2tag_cM = atof(argv[7]);
		double N_e = atof(argv[8]);
		double allele_eps = atof(argv[9]);
		char* fore_back_output_fp = argv[10];

		calculate_haplotype_emission_probabilities_per_reference_Full_States(reference_tag_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			recombination_rate_dir,
			min_tar2tag_cM,
			N_e, allele_eps,
			fore_back_output_fp);
	} // -calculate_haplotype_emission_probabilities_per_reference_Full_States option.


	else if (t_string::compare_strings(argv[1], "-get_R2_per_imputed_genotypes"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Imputed genotypes matrix file path] [Imputed genotypes sample ids path] \
[Known genotypes matrix bed file path] [Known genotypes samples path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* imputed_genotypes_fp = argv[2];
		char* imputed_sample_ids_list_fp = argv[3];
		char* known_genotypes_fp = argv[4];
		char* known_sample_ids_list_fp = argv[5];

		get_R2_per_imputed_genotypes(imputed_genotypes_fp, imputed_sample_ids_list_fp, known_genotypes_fp, known_sample_ids_list_fp);
	}
	else if (t_string::compare_strings(argv[1], "-get_query_haplotype_frequency_per_reference"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Query haplocoded matbed] [Query smaple ids list] [Reference haplocoded matbed] [Reference sample ids list] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_haplocoded_geno_fp = argv[2];
		char* query_sample_ids_fp = argv[3];
		char* ref_haplocoded_geno_fp = argv[4];
		char* ref_sample_ids_fp = argv[5];
		char* op_prefix = argv[6];

		get_query_haplotype_frequency_per_reference(query_haplocoded_geno_fp, query_sample_ids_fp,
			ref_haplocoded_geno_fp, ref_sample_ids_fp,
			op_prefix);
	} // -get_query_haplotype_frequency_per_reference option.
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack"))
	{
		if (argc != 15)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Reference sample list file] \
[Recomb. maps directory] \
[Variant start index] [Variant end index] \
[kmer vicinity length] \
[kmer concordance log weight (Higher means more concordant consecutive kmers)] \
[N_e] [# query subjects 2 decode] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_sample_ids_fp = argv[6];
		char* recombination_rate_dir = argv[7];
		int var_start_i = atoi(argv[8]);
		int var_end_i = atoi(argv[9]);
		int n_kmer_vicinity_vars = atoi(argv[10]);
		double kmer_concordance_log_weight = atof(argv[11]);
		double N_e = atof(argv[12]);
		int n_query_subj_2_decode = atoi(argv[13]);
		char* op_prefix = argv[14];

		decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
			ref_original_haplocoded_geno_fp, ref_sample_ids_fp,
			recombination_rate_dir,
			var_start_i, var_end_i,
			n_kmer_vicinity_vars,
			kmer_concordance_log_weight,
			N_e,
			n_query_subj_2_decode,
			op_prefix);
	} // -decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack option.
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_MT"))
	{
		if (argc != 16)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Reference sample list file] \
[Recomb. maps directory] \
[Variant start index] [Variant end index] \
[kmer vicinity length] \
[kmer concordance log weight (Higher means more concordant consecutive kmers)] \
[N_e] [# query subjects 2 decode] [# threads to use] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_sample_ids_fp = argv[6];
		char* recombination_rate_dir = argv[7];
		int var_start_i = atoi(argv[8]);
		int var_end_i = atoi(argv[9]);
		int n_kmer_vicinity_vars = atoi(argv[10]);
		double kmer_concordance_log_weight = atof(argv[11]);
		double N_e = atof(argv[12]);
		int n_query_subj_2_decode = atoi(argv[13]);
		int n_threads = atoi(argv[14]);
		char* op_prefix = argv[15];

		decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_MT(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
			ref_original_haplocoded_geno_fp, ref_sample_ids_fp,
			recombination_rate_dir,
			var_start_i, var_end_i,
			n_kmer_vicinity_vars,
			kmer_concordance_log_weight,
			N_e,
			n_query_subj_2_decode, n_threads,
			op_prefix);
	} // -decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching option.
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi"))
	{
		if (argc != 15)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Reference sample list file] \
[Recomb. maps directory] \
[Variant start index] [Variant end index] \
[kmer vicinity length] \
[kmer concordance log weight (Higher means more concordant consecutive kmers)] \
[N_e] [# query subjects 2 decode] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_sample_ids_fp = argv[6];
		char* recombination_rate_dir = argv[7];
		int var_start_i = atoi(argv[8]);
		int var_end_i = atoi(argv[9]);
		int n_kmer_vicinity_vars = atoi(argv[10]);
		double kmer_concordance_log_weight = atof(argv[11]);
		double N_e = atof(argv[12]);
		int n_query_subj_2_decode = atoi(argv[13]);
		char* op_prefix = argv[14];

		decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
			ref_original_haplocoded_geno_fp, ref_sample_ids_fp,
			recombination_rate_dir,
			var_start_i, var_end_i,
			n_kmer_vicinity_vars,
			kmer_concordance_log_weight,
			N_e,
			n_query_subj_2_decode,
			op_prefix);
	} // -decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching option.
	else if (t_string::compare_strings(argv[1], "-get_proxization1_2_proxization_mapped_haplotype_frequency_correlation"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Query panel matbed] [Proxy1 matbed] [Proxy2 matbed] [Query sample list file] [# vicinity] [Output file prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* proxy1_haplocoded_geno_signal_fp = argv[3];
		char* proxy2_haplocoded_geno_signal_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int l_hap_win = atoi(argv[6]);
		char* op_prefix = argv[7];

		// This option correlates the frequencies of haplotypes between two different proxies of the same reference panel.
		get_proxization1_2_proxization_mapped_haplotype_frequency_correlation(haplocoded_geno_signal_fp,
			proxy1_haplocoded_geno_signal_fp,
			proxy2_haplocoded_geno_signal_fp,
			sample_ids_list_fp,
			l_hap_win,
			op_prefix);
	} // -get_proxization1_2_proxization_mapped_haplotype_frequency_correlation option.
	else if (t_string::compare_strings(argv[1], "-random_per_sample_hap_switch"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Query panel matbed] [Query sample list file] \
[Switch probability] [Output matbed file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		double switch_prob = atof(argv[4]);
		char* op_matbed_fp = argv[5];

		random_per_sample_hap_switch(haplocoded_geno_signal_fp, sample_ids_list_fp, switch_prob, op_matbed_fp);
	} // -random_per_sample_hap_switch option.
	else if (t_string::compare_strings(argv[1], "-linking_attack_per_haplotype_frequency_signatures_per_ref_panel"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Query panel matbed] [Query sample list file] \
[Reference panel matbed] [Reference sample list file] \
[Proxized panel matbed] [Proxized sample list file] \
[# variants per window] [# variants per step] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* target_query_panel_matbed_fp = argv[2];
		char* target_query_sample_list_fp = argv[3];
		char* ref_panel_matbed_fp = argv[4];
		char* ref_panel_sample_list_fp = argv[5];
		char* target_proxy_panel_matbed_fp = argv[6];
		char* target_proxy_panel_sample_list_fp = argv[7];
		int n_vicinity = atoi(argv[8]);
		int n_step = atoi(argv[9]);
		char* op_fp = argv[10];

		linking_attack_per_haplotype_frequency_signatures_per_ref_panel(target_query_panel_matbed_fp, target_query_sample_list_fp,
			ref_panel_matbed_fp, ref_panel_sample_list_fp,
			target_proxy_panel_matbed_fp, target_proxy_panel_sample_list_fp,
			n_vicinity,
			n_step,
			op_fp);
	} // -linking_attack_per_haplotype_frequency_signatures_per_ref_panel
	else if (t_string::compare_strings(argv[1], "-resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling"))
	{
		if (argc != 17)
		{
			fprintf(stderr, "USAGE: %s %s \n\
[Haplocoded tag genotype signals BED file path (HC)] \n\
[Haplocoded target genotype signals BED file path (HC)] \n\
[Sample ids list file path] \n\
[3-column Recombination Rate files dir] \n\
[Resampled size] \n\
[N_e_2_n_ref_haplotypes (Ne / N_haplotypes)] \n\
[Minimum cM different between recomb. points] \n\
[Allele error prob (10e-4)] \n\
[Max segment length cutoff in bps] \n\
[Max segment length cutoff in cM] \n\
[Max segment length cutoff in # vars] \n\
[# threads] \n\
[start position] [end position]\n\
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_genotype_matrix_fp = argv[2];
		char* haplocoded_target_genotype_matrix_fp = argv[3];
		char* sample_ids_list_fp = argv[4];
		char* recombination_rate_dir = argv[5];
		int n_resampled_sample_size = atof(argv[6]);
		double N_e_2_n_ref_haplotypes = atof(argv[7]);
		double min_cM_delta_per_recomb = atof(argv[8]);
		double allele_error_prob = atof(argv[9]);
		double segment_length_cutoff_bp = atof(argv[10]);
		double segment_length_cutoff_cM = atof(argv[11]);
		double segment_length_cutoff_in_var_number = atof(argv[12]);
		int n_threads = atoi(argv[13]);
		int start_posn = atoi(argv[14]);
		int end_posn = atoi(argv[15]);
		char* op_fp = argv[16];


		resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb_tag_anchored_target_resampling(haplocoded_genotype_matrix_fp,
			haplocoded_target_genotype_matrix_fp,
			sample_ids_list_fp,
			recombination_rate_dir,
			n_resampled_sample_size,
			N_e_2_n_ref_haplotypes,
			min_cM_delta_per_recomb,
			allele_error_prob,
			segment_length_cutoff_bp,
			segment_length_cutoff_cM,
			segment_length_cutoff_in_var_number,
			NULL,
			n_threads,
			start_posn, end_posn,
			op_fp);
	} // -resample_phased_haplotypes_per_AF_priors_per_recombination_rates_multithreaded option.
	else if (t_string::compare_strings(argv[1], "-resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb"))
	{
		if (argc != 15)
		{
			fprintf(stderr, "USAGE: %s -resample_phased_haplotypes_per_AF_priors_per_recombination_rates_multithreaded \n\
[Haplocoded genotype signals BED file path (HC)] \n\
[Sample ids list file path] \n\
[3-column Recombination Rate files dir] \n\
[Upsampling rate] \n\
[N_e_2_n_ref_haplotypes (Ne / N_haplotypes)] \n\
[Allele error prob (10e-4)] \n\
[# threads] \n\
[start position] [end position]\n\
[Output file prefix]\n", argv[0]);
			exit(1);
		}

		char* haplocoded_genotype_matrix_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* recombination_rate_dir = argv[4];
		int n_resampled_sample_size = atof(argv[5]);
		double N_e_2_n_ref_haplotypes = atof(argv[6]);
		double allele_error_prob = atof(argv[7]);
		double segment_length_cutoff_bp = atof(argv[8]);
		double segment_length_cutoff_cM = atof(argv[9]);
		double segment_length_cutoff_in_var_number = atof(argv[10]);
		int n_threads = atoi(argv[11]);
		int start_posn = atoi(argv[12]);
		int end_posn = atoi(argv[13]);
		char* op_prefix = argv[14];

		resample_phased_haplotypes_per_recombination_rates_length_cutoff_multithreaded_save_recomb(haplocoded_genotype_matrix_fp,
			sample_ids_list_fp,
			recombination_rate_dir,
			n_resampled_sample_size,
			N_e_2_n_ref_haplotypes,
			allele_error_prob,
			segment_length_cutoff_bp,
			segment_length_cutoff_cM,
			segment_length_cutoff_in_var_number,
			NULL,
			n_threads,
			start_posn, end_posn,
			op_prefix);
	} // -resample_phased_haplotypes_per_AF_priors_per_recombination_rates_multithreaded option.
	else if (t_string::compare_strings(argv[1], "-linking_attack_per_haplotype_frequency_signatures"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Reference1 matbed file] [Reference1 sample id's] \
[Reference2 matbed file] [Reference2 sample id's] \
[# variants per window] [# variants per step] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* ref1_panel_matbed_fp = argv[2];
		char* ref1_sample_list_fp = argv[3];
		char* ref2_panel_matbed_fp = argv[4];
		char* ref2_sample_list_fp = argv[5];
		int n_vars_per_win = atoi(argv[6]);
		int n_step_vars = atoi(argv[7]);
		char* op_fp = argv[8];

		linking_attack_per_haplotype_frequency_signatures(ref1_panel_matbed_fp, ref1_sample_list_fp,
			ref2_panel_matbed_fp, ref2_sample_list_fp,
			n_vars_per_win,
			n_step_vars,
			op_fp);
	} // -linking_attack_per_haplotype_frequency_signatures option.
	else if (t_string::compare_strings(argv[1], "-locally_permute_indices"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s %s [# elements] [# vicinity]\n", argv[0], argv[1]);
			exit(1);
		}

		int n_elements = atoi(argv[2]);
		int n_vicinity = atoi(argv[3]);
		t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());
		vector<int>* shuffled_i = locally_permute_indices(n_elements, n_vicinity, rng);

		if((int)shuffled_i->size() != n_elements)
		{
			fprintf(stderr, "Sanity check failed: %d vs %d elements\n", (int)shuffled_i->size(), n_elements);
			exit(1);
		}

		vector<int>* sorted_shuffled_i = new vector<int>();
		for(int i = 0; i < (int)shuffled_i->size(); i++)
		{
			fprintf(stderr, "%d: %d\n", i, shuffled_i->at(i));
			sorted_shuffled_i->push_back(shuffled_i->at(i));
		} // i loop.

		sort(sorted_shuffled_i->begin(), sorted_shuffled_i->end());

		for (int i = 0; i < (int)shuffled_i->size(); i++)
		{
			if(i != sorted_shuffled_i->at(i))
			{
				fprintf(stderr, "Sanity check failed: %d, %d\n",
						i, sorted_shuffled_i->at(i));

				exit(1);
			}
		} // i loop.
	} // -locally_permute_indices option.
	else if (t_string::compare_strings(argv[1], "-generic_map_genotype_regs_per_src_dest_regs"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Target genotype signal matbed] [Sample list file] [Source coordinates BED file] [Destination coordinates BED file] [# threads] [Output panel file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* target_geno_sig_regs_matbed_fp = argv[2];
		char* sample_list_fp = argv[3];
		char* source_regs_BED_fp = argv[4];
		char* dest_regs_BED_fp = argv[5];
		int n_threads = atoi(argv[6]);
		char* op_fp = argv[7];

		vector<t_annot_region*>* target_geno_sig_regs = load_variant_signal_regions_wrapper(target_geno_sig_regs_matbed_fp, sample_list_fp);
		vector<char*>* sample_ids = buffer_file(sample_list_fp);

		fprintf(stderr, "Loaded %d genotype signal regions for %d subjects.\n", (int)target_geno_sig_regs->size(), (int)sample_ids->size());

		// Map the coordinates.
		generic_map_regions(target_geno_sig_regs, source_regs_BED_fp, dest_regs_BED_fp);
		binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(target_geno_sig_regs, sample_ids, true, n_threads, op_fp);
	} // -generic_map_regions_per_src_des_regs option.
	else if (t_string::compare_strings(argv[1], "-generic_map_regions_per_src_dest_regs"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Target regions BED file] [Source coordinates BED file] [Destination coordinates BED file] [Output directory]\n", argv[0], argv[1]);

			exit(1);
		}

		char* target_regs_BED_fp = argv[2];
		char* source_regs_BED_fp = argv[3];
		char* dest_regs_BED_fp = argv[4];
		char* op_fp = argv[5];

		vector<t_annot_region*>* target_regs = load_BED(target_regs_BED_fp);
		generic_map_regions(target_regs, source_regs_BED_fp, dest_regs_BED_fp);

		dump_BED(op_fp, target_regs);
	} // -generic_map_regions_per_src_des_regs option.
	else if (t_string::compare_strings(argv[1], "-anonymize_tag_target_genetic_map_coords_per_query_ref_variants"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Tag Variant regions BED file] [Tag/Target Variant regions BED file] [Genetic map file] [cM noise S.D.] [Anonymized chromosome length] [Output directory]\n", argv[0], argv[1]);

			exit(1);
		}

		char* tag_variants_BED_fp = argv[2];
		char* tag_target_variants_BED_fp = argv[3];
		char* genetic_map_file = argv[4];
		double cM_noise_SD = atof(argv[5]);
		unsigned int max_base_posn = (unsigned int)(atoi(argv[6]));
		char* op_dir = argv[7];

		anonymize_tag_target_genetic_map_coordinates_per_ref_query_variants(tag_variants_BED_fp, tag_target_variants_BED_fp,
																			genetic_map_file, cM_noise_SD, max_base_posn, op_dir);
	} // -anonymize_tag_target_genetic_map_coords

	else if (t_string::compare_strings(argv[1], "-anonymize_tag_target_genetic_map_coords"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Tag Variant regions BED file] [Target Variant regions BED file] [Genetic map file] [cM noise S.D.] [Anonymized chromosome length] [Output directory]\n", argv[0], argv[1]);

			exit(1);
		}

		char* tag_variants_BED_fp = argv[2];
		char* target_variants_BED_fp = argv[3];
		char* genetic_map_file = argv[4];
		double cM_noise_SD = atof(argv[5]);
		unsigned int max_base_posn = (unsigned int)(atoi(argv[6]));
		char* op_dir = argv[7];

		anonymize_tag_target_genetic_map_coordinates(tag_variants_BED_fp, target_variants_BED_fp,
														genetic_map_file, cM_noise_SD, max_base_posn, op_dir);
	} // -anonymize_tag_target_genetic_map_coords
	else if (t_string::compare_strings(argv[1], "-anonymize_genetic_distances"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Variant regions BED file] [Genetic map file] [cM noise S.D.] [Output genetic map file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* variants_BED_fp = argv[2];
		char* genetic_map_file = argv[3];
		double cM_noise_SD = atof(argv[4]);
		char* op_fp = argv[5];

		anonymize_genetic_distances(variants_BED_fp, genetic_map_file, cM_noise_SD, op_fp);
	} // -anonymize_genetic_distances option.
	else if (t_string::compare_strings(argv[1], "-permute_proxize_genotype_signal_regions"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Panel matbed file] [Sample ids list file] [Permuted proxizing mapping BED] [# threads] [Output matbed file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* permute_proxy_mappings_BED_fp = argv[4];
		int n_threads = atoi(argv[5]);
		char* proxized_haplocoded_var_regs_op_fp = argv[6];

		proxize_variants_per_locality_permutation(haplocoded_geno_signal_fp, sample_ids_list_fp,
			permute_proxy_mappings_BED_fp,
			n_threads,
			proxized_haplocoded_var_regs_op_fp);
	} // -permute_proxize_genotype_signal_regions option.
	else if (t_string::compare_strings(argv[1], "-generate_permute_proxizing_parameters"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [BED file] [# vicinity regs] [Local permutation probability] [Inversion prob.] [Output BED]\n", argv[0], argv[1]);

			exit(1);
		}

		char* var_regs_BED_fp = argv[2];
		int n_vicinity = atoi(argv[3]);
		double local_permute_prob = atof(argv[4]);
		double geno_inv_prob = atof(argv[5]);
		char* op_BED_fp = argv[6];

		generate_permute_proxizing_parameters(var_regs_BED_fp, n_vicinity, local_permute_prob, geno_inv_prob, op_BED_fp);
	} // -sliding_window_permute_regions option.
	else if (t_string::compare_strings(argv[1], "-simple_decompose_untyped_variants_multithreaded"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Typed haplocoded geno matbed] [Untyped haplocoded geno matbed] [Panel sample list file] \
[Min Alternate Allele Frequency to partition] [Shuffle partitioned variants? (0/1)] [# threads] [Output file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* typed_var_geno_sig_regs_fp = argv[2];
		char* untyped_var_geno_sig_regs_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		double min_AAF_per_decomp_var = atof(argv[5]);
		bool shuffle_decomp_variants = (argv[6][0] == '1');
		int n_threads = atoi(argv[7]);
		char* op_prefix = argv[8];

		simple_decompose_untyped_variants_multithreaded(typed_var_geno_sig_regs_fp, untyped_var_geno_sig_regs_fp, 
														panel_sample_list_fp, min_AAF_per_decomp_var,
														shuffle_decomp_variants, n_threads, op_prefix);
	} // -simple_decompose_untyped_variants_multithreaded
	else if (t_string::compare_strings(argv[1], "-simple_decompose_untyped_variants"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Typed haplocoded geno matbed] [Untyped haplocoded geno matbed] [Panel sample list file] \
[Min Alternate Allele Frequency to partition] [Shuffle partitioned variants? (0/1)] [# threads] [Output file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* typed_var_geno_sig_regs_fp = argv[2];
		char* untyped_var_geno_sig_regs_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		double min_AAF_per_decomp_var = atof(argv[5]);
		bool shuffle_decomp_variants = (argv[6][0] == '1');
		int n_threads = atoi(argv[7]);
		char* op_prefix = argv[8];

		simple_decompose_untyped_variants(typed_var_geno_sig_regs_fp, untyped_var_geno_sig_regs_fp, panel_sample_list_fp, min_AAF_per_decomp_var, 
			shuffle_decomp_variants, n_threads, op_prefix);
	} // -simple_decompose_untyped_variants
	else if (t_string::compare_strings(argv[1], "-combine_BEAGLE_imputed_decomposed_genotype_probabilities_multithreaded"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [BEAGLE imputed vcf file path] \
[Variant decomposed variants interval file] [Panel sample list file] [# threads] [Output prefix]\n", argv[0], argv[1]);

			exit(1);
		}

		char* beagle_imputed_vcf_fp = argv[2];
		char* variant_decomposed_intervals_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		int n_threads = atoi(argv[5]);
		char* recomposed_geno_panel_matbed = argv[6];
		combine_BEAGLE_imputed_decomposed_genotype_probabilities_multithreaded(beagle_imputed_vcf_fp,
			variant_decomposed_intervals_fp,
			panel_sample_list_fp,
			n_threads,
			recomposed_geno_panel_matbed);
	}
	else if (t_string::compare_strings(argv[1], "-combine_BEAGLE_imputed_decomposed_genotype_probabilities"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [BEAGLE imputed vcf file path] \
[Variant decomposed variants interval file] [Panel sample list file] [Output prefix]\n", argv[0], argv[1]);

			exit(1);
		}

		char* beagle_imputed_vcf_fp = argv[2];
		char* untyped_variant_decomposed_intervals_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		char* recomposed_geno_panel_matbed = argv[5];

		combine_BEAGLE_imputed_decomposed_genotype_probabilities(beagle_imputed_vcf_fp,
			untyped_variant_decomposed_intervals_fp,
			panel_sample_list_fp,
			recomposed_geno_panel_matbed);
	} // -combine_BEAGLE_imputed_decomposed_genotype_probabilities option.
	else if (t_string::compare_strings(argv[1], "-get_untyped_variant_LD_statistics"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Typed haplocoded geno matbed] [Untyped haplocoded geno matbed] [Panel sample list file] [# tag-tag blocks to process] [Output file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* typed_var_geno_sig_regs_fp = argv[2];
		char* untyped_var_geno_sig_regs_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		int n_blocks_2_process = atoi(argv[5]);
		char* op_prefix = argv[6];

		get_untyped_variant_LD_statistics(typed_var_geno_sig_regs_fp, untyped_var_geno_sig_regs_fp, panel_sample_list_fp, n_blocks_2_process, op_prefix);
	} // -get_consecutive_block_variant_correlations
	else if (t_string::compare_strings(argv[1], "-get_consecutive_block_variant_correlations"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Panel haplocoded geno matbed] [Panel sample list file] [Block length] [Step size (in vars)] [Output file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* panel_target_geno_sig_regs_fp = argv[2];
		char* panel_sample_list_fp = argv[3];
		int l_block = atoi(argv[4]);
		int l_step = atoi(argv[5]);
		char* op_fp = argv[6];

		get_consecutive_block_variant_correlations(panel_target_geno_sig_regs_fp, panel_sample_list_fp, l_block, l_step, op_fp);
	} // -get_consecutive_block_variant_correlations option.
	else if (t_string::compare_strings(argv[1], "-get_consecutive_variant_correlations"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "USAGE: %s %s [Panel haplocoded geno matbed] [Panel sample list file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* panel_target_geno_sig_regs_fp = argv[2];
		char* panel_sample_list_fp = argv[3];
		get_consecutive_variant_correlations(panel_target_geno_sig_regs_fp, panel_sample_list_fp);
	} // -get_consecutive_variant_correlations option.
	else if (t_string::compare_strings(argv[1], "-calculate_proxy2clear_var2var_correlation_stats"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Proxized panel haplocoded matbed] [Proxy sample list] [Original panel haplocoded matbed] [Orig. sample list]\n", argv[0], argv[1]);

			exit(1);
		}

		char* proxized_panel_matbed = argv[2];
		char* proxized_panel_sample_list_fp = argv[3];
		char* cleartext_panel_matbed = argv[4];
		char* cleartext_panel_sample_list_fp = argv[5];

		calculate_proxy2clear_var2var_correlation_stats(proxized_panel_matbed, proxized_panel_sample_list_fp,
														cleartext_panel_matbed, cleartext_panel_sample_list_fp);
	} // -calculate_proxy2clear_var2var_correlation_stats option.
	else if (t_string::compare_strings(argv[1], "-decode_untyped_variant_reference_panel_per_target_permutation"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Recoded Reference panel haplocoded target matbed] \
[Proxy-2-target mapping BED file] [Panel sample list file] [Output decoded reference panel matbed file]\n", argv[0], argv[1]);

			exit(1);
		}

		char* recoded_ref_panel_target_haplocoded_matbed_fp = argv[2];
		char* proxy_2_target_mapping_BED_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		char* op_matbed_fp = argv[5];

		decode_untyped_variant_reference_panel_per_target_permutation(recoded_ref_panel_target_haplocoded_matbed_fp,
			proxy_2_target_mapping_BED_fp,
			panel_sample_list_fp,
			op_matbed_fp);
	}
	else if (t_string::compare_strings(argv[1], "-recode_untyped_variant_reference_panel_per_target_permutation"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Ref. panel haplocoded tag matbed] [Ref. panel haplocoded target matbed] \
[Panel sample list file] [Untyped variant permutation vicinity size] [Untyped genotype allele switch probability] [Output prefix]\n", argv[0], argv[1]);

			exit(1);
		}

		char* ref_panel_tag_haplocoded_matbed_fp = argv[2];
		char* ref_panel_target_haplocoded_matbed_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		int untyped_var_perm_n_vicinity = atoi(argv[5]);
		double untyped_var_allele_switch_prob = atof(argv[6]);
		char* op_prefix = argv[7];

		recode_untyped_variant_reference_panel_per_target_permutation(ref_panel_tag_haplocoded_matbed_fp,
			ref_panel_target_haplocoded_matbed_fp,
			panel_sample_list_fp, 
			untyped_var_perm_n_vicinity,
			untyped_var_allele_switch_prob,
			op_prefix);
	} // -recode_untyped_variant_reference_panel_per_target_permutation 
	else if (t_string::compare_strings(argv[1], "-check_singular_allele1_2_haplotype_matrix_right_inversion"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [# vars] [# haplotypes] [allele1 prob.]\n", argv[0], argv[1]);

			exit(1);
		}

		int n_block_target_vars = atoi(argv[2]);
		int n_panel_haplotypes = atoi(argv[3]);
		double allele1_prob = atof(argv[4]);

		t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

		double** allele1_2_hap = get_random_allele1_2_hap_matrix(rng, n_block_target_vars, n_panel_haplotypes, allele1_prob);
		fprintf(stderr, "Original matrix:\n");
		print_matrix(allele1_2_hap, n_block_target_vars, n_panel_haplotypes);
	
		vector<int>* independent_var_row_i = matrix_get_independent_row_indices(allele1_2_hap, n_block_target_vars, n_panel_haplotypes);

		fprintf(stderr, "Found %d independent variants:\n", (int)independent_var_row_i->size());

		int n_ind_rows = (int)independent_var_row_i->size();
		double** allele1_2_hap_ind = allocate_matrix(n_ind_rows, n_panel_haplotypes);
		int cur_ind_row_i = 0;
		for (int i_var = 0; i_var < (int)independent_var_row_i->size(); i_var++)
		{
			fprintf(stderr, "%d\n", independent_var_row_i->at(i_var));

			int cur_row_i = independent_var_row_i->at(i_var);
			for (int i_hap = 0; i_hap < n_panel_haplotypes; i_hap++)
			{
				allele1_2_hap_ind[cur_ind_row_i][i_hap] = allele1_2_hap[cur_row_i][i_hap];
			} // i_hap loop.

			cur_ind_row_i++;
		} // i_var loop.

		fprintf(stderr, "Copied %d independent variants (rows):\n", cur_ind_row_i);
		print_matrix(allele1_2_hap_ind, n_ind_rows, n_panel_haplotypes);

		// Do inverse operation:
		double** Bt = transpose_matrix(allele1_2_hap_ind, n_ind_rows, n_panel_haplotypes, NULL);
		double** BBt = matrix_multiply(allele1_2_hap_ind, n_ind_rows, n_panel_haplotypes, Bt, n_panel_haplotypes, n_ind_rows, NULL);

		fprintf(stderr, "BBt:\n");
		print_matrix(BBt, n_ind_rows, n_ind_rows);

		fprintf(stderr, "BBt (linear):\n");
		for (int i_row = 0; i_row < n_ind_rows; i_row++)
		{
			for (int i_col = 0; i_col < n_ind_rows; i_col++)
			{
				fprintf(stderr, "%.0f,", BBt[i_row][i_col]);
			} // i_row
		} // i_col

		fprintf(stderr, "\nInverting BBt:\n");
		double** inv_BBt = invert_matrix_GJ(BBt, n_ind_rows, n_ind_rows, NULL);

		// Test inversion.
		double** BBt_inv_BBt = matrix_multiply(BBt, n_ind_rows, n_ind_rows, inv_BBt, n_ind_rows, n_ind_rows, NULL);
		double ident_eps = pow(10, -3);
		if (!matrix_is_identity(BBt_inv_BBt, n_ind_rows, n_ind_rows, ident_eps))
		{
			fprintf(stderr, "Sanity check failed: Inversion error..\n");
			exit(1);
		}

		// Now test the right inverse.
		fprintf(stderr, "Calculating right inverse..\n");
		double** invr_B = matrix_multiply(Bt, n_panel_haplotypes, n_ind_rows, inv_BBt, n_ind_rows, n_ind_rows, NULL);

		fprintf(stderr, "Testing right inverse..\n");
		double** ident_check_mat = matrix_multiply(allele1_2_hap_ind, n_ind_rows, n_panel_haplotypes, invr_B, n_panel_haplotypes, n_ind_rows, NULL);
		if (!matrix_is_identity(ident_check_mat, n_ind_rows, n_ind_rows, ident_eps))
		{
			fprintf(stderr, "Right inverse is not correctly computed.\n");
			exit(1);
		}

		fprintf(stderr, "Success..\n");

		exit(1);
	} // -check_singular_allele1_2_haplotype_matrix_right_inversion option.
	else if (t_string::compare_strings(argv[1], "-generate_untyped_variant_recoded_reference_panel"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Ref. panel haplocoded tag matbed] [Ref. panel haplocoded target matbed] \
[Panel sample list file] [Output Directory] [Output prefix]\n", argv[0], argv[1]);

			exit(1);
		}

		char* ref_panel_tag_haplocoded_matbed_fp = argv[2];
		char* ref_panel_target_haplocoded_matbed_fp = argv[3];
		char* panel_sample_list_fp = argv[4];
		char* op_dir = argv[5];
		char* op_prefix = argv[6];

		generate_untyped_variant_recoded_reference_panel(ref_panel_tag_haplocoded_matbed_fp, ref_panel_target_haplocoded_matbed_fp, panel_sample_list_fp,
			op_dir, op_prefix);
	} // -generate_untyped_variant_recoded_reference_panel option.
	else if (t_string::compare_strings(argv[1], "-calculate_proxy2clear_pairwise_distance_stats"))
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s %s [Proxy panel matbed file] [Proxy panel sample list file] \
[Cleartext panel matbed file] [Cleartext sample list file] \
[Cleartext windowizing win. length] [Proxy windowizing win. length] [Length of correlation window] \
[Output file path]\n", argv[0], argv[1]);

			exit(1);
		}

		char* proxized_panel_matbed = argv[2];
		char* proxized_panel_sample_list_fp = argv[3];
		char* cleartext_panel_matbed = argv[4];
		char* cleartext_panel_sample_list_fp = argv[5];
		int l_cleartext_windowizing_win = atoi(argv[6]);
		int l_proxy_windowizing_win = atoi(argv[7]);
		int l_corr_win = atoi(argv[8]);
		char* op_prefix = argv[9];

		calculate_proxy2clear_pairwise_distance_stats(proxized_panel_matbed, proxized_panel_sample_list_fp,
												cleartext_panel_matbed, cleartext_panel_sample_list_fp, 
												l_cleartext_windowizing_win, l_proxy_windowizing_win,
												l_corr_win,
												op_prefix);
	} // -calculate_Homer_LRT_statistics_on_proxized_panels option.
	else if (t_string::compare_strings(argv[1], "-calculate_windowed_stats_Homer_t_statistics_on_proxized_panels"))
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s %s [Original query panel] [Query sample list] \
[Reference panel] [Reference sample list] \
[Database panel] [Database sample list] \
[# vicinity variants] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* target_query_panel_matbed = argv[2];
		char* target_query_panel_sample_list_fp = argv[3];
		char* ref_AF_panel_matbed = argv[4];
		char* ref_AF_panel_sample_list_fp = argv[5];
		char* target_AF_database_panel_matbed = argv[6];
		char* target_AF_database_panel_sample_list_fp = argv[7];
		int n_vicinity = atoi(argv[8]);
		char* op_prefix = argv[9];

		calculate_windowed_Homer_t_statistics_on_proxized_panels(target_query_panel_matbed, target_query_panel_sample_list_fp,
			ref_AF_panel_matbed, ref_AF_panel_sample_list_fp, // How can an attacker have access to proxized reference panel???
			target_AF_database_panel_matbed, target_AF_database_panel_sample_list_fp,
			n_vicinity,
			op_prefix);
	}
	else if (t_string::compare_strings(argv[1], "-pool_summarize_Homer_t_statistics_per_query"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Per query per variant Homer-t Dyij stats file] [Query sample list] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* per_query_per_variant_files_list_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* summarized_results_op_fp = argv[4];

		pool_summarize_Homer_t_statistics_per_query(per_query_per_variant_files_list_fp, sample_ids_list_fp, summarized_results_op_fp);
	} // -pool_summarize_Homer_t_statistics_per_query option.
	else if (t_string::compare_strings(argv[1], "-calculate_Homer_t_statistics_on_proxized_panels"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Original query panel] [Query sample list] \
[Reference panel] [Reference sample list] \
[Database panel] [Database sample list] \
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_panel_matbed = argv[2];
		char* query_panel_sample_list_fp = argv[3];
		char* ref_panel_matbed = argv[4];
		char* ref_panel_sample_list_fp = argv[5];
		char* database_panel_matbed = argv[6];
		char* database_panel_sample_list_fp = argv[7];
		char* op_prefix = argv[8];

		calculate_Homer_t_statistics_on_proxized_panels(query_panel_matbed, query_panel_sample_list_fp,
															ref_panel_matbed, ref_panel_sample_list_fp, 
															database_panel_matbed, database_panel_sample_list_fp,
															op_prefix);
	} // -calculate_Homer_LRT_statistics_on_proxized_panels option.
	else if (t_string::compare_strings(argv[1], "-pool_summarize_Sankararaman_LRT_statistics_per_query"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [Per query per variant Homer-t Dyij stats file] [Query sample list] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* per_query_per_variant_files_list_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* summarized_results_op_fp = argv[4];

		pool_summarize_Sankararaman_LRT_statistics_per_query(per_query_per_variant_files_list_fp, sample_ids_list_fp, summarized_results_op_fp);
	} // -pool_summarize_Homer_t_statistics_per_query option.
	else if (t_string::compare_strings(argv[1], "-write_securegenome_input_files"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Original query panel] [Query sample list] \
[Reference panel] [Reference sample list] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_panel_matbed = argv[2];
		char* query_panel_sample_list_fp = argv[3];
		char* ref_panel_matbed = argv[4];
		char* ref_panel_sample_list_fp = argv[5];
		char* op_prefix = argv[6];

		write_securegenome_input_files(query_panel_matbed, query_panel_sample_list_fp,
			ref_panel_matbed, ref_panel_sample_list_fp,
			op_prefix);
	} // -write_securegenome_input_files option.
	else if (t_string::compare_strings(argv[1], "-calculate_Sankararaman_LRT_statistics_on_proxized_panels"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Original query panel] [Query sample list] \
[Reference panel] [Reference sample list] \
[Database panel] [Database sample list] \
[MAF cutoff] \
[Min var-var dist] \
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_panel_matbed = argv[2];
		char* query_panel_sample_list_fp = argv[3];
		char* ref_panel_matbed = argv[4];
		char* ref_panel_sample_list_fp = argv[5];
		char* database_panel_matbed = argv[6];
		char* database_panel_sample_list_fp = argv[7];
		double maf_cutoff = atof(argv[8]);
		int min_var2var_dist = atoi(argv[9]);
		char* op_prefix = argv[10];

		calculate_Sankararaman_LRT_statistics_on_proxized_panels(query_panel_matbed, query_panel_sample_list_fp,
			ref_panel_matbed, ref_panel_sample_list_fp,
			database_panel_matbed, database_panel_sample_list_fp,
			maf_cutoff,
			min_var2var_dist,
			op_prefix);
	} // -calculate_Homer_LRT_statistics_on_proxized_panels option.
	else if (t_string::compare_strings(argv[1], "-modify_per_site_mixing_parameters"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Parameters file] [# 1st to mod.] [# 2nd to mod.] [# 3rd to mod.] [modif params file]\n", argv[0], argv[1]);
			exit(1);
		}

		char* parameter_op_fp = argv[2];
		int n_first_degree_modif = atoi(argv[3]);
		int n_sec_degree_modif = atoi(argv[4]);
		int n_third_degree_modif = atoi(argv[5]);
		char* modif_parameter_op_fp = argv[6];

		fprintf(stderr, "Modifying %d first degree and %d second degree parameters from per-variant parameter file %s and saving to %s\n", 
				n_first_degree_modif, n_sec_degree_modif, parameter_op_fp, modif_parameter_op_fp);

		fprintf(stderr, "Loading mixing parameters per variant..\n");
		int loaded_n_vicinity = 0;
		int loaded_avg_geno_mod = 0;
		vector<t_annot_region*>* loaded_per_var_proxy_params = load_per_variant_site_mixing_parameters(parameter_op_fp, loaded_n_vicinity, loaded_avg_geno_mod);

		fprintf(stderr, "Loaded proxization parameters for %d variants over [-%d, +%d] vicinity using %d-modular arithmetic.\n",
			(int)loaded_per_var_proxy_params->size(),
			loaded_n_vicinity, loaded_n_vicinity, loaded_avg_geno_mod);

		int n_modif_params = 0;
		for (int i_var = 0; i_var < (int)loaded_per_var_proxy_params->size(); i_var++)
		{
			void** cur_var_loaded_reg_info = (void**)(loaded_per_var_proxy_params->at(i_var)->data);
			int* loaded_per_var_weights = (int*)(cur_var_loaded_reg_info[2]);
			int** loaded_per_var2var_weights = (int**)(cur_var_loaded_reg_info[3]);
			int*** loaded_per_var2var2var_weights = (int***)(cur_var_loaded_reg_info[4]);

			int cur_var_n_modif_first_deg = 0;
			int cur_var_n_modif_sec_deg = 0;
			int cur_var_n_modif_third_deg = 0;

			for (int vic_i = 0; vic_i < 2 * loaded_n_vicinity + 1; vic_i++)
			{
				if (cur_var_n_modif_first_deg >= n_first_degree_modif)
				{
					break;
				}

				if (loaded_per_var_weights[vic_i] > 0)
				{
					loaded_per_var_weights[vic_i] = 0;
					cur_var_n_modif_first_deg++;
					n_modif_params++;
				}
			} // vic_i loop.

			for (int vic_i = 0; vic_i < 2 * loaded_n_vicinity + 1; vic_i++)
			{
				for (int vic_j = 0; vic_j < 2 * loaded_n_vicinity + 1; vic_j++)
				{
					if (cur_var_n_modif_sec_deg >= n_sec_degree_modif)
					{
						break;
					}

					if (loaded_per_var2var_weights[vic_i][vic_j] > 0)
					{
						n_modif_params++;
						loaded_per_var2var_weights[vic_i][vic_j] = 0;
						cur_var_n_modif_sec_deg++;
					}
				} // vic_j loop.
			} // vic_i loop.

			for (int vic_i = 0; vic_i < 2 * loaded_n_vicinity + 1; vic_i++)
			{
				for (int vic_j = 0; vic_j < vic_i; vic_j++)
				{
					for (int vic_k = 0; vic_k < vic_j; vic_k++)
					{
						if (cur_var_n_modif_third_deg >= n_third_degree_modif)
						{
							break;
						}

						if (loaded_per_var2var2var_weights[vic_i][vic_j][vic_k] > 0)
						{
							n_modif_params++;
							loaded_per_var2var2var_weights[vic_i][vic_j][vic_k] = 0;
							cur_var_n_modif_third_deg++;
						}
					}
				} // vic_j loop.
			} // vic_i loop.
		} // i_var loop.

		fprintf(stderr, "Modified %d parameters on %d variants.\n", n_modif_params, (int)loaded_per_var_proxy_params->size());

		save_per_variant_site_mixing_parameters(loaded_per_var_proxy_params, loaded_n_vicinity, loaded_avg_geno_mod, modif_parameter_op_fp);
	} // -modify_per_site_mixing_parameters option.
	else if (t_string::compare_strings(argv[1], "-generate_save_per_site_mixing_parameters_LD_aware"))
	{
		if (argc != 13)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file (BED)] \
[Vicinity size] \
[Variant Weight Prob.] [Var-var interaction prob.] [Var-Var-Var interaction prob.] [Weight inversion Prob.] \
[Allele Coding Modulus] [Norm. N_e] [Min. # of parameters per variant] [Genetic distance map directory] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* target_vars_BED_fp = argv[2];
		int n_vicinity = atoi(argv[3]);
		double per_var_weight_prob = atof(argv[4]);
		double per_var2var_interaction_prob = atof(argv[5]);
		double per_var2var2var_interaction_prob = atof(argv[6]);
		double weight_inversion_prob = atof(argv[7]);
		int avg_geno_mod = atoi(argv[8]);
		double normalized_N_e = atof(argv[9]);
		int min_n_params_per_var = atoi(argv[10]);
		char* recombination_rate_dir = argv[11];
		char* parameter_op_fp = argv[12];

		t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

		vector<t_annot_region*>* haplocoded_var_regs = load_BED(target_vars_BED_fp);

		generate_save_per_site_mixing_parameters_LD_aware(haplocoded_var_regs,
			rng,
			n_vicinity,
			per_var_weight_prob, 
			weight_inversion_prob,
			per_var2var_interaction_prob,
			per_var2var2var_interaction_prob,
			avg_geno_mod,
			normalized_N_e,
			min_n_params_per_var,
			recombination_rate_dir,
			parameter_op_fp);
	} // -generate_save_per_site_mixing_parameters_LD_aware option.
	else if (t_string::compare_strings(argv[1], "-generate_save_per_site_mixing_parameters"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Vicinity size] [Variant Weight Prob.] [Var-var interaction prob.] [Allele Coding Modulus] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		int n_vicinity = atoi(argv[4]);
		double per_var_weight_prob = atof(argv[5]);
		double per_var2var_interaction_prob = atof(argv[6]);
		int avg_geno_mod = atoi(argv[7]);
		char* parameter_op_fp = argv[8];

		t_rng* rng = new t_rng(t_seed_manager::seed_me_getrandom());

		vector<t_annot_region*>* haplocoded_var_regs = load_variant_signal_regions_wrapper(haplocoded_geno_signal_fp, sample_ids_list_fp);
		//vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		generate_save_per_site_mixing_parameters(haplocoded_var_regs, 
												rng, 
												n_vicinity, per_var_weight_prob, per_var2var_interaction_prob, avg_geno_mod, 
												parameter_op_fp);

		// Now load and test.
		fprintf(stderr, "Loading mixing parameters per variant..\n");
		int loaded_n_vicinity = 0;
		int loaded_avg_geno_mod = 0;
		vector<t_annot_region*>* loaded_per_var_proxy_params = load_per_variant_site_mixing_parameters(parameter_op_fp, loaded_n_vicinity, loaded_avg_geno_mod);

		fprintf(stderr, "Loaded proxization parameters for %d variants over [-%d, +%d] vicinity using %d-modular arithmetic.\n", 
			(int)loaded_per_var_proxy_params->size(),
			loaded_n_vicinity, loaded_n_vicinity, loaded_avg_geno_mod);

		for (int i_var = 0; i_var < 10; i_var++)
		{
			fprintf(stderr, "%s:%d-%d loaded parameters:\n", 
				loaded_per_var_proxy_params->at(i_var)->chrom, 
				loaded_per_var_proxy_params->at(i_var)->start, 
				loaded_per_var_proxy_params->at(i_var)->end);

			void** reg_info = (void**)(loaded_per_var_proxy_params->at(i_var)->data);
			int* per_var_weight = (int*)(reg_info[2]);
			int** per_var2var_interaction_weight = (int**)(reg_info[3]);

			fprintf(stderr, "Per variant weights:\n");
			for(int vic_i = 0; vic_i < 2*loaded_n_vicinity+1; vic_i++)
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
	} // -generate_save_per_site_mixing_parameters option.
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Reference sample list file] [Variant start index] [Variant end index] [Half vicinity size (i.e., k)] [Output prefix]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_sample_ids_fp = argv[6];
		int var_start_i = atoi(argv[7]);
		int var_end_i = atoi(argv[8]);
		int n_vicinity = atoi(argv[9]);
		char* op_prefix = argv[10];

		decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, 
			query_sample_ids_fp,
			ref_original_haplocoded_geno_fp, ref_sample_ids_fp,
			var_start_i, var_end_i,
			n_vicinity,
			op_prefix);
	} // -decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM option.
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Reference sample list file] [Half vicinity size (i.e., k)]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_sample_ids_fp = argv[6];
		int n_vicinity = atoi(argv[7]);

		decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
			ref_original_haplocoded_geno_fp, ref_sample_ids_fp,
			n_vicinity);
	}
	else if (t_string::compare_strings(argv[1], "-decode_site_alleles_per_proxized_reference"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Query original haplocoded matbed file] [Query proxied haplocoded matbed file] [Query sample list file] \
[Reference original haplocoded matbed file] [Ref. proxized haplocoded matbed file] [Reference sample list file] [Half vicinity size (i.e., k)]\n", argv[0], argv[1]);
			exit(1);
		}

		char* query_original_haplocoded_geno_fp = argv[2];
		char* query_proxized_haplocoded_geno_fp = argv[3];
		char* query_sample_ids_fp = argv[4];
		char* ref_original_haplocoded_geno_fp = argv[5];
		char* ref_proxized_haplocoded_geno_fp = argv[6];
		char* ref_sample_ids_fp = argv[7];
		int n_vicinity = atoi(argv[8]);

		decode_site_alleles_per_proxized_reference(query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
													ref_original_haplocoded_geno_fp, ref_proxized_haplocoded_geno_fp, ref_sample_ids_fp,
													n_vicinity);
	} // -decode_site_alleles_per_proxized_reference option.
	else if (t_string::compare_strings(argv[1], "-proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Per variant vicinity weight params file path] [Allele error probability] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* per_var_weight_params_fp = argv[4];
		double allele_err_eps = atof(argv[5]);
		char* proxized_haplocoded_var_regs_op_fp = argv[6];

		proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(haplocoded_geno_signal_fp, sample_ids_list_fp,
			per_var_weight_params_fp,
			allele_err_eps,
			proxized_haplocoded_var_regs_op_fp);
	} // -proxize_variants_per_vicinity_non_linear_modular_average option.
	else if (t_string::compare_strings(argv[1], "-MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Per variant vicinity weight params file path] [Allele error probability] [# threads] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* per_var_weight_params_fp = argv[4];
		double allele_err_eps = atof(argv[5]);
		int n_threads = atoi(argv[6]);
		char* proxized_haplocoded_var_regs_op_fp = argv[7];

		MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters(haplocoded_geno_signal_fp, sample_ids_list_fp,
			per_var_weight_params_fp,
			allele_err_eps,
			n_threads,
			proxized_haplocoded_var_regs_op_fp);
	} // -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters
	else if (t_string::compare_strings(argv[1], "-MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Per variant vicinity weight params file path] [Allele error probability] [# threads] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* per_var_weight_params_fp = argv[4];
		double allele_err_eps = atof(argv[5]);
		int n_threads = atoi(argv[6]);
		char* proxized_haplocoded_var_regs_op_fp = argv[7];

		MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer(haplocoded_geno_signal_fp, sample_ids_list_fp,
			per_var_weight_params_fp,
			allele_err_eps,
			n_threads,
			proxized_haplocoded_var_regs_op_fp);
	} // -MT_proxize_variants_per_vicinity_non_linear_modular_average_per_var_custom_filters_no_buffer
	else if (t_string::compare_strings(argv[1], "-proxize_variants_per_vicinity_non_linear_modular_average_uniform"))
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Vicinity size] [Weight params file path] [Variant Weight Prob.] [Var-var interaction prob.] [Allele Coding Modulus] [Allele error probability] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		int n_vicinity = atoi(argv[4]);
		char* weight_params_fp = argv[5];
		double per_var_weight_prob = atof(argv[6]);
		double per_var2var_interaction_prob = atof(argv[7]);
		int avg_geno_mod = atoi(argv[8]);
		double allele_err_eps = atof(argv[9]);
		char* proxized_haplocoded_var_regs_op_fp = argv[10];

		proxize_variants_per_vicinity_non_linear_modular_average_uniform(haplocoded_geno_signal_fp, sample_ids_list_fp,
			n_vicinity,
			weight_params_fp,
			per_var_weight_prob,
			per_var2var_interaction_prob,
			avg_geno_mod,
			allele_err_eps,
			proxized_haplocoded_var_regs_op_fp);
	} // -proxize_variants_per_vicinity_non_linear_modular_average option.
	else if (t_string::compare_strings(argv[1], "-proxize_variants_per_vicinity_modular_average"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s %s [Haplocoded var regs file] [Sample list file] \
[Vicinity size] [Variant Weight] [Allele Coding Modulus] [Allele error probability] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* haplocoded_geno_signal_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		int n_vicinity = atoi(argv[4]);
		int per_var_weight = atoi(argv[5]);
		int avg_geno_mod = atoi(argv[6]);
		double allele_err_eps = atof(argv[7]);
		char* proxized_haplocoded_var_regs_op_fp = argv[8];

		proxize_variants_per_vicinity_modular_average(haplocoded_geno_signal_fp, sample_ids_list_fp,
			n_vicinity, per_var_weight,
			avg_geno_mod,
			allele_err_eps,
			proxized_haplocoded_var_regs_op_fp);
	} // -proxize_variants_per_vicinity_modular_average option.
	else if (t_string::compare_strings(argv[1], "-get_hapfreq_confusion_matrix"))
	{
		if (argc != 10)
		{
			fprintf(stderr, "USAGE: %s %s [Reference1 matbed file] [Reference1 sample id's] \
[Reference2 matbed file] [Reference2 sample id's] \
[# variants per window] \
[Min. HF ref1 hap. freq. to compare] [Max. ref2 hap. freq. to compare] \
[Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* ref1_panel_matbed_fp = argv[2];
		char* ref1_sample_list_fp = argv[3];
		char* ref2_panel_matbed_fp = argv[4];
		char* ref2_sample_list_fp = argv[5];
		int n_vars_per_win = atoi(argv[6]);
		double min_HF_ref1_freq_2_quantify = atof(argv[7]);
		double max_LF_ref2_freq_2_quantify = atof(argv[8]);
		char* op_fp = argv[9];

		get_hapfreq_confusion_matrix(ref1_panel_matbed_fp, ref1_sample_list_fp,
			ref2_panel_matbed_fp, ref2_sample_list_fp,
			n_vars_per_win,
			min_HF_ref1_freq_2_quantify,
			max_LF_ref2_freq_2_quantify,
			op_fp);
	} // -haplotype_frequency_attack option.
	else if (t_string::compare_strings(argv[1], "-compare_per_win_haplotypes"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "USAGE: %s %s [Reference1 matbed file] [Reference1 sample id's] \
[Reference2 matbed file] [Reference2 sample id's] \
[# variants per window] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* ref1_panel_matbed_fp = argv[2];
		char* ref1_sample_list_fp = argv[3];
		char* ref2_panel_matbed_fp = argv[4];
		char* ref2_sample_list_fp = argv[5];
		int n_vars_per_win = atoi(argv[6]);
		char* op_fp = argv[7];

		compare_per_win_haplotypes(ref1_panel_matbed_fp, ref1_sample_list_fp, 
									ref2_panel_matbed_fp, ref2_sample_list_fp, 
									n_vars_per_win, op_fp);
	} // -compare_per_win_haplotypes option.
	else if (t_string::compare_strings(argv[1], "-get_per_win_haplotype_frequencies_per_reference"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Reference matbed file] [Referene sample id's] [# variants per window] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* ref_panel_matbed_fp = argv[2];
		char* ref_sample_list_fp = argv[3];
		int n_vars_per_win = atoi(argv[4]);
		char* op_fp = argv[5];

		get_per_win_haplotype_frequencies_per_reference(ref_panel_matbed_fp, ref_sample_list_fp, n_vars_per_win, op_fp);
	} // get_per_window_sampling_probability_per_reference option.
	else if (t_string::compare_strings(argv[1], "-convert_BEAGLE_imputed_meta_panel_2_matbed"))
	{
		if (argc != 6)
		{
			fprintf(stderr ,"USAGE: %s %s [VCF file path] [Sample id list fp] [Imputed genotype randomization weight] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* vcf_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		double imp_geno_rand_weight = atof(argv[4]);
		char* op_fp = argv[5];

		convert_BEAGLE_imputed_meta_panel_2_matbed(vcf_fp, sample_ids_list_fp, imp_geno_rand_weight, op_fp);
	} // -convert_BEAGLE_imputed_meta_panel_2_matbed option.
	else if (t_string::compare_strings(argv[1], "-summarize_sampled_segments_per_resampled_haplotype_info"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s %s [Resampling haplotype info sigbed file path] [Resampling sample list file] \
[Generating sample identifiers list] [Output file path]\n", argv[0], argv[1]);
			exit(1);
		}

		char* resampling_hap_info_sigbed_fp = argv[2];
		char* resampling_sample_list_fp = argv[3];
		char* generating_sample_list_fp = argv[4];
		char* op_fp = argv[5];

		summarize_sampled_segments_per_resampled_haplotype_info(resampling_hap_info_sigbed_fp, resampling_sample_list_fp, generating_sample_list_fp, op_fp);
	} // summarize_sampled_segments_per_resampled_haplotype_info option.	
	
	clock_t end_c = clock();
	time_t end_t = time(NULL);
	fprintf(stderr, "%s finished option \"%s\" in %d (%d) seconds.\n", argv[0], argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));

	FILE* f_beacon = open_f("beacon.log", "a");
	fprintf(f_beacon, "%s finished option \"%s\" in %d (%d) seconds.\n", argv[0], argv[1], (int)(end_t - start_t), (int)((end_c - start_c) / CLOCKS_PER_SEC));
	close_f(f_beacon, NULL);
}