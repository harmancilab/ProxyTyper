#ifndef __VCF__
#define __VCF__

#include <vector>
using namespace std;

struct t_annot_region;
struct t_restr_annot_region_list;

class t_string;

void extract_hapgen2_files_per_single_reference(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp);

void convert_HAP_LEGEND_2_haplocoded_matbed(char* hap_fp, char* legend_fp, char* chrom_str, bool add_AF_info_2_id, char* output_matbed_fp);

void convert_GEN_2_matbed(char* gen_fp, char* chrom_str, char* output_matbed_fp);

void extract_hapgen2_files(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp);

void replace_variant_alleles_per_genotypes_signals(char* geno_regs_BED_fp, char* sample_ids_list_fp, char* new_alleles_BED_fp, char* op_fp);

bool is_haplocoded_allele_het(char haplocoded_geno);

void add_noise_to_geno_signal_regions(char* geno_sig_regs_fp, char* sample_ids_list_fp, double allele_error_prob, char* op_fp);

void concat_variant_wide_genotype_signal_regions(char* matbed1_fp, char* sample1_list_fp, char* matbed2_fp, char* sample2_list_fp, bool remove_unique_coord_vars_flag, char* op_fp);

void concat_variant_wide_genotype_signal_regions_per_geno_sig_list(char* matbed_file_list_fp, char* sample_list_fp, bool remove_unique_coord_vars_flag, char* op_fp);

void uniquefy_genotype_signal_regions(char* geno_reg_sig_file, char* sample_list_fp, char* op_fp);

void load_dump_allele_switching_events_in_haplotypes(char* haplo_swtiching_geno_matrix_fp, char* sample_ids_list_fp);

int get_max_genotype_value(vector<t_annot_region*>* regs, vector<char*>* sample_ids);

int get_genotype_per_haplocoded_genotype(char haplocoded_geno);
int get_allele_per_haplotype(char geno, int hap_i);
int get_allele_haploswitch_per_haplotype(char geno, int hap_i, char& haplo_switch);

void dump_haplo_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, int haplo_2_extract, const char* op_fp);

vector<t_annot_region*>* load_recombination_rate_regions_per_map_file(char* recombination_rate_dir);

double get_recomb_prob_per_delta_cM(double delta_cM);

double get_recomb_rate_difference_per_regions(t_restr_annot_region_list* restr_recomb_rate_regs, t_annot_region* reg1, t_annot_region* reg2);

void filter_redundant_variants_per_genotype_distance(char* genotype_matbed_fp, char* sample_ids_list_fp,
	double max_geno_distance, int max_index_difference,
	char* unique_regs_BED_op_fp);

// Save/Load regions+subjects+binary_matrix_genotypes data for loading genotypes into python/R, etc.
void binarize_variant_genotype_signal_regions_per_matrix_subjects_regions(vector<t_annot_region*>* genotype_signal_regions, vector<char*>* geno_sample_ids,
	const bool compress_geno_matrix,
	const char* op_prefix);

vector<t_annot_region*>* load_variant_genotype_signal_regions_per_matrix_subjects_regions(const char* op_prefix, vector<char*>* geno_sample_ids);

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids);
vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids);
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp);
void binarize_variant_signal_regions_wrapper(vector<t_annot_region*>* roi_regs_w_signals, vector<char*>* sample_ids, char* op_fp);
void binarize_variant_genotype_signal_regions_per_matrix_subjects_regions_multithreaded(vector<t_annot_region*>* genotype_signal_regions, vector<char*>* geno_sample_ids,
	const bool compress_geno_matrix,
	int n_threads,
	const char* op_prefix);

unsigned char* load_MatSubjReg_unstructured(char* matsubjreg_prefix, vector<char*>* subj_ids, vector<t_annot_region*>* var_regs);
void save_MatSubjReg_unstructured(vector<t_annot_region*>* var_regs, vector<char*>* subj_ids, const unsigned char* geno_mat, bool compress, char* op_prefix);
vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp);


void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_regs_only, const char* op_fp);

void extract_genotype_signals_per_chr_ids(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* chr_ids_list_fp, char* op_fp);
void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp);
void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp);

void convert_haplocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, bool phase_op, char* op_VCF_fp);

void convert_genocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* op_VCF_fp);

void convert_haplocoded_signal_regs_2_VCF_w_header(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* header_file, bool phase_op, char* op_VCF_fp);

void convert_haplocoded_signal_regs_2_VCF_w_header_multithreaded(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* header_file, bool phase_op, int n_threads, char* op_VCF_fp);

char vcf_entry2geno(const char* vcf_entry, bool haplocoded);
void geno2vcf_entry(const char geno, char* vcf_geno_str, bool haplocoded);

void extract_genotype_signals_per_VCF_memoptimized(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* var_regions_BED_fp,			// Regions to focus on while extracting.
	char* chr_id_2_process,				// Chromosome to process.
	char* bin_seq_dir,					// Sequences directory.
	bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
	bool match_region_names_flag,		// This is a flag to enforce matching of region names.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	char* op_fp);						// Output file path.

void extract_genotype_signals_per_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* var_regions_BED_fp,			// Regions to focus on while extracting.
	char* chr_id_2_process,				// Chromosome to process.
	char* bin_seq_dir,					// Sequences directory.
	bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
	bool match_region_names_flag,		// This is a flag to enforce matching of region names.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	char* op_fp);						// Output file path.

void generate_random_sampled_genotype_matrix(char* genotype_regs_BED_fp, char* sample_ids_list_fp,
	int n_samples_2_generate,
	char* op_fp);

void merge_samples_per_regions_list(char* regions_list_fp,
	char* samples_list_fp,
	bool match_region_names,
	char* regions_op_fp,
	char* samples_op_fp);

void merge_samples_per_matching_genotype_signal_regions(char* sample1_regs_fp, char* sample1_list_fp,
														char* sample2_regs_fp, char* sample2_list_fp, bool match_region_names, char* op_fp);

class t_rng;
struct t_annot_region;

char* copy_nucs(char* seq);
int* copy_bps(char* seq, int* bps);
void get_mutation_matrices_per_DAF(char* snp_vcf_fp, double min_daf, double max_daf);
char mutate_nuc_per_transition_transversion(char nuc);
char mutate_nuc(char nuc, t_rng* rng);
char mutate_nuc_2_GC(char nuc, t_rng* rng);
bool mutate_GC_2_AU(char* seq);
char* get_random_seq(int l);


#endif // __VCF__