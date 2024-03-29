#ifndef __VCF__
#define __VCF__

#include <vector>
using namespace std;

struct t_annot_region;
struct t_restr_annot_region_list;

class t_string;

struct t_genotype
{
	char mat_allele;
	char pat_allele;
	int quality;
	bool phased;
};

/*
TODO::All variation infor must be encapsulated to t_variant_info structure that is pointed to by the annotation_info member.
*/
struct t_vcf_info
{
	int ref_pos;
	char* info_str;
	char* ref_allele_str;
	char* alt_allele_str;

	// Genotype information.
	char* genotype_format;
	vector<t_genotype*>* genotypes;
};

vector<t_annot_region*>* load_recombination_rates(char* recombination_rate_fp);

double get_avg_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs);

double get_cumulative_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs);

double get_cumulative_recomb_rate_per_variant(t_annot_region* var_reg, vector<t_annot_region*>* recomb_regs);

void harmonize_dbGAP_GenoPheno_files(char* vcf_fp, char* vcf_samples_list_fp, char* genome_dir, char* recoded_vcf_fp);

void harmonize_dbGAP_GenoPheno_files_forward_only_check(char* vcf_fp, char* vcf_samples_list_fp, char* genome_dir, char* recoded_vcf_fp);


void extract_hapgen2_files(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp);


void replace_variant_alleles_per_genotypes_signals(char* geno_regs_BED_fp, char* sample_ids_list_fp, char* new_alleles_BED_fp, char* op_fp);

bool is_haplocoded_allele_het(char haplocoded_geno);

void concat_variant_wide_genotype_signal_regions(char* matbed1_fp, char* sample1_list_fp, char* matbed2_fp, char* sample2_list_fp, bool remove_unique_coord_vars_flag, char* op_fp);

void load_dump_allele_switching_events_in_haplotypes(char* haplo_swtiching_geno_matrix_fp, char* sample_ids_list_fp);

int get_max_genotype_value(vector<t_annot_region*>* regs, vector<char*>* sample_ids);

int get_genotype_per_haplocoded_genotype(char haplocoded_geno);
int get_allele_per_haplotype(char geno, int hap_i);
int get_allele_haploswitch_per_haplotype(char geno, int hap_i, char& haplo_switch);

void dump_haplo_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, int haplo_2_extract, const char* op_fp);

vector<t_annot_region*>* load_recombination_rate_regions_per_map_file(char* recombination_rate_dir);

double get_recomb_prob_per_delta_cM(double delta_cM);

double get_recomb_rate_difference_per_regions(t_restr_annot_region_list* restr_recomb_rate_regs, t_annot_region* reg1, t_annot_region* reg2);

struct t_TF_motif_effect_info
{
	// The effected TF motif region.
	char* ref_motif;
	char* alt_motif;
	double alt_2_ref_binding_affinity_diff;
};

struct t_miRNA_binding_site_effect_info
{
	// Affected miRNA binding site.
	char* ref_miRNA_binding_site;
	char* alt_miRNA_binding_site;

	double alt_2_ref_binding_affinity_diff;
};

// This effect info marks effect of a variants on a specific element.
struct t_variant_effect_info
{
	// Each effect has a type.
	int effect_type;

	// Each effect associates with a functional element type.
	int effected_element_type;

	// Each effect overlaps and effects a genomic element.
	t_annot_region* effected_element_region;

	// Each effect has a element specific effect information.
	void* element_specific_effect_info;
};

struct t_variant_info
{
	// Every variant has a type.
	int var_type;

	// Every variant has a ref and alt allele.
	char* ref_allele;
	char* alt_allele;
	char* neighbor_seq;

	// Every variant can be assigned a list of effects.
	vector<t_variant_effect_info*>* variant_effects;

	// Every variant can have its own data structure.
	void* var_data;
};

char** get_variant_effect_name_strings_list();

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids);
vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids);
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp);

vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp);

void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_regs_only, const char* op_fp);

void extract_genotype_signals_per_chr_ids(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* chr_ids_list_fp, char* op_fp);
void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp);
void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp);

char get_nuc_per_transcript_i(t_annot_region* transcript, int nuc_trans_i);

void convert_haplocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, bool phase_op, char* op_VCF_fp);

void convert_genocoded_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, char* op_VCF_fp);

// This function adds variant effect information for a list of candidate variants.
void add_variant_effect_information(vector<t_annot_region*>* candidate_var_regs, 
									vector<t_annot_region*>* gencode_transcript_regs,
									vector<t_annot_region*>* gencode_cds_regs,
									vector<t_annot_region*>* gencode_5p_UTR_regs,
									vector<t_annot_region*>* gencode_3p_UTR_regs,
									vector<t_annot_region*>* gencode_promoter_regs,
									vector<t_annot_region*>* gencode_intronic_regs,
									vector<t_annot_region*>* encode_motif_regs);

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

vector<t_annot_region*>* load_VCF_regions(char* vcf_fp, bool parse_genotypes);

// All the information on the line.
struct t_vcf_entry
{
	// Chromosome identifier.
	char* chrom; // This is an id, not a number. 

	// The reference index position.
	int ref_pos;

	// identifier for this vcf entry.
	char* id;

	// All the info string for this entry.
	char* info_str;

	// Genotype format string.
	char* genotype_format;

	// The reference and alternative alleles.
	t_string* ref_allele;
	vector<t_string*>* alt_alleles;

	// These are the paternal and maternal alleles for this variant that are annotated in the file.
	vector<t_genotype*>* genotypes;
};


class t_rng;
struct t_annot_region;
struct t_vcf_entry;

char* copy_nucs(char* seq);
int* copy_bps(char* seq, int* bps);
void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count);
void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char* ancestral_allele,
	char* derived_allele,
	int& alternate_allele_count,
	int& total_allele_count);
bool verify_snp(t_vcf_entry* snp, char* ncrna_nucs, char ncrna_strand, int snp_i_seq);
void parse_variants_per_info_string(char* info_str,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count);
void get_mutation_matrices_per_DAF(char* snp_vcf_fp, double min_daf, double max_daf);
char mutate_nuc_per_transition_transversion(char nuc);
char mutate_nuc(char nuc, t_rng* rng);
char mutate_nuc_2_GC(char nuc, t_rng* rng);
bool mutate_GC_2_AU(char* seq);
bool get_DAF_per_snp_region(t_annot_region* snp_region, double& daf);
void length_based_Duplex_deltaG_simulation();
char* get_random_seq(int l);
//
//struct t_genotype
//{
//	char mat_allele;
//	char pat_allele;
//	int quality;
//	bool phased;
//};
//
//struct t_vcf_info
//{
//	int ref_pos;
//	char* info_str;
//	char* ref_allele_str;
//	char* alt_allele_str;
//
//	// Genotype information.
//	char* genotype_format;
//	vector<t_genotype*>* genotypes;
//};

vector<t_annot_region*>* load_VCF_regions(char* vcf_fp, bool parse_genotypes);

//// All the information on the line.
//struct t_vcf_entry
//{
//	// Chromosome identifier.
//	char* chrom; // This is an id, not a number. 
//
//				 // The reference index position.
//	int ref_pos;
//
//	// identifier for this vcf entry.
//	char* id;
//
//	// All the info string for this entry.
//	char* info_str;
//
//	// Genotype format string.
//	char* genotype_format;
//
//	// The reference and alternative alleles.
//	t_string* ref_allele;
//	vector<t_string*>* alt_alleles;
//
//	// These are the paternal and maternal alleles for this variant that are annotated in the file.
//	vector<t_genotype*>* genotypes;
//};

class t_VCF
{
public:
	t_VCF(char* fp);
	~t_VCF();

	vector<t_vcf_entry*>* entries;

	// Following correspond to each other.
	vector<vector<t_vcf_entry*>*>* entries_per_chromosome;
	vector<char*>* chrom_ids;

	// Dump BED with extra information on the variants.
	void dump_VCF_BED(char* bed_fp);

	void dump_paternal_alternate_sites_in_ref_coords(char* bed_fp);
	void dump_maternal_alternate_sites_in_ref_coords(char* bed_fp);
	void dump_het_variants_in_ref_coords(char* bed_fp);

	// Parse the info string for this entry and return the information.
	char* get_info_per_entry(t_vcf_entry* vcf_entry, char* info_id);
};

#endif // __VCF__