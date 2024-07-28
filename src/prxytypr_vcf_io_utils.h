#ifndef __VCF_IO_UTILS__
#define __VCF_IO_UTILS__


void extract_genotype_signals_per_VCF_no_buffer_multithreaded(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	bool add_AF_info_2_id,
	int n_threads,
	char* op_prefix);


#endif // __VCF_IO_UTILS__