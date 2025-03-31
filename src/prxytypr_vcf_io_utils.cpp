#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <ctype.h>
#include "prxytypr_utils.h"
#include "prxytypr_vector_macros.h"
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
#include "prxytypr_proxytyper.h"
#include "prxytypr_ansi_mutex.h"
#include <queue>
#include <zlib.h>

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

const bool __DUMP_VCF_IO_UTILS_MSGS__ = false;

void strcpy_protected(char* dest, size_t l_dest, char* src)
{
	if (strlen(src) >= l_dest)
	{
		fprintf(stderr, "Cannot copy %s\n", src);
		exit(1);
	}

	strcpy(dest, src);
}

static void* thread_callback_extract_genotype_signals_per_VCF_no_buffer_multithreaded(void* __thread_info_ptr)
{
	void** cur_thread_info = (void**)(__thread_info_ptr);
	int* int_vals = (int*)(cur_thread_info[0]);
	int i_thread = int_vals[0];
	//int n_threads = int_vals[1];
	int_vals[2] = 0; 
	int haplotype_specific_encoding = int_vals[3];
	int add_AF_info_2_id = int_vals[4];

	queue<char*>* cur_thread_q = (queue<char*>*)(cur_thread_info[2]);

	t_ansi_mutex* cur_thread_mutex = (t_ansi_mutex*)(cur_thread_info[3]);

	vector<char*>* vcf_sample_ids = (vector<char*>*)(cur_thread_info[4]);
	
	// Allocate the genotype vector.
	char* cur_var_geno_sig = new char[(int)vcf_sample_ids->size() + 2];

	int LINE_BUFFER_L_BUFF = 1000 * 1000;
	char* buff = new char[LINE_BUFFER_L_BUFF];

	char variants_BED_fp[1000];
	sprintf(variants_BED_fp, "vcf_variants_%d.bed", i_thread);
	char geno_matrix_fp[1000];
	sprintf(geno_matrix_fp, "vcf_genotypes_%d.matrix.gz", i_thread);
	FILE* f_vars_BED = open_f(variants_BED_fp, "w");
	//FILE* f_geno_matrix = open_f(geno_matrix_fp, "wb");

	gzFile f_geno_matrix = gzopen(geno_matrix_fp, "wb");
	if (!f_geno_matrix)
	{
		fprintf(stderr, "Failed to open VCF file: %s (Make sure it is gzipped?)\n", geno_matrix_fp);
		exit(1);
	}

	// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
	size_t INFO_FIELD_L_BUFF = 1000 * 1000;
	size_t VAR_COL_L_BUFF = 10000;
	char* chrom = new char[VAR_COL_L_BUFF];
	char* posn_str = new char[VAR_COL_L_BUFF];
	char* id = new char[VAR_COL_L_BUFF];
	char* ref = new char[VAR_COL_L_BUFF];
	char* alt = new char[VAR_COL_L_BUFF];
	char* qual = new char[VAR_COL_L_BUFF];
	char* filter = new char[VAR_COL_L_BUFF];
	char* info = new char[INFO_FIELD_L_BUFF];
	char* format = new char[VAR_COL_L_BUFF];

	int n_err_msgs = 0;
	while (1)
	{
		// Read the next line.
		char* cur_line = NULL;
		bool read_null = false;

		// Wait until we can lock the mutex.
		while (!cur_thread_mutex->try_lock_mutex()) {}

		// If the q is not empty, load the line and pop it.
		bool q_empty = false;
		if (!cur_thread_q->empty())
		{
			//if (i_thread == 0 &&
			//	int_vals[2] % 1000 == 0)
			//{
			//	fprintf(stderr, "@ %d. variant [%d variants in buffer]\r", int_vals[2], vecsize(cur_thread_q));
			//}
			cur_line = cur_thread_q->front();
			cur_thread_q->pop();

			// If we popped a null value, this signals end of file.
			if (cur_line == NULL)
			{
				read_null = true;
			}
		}
		else
		{
			q_empty = true;
		}

		// Release the mutex
		cur_thread_mutex->release_mutex();

		// If the q was empty, read again.
		if (q_empty == true)
		{
			continue;
		}

		// If the popped line is null, break.
		if (read_null)
		{
			break;
		}

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(chrom, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(posn_str, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(id, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(ref, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(alt, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(qual, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(filter, VAR_COL_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(info, INFO_FIELD_L_BUFF, buff);
		t_string::get_next_token(cur_line, buff, LINE_BUFFER_L_BUFF, "\t", i_cur_char);
		strcpy_protected(format, VAR_COL_L_BUFF, buff);

		//////////////////////////////////////////////////////////////////////////////////////////
		// Do the filtering here: Take SNPs and valid mutations.
		if (t_string::string_length(ref) != 1 ||
			t_string::string_length(alt) != 1)
		{
			delete[] cur_line;
			continue;
		}

		if (toupper(ref[0]) != 'A' && 
			toupper(ref[0]) != 'C' && 
			toupper(ref[0]) != 'G' && 
			toupper(ref[0]) != 'T')
		{
			delete[] cur_line;
			continue;
		}

		if (toupper(alt[0]) != 'A' && 
			toupper(alt[0]) != 'C' &&
			toupper(alt[0]) != 'G' &&
			toupper(alt[0]) != 'T')
		{
			delete[] cur_line;
			continue;
		}
		//////////////////////////////////////////////////////////////////////////////////////////

		int cur_var_posn = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
		int cur_var_end = cur_var_posn + t_string::string_length(ref) - 1;

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

		if (__DUMP_VCF_IO_UTILS_MSGS__)
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
		memset(cur_var_geno_sig, 0xff, sizeof(char) * vecsize_t(vcf_sample_ids));

		bool correctly_parsed_all_genotypes = true;
		double AAC_cnt = 0;
		double n_total_valid_haps = 0;
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
				if (n_err_msgs < 100)
				{
					fprintf(stderr, "Failed to read the genotype correctly for entry %d in line: %s:%s: %s. sample: %s (Potentially multiallelic)\n", i_s, chrom, posn_str, id, buff);
					n_err_msgs++;
				}
				//exit(1);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[1] == '/' && buff[0] != '.' && buff[2] != '.')
				{
					if (n_err_msgs < 100)
					{
						fprintf(stderr, "***WARNING::There are unphased genotypes although phased option is selected..***\n");
						n_err_msgs++;
					}
				}

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

					AAC_cnt += (geno0_val + geno1_val);
					n_total_valid_haps+=2;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;

					AAC_cnt += (geno0_val + geno1_val);
					n_total_valid_haps += 2;
				}
			} // phased panel block.
			else
			{
				// Unphased panel processing.
				if (buff[1] == '|' && buff[0] != '.' && buff[2] != '.')
				{
					if (n_err_msgs < 100)
					{
						fprintf(stderr, "***WARNING::There are phased genotypes although unphased option is selected..***\n");
						n_err_msgs++;
					}
				}

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

				if (cur_var_geno_sig[i_s] >= 0)
				{
					AAC_cnt += cur_var_geno_sig[i_s];
					n_total_valid_haps += 2;

				}
			} // unphased panel block.

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
			// If the variant identifier is not found, replace it with a generic one.
			if (id[0] == '.')
			{
				sprintf(id, "Thr%dVar%d", i_thread, int_vals[2]);
			}

			// Update the id with frequency.
			double AAF = 0;
			if (n_total_valid_haps > 0)
			{
				AAF = AAC_cnt / n_total_valid_haps;
			}

			// Use full id?
			char full_id[1000];
			if (add_AF_info_2_id == 1)
			{
				t_string::replace_avoid_list(id, ":*-_$%@!~ ()[]|#%^&?", '.');
				sprintf(full_id, "%s_%s_%s_%.5f", id, ref, alt, AAF);
			}
			else
			{
				strcpy(full_id, id);
			}
			
			// Write the bed file.
			fprintf(f_vars_BED, "%s\t%d\t%d\t%s\t.\t+\n", chrom,
				translate_coord(cur_var_posn, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(cur_var_end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				full_id);

			// Write the genome matrix.
			//fwrite(cur_var_geno_sig, sizeof(char), vecsize_t(vcf_sample_ids), f_geno_matrix);
			//int bytes_written = gzwrite(f_cur_thread_proxized_geno_mat, cur_subj_var_i_proxized_geno_sig, vecsize(sample_ids));
			int bytes_written = gzwrite(f_geno_matrix, cur_var_geno_sig, vecsize_t(vcf_sample_ids) * sizeof(char));
			if (bytes_written != vecsize(vcf_sample_ids))
			{
				fprintf(stderr, "Could not write genotype vector..\n");
				exit(1);
			}

			// Update the number of lines processed by the thread.
			int_vals[2]++;
		} // genotype parsing check.

		delete[] cur_line;
	} // buffer reading loop.

	fprintf(stderr, "Thread-%d completed..          \r", i_thread);
	close_f(f_vars_BED, NULL);
	//close_f(f_geno_matrix, geno_matrix_fp);
	gzclose(f_geno_matrix);

	return(NULL);
} // thread_callback_extract_genotype_signals_per_VCF_no_buffer_multithreaded function.

// Extract VCF and save to matrix.
void extract_genotype_signals_per_VCF_no_buffer_multithreaded(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
	bool add_AF_info_2_id,
	int n_threads,
	char* op_prefix)						// Output file path.
{
	auto vcf_import_start_chrono = std::chrono::high_resolution_clock::now();
	fprintf(stderr, "%d-thread parsing %s:\n\
Subjects: %s\n\
haplotype_specific_coding: %d\n\
add_AF_info_2_id: %d\n\
op_prefix: %s\n", n_threads, vcf_fp, vcf_sample_ids_list_fp, (int)haplotype_specific_encoding, (int)add_AF_info_2_id, op_prefix);

	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(1);
	}

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", (int)vcf_sample_ids->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}

	vector<queue<char*>*>* per_thread_VCF_buffer = new vector<queue<char*>*>();
	vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
	vector<t_ansi_mutex*>* per_thread_mutexes = new vector<t_ansi_mutex*>();
	vector<int*>* per_thread_int_vals = new vector<int*>();
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		void** cur_thread_info = new void* [10];
		int* int_vals = new int[10];
		int_vals[0] = i_thread;
		int_vals[1] = n_threads;
		int_vals[2] = 0; // # of VCF lines processed.
		int_vals[3] = (int)haplotype_specific_encoding;
		int_vals[4] = (int)add_AF_info_2_id;
		cur_thread_info[0] = int_vals;
		per_thread_int_vals->push_back(int_vals);

		double* dbl_vals = new double[10];
		cur_thread_info[1] = dbl_vals;

		queue<char*>* cur_thread_q = new queue<char*>();
		cur_thread_info[2] = cur_thread_q;
		per_thread_VCF_buffer->push_back(cur_thread_q);

		t_ansi_mutex* cur_thread_mutex = new t_ansi_mutex();
		cur_thread_info[3] = cur_thread_mutex;
		per_thread_mutexes->push_back(cur_thread_mutex);

		cur_thread_info[4] = vcf_sample_ids;

		fprintf(stderr, "Starting thread-%d..         \r", i_thread);
		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_extract_genotype_signals_per_VCF_no_buffer_multithreaded, cur_thread_info);
		cur_thread->run_thread();
		threads->push_back(cur_thread);
	} // i_thread loop.

	// Start loading VCF.
	//FILE* f_vcf = open_f(vcf_fp, "r");
	gzFile vcf_gzFile = gzopen(vcf_fp, "rb");
	if (!vcf_gzFile)
	{
		fprintf(stderr, "Failed to open VCF file: %s (Make sure it is gzipped?)\n", vcf_fp);
		exit(1);
	}

	//char* buff = new char[100 * 1000];
	int n_processed_regions = 0;
	int n_assigned_regions = 0;

	size_t l_line_buffer = 10 * 1000 * 1000;
	char* line_buffer = new char[l_line_buffer];

	while (gzgets(vcf_gzFile, line_buffer, l_line_buffer))
	{
		int l_str = strlen(line_buffer);
		char* cur_line = new char[l_str + 1];
		memcpy(cur_line, line_buffer, l_str);
		cur_line[l_str - 1] = 0;

		//char* cur_line = getline(f_vcf);
		//if (cur_line == NULL)
		//{
		//	break;
		//}

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

		// Add the line to the corresponding thread.
		int thread_i_2_process = n_processed_regions % n_threads;

		// Wait to lock the mutex.
		while (!per_thread_mutexes->at(thread_i_2_process)->try_lock_mutex()) {}

		// Mutex is locked, add to the list, then release.
		per_thread_VCF_buffer->at(thread_i_2_process)->push(cur_line);

		per_thread_mutexes->at(thread_i_2_process)->release_mutex();
	} // vcf file reading loop.

	gzclose(vcf_gzFile);

	// Finish the queues.
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		while (!per_thread_mutexes->at(thread_i)->try_lock_mutex()) {}

		// Mutex is locked, add to the list, then release.
		per_thread_VCF_buffer->at(thread_i)->push(NULL);

		per_thread_mutexes->at(thread_i)->release_mutex();
	} // thread_i loop.

	vector<char*>* bed_files = new vector<char*>();
	vector<char*>* geno_matrix_files = new vector<char*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		threads->at(thread_i)->wait_thread();
		fprintf(stderr, "Thread-%d completed [%d variants..]        \r", thread_i, per_thread_int_vals->at(thread_i)[2]);

		char variants_BED_fp[1000];
		sprintf(variants_BED_fp, "vcf_variants_%d.bed", thread_i);
		char geno_matrix_fp[1000];
		sprintf(geno_matrix_fp, "vcf_genotypes_%d.matrix.gz", thread_i);

		if (!check_file(geno_matrix_fp) ||
			!check_file(variants_BED_fp))
		{
			fprintf(stderr, "Could not find genotype matrix files: %s, %s\n", variants_BED_fp, geno_matrix_fp);
			exit(1);
		}

		bed_files->push_back(t_string::copy_me_str(variants_BED_fp));
		geno_matrix_files->push_back(t_string::copy_me_str(geno_matrix_fp));
	} // thread_i loop.

	fprintf(stderr, "Pooling %d files..                            \n", vecsize(bed_files));
	// Pool.
	char geno_matrix_fp[1000];
	sprintf(geno_matrix_fp, "%s_genotypes.matrix.gz", op_prefix);
	concatenateGzipFiles(geno_matrix_fp, geno_matrix_files);

	// Pool the per thread matrices.
	char variants_BED_fp[1000];
	sprintf(variants_BED_fp, "%s_variants.bed", op_prefix);
	concatenateGzipFiles(variants_BED_fp, bed_files);

	char subject_ids_list_fp[1000];
	sprintf(subject_ids_list_fp, "%s_subjects.list", op_prefix);
	FILE* f_subjects = open_f(subject_ids_list_fp, "w");
	for (int i_s = 0; i_s < vecsize(vcf_sample_ids); i_s++)
	{
		fprintf(f_subjects, "%s\n", vcf_sample_ids->at(i_s));
	} // i_s loop.
	close_f(f_subjects, NULL);

	fprintf(stderr, "Cleaning up..\n");
	for(int i_thread = 0; i_thread < vecsize(bed_files); i_thread++)
	{
		delete_file(bed_files->at(i_thread));
		delete_file(geno_matrix_files->at(i_thread));
	} // thread_i loop.

	auto vcf_import_end_chrono = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> vcf_import_duration = vcf_import_end_chrono - vcf_import_start_chrono;
	fprintf(stderr, "VCF import finished in %.3f seconds..\n", vcf_import_duration.count());
} // extract_genotype_signals_per_VCF_no_buffer_multithreaded


