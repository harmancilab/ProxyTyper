#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "prxytypr_utils.h"
#include "prxytypr_ansi_string.h"
#include "prxytypr_rng.h"
#include "prxytypr_vector_macros.h"
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

#include <vector>
#include <algorithm>
#include<functional>

using namespace std;

const bool __DUMP_KMER_TRACKING_MSGS__ = false;

//struct t_kmer
//{
//	int kmer_length;
//	char* kmer;
//};

void get_allele_error_per_decoded_panel_known_panel(char* decoded_panel_matbed_fp, char* known_panel_matbed_fp, char* sample_ids_list_fp, char* op_fp)
{
	vector<t_annot_region*>* decoded_geno_regs = load_variant_signal_regions_wrapper(decoded_panel_matbed_fp, sample_ids_list_fp);
	sort(decoded_geno_regs->begin(), decoded_geno_regs->end(), sort_regions);
	vector<t_annot_region*>* known_geno_regs = load_variant_signal_regions_wrapper(known_panel_matbed_fp, sample_ids_list_fp);
	sort(known_geno_regs->begin(), known_geno_regs->end(), sort_regions);

	if (decoded_geno_regs->size() != known_geno_regs->size())
	{
		fprintf(stderr, "Genotype regions are not the same size: %d/%d\n", 
			vecsize(decoded_geno_regs), vecsize(known_geno_regs));

		exit(1);
	}

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_s = 0; i_s < vecsize(sample_ids); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			int n_errors = 0;

			for (int i_var = 0; i_var < vecsize(decoded_geno_regs); i_var++)
			{
				if (decoded_geno_regs->at(i_var)->start != known_geno_regs->at(i_var)->start)
				{
					fprintf(stderr, "Variant coordiantes are not matching between known and decoded variants.\n");
					exit(1);
				}

				void** decoded_var_info = (void**)(decoded_geno_regs->at(i_var)->data);
				char* decoded_geno = (char* )(decoded_var_info[0]);
				void** known_var_info = (void**)(known_geno_regs->at(i_var)->data);
				char* known_geno = (char*)(known_var_info[0]);

				char decoded_allele = get_allele_per_haplotype(decoded_geno[i_s], i_hap);
				char known_alelle = get_allele_per_haplotype(known_geno[i_s], i_hap);

				if (decoded_allele != known_alelle)
				{
					n_errors++;
				}
			} // 

			fprintf(f_op, "%s\t%d\t%d\t%d\n", sample_ids->at(i_s), i_hap, n_errors, vecsize(decoded_geno_regs));
		} // i_hap loop.
	} // i_s loop.

	close_f(f_op, NULL);
}

void delete_kmer(t_kmer* kmer)
{
	delete[] kmer->kmer;
	delete kmer;
}

void delete_kmers(vector<t_kmer*>* kmers)
{
	for (size_t i_k = 0; i_k < kmers->size(); i_k++)
	{
		delete_kmer(kmers->at(i_k));
	} // i_k loop.

	delete kmers;
}

void dump_kmer(t_kmer* kmer)
{
	fprintf(stderr, "kmer is: ");
	for (int i = 0; i < kmer->kmer_length; i++)
	{
		fprintf(stderr, "%d", kmer->kmer[i]);
	} // i loop.

	fprintf(stderr, "\n");
}

void dump_kmers(vector<t_kmer*>* kmers)
{
	for (size_t i_k = 0; i_k < kmers->size(); i_k++)
	{
		fprintf(stderr, "k-mer %d: ", (int)i_k);
		for (int i = 0; i < kmers->at(i_k)->kmer_length; i++)
		{
			fprintf(stderr, "%d", kmers->at(i_k)->kmer[i]);
		} // i loop.

		fprintf(stderr, "\n");
	} // i_k loop.
}

bool compare_kmers(t_kmer* kmer1, t_kmer* kmer2)
{
	for (int i = 0; i < kmer1->kmer_length; i++)
	{
		if (kmer1->kmer[i] != kmer2->kmer[i])
		{
			return(false);
		}
	}
	return(true);
}

bool sort_kmers(t_kmer* kmer1, t_kmer* kmer2)
{
	for (int i = 0; i < kmer1->kmer_length; i++)
	{
		if (kmer1->kmer[i] != kmer2->kmer[i])
		{
			return(kmer1->kmer[i] < kmer2->kmer[i]);
		}
	} // i loop.

	return(false);
}

t_kmer* copy_kmer(char* kmer_array, int kmer_length)
{
	t_kmer* kmer_copy = new t_kmer();
	kmer_copy->kmer_length = kmer_length;
	kmer_copy->kmer = new char[kmer_length + 2];

	for (int i = 0; i < kmer_length; i++)
	{
		kmer_copy->kmer[i] = kmer_array[i];
	} // i loop.

	return(kmer_copy);
} // copy_kmer function.

t_kmer* copy_kmer(t_kmer* kmer)
{
	t_kmer* kmer_copy = new t_kmer();
	kmer_copy->kmer = new char[kmer->kmer_length + 2];
	kmer_copy->kmer_length = kmer->kmer_length;

	for (int i = 0; i < kmer->kmer_length; i++)
	{
		kmer_copy->kmer[i] = kmer->kmer[i];
	} // i loop.

	return(kmer_copy);
} // copy_kmer function.

vector<t_kmer*>* copy_kmers(vector<t_kmer*>* kmers_list)
{
	vector<t_kmer*>* copied_kmers = new vector<t_kmer*>();
	for (size_t i = 0; i < kmers_list->size(); i++)
	{
		t_kmer* kmer_copy = copy_kmer(kmers_list->at(i));
		copied_kmers->push_back(kmer_copy);
	} // i loop.

	return(copied_kmers);
}

bool shift_kmers(vector<t_kmer*>* kmers_list, int n_nucs_2_shift_left)
{
	for (int i = 0; i < (int)kmers_list->size(); i++)
	{
		char* new_kmer = new char[kmers_list->at(i)->kmer_length + 1];
		memset(new_kmer, 0, sizeof(char) * (kmers_list->at(i)->kmer_length + 1));

		for (int nuc_i = 0; nuc_i < kmers_list->at(i)->kmer_length - n_nucs_2_shift_left; nuc_i++)
		{
			new_kmer[nuc_i] = kmers_list->at(i)->kmer[nuc_i + n_nucs_2_shift_left];
		} // i loop.

		delete[] kmers_list->at(i)->kmer;
		kmers_list->at(i)->kmer = new_kmer;
	} // i loop.

	return true;
}

vector<t_kmer*>* get_unique_kmers(vector<t_kmer*>* kmers_list)
{
	vector<t_kmer*>* unique_kmers = new vector<t_kmer*>();

	if (kmers_list->size() == 0)
	{
		return(unique_kmers);
	}

	// Sort the kmers.
	sort(kmers_list->begin(), kmers_list->end(), sort_kmers);

	unique_kmers->push_back(copy_kmer(kmers_list->at(0)));
	for (size_t i_kmer = 1; i_kmer < kmers_list->size(); i_kmer++)
	{
		if (!compare_kmers(kmers_list->at(i_kmer), unique_kmers->back()))
		{
			unique_kmers->push_back(copy_kmer(kmers_list->at(i_kmer)));
		}
	} // i_kmer loop.

	return(unique_kmers);
}
//
//void get_unique_kmers_w_counts_unsorted_indices(vector<t_kmer*>* par_kmers_list, 
//	vector<t_kmer*>* unique_kmers, vector<int>* unique_kmer_counts,
//	vector<int>* unique_kmer_i_per_orig_kmer_i)
//{
//	//vector<t_kmer*>* unique_kmers = new vector<t_kmer*>();
//
//	// Copy the kmers.
//	vector<t_kmer*>* kmers_list = new vector<t_kmer*>();
//	kmers_list->insert(kmers_list->end(), par_kmers_list->begin(), par_kmers_list->end());
//
//	if (kmers_list->size() == 0)
//	{
//		return;
//	}
//
//	// Sort the kmers.
//	sort(kmers_list->begin(), kmers_list->end(), sort_kmers);
//
//	t_kmer* cur_unique_kmer = kmers_list->at(0);
//	int n_kmer_per_cur_unique_kmer = 1;
//	for (size_t i_kmer = 1; i_kmer < kmers_list->size(); i_kmer++)
//	{
//		if (!compare_kmers(kmers_list->at(i_kmer), cur_unique_kmer))
//		{
//			// Copy the current unique kmer.
//			unique_kmers->push_back(copy_kmer(cur_unique_kmer));
//			unique_kmer_counts->push_back(n_kmer_per_cur_unique_kmer);
//			unique_kmer_i_per_orig_kmer_i->push_back(unique_kmers->size() - 1);
//
//			cur_unique_kmer = kmers_list->at(i_kmer);
//			n_kmer_per_cur_unique_kmer = 1;
//		}
//		else
//		{
//			// Add the pointer to the unique kmer for this original kmer.
//			unique_kmer_i_per_orig_kmer_i->push_back(unique_kmers->size() - 1);
//			n_kmer_per_cur_unique_kmer++;
//		}
//	} // i_kmer loop.
//
//	// Add the last kmer.
//	if (n_kmer_per_cur_unique_kmer > 0)
//	{
//		unique_kmers->push_back(copy_kmer(cur_unique_kmer));
//		unique_kmer_counts->push_back(n_kmer_per_cur_unique_kmer);
//		unique_kmer_i_per_orig_kmer_i->push_back(unique_kmers->size() - 1);
//	}
//
//	delete kmers_list;
//}


void get_unique_kmers_w_counts(vector<t_kmer*>* par_kmers_list, vector<t_kmer*>* unique_kmers, vector<int>* unique_kmer_counts)
{
	//vector<t_kmer*>* unique_kmers = new vector<t_kmer*>();

	// Copy the kmers.
	vector<t_kmer*>* kmers_list = new vector<t_kmer*>();
	kmers_list->insert(kmers_list->end(), par_kmers_list->begin(), par_kmers_list->end());

	if (kmers_list->size() == 0)
	{
		return;
	}

	// Sort the kmers.
	sort(kmers_list->begin(), kmers_list->end(), sort_kmers);

	t_kmer* cur_unique_kmer = kmers_list->at(0);
	int n_kmer_per_cur_unique_kmer = 1;
	for (size_t i_kmer = 1; i_kmer < kmers_list->size(); i_kmer++)
	{
		if (!compare_kmers(kmers_list->at(i_kmer), cur_unique_kmer))
		{
			// Copy the current unique kmer.
			unique_kmers->push_back(copy_kmer(cur_unique_kmer));
			unique_kmer_counts->push_back(n_kmer_per_cur_unique_kmer);

			cur_unique_kmer = kmers_list->at(i_kmer);
			n_kmer_per_cur_unique_kmer = 1;
		}
		else
		{
			n_kmer_per_cur_unique_kmer++;
		}
	} // i_kmer loop.

	// Add the last kmer.
	if (n_kmer_per_cur_unique_kmer > 0)
	{
		unique_kmers->push_back(copy_kmer(cur_unique_kmer));
		unique_kmer_counts->push_back(n_kmer_per_cur_unique_kmer);
	}

	delete kmers_list;
}


t_kmer* allocate_kmer(int kmer_length)
{
	t_kmer* new_kmer = new t_kmer();
	new_kmer->kmer_length = kmer_length;
	new_kmer->kmer = new char[kmer_length + 1];

	// Reset to a 0 length string.
	memset(new_kmer->kmer, 0, sizeof(char) * (kmer_length + 1));

	return(new_kmer);
}

vector<t_kmer*>* expand_kmers_by_last_nuc(vector<t_kmer*>* kmers_list)
{
	vector<t_kmer*>* expanded_kmers = new vector<t_kmer*>();

	// First, copy the kmers.
	vector<t_kmer*>* kmers0 = copy_kmers(kmers_list);
	vector<t_kmer*>* kmers1 = copy_kmers(kmers_list);

	// Replace the last nucleotide accordingly.
	for (int i = 0; i < (int)kmers_list->size(); i++)
	{
		int kmer_length = kmers_list->at(i)->kmer_length;
		kmers0->at(i)->kmer[kmer_length - 1] = 0;
		kmers1->at(i)->kmer[kmer_length - 1] = 1;
	} // i loop.

	// Add the expanded kmers to the list; this doubles the size of kmers.
	expanded_kmers->insert(expanded_kmers->end(), kmers0->begin(), kmers0->end());
	expanded_kmers->insert(expanded_kmers->end(), kmers1->begin(), kmers1->end());

	return(expanded_kmers);
}

vector<t_kmer*>* extract_kmers_per_haplotype(vector<char*>* per_ind_haplotypes, int i_win_start, int i_win_end)
{
	vector<t_kmer*>* kmers = new vector<t_kmer*>();

	int kmer_length = i_win_end - i_win_start + 1;

	// For every negative index, assign 0 to kmer.
	for (int i_s = 0; i_s < (int)per_ind_haplotypes->size(); i_s++)
	{
		t_kmer* new_kmer = allocate_kmer(kmer_length);
		new_kmer->kmer_length = kmer_length;

		for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
		{
			if (nuc_i >= 0)
			{
				// Copy the kmer.
				new_kmer->kmer[nuc_i - i_win_start] = (per_ind_haplotypes->at(i_s))[nuc_i];
			}
		} // nuc_i loop.

		kmers->push_back(new_kmer);
	} // i loop.

	return(kmers);
}

vector<t_kmer*>* extract_kmers_per_haplotype(vector<char**>* per_ind_haplotypes, int i_win_start, int i_win_end)
{
	vector<t_kmer*>* kmers = new vector<t_kmer*>();

	int kmer_length = i_win_end - i_win_start + 1;

	// For every negative index, assign 0 to kmer.
	for (size_t i_s = 0; i_s < per_ind_haplotypes->size(); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			t_kmer* new_kmer = allocate_kmer(kmer_length);
			new_kmer->kmer_length = kmer_length;

			for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
			{
				if (nuc_i >= 0)
				{
					// Copy the kmer.
					new_kmer->kmer[nuc_i - i_win_start] = (per_ind_haplotypes->at(i_s)[i_hap])[nuc_i];
				}
			} // nuc_i loop.

			kmers->push_back(new_kmer);
		} // i_hap loop.
	} // i loop.

	return(kmers);
}

void decode_site_alleles_per_proxized_reference(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_proxized_haplocoded_geno_fp, char* ref_sample_ids_fp,
	int n_vicinity)
{
	fprintf(stderr, "Decoding proxized alleles of query site using reference panel.\n");

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_proxized_geno_var_regs = load_variant_signal_regions_wrapper(ref_proxized_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	fprintf(stderr, "Loaded %d (%d), %d/%d (%d) variants for query-proxy, ref-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_proxized_geno_var_regs->size(), (int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_ref_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_proxized_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	// We track a list of candidate kmers at each position.
	FILE* f_dec_stats = open_f("decoded_stats.txt", "w");
	FILE* f_per_sample_dec_stats = open_f("per_sample_decoding_stats.txt", "w");
	for (int i_s = 0; i_s < (int)query_sample_ids->size(); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			fprintf(stderr, "Decoding query sample %s[%d]\n", query_sample_ids->at(i_s), i_hap);

			int n_known_query_mismatching_posns = 0;
			int n_no_path_sites = 0;

			vector<t_kmer*>* current_tracked_candidate_kmers = new vector<t_kmer*>();
			for (int i_var = 0; i_var < (int)query_proxized_geno_var_regs->size(); i_var++)
			{
				int i_win_start = i_var - n_vicinity;
				int i_win_end = i_var + n_vicinity;

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				// Extract the current window kmers within reference panel.
				vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
				}
				//dump_kmers(cur_win_all_ref_kmers);

				vector<t_kmer*>* cur_win_unique_all_ref_kmers = get_unique_kmers(cur_win_all_ref_kmers);
				int n_unique_ref_kmers = (int)cur_win_unique_all_ref_kmers->size();

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d unique all kmers @ [%d, %d]\n", (int)cur_win_unique_all_ref_kmers->size(), i_win_start, i_win_end);
				}

				delete_kmers(cur_win_all_ref_kmers);
				delete_kmers(cur_win_unique_all_ref_kmers);

				// This allele is important to select the list of haplotypes with matching proxy-allele in the reference panel.
				char cur_query_proxy_allele = per_query_subj_proxy_haplotypes->at(i_s)[i_hap][i_var];

				// Find the ref subjects with matching proxy-allele to the query.
				// Following loop matches the indexing order for building the haplotypes, i.e., [i_s][i_hap][i_var].
				vector<char*>* cur_win_candidate_full_ref_haplotypes = new vector<char*>();
				double tot_ref_panel_proxy_AF = 0;
				for (size_t i_r_s = 0; i_r_s < ref_sample_ids->size(); i_r_s++)
				{
					for (int i_r_hap = 0; i_r_hap < 2; i_r_hap++)
					{
						char cur_ref_proxy_allele = per_ref_subj_proxy_haplotypes->at(i_r_s)[i_r_hap][i_var];

						// If the query proxy allele matches the reference proxy allele, add the current reference kmer to the candidate reference kmers list.
						if (cur_query_proxy_allele == cur_ref_proxy_allele)
						{
							// Add this to the list of reference haplotypes (not kmers, yet).
							cur_win_candidate_full_ref_haplotypes->push_back(per_ref_subj_orig_haplotypes->at(i_r_s)[i_r_hap]);
						}

						tot_ref_panel_proxy_AF += cur_ref_proxy_allele;
					} // i_r_hap loop.
				} // i_r_s loop.

				// For the haplotypes with matching proxy allele between ref and query, get the original kmers for the reference panel.
				vector<t_kmer*>* cur_win_candidate_ref_kmers = extract_kmers_per_haplotype(cur_win_candidate_full_ref_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					// This is for sanity checks.
					fprintf(stderr, "%s:%d (%s) => Ref. Panel AF: %.4f\n",
						query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
						query_proxized_geno_var_regs->at(i_var)->name,
						tot_ref_panel_proxy_AF / ((int)ref_sample_ids->size() * 2));
				}

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d candidate original ref. haplotypes matching to query proxy allele %d:\n", (int)cur_win_candidate_full_ref_haplotypes->size(), cur_query_proxy_allele);
				}

				//dump_kmers(cur_win_candidate_ref_kmers);
				delete cur_win_candidate_full_ref_haplotypes;

				// If there are no reference kmers, we will simply not update the tracked kmers, this should simply re-start the whole process in the next variant.
				if (cur_win_candidate_ref_kmers->size() == 0)
				{
					if (__DUMP_KMER_TRACKING_MSGS__)
					{
						fprintf(stderr, "Could not find any more candidate haplotypes; Re-initializing..\n");
					}

					n_no_path_sites++;
					//getc(stdin);
				}

				// Now we have the candidates that match to the current query's proxy allele at the position.
				vector<t_kmer*>* cur_win_unique_candidate_ref_kmers = get_unique_kmers(cur_win_candidate_ref_kmers);
				delete_kmers(cur_win_candidate_ref_kmers);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d unique ref. candidate kmers with proxy allele %d.\n", (int)cur_win_unique_candidate_ref_kmers->size(), cur_query_proxy_allele);
					//dump_kmers(cur_win_unique_candidate_ref_kmers);
				}

				// Expand the candidate kmers from previous position after shifting.
				if (current_tracked_candidate_kmers == NULL ||
					current_tracked_candidate_kmers->size() == 0)
				{
					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Copying first set of k-mers from the unique candidates.\n"); }

					current_tracked_candidate_kmers = copy_kmers(cur_win_unique_candidate_ref_kmers);
				}
				else
				{
					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Shifting and expanding %d tracked k-mers.\n", (int)current_tracked_candidate_kmers->size()); }

					shift_kmers(current_tracked_candidate_kmers, 1);

					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Expanding %d tracked k-mers.\n", (int)current_tracked_candidate_kmers->size()); }

					vector<t_kmer*>* expanded_tracked_kmers = expand_kmers_by_last_nuc(current_tracked_candidate_kmers);

					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Expanded to %d kmers..\n", (int)expanded_tracked_kmers->size()); dump_kmers(expanded_tracked_kmers); }

					// Compare and filter the expanded candidate kmers with the filtered compatible kmers.
					vector<t_kmer*>* filtered_expanded_tracked_kmers = new vector<t_kmer*>();
					for (int i_kmer1 = 0; i_kmer1 < (int)expanded_tracked_kmers->size(); i_kmer1++)
					{
						for (int i_kmer2 = 0; i_kmer2 < (int)cur_win_unique_candidate_ref_kmers->size(); i_kmer2++)
						{
							if (compare_kmers(expanded_tracked_kmers->at(i_kmer1), cur_win_unique_candidate_ref_kmers->at(i_kmer2)))
							{
								filtered_expanded_tracked_kmers->push_back(expanded_tracked_kmers->at(i_kmer1));
								break;
							}
						} // i_kmer2 loop.
					} // i_kmer1 loop.

					if (filtered_expanded_tracked_kmers->size() == 0)
					{
						n_no_path_sites++;
						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "Could not track any more k-mers.\n");
						}
						//exit(1);
					}

					// We have the new updated set of kmers; update the current list of kmers.
					//delete_kmers(current_tracked_candidate_kmers);

					// Update the current tracked candidate kmers.
					if (current_tracked_candidate_kmers == NULL)
					{
						fprintf(stderr, "Sanity check failed..\n");
						exit(1);
					}
					else
					{
						delete_kmers(current_tracked_candidate_kmers);
					}

					current_tracked_candidate_kmers = get_unique_kmers(filtered_expanded_tracked_kmers);
				}

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%s[%d]::%s:%d => %d kmers are being tracked.\n",
						query_sample_ids->at(i_s), i_hap,
						query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
						(int)current_tracked_candidate_kmers->size());

					dump_kmers(current_tracked_candidate_kmers);
				}

				// Copy the known k-mer at this position:
				char cur_subj_haplotype[100];
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					if (nuc_i >= 0)
					{
						cur_subj_haplotype[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(i_s)[i_hap][nuc_i];
					}
					else
					{
						cur_subj_haplotype[nuc_i - i_win_start] = 0;
					}
				} // nuc_i loop.

				t_kmer* cur_subj_kmer = copy_kmer(cur_subj_haplotype, i_win_end - i_win_start + 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Known original query haplotype (Private):\n");
					dump_kmer(cur_subj_kmer);
				}

				bool found_kmer = false;
				for (size_t i_kmer = 0; i_kmer < current_tracked_candidate_kmers->size(); i_kmer++)
				{
					if (compare_kmers(cur_subj_kmer, current_tracked_candidate_kmers->at(i_kmer)))
					{
						found_kmer = true;
						break;
					}
				} // i_kmer loop.

				if (!found_kmer)
				{
					if (__DUMP_KMER_TRACKING_MSGS__)
					{
						fprintf(stderr, "!Could not find kmer!\n");
					}
					//getc(stdin);
					n_known_query_mismatching_posns++;
				}

				fprintf(f_per_sample_dec_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n",
					query_sample_ids->at(i_s), i_hap,
					query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
					(int)current_tracked_candidate_kmers->size(), n_unique_ref_kmers,
					n_known_query_mismatching_posns);
			} // i_var loop.

			fprintf(f_dec_stats, "%d\t%d\t%d\t%d\n", i_s, i_hap, n_known_query_mismatching_posns, n_no_path_sites);
			fflush(f_dec_stats);
		} // i_hap loop.
	} // i_s loop.

	close_f(f_dec_stats, NULL);
	close_f(f_per_sample_dec_stats, NULL);
} // decode_site_alleles_per_proxized_reference function.




void decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
																					char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
																					int n_vicinity)
{
	fprintf(stderr, "Decoding proxized alleles of query site using reference panel.\n");

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	//vector<t_annot_region*>* ref_proxized_geno_var_regs = load_variant_signal_regions_wrapper(ref_proxized_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	// We track a list of candidate kmers at each position.
	FILE* f_dec_stats = open_f("decoded_stats.txt", "w");
	FILE* f_per_sample_dec_stats = open_f("per_sample_decoding_stats.txt", "w");

	// The basic idea behind this is:
	// 1. Find all of the unique reference k-mers and frequencies in the current window.
	// 2. Extract the target individual's proxy haplotype and its frequency in the current window.
	// 3. Identify the candidate reference haplotypes match to in frequency to the target proxy haplotype.
	// 4. Shift-and-expand the tracked haplotypes from previous window.
	// 5. Overlap the expanded haplotypes with the current window's candidate haplotypes.
	// 6. Set the overlapping haplotypes as the new tracked haplotypes.
	// 
	for (int i_s = 0; i_s < (int)query_sample_ids->size(); i_s++)
	{
		for (int i_hap = 0; i_hap < 2; i_hap++)
		{
			fprintf(stderr, "Decoding query sample %s[%d]\n", query_sample_ids->at(i_s), i_hap);

			int n_known_query_mismatching_posns = 0;
			int n_no_path_sites = 0;

			// These are the kmers that are being tracked for the query panel.
			vector<t_kmer*>* current_tracked_candidate_kmers = new vector<t_kmer*>();

			// We are assuming that the variants are correctly aligned between reference and query.
			for (int i_var = 0; i_var < (int)query_proxized_geno_var_regs->size(); i_var++)
			{
				// Set the window starts and ends.
				int i_win_start = MAX(0, i_var - n_vicinity);
				int i_win_end = MIN(i_var + n_vicinity, vecsize(query_proxized_geno_var_regs)-1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				// Extract the current window's kmers within reference panel.
				vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

				// Extract the current window's kmer for the query panel.
				vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
				}
				//dump_kmers(cur_win_all_ref_kmers);

				// Get the window reference k-mers
				vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_ref_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_cnts);

				// Get the query proxy k-mers.
				vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

				int n_unique_ref_kmers = (int)cur_win_unique_ref_kmers->size();

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d unique all kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
				}

				// Assign the frequency to this query kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(i_s)[i_hap];
				t_kmer* cur_query_subj_proxy_kmer = cur_win_all_query_proxy_kmers->at(2 * i_s + i_hap);

				for(int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				}

				double cur_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
						break;
					}
				} // i_kmer loop.
				
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0)
				{
					fprintf(stderr, "Check failed!\n");
					exit(1);
				}

				// Following adds all of the reference kmers that have same or smaller frequency.
				vector<t_kmer*>* cur_win_unique_candidate_ref_kmers = new vector<t_kmer*>();
				for (int j_kmer = 0; j_kmer < vecsize(cur_win_unique_ref_kmers); j_kmer++)
				{
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_cnts->at(j_kmer) / (2 * vecsize(ref_sample_ids));

					// This is the main check that adds the candidates at this position by matching the frequency of the 
					// query proxy haplotype to the frequency distribution of the reference haplotype list.
					// Any reference haplotype that satisfies following is compared with the shift-expanded kmers from previous window.
					//if (fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq) < 0.01)
					//if (cur_ref_kmer_freq >= cur_query_subj_proxy_kmer_freq)
					if (cur_ref_kmer_freq <= cur_query_subj_proxy_kmer_freq)
					{
						cur_win_unique_candidate_ref_kmers->push_back(cur_win_unique_ref_kmers->at(j_kmer));

						if (cur_ref_kmer_freq == 0)
						{
							fprintf(stderr, "Check failed!\n");
							exit(1);
						}

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "Adding ref candidate: Proxy-hap Freq: %.3f ;; Ref-hap %d/%d Freq: %.3f\n",
								cur_query_subj_proxy_kmer_freq,
								j_kmer, vecsize(cur_win_unique_ref_kmers),
								cur_ref_kmer_freq);
						}
					}
				} // j_kmer loop.

				//if (__DUMP_KMER_TRACKING_MSGS__)
				//{
				//	// This is for sanity checks.
				//	fprintf(stderr, "%s:%d (%s)\n",
				//		query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
				//		query_proxized_geno_var_regs->at(i_var)->name);
				//}

				/*if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d candidate original ref. haplotypes matching to query proxy allele %d:\n", (int)cur_win_candidate_full_ref_haplotypes->size(), cur_query_proxy_allele);
				}*/

				//dump_kmers(cur_win_candidate_ref_kmers);
				//delete cur_win_candidate_full_ref_haplotypes;

				// If there are no reference kmers, we will simply not update the tracked kmers, this should simply re-start the whole process in the next variant.
				if (cur_win_unique_candidate_ref_kmers->size() == 0)
				{
					if (__DUMP_KMER_TRACKING_MSGS__)
					{
						fprintf(stderr, "Could not find any more candidate haplotypes; Re-initializing..\n");
					}

					n_no_path_sites++;
					//getc(stdin);
				}

				//if (__DUMP_KMER_TRACKING_MSGS__)
				//{
				//	fprintf(stderr, "Extracted %d unique ref. candidate kmers with proxy allele %d.\n", (int)cur_win_unique_candidate_ref_kmers->size(), cur_query_proxy_allele);
				//	//dump_kmers(cur_win_unique_candidate_ref_kmers);
				//}

				// Expand the candidate kmers from previous position after shifting.
				if (current_tracked_candidate_kmers == NULL ||
					current_tracked_candidate_kmers->size() == 0)
				{
					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Copying first set of k-mers from the unique candidates.\n"); }

					current_tracked_candidate_kmers = copy_kmers(cur_win_unique_candidate_ref_kmers);
				}
				else
				{
					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Shifting and expanding %d tracked k-mers.\n", (int)current_tracked_candidate_kmers->size()); }

					shift_kmers(current_tracked_candidate_kmers, 1);

					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Expanding %d tracked k-mers.\n", (int)current_tracked_candidate_kmers->size()); }

					vector<t_kmer*>* expanded_tracked_kmers = expand_kmers_by_last_nuc(current_tracked_candidate_kmers);

					if (__DUMP_KMER_TRACKING_MSGS__) { fprintf(stderr, "Expanded to %d kmers..\n", (int)expanded_tracked_kmers->size()); dump_kmers(expanded_tracked_kmers); }

					// Compare and filter the expanded candidate kmers with the filtered compatible kmers.
					vector<t_kmer*>* filtered_expanded_tracked_kmers = new vector<t_kmer*>();
					for (int i_kmer1 = 0; i_kmer1 < (int)expanded_tracked_kmers->size(); i_kmer1++)
					{
						for (int i_kmer2 = 0; i_kmer2 < (int)cur_win_unique_candidate_ref_kmers->size(); i_kmer2++)
						{
							if (compare_kmers(expanded_tracked_kmers->at(i_kmer1), cur_win_unique_candidate_ref_kmers->at(i_kmer2)))
							{
								filtered_expanded_tracked_kmers->push_back(expanded_tracked_kmers->at(i_kmer1));
								break;
							}
						} // i_kmer2 loop.
					} // i_kmer1 loop.

					if (filtered_expanded_tracked_kmers->size() == 0)
					{
						n_no_path_sites++;
						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "Could not track any more k-mers.\n");
						}
						//exit(1);
					}

					// We have the new updated set of kmers; update the current list of kmers.
					//delete_kmers(current_tracked_candidate_kmers);

					// Update the current tracked candidate kmers.
					if (current_tracked_candidate_kmers == NULL)
					{
						fprintf(stderr, "Sanity check failed..\n");
						exit(1);
					}
					else
					{
						delete_kmers(current_tracked_candidate_kmers);
					}

					current_tracked_candidate_kmers = get_unique_kmers(filtered_expanded_tracked_kmers);
				}

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%s[%d]::%s:%d => %d kmers are being tracked.\n",
						query_sample_ids->at(i_s), i_hap,
						query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
						(int)current_tracked_candidate_kmers->size());

					dump_kmers(current_tracked_candidate_kmers);
				}

				// Copy the known k-mer at this position:
				char cur_subj_haplotype[100];
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					if (nuc_i >= 0)
					{
						cur_subj_haplotype[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(i_s)[i_hap][nuc_i];
					}
					else
					{
						cur_subj_haplotype[nuc_i - i_win_start] = 0;
					}
				} // nuc_i loop.

				t_kmer* cur_subj_kmer = copy_kmer(cur_subj_haplotype, i_win_end - i_win_start + 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Known original query haplotype (Private):\n");
					dump_kmer(cur_subj_kmer);
				}

				bool found_kmer = false;
				for (size_t i_kmer = 0; i_kmer < current_tracked_candidate_kmers->size(); i_kmer++)
				{
					if (compare_kmers(cur_subj_kmer, current_tracked_candidate_kmers->at(i_kmer)))
					{
						found_kmer = true;
						break;
					}
				} // i_kmer loop.

				if (!found_kmer)
				{
					if (__DUMP_KMER_TRACKING_MSGS__)
					{
						fprintf(stderr, "!Could not find kmer!\n");
					}
					//getc(stdin);
					n_known_query_mismatching_posns++;
				}

				fprintf(f_per_sample_dec_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n",
					query_sample_ids->at(i_s), i_hap,
					query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
					(int)current_tracked_candidate_kmers->size(), n_unique_ref_kmers,
					n_known_query_mismatching_posns);
			} // i_var loop.

			fprintf(f_dec_stats, "%d\t%d\t%d\t%d\n", i_s, i_hap, n_known_query_mismatching_posns, n_no_path_sites);
			fflush(f_dec_stats);
		} // i_hap loop.
	} // i_s loop.

	close_f(f_dec_stats, NULL);
	close_f(f_per_sample_dec_stats, NULL);
} // decode_site_alleles_per_proxized_reference function.

double get_shifted_kmer_distance(t_kmer* left_kmer, t_kmer* right_kmer)
{
	double dist_val = 0;
	for (int i = 1; i < left_kmer->kmer_length; i++)
	{
		if (left_kmer->kmer[i] != right_kmer->kmer[i - 1])
		{
			dist_val++;
		}
	} // i loop.

	return(dist_val);
}

bool are_kmers_shift_compatible(t_kmer* left_kmer, t_kmer* right_kmer)
{
	for (int i = 1; i < left_kmer->kmer_length; i++)
	{
		if (left_kmer->kmer[i] != right_kmer->kmer[i - 1])
		{
			return(false);
		}
	} // i loop.

	return(true);
}

void decode_site_alleles_per_proxized_reference_2_hapfreq_ref_histogram_matching_HMM(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	int var_start_i, int var_end_i,
	int n_vicinity,
	char* op_prefix)
{
	fprintf(stderr, "Decoding proxized alleles of query site using histogram matching HMM:\n\
Query Original: %s\n\
Query Proxy: %s\n\
Query Sample list: %s\n\
Ref. Original: %s\n\
Ref. Sample list: %s\n\
# vic.: %d\n\
Variant range: [%d-%d]\n", query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
		ref_original_haplocoded_geno_fp, ref_sample_ids_fp, n_vicinity,
		var_start_i, var_end_i);

	// If this is set, the emission probabilities for each proxy kmer will be normalized by the total among query kmers at each window.
	bool USE_NORMALIZED_EMIT_PROB = false;

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	sort(query_proxized_geno_var_regs->begin(), query_proxized_geno_var_regs->end(), sort_regions);
	sort(query_original_geno_var_regs->begin(), query_original_geno_var_regs->end(), sort_regions);
	sort(ref_orig_geno_var_regs->begin(), ref_orig_geno_var_regs->end(), sort_regions);

	int max_query_proxy_geno = get_max_genotype_value(query_proxized_geno_var_regs, query_sample_ids);
	int max_query_original_geno = get_max_genotype_value(query_original_geno_var_regs, query_sample_ids);
	int max_ref_original_geno = get_max_genotype_value(ref_orig_geno_var_regs, ref_sample_ids);

	if (max_query_proxy_geno != 3 ||
		max_query_original_geno != 3 ||
		max_ref_original_geno != 3)
	{
		fprintf(stderr, "One of the panels is not haplocoded..\n");
		exit(1);
	}

	if (query_proxized_geno_var_regs->size() != query_original_geno_var_regs->size() ||
		query_proxized_geno_var_regs->size() != ref_orig_geno_var_regs->size())
	{
		fprintf(stderr, "Variant counts are not the same: %d, %d, %d\n",
				vecsize(query_proxized_geno_var_regs), vecsize(query_original_geno_var_regs),
				vecsize(ref_orig_geno_var_regs));

		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	// Check & fix variant range.
	if (var_start_i < 0 ||
		var_start_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_start_i = 0;
	}

	if (var_end_i < 0 || 
		var_end_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_end_i = vecsize(query_proxized_geno_var_regs);
	}

	if (var_end_i < (var_start_i + 100))
	{
		fprintf(stderr, "Invalid range: %d-%d\n", var_start_i, var_end_i);
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	// We track a list of candidate kmers at each position.
	char dec_summary_stats_fp[1000];
	sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary.txt", op_prefix);
	FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	char per_sample_per_var_dec_stats_fp[1000];
	sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats.txt", op_prefix);
	FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Go over each proxy query subject.
	for (int proxy_query_i_s = 0; proxy_query_i_s < (int)query_sample_ids->size(); proxy_query_i_s++)
	{
		// Process each haplotype, we assume these are perfectly phased, which can be accessed in task specific results.
		for (int proxy_query_i_hap = 0; proxy_query_i_hap < 2; proxy_query_i_hap++)
		{
			fprintf(stderr, "Decoding query sample %s[%d]\n", query_sample_ids->at(proxy_query_i_s), proxy_query_i_hap);

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Following store the window (or variant specific information)
			//
			// This is the max log-probability score for each unique kmer in the window.
			double** per_win_per_ref_kmer_log_cumul_probs = new double*[query_proxized_geno_var_regs->size()];
			// # of unique ref. kmers in the window.
			int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];
			// List of unique ref. kmers per window.
			vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
			// The best scoring unique kmer index in the previous window for each current window.
			int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
			// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
			double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];
			// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
			vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];

			// Following loop goes over all of the windows and fills up the DP arrays including scores and the backtracking array.
			// We are assuming that the variants are correctly aligned between reference and query.
			for (int i_var = var_start_i; i_var < var_end_i; i_var++)
			{
				// Set the window starts and ends.
				int i_win_start = MAX(0, i_var - n_vicinity);
				int i_win_end = MIN(i_var + n_vicinity, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////////////////////////////////
				// Extract the current window's kmers within reference panel.
				vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

				// Extract the current window's kmer for the query panel.
				vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d reference kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
				}
				//dump_kmers(cur_win_all_ref_kmers);

				/////////////////////////////////////////////////////////////////////////////////
				// Get the reference haplotype frequencies.
				// Get the window reference k-mers
				vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_ref_kmer_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_kmer_cnts);

				delete_kmers(cur_win_all_ref_kmers);

				int n_unique_ref_haps = cur_win_unique_ref_kmer_cnts->size();

				// These keep state of ref kmer scores, backtracks, etc.
				per_win_per_unique_ref_kmer_cnts[i_var] = cur_win_unique_ref_kmer_cnts;
				per_win_n_ref_kmers[i_var] = n_unique_ref_haps;
				per_win_per_ref_kmer_log_cumul_probs[i_var] = new double[n_unique_ref_haps];
				per_win_unique_ref_kmers[i_var] = cur_win_unique_ref_kmers;
				per_win_backtracking_kmer_i[i_var] = new int[n_unique_ref_haps];

				// Initialize the score and backtracking haplotype index in previous window.
				for (int i_kmer = 0; i_kmer < n_unique_ref_haps; i_kmer++)
				{
					per_win_per_ref_kmer_log_cumul_probs[i_var][i_kmer] = xlog(0.0);
					per_win_backtracking_kmer_i[i_var][i_kmer] = -1;
				} // i_kmer loop.
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Get the query proxy k-mers.
				vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d unique reference kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
					fprintf(stderr, "%d unique proxy query kmers @ [%d, %d]\n", (int)cur_win_unique_query_proxy_kmers->size(), i_win_start, i_win_end);
				}

				// Get the current proxy haplotype of the query subject and extract the kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(proxy_query_i_s)[proxy_query_i_hap];
				t_kmer* cur_query_subj_proxy_kmer = copy_kmer(cur_win_all_query_proxy_kmers->at(2 * proxy_query_i_s + proxy_query_i_hap));
				delete_kmers(cur_win_all_query_proxy_kmers);

				// This is a sanity check for the extracted proxy kmer.
				for (int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				}				

				// Search for the query proxy kmer in unique kmers for the window.
				double cur_query_subj_proxy_kmer_cnt = 0.0;
				double total_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
						//break;
					}

					total_query_subj_proxy_kmer_cnt += cur_win_unique_query_proxy_cnts->at(i_kmer);
				} // i_kmer loop.
				delete_kmer(cur_query_subj_proxy_kmer);
				delete_kmers(cur_win_unique_query_proxy_kmers);
				delete(cur_win_unique_query_proxy_cnts);

				// This is the frequency of proxy haplotype at this window that we will match to from the reference panel.
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0 ||
					cur_query_subj_proxy_kmer_freq > 1.0)
				{
					fprintf(stderr, "Sanity Check failed: Found 0 frequency for the query proxy haplotype.\n");
					exit(1);
				}

				// Sanity check on the total count of the query kmer count.
				if (total_query_subj_proxy_kmer_cnt != (2 * vecsize(query_sample_ids)))
				{
					fprintf(stderr, "Sanity Check failed: Total query ref kmer counts does not add up to 1: %.3f/%d",
							total_query_subj_proxy_kmer_cnt, 2 * vecsize(query_sample_ids));

					exit(1);
				}

				// Save the frequency of this haplotype.
				per_win_query_proxy_kmer_freq[i_var] = cur_query_subj_proxy_kmer_freq;
				/////////////////////////////////////////////////////////////////////////////////

				// This is the normalizing factor for the emission probabilities; ensures 
				// that the total emission probabilities of query kmers from reference kmers are equal to 1.
				double query_emit_log_normalizing_factor = xlog(0.0);
				double total_ref_unique_kmer_prob = 0; // This is for sanity check.
				for (int cur_kmer = 0; cur_kmer < vecsize(cur_win_unique_ref_kmers); cur_kmer++)
				{
					// Get the kmer frequency for this kmer.
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer) / (2 * vecsize(ref_sample_ids));
					if (cur_ref_kmer_freq == 0)
					{
						fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
						exit(1);
					}

					total_ref_unique_kmer_prob += cur_ref_kmer_freq;

					double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
					double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

					query_emit_log_normalizing_factor = xlog_sum(query_emit_log_normalizing_factor, proxy_kmer_emit_log_prob_per_cur_kmer);
				} // cur_kmer loop.

				// Make sure that total kmer counts equal to 1.
				if (fabs(total_ref_unique_kmer_prob - 1.0) > 0.001)
				{
					fprintf(stderr, "Sanity check failed on total unique ref prob: %.4f\n", total_ref_unique_kmer_prob);
					exit(1);
				}

				/////////////////////////////////////////////////////////////////////////////////
				// Following resets the normalizing factor to 1. 
				// If this is uncommented, the normalizing probability check below should also be commented out.
				if (!USE_NORMALIZED_EMIT_PROB)
				{
					query_emit_log_normalizing_factor = xlog(1.0);
				}
				
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Calculate DP arrays.
				// This is the main loop that updates scores and the backtracking array.
				if (i_var > var_start_i)
				{
					// Go over all current haplotypes and update scores by emission/transition on previous window.
					double total_cur_proxy_kmer_emit_log_prob_by_ref_kmers = xlog(0.0);
					for (int cur_kmer = 0; cur_kmer < vecsize(cur_win_unique_ref_kmers); cur_kmer++)
					{
						// The initial value of the cumulative prob must be exactly 0.
						if (per_win_per_ref_kmer_log_cumul_probs[i_var][cur_kmer] != xlog(0.0))
						{
							fprintf(stderr, "Sanity check failed: Could not init probs.\n");
							exit(1);
						}

						// Get the kmer frequency for this kmer.
						double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer) / (2 * vecsize(ref_sample_ids));
						if (cur_ref_kmer_freq == 0 || cur_ref_kmer_freq > 1.0)
						{
							fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is not valid: %.4f\n", cur_ref_kmer_freq);
							exit(1);
						}

						// This scores the proxy kmer emission by the current ref kmer: exp(-abs(delta hap freq.)).
						// Note that this does not rely on the previous window.
						//double proxy_kmer_emit_prob_per_cur_kmer = xlog(cur_ref_kmer_freq) - 1 * fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						//double proxy_kmer_emit_prob_per_cur_kmer = -1 * fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

						// This is the normalized kmer emission prob score.
						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_div(proxy_kmer_emit_log_prob_per_cur_kmer, query_emit_log_normalizing_factor);

						// This is for sanity check.
						total_cur_proxy_kmer_emit_log_prob_by_ref_kmers = xlog_sum(total_cur_proxy_kmer_emit_log_prob_by_ref_kmers, norm_proxy_kmer_emit_log_prob_per_cur_kmer);

						// Go over all of the previous window's haplotypes.
						for (int prev_kmer = 0; prev_kmer < per_win_n_ref_kmers[i_var - 1]; prev_kmer++)
						{
							// Get transition prob: 1/0 based on kmer overlap; transitions can include genetic distances? We need to keep track of 
							// all haplotypes individually for this.
							if (are_kmers_shift_compatible(per_win_unique_ref_kmers[i_var - 1]->at(prev_kmer), cur_win_unique_ref_kmers->at(cur_kmer)))
							{
								// Update the probability with the current previous haplotype if it leads to higher emission probability.
								double cur_prev2cur_trans_emit_log_prob = xlog_mul(per_win_per_ref_kmer_log_cumul_probs[i_var - 1][prev_kmer], norm_proxy_kmer_emit_log_prob_per_cur_kmer);
								if (per_win_per_ref_kmer_log_cumul_probs[i_var][cur_kmer] < cur_prev2cur_trans_emit_log_prob)
								{
									per_win_per_ref_kmer_log_cumul_probs[i_var][cur_kmer] = cur_prev2cur_trans_emit_log_prob;
									per_win_backtracking_kmer_i[i_var][cur_kmer] = prev_kmer;
								}
							} // transition check.
						} // prev_kmer loop.

						//fprintf(stderr, "i_var: %d: %d. haplotype: Score: %.5f (%d)\n", );
					} // cur_kmer loop.

					if (USE_NORMALIZED_EMIT_PROB && fabs(exp(total_cur_proxy_kmer_emit_log_prob_by_ref_kmers) - 1.0) > 0.001)
					{
						fprintf(stderr, "Sanity check failed: Total normalized emission probability is not 1: %.4f\n", exp(total_cur_proxy_kmer_emit_log_prob_by_ref_kmers));
						exit(1);
					}
				}
				else
				{
					// Processing of first window.
					for (int cur_kmer = 0; cur_kmer < vecsize(cur_win_unique_ref_kmers); cur_kmer++)
					{
						double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer) / (2 * vecsize(ref_sample_ids));

						if (cur_ref_kmer_freq == 0)
						{
							fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
							exit(1);
						}

						// This scores the proxy kmer emission by the current ref kmer.
						//double proxy_kmer_emit_prob_per_cur_kmer = xlog(cur_ref_kmer_freq) - 1 * fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						//double proxy_kmer_emit_prob_per_cur_kmer = - 1 * fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;
						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_div(proxy_kmer_emit_log_prob_per_cur_kmer, query_emit_log_normalizing_factor);

						per_win_per_ref_kmer_log_cumul_probs[i_var][cur_kmer] = norm_proxy_kmer_emit_log_prob_per_cur_kmer;
					} // cur_kmer loop.
				} // i_var>0 check.
				/////////////////////////////////////////////////////////////////////////////////
			} // i_var loop.

			/////////////////////////////////////////////////////////////////////////////////
			// Backtrack the scores; Find the top scoring haplotype at the last variant.
			int top_scoring_kmer_i = -1;
			double top_scoring_kmer_score = xlog(0.0);
			for (int i_kmer = 0; i_kmer < per_win_n_ref_kmers[var_end_i - 1]; i_kmer++)
			{
				if (top_scoring_kmer_score < per_win_per_ref_kmer_log_cumul_probs[var_end_i - 1][i_kmer])
				{
					top_scoring_kmer_i = i_kmer;
					top_scoring_kmer_score = per_win_per_ref_kmer_log_cumul_probs[var_end_i - 1][i_kmer];
				}
			} // i_kmer loop.

			double total_max_proxy_score_per_full_hap = top_scoring_kmer_score;

			// Start backtracking.
			int n_errs = 0;
			for (int i_var = var_end_i - 1; i_var >= var_start_i; i_var--)
			{
				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Backtracking @ var_i: %d ref kmer %d (%.5f)\n", i_var, top_scoring_kmer_i, per_win_per_ref_kmer_log_cumul_probs[i_var][top_scoring_kmer_i]);
				}

				int i_win_start = MAX(0, i_var - n_vicinity);
				int i_win_end = MIN(i_var + n_vicinity, vecsize(query_proxized_geno_var_regs) - 1);

				// Copy the original query subject's haplotype for the current query haplotype.
				char cur_win_query_subj_orig_kmer_nucs[100];
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_query_subj_orig_kmer_nucs[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(proxy_query_i_s)[proxy_query_i_hap][nuc_i];
				} // nuc_i loop.

				// Generate a kmer from the allele sequence.
				t_kmer* cur_win_query_subj_orig_kmer = copy_kmer(cur_win_query_subj_orig_kmer_nucs, i_win_end - i_win_start + 1);

				//////////////
				// This is for sanity check.
				vector<t_kmer*>* cur_win_all_query_orig_kmers = extract_kmers_per_haplotype(per_query_subj_orig_haplotypes, i_win_start, i_win_end);
				if (!compare_kmers(cur_win_query_subj_orig_kmer, cur_win_all_query_orig_kmers->at(2 * proxy_query_i_s + proxy_query_i_hap)))
				{
					fprintf(stderr, "Sanity check failed: Original k-mer mismatch @ backtrack.\n");
					exit(1);
				}
				delete_kmers(cur_win_all_query_orig_kmers);
				//////////////

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Var %d: Known original query haplotype (Private) [Proxy. freq: %.4f]:\n", i_var, per_win_query_proxy_kmer_freq[i_var]);
					dump_kmer(cur_win_query_subj_orig_kmer);

					double backtracked_ref_hap_AF = (double)(per_win_per_unique_ref_kmer_cnts[i_var]->at(top_scoring_kmer_i)) / (2 * vecsize(ref_sample_ids));
					fprintf(stderr, "Var %d: Backtracked reference kmer (%.4f):\n", i_var, backtracked_ref_hap_AF);
					dump_kmer(per_win_unique_ref_kmers[i_var]->at(top_scoring_kmer_i));
				}

				// Compare the max scoring ref kmer and the original query kmer.
				bool found_kmer = false;
				if (compare_kmers(cur_win_query_subj_orig_kmer, per_win_unique_ref_kmers[i_var]->at(top_scoring_kmer_i)))
				{
					found_kmer = true;
				}
				else
				{
					n_errs++;
				}
				delete_kmer(cur_win_query_subj_orig_kmer);				

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					if (!found_kmer)
					{
						fprintf(stderr, "ERROR @ var_i=%d\n", i_var);
					}
				}

				double backtracked_ref_hap_AF = (double)(per_win_per_unique_ref_kmer_cnts[i_var]->at(top_scoring_kmer_i)) / (2 * vecsize(ref_sample_ids));

				fprintf(f_per_sample_per_var_dec_stats, "%s\t%d\t%s\t%d\t%d\t%.4f\t%.4f\t%d\n",
					query_sample_ids->at(proxy_query_i_s), proxy_query_i_hap,
					query_proxized_geno_var_regs->at(i_var)->chrom, query_proxized_geno_var_regs->at(i_var)->start,
					per_win_n_ref_kmers[i_var], 
					per_win_query_proxy_kmer_freq[i_var], backtracked_ref_hap_AF, (int)found_kmer);

				// Update the state.
				if (i_var > var_start_i)
				{
					// Note that we update top_scoring_hap_i here, which refers to the previous variant's kmers now.
					top_scoring_kmer_i = per_win_backtracking_kmer_i[i_var][top_scoring_kmer_i];					

					if (top_scoring_kmer_i == -1 ||
						top_scoring_kmer_i >= per_win_n_ref_kmers[i_var-1])
					{
						fprintf(stderr, "Sanity check failed: Cannot backtrack @ var_i=%d; %d/%d\n", i_var, top_scoring_kmer_i, per_win_n_ref_kmers[i_var - 1]);
						exit(1);
					}
				}
				else
				{
					break;
				}
			} // i_var loop.

			fprintf(stderr, "%s\t%d\t%d\t%d\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_i_hap, n_errs, total_max_proxy_score_per_full_hap);
			fprintf(f_dec_stats, "%d\t%d\t%d\t%.5f\n", proxy_query_i_s, proxy_query_i_hap, n_errs, total_max_proxy_score_per_full_hap);
			fflush(f_dec_stats);

			// Free memory, this is necessary for long runs.
			//for(int i_var = 0; i_var < vecsize(query_proxized_geno_var_regs); i_var++)
			for (int i_var = var_start_i; i_var < var_end_i; i_var++)
			{
				// This is the max log-probability score for each unique kmer in the window.
				//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
				delete[] per_win_per_ref_kmer_log_cumul_probs[i_var];

				// # of unique ref. kmers in the window.
				//int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];

				// List of unique ref. kmers per window.
				//vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
				delete_kmers(per_win_unique_ref_kmers[i_var]);

				// The best scoring unique kmer index in the previous window for each current window.
				//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
				delete[] per_win_backtracking_kmer_i[i_var];

				// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
				//double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];

				// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
				//vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
				delete(per_win_per_unique_ref_kmer_cnts[i_var]);
			} // i_var loop.

			delete[] per_win_unique_ref_kmers;
			delete[] per_win_per_unique_ref_kmer_cnts;
			delete[] per_win_query_proxy_kmer_freq;
			delete[] per_win_backtracking_kmer_i;
			delete[] per_win_n_ref_kmers;
			delete[] per_win_per_ref_kmer_log_cumul_probs;
		} // i_hap loop.
	} // i_s loop.

	close_f(f_dec_stats, NULL);
	close_f(f_per_sample_per_var_dec_stats, NULL);
} // decode_site_alleles_per_proxized_reference function.

double get_self(double a)
{
	return(a);
}

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_Buffered(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	char* op_prefix)
{
	fprintf(stderr, "Decoding proxized alleles of query site using histogram matching HMM with Viterbi:\n\
Query Original: %s\n\
Query Proxy: %s\n\
Query Sample list: %s\n\
Ref. Original: %s\n\
Ref. Sample list: %s\n\
# vic.: %d\n\
Variant range: [%d-%d]\n\
k-mer vic: %d\n\
k-mer conc. log weight: %.3f\n\
N_e: %.4f\n\
# Query Subjects 2 Decode: %d\n", query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
ref_original_haplocoded_geno_fp, ref_sample_ids_fp, n_kmer_vicinity_vars,
var_start_i, var_end_i, n_kmer_vicinity_vars, kmer_concordance_log_weight, N_e, n_query_subjects_2_decode);

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	sort(query_proxized_geno_var_regs->begin(), query_proxized_geno_var_regs->end(), sort_regions);
	sort(query_original_geno_var_regs->begin(), query_original_geno_var_regs->end(), sort_regions);
	sort(ref_orig_geno_var_regs->begin(), ref_orig_geno_var_regs->end(), sort_regions);

	int max_query_proxy_geno = get_max_genotype_value(query_proxized_geno_var_regs, query_sample_ids);
	int max_query_original_geno = get_max_genotype_value(query_original_geno_var_regs, query_sample_ids);
	int max_ref_original_geno = get_max_genotype_value(ref_orig_geno_var_regs, ref_sample_ids);

	// Input checks.
	if (var_start_i <= n_kmer_vicinity_vars ||
		var_end_i >= (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars - 1))
	{
		fprintf(stderr, "Make sure variant start-end indices are within the kmer vicinity of variant ends: [%d-%d] vs [%d-%d]\n",
			var_start_i, var_end_i,
			n_kmer_vicinity_vars, (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars));
		exit(1);
	}

	if (max_query_proxy_geno != 3 ||
		max_query_original_geno != 3 ||
		max_ref_original_geno != 3)
	{
		fprintf(stderr, "One of the panels is not haplocoded..\n");
		exit(1);
	}

	if (query_proxized_geno_var_regs->size() != query_original_geno_var_regs->size() ||
		query_proxized_geno_var_regs->size() != ref_orig_geno_var_regs->size())
	{
		fprintf(stderr, "Variant counts are not the same: %d, %d, %d\n",
			vecsize(query_proxized_geno_var_regs), vecsize(query_original_geno_var_regs),
			vecsize(ref_orig_geno_var_regs));

		exit(1);
	}

	// Check to make sure sorted variant coordinates match.
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		if (ref_orig_geno_var_regs->at(i_var)->start != query_proxized_geno_var_regs->at(i_var)->start ||
			ref_orig_geno_var_regs->at(i_var)->start != query_original_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: The variant coordinates are not matching among panels.\n");
			exit(1);
		}
	} // i_var loop.

	//////////////////////////////////////////////////////////////////////////////////////////
	// Check & fix variant range.
	if (var_start_i < 1 ||
		var_start_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_start_i = 1;
	}

	if (var_end_i < 0 ||
		var_end_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_end_i = vecsize(query_proxized_geno_var_regs);
	}

	if (var_end_i < (var_start_i + 100))
	{
		fprintf(stderr, "Invalid range: %d-%d\n", var_start_i, var_end_i);
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	//////////////////////////////////////////////////////////////////////////////////////////
	// Assign the maximum frequency alleles to all reference variants.
	fprintf(stderr, "Assigning maximum frequency alleles to reference regions.\n");
	char max_AFs_bed_fp[1000];
	sprintf(max_AFs_bed_fp, "%s_max_AF_alleles.bed", op_prefix);
	FILE* f_max_AF_bed = open_f(max_AFs_bed_fp, "w");
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		void** cur_ref_var_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
		char* cur_ref_geno_sig = (char*)(cur_ref_var_info[0]);

		double total_AA = 0;
		for (int i_s = 0; i_s < vecsize(ref_sample_ids); i_s++)
		{
			total_AA += get_genotype_per_haplocoded_genotype(cur_ref_geno_sig[i_s]);
		} // i_s loop.

		double AAF = total_AA / (2 * vecsize(ref_sample_ids));

		if (AAF > 0.5)
		{
			ref_orig_geno_var_regs->at(i_var)->score = 1;
		}
		else
		{
			ref_orig_geno_var_regs->at(i_var)->score = 0;
		}

		fprintf(f_max_AF_bed, "%s\t%d\t%d\t%s\t%d\t+\n",
			ref_orig_geno_var_regs->at(i_var)->chrom,
			translate_coord(ref_orig_geno_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(ref_orig_geno_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			ref_orig_geno_var_regs->at(i_var)->name,
			ref_orig_geno_var_regs->at(i_var)->score);
	} // i_var loop.
	close_f(f_max_AF_bed, NULL);
	//////////////////////////////////////////////////////////////////////////////////////////

	// Set the decoded allele information for each query individual.
	fprintf(stderr, "Setting up the decoded genotype signals to original genotype signal regions.\n");
	for (int i_var = 0; i_var < vecsize(query_original_geno_var_regs); i_var++)
	{
		void** cur_var_info = (void**)(query_original_geno_var_regs->at(i_var)->data);

		char* decoded_geno_sig = new char[query_sample_ids->size() * 2];
		memset(decoded_geno_sig, 0, sizeof(char) * (query_sample_ids->size() * 2));

		void** new_var_info = new void* [10];
		new_var_info[0] = cur_var_info[0];
		new_var_info[1] = decoded_geno_sig;

		query_original_geno_var_regs->at(i_var)->data = new_var_info;
	} // i_var loop.
	//////////////////////////////////////////////////////////////////////////////////////////

	//fprintf(stderr, "Running Viterbi on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
	char cur_chr_recombination_rate_fp[1000];
	sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, ref_orig_geno_var_regs->at(0)->chrom);
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		exit(1);
	}

	// Add one element to the beginning of the regions for the current chromosome.
	// Assign the recomb rates.
	fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
	for (int i_reg = 0; i_reg < vecsize(ref_orig_geno_var_regs); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(ref_orig_geno_var_regs->at(i_reg), cur_chrom_recomb_regs);
		ref_orig_geno_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.

	int n_states = 2 * vecsize(ref_sample_ids);
	int n_ref_haplotypes = n_states;

	// We track a list of candidate kmers at each position.
	char dec_summary_stats_fp[1000];
	sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary.txt", op_prefix);
	FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	char per_sample_per_var_dec_stats_fp[1000];
	sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats.txt", op_prefix);
	FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Allocate arrays for storing scores and backtracking info.
	int n_vars = vecsize(ref_orig_geno_var_regs);

	// Allocate and initialize the forward/backward arrays.
	double*** ML_scores_per_hap = new double** [n_vars + 2];
	int*** ML_prev_state_per_hap = new int** [n_vars + 2];
	for (int var_i = 0; var_i <= n_vars + 1; var_i++)
	{
		ML_scores_per_hap[var_i] = new double* [2];
		ML_prev_state_per_hap[var_i] = new int* [2];

		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			ML_scores_per_hap[var_i][proxy_query_hap_i] = new double[n_states + 2];
			ML_prev_state_per_hap[var_i][proxy_query_hap_i] = new int[n_states + 2];
		} // proxy_query_hap_i loop.
	} // i loop.

	if (n_query_subjects_2_decode > vecsize(query_sample_ids) ||
		n_query_subjects_2_decode <= 0)
	{
		n_query_subjects_2_decode = vecsize(query_sample_ids);
	}






	// This is the max log-probability score for each unique kmer in the window.
	//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
	// # of unique ref. kmers in the window.
	int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];
	// List of unique ref. kmers per window.
	vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
	// This is the list of all ref kmers in the current window, we need this for keeping track of haplotype sequence overlaps.
	vector<t_kmer*>** per_win_all_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
	// The best scoring unique kmer index in the previous window for each current window.
	//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
	// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
	double*** per_win_per_query_proxy_per_hap_kmer_freq = new double** [query_proxized_geno_var_regs->size()];
	// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
	vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
	// This is the list of indices for each original ref kmer to map to unique kmers.
	vector<int>** per_win_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>*[query_proxized_geno_var_regs->size()];
	// This is the query hapfreq emission log probability for each unsorted reference kmer.
	vector<double>**** per_win_per_query_per_hap_per_ref_hap_query_emission_prob = new vector<double>***[query_proxized_geno_var_regs->size()];

	// Process all variants.
	fprintf(stderr, "Assigning kmer frequency and emission statistics..\n");
	for (int var_i = var_start_i - 1; var_i <= var_end_i; var_i++)
	{
		// Extract the k-mers on this window.
		// Set the window starts and ends.
		int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
		int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "================================\n");
			fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Extract the current window's kmers within reference panel.
		vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

		// Extract the current window's kmer for the query panel.
		vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "Extracted %d reference kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
		}
		//dump_kmers(cur_win_all_ref_kmers);

		/////////////////////////////////////////////////////////////////////////////////
		// Get the reference haplotype frequencies.
		// Get the window reference k-mers
		vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_ref_kmer_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_kmer_cnts);

		// This is necessary to map the reference k-mers to their probabilities.
		vector<int>* cur_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>();
		for (int i_kmer = 0; i_kmer < vecsize(cur_win_all_ref_kmers); i_kmer++)
		{
			for (int i_unique_kmer = 0; i_unique_kmer < vecsize(cur_win_unique_ref_kmers); i_unique_kmer++)
			{
				if (compare_kmers(cur_win_all_ref_kmers->at(i_kmer), cur_win_unique_ref_kmers->at(i_unique_kmer)))
				{
					cur_unique_ref_kmer_i_per_orig_ref_kmer_i->push_back(i_unique_kmer);
					break;
				}
			} // i_unique_kmer
		} // i_kmer loop.

		//delete_kmers(cur_win_all_ref_kmers);

		int n_unique_ref_haps = cur_win_unique_ref_kmer_cnts->size();

		// These keep state of ref kmer scores, backtracks, etc.
		per_win_per_unique_ref_kmer_cnts[var_i] = cur_win_unique_ref_kmer_cnts;
		per_win_n_ref_kmers[var_i] = n_unique_ref_haps;
		per_win_unique_ref_kmers[var_i] = cur_win_unique_ref_kmers;
		per_win_all_ref_kmers[var_i] = cur_win_all_ref_kmers;
		per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i] = cur_unique_ref_kmer_i_per_orig_ref_kmer_i;

		/////////////////////////////////////////////////////////////////////////////////
		// Get the query proxy k-mers.
		vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "%d unique reference kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
			fprintf(stderr, "%d unique proxy query kmers @ [%d, %d]\n", (int)cur_win_unique_query_proxy_kmers->size(), i_win_start, i_win_end);
		}

		per_win_per_query_proxy_per_hap_kmer_freq[var_i] = new double* [n_query_subjects_2_decode + 2];
		per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i] = new vector<double>**[n_query_subjects_2_decode];
		for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
		{
			per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s] = new double[2];
			per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s] = new vector<double>*[2];

			for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
			{
				// Get the current proxy haplotype of the query subject and extract the kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i];
				t_kmer* cur_query_subj_proxy_kmer = copy_kmer(cur_win_all_query_proxy_kmers->at(2 * proxy_query_i_s + proxy_query_hap_i));

				// This is a sanity check for the extracted proxy kmer.
				for (int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				}

				// Search for the query proxy kmer in unique kmers for the window.
				double cur_query_subj_proxy_kmer_cnt = 0.0;
				double total_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
						//break;
					}

					total_query_subj_proxy_kmer_cnt += cur_win_unique_query_proxy_cnts->at(i_kmer);
				} // i_kmer loop.
				delete_kmer(cur_query_subj_proxy_kmer);
				//delete_kmers(cur_win_unique_query_proxy_kmers);
				//delete(cur_win_unique_query_proxy_cnts);

				// This is the frequency of proxy haplotype at this window that we will match to from the reference panel.
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0 ||
					cur_query_subj_proxy_kmer_freq > 1.0)
				{
					fprintf(stderr, "Sanity Check failed: Found 0 frequency for the query proxy haplotype.\n");
					exit(1);
				}

				// Sanity check on the total count of the query kmer count.
				if (total_query_subj_proxy_kmer_cnt != (2 * vecsize(query_sample_ids)))
				{
					fprintf(stderr, "Sanity Check failed: Total query ref kmer counts does not add up to 1: %.3f/%d",
						total_query_subj_proxy_kmer_cnt, 2 * vecsize(query_sample_ids));

					exit(1);
				}

				// Save the frequency of this haplotype.
				//per_win_query_proxy_kmer_freq[var_i] = cur_query_subj_proxy_kmer_freq;
				per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s][proxy_query_hap_i] = cur_query_subj_proxy_kmer_freq;
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Calculate the normalized proxy query kmer emission probability for all of the reference kmers at this window.
				double query_emit_log_normalizing_factor = xlog(0.0);
				double total_ref_unique_kmer_prob = 0; // This is for sanity check.

				// This stores the list of emission probabilities for this proxy query kmer at the respective ref haplotype.
				vector<double>* per_state_query_emission_prob = new vector<double>();
				per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i] = per_state_query_emission_prob;
				//per_win_query_per_ref_hap_query_emission_prob[var_i] = per_state_query_emission_prob;

				for (int cur_ref_kmer_i = 0; cur_ref_kmer_i < n_ref_haplotypes; cur_ref_kmer_i++)
				{
					// First get the unique kmer index for this kmer.
					int cur_kmer_unique_kmer_i = cur_unique_ref_kmer_i_per_orig_ref_kmer_i->at(cur_ref_kmer_i);
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer_unique_kmer_i) / (2 * vecsize(ref_sample_ids));
					if (cur_ref_kmer_freq == 0)
					{
						fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
						exit(1);
					}

					if (!compare_kmers(cur_win_all_ref_kmers->at(cur_ref_kmer_i), cur_win_unique_ref_kmers->at(cur_kmer_unique_kmer_i)))
					{
						fprintf(stderr, "Sanity check failed: All-2-unique mapped k-mer does not match unique kmer.\n");
						exit(1);
					}

					total_ref_unique_kmer_prob += cur_ref_kmer_freq;

					double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
					double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

					// Add this emission probability to the list of emission probabilities.
					per_state_query_emission_prob->push_back(proxy_kmer_emit_log_prob_per_cur_kmer);

					// Update the normalizing factor.
					query_emit_log_normalizing_factor = xlog_sum(query_emit_log_normalizing_factor, proxy_kmer_emit_log_prob_per_cur_kmer);
				} // cur_kmer loop.

				// Normalize all of the emission probabilities.
				for (int cur_kmer = 0; cur_kmer < n_ref_haplotypes; cur_kmer++)
				{
					per_state_query_emission_prob->at(cur_kmer) = xlog_div(per_state_query_emission_prob->at(cur_kmer), query_emit_log_normalizing_factor);
				} // cur_kmer loop.

				// Note that this check does not hold any more since we are working over all of the haplotypes.
				//// Make sure that total kmer counts equal to 1.
				//if (fabs(total_ref_unique_kmer_prob - 1.0) > 0.001)
				//{
				//	fprintf(stderr, "Sanity check failed on total unique ref prob: %.4f\n", total_ref_unique_kmer_prob);
				//	exit(1);
				//}
				/////////////////////////////////////////////////////////////////////////////////
			}
		} // proxy_query_i_s option.

		delete_kmers(cur_win_all_query_proxy_kmers);
		delete_kmers(cur_win_unique_query_proxy_kmers);
	} // var_i loop.

	//////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Decoding %d/%d subjects.\n", n_query_subjects_2_decode, vecsize(query_sample_ids));

	double MIN_KMER_CONC_PROB = pow(10, -4);

	// Start looping over all individuals.
	for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
	{
		// Use all variants on the chromosome.
		vector<t_annot_region*>* cur_win_var_regs = ref_orig_geno_var_regs;

		// Start recursing over the variants.
		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			// Initialize all scores and backtracking states.
			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				// INitialize all the scores for all variants.
				memset(ML_scores_per_hap[var_i][proxy_query_hap_i], 0, sizeof(double) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = xlog(0.0);
				} // cur_state loop.

				// INitialize prev stats for all variants.
				memset(ML_prev_state_per_hap[var_i][proxy_query_hap_i], 0, sizeof(int) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = -1;
				} // cur_state loop.
			} // i loop.

			// Initialize the state probabilities for both haplotypes.
			// Initialize the first variant's scores.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				ML_scores_per_hap[var_start_i - 1][proxy_query_hap_i][state_i] = xlog((double)1.0 / n_ref_haplotypes);
			} // state_i loop.

			// This is the max log-probability score for each unique kmer in the window.
			//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
			// # of unique ref. kmers in the window.
			// Process all variants.
			for (int var_i = var_start_i; var_i < var_end_i; var_i++)
			{
				if (var_i % 10 == 0)
				{
					t_string::print_padded_string(stderr, '\r', 100, "Viterbi: sample_i: %d/%d: var_i: %d", proxy_query_i_s, vecsize(query_sample_ids), var_i);
				}

				// Extract the k-mers on this window.
				// Set the window starts and ends.
				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////
				// Pre-compute the transition probabilities:
				double prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				if (var_i > 0)
				{
					prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				}

				double cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
				if (var_i < n_vars)
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

				if (n_ref_haplotypes != vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]));
					exit(1);
				}

				if (n_ref_haplotypes != (2 * vecsize(ref_sample_ids)))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]));
					exit(1);
				}

				// Loop over all the states.
				double total_normalized_kmer_emission_log_prob = xlog(0.0); // This is for sanity check.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					// Recurse over the previous states.
					ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog(0.0);

					// This is a sanity check value on the emission probability of query kmer by ref kmers.
					total_normalized_kmer_emission_log_prob = xlog_sum(total_normalized_kmer_emission_log_prob, per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i]->at(state_i));

					for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
					{
						// Calculate the kmer concordance penalty.
						double kmer_concordance_log_prob = xlog(1.0);

						/////////////////
						// Get the kmer concordance penalty score for the previous state.
						double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i));

						kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

						// Exclude very low probability paths that add penalty of 10^-4.
						if (kmer_concordance_log_prob < log(MIN_KMER_CONC_PROB))
						{
							continue;
						}

						// We should not need this check any more since we buffered previous window's kmers.
						if (shifted_kmer_dist == 0 &&
							!are_kmers_shift_compatible(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i)))
						{
							fprintf(stderr, "var %d: %d -> %d is not k-mer compatible with 0 shifted distance.\n", var_i, prev_state_i, state_i);
							exit(1);
						}

						// We are currenty looping over all of the haplotypes (not unique kmers), we do not have exact frequency info for this yet.
						// We need to extract the frequency of the reference kmer at this position; however, we do not have access to this info, yet.
						int cur_unique_ref_kmer_i = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]->at(state_i);
						double cur_ref_kmer_freq = (double)per_win_per_unique_ref_kmer_cnts[var_i]->at(cur_unique_ref_kmer_i) / n_ref_haplotypes;
						if (cur_ref_kmer_freq == 0 || cur_ref_kmer_freq > 1.0)
						{
							fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is not valid: %.4f\n", cur_ref_kmer_freq);
							exit(1);
						}

						// Sanity check on the identified unique ref kmer index.
						for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
						{
							if (per_win_unique_ref_kmers[var_i]->at(cur_unique_ref_kmer_i)->kmer[nuc_i - i_win_start] != per_win_all_ref_kmers[var_i]->at(state_i)->kmer[nuc_i - i_win_start])
							{
								fprintf(stderr, "Ref kmer from unique kmers do not match all ref kmers.\n");
								exit(1);
							}
						} // nuc_i loop.
						//////////////////////////////////////////////////////////////////////////////////////////

						// Now update transition probability: This is added to Li-Stephens transition probability.
						//double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						//double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

						//// This is the normalized kmer emission prob score.
						//double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_div(proxy_kmer_emit_log_prob_per_cur_kmer, query_emit_log_normalizing_factor);

						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i]->at(state_i);
						//double norm_proxy_kmer_emit_log_prob_per_cur_kmer = per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i);

						//if (fabs(test_emit_val - norm_proxy_kmer_emit_log_prob_per_cur_kmer) > 0.001)
						//{
						//	fprintf(stderr, "Sanity check failed: Emission prob mismatch: %.4f vs %.4f\n",
						//		test_emit_val, norm_proxy_kmer_emit_log_prob_per_cur_kmer);
						//	exit(1);
						//}

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
						// Set the transition probabilities.
						double trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, kmer_concordance_log_prob);
						if (state_i == prev_state_i)
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(self_prob));
						}
						else
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(other_prob));
						}

						/*fprintf(stderr, "%d: %d->%d: trans prob: %.4fx %.4f / %.4f ;; %.4f\n", var_i,
							state_i, prev_state_i,
							exp(norm_proxy_kmer_emit_log_prob_per_cur_kmer), self_prob, other_prob,
							trans_prob);*/

							// Add the scaler for this position.
						ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = MAX(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i],
							xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i]));

						if (xlog_comp(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i], xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i])))
						{
							ML_prev_state_per_hap[var_i][proxy_query_hap_i][state_i] = prev_state_i;
						}

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "ML[%d][%d]: %.5f: P_trans_emit(%d->%d): %.5f\n",
								var_i, state_i, ML_scores_per_hap[var_i][proxy_query_hap_i][state_i],
								prev_state_i, state_i, trans_emit_prob);
						}
					} // prev_state_i loop.

					//fprintf(stderr, "Var %d: %d non-compatible kmers with previous haps.\n", var_i, state_i, n_non_compatible_kmers);
				} // state_i loop.

				// Following checks to make sure that we used a valid probability distribution for the emission of kmer's based on
				// haplotype frequencies.
				if (fabs(exp(total_normalized_kmer_emission_log_prob) - 1.0) > 0.001)
				{
					fprintf(stderr, "Sanity check failed: Total kmer emission freq is not valid: %.4f\n", total_normalized_kmer_emission_log_prob);
					exit(1);
				}

			} // var_i loop.

			// Compute the total log forward and backward probabilities.
			int top_scoring_hap_i = -1;
			double top_scoring_kmer_score = xlog(0.0);
			double least_scoring_kmer_score = xlog(1.0);

			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				top_scoring_kmer_score = MAX(top_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
				if (top_scoring_kmer_score == ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i])
				{
					top_scoring_hap_i = state_i;
				}

				least_scoring_kmer_score = MIN(least_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
			} // state_i loop.

			double total_max_proxy_score_per_full_hap = top_scoring_kmer_score;

			fprintf(stderr, "Haplotype %d total probabilities: ML=%.5f @ state %d\n", proxy_query_hap_i, total_max_proxy_score_per_full_hap, top_scoring_hap_i);

			// Trace the states back to identify the optimal path.
			fprintf(stderr, "Tracing back the Viterbi path for each haplotype.\n");
			//int** per_hap_per_variant_viterbi_haplotype = new int* [2];

			int cur_state = top_scoring_hap_i;
			int* per_variant_viterbi_haplotype = new int[n_vars + 2];
			int n_errs = 0;
			int n_max_AF_errs = 0;
			for (int var_i = var_end_i - 1; var_i >= var_start_i; var_i--)
			{
				fprintf(stderr, "=====================================\n");

				// Set the optimal state for the current index.
				per_variant_viterbi_haplotype[var_i] = cur_state;

				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				// Copy the original query subject's haplotype for the current query haplotype.
				char cur_win_query_subj_orig_kmer_nucs[100];
				memset(cur_win_query_subj_orig_kmer_nucs, 0, 100);
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_query_subj_orig_kmer_nucs[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i][nuc_i];
				} // nuc_i loop.

				char cur_win_ref_kmer_nucs[100];
				memset(cur_win_ref_kmer_nucs, 0, 100);
				int ref_i_s = cur_state / 2;
				int ref_i_hap = cur_state % 2;
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_ref_kmer_nucs[nuc_i - i_win_start] = per_ref_subj_orig_haplotypes->at(ref_i_s)[ref_i_hap][nuc_i];
				} // nuc_i loop.

				int l_kmer = i_win_end - i_win_start + 1;
				t_kmer* orig_query_kmer = copy_kmer(cur_win_query_subj_orig_kmer_nucs, l_kmer);
				t_kmer* ref_kmer = copy_kmer(cur_win_ref_kmer_nucs, l_kmer);

				// Get the query proxy and ref. kmer frequencies, these should be pretty much aligned.
				double query_proxy_kmer_freq = per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s][proxy_query_hap_i];

				int cur_ref_kmer_uniq_kmer_i = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]->at(cur_state);
				double backtracked_ref_hap_AF = (double)(per_win_per_unique_ref_kmer_cnts[var_i]->at(cur_ref_kmer_uniq_kmer_i)) / (2 * vecsize(ref_sample_ids));

				// Write.
				fprintf(stderr, "Variant %d: \n", var_i);

				fprintf(stderr, "Backtracked original k-mer @ var_i=%d [Query Proxy kmer freq: %.4f]\n", var_i, query_proxy_kmer_freq);
				dump_kmer(orig_query_kmer);
				fprintf(stderr, "Backtracked reference k-mer @ var_i=%d [Reference kmer freq: %.4f] (state: %d)\n", var_i, backtracked_ref_hap_AF, cur_state);
				dump_kmer(ref_kmer);

				bool found_kmer = false;
				//if (compare_kmers(ref_kmer, orig_query_kmer))
				if (ref_kmer->kmer[n_kmer_vicinity_vars] == orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					found_kmer = true;
				}
				else
				{
					n_errs++;
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////
				// Set the allele to this variant.
				char decoded_allele = ref_kmer->kmer[n_kmer_vicinity_vars];

				void** query_orig_reg_info = (void**)(query_original_geno_var_regs->at(var_i)->data);
				char* decoded_query_geno_sig = (char*)(query_orig_reg_info[1]);

				// Assign the decoded allele.
				int cur_geno = decoded_query_geno_sig[proxy_query_i_s];
				cur_geno = cur_geno | (decoded_allele << proxy_query_hap_i);
				decoded_query_geno_sig[proxy_query_i_s] = cur_geno;
				///////////////////////////////////////////////////////////////////////////////////////////////////

				// Check MAF AF errors.
				int max_AF_allele = ref_orig_geno_var_regs->at(var_i)->score;

				bool found_max_AF_allele = false;
				if (max_AF_allele != orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					n_max_AF_errs++;
				}
				else
				{
					found_max_AF_allele = true;
				}

				// Write the per subject info.
				fprintf(f_per_sample_per_var_dec_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n",
					query_sample_ids->at(proxy_query_i_s), proxy_query_hap_i,
					query_proxized_geno_var_regs->at(var_i)->chrom, query_proxized_geno_var_regs->at(var_i)->start,
					per_win_n_ref_kmers[var_i],
					cur_state,
					query_proxy_kmer_freq, backtracked_ref_hap_AF, (int)found_kmer, (int)found_max_AF_allele);
				fflush(f_per_sample_per_var_dec_stats);

				// Get the new state for the previous variant.
				cur_state = ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state];

				// Note that we have effectively moved one variant back.
				if (cur_state == -1)
				{
					if (var_i >= var_start_i)
					{
						fprintf(stderr, "Sanity check failed: Backtracking stopped before starting index.\n");
						exit(1);
					}

					fprintf(stderr, "Breaking @ variant %d.\n", var_i);
					break;
				}
				else
				{
					// Make sure that this state is valid.
					if (var_i > var_start_i &&
						cur_state > vecsize(per_win_all_ref_kmers[var_i - 1]))
					{
						fprintf(stderr, "Sanity check failed: Invalid backtracked state index @ var_i=%d.\n", var_i);
						exit(1);
					}
				}
			} // var_i loop.

			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fprintf(f_dec_stats, "%d\t%d\t%d\t%d\t%.5f\t%.5f\n", proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fflush(f_dec_stats);
		} // query_hap_i loop.

		// Save the genotypes for the decoded alleles on this subject.
		vector<t_annot_region*>* cur_sample_decoded_geno_regs = new vector<t_annot_region*>();
		for (int i_var = var_start_i; i_var < var_end_i; i_var++)
		{
			void** cur_var_reg_info = (void**)(query_original_geno_var_regs->at(i_var)->data);
			char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

			void** decoded_var_reg_info = new void* [3];
			char* decoded_geno_sig = new char[5];
			decoded_geno_sig[0] = var_reg_decoded_geno_sig[proxy_query_i_s];
			decoded_var_reg_info[0] = decoded_geno_sig;

			t_annot_region* decoded_geno_reg = duplicate_region(query_original_geno_var_regs->at(i_var));
			decoded_geno_reg->data = decoded_var_reg_info;
			cur_sample_decoded_geno_regs->push_back(decoded_geno_reg);
		} // i_var loop.

		// Set the sample id's and file name.
		vector<char*>* decoded_query_sample_ids = new vector<char*>();
		decoded_query_sample_ids->push_back(query_sample_ids->at(proxy_query_i_s));
		char decoded_sample_fp[1000];
		sprintf(decoded_sample_fp, "decoded_sample_list_%d.txt", proxy_query_i_s);
		FILE* f_decoded_sample_list = open_f(decoded_sample_fp, "w");
		fprintf(f_decoded_sample_list, "%s\n", query_sample_ids->at(proxy_query_i_s));
		close_f(f_decoded_sample_list, decoded_sample_fp);

		fprintf(stderr, "Saving decoded genotypes for the subject..");
		char decoded_query_geno_fp[1000];
		sprintf(decoded_query_geno_fp, "decoded_sample_%d.matbed.gz", proxy_query_i_s);
		binarize_variant_genotype_signal_regions(cur_sample_decoded_geno_regs, NULL, decoded_query_sample_ids, decoded_query_geno_fp);

		// Free memory for the decoded genotype regions on this subject.
		delete_annot_regions(cur_sample_decoded_geno_regs);
	} // proxy_query_i_s loop.

	// Save the decoded genotypes: Replace the decoded genotypes with actual genotypes.
	vector<t_annot_region*>* decoded_geno_regs = new vector<t_annot_region*>();
	for (int i_var = var_start_i; i_var < var_end_i; i_var++)
	{
		void** cur_var_reg_info = (void**)(query_original_geno_var_regs->at(i_var)->data);
		//char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

		// Replace the original variant signal to decoded genotype signal.
		cur_var_reg_info[0] = cur_var_reg_info[1];

		t_annot_region* decoded_geno_reg = duplicate_region(query_original_geno_var_regs->at(i_var));
		decoded_geno_reg->data = cur_var_reg_info;
		decoded_geno_regs->push_back(decoded_geno_reg);
	} // i_var loop.

	fprintf(stderr, "Saving all decoded genotypes for the subject..");
	char decoded_query_geno_fp[1000];
	sprintf(decoded_query_geno_fp, "decoded_all_sample.matbed.gz");
	binarize_variant_genotype_signal_regions(decoded_geno_regs, NULL, query_sample_ids, decoded_query_geno_fp);

	fprintf(stderr, "Done!\n");
} // decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_Buffered function.













static void* thread_callback_Viterbi_Decoder(void* __thread_info_ptr_ptr)
{
	void** __thread_info_ptr = (void**)(__thread_info_ptr_ptr);

	//////////////////////////////////////////////////////////////////////////////////////////
	int* int_pars = (int*)(__thread_info_ptr[0]);
	int thread_i = int_pars[0];
	int n_threads = int_pars[1];
	int n_query_subjects_2_decode = int_pars[2];
	int n_ref_haplotypes = int_pars[3];
	int var_start_i = int_pars[4];
	int var_end_i = int_pars[5];
	int n_kmer_vicinity_vars = int_pars[6];

	//////////////////////////////////////////////////////////////////////////////////////////
	double* dbl_pars = (double*)(__thread_info_ptr[1]);
	double N_e = dbl_pars[0];
	double kmer_concordance_log_weight = dbl_pars[1];
	//////////////////////////////////////////////////////////////////////////////////////////

	bool SAVE_IND_GENO_FILES = false;

	vector<t_annot_region*>* ref_orig_geno_var_regs = (vector<t_annot_region*>*)(__thread_info_ptr[2]);
	vector<t_annot_region*>* query_proxized_geno_var_regs = (vector<t_annot_region*>*)(__thread_info_ptr[3]);
	//vector<t_annot_region*>* query_original_geno_var_regs = (vector<t_annot_region*>*)(__thread_info_ptr[4]);

	vector<char*>* ref_sample_ids = (vector<char*>*)(__thread_info_ptr[5]);
	vector<char*>* query_sample_ids = (vector<char*>*)(__thread_info_ptr[6]);

	vector<char**>* per_query_subj_orig_haplotypes = (vector<char**>*)(__thread_info_ptr[7]);
	vector<char**>* per_ref_subj_orig_haplotypes = (vector<char**>*)(__thread_info_ptr[8]);

	int* per_win_n_ref_kmers = (int*)(__thread_info_ptr[9]);
	vector<t_kmer*>** per_win_unique_ref_kmers = (vector<t_kmer*>**)(__thread_info_ptr[10]);
	vector<t_kmer*>** per_win_all_ref_kmers = (vector<t_kmer*>**)(__thread_info_ptr[11]);
	double*** per_win_per_query_proxy_per_hap_kmer_freq = (double***)(__thread_info_ptr[12]);
	vector<int>** per_win_per_unique_ref_kmer_cnts = (vector<int>**)(__thread_info_ptr[13]);
	vector<int>** per_win_unique_ref_kmer_i_per_orig_ref_kmer_i = (vector<int>**)(__thread_info_ptr[14]);
	vector<double>**** per_win_per_query_per_hap_per_ref_hap_query_emission_prob = (vector<double>****)(__thread_info_ptr[15]);

	int*** per_win_per_ref_hap_transible_haps = (int***)(__thread_info_ptr[16]);

	char* op_prefix = (char*)(__thread_info_ptr[17]);

	int n_states = n_ref_haplotypes;
	int n_vars = vecsize(ref_orig_geno_var_regs);

	//////////////////////////////////////////////////////////////////////////////////////////
	t_string::print_padded_string(stderr, '\r', 100, "Thread %d: Decoding %d/%d subjects.                  \r", thread_i, n_query_subjects_2_decode, vecsize(query_sample_ids));

	double MIN_KMER_CONC_PROB = pow(10, -4);

	// Allocate and initialize the forward/backward arrays.
	double*** ML_scores_per_hap = new double** [n_vars + 2];
	int*** ML_prev_state_per_hap = new int** [n_vars + 2];

	// Allocate memory only for the variants of interest; this allocates a lot of memory.
	for (int var_i = (var_start_i-n_kmer_vicinity_vars); var_i <= (var_end_i+n_kmer_vicinity_vars); var_i++)
	{
		ML_scores_per_hap[var_i] = new double* [2];
		ML_prev_state_per_hap[var_i] = new int* [2];

		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			ML_scores_per_hap[var_i][proxy_query_hap_i] = new double[n_states + 2];
			ML_prev_state_per_hap[var_i][proxy_query_hap_i] = new int[n_states + 2];
		} // proxy_query_hap_i loop.
	} // i loop.

	// We track a list of candidate kmers at each position.
	char dec_summary_stats_fp[1000];
	sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary_thread_%d.txt", op_prefix, thread_i);
	FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	char per_sample_per_var_dec_stats_fp[1000];
	sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats_thread_%d.txt", op_prefix, thread_i);
	FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Start looping over all individuals.
	for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
	{
		if (proxy_query_i_s % n_threads != thread_i)
		{
			continue;
		}

		// Use all variants on the chromosome.
		vector<t_annot_region*>* cur_win_var_regs = ref_orig_geno_var_regs;

		// Start recursing over the variants.
		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			// Initialize all scores and backtracking states.
			for (int var_i = var_start_i-1; var_i <= var_end_i; var_i++)
			{
				// INitialize all the scores for all variants.
				memset(ML_scores_per_hap[var_i][proxy_query_hap_i], 0, sizeof(double) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = xlog(0.0);
				} // cur_state loop.

				// INitialize prev stats for all variants.
				memset(ML_prev_state_per_hap[var_i][proxy_query_hap_i], 0, sizeof(int) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = -1;
				} // cur_state loop.
			} // i loop.

			// Initialize the state probabilities for both haplotypes.
			// Initialize the first variant's scores.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				ML_scores_per_hap[var_start_i - 1][proxy_query_hap_i][state_i] = xlog((double)1.0 / n_ref_haplotypes);
			} // state_i loop.

			// This is the max log-probability score for each unique kmer in the window.
			//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
			// # of unique ref. kmers in the window.
			// Process all variants.
			for (int var_i = var_start_i; var_i < var_end_i; var_i++)
			{
				if (var_i % 10 == 0)
				{
					t_string::print_padded_string(stderr, '\r', 100, "Viterbi Thread %d: sample_i: %d/%d: var_i: %d", thread_i, proxy_query_i_s, vecsize(query_sample_ids), var_i);
				}

				// Extract the k-mers on this window.
				// Set the window starts and ends.
				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////
				// Pre-compute the transition probabilities:
				double prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				if (var_i > 0)
				{
					prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				}

				double cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
				if (var_i < n_vars)
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

				if (n_ref_haplotypes != vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]));
					exit(1);
				}

				if (n_ref_haplotypes != (2 * vecsize(ref_sample_ids)))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]));
					exit(1);
				}

				// Loop over all the states.
				double total_normalized_kmer_emission_log_prob = xlog(0.0); // This is for sanity check.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					// Recurse over the previous states.
					ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog(0.0);

					// This is a sanity check value on the emission probability of query kmer by ref kmers.
					total_normalized_kmer_emission_log_prob = xlog_sum(total_normalized_kmer_emission_log_prob, per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i]->at(state_i));

					
					//for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
					for (int prev_state_i_transible_i = 0; prev_state_i_transible_i < n_ref_haplotypes; prev_state_i_transible_i++)
					{
						int prev_state_i = per_win_per_ref_hap_transible_haps[var_i][state_i][prev_state_i_transible_i];

						// Are we at the end of the transible previous states for this state?
						if (prev_state_i == -1)
						{
							break;
						}

						// Calculate the kmer concordance penalty.
						double kmer_concordance_log_prob = xlog(1.0);

						/////////////////
						// Get the kmer concordance penalty score for the previous state.
						double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i));

						kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

						// Exclude very low probability paths that add penalty of 10^-4.
						if (kmer_concordance_log_prob < log(MIN_KMER_CONC_PROB))
						{
							continue;
						}

						// We should not need this check any more since we buffered previous window's kmers.
						if (shifted_kmer_dist == 0 &&
							!are_kmers_shift_compatible(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i)))
						{
							fprintf(stderr, "var %d: %d -> %d is not k-mer compatible with 0 shifted distance.\n", var_i, prev_state_i, state_i);
							exit(1);
						}

						// We are currenty looping over all of the haplotypes (not unique kmers), we do not have exact frequency info for this yet.
						// We need to extract the frequency of the reference kmer at this position; however, we do not have access to this info, yet.
						int cur_unique_ref_kmer_i = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]->at(state_i);
						double cur_ref_kmer_freq = (double)per_win_per_unique_ref_kmer_cnts[var_i]->at(cur_unique_ref_kmer_i) / n_ref_haplotypes;
						if (cur_ref_kmer_freq == 0 || cur_ref_kmer_freq > 1.0)
						{
							fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is not valid: %.4f\n", cur_ref_kmer_freq);
							exit(1);
						}

						// Sanity check on the identified unique ref kmer index.
						for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
						{
							if (per_win_unique_ref_kmers[var_i]->at(cur_unique_ref_kmer_i)->kmer[nuc_i - i_win_start] != per_win_all_ref_kmers[var_i]->at(state_i)->kmer[nuc_i - i_win_start])
							{
								fprintf(stderr, "Ref kmer from unique kmers do not match all ref kmers.\n");
								exit(1);
							}
						} // nuc_i loop.
						//////////////////////////////////////////////////////////////////////////////////////////

						// Now update transition probability: This is added to Li-Stephens transition probability.
						//double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						//double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

						//// This is the normalized kmer emission prob score.
						//double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_div(proxy_kmer_emit_log_prob_per_cur_kmer, query_emit_log_normalizing_factor);

						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i]->at(state_i);
						//double norm_proxy_kmer_emit_log_prob_per_cur_kmer = per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i);

						//if (fabs(test_emit_val - norm_proxy_kmer_emit_log_prob_per_cur_kmer) > 0.001)
						//{
						//	fprintf(stderr, "Sanity check failed: Emission prob mismatch: %.4f vs %.4f\n",
						//		test_emit_val, norm_proxy_kmer_emit_log_prob_per_cur_kmer);
						//	exit(1);
						//}

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
						// Set the transition probabilities.
						double trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, kmer_concordance_log_prob);
						if (state_i == prev_state_i)
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(self_prob));
						}
						else
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(other_prob));
						}

						/*fprintf(stderr, "%d: %d->%d: trans prob: %.4fx %.4f / %.4f ;; %.4f\n", var_i,
							state_i, prev_state_i,
							exp(norm_proxy_kmer_emit_log_prob_per_cur_kmer), self_prob, other_prob,
							trans_prob);*/

							// Add the scaler for this position.
						ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = MAX(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i],
							xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i]));

						if (xlog_comp(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i], xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i])))
						{
							ML_prev_state_per_hap[var_i][proxy_query_hap_i][state_i] = prev_state_i;
						}

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "ML[%d][%d]: %.5f: P_trans_emit(%d->%d): %.5f\n",
								var_i, state_i, ML_scores_per_hap[var_i][proxy_query_hap_i][state_i],
								prev_state_i, state_i, trans_emit_prob);
						}
					} // prev_state_i loop.

					//fprintf(stderr, "Var %d: %d non-compatible kmers with previous haps.\n", var_i, state_i, n_non_compatible_kmers);
				} // state_i loop.

				// Following checks to make sure that we used a valid probability distribution for the emission of kmer's based on
				// haplotype frequencies.
				if (fabs(exp(total_normalized_kmer_emission_log_prob) - 1.0) > 0.001)
				{
					fprintf(stderr, "Sanity check failed: Total kmer emission freq is not valid: %.4f\n", total_normalized_kmer_emission_log_prob);
					exit(1);
				}

			} // var_i loop.

			// Compute the total log forward and backward probabilities.
			int top_scoring_hap_i = -1;
			double top_scoring_kmer_score = xlog(0.0);
			double least_scoring_kmer_score = xlog(1.0);

			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				top_scoring_kmer_score = MAX(top_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
				if (top_scoring_kmer_score == ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i])
				{
					top_scoring_hap_i = state_i;
				}

				least_scoring_kmer_score = MIN(least_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
			} // state_i loop.

			double total_max_proxy_score_per_full_hap = top_scoring_kmer_score;

			fprintf(stderr, "Haplotype %d total probabilities: ML=%.5f @ state %d\n", proxy_query_hap_i, total_max_proxy_score_per_full_hap, top_scoring_hap_i);

			// Trace the states back to identify the optimal path.
			fprintf(stderr, "Tracing back the Viterbi path for each haplotype.\n");
			//int** per_hap_per_variant_viterbi_haplotype = new int* [2];

			int cur_state = top_scoring_hap_i;
			int* per_variant_viterbi_haplotype = new int[n_vars + 2];
			int n_errs = 0;
			int n_max_AF_errs = 0;
			for (int var_i = var_end_i - 1; var_i >= var_start_i; var_i--)
			{
				// Set the optimal state for the current index.
				per_variant_viterbi_haplotype[var_i] = cur_state;

				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				// Copy the original query subject's haplotype for the current query haplotype.
				char cur_win_query_subj_orig_kmer_nucs[100];
				memset(cur_win_query_subj_orig_kmer_nucs, 0, 100);
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_query_subj_orig_kmer_nucs[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i][nuc_i];
				} // nuc_i loop.

				char cur_win_ref_kmer_nucs[100];
				memset(cur_win_ref_kmer_nucs, 0, 100);
				int ref_i_s = cur_state / 2;
				int ref_i_hap = cur_state % 2;
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_ref_kmer_nucs[nuc_i - i_win_start] = per_ref_subj_orig_haplotypes->at(ref_i_s)[ref_i_hap][nuc_i];
				} // nuc_i loop.

				int l_kmer = i_win_end - i_win_start + 1;
				t_kmer* orig_query_kmer = copy_kmer(cur_win_query_subj_orig_kmer_nucs, l_kmer);
				t_kmer* ref_kmer = copy_kmer(cur_win_ref_kmer_nucs, l_kmer);

				// Get the query proxy and ref. kmer frequencies, these should be pretty much aligned.
				double query_proxy_kmer_freq = per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s][proxy_query_hap_i];

				int cur_ref_kmer_uniq_kmer_i = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]->at(cur_state);
				double backtracked_ref_hap_AF = (double)(per_win_per_unique_ref_kmer_cnts[var_i]->at(cur_ref_kmer_uniq_kmer_i)) / (2 * vecsize(ref_sample_ids));

				// Write.
				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "=====================================\n");

					fprintf(stderr, "Variant %d: \n", var_i);

					fprintf(stderr, "Backtracked original k-mer @ var_i=%d [Query Proxy kmer freq: %.4f]\n", var_i, query_proxy_kmer_freq);
					dump_kmer(orig_query_kmer);
					fprintf(stderr, "Backtracked reference k-mer @ var_i=%d [Reference kmer freq: %.4f] (state: %d)\n", var_i, backtracked_ref_hap_AF, cur_state);
					dump_kmer(ref_kmer);
				}

				bool found_kmer = false;
				//if (compare_kmers(ref_kmer, orig_query_kmer))
				if (ref_kmer->kmer[n_kmer_vicinity_vars] == orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					found_kmer = true;
				}
				else
				{
					n_errs++;
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////
				// Set the allele to this variant.
				char decoded_allele = ref_kmer->kmer[n_kmer_vicinity_vars];

				void** ref_orig_reg_info = (void**)(ref_orig_geno_var_regs->at(var_i)->data);
				char* decoded_query_geno_sig = (char*)(ref_orig_reg_info[1]);

				// Assign the decoded allele.
				int cur_geno = decoded_query_geno_sig[proxy_query_i_s];
				cur_geno = cur_geno | (decoded_allele << proxy_query_hap_i);
				decoded_query_geno_sig[proxy_query_i_s] = cur_geno;
				///////////////////////////////////////////////////////////////////////////////////////////////////

				// Check MAF AF errors.
				int max_AF_allele = ref_orig_geno_var_regs->at(var_i)->score;

				bool found_max_AF_allele = false;
				if (max_AF_allele != orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					n_max_AF_errs++;
				}
				else
				{
					found_max_AF_allele = true;
				}

				// Write the per subject info.
				fprintf(f_per_sample_per_var_dec_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n",
					query_sample_ids->at(proxy_query_i_s), proxy_query_hap_i,
					ref_orig_geno_var_regs->at(var_i)->chrom, ref_orig_geno_var_regs->at(var_i)->start,
					per_win_n_ref_kmers[var_i],
					cur_state,
					query_proxy_kmer_freq, backtracked_ref_hap_AF, (int)found_kmer, (int)found_max_AF_allele);
				fflush(f_per_sample_per_var_dec_stats);

				// Get the new state for the previous variant.
				cur_state = ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state];

				// Note that we have effectively moved one variant back.
				if (cur_state == -1)
				{
					if (var_i >= var_start_i)
					{
						fprintf(stderr, "Sanity check failed: Backtracking stopped before starting index.\n");
						exit(1);
					}

					fprintf(stderr, "Breaking @ variant %d.\n", var_i);
					break;
				}
				else
				{
					// Make sure that this state is valid.
					if (var_i > var_start_i &&
						cur_state > vecsize(per_win_all_ref_kmers[var_i - 1]))
					{
						fprintf(stderr, "Sanity check failed: Invalid backtracked state index @ var_i=%d.\n", var_i);
						exit(1);
					}
				}
			} // var_i loop.

			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fprintf(f_dec_stats, "%d\t%d\t%d\t%d\t%.5f\t%.5f\n", proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fflush(f_dec_stats);
		} // query_hap_i loop.

		if (SAVE_IND_GENO_FILES)
		{
			// Save the genotypes for the decoded alleles on this subject.
			vector<t_annot_region*>* cur_sample_decoded_geno_regs = new vector<t_annot_region*>();
			for (int i_var = var_start_i; i_var < var_end_i; i_var++)
			{
				void** cur_var_reg_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
				char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

				void** decoded_var_reg_info = new void* [3];
				char* decoded_geno_sig = new char[5];
				decoded_geno_sig[0] = var_reg_decoded_geno_sig[proxy_query_i_s];
				decoded_var_reg_info[0] = decoded_geno_sig;

				t_annot_region* decoded_geno_reg = duplicate_region(ref_orig_geno_var_regs->at(i_var));
				decoded_geno_reg->data = decoded_var_reg_info;
				cur_sample_decoded_geno_regs->push_back(decoded_geno_reg);
			} // i_var loop.

			// Set the sample id's and file name.
			vector<char*>* decoded_query_sample_ids = new vector<char*>();
			decoded_query_sample_ids->push_back(query_sample_ids->at(proxy_query_i_s));
			char decoded_sample_fp[1000];
			sprintf(decoded_sample_fp, "decoded_sample_list_%d.txt", proxy_query_i_s);
			FILE* f_decoded_sample_list = open_f(decoded_sample_fp, "w");
			fprintf(f_decoded_sample_list, "%s\n", query_sample_ids->at(proxy_query_i_s));
			close_f(f_decoded_sample_list, decoded_sample_fp);

			fprintf(stderr, "Saving decoded genotypes for the subject..");
			char decoded_query_geno_fp[1000];
			sprintf(decoded_query_geno_fp, "decoded_sample_%d.matbed.gz", proxy_query_i_s);
			binarize_variant_genotype_signal_regions(cur_sample_decoded_geno_regs, NULL, decoded_query_sample_ids, decoded_query_geno_fp);

			// Free memory for the decoded genotype regions on this subject.
			delete_annot_regions(cur_sample_decoded_geno_regs);
		} // individual genotype saving check.
	} // proxy_query_i_s loop.

	close_f(f_dec_stats, dec_summary_stats_fp);
	close_f(f_per_sample_per_var_dec_stats, per_sample_per_var_dec_stats_fp);

	return NULL;
} // thread_callback_Viterbi_Decoder function.

static void* thread_callback_Transible_Calculator(void* __thread_info_ptr_ptr)
{
	void** __thread_info_ptr = (void**)(__thread_info_ptr_ptr);

	//////////////////////////////////////////////////////////////////////////////////////////
	int* int_pars = (int*)(__thread_info_ptr[0]);
	int thread_i = int_pars[0];
	int n_threads = int_pars[1];
	//int n_query_subjects_2_decode = int_pars[2];
	int n_ref_haplotypes = int_pars[3];
	int var_start_i = int_pars[4];
	int var_end_i = int_pars[5];
	//int n_kmer_vicinity_vars = int_pars[6];

	//////////////////////////////////////////////////////////////////////////////////////////
	double* dbl_pars = (double*)(__thread_info_ptr[1]);
	double MIN_KMER_CONC_PROB = dbl_pars[0];
	double kmer_concordance_log_weight = dbl_pars[1];

	//////////////////////////////////////////////////////////////////////////////////////////
	int*** per_win_per_ref_hap_transible_haps = (int***)(__thread_info_ptr[2]);
	//vector<t_annot_region*>* query_original_geno_var_regs = (vector<t_annot_region*>*)(__thread_info_ptr[3]);
	vector<t_kmer*>** per_win_all_ref_kmers = (vector<t_kmer*>**)(__thread_info_ptr[4]);

	//double MIN_KMER_CONC_PROB = xlog(0.0001);

	for (int var_i = var_start_i; var_i <= var_end_i; var_i++)
	{
		if (var_i % n_threads != thread_i)
		{
			continue;
		}

		if (var_i % 100 == 0)
		{
			fprintf(stderr, "Thread %d: Setting transibles @ variant %d        \r", thread_i, var_i);
		}

		per_win_per_ref_hap_transible_haps[var_i] = new int* [n_ref_haplotypes + 2];
		for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
		{
			int* cur_state_transible_prev_state = new int[n_ref_haplotypes + 2];
			per_win_per_ref_hap_transible_haps[var_i][state_i] = cur_state_transible_prev_state;

			int n_cur_transibles = 0;
			for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
			{
				// Calculate the kmer concordance penalty.
				double kmer_concordance_log_prob = xlog(1.0);

				// Get the kmer concordance penalty score for the previous state.
				double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i));

				kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

				// Exclude very low probability paths that add penalty of 10^-4.
				if (kmer_concordance_log_prob >= log(MIN_KMER_CONC_PROB))
				{
					cur_state_transible_prev_state[n_cur_transibles] = prev_state_i;
					n_cur_transibles++;
				}
			} // prev_state_i loop.

			// Set the final transible to -1 to indicate end.
			cur_state_transible_prev_state[n_cur_transibles] = -1;
		} // state_i loop.
	} // var_i loop.

	return NULL;
} // thread_callback_Transible_Calculator function.

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi_MT(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	int n_threads,
	char* op_prefix)
{
	fprintf(stderr, "Decoding proxized alleles of query site using histogram matching HMM with Viterbi-MT:\n\
Query Original: %s\n\
Query Proxy: %s\n\
Query Sample list: %s\n\
Ref. Original: %s\n\
Ref. Sample list: %s\n\
# vic.: %d\n\
Variant range: [%d-%d]\n\
k-mer vic: %d\n\
k-mer conc. log weight: %.3f\n\
N_e: %.4f\n\
# threads: %d\n\
# Query Subjects 2 Decode: %d\n", query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
ref_original_haplocoded_geno_fp, ref_sample_ids_fp, n_kmer_vicinity_vars,
var_start_i, var_end_i, n_kmer_vicinity_vars, kmer_concordance_log_weight, N_e, n_threads, n_query_subjects_2_decode);

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	sort(query_proxized_geno_var_regs->begin(), query_proxized_geno_var_regs->end(), sort_regions);
	sort(query_original_geno_var_regs->begin(), query_original_geno_var_regs->end(), sort_regions);
	sort(ref_orig_geno_var_regs->begin(), ref_orig_geno_var_regs->end(), sort_regions);

	int max_query_proxy_geno = get_max_genotype_value(query_proxized_geno_var_regs, query_sample_ids);
	int max_query_original_geno = get_max_genotype_value(query_original_geno_var_regs, query_sample_ids);
	int max_ref_original_geno = get_max_genotype_value(ref_orig_geno_var_regs, ref_sample_ids);

	// Input checks.
	if (var_start_i <= n_kmer_vicinity_vars ||
		var_end_i >= (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars - 1))
	{
		fprintf(stderr, "Make sure variant start-end indices are within the kmer vicinity of variant ends: [%d-%d] vs [%d-%d]\n",
			var_start_i, var_end_i,
			n_kmer_vicinity_vars, (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars));
		exit(1);
	}

	if (max_query_proxy_geno != 3 ||
		max_query_original_geno != 3 ||
		max_ref_original_geno != 3)
	{
		fprintf(stderr, "One of the panels is not haplocoded..\n");
		exit(1);
	}

	if (query_proxized_geno_var_regs->size() != query_original_geno_var_regs->size() ||
		query_proxized_geno_var_regs->size() != ref_orig_geno_var_regs->size())
	{
		fprintf(stderr, "Variant counts are not the same: %d, %d, %d\n",
			vecsize(query_proxized_geno_var_regs), vecsize(query_original_geno_var_regs),
			vecsize(ref_orig_geno_var_regs));

		exit(1);
	}

	// Check to make sure sorted variant coordinates match.
	// Note that this does not enforce a matching for the reference any more.
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		/*if (ref_orig_geno_var_regs->at(i_var)->start != query_proxized_geno_var_regs->at(i_var)->start ||
			ref_orig_geno_var_regs->at(i_var)->start != query_original_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: The variant coordinates are not matching among panels.\n");
			exit(1);
		}*/

		if (query_original_geno_var_regs->at(i_var)->start != query_proxized_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: The variant coordinates are not matching among panels.\n");
			exit(1);
		}
	} // i_var loop.

	//////////////////////////////////////////////////////////////////////////////////////////
	// Check & fix variant range.
	if (var_start_i < 1 ||
		var_start_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_start_i = 1;
	}

	if (var_end_i < 0 ||
		var_end_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_end_i = vecsize(query_proxized_geno_var_regs);
	}

	if (var_end_i < (var_start_i + 100))
	{
		fprintf(stderr, "Invalid range: %d-%d\n", var_start_i, var_end_i);
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	//////////////////////////////////////////////////////////////////////////////////////////
	// Assign the maximum frequency alleles to all reference variants.
	fprintf(stderr, "Assigning maximum frequency alleles to reference regions.\n");
	char max_AFs_bed_fp[1000];
	sprintf(max_AFs_bed_fp, "%s_max_AF_alleles.bed", op_prefix);
	FILE* f_max_AF_bed = open_f(max_AFs_bed_fp, "w");
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		void** cur_ref_var_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
		char* cur_ref_geno_sig = (char*)(cur_ref_var_info[0]);

		double total_AA = 0;
		for (int i_s = 0; i_s < vecsize(ref_sample_ids); i_s++)
		{
			total_AA += get_genotype_per_haplocoded_genotype(cur_ref_geno_sig[i_s]);
		} // i_s loop.

		double AAF = total_AA / (2 * vecsize(ref_sample_ids));

		if (AAF > 0.5)
		{
			ref_orig_geno_var_regs->at(i_var)->score = 1;
		}
		else
		{
			ref_orig_geno_var_regs->at(i_var)->score = 0;
		}

		fprintf(f_max_AF_bed, "%s\t%d\t%d\t%s\t%d\t+\n",
			ref_orig_geno_var_regs->at(i_var)->chrom,
			translate_coord(ref_orig_geno_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(ref_orig_geno_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			ref_orig_geno_var_regs->at(i_var)->name,
			ref_orig_geno_var_regs->at(i_var)->score);
	} // i_var loop.
	close_f(f_max_AF_bed, NULL);
	//////////////////////////////////////////////////////////////////////////////////////////

	// Set the decoded allele information for each query individual.
	fprintf(stderr, "Setting up the decoded genotype signals to original genotype signal regions.\n");
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		void** cur_var_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);

		char* decoded_geno_sig = new char[query_sample_ids->size() * 2];
		memset(decoded_geno_sig, 0, sizeof(char) * (query_sample_ids->size() * 2));

		void** new_var_info = new void* [10];
		new_var_info[0] = cur_var_info[0];
		new_var_info[1] = decoded_geno_sig;

		ref_orig_geno_var_regs->at(i_var)->data = new_var_info;
	} // i_var loop.
	//////////////////////////////////////////////////////////////////////////////////////////

	//fprintf(stderr, "Running Viterbi on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
	char cur_chr_recombination_rate_fp[1000];
	sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, ref_orig_geno_var_regs->at(0)->chrom);
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		exit(1);
	}

	// Add one element to the beginning of the regions for the current chromosome.
	// Assign the recomb rates.
	fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
	for (int i_reg = 0; i_reg < vecsize(ref_orig_geno_var_regs); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(ref_orig_geno_var_regs->at(i_reg), cur_chrom_recomb_regs);
		ref_orig_geno_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.

	int n_states = 2 * vecsize(ref_sample_ids);
	int n_ref_haplotypes = n_states;

	//// We track a list of candidate kmers at each position.
	//char dec_summary_stats_fp[1000];
	//sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary.txt", op_prefix);
	//FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	//char per_sample_per_var_dec_stats_fp[1000];
	//sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats.txt", op_prefix);
	//FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Allocate arrays for storing scores and backtracking info.
	//int n_vars = vecsize(ref_orig_geno_var_regs);

	if (n_query_subjects_2_decode > vecsize(query_sample_ids) ||
		n_query_subjects_2_decode <= 0)
	{
		n_query_subjects_2_decode = vecsize(query_sample_ids);
	}

	// This is the max log-probability score for each unique kmer in the window.
	//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
	// # of unique ref. kmers in the window.
	int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];
	// List of unique ref. kmers per window.
	vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
	// This is the list of all ref kmers in the current window, we need this for keeping track of haplotype sequence overlaps.
	vector<t_kmer*>** per_win_all_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
	// The best scoring unique kmer index in the previous window for each current window.
	//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
	// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
	double*** per_win_per_query_proxy_per_hap_kmer_freq = new double**[query_proxized_geno_var_regs->size()];
	// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
	vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
	// This is the list of indices for each original ref kmer to map to unique kmers.
	vector<int>** per_win_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>*[query_proxized_geno_var_regs->size()];
	// This is the query hapfreq emission log probability for each unsorted reference kmer.
	vector<double>**** per_win_per_query_per_hap_per_ref_hap_query_emission_prob = new vector<double>***[query_proxized_geno_var_regs->size()];

	// Process all variants.
	fprintf(stderr, "Assigning kmer frequency and emission statistics..\n");
	for (int var_i = var_start_i - 1; var_i <= var_end_i; var_i++)
	{
		// Extract the k-mers on this window.
		// Set the window starts and ends.
		int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
		int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "================================\n");
			fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
		}

		/////////////////////////////////////////////////////////////////////////////////
		// Extract the current window's kmers within reference panel.
		vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

		// Extract the current window's kmer for the query panel.
		vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "Extracted %d reference kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
		}
		//dump_kmers(cur_win_all_ref_kmers);

		/////////////////////////////////////////////////////////////////////////////////
		// Get the reference haplotype frequencies.
		// Get the window reference k-mers
		vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_ref_kmer_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_kmer_cnts);

		// This is necessary to map the reference k-mers to their probabilities.
		vector<int>* cur_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>();
		for (int i_kmer = 0; i_kmer < vecsize(cur_win_all_ref_kmers); i_kmer++)
		{
			for (int i_unique_kmer = 0; i_unique_kmer < vecsize(cur_win_unique_ref_kmers); i_unique_kmer++)
			{
				if (compare_kmers(cur_win_all_ref_kmers->at(i_kmer), cur_win_unique_ref_kmers->at(i_unique_kmer)))
				{
					cur_unique_ref_kmer_i_per_orig_ref_kmer_i->push_back(i_unique_kmer);
					break;
				}
			} // i_unique_kmer
		} // i_kmer loop.

		//delete_kmers(cur_win_all_ref_kmers);

		int n_unique_ref_haps = cur_win_unique_ref_kmer_cnts->size();

		// These keep state of ref kmer scores, backtracks, etc.
		per_win_per_unique_ref_kmer_cnts[var_i] = cur_win_unique_ref_kmer_cnts;
		per_win_n_ref_kmers[var_i] = n_unique_ref_haps;
		per_win_unique_ref_kmers[var_i] = cur_win_unique_ref_kmers;
		per_win_all_ref_kmers[var_i] = cur_win_all_ref_kmers;
		per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i] = cur_unique_ref_kmer_i_per_orig_ref_kmer_i;

		/////////////////////////////////////////////////////////////////////////////////
		// Get the query proxy k-mers.
		vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
		vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
		get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

		if (__DUMP_KMER_TRACKING_MSGS__)
		{
			fprintf(stderr, "%d unique reference kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
			fprintf(stderr, "%d unique proxy query kmers @ [%d, %d]\n", (int)cur_win_unique_query_proxy_kmers->size(), i_win_start, i_win_end);
		}

		per_win_per_query_proxy_per_hap_kmer_freq[var_i] = new double* [n_query_subjects_2_decode + 2];
		per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i] = new vector<double>**[n_query_subjects_2_decode];
		for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
		{
			per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s] = new double[2];
			per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s] = new vector<double>*[2];

			for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
			{
				// Get the current proxy haplotype of the query subject and extract the kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i];
				t_kmer* cur_query_subj_proxy_kmer = copy_kmer(cur_win_all_query_proxy_kmers->at(2 * proxy_query_i_s + proxy_query_hap_i));

				// This is a sanity check for the extracted proxy kmer.
				for (int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				}

				// Search for the query proxy kmer in unique kmers for the window.
				double cur_query_subj_proxy_kmer_cnt = 0.0;
				double total_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
						//break;
					}

					total_query_subj_proxy_kmer_cnt += cur_win_unique_query_proxy_cnts->at(i_kmer);
				} // i_kmer loop.
				delete_kmer(cur_query_subj_proxy_kmer);
				//delete_kmers(cur_win_unique_query_proxy_kmers);
				//delete(cur_win_unique_query_proxy_cnts);

				// This is the frequency of proxy haplotype at this window that we will match to from the reference panel.
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0 ||
					cur_query_subj_proxy_kmer_freq > 1.0)
				{
					fprintf(stderr, "Sanity Check failed: Found 0 frequency for the query proxy haplotype.\n");
					exit(1);
				}

				// Sanity check on the total count of the query kmer count.
				if (total_query_subj_proxy_kmer_cnt != (2 * vecsize(query_sample_ids)))
				{
					fprintf(stderr, "Sanity Check failed: Total query ref kmer counts does not add up to 1: %.3f/%d",
						total_query_subj_proxy_kmer_cnt, 2 * vecsize(query_sample_ids));

					exit(1);
				}

				// Save the frequency of this haplotype.
				//per_win_query_proxy_kmer_freq[var_i] = cur_query_subj_proxy_kmer_freq;
				per_win_per_query_proxy_per_hap_kmer_freq[var_i][proxy_query_i_s][proxy_query_hap_i] = cur_query_subj_proxy_kmer_freq;
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Calculate the normalized proxy query kmer emission probability for all of the reference kmers at this window.
				double query_emit_log_normalizing_factor = xlog(0.0);
				double total_ref_unique_kmer_prob = 0; // This is for sanity check.

				// This stores the list of emission probabilities for this proxy query kmer at the respective ref haplotype.
				vector<double>* per_state_query_emission_prob = new vector<double>();
				per_win_per_query_per_hap_per_ref_hap_query_emission_prob[var_i][proxy_query_i_s][proxy_query_hap_i] = per_state_query_emission_prob;
				//per_win_query_per_ref_hap_query_emission_prob[var_i] = per_state_query_emission_prob;

				for (int cur_ref_kmer_i = 0; cur_ref_kmer_i < n_ref_haplotypes; cur_ref_kmer_i++)
				{
					// First get the unique kmer index for this kmer.
					int cur_kmer_unique_kmer_i = cur_unique_ref_kmer_i_per_orig_ref_kmer_i->at(cur_ref_kmer_i);
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer_unique_kmer_i) / (2 * vecsize(ref_sample_ids));
					if (cur_ref_kmer_freq == 0)
					{
						fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
						exit(1);
					}

					if (!compare_kmers(cur_win_all_ref_kmers->at(cur_ref_kmer_i), cur_win_unique_ref_kmers->at(cur_kmer_unique_kmer_i)))
					{
						fprintf(stderr, "Sanity check failed: All-2-unique mapped k-mer does not match unique kmer.\n");
						exit(1);
					}

					total_ref_unique_kmer_prob += cur_ref_kmer_freq;

					double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
					double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

					// Add this emission probability to the list of emission probabilities.
					per_state_query_emission_prob->push_back(proxy_kmer_emit_log_prob_per_cur_kmer);

					// Update the normalizing factor.
					query_emit_log_normalizing_factor = xlog_sum(query_emit_log_normalizing_factor, proxy_kmer_emit_log_prob_per_cur_kmer);
				} // cur_kmer loop.

				// Normalize all of the emission probabilities.
				for (int cur_kmer = 0; cur_kmer < n_ref_haplotypes; cur_kmer++)
				{
					per_state_query_emission_prob->at(cur_kmer) = xlog_div(per_state_query_emission_prob->at(cur_kmer), query_emit_log_normalizing_factor);
				} // cur_kmer loop.

				// Note that this check does not hold any more since we are working over all of the haplotypes.
				//// Make sure that total kmer counts equal to 1.
				//if (fabs(total_ref_unique_kmer_prob - 1.0) > 0.001)
				//{
				//	fprintf(stderr, "Sanity check failed on total unique ref prob: %.4f\n", total_ref_unique_kmer_prob);
				//	exit(1);
				//}
				/////////////////////////////////////////////////////////////////////////////////
			}
		} // proxy_query_i_s option.

		delete_kmers(cur_win_all_query_proxy_kmers);
		delete_kmers(cur_win_unique_query_proxy_kmers);
	} // var_i loop.

	double MIN_KMER_CONC_PROB = pow(10, -4);

	fprintf(stderr, "Setting transibles array.\n");
	int*** per_win_per_ref_hap_transible_haps = new int** [vecsize(query_original_geno_var_regs) + 2];
	vector<t_ansi_thread*>* transible_comp_threads = new vector<t_ansi_thread*>();
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		void** __thread_info_ptr = new void* [50];

		//////////////////////////////////////////////////////////////////////////////////////////
		//int* int_pars = (int*)(__thread_info_ptr[0]);
		int* int_pars = new int[20];
		int_pars[0] = i_thread;
		int_pars[1] = n_threads;
		int_pars[2] = n_query_subjects_2_decode;
		int_pars[3] = n_ref_haplotypes;
		int_pars[4] = var_start_i;
		int_pars[5] = var_end_i;
		int_pars[6] = n_kmer_vicinity_vars;
		__thread_info_ptr[0] = int_pars;

		//////////////////////////////////////////////////////////////////////////////////////////
		double* dbl_pars = new double[10];
		dbl_pars[0] = MIN_KMER_CONC_PROB;
		dbl_pars[1] = kmer_concordance_log_weight;
		__thread_info_ptr[1] = dbl_pars;

		//////////////////////////////////////////////////////////////////////////////////////////
		__thread_info_ptr[2] = per_win_per_ref_hap_transible_haps;
		__thread_info_ptr[3] = query_original_geno_var_regs;
		__thread_info_ptr[4] = per_win_all_ref_kmers;

		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_Transible_Calculator, __thread_info_ptr);
		cur_thread->run_thread();
		transible_comp_threads->push_back(cur_thread);
	} // i_thread loop.

	fprintf(stderr, "Waiting for threads..\n");
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		transible_comp_threads->at(i_thread)->wait_thread();
	} // i_thread loop.
	fprintf(stderr, "Threads finished..\n");

	// Start the threads.
	vector<t_ansi_thread*>* threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting thread %d..         \r", thread_i);
		void** cur_thread_info_ptr = new void*[20];

		//////////////////////////////////////////////////////////////////////////////////////////
		int* int_pars = new int[20];
		cur_thread_info_ptr[0] = int_pars;

		//int* int_pars = (int*)(__thread_info_ptr[0]);
		int_pars[0] = thread_i;
		int_pars[1] = n_threads;
		int_pars[2] = n_query_subjects_2_decode;
		int_pars[3] = n_ref_haplotypes;
		int_pars[4] = var_start_i;
		int_pars[5] = var_end_i;
		int_pars[6] = n_kmer_vicinity_vars;

		//////////////////////////////////////////////////////////////////////////////////////////
		double* dbl_pars = new double[20];
		cur_thread_info_ptr[1] = dbl_pars;
		dbl_pars[0] = N_e;
		dbl_pars[1] = kmer_concordance_log_weight;

		//////////////////////////////////////////////////////////////////////////////////////////

		cur_thread_info_ptr[2] = ref_orig_geno_var_regs;
		cur_thread_info_ptr[3] = query_proxized_geno_var_regs;
		cur_thread_info_ptr[4] = query_original_geno_var_regs;

		cur_thread_info_ptr[5] = ref_sample_ids;
		cur_thread_info_ptr[6] = query_sample_ids;

		cur_thread_info_ptr[7] = per_query_subj_orig_haplotypes;
		cur_thread_info_ptr[8] = per_ref_subj_orig_haplotypes;

		cur_thread_info_ptr[9] = per_win_n_ref_kmers;
		cur_thread_info_ptr[10] = per_win_unique_ref_kmers;
		cur_thread_info_ptr[11] = per_win_all_ref_kmers;
		cur_thread_info_ptr[12] = per_win_per_query_proxy_per_hap_kmer_freq;
		cur_thread_info_ptr[13] = per_win_per_unique_ref_kmer_cnts;
		cur_thread_info_ptr[14] = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i;
		cur_thread_info_ptr[15] = per_win_per_query_per_hap_per_ref_hap_query_emission_prob;
		cur_thread_info_ptr[16] = per_win_per_ref_hap_transible_haps;

		cur_thread_info_ptr[17] = op_prefix;

		t_ansi_thread* cur_thread = new t_ansi_thread(thread_callback_Viterbi_Decoder, cur_thread_info_ptr);

		cur_thread->run_thread();
		fprintf(stderr, "Started thread %d..         \r", thread_i);

		threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "\nWaiting for threads to finish.\n");
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		threads->at(i_thread)->wait_thread();
	} // i_thread loop.

	// Save the decoded genotypes: Replace the decoded genotypes with actual genotypes.
	vector<t_annot_region*>* decoded_geno_regs = new vector<t_annot_region*>();
	for (int i_var = var_start_i; i_var < var_end_i; i_var++)
	{
		void** cur_var_reg_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
		//char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

		// Replace the original variant signal to decoded genotype signal.
		cur_var_reg_info[0] = cur_var_reg_info[1];

		t_annot_region* decoded_geno_reg = duplicate_region(ref_orig_geno_var_regs->at(i_var));
		decoded_geno_reg->data = cur_var_reg_info;
		decoded_geno_regs->push_back(decoded_geno_reg);
	} // i_var loop.

	fprintf(stderr, "Saving all decoded genotypes for the subject..");
	char decoded_query_geno_fp[1000];
	sprintf(decoded_query_geno_fp, "%s_decoded_ref_genotypes.matbed.gz", op_prefix);
	binarize_variant_genotype_signal_regions(decoded_geno_regs, NULL, query_sample_ids, decoded_query_geno_fp);

	fprintf(stderr, "Done!\n");
} // decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi function.

void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	char* op_prefix)
{
	fprintf(stderr, "Decoding proxized alleles of query site using histogram matching HMM with Viterbi:\n\
Query Original: %s\n\
Query Proxy: %s\n\
Query Sample list: %s\n\
Ref. Original: %s\n\
Ref. Sample list: %s\n\
# vic.: %d\n\
Variant range: [%d-%d]\n\
k-mer vic: %d\n\
k-mer conc. log weight: %.3f\n\
N_e: %.4f\n\
# Query Subjects 2 Decode: %d\n", query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
ref_original_haplocoded_geno_fp, ref_sample_ids_fp, n_kmer_vicinity_vars,
var_start_i, var_end_i, n_kmer_vicinity_vars, kmer_concordance_log_weight, N_e, n_query_subjects_2_decode);

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	sort(query_proxized_geno_var_regs->begin(), query_proxized_geno_var_regs->end(), sort_regions);
	sort(query_original_geno_var_regs->begin(), query_original_geno_var_regs->end(), sort_regions);
	sort(ref_orig_geno_var_regs->begin(), ref_orig_geno_var_regs->end(), sort_regions);

	int max_query_proxy_geno = get_max_genotype_value(query_proxized_geno_var_regs, query_sample_ids);
	int max_query_original_geno = get_max_genotype_value(query_original_geno_var_regs, query_sample_ids);
	int max_ref_original_geno = get_max_genotype_value(ref_orig_geno_var_regs, ref_sample_ids);

	// Input checks.
	if (var_start_i <= n_kmer_vicinity_vars ||
		var_end_i >= (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars - 1))
	{
		fprintf(stderr, "Make sure variant start-end indices are within the kmer vicinity of variant ends: [%d-%d] vs [%d-%d]\n",
			var_start_i, var_end_i,
			n_kmer_vicinity_vars, (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars));
		exit(1);
	}

	if (max_query_proxy_geno != 3 ||
		max_query_original_geno != 3 ||
		max_ref_original_geno != 3)
	{
		fprintf(stderr, "One of the panels is not haplocoded..\n");
		exit(1);
	}

	if (query_proxized_geno_var_regs->size() != query_original_geno_var_regs->size() ||
		query_proxized_geno_var_regs->size() != ref_orig_geno_var_regs->size())
	{
		fprintf(stderr, "Variant counts are not the same: %d, %d, %d\n",
			vecsize(query_proxized_geno_var_regs), vecsize(query_original_geno_var_regs),
			vecsize(ref_orig_geno_var_regs));

		exit(1);
	}

	// Check to make sure sorted variant coordinates match.
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		if (ref_orig_geno_var_regs->at(i_var)->start != query_proxized_geno_var_regs->at(i_var)->start ||
			ref_orig_geno_var_regs->at(i_var)->start != query_original_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: The variant coordinates are not matching among panels.\n");
			exit(1);
		}
	} // i_var loop.

	//////////////////////////////////////////////////////////////////////////////////////////
	// Check & fix variant range.
	if (var_start_i < 1 ||
		var_start_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_start_i = 1;
	}

	if (var_end_i < 0 ||
		var_end_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_end_i = vecsize(query_proxized_geno_var_regs);
	}

	if (var_end_i < (var_start_i + 100))
	{
		fprintf(stderr, "Invalid range: %d-%d\n", var_start_i, var_end_i);
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	//////////////////////////////////////////////////////////////////////////////////////////
	// Assign the maximum frequency alleles to all reference variants.
	fprintf(stderr, "Assigning maximum frequency alleles to reference regions.\n");
	char max_AFs_bed_fp[1000];
	sprintf(max_AFs_bed_fp, "%s_max_AF_alleles.bed", op_prefix);
	FILE* f_max_AF_bed = open_f(max_AFs_bed_fp, "w");
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		void** cur_ref_var_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
		char* cur_ref_geno_sig = (char*)(cur_ref_var_info[0]);

		double total_AA = 0;
		for (int i_s = 0; i_s < vecsize(ref_sample_ids); i_s++)
		{
			total_AA += get_genotype_per_haplocoded_genotype(cur_ref_geno_sig[i_s]);
		} // i_s loop.

		double AAF = total_AA / (2 * vecsize(ref_sample_ids));

		if (AAF > 0.5)
		{
			ref_orig_geno_var_regs->at(i_var)->score = 1;
		}
		else
		{
			ref_orig_geno_var_regs->at(i_var)->score = 0;
		}

		fprintf(f_max_AF_bed, "%s\t%d\t%d\t%s\t%d\t+\n",
			ref_orig_geno_var_regs->at(i_var)->chrom,
			translate_coord(ref_orig_geno_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(ref_orig_geno_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			ref_orig_geno_var_regs->at(i_var)->name,
			ref_orig_geno_var_regs->at(i_var)->score);
	} // i_var loop.
	close_f(f_max_AF_bed, NULL);
	//////////////////////////////////////////////////////////////////////////////////////////

	// Set the decoded allele information for each query individual.
	fprintf(stderr, "Setting up the decoded genotype signals to original genotype signal regions.\n");
	for (int i_var = 0; i_var < vecsize(query_original_geno_var_regs); i_var++)
	{
		void** cur_var_info = (void**)(query_original_geno_var_regs->at(i_var)->data);

		char* decoded_geno_sig = new char[query_sample_ids->size() * 2];
		memset(decoded_geno_sig, 0, sizeof(char)* (query_sample_ids->size() * 2));

		void** new_var_info = new void* [10];
		new_var_info[0] = cur_var_info[0];
		new_var_info[1] = decoded_geno_sig;

		query_original_geno_var_regs->at(i_var)->data = new_var_info;
	} // i_var loop.
	//////////////////////////////////////////////////////////////////////////////////////////

	//fprintf(stderr, "Running Viterbi on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
	char cur_chr_recombination_rate_fp[1000];
	sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, ref_orig_geno_var_regs->at(0)->chrom);
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		exit(1);
	}

	// Add one element to the beginning of the regions for the current chromosome.
	// Assign the recomb rates.
	fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
	for (int i_reg = 0; i_reg < vecsize(ref_orig_geno_var_regs); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(ref_orig_geno_var_regs->at(i_reg), cur_chrom_recomb_regs);
		ref_orig_geno_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.

	int n_states = 2 * vecsize(ref_sample_ids);
	int n_ref_haplotypes = n_states;

	// We track a list of candidate kmers at each position.
	char dec_summary_stats_fp[1000];
	sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary.txt", op_prefix);
	FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	char per_sample_per_var_dec_stats_fp[1000];
	sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats.txt", op_prefix);
	FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Allocate arrays for storing scores and backtracking info.
	int n_vars = vecsize(ref_orig_geno_var_regs);

	// Allocate and initialize the forward/backward arrays.
	double*** ML_scores_per_hap = new double** [n_vars + 2];
	int*** ML_prev_state_per_hap = new int** [n_vars + 2];
	for (int var_i = 0; var_i <= n_vars + 1; var_i++)
	{
		ML_scores_per_hap[var_i] = new double* [2];
		ML_prev_state_per_hap[var_i] = new int* [2];

		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			ML_scores_per_hap[var_i][proxy_query_hap_i] = new double[n_states + 2];
			ML_prev_state_per_hap[var_i][proxy_query_hap_i] = new int[n_states + 2];
		} // proxy_query_hap_i loop.
	} // i loop.

	if (n_query_subjects_2_decode > vecsize(query_sample_ids) ||
		n_query_subjects_2_decode <= 0)
	{
		n_query_subjects_2_decode = vecsize(query_sample_ids);
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	fprintf(stderr, "Decoding %d/%d subjects.\n", n_query_subjects_2_decode, vecsize(query_sample_ids));

	double MIN_KMER_CONC_PROB = pow(10, -4);

	// Start looping over all individuals.
	for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
	{
		// Use all variants on the chromosome.
		vector<t_annot_region*>* cur_win_var_regs = ref_orig_geno_var_regs;

		// Start recursing over the variants.
		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			// Initialize all scores and backtracking states.
			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				// INitialize all the scores for all variants.
				memset(ML_scores_per_hap[var_i][proxy_query_hap_i], 0, sizeof(double) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = xlog(0.0);
				} // cur_state loop.

				// INitialize prev stats for all variants.
				memset(ML_prev_state_per_hap[var_i][proxy_query_hap_i], 0, sizeof(int) * (n_states + 1));
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state] = -1;
				} // cur_state loop.
			} // i loop.

			// Initialize the state probabilities for both haplotypes.
			// Initialize the first variant's scores.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				ML_scores_per_hap[var_start_i - 1][proxy_query_hap_i][state_i] = xlog((double)1.0 / n_ref_haplotypes);
			} // state_i loop.

			// This is the max log-probability score for each unique kmer in the window.
			//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
			// # of unique ref. kmers in the window.
			int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];
			// List of unique ref. kmers per window.
			vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
			// This is the list of all ref kmers in the current window, we need this for keeping track of haplotype sequence overlaps.
			vector<t_kmer*>** per_win_all_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
			// The best scoring unique kmer index in the previous window for each current window.
			//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
			// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
			double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];
			// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
			vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
			// This is the list of indices for each original ref kmer to map to unique kmers.
			vector<int>** per_win_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>*[query_proxized_geno_var_regs->size()];
			// This is the query hapfreq emission log probability for each unsorted reference kmer.
			vector<double>** per_win_query_per_ref_hap_query_emission_prob = new vector<double>*[query_proxized_geno_var_regs->size()];

			// Process all variants.
			for (int var_i = var_start_i - 1; var_i <= var_end_i; var_i++)
			{
				if (var_i % 10 == 0)
				{
					fprintf(stderr, "Viterbi: sample_i: %d/%d: var_i: %d         \r", proxy_query_i_s, vecsize(query_sample_ids), var_i);
				}

				// Extract the k-mers on this window.
				// Set the window starts and ends.
				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////////////////////////////////
				// Extract the current window's kmers within reference panel.
				vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

				// Extract the current window's kmer for the query panel.
				vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d reference kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
				}
				//dump_kmers(cur_win_all_ref_kmers);

				/////////////////////////////////////////////////////////////////////////////////
				// Get the reference haplotype frequencies.
				// Get the window reference k-mers
				vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_ref_kmer_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_kmer_cnts);

				// This is necessary to map the reference k-mers to their probabilities.
				vector<int>* cur_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>();
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_all_ref_kmers); i_kmer++)
				{
					for (int i_unique_kmer = 0; i_unique_kmer < vecsize(cur_win_unique_ref_kmers); i_unique_kmer++)
					{
						if (compare_kmers(cur_win_all_ref_kmers->at(i_kmer), cur_win_unique_ref_kmers->at(i_unique_kmer)))
						{
							cur_unique_ref_kmer_i_per_orig_ref_kmer_i->push_back(i_unique_kmer);
							break;
						}
					} // i_unique_kmer
				} // i_kmer loop.

				//delete_kmers(cur_win_all_ref_kmers);

				int n_unique_ref_haps = cur_win_unique_ref_kmer_cnts->size();

				// These keep state of ref kmer scores, backtracks, etc.
				per_win_per_unique_ref_kmer_cnts[var_i] = cur_win_unique_ref_kmer_cnts;
				per_win_n_ref_kmers[var_i] = n_unique_ref_haps;
				//per_win_per_ref_kmer_log_cumul_probs[var_i] = new double[n_unique_ref_haps];
				//per_win_backtracking_kmer_i[var_i] = new int[n_unique_ref_haps];
				per_win_unique_ref_kmers[var_i] = cur_win_unique_ref_kmers;
				per_win_all_ref_kmers[var_i] = cur_win_all_ref_kmers;
				per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i] = cur_unique_ref_kmer_i_per_orig_ref_kmer_i;

				//// Initialize the score and backtracking haplotype index in previous window.
				//for (int i_kmer = 0; i_kmer < n_unique_ref_haps; i_kmer++)
				//{
				//	//per_win_per_ref_kmer_log_cumul_probs[var_i][i_kmer] = xlog(0.0);
				//	//per_win_backtracking_kmer_i[var_i][i_kmer] = -1;
				//} // i_kmer loop.
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Get the query proxy k-mers.
				vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d unique reference kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
					fprintf(stderr, "%d unique proxy query kmers @ [%d, %d]\n", (int)cur_win_unique_query_proxy_kmers->size(), i_win_start, i_win_end);
				}

				// Get the current proxy haplotype of the query subject and extract the kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i];
				t_kmer* cur_query_subj_proxy_kmer = copy_kmer(cur_win_all_query_proxy_kmers->at(2 * proxy_query_i_s + proxy_query_hap_i));
				delete_kmers(cur_win_all_query_proxy_kmers);

				// This is a sanity check for the extracted proxy kmer.
				for (int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				}

				// Search for the query proxy kmer in unique kmers for the window.
				double cur_query_subj_proxy_kmer_cnt = 0.0;
				double total_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
						//break;
					}

					total_query_subj_proxy_kmer_cnt += cur_win_unique_query_proxy_cnts->at(i_kmer);
				} // i_kmer loop.
				delete_kmer(cur_query_subj_proxy_kmer);
				delete_kmers(cur_win_unique_query_proxy_kmers);
				delete(cur_win_unique_query_proxy_cnts);

				// This is the frequency of proxy haplotype at this window that we will match to from the reference panel.
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0 ||
					cur_query_subj_proxy_kmer_freq > 1.0)
				{
					fprintf(stderr, "Sanity Check failed: Found 0 frequency for the query proxy haplotype.\n");
					exit(1);
				}

				// Sanity check on the total count of the query kmer count.
				if (total_query_subj_proxy_kmer_cnt != (2 * vecsize(query_sample_ids)))
				{
					fprintf(stderr, "Sanity Check failed: Total query ref kmer counts does not add up to 1: %.3f/%d",
						total_query_subj_proxy_kmer_cnt, 2 * vecsize(query_sample_ids));

					exit(1);
				}

				// Save the frequency of this haplotype.
				per_win_query_proxy_kmer_freq[var_i] = cur_query_subj_proxy_kmer_freq;
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Calculate the normalized proxy query kmer emission probability for all of the reference kmers at this window.
				double query_emit_log_normalizing_factor = xlog(0.0);
				double total_ref_unique_kmer_prob = 0; // This is for sanity check.

				// This stores the list of emission probabilities for this proxy query kmer at the respective ref haplotype.
				vector<double>* per_state_query_emission_prob = new vector<double>();
				per_win_query_per_ref_hap_query_emission_prob[var_i] = per_state_query_emission_prob;

				for (int cur_ref_kmer_i = 0; cur_ref_kmer_i < n_ref_haplotypes; cur_ref_kmer_i++)
				{
					// First get the unique kmer index for this kmer.
					int cur_kmer_unique_kmer_i = cur_unique_ref_kmer_i_per_orig_ref_kmer_i->at(cur_ref_kmer_i);
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer_unique_kmer_i) / (2 * vecsize(ref_sample_ids));
					if (cur_ref_kmer_freq == 0)
					{
						fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
						exit(1);
					}

					if (!compare_kmers(cur_win_all_ref_kmers->at(cur_ref_kmer_i), cur_win_unique_ref_kmers->at(cur_kmer_unique_kmer_i)))
					{
						fprintf(stderr, "Sanity check failed: All-2-unique mapped k-mer does not match unique kmer.\n");
						exit(1);
					}

					total_ref_unique_kmer_prob += cur_ref_kmer_freq;

					double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
					double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

					// Add this emission probability to the list of emission probabilities.
					per_state_query_emission_prob->push_back(proxy_kmer_emit_log_prob_per_cur_kmer);

					// Update the normalizing factor.
					query_emit_log_normalizing_factor = xlog_sum(query_emit_log_normalizing_factor, proxy_kmer_emit_log_prob_per_cur_kmer);
				} // cur_kmer loop.

				// Normalize all of the emission probabilities.
				for (int cur_kmer = 0; cur_kmer < n_ref_haplotypes; cur_kmer++)
				{
					per_state_query_emission_prob->at(cur_kmer) = xlog_div(per_state_query_emission_prob->at(cur_kmer), query_emit_log_normalizing_factor);
				} // cur_kmer loop.

				// Note that this check does not hold any more since we are working over all of the haplotypes.
				//// Make sure that total kmer counts equal to 1.
				//if (fabs(total_ref_unique_kmer_prob - 1.0) > 0.001)
				//{
				//	fprintf(stderr, "Sanity check failed on total unique ref prob: %.4f\n", total_ref_unique_kmer_prob);
				//	exit(1);
				//}
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////
				// Pre-compute the transition probabilities:
				double prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				if (var_i > 0)
				{
					prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
				}

				double cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
				if (var_i < n_vars)
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

				if (n_ref_haplotypes != vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i));
					exit(1);
				}

				if (n_ref_haplotypes != (2 * vecsize(ref_sample_ids)))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i));
					exit(1);
				}

				////////////////////////////////////////////////////////////////////////////////////
				// Check boundary condition.
				if (var_i < var_start_i)
				{
					continue;
				}
				////////////////////////////////////////////////////////////////////////////////////

				// Loop over all the states.
				double total_normalized_kmer_emission_log_prob = xlog(0.0); // This is for sanity check.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					// Recurse over the previous states.
					ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog(0.0);

					// This is a sanity check value on the emission probability of query kmer by ref kmers.
					total_normalized_kmer_emission_log_prob = xlog_sum(total_normalized_kmer_emission_log_prob, per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i));

					for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
					{
						// Calculate the kmer concordance penalty.
						double kmer_concordance_log_prob = xlog(1.0);

						/////////////////
						// Get the kmer concordance penalty score for the previous state.
						double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i));

						kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

						// Exclude very low probability paths that add penalty of 10^-4.
						if (kmer_concordance_log_prob < log(MIN_KMER_CONC_PROB))
						{
							continue;
						}

						// We should not need this check any more since we buffered previous window's kmers.
						if(shifted_kmer_dist == 0 &&
							!are_kmers_shift_compatible(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i)))
						{
							fprintf(stderr, "var %d: %d -> %d is not k-mer compatible with 0 shifted distance.\n", var_i, prev_state_i, state_i);
							exit(1);
						}

						// We are currenty looping over all of the haplotypes (not unique kmers), we do not have exact frequency info for this yet.
						// We need to extract the frequency of the reference kmer at this position; however, we do not have access to this info, yet.
						int cur_unique_ref_kmer_i = cur_unique_ref_kmer_i_per_orig_ref_kmer_i->at(state_i);
						double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_unique_ref_kmer_i) / n_ref_haplotypes;
						if (cur_ref_kmer_freq == 0 || cur_ref_kmer_freq > 1.0)
						{
							fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is not valid: %.4f\n", cur_ref_kmer_freq);
							exit(1);
						}

						// Sanity check on the identified unique ref kmer index.
						for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
						{
							if (cur_win_unique_ref_kmers->at(cur_unique_ref_kmer_i)->kmer[nuc_i - i_win_start] != cur_win_all_ref_kmers->at(state_i)->kmer[nuc_i - i_win_start])
							{
								fprintf(stderr, "Ref kmer from unique kmers do not match all ref kmers.\n");
								exit(1);
							}
						} // nuc_i loop.
						//////////////////////////////////////////////////////////////////////////////////////////

						// Now update transition probability: This is added to Li-Stephens transition probability.
						double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
						double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

						// This is the normalized kmer emission prob score.
						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_div(proxy_kmer_emit_log_prob_per_cur_kmer, query_emit_log_normalizing_factor);

						double test_emit_val = per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i);

						if (fabs(test_emit_val - norm_proxy_kmer_emit_log_prob_per_cur_kmer) > 0.001)
						{
							fprintf(stderr, "Sanity check failed: Emission prob mismatch: %.4f vs %.4f\n", 
								test_emit_val, norm_proxy_kmer_emit_log_prob_per_cur_kmer);
							exit(1);
						}

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
						// Set the transition probabilities.
						double trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, kmer_concordance_log_prob);
						if (state_i == prev_state_i)
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(self_prob));
						}
						else
						{
							trans_emit_prob = xlog_mul(trans_emit_prob, xlog(other_prob));
						}

						/*fprintf(stderr, "%d: %d->%d: trans prob: %.4fx %.4f / %.4f ;; %.4f\n", var_i,
							state_i, prev_state_i,
							exp(norm_proxy_kmer_emit_log_prob_per_cur_kmer), self_prob, other_prob,
							trans_prob);*/

						// Add the scaler for this position.
						ML_scores_per_hap[var_i][proxy_query_hap_i][state_i] = MAX(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i], 
																					xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i]));

						if (xlog_comp(ML_scores_per_hap[var_i][proxy_query_hap_i][state_i], xlog_mul(trans_emit_prob, ML_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i])))
						{
							ML_prev_state_per_hap[var_i][proxy_query_hap_i][state_i] = prev_state_i;
						}

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "ML[%d][%d]: %.5f: P_trans_emit(%d->%d): %.5f\n",
								var_i, state_i, ML_scores_per_hap[var_i][proxy_query_hap_i][state_i],
								prev_state_i, state_i, trans_emit_prob);
						}
					} // prev_state_i loop.

					//fprintf(stderr, "Var %d: %d non-compatible kmers with previous haps.\n", var_i, state_i, n_non_compatible_kmers);
				} // state_i loop.


				// Following checks to make sure that we used a valid probability distribution for the emission of kmer's based on
				// haplotype frequencies.
				if (fabs(exp(total_normalized_kmer_emission_log_prob) - 1.0) > 0.001)
				{
					fprintf(stderr, "Sanity check failed: Total kmer emission freq is not valud: %.4f\n", total_normalized_kmer_emission_log_prob);
					exit(1);
				}
			} // var_i loop.


			// Compute the total log forward and backward probabilities.
			int top_scoring_hap_i = -1;
			double top_scoring_kmer_score = xlog(0.0);
			double least_scoring_kmer_score = xlog(1.0);
			
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				top_scoring_kmer_score = MAX(top_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
				if (top_scoring_kmer_score == ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i])
				{
					top_scoring_hap_i = state_i;
				}

				least_scoring_kmer_score = MIN(least_scoring_kmer_score, ML_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]);
			} // state_i loop.

			double total_max_proxy_score_per_full_hap = top_scoring_kmer_score;

			fprintf(stderr, "Haplotype %d total probabilities: ML=%.5f @ state %d\n", proxy_query_hap_i, total_max_proxy_score_per_full_hap, top_scoring_hap_i);

			// Trace the states back to identify the optimal path.
			fprintf(stderr, "Tracing back the Viterbi path for each haplotype.\n");
			//int** per_hap_per_variant_viterbi_haplotype = new int* [2];

			int cur_state = top_scoring_hap_i;
			int* per_variant_viterbi_haplotype = new int[n_vars + 2];
			int n_errs = 0;
			int n_max_AF_errs = 0;
			for (int var_i = var_end_i-1; var_i >= var_start_i; var_i--)
			{
				fprintf(stderr, "=====================================\n");

				// Set the optimal state for the current index.
				per_variant_viterbi_haplotype[var_i] = cur_state;

				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				// Copy the original query subject's haplotype for the current query haplotype.
				char cur_win_query_subj_orig_kmer_nucs[100];
				memset(cur_win_query_subj_orig_kmer_nucs, 0, 100);
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_query_subj_orig_kmer_nucs[nuc_i - i_win_start] = per_query_subj_orig_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i][nuc_i];
				} // nuc_i loop.

				char cur_win_ref_kmer_nucs[100];
				memset(cur_win_ref_kmer_nucs, 0, 100);
				int ref_i_s = cur_state / 2;
				int ref_i_hap = cur_state % 2;
				for (int nuc_i = i_win_start; nuc_i <= i_win_end; nuc_i++)
				{
					cur_win_ref_kmer_nucs[nuc_i - i_win_start] = per_ref_subj_orig_haplotypes->at(ref_i_s)[ref_i_hap][nuc_i];
				} // nuc_i loop.

				int l_kmer = i_win_end - i_win_start + 1;
				t_kmer* orig_query_kmer = copy_kmer(cur_win_query_subj_orig_kmer_nucs, l_kmer);
				t_kmer* ref_kmer = copy_kmer(cur_win_ref_kmer_nucs, l_kmer);

				// Get the query proxy and ref. kmer frequencies, these should be pretty much aligned.
				double query_proxy_kmer_freq = per_win_query_proxy_kmer_freq[var_i];

				int cur_ref_kmer_uniq_kmer_i = per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i]->at(cur_state);
				double backtracked_ref_hap_AF = (double)(per_win_per_unique_ref_kmer_cnts[var_i]->at(cur_ref_kmer_uniq_kmer_i)) / (2 * vecsize(ref_sample_ids));

				// Write.
				fprintf(stderr, "Variant %d: \n", var_i);

				fprintf(stderr, "Backtracked original k-mer @ var_i=%d [Query Proxy kmer freq: %.4f]\n", var_i, query_proxy_kmer_freq);
				dump_kmer(orig_query_kmer);
				fprintf(stderr, "Backtracked reference k-mer @ var_i=%d [Reference kmer freq: %.4f] (state: %d)\n", var_i, backtracked_ref_hap_AF, cur_state);
				dump_kmer(ref_kmer);

				bool found_kmer = false;
				//if (compare_kmers(ref_kmer, orig_query_kmer))
				if (ref_kmer->kmer[n_kmer_vicinity_vars] == orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					found_kmer = true;
				}
				else
				{
					n_errs++;
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////
				// Set the allele to this variant.
				char decoded_allele = ref_kmer->kmer[n_kmer_vicinity_vars];

				void** query_orig_reg_info = (void**)(query_original_geno_var_regs->at(var_i)->data);
				char* decoded_query_geno_sig = (char*)(query_orig_reg_info[1]);

				// Assign the decoded allele.
				int cur_geno = decoded_query_geno_sig[proxy_query_i_s];
				cur_geno = cur_geno | (decoded_allele << proxy_query_hap_i);
				decoded_query_geno_sig[proxy_query_i_s] = cur_geno;
				///////////////////////////////////////////////////////////////////////////////////////////////////

				// Check MAF AF errors.
				int max_AF_allele = ref_orig_geno_var_regs->at(var_i)->score;

				bool found_max_AF_allele = false;
				if (max_AF_allele != orig_query_kmer->kmer[n_kmer_vicinity_vars])
				{
					n_max_AF_errs++;
				}
				else
				{
					found_max_AF_allele = true;
				}

				// Write the per subject info.
				fprintf(f_per_sample_per_var_dec_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n",
						query_sample_ids->at(proxy_query_i_s), proxy_query_hap_i,
						query_proxized_geno_var_regs->at(var_i)->chrom, query_proxized_geno_var_regs->at(var_i)->start,
						per_win_n_ref_kmers[var_i],
						cur_state,
						per_win_query_proxy_kmer_freq[var_i], backtracked_ref_hap_AF, (int)found_kmer, (int)found_max_AF_allele);
				fflush(f_per_sample_per_var_dec_stats);

				// Get the new state for the previous variant.
				cur_state = ML_prev_state_per_hap[var_i][proxy_query_hap_i][cur_state];

				// Note that we have effectively moved one variant back.
				if (cur_state == -1)
				{
					if (var_i >= var_start_i)
					{
						fprintf(stderr, "Sanity check failed: Backtracking stopped before starting index.\n");
						exit(1);
					}

					fprintf(stderr, "Breaking @ variant %d.\n", var_i);
					break;
				}
				else
				{
					// Make sure that this state is valid.
					if (var_i > var_start_i &&
						cur_state > vecsize(per_win_all_ref_kmers[var_i - 1]))
					{
						fprintf(stderr, "Sanity check failed: Invalid backtracked state index @ var_i=%d.\n", var_i);
						exit(1);
					}
				}
			} // var_i loop.

			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fprintf(f_dec_stats, "%d\t%d\t%d\t%d\t%.5f\t%.5f\n", proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_max_proxy_score_per_full_hap, least_scoring_kmer_score);
			fflush(f_dec_stats);

			// Free memory, this is necessary for long runs.
			for (int i_var = var_start_i; i_var < var_end_i; i_var++)
			{
				// This is the max log-probability score for each unique kmer in the window.
				//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
				//delete[] per_win_per_ref_kmer_log_cumul_probs[i_var];

				// # of unique ref. kmers in the window.
				//int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];

				// List of unique ref. kmers per window.
				//vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
				delete_kmers(per_win_unique_ref_kmers[i_var]);
				delete_kmers(per_win_all_ref_kmers[i_var]);
				delete(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[i_var]);

				// The best scoring unique kmer index in the previous window for each current window.
				//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
				//delete[] per_win_backtracking_kmer_i[i_var];

				// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
				//double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];

				// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
				//vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
				delete(per_win_per_unique_ref_kmer_cnts[i_var]);

				delete(per_win_query_per_ref_hap_query_emission_prob[i_var]);
			} // i_var loop.

			delete[] per_win_unique_ref_kmers;
			delete[] per_win_per_unique_ref_kmer_cnts;
			delete[] per_win_query_proxy_kmer_freq;
			delete[] per_win_unique_ref_kmer_i_per_orig_ref_kmer_i;
			delete[] per_win_n_ref_kmers;
			delete[] per_win_query_per_ref_hap_query_emission_prob;
			delete[] per_win_all_ref_kmers;
		} // query_hap_i loop.

		// Save the genotypes for the decoded alleles on this subject.
		vector<t_annot_region*>* cur_sample_decoded_geno_regs = new vector<t_annot_region*>();
		for (int i_var = var_start_i; i_var < var_end_i; i_var++)
		{
			void** cur_var_reg_info = (void** )(query_original_geno_var_regs->at(i_var)->data);
			char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

			void** decoded_var_reg_info = new void* [3];
			char* decoded_geno_sig = new char[5];
			decoded_geno_sig[0] = var_reg_decoded_geno_sig[proxy_query_i_s];
			decoded_var_reg_info[0] = decoded_geno_sig;

			t_annot_region* decoded_geno_reg = duplicate_region(query_original_geno_var_regs->at(i_var));
			decoded_geno_reg->data = decoded_var_reg_info;
			cur_sample_decoded_geno_regs->push_back(decoded_geno_reg);
		} // i_var loop.

		// Set the sample id's and file name.
		vector<char*>* decoded_query_sample_ids = new vector<char*>();
		decoded_query_sample_ids->push_back(query_sample_ids->at(proxy_query_i_s));
		char decoded_sample_fp[1000];
		sprintf(decoded_sample_fp, "decoded_sample_list_%d.txt", proxy_query_i_s);
		FILE* f_decoded_sample_list = open_f(decoded_sample_fp, "w");
		fprintf(f_decoded_sample_list, "%s\n", query_sample_ids->at(proxy_query_i_s));
		close_f(f_decoded_sample_list, decoded_sample_fp);

		fprintf(stderr, "Saving decoded genotypes for the subject..");
		char decoded_query_geno_fp[1000];
		sprintf(decoded_query_geno_fp, "decoded_sample_%d.matbed.gz", proxy_query_i_s);
		binarize_variant_genotype_signal_regions(cur_sample_decoded_geno_regs, NULL, decoded_query_sample_ids, decoded_query_geno_fp);

		// Free memory for the decoded genotype regions on this subject.
		delete_annot_regions(cur_sample_decoded_geno_regs);
	} // proxy_query_i_s loop.

	// Save the decoded genotypes: Replace the decoded genotypes with actual genotypes.
	vector<t_annot_region*>* decoded_geno_regs = new vector<t_annot_region*>();
	for (int i_var = var_start_i; i_var < var_end_i; i_var++)
	{
		void** cur_var_reg_info = (void**)(query_original_geno_var_regs->at(i_var)->data);
		//char* var_reg_decoded_geno_sig = (char*)(cur_var_reg_info[1]);

		// Replace the original variant signal to decoded genotype signal.
		cur_var_reg_info[0] = cur_var_reg_info[1];

		t_annot_region* decoded_geno_reg = duplicate_region(query_original_geno_var_regs->at(i_var));
		decoded_geno_reg->data = cur_var_reg_info;
		decoded_geno_regs->push_back(decoded_geno_reg);
	} // i_var loop.

	fprintf(stderr, "Saving all decoded genotypes for the subject..");
	char decoded_query_geno_fp[1000];
	sprintf(decoded_query_geno_fp, "decoded_all_sample.matbed.gz");
	binarize_variant_genotype_signal_regions(decoded_geno_regs, NULL, query_sample_ids, decoded_query_geno_fp);

	fprintf(stderr, "Done!\n");
} // decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_Viterbi function.


void decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack(char* query_original_haplocoded_geno_fp, char* query_proxized_haplocoded_geno_fp, char* query_sample_ids_fp,
	char* ref_original_haplocoded_geno_fp, char* ref_sample_ids_fp,
	char* recombination_rate_dir,
	int var_start_i, int var_end_i,
	int n_kmer_vicinity_vars,
	double kmer_concordance_log_weight,
	double N_e,
	int n_query_subjects_2_decode,
	char* op_prefix)
{
	fprintf(stderr, "Decoding proxized alleles of query site using histogram matching HMM with fore-back:\n\
Query Original: %s\n\
Query Proxy: %s\n\
Query Sample list: %s\n\
Ref. Original: %s\n\
Ref. Sample list: %s\n\
# vic.: %d\n\
Variant range: [%d-%d]\n\
kmer_concordance_log_weight: %.4f\n\
N_e: %.4f\n\
# Query Subjects 2 Decode: %d", query_original_haplocoded_geno_fp, query_proxized_haplocoded_geno_fp, query_sample_ids_fp,
ref_original_haplocoded_geno_fp, ref_sample_ids_fp, n_kmer_vicinity_vars,
var_start_i, var_end_i, kmer_concordance_log_weight, N_e, n_query_subjects_2_decode);

	vector<t_annot_region*>* query_proxized_geno_var_regs = load_variant_signal_regions_wrapper(query_proxized_haplocoded_geno_fp, query_sample_ids_fp);
	vector<t_annot_region*>* query_original_geno_var_regs = load_variant_signal_regions_wrapper(query_original_haplocoded_geno_fp, query_sample_ids_fp);
	vector<char*>* query_sample_ids = buffer_file(query_sample_ids_fp);
	vector<t_annot_region*>* ref_orig_geno_var_regs = load_variant_signal_regions_wrapper(ref_original_haplocoded_geno_fp, ref_sample_ids_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_fp);

	sort(query_proxized_geno_var_regs->begin(), query_proxized_geno_var_regs->end(), sort_regions);
	sort(query_original_geno_var_regs->begin(), query_original_geno_var_regs->end(), sort_regions);
	sort(ref_orig_geno_var_regs->begin(), ref_orig_geno_var_regs->end(), sort_regions);

	int max_query_proxy_geno = get_max_genotype_value(query_proxized_geno_var_regs, query_sample_ids);
	int max_query_original_geno = get_max_genotype_value(query_original_geno_var_regs, query_sample_ids);
	int max_ref_original_geno = get_max_genotype_value(ref_orig_geno_var_regs, ref_sample_ids);

	// Input checks.
	if (var_start_i <= n_kmer_vicinity_vars ||
		var_end_i >= (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars - 1))
	{
		fprintf(stderr, "Make sure variant start-end indices are within the kmer vicinity of variant ends: [%d-%d] vs [%d-%d]\n", 
				var_start_i, var_end_i, 
				n_kmer_vicinity_vars, (vecsize(query_proxized_geno_var_regs) - n_kmer_vicinity_vars));
		exit(1);
	}

	if (max_query_proxy_geno != 3 ||
		max_query_original_geno != 3 ||
		max_ref_original_geno != 3)
	{
		fprintf(stderr, "One of the panels is not haplocoded..\n");
		exit(1);
	}

	if (query_proxized_geno_var_regs->size() != query_original_geno_var_regs->size() ||
		query_proxized_geno_var_regs->size() != ref_orig_geno_var_regs->size())
	{
		fprintf(stderr, "Sanity check failed: Variant counts are not the same: %d, %d, %d\n",
			vecsize(query_proxized_geno_var_regs), vecsize(query_original_geno_var_regs),
			vecsize(ref_orig_geno_var_regs));

		exit(1);
	}

	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		if (ref_orig_geno_var_regs->at(i_var)->start != query_proxized_geno_var_regs->at(i_var)->start ||
			ref_orig_geno_var_regs->at(i_var)->start != query_original_geno_var_regs->at(i_var)->start)
		{
			fprintf(stderr, "Sanity check failed: The variant coordinates are not matching among panels.\n");
			exit(1);
		}
	} // i_var loop.

	//////////////////////////////////////////////////////////////////////////////////////////
	// Check & fix variant range.
	if (var_start_i < 1 ||
		var_start_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_start_i = 1;
	}

	if (var_end_i < 0 ||
		var_end_i >= vecsize(query_proxized_geno_var_regs))
	{
		var_end_i = vecsize(query_proxized_geno_var_regs);
	}

	if (var_end_i < (var_start_i + 100))
	{
		fprintf(stderr, "Invalid range: %d-%d\n", var_start_i, var_end_i);
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////

	fprintf(stderr, "Loaded %d (%d), %d (%d) variants for query-proxy, ref-original panels.\n",
		(int)query_proxized_geno_var_regs->size(), (int)query_sample_ids->size(),
		(int)ref_orig_geno_var_regs->size(), (int)ref_sample_ids->size());

	fprintf(stderr, "Extracting per subject haplotypes.\n");
	vector<char**>* per_ref_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(ref_orig_geno_var_regs, ref_sample_ids);
	vector<char**>* per_query_subj_orig_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_original_geno_var_regs, query_sample_ids);
	vector<char**>* per_query_subj_proxy_haplotypes = get_per_subject_haplotypes_per_haplocoded_var_regs(query_proxized_geno_var_regs, query_sample_ids);

	//////////////////////////////////////////////////////////////////////////////////////////
	// Assign the maximum frequency alleles to all reference variants.
	fprintf(stderr, "Assigning maximum frequency alleles to reference regions.\n");
	char max_AFs_bed_fp[1000];
	sprintf(max_AFs_bed_fp, "%s_max_AF_alleles.bed", op_prefix);
	FILE* f_max_AF_bed = open_f(max_AFs_bed_fp, "w");
	for (int i_var = 0; i_var < vecsize(ref_orig_geno_var_regs); i_var++)
	{
		void** cur_ref_var_info = (void**)(ref_orig_geno_var_regs->at(i_var)->data);
		char* cur_ref_geno_sig = (char*)(cur_ref_var_info[0]);

		double total_AA = 0;
		for (int i_s = 0; i_s < vecsize(ref_sample_ids); i_s++)
		{
			total_AA += get_genotype_per_haplocoded_genotype(cur_ref_geno_sig[i_s]);
		} // i_s loop.

		double AAF = total_AA / (2 * vecsize(ref_sample_ids));

		if (AAF > 0.5)
		{
			ref_orig_geno_var_regs->at(i_var)->score = 1;
		}
		else
		{
			ref_orig_geno_var_regs->at(i_var)->score = 0;
		}

		fprintf(f_max_AF_bed, "%s\t%d\t%d\t%s\t%d\t+\n",
				ref_orig_geno_var_regs->at(i_var)->chrom,
				translate_coord(ref_orig_geno_var_regs->at(i_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(ref_orig_geno_var_regs->at(i_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				ref_orig_geno_var_regs->at(i_var)->name,
				ref_orig_geno_var_regs->at(i_var)->score);
	} // i_var loop.
	close_f(f_max_AF_bed, NULL);
	//////////////////////////////////////////////////////////////////////////////////////////
	
	char cur_chr_recombination_rate_fp[1000];
	sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, ref_orig_geno_var_regs->at(0)->chrom);
	vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
	if (cur_chrom_recomb_regs == NULL)
	{
		fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
		exit(1);
	}

	// Add one element to the beginning of the regions for the current chromosome.
	// Assign the recomb rates.
	fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", vecsize(cur_chrom_recomb_regs));
	for (int i_reg = 0; i_reg < vecsize(ref_orig_geno_var_regs); i_reg++)
	{
		double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(ref_orig_geno_var_regs->at(i_reg), cur_chrom_recomb_regs);
		ref_orig_geno_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
	} // i_reg loop.

	int n_states = 2 * vecsize(ref_sample_ids);
	int n_ref_haplotypes = n_states;

	// We track a list of candidate kmers at each position.
	char dec_summary_stats_fp[1000];
	sprintf(dec_summary_stats_fp, "%s_per_haplotype_decoding_summary.txt", op_prefix);
	FILE* f_dec_stats = open_f(dec_summary_stats_fp, "w");

	char per_sample_per_var_dec_stats_fp[1000];
	sprintf(per_sample_per_var_dec_stats_fp, "%s_per_per_sample_decoding_stats.txt", op_prefix);
	FILE* f_per_sample_per_var_dec_stats = open_f(per_sample_per_var_dec_stats_fp, "w");

	// Allocate arrays for storing scores and backtracking info.
	int n_vars = vecsize(ref_orig_geno_var_regs);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Foreback arrays are allocated only once and reused. Make sure inits are done correctly.
	// Allocate and initialize the forward/backward arrays.
	double*** back_scores_per_hap = new double** [n_vars + 2];
	double*** fore_scores_per_hap = new double** [n_vars + 2];

	for (int var_i = 0; var_i <= n_vars + 1; var_i++)
	{
		back_scores_per_hap[var_i] = new double* [2];
		fore_scores_per_hap[var_i] = new double* [2];

		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			fore_scores_per_hap[var_i][proxy_query_hap_i] = new double[n_states + 2];
			back_scores_per_hap[var_i][proxy_query_hap_i] = new double[n_states + 2];
		} // proxy_query_hap_i loop.
	} // i loop.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Start looping over all individuals.
	if (n_query_subjects_2_decode > vecsize(query_sample_ids) ||
		n_query_subjects_2_decode <= 0)
	{
		n_query_subjects_2_decode = vecsize(query_sample_ids);
	}

	fprintf(stderr, "Decoding %d/%d subjects.\n", n_query_subjects_2_decode, vecsize(query_sample_ids));

	for (int proxy_query_i_s = 0; proxy_query_i_s < n_query_subjects_2_decode; proxy_query_i_s++)
	{
		// Use all variants on the chromosome.
		vector<t_annot_region*>* cur_win_var_regs = ref_orig_geno_var_regs;

		// Start recursing over the variants.
		for (int proxy_query_hap_i = 0; proxy_query_hap_i < 2; proxy_query_hap_i++)
		{
			// Initialize all scores and backtracking states.
			for (int var_i = 0; var_i <= n_vars + 1; var_i++)
			{
				// INitialize all the scores for all variants.
				for (int cur_state = 0; cur_state < n_ref_haplotypes; cur_state++)
				{
					fore_scores_per_hap[var_i][proxy_query_hap_i][cur_state] = xlog(0.0);
					back_scores_per_hap[var_i][proxy_query_hap_i][cur_state] = xlog(0.0);
				} // cur_state loop.
			} // i loop.

			// Initialize the boundary conditions for fore and back variables. These are set only once and used at the boundaries.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				fore_scores_per_hap[var_start_i - 1][proxy_query_hap_i][state_i] = xlog((double)1.0 / n_ref_haplotypes);
				back_scores_per_hap[var_end_i][proxy_query_hap_i][state_i] = xlog((double)1.0 / n_ref_haplotypes);
			} // state_i loop.

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Following are the kmer statistics and information in each window.
			// These are allocated here and initialized while calculating fore array and reused in back array.
			// # of unique ref. kmers in the window.
			int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];
			// List of unique ref. kmers per window.
			vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
			// This is the list of all unsorted ref kmers in the current window, we need this for keeping track of haplotype sequence overlaps.
			vector<t_kmer*>** per_win_all_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
			// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window for the current query's kmer on the haplotype.
			double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];
			// The number of ref unique kmer counts per unique kmer at each window; this is used for calculating kmer frequencies in the window.
			vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
			// This is the list of indices for each original (unsorted) ref kmer to map to unique kmers.
			// This is necessary for easy mapping of the frequencies for unsorted reference kmers. This is needed since each state is one haplotype (not unique).
			vector<int>** per_win_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>*[query_proxized_geno_var_regs->size()];
			// This is the query hapfreq emission log probability for each unsorted reference kmer.
			vector<double>** per_win_query_per_ref_hap_query_emission_prob = new vector<double>*[query_proxized_geno_var_regs->size()];

			// Process all variants: Note that this loop sets up the kmer statistics for boundaries as well.
			// The boundaries are only processed to setup the kmer statistics.
			for (int var_i = (var_start_i-1); var_i <= var_end_i; var_i++)
			{
				if (var_i % 10 == 0)
				{
					fprintf(stderr, "Forward: sample_i: %d/%d: var_i: %d         \r", proxy_query_i_s, vecsize(query_sample_ids), var_i);
				}

				// Extract the k-mers on this window.
				// Set the window starts and ends.
				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////////////////////////////////
				// Extract the current window's kmers within reference panel.
				vector<t_kmer*>* cur_win_all_ref_kmers = extract_kmers_per_haplotype(per_ref_subj_orig_haplotypes, i_win_start, i_win_end);

				// Extract the current window's kmer for the query panel.
				vector<t_kmer*>* cur_win_all_query_proxy_kmers = extract_kmers_per_haplotype(per_query_subj_proxy_haplotypes, i_win_start, i_win_end);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "Extracted %d reference kmers @ [%d, %d]\n", (int)cur_win_all_ref_kmers->size(), i_win_start, i_win_end);
					//dump_kmers(cur_win_all_ref_kmers);
				}
				
				/////////////////////////////////////////////////////////////////////////////////
				// Get the reference haplotype frequencies.
				// Get the window reference k-mers
				vector<t_kmer*>* cur_win_unique_ref_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_ref_kmer_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_ref_kmers, cur_win_unique_ref_kmers, cur_win_unique_ref_kmer_cnts);

				// This is necessary to map the reference k-mers to their probabilities.
				vector<int>* cur_unique_ref_kmer_i_per_orig_ref_kmer_i = new vector<int>();
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_all_ref_kmers); i_kmer++)
				{
					for (int i_unique_kmer = 0; i_unique_kmer < vecsize(cur_win_unique_ref_kmers); i_unique_kmer++)
					{
						if (compare_kmers(cur_win_all_ref_kmers->at(i_kmer), cur_win_unique_ref_kmers->at(i_unique_kmer)))
						{
							cur_unique_ref_kmer_i_per_orig_ref_kmer_i->push_back(i_unique_kmer);
							break;
						}
					} // i_unique_kmer
				} // i_kmer loop.

				int n_unique_ref_haps = cur_win_unique_ref_kmer_cnts->size();

				// These keep state of ref kmer scores, backtracks, etc.
				per_win_per_unique_ref_kmer_cnts[var_i] = cur_win_unique_ref_kmer_cnts;
				per_win_n_ref_kmers[var_i] = n_unique_ref_haps;
				per_win_unique_ref_kmers[var_i] = cur_win_unique_ref_kmers;
				per_win_all_ref_kmers[var_i] = cur_win_all_ref_kmers;
				per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[var_i] = cur_unique_ref_kmer_i_per_orig_ref_kmer_i;

				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Get the query proxy k-mers.
				vector<t_kmer*>* cur_win_unique_query_proxy_kmers = new vector<t_kmer*>();
				vector<int>* cur_win_unique_query_proxy_cnts = new vector<int>();
				get_unique_kmers_w_counts(cur_win_all_query_proxy_kmers, cur_win_unique_query_proxy_kmers, cur_win_unique_query_proxy_cnts);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "%d unique reference kmers @ [%d, %d]\n", (int)cur_win_unique_ref_kmers->size(), i_win_start, i_win_end);
					fprintf(stderr, "%d unique proxy query kmers @ [%d, %d]\n", (int)cur_win_unique_query_proxy_kmers->size(), i_win_start, i_win_end);
				}

				// Get the current proxy haplotype of the query subject and extract the kmer.
				char* cur_query_subj_proxy_hap = per_query_subj_proxy_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i];
				t_kmer* cur_query_subj_proxy_kmer = copy_kmer(cur_win_all_query_proxy_kmers->at(2 * proxy_query_i_s + proxy_query_hap_i));

				// Note that the proxy query kmer does not mean anything to the adversary, she just needs the frequency of it within the proxy panel.
				delete_kmers(cur_win_all_query_proxy_kmers);

				// This is a sanity check for the extracted proxy kmer.
				for (int check_var_i = i_win_start; check_var_i <= i_win_end; check_var_i++)
				{
					if (cur_query_subj_proxy_hap[check_var_i] != cur_query_subj_proxy_kmer->kmer[check_var_i - i_win_start])
					{
						fprintf(stderr, "Sanity check failed while getting the kmer for subject.\n");

						exit(1);
					}
				} // check_var_i loop.

				// Search for the query proxy kmer in unique kmers for the window.
				double cur_query_subj_proxy_kmer_cnt = 0.0;
				double total_query_subj_proxy_kmer_cnt = 0.0;
				for (int i_kmer = 0; i_kmer < vecsize(cur_win_unique_query_proxy_kmers); i_kmer++)
				{
					if (compare_kmers(cur_win_unique_query_proxy_kmers->at(i_kmer), cur_query_subj_proxy_kmer))
					{
						cur_query_subj_proxy_kmer_cnt = cur_win_unique_query_proxy_cnts->at(i_kmer);
					} // current unique proxy kmer check for the current subject's query kmer.

					total_query_subj_proxy_kmer_cnt += cur_win_unique_query_proxy_cnts->at(i_kmer);
				} // i_kmer loop.
				delete_kmer(cur_query_subj_proxy_kmer);
				delete_kmers(cur_win_unique_query_proxy_kmers);
				delete(cur_win_unique_query_proxy_cnts);

				// This is the frequency of proxy haplotype at this window that we will match to from the reference panel.
				double cur_query_subj_proxy_kmer_freq = cur_query_subj_proxy_kmer_cnt / (2 * vecsize(query_sample_ids));
				if (cur_query_subj_proxy_kmer_freq == 0 ||
					cur_query_subj_proxy_kmer_freq > 1.0)
				{
					fprintf(stderr, "Sanity Check failed: Found 0 frequency for the query proxy haplotype.\n");
					exit(1);
				}

				// Sanity check on the total count of the query kmer count.
				if (total_query_subj_proxy_kmer_cnt != (2 * vecsize(query_sample_ids)))
				{
					fprintf(stderr, "Sanity Check failed: Total query ref kmer counts does not add up to 1: %.3f/%d",
						total_query_subj_proxy_kmer_cnt, 2 * vecsize(query_sample_ids));

					exit(1);
				}

				// Save the frequency of query's kmer within its own cohort.
				per_win_query_proxy_kmer_freq[var_i] = cur_query_subj_proxy_kmer_freq;
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////////////////////////////////
				// Calculate the normalized proxy query kmer emission probability for all of the reference kmers at this window.
				double query_emit_log_normalizing_factor = xlog(0.0);
				double total_ref_unique_kmer_prob = 0; // This is for sanity check.

				// This stores the list of emission probabilities for this proxy query kmer at the respective ref haplotype.
				vector<double>* per_state_query_emission_prob = new vector<double>();
				per_win_query_per_ref_hap_query_emission_prob[var_i] = per_state_query_emission_prob;

				for (int cur_ref_kmer_i = 0; cur_ref_kmer_i < n_ref_haplotypes; cur_ref_kmer_i++)
				{
					// First get the unique kmer index for this kmer.
					int cur_kmer_unique_kmer_i = cur_unique_ref_kmer_i_per_orig_ref_kmer_i->at(cur_ref_kmer_i);
					double cur_ref_kmer_freq = (double)cur_win_unique_ref_kmer_cnts->at(cur_kmer_unique_kmer_i) / (2 * vecsize(ref_sample_ids));
					if (cur_ref_kmer_freq == 0)
					{
						fprintf(stderr, "Sanity check failed: Unique reference k-mer frequency is 0.\n");
						exit(1);
					}

					if (!compare_kmers(cur_win_all_ref_kmers->at(cur_ref_kmer_i), cur_win_unique_ref_kmers->at(cur_kmer_unique_kmer_i)))
					{
						fprintf(stderr, "Sanity check failed: All-2-unique mapped k-mer does not match unique kmer.\n");
						exit(1);
					}

					total_ref_unique_kmer_prob += cur_ref_kmer_freq;

					double abs_delta_hapfreq = fabs(cur_ref_kmer_freq - cur_query_subj_proxy_kmer_freq);
					double proxy_kmer_emit_log_prob_per_cur_kmer = -1 * abs_delta_hapfreq;

					// Add this emission probability to the list of emission probabilities.
					per_state_query_emission_prob->push_back(proxy_kmer_emit_log_prob_per_cur_kmer);

					// Update the normalizing factor.
					query_emit_log_normalizing_factor = xlog_sum(query_emit_log_normalizing_factor, proxy_kmer_emit_log_prob_per_cur_kmer);
				} // cur_kmer loop.

				// Normalize all of the emission probabilities.
				for (int cur_kmer = 0; cur_kmer < n_ref_haplotypes; cur_kmer++)
				{
					per_state_query_emission_prob->at(cur_kmer) = xlog_div(per_state_query_emission_prob->at(cur_kmer), query_emit_log_normalizing_factor);
				} // cur_kmer loop.

				// Note that this check does not hold any more since we are working over all of the haplotypes.
				//// Make sure that total kmer counts equal to 1.
				//if (fabs(total_ref_unique_kmer_prob - 1.0) > 0.001)
				//{
				//	fprintf(stderr, "Sanity check failed on total unique ref prob: %.4f\n", total_ref_unique_kmer_prob);
				//	exit(1);
				//}
				/////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////////////
				// Pre-compute the transition probabilities: This is calculated once before all states since it only depends on genetic distance.
				double prev_var_cM = cur_win_var_regs->at(var_i -1)->dbl_score;

				double cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;

				double r_m = fabs(cur_var_cM - prev_var_cM);
				double rho_m = 4 * N_e * r_m;

				double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

				double other_prob = tau_m / n_ref_haplotypes;
				double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

				//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
				/////////////////////////////////////////////////////

				if (n_ref_haplotypes != vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i));
					exit(1);
				}

				if (n_ref_haplotypes != (2 * vecsize(ref_sample_ids)))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i));
					exit(1);
				}

				////////////////////////////////////////////////////////////////////////////////////////////////
				// This is an important check, we do up to here at boundaries to setup the kmer frequency statistics.
				// Do not calculate the arrays for boundary points.
				if (var_i <= var_start_i - 1 ||
					var_i >= var_end_i)
				{
					continue;
				}
				////////////////////////////////////////////////////////////////////////////////////////////////

				// Start calculating fore array.
				double total_normalized_kmer_emission_log_prob = xlog(0.0);
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					// Initialize the foreward score at this variant and state (reference haplotype)
					fore_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog(0.0);

					// This is a sanity check value on the emission probability of query kmer by ref kmers.
					total_normalized_kmer_emission_log_prob = xlog_sum(total_normalized_kmer_emission_log_prob, per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i));

					// We have the emission of query proxy from current state, now calculate the transitions.					
					for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
					{
						// Calculate the kmer concordance penalty.
						double kmer_concordance_log_prob = xlog(1.0);

						/////////////////
						// Get the kmer concordance penalty score for the previous state.
						double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i));

						kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

						// Sanity check on concordant kmers; if they are concordant, the distance must be 0.
						if (are_kmers_shift_compatible(per_win_all_ref_kmers[var_i - 1]->at(prev_state_i), per_win_all_ref_kmers[var_i]->at(state_i)) &&
							shifted_kmer_dist > 0)
						{
							fprintf(stderr, "Sanity check failed: Compatability failure.\n");
							exit(1);
						}
						/////////////////

						/////////////////
						// Set the emission probability at this window/state and add the kmer concordance penalty.
						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_mul(per_win_query_per_ref_hap_query_emission_prob[var_i]->at(state_i), kmer_concordance_log_prob);

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
						// Set the transition probabilities.
						double trans_emit_prob = xlog(1.0);
						if (state_i == prev_state_i)
						{
							trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, xlog(self_prob));
						}
						else
						{
							trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, xlog(other_prob));
						}

						/*fprintf(stderr, "%d: %d->%d: trans prob: %.4fx %.4f / %.4f ;; %.4f\n", var_i,
							state_i, prev_state_i,
							exp(norm_proxy_kmer_emit_log_prob_per_cur_kmer), self_prob, other_prob,
							trans_prob);*/

							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Always emits what is on the haplotype.

						// Update fore array score.
						fore_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog_sum(fore_scores_per_hap[var_i][proxy_query_hap_i][state_i],
																xlog_mul(trans_emit_prob, fore_scores_per_hap[var_i - 1][proxy_query_hap_i][prev_state_i]));

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "fore[%d][%d]: %.5f: P_trans_emit(%d->%d): %.5f\n",
								var_i, state_i, fore_scores_per_hap[var_i][proxy_query_hap_i][state_i],
								prev_state_i, state_i, trans_emit_prob);
						}
					} // prev_state_i loop.

					//fprintf(stderr, "Var %d: %d non-compatible kmers with previous haps.\n", var_i, state_i, n_non_compatible_kmers);
				} // state_i loop.

				// Following checks to make sure that we used a valid probability distribution for the emission of kmer's based on
				// haplotype frequencies.
				if (fabs(exp(total_normalized_kmer_emission_log_prob) - 1.0) > 0.001)
				{
					fprintf(stderr, "Sanity check failed: Total kmer emission freq is not valud: %.4f\n", total_normalized_kmer_emission_log_prob);
					exit(1);
				}
			} // foreward var_i loop.

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Following is the backward array calculation.
			// At this point, we setup the necessary k-mer statistics for all windows (i.e. variant positions) including emission probabilities, kmers, etc.
			// Calculate the backward array.
			for (int var_i = var_end_i-1; var_i >= var_start_i; var_i--)
			{
				if (var_i % 10 == 0)
				{
					fprintf(stderr, "Backward: sample_i: %d/%d: var_i: %d         \r", proxy_query_i_s, vecsize(query_sample_ids), var_i);
				}

				// Extract the k-mers on this window.
				// Set the window starts and ends.
				int i_win_start = MAX(0, var_i - n_kmer_vicinity_vars);
				int i_win_end = MIN(var_i + n_kmer_vicinity_vars, vecsize(query_proxized_geno_var_regs) - 1);

				if (__DUMP_KMER_TRACKING_MSGS__)
				{
					fprintf(stderr, "================================\n");
					fprintf(stderr, "Extracting the original reference kmers @ [%d, %d]\n", i_win_start, i_win_end);
				}

				/////////////////////////////////////////////////////
				// Pre-compute the transition probability at this variant.
				double cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;

				// Note that for backward variable we do a transition from current to next.
				double next_var_cM = cur_win_var_regs->at(var_i+1)->dbl_score;

				double r_m = fabs(cur_var_cM - next_var_cM);
				double rho_m = 4 * N_e * r_m;

				double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

				double other_prob = tau_m / n_ref_haplotypes;
				double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

				//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, next_var_cM, cur_var_cM);
				/////////////////////////////////////////////////////

				//if (n_ref_haplotypes != vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i))
				//{
				//	fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n", n_ref_haplotypes, vecsize(cur_unique_ref_kmer_i_per_orig_ref_kmer_i));
				//	exit(1);
				//}

				if (n_ref_haplotypes != (2 * vecsize(ref_sample_ids)))
				{
					fprintf(stderr, "Sanity check failed on ref haplotype counts: %d/%d\n",
						n_ref_haplotypes, (2 * vecsize(ref_sample_ids)));
					exit(1);
				}

				// Loop over all the states.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					// Recurse over the previous states.
					back_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog(0.0);

					for (int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
					{
						/////////////////
						// Set the kmer concordance penalty.
						double kmer_concordance_log_prob = xlog(1.0);

						// The check is necessary since back array is calculated with emission of the next window's query proxy kmer.
						// Make sure we go from left to right in concordance assignment.
						double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_i]->at(state_i), per_win_all_ref_kmers[var_i + 1]->at(next_state_i));

						kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_i]->at(state_i)->kmer_length);

						//if (!are_kmers_shift_compatible(cur_win_all_ref_kmers->at(state_i), per_win_all_ref_kmers[var_i + 1]->at(next_state_i)))
						//{
						//	//fprintf(stderr, "%d -> %d is not k-mer compatible\n", prev_state_i, state_i);
						//	n_non_compatible_kmers++;
						//	continue;
						//}
						/////////////////

						/////////////////
						// Set the normalized emission prob: This is a little tricky for backward array because we emit in the next window/state
						// We always have the kmer concordance. If the next window is the boundary, do not use emission, use only kmer concordance.
						double norm_proxy_kmer_emit_log_prob_per_cur_kmer = kmer_concordance_log_prob;
						/////////////////

						/////////////////
						// Get the emission probability by next reference kmer but make sure to handle the boundary condition, we need to handle this since
						// we only want to look at the kmer concordance at the boundaries, not emissions.
						if (var_i < (var_end_i - 1))
						{
							// Use the same emission probability that we used for fore array.
							// In the backward array, we are emitting the next kmer at the next state.
							norm_proxy_kmer_emit_log_prob_per_cur_kmer = xlog_mul(per_win_query_per_ref_hap_query_emission_prob[var_i + 1]->at(next_state_i),
																					kmer_concordance_log_prob);

							if (norm_proxy_kmer_emit_log_prob_per_cur_kmer == xlog(0.0))
							{
								fprintf(stderr, "Sanity check failed; kmer emission at next window is 0.\n");
								exit(1);
							}
						}
						/////////////////

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
						// Set the transition probabilities.
						double trans_emit_prob = xlog(1.0);
						if (state_i == next_state_i)
						{
							trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, xlog(self_prob));
						}
						else
						{
							trans_emit_prob = xlog_mul(norm_proxy_kmer_emit_log_prob_per_cur_kmer, xlog(other_prob));
						}

						/*fprintf(stderr, "%d: %d->%d: trans prob: %.4fx %.4f / %.4f ;; %.4f\n", var_i,
							state_i, prev_state_i,
							exp(norm_proxy_kmer_emit_log_prob_per_cur_kmer), self_prob, other_prob,
							trans_prob);*/

						// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

						back_scores_per_hap[var_i][proxy_query_hap_i][state_i] = xlog_sum(back_scores_per_hap[var_i][proxy_query_hap_i][state_i], 
																						xlog_mul(trans_emit_prob, back_scores_per_hap[var_i + 1][proxy_query_hap_i][next_state_i]));

						if (__DUMP_KMER_TRACKING_MSGS__)
						{
							fprintf(stderr, "back[%d][%d]: %.5f: P_trans_emit(%d->%d): %.5f\n",
								var_i, state_i, back_scores_per_hap[var_i][proxy_query_hap_i][state_i],
								state_i, next_state_i, trans_emit_prob);
						}
					} // next_state_i loop.

					//fprintf(stderr, "Var %d: %d non-compatible kmers with previous haps.\n", var_i, state_i, n_non_compatible_kmers);
				} // state_i loop.
			} // backward variable var_i loop.

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Fore-back is finished, assign the probabilities.
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// Compute the total log forward and backward probabilities.
			double total_foreback_score_per_fore = xlog(0.0);
			double total_foreback_score_per_back = xlog(0.0);

			// Calculate total probabilities for sanity checkking.
			// This is a little tricky since we look at kmer concordance at the boundaries.
			// The boundaries also include the hap-hap transitions from genetic map, which means we need to weight the kmer concordance with hap-hap transition.
			for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
			{
				// Fore array: Look at the last boundary window (index: var_end_i) and loop over all possible states at the boundary.
				// Make sure to include the hap-hap transitions and the kmer concordance from var_end_i-1 to var_end.
				for (int last_state_i = 0; last_state_i < n_ref_haplotypes; last_state_i++)
				{
					/////////////////////////
					// Pre-compute the transition probability at this variant.
					double cur_var_cM = cur_win_var_regs->at(var_end_i-1)->dbl_score;

					// Note that for backward variable we do a transition from current to next.
					double next_var_cM = cur_win_var_regs->at(var_end_i)->dbl_score;

					double r_m = fabs(cur_var_cM - next_var_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

					double state_trans_log_prob = xlog(other_prob);
					if (state_i == last_state_i)
					{
						state_trans_log_prob = xlog(self_prob);
					}

					/////////////////////////
					// Get k-mer transition score.
					double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_end_i-1]->at(state_i), per_win_all_ref_kmers[var_end_i]->at(last_state_i));
					double kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_end_i]->at(state_i)->kmer_length);

					/////////////////////////
					double state_kmer_transition_prob = xlog_mul(kmer_concordance_log_prob, state_trans_log_prob);

					/////////////////////////
					total_foreback_score_per_fore = xlog_sum(total_foreback_score_per_fore, 
															xlog_mul(state_kmer_transition_prob, fore_scores_per_hap[var_end_i - 1][proxy_query_hap_i][state_i]));
				} // last_state_i loop.
				
				// Back array: Loop over all the haplotypes at the first boundary window (index: var_start_i-1).
				// Account for (1) hap-hap transitions, (2) kmer concordance, (3) Emission of the first window.
				for (int first_state_i = 0; first_state_i < n_ref_haplotypes; first_state_i++)
				{
					/////////////////////////
					// Pre-compute the transition probability at this variant.
					double cur_var_cM = cur_win_var_regs->at(var_start_i - 1)->dbl_score;

					// Note that for backward variable we do a transition from current to next.
					double next_var_cM = cur_win_var_regs->at(var_start_i)->dbl_score;

					double r_m = fabs(cur_var_cM - next_var_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

					double state_trans_log_prob = xlog(other_prob);
					if (state_i == first_state_i)
					{
						state_trans_log_prob = xlog(self_prob);
					}

					/////////////////////////
					// Get k-mer transition score.
					double shifted_kmer_dist = get_shifted_kmer_distance(per_win_all_ref_kmers[var_start_i-1]->at(first_state_i), per_win_all_ref_kmers[var_start_i]->at(state_i));
					double kmer_concordance_log_prob = -1 * kmer_concordance_log_weight * shifted_kmer_dist / (per_win_all_ref_kmers[var_start_i-1]->at(state_i)->kmer_length);

					// Emission probability comes from the buffered emission probability at the first window for the current state.
					// Note that the emissions are not necessary for fore array since we emit the current window in fore array but back array emits next window.
					double emission_concordance_log_prob = xlog_mul(kmer_concordance_log_prob, per_win_query_per_ref_hap_query_emission_prob[var_start_i]->at(state_i));

					/////////////////////////
					double state_kmer_transition_prob = xlog_mul(emission_concordance_log_prob, state_trans_log_prob);

					/////////////////////////
					total_foreback_score_per_back = xlog_sum(total_foreback_score_per_back,
															xlog_mul(state_kmer_transition_prob,
																	back_scores_per_hap[var_start_i][proxy_query_hap_i][state_i]));
				} // first_state_i loop.
			} // state_i loop.

			//double total_max_proxy_score_per_full_hap = top_scoring_kmer_score;

			fprintf(stderr, "Haplotype %d total probabilities: Fore=%.5f; Back=%.5f\n", 
					proxy_query_hap_i, total_foreback_score_per_fore, total_foreback_score_per_back);

			// Calculate allele probability at each variant: Use the fore-back arrays at each position.
			int n_errs = 0;
			int n_max_AF_errs = 0;
			for (int var_i = var_start_i; var_i < var_end_i; var_i++)
			{
				double per_allele_prob[2];
				per_allele_prob[0] = xlog(0.0);
				per_allele_prob[1] = xlog(0.0);

				// Go over all states and assign the per allele probabilities for the current window's center variant.
				for (int i_state = 0; i_state < n_ref_haplotypes; i_state++)
				{
					double cur_state_foreback_score = xlog_mul(fore_scores_per_hap[var_i][proxy_query_hap_i][i_state],
																back_scores_per_hap[var_i][proxy_query_hap_i][i_state]);

					// Get the reference subject index and haplotype index.
					int ref_i_s = i_state / 2;
					int ref_i_hap = i_state % 2;

					// Get the variant allele from reference haplotype.
					char cur_state_var_allele = per_ref_subj_orig_haplotypes->at(ref_i_s)[ref_i_hap][var_i];

					// Update the allele probability at this variant's corresponding allele.
					per_allele_prob[(int)cur_state_var_allele] = xlog_sum(per_allele_prob[(int)cur_state_var_allele], cur_state_foreback_score);
				} // i_state loop.

				// Select the map allele.
				char map_foreback_allele = 0;
				if (per_allele_prob[1] > per_allele_prob[0])
				{
					map_foreback_allele = 1;
				}

				bool found_kmer = false;
				char cur_subj_orig_allele = per_query_subj_orig_haplotypes->at(proxy_query_i_s)[proxy_query_hap_i][var_i];
				if (map_foreback_allele != cur_subj_orig_allele)
				{
					n_errs++;
				}
				else
				{
					found_kmer = true;
				}

				int max_AF_allele = ref_orig_geno_var_regs->at(var_i)->score;

				bool found_max_AF_allele = false;
				if (max_AF_allele != cur_subj_orig_allele)
				{
					n_max_AF_errs++;
				}
				else
				{
					found_max_AF_allele = true;
				}

				// Write the per subject info: Note that the backtracking does not work here, i.e., we do not have a way to track the AF of kmers in ForeBack
				// This is because ForeBack works at the variant level using per variant allele posterior probabilities. It does not do a complete
				// traceback of the haplotype.
				double backtracked_ref_hap_AF = 0.0;

				fprintf(f_per_sample_per_var_dec_stats, "%s\t%d\t%s\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n",
					query_sample_ids->at(proxy_query_i_s), proxy_query_hap_i,
					query_proxized_geno_var_regs->at(var_i)->chrom, query_proxized_geno_var_regs->at(var_i)->start,
					per_win_n_ref_kmers[var_i],
					per_win_query_proxy_kmer_freq[var_i], backtracked_ref_hap_AF, (int)found_kmer, (int)found_max_AF_allele);
				fflush(f_per_sample_per_var_dec_stats);
			} // var_i loop.

			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_foreback_score_per_fore, total_foreback_score_per_back);
			fprintf(f_dec_stats, "%s\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n", op_prefix, proxy_query_i_s, proxy_query_hap_i, n_errs, n_max_AF_errs, total_foreback_score_per_fore, total_foreback_score_per_back);
			fflush(f_dec_stats);

			// Free memory, this is necessary for long runs.
			for (int i_var = var_start_i; i_var < var_end_i; i_var++)
			{
				// This is the max log-probability score for each unique kmer in the window.
				//double** per_win_per_ref_kmer_log_cumul_probs = new double* [query_proxized_geno_var_regs->size()];
				//delete[] per_win_per_ref_kmer_log_cumul_probs[i_var];

				// # of unique ref. kmers in the window.
				//int* per_win_n_ref_kmers = new int[query_proxized_geno_var_regs->size()];

				// List of unique ref. kmers per window.
				//vector<t_kmer*>** per_win_unique_ref_kmers = new vector<t_kmer*>*[query_proxized_geno_var_regs->size()];
				delete_kmers(per_win_unique_ref_kmers[i_var]);
				delete_kmers(per_win_all_ref_kmers[i_var]);
				delete(per_win_unique_ref_kmer_i_per_orig_ref_kmer_i[i_var]);

				// The best scoring unique kmer index in the previous window for each current window.
				//int** per_win_backtracking_kmer_i = new int* [query_proxized_geno_var_regs->size()];
				//delete[] per_win_backtracking_kmer_i[i_var];

				// The frequency of the proxy kmer at each window. This is query sample specific, i.e., one value at each window.
				//double* per_win_query_proxy_kmer_freq = new double[query_proxized_geno_var_regs->size()];

				// The number of ref unique kmer counts at each window; this is used for calculating kmer frequencies in the window.
				//vector<int>** per_win_per_unique_ref_kmer_cnts = new vector<int>*[query_proxized_geno_var_regs->size()];
				delete(per_win_per_unique_ref_kmer_cnts[i_var]);

				delete(per_win_query_per_ref_hap_query_emission_prob[i_var]);
			} // i_var loop.

			delete[] per_win_unique_ref_kmers;
			delete[] per_win_per_unique_ref_kmer_cnts;
			delete[] per_win_query_proxy_kmer_freq;
			delete[] per_win_unique_ref_kmer_i_per_orig_ref_kmer_i;
			delete[] per_win_n_ref_kmers;
			delete[] per_win_query_per_ref_hap_query_emission_prob;
			delete[] per_win_all_ref_kmers;
		} // query_hap_i loop.
	} // test_sample_i loop.

	close_f(f_dec_stats, dec_summary_stats_fp);
	close_f(f_per_sample_per_var_dec_stats, per_sample_per_var_dec_stats_fp);

	fprintf(stderr, "Done!\n");
} // decode_site_alleles_per_proxized_reference_per_LiStephens_Histogram_matching_ForeBack function.