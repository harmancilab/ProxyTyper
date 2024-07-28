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
#include "prxytypr_proxytyper.h"

#include <vector>
#include <algorithm>
#include<functional>

using namespace std;

const bool __DUMP_VARIANT_MATRIX_FACTORIZATION__ = false;

bool is_all_zero_row(double* row, int numCols, double zero_check_eps)
{
	for (int i_col = 0; i_col < numCols; i_col++)
	{
		if (fabs(row[i_col]) > zero_check_eps)
		{
			return(false);
		}
	} // i_col loop.

	return(true);
}

int get_pivot_row_i(double* row, int numCols, double eps)
{
	for (int i_col = 0; i_col < numCols; i_col++)
	{
		if (fabs(row[i_col]) > eps)
		{
			return(i_col);
		}
	} // i_col loop.

	return(-1);
}

bool sort_two_rows_per_increasing_pivot_index(void** row_info1, void** row_info2)
{
	//double* row1 = (double*)(row_info1[0]);
	//double* row2 = (double*)(row_info2[0]);

	int* row1_pivot_i_ptr = (int*)(row_info1[2]);
	int* row2_pivot_i_ptr = (int*)(row_info2[2]);

	return(row1_pivot_i_ptr[0] < row2_pivot_i_ptr[0]);
}

bool sort_two_rows_per_decreasing_pivot_index(void** row_info1, void** row_info2)
{
	//double* row1 = (double*)(row_info1[0]);
	//double* row2 = (double*)(row_info2[0]);

	int* row1_pivot_i_ptr = (int*)(row_info1[2]);
	int* row2_pivot_i_ptr = (int*)(row_info2[2]);

	return(row1_pivot_i_ptr[0] > row2_pivot_i_ptr[0]);
}

double** sort_matrix_rows_per_pivot_index(double** matrix, int numRows, int numCols, double zero_eps, int* sorted_row_i, bool sort_increasing)
{
	vector<void**>* rows_info_vector = new vector<void**>();

	for (int i_row = 0; i_row < numRows; i_row++)
	{
		void** cur_row_info = new void* [3];
		cur_row_info[0] = matrix[i_row];

		int cur_row_pivot_i = get_pivot_row_i(matrix[i_row], numCols, zero_eps);

		int* row_i_ptr = new int[2];
		row_i_ptr[0] = i_row;
		cur_row_info[1] = row_i_ptr;

		int* pivot_i_ptr = new int[2];
		pivot_i_ptr[0] = cur_row_pivot_i;
		cur_row_info[2] = pivot_i_ptr;

		rows_info_vector->push_back(cur_row_info);
	} // i_row loop.

	if (sort_increasing)
	{
		sort(rows_info_vector->begin(), rows_info_vector->end(), sort_two_rows_per_increasing_pivot_index);
	}
	else
	{
		sort(rows_info_vector->begin(), rows_info_vector->end(), sort_two_rows_per_decreasing_pivot_index);
	}


	double** sorted_matrix = allocate_matrix(numRows, numCols);
	for (int i_row = 0; i_row < numRows; i_row++)
	{
		void** cur_row_info = rows_info_vector->at(i_row);
		int* cur_row_i_ptr = (int*)(cur_row_info[1]);
		sorted_row_i[i_row] = cur_row_i_ptr[0];

		sorted_matrix[i_row] = (double*)(cur_row_info[0]);
	} // i_row loop.

	return(sorted_matrix);
}

// This function returns the list of independent rows from a matrix.
vector<int>* matrix_get_independent_row_indices(double** matrix, int numRows, int numCols)
{
	/*double** tempMatrix = (double**)malloc(numRows * sizeof(double*));
	for (int i = 0; i < numRows; i++) {
		tempMatrix[i] = (double*)malloc(numCols * sizeof(double));
		for (int j = 0; j < numCols; j++) {
			tempMatrix[i][j] = matrix[i][j];
		}
	}*/

	// 
	double zero_check_eps = pow(10, -3);

	double** tempMatrix = allocate_copy_matrix(matrix, numRows, numCols);

	double** trans_tempMatrix = transpose_matrix(tempMatrix, numRows, numCols, NULL);

	// Sort the columns to push first col first.
	int* sorted_col_indices = new int[numCols];
	double** col_sorted_trans_tempMatrix = sort_matrix_rows_per_pivot_index(trans_tempMatrix, numCols, numRows, zero_check_eps, sorted_col_indices, false);

	if (__DUMP_VARIANT_MATRIX_FACTORIZATION__)
	{
		fprintf(stderr, "Col sorted trans matrix:\n");
		print_matrix(col_sorted_trans_tempMatrix, numCols, numRows);
	}

	double** trans_col_sorted_trans_tempMatrix = transpose_matrix(col_sorted_trans_tempMatrix, numCols, numRows, NULL);

	int* sorted_row_indices = new int[numRows];
	//double** sorted_matrix = sort_matrix_rows_per_pivot_index(tempMatrix, numRows, numCols, sorted_row_indices, true);
	double** sorted_matrix = sort_matrix_rows_per_pivot_index(trans_col_sorted_trans_tempMatrix, numRows, numCols, zero_check_eps, sorted_row_indices, true);

	if (__DUMP_VARIANT_MATRIX_FACTORIZATION__)
	{
		fprintf(stderr, "Final Row/Col Pivot sorted matrix:\n");
		for (int row_i = 0; row_i < numRows; row_i++)
		{
			fprintf(stderr, "%d:", sorted_row_indices[row_i]);
			for (int col_i = 0; col_i < numCols; col_i++)
			{
				fprintf(stderr, "\t%.3f", sorted_matrix[row_i][col_i]);
			}
			fprintf(stderr, "\n");
		}
	}
	//print_matrix(sorted_matrix, numRows, numCols);

	//exit();

	tempMatrix = sorted_matrix;

	int* row_indices = sorted_row_indices;
	//int* row_indices = new int[numRows];
	//for (int i = 0; i < numRows; i++)
	//{
	//	row_indices[i] = i;
	//} // i loop.


	// This loop goes over all the rows once for the whole process.
	for (int i = 0; i < numRows; i++)
	{
		// Find pivot element: Looks for non-zero element at the i^th column for all the remaining rows.
		// It looks for the first row with a non-zero value at diagonal i'th column, which will be swapped to i^th row.
		// If there is no row like this, it moves forward.
		// This misses the fact that 
		int pivotRow = i;
		while (pivotRow < numRows)
		{
			if (is_all_zero_row(tempMatrix[pivotRow], numCols, zero_check_eps))
			{
				pivotRow++;
			}
			else if (i < numCols)
			{
				if (fabs(tempMatrix[pivotRow][i]) < zero_check_eps)
				{
					pivotRow++;
				}
				else
				{
					// The basic idea is that before we hit the out-of-bounds, we should find a pivot or all rows must be zero as check by above.
					break;
				}
			}
			else if (i >= numCols)
			{
				// if this is not an all_zero_row, there is a problem with this matrix.
				fprintf(stderr, "Sanity check error..\n");
				exit(1);
			}
		}

		if (__DUMP_VARIANT_MATRIX_FACTORIZATION__)
		{
			fprintf(stderr, "i:%d/%d, pivotRow: %d\n", i, numRows, pivotRow);
		}

		if (pivotRow == numRows)
		{
			// All remaining rows are zero in this column, move to next column
			continue;
		}

		if (pivotRow != i)
		{
			// Swap rows to bring pivot element to diagonal
			double* temp = tempMatrix[i];
			tempMatrix[i] = tempMatrix[pivotRow];
			tempMatrix[pivotRow] = temp;

			// Also swap the indices for these rows.
			int temp_row_i = row_indices[i];
			row_indices[i] = row_indices[pivotRow];
			row_indices[pivotRow] = temp_row_i;
		}

		// Eliminate column below pivot element
		// This is different from original procedure.
		//for (int j = i + 1; j < numRows; j++)
		for (int j = 0; j < numRows; j++)
		{
			if (j != i)
			{
				double factor = tempMatrix[j][i] / tempMatrix[i][i];

				if (tempMatrix[i][i] == 0)
				{
					fprintf(stderr, "Sanity check failed: Error; diag is 0\n");
					exit(1);
				}

				//for (int k = i; k < numCols; k++) 
				for (int k = 0; k < numCols; k++)
				{
					tempMatrix[j][k] -= factor * tempMatrix[i][k];

					// Following enforces small values to 0 to make sure we do not have numerical issues.
					if (fabs(tempMatrix[j][k]) < zero_check_eps)
					{
						tempMatrix[j][k] = 0;
					}
				}
			}
		} // j loop.

		if (__DUMP_VARIANT_MATRIX_FACTORIZATION__)
		{
			fprintf(stderr, "i: %d; temp_matrix:\n", i);
			for (int row_i = 0; row_i < numRows; row_i++)
			{
				fprintf(stderr, "%d:", row_indices[row_i]);
				for (int col_i = 0; col_i < numCols; col_i++)
				{
					fprintf(stderr, "\t%.3f", tempMatrix[row_i][col_i]);
				}
				fprintf(stderr, "\n");
			}
		}
		//print_matrix(tempMatrix, numRows, numCols);
	} // i loop.

	// Find independent rows
	//int* independentRows = (int*)malloc(numRows * sizeof(int));
	//int* independentRows = new int[numRows];
	vector<int>* independentRows = new vector<int>();

	//int numIndependentRows = 0;
	for (int i = 0; i < numRows; i++) {
		int allZeroes = 1;

		//for (int j = 0; j < numCols - 1; j++)
		for (int j = 0; j < numCols; j++)
		{
			//if (tempMatrix[i][j] != 0) 
			if (fabs(tempMatrix[i][j]) > zero_check_eps)
			{
				allZeroes = 0;
				break;
			}
		}

		//if (!allZeroes || tempMatrix[i][numCols - 1] != 0) 
		if (!allZeroes)
		{
			independentRows->push_back(row_indices[i]);
			//independentRows[numIndependentRows] = i;
			//numIndependentRows++;
		}
	}

	delete_matrix(tempMatrix, numRows, numCols);

	//fprintf(stderr, "%d independent rows\n", numIndependentRows);

	//return (int*)realloc(independentRows, numIndependentRows * sizeof(int));
	return(independentRows);
}


void process_regions_per_callback(vector<t_annot_region*>* regions, function<void(t_annot_region*)>& callback)
{
	for (size_t i_reg = 0; i_reg < regions->size(); i_reg++)
	{
		callback(regions->at(i_reg));
	} // i_reg loop.
}

bool matrix_is_identity(double** ident_check_mat, int nrow, int ncol, double ident_check_eps)
{
	for (int row = 0; row < nrow; row++)
	{
		for (int col = 0; col < ncol; col++)
		{
			if (row == col)
			{
				if (fabs(ident_check_mat[row][col] - 1.0) > ident_check_eps)
				{
					fprintf(stderr, "M[%d][%d]=%.5f\n", row, col, ident_check_mat[row][col]);
					return(false);
				}
			}
			else
			{
				if (fabs(ident_check_mat[row][col]) > ident_check_eps)
				{
					fprintf(stderr, "M[%d][%d]=%.5f\n", row, col, ident_check_mat[row][col]);
					return(false);
				}
			}
		} // col loop.
	} // row loop.

	return(true);
}

double** get_random_allele1_2_hap_matrix(t_rng* rng, int n_block_target_vars, int n_panel_haplotypes, double allele1_prob)
{
	double** rand_allele1_2_hap_matrix = allocate_matrix(n_block_target_vars, n_panel_haplotypes);

	// Total random hap matrix, this is a dense matrix.
	for (int i_var = 0; i_var < n_block_target_vars; i_var++)
	{
		for (int i_hap = 0; i_hap < n_panel_haplotypes; i_hap++)
		{
			double cur_val = (rng->random_double_ran3() < allele1_prob) ? (1.0) : (0.0);
			rand_allele1_2_hap_matrix[i_var][i_hap] = cur_val;
		} // i_hap loop. 
	} // i_var loop.

	return(rand_allele1_2_hap_matrix);

	//int n_haps_in_list = MAX(n_block_target_vars, n_panel_haplotypes);

	//// Add each variant to at least to 1 haplotype.
	//vector<int>* hap_i_list = new vector<int>();
	////int cur_i_hap = 0;
	//for (int i_hap = 0; i_hap < n_haps_in_list; i_hap++)
	//{
	//	// Select the haplotype:
	//	int cur_i_hap = i_hap % n_panel_haplotypes;
	//	hap_i_list->push_back(cur_i_hap);
	//} // i_var loop.

	//// Shuffle haplotypes to be distributed to variants.
	//vector<int>* rand_i = rng->permute_indices(n_haps_in_list, n_haps_in_list);
	//int i_var = 0;
	//for (int i_hap = 0; i_hap < n_haps_in_list; i_hap++)
	//{
	//	int cur_i_var = i_hap % n_block_target_vars;
	//	int cur_i_hap = hap_i_list->at(rand_i->at(i_hap));
	//	rand_allele1_2_hap_matrix[cur_i_var][cur_i_hap] = 1;
	//} // i_var loop.
}


void generate_untyped_variant_recoded_reference_panel(char* ref_panel_tag_haplocoded_matbed_fp, char* ref_panel_target_haplocoded_matbed_fp, char* panel_sample_list_fp,
	char* op_dir, char* op_prefix)
{
	fprintf(stderr, "THIS FUNCTION IS NOT USED ANY MORE SINCE WE NEED THE LEFT INVERSE OF B THAT IS VERY LARGE IN DIMENSIONS.\n");
	fprintf(stderr, "THIS FUNCTION IS NOT USED ANY MORE SINCE WE NEED THE LEFT INVERSE OF B THAT IS VERY LARGE IN DIMENSIONS.\n");
	fprintf(stderr, "THIS FUNCTION IS NOT USED ANY MORE SINCE WE NEED THE LEFT INVERSE OF B THAT IS VERY LARGE IN DIMENSIONS.\n");
	fprintf(stderr, "THIS FUNCTION IS NOT USED ANY MORE SINCE WE NEED THE LEFT INVERSE OF B THAT IS VERY LARGE IN DIMENSIONS.\n");
	fprintf(stderr, "THIS FUNCTION IS NOT USED ANY MORE SINCE WE NEED THE LEFT INVERSE OF B THAT IS VERY LARGE IN DIMENSIONS.\n");

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

	for (size_t i_tag = 0; i_tag < panel_tag_geno_sig_regs->size(); i_tag++)
	{
		panel_tag_geno_sig_regs->at(i_tag)->score = 0;
	} // i_tag loop.

	for (size_t i_target = 0; i_target < panel_target_geno_sig_regs->size(); i_target++)
	{
		panel_target_geno_sig_regs->at(i_target)->score = 1;
	} // i_target loop.

	//process_regions_per_callback(panel_tag_geno_sig_regs, [](t_annot_region* reg) {reg->score = 0; });
	//process_regions_per_callback(panel_target_geno_sig_regs, [](t_annot_region* reg) {reg->score = 1; });

	// Add tag and target regions to the list of all regions.
	vector<t_annot_region*>* tag_target_geno_sig_regs = new vector<t_annot_region*>();
	tag_target_geno_sig_regs->insert(tag_target_geno_sig_regs->end(), panel_tag_geno_sig_regs->begin(), panel_tag_geno_sig_regs->end());
	tag_target_geno_sig_regs->insert(tag_target_geno_sig_regs->end(), panel_target_geno_sig_regs->begin(), panel_target_geno_sig_regs->end());

	// Sort and restructure regions.
	t_restr_annot_region_list* restr_tag_target_geno_sig_regs = restructure_annot_regions(tag_target_geno_sig_regs);

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	// We have the allele1-2-haplotype mapping matrix, A; A=TxB
	// For this, we selected B as a random permutation matrix -> This corresponds to censoring.
	// We can also select a non-singular B matrix -> then get T by AxinvR(B).
	// Number of untyped: n ; number of haplotypes: s ; A->[nxs], B->[n'xs], T->[nxn'] 
	// Allele-1 probabilities: Phap(1)->[sx1] (fore-back probs for each hap); P(1)->[nx1] where P(1)=AxPhap(1)
	// We do imputation on variants defined by B. So we will have [n'x1] probabilities for one allele; which is calculated as BxPhap(1)
	// We will use TxP(1) to map the final probabilities.
	for (size_t i_chr = 0; i_chr < restr_tag_target_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_tag_target_regs = restr_tag_target_geno_sig_regs->regions_per_chrom[i_chr];

		fprintf(stderr, "Processing %d tag/targets on %s\n", (int)cur_chr_tag_target_regs->size(),
			restr_tag_target_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_block_target_regs = new vector<t_annot_region*>();
		t_annot_region* cur_block_starting_tag_reg = NULL;
		int untyped_block_i = 0;
		for (size_t i_reg = 0; i_reg < cur_chr_tag_target_regs->size(); i_reg++)
		{
			// Find the first tag region.
			if (cur_chr_tag_target_regs->at(i_reg)->score == 0 &&
				cur_block_target_regs->size() > 0)
			{
				if (cur_block_starting_tag_reg == NULL)
				{
					fprintf(stderr, "Found first tag variant: %s:%d-%d::%s\n", cur_chr_tag_target_regs->at(i_reg)->chrom,
						cur_chr_tag_target_regs->at(i_reg)->start, cur_chr_tag_target_regs->at(i_reg)->end,
						cur_chr_tag_target_regs->at(i_reg)->name);
				}
				else
				{
					fprintf(stderr, "Found block: %s:%d-%d [%s] -- %s:%d-%d [%s]\n",
						cur_block_starting_tag_reg->chrom,
						cur_block_starting_tag_reg->start, cur_block_starting_tag_reg->end,
						cur_block_starting_tag_reg->name,
						cur_chr_tag_target_regs->at(i_reg)->chrom,
						cur_chr_tag_target_regs->at(i_reg)->start, cur_chr_tag_target_regs->at(i_reg)->end,
						cur_chr_tag_target_regs->at(i_reg)->name);
				}

				// Process the targets in the current block.
				fprintf(stderr, "Processing %d target regions.\n", (int)cur_block_target_regs->size());

				int n_block_target_vars = (int)cur_block_target_regs->size();

				// Extract the allele1-2-haplotypes matrix.
				double** A = allocate_matrix(n_block_target_vars, n_panel_haplotypes);

				for (int i_target_var = 0; i_target_var < (int)cur_block_target_regs->size(); i_target_var++)
				{
					void** cur_target_reg_info = (void**)(cur_block_target_regs->at(i_target_var)->data);
					char* cur_target_var_geno_sig = (char*)(cur_target_reg_info[0]);

					int panel_hap_i = 0;
					for (int i_s = 0; i_s < (int)panel_sample_ids->size(); i_s++)
					{
						for (int i_hap = 0; i_hap < 2; i_hap++)
						{
							A[i_target_var][panel_hap_i] = get_allele_per_haplotype(cur_target_var_geno_sig[i_s], i_hap);

							panel_hap_i++;
						} // i_hap loop.
					} // i_s loop.
				} // i_target_var loop.			

				vector<int>* independent_row_A_i = matrix_get_independent_row_indices(A, n_block_target_vars, n_panel_haplotypes);

				fprintf(stderr, "Extracted allele1_2_hap matrix (%d/%d independent variants).\n",
					(int)independent_row_A_i->size(), n_block_target_vars);

				// Select a random B binary matrix that is sparse.
				double** raw_B = get_random_allele1_2_hap_matrix(rng, n_block_target_vars, n_panel_haplotypes, 0.3);
				vector<int>* independent_row_B_i = matrix_get_independent_row_indices(raw_B, n_block_target_vars, n_panel_haplotypes);

				// Remove the linearly dependent rows from raw B matrix.
				int n_B_rows = (int)independent_row_B_i->size();
				fprintf(stderr, "Generated %d independent rows for B.\n", n_B_rows);
				double** B = allocate_matrix(n_B_rows, n_panel_haplotypes);
				int cur_ind_row_i = 0;
				for (int i_var = 0; i_var < (int)independent_row_B_i->size(); i_var++)
				{
					//fprintf(stderr, "%d\n", independent_row_B_i->at(i_var));

					int cur_row_i = independent_row_B_i->at(i_var);
					for (int i_hap = 0; i_hap < n_panel_haplotypes; i_hap++)
					{
						B[cur_ind_row_i][i_hap] = raw_B[cur_row_i][i_hap];
					} // i_hap loop.

					cur_ind_row_i++;
				} // i_var loop.

				double** Bt = transpose_matrix(B, n_B_rows, n_panel_haplotypes, NULL);

				fprintf(stderr, "Calculated B and Bt.\n");

				double** BBt = matrix_multiply(B, n_B_rows, n_panel_haplotypes, Bt, n_panel_haplotypes, n_B_rows, NULL);
				double** invBBt = invert_matrix_GJ(BBt, n_B_rows, n_B_rows, NULL);
				double** ident_check_mat_BBt = matrix_multiply(BBt, n_B_rows, n_B_rows, invBBt, n_B_rows, n_B_rows, NULL);
				double ident_check_eps = pow(10, -6);
				if (!matrix_is_identity(ident_check_mat_BBt, n_B_rows, n_B_rows, ident_check_eps))
				{
					fprintf(stderr, "BBt inverse is not correctly computed.\n");
					exit(1);
				}

				fprintf(stderr, "Calculated invBBt.\n");

				// Calculate T using right pseudoinverse of B: Br=B'(BB')^-1; T=AxBr
				// For simplicity, assume n'=n for now.
				// Note that Br has a large calculation; B' -> sxn ; BB' -> nxn ; thus Br is sxn
				// So T=AxBr is nxn.
				double** invr_B = matrix_multiply(Bt, n_panel_haplotypes, n_B_rows, invBBt, n_B_rows, n_B_rows, NULL);

				// Do a sanity check on the right inverse.
				double** ident_check_mat = matrix_multiply(B, n_B_rows, n_panel_haplotypes, invr_B, n_panel_haplotypes, n_B_rows, NULL);
				//double ident_check_eps = pow(10, -6);
				if (!matrix_is_identity(ident_check_mat, n_B_rows, n_B_rows, ident_check_eps))
				{
					fprintf(stderr, "Right inverse is not correctly computed.\n");
					exit(1);
				}

				//fprintf(stderr, "Looks good..\n"); exit(1);

				fprintf(stderr, "Calculated right inverse of B.\n");

				double** T = matrix_multiply(A, n_block_target_vars, n_panel_haplotypes, invr_B, n_panel_haplotypes, n_B_rows, NULL);

				fprintf(stderr, "Calculated T, saving.\n");

				// Save T to be used later.
				char T_matrix_fp[1000];
				//sprintf(T_matrix_fp, "%s_%d_T.mat.gz", op_dir, cur_block_starting_tag_reg->chrom, cur_block_starting_tag_reg->start);
				sprintf(T_matrix_fp, "%s/%s_untyped_block_%d_T.mat.gz", op_dir, op_prefix, untyped_block_i);

				fprintf(stderr, "Calculated T, saving to %s\n", T_matrix_fp);
				save_matrix_binary(T, n_block_target_vars, n_B_rows, T_matrix_fp);

				fprintf(stderr, "Setting the proxized target variants using B.\n");

				//// Generate the new proxy target variant genotypes: Use B to generate the signals: Var-2-hap mapping.
				//for (int i_target_var = 0; i_target_var < cur_block_target_regs->size(); i_target_var++)
				//{
				//	void** cur_target_reg_info = (void**)(cur_block_target_regs->at(i_target_var)->data);
				//	char* cur_target_var_geno_sig = (char*)(cur_target_reg_info[0]);

				//	char* cur_target_var_proxy_geno_sig = new char[panel_sample_ids->size() + 1];
				//	memset(cur_target_var_proxy_geno_sig, 0, sizeof(char) * (panel_sample_ids->size() + 1));

				//	//char* cur_target_var_proxy_geno_sig = (char*)(cur_target_reg_info[1]);

				//	for (int i_hap = 0; i_hap < n_panel_haplotypes; i_hap++)
				//	{
				//		int i_s = i_hap / 2;
				//		int cur_i_s_hap_i = i_hap % 2;
				//		char cur_geno = cur_target_var_proxy_geno_sig[i_s];

				//		if (i_s >= panel_sample_ids->size())
				//		{
				//			fprintf(stderr, "Sanity check failed: Sample index out of bounds.\n");
				//			exit(1);
				//		}

				//		// Set B to the correct allele of the genotype.
				//		char cur_B_int_val = (char)(B[i_target_var][i_hap]);
				//		cur_geno = cur_geno | (cur_B_int_val << cur_i_s_hap_i);

				//		// Update the genotype.
				//		cur_target_var_proxy_geno_sig[i_s] = cur_geno;
				//	} // i_hap loop.
				//} // i_target_var loop.

				// Clear the target regions.
				cur_block_target_regs->clear();

				untyped_block_i++;

				cur_block_starting_tag_reg = cur_chr_tag_target_regs->at(i_reg);
			}
			else if (cur_chr_tag_target_regs->at(i_reg)->score == 1)
			{
				cur_block_target_regs->push_back(cur_chr_tag_target_regs->at(i_reg));
			}
		} // i_reg loop.
	} // i_chr loop.

	// Replace the signal and save.
}