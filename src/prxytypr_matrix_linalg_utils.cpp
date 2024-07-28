#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_utils.h"
#include "prxytypr_histogram.h"
#include "prxytypr_ansi_string.h"
//#include "prxytypr_matrix.h"
#include <math.h>
#include <string.h>

using namespace std;

bool __DUMP_LINALG_MESSAGES__ = false;

/*
Simple matrix library for processing matrices for development purposes. 
!!! All indices start with 0 (not 1 as in t_matrix class). !!! 
!!! All indices start with 0 (not 1 as in t_matrix class). !!! 
!!! All indices start with 0 (not 1 as in t_matrix class). !!! 
*/

double** allocate_matrix(int nrow, int ncol)
{
	double** matrix = new double*[nrow];
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		matrix[i_row] = new double[ncol + 1];
		memset(matrix[i_row], 0, sizeof(double) * ncol);
	} // i_row loop.

	return(matrix);
}

void delete_matrix(double** matrix, int n_row, int n_col)
{
	for (int i_row = 0; i_row < n_row; i_row++)
	{
		delete[] matrix[i_row];
	} // i_row loop.

	delete[] matrix;
}

double** allocate_matrix(int nrow, int ncol, double init_val)
{
	double** matrix = new double*[nrow];
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		matrix[i_row] = new double[ncol + 1];
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			matrix[i_row][i_col] = init_val;
		}
	} // i_row loop.

	return(matrix);
}

double** transpose_matrix(double** matrix, int nrow, int ncol, double** res_matrix)
{
	double** trans_matrix = res_matrix;
	if (trans_matrix == NULL)
	{
		trans_matrix = allocate_matrix(ncol, nrow);
	}

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			trans_matrix[i_col][i_row] = matrix[i_row][i_col];
		} // i_col loop.
	} // i_row loop.

	return(trans_matrix);
}

double** get_matrix_per_col_vector(vector<double>* row_vector, double** res_matrix)
{
	double** col_matrix = res_matrix;
	int nrow = (int)row_vector->size();
	if (col_matrix == NULL)
	{
		col_matrix = allocate_matrix(nrow, 1);
	}

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		col_matrix[i_row][0] = row_vector->at(i_row);
	} // i_row loop.

	return(col_matrix);
}

double** get_matrix_per_col_vector(double* row_vector, int nrow, double** res_matrix)
{
	double** col_matrix = res_matrix;
	if (col_matrix == NULL)
	{
		col_matrix = allocate_matrix(nrow, 1);
	}
	
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		col_matrix[i_row][0] = row_vector[i_row];
	}

	return(col_matrix);
}

double** get_matrix_per_row_vector_list(vector<double*>* row_vector_list, int n_col, double** res_matrix)
{
	int n_row = (int)row_vector_list->size();
	double** matrix = res_matrix;
	if (matrix == NULL)
	{
		matrix = allocate_matrix(n_row, n_col);
	}	
	
	for (int i_row = 0; i_row < n_row; i_row++)
	{
		for (int i_col = 0; i_col < n_col; i_col++)
		{
			matrix[i_row][i_col] = row_vector_list->at(i_row)[i_col];
		} // i_col loop.
	} // i_row loop.

	return(matrix);
}

double** allocate_copy_matrix(double** matrix, int nrow, int ncol)
{
	double** matrix_c = allocate_matrix(nrow, ncol);
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			matrix_c[i_row][i_col] = matrix[i_row][i_col];
		} // i_col loop.
	} // i_row loop.

	return(matrix_c);
}


void print_matrix(double** matrix, int nrow, int ncol)
{
	for (int i = 0; i < nrow; i++)
	{
		fprintf(stderr, "%d:", i);

		for (int j = 0; j < ncol; j++)
		{
			fprintf(stderr, "\t%.3f", matrix[i][j]);
		} // j loop.

		fprintf(stderr, "\n");
	} // i loop.
}

void copy_matrix(double** src_matrix, int src_nrow, int src_ncol,
	double** dest_matrix, int dest_nrow, int dest_ncol)
{
	if (src_nrow > dest_nrow ||
		src_ncol > dest_ncol)
	{
		fprintf(stderr, "Non-conformant dimensions for matrix copy: %d, %d; %d, %d", src_nrow, src_ncol, dest_nrow, dest_ncol);
		exit(1);
	}

	for (int i_row = 0; i_row < src_nrow; i_row++)
	{
		for (int i_col = 0; i_col < src_ncol; i_col++)
		{
			dest_matrix[i_row][i_col] = src_matrix[i_row][i_col];
		} // i_col loop.
	} // i_row loop.
}

double** get_submatrix(double** src_matrix, 
	int src_row_i_2_start, int src_col_i_2_start, 
	int src_nrow, int src_ncol,
	int nrows_2_copy, int ncols_2_copy, double** dest_matrix_val)
{
	double** dest_matrix = dest_matrix_val;
	if(dest_matrix == NULL)
	{
		dest_matrix = allocate_matrix(nrows_2_copy, ncols_2_copy);
	}
	
	for (int i_row = 0; i_row < nrows_2_copy; i_row++)
	{
		for (int i_col = 0; i_col < ncols_2_copy; i_col++)
		{
			int i_row_2_copy = src_row_i_2_start + i_row;
			int i_col_2_copy = src_col_i_2_start + i_col;

			dest_matrix[i_row][i_col] = src_matrix[i_row_2_copy][i_col_2_copy];
		} // i_col loop.
	} // i_row loop.

	return(dest_matrix);
}

// https://www.codesansar.com/numerical-methods/matrix-inverse-using-gauss-jordan-cpp-output.htm
double** invert_matrix_GJ(double** matrix, int nrow, int ncol, double** res_matrix)
{
	if (ncol != nrow)
	{
		fprintf(stderr, "Cannot invert a non-square matrix.\n");
		exit(1);
	}

	if (__DUMP_LINALG_MESSAGES__)
	{
		fprintf(stderr, "Inverting:\n");
		print_matrix(matrix, nrow, ncol);
	}

	double** a = allocate_matrix(nrow, 2 * ncol);

	//float a[SIZE][SIZE], x[SIZE], ratio;
	double ratio;
	copy_matrix(matrix, nrow, ncol, a, nrow, 2 * ncol);
	int i, j, k, n;
	n = nrow;

	///* Setting precision and writing floating point values in fixed-point notation. */
	//cout << setprecision(3) << fixed;

	///* Inputs */
	///* 1. Reading order of matrix */
	//cout << "Enter order of matrix: ";
	//cin >> n;

	///* 2. Reading Matrix */
	//cout << "Enter coefficients of Matrix: " << endl;
	//for (i = 1;i <= n;i++)
	//{
	//	for (j = 1;j <= n;j++)
	//	{
	//		cout << "a[" << i << "]" << j << "]= ";
	//		cin >> a[i][j];
	//	}
	//}

	/* Augmenting Identity Matrix of Order n */
	for (i = 0;i < n;i++)
	{
		for (j = 0;j < n;j++)
		{
			if (i == j)
			{
				a[i][j + n] = 1;
			}
			else
			{
				a[i][j + n] = 0;
			}
		}
	}

	if (__DUMP_LINALG_MESSAGES__)
	{
		fprintf(stderr, "Copied, augmented:\n");
		print_matrix(a, nrow, 2 * nrow);
	}

	/* Applying Gauss Jordan Elimination */
	for (i = 0; i < n; i++)
	{
		if (a[i][i] == 0.0)
		{
			//cout << "Mathematical Error!";
			fprintf(stderr, "Mathematical error..\n");
			exit(1);
		}
		for (j = 0; j < n; j++)
		{
			if (i != j)
			{
				ratio = a[j][i] / a[i][i];

				for (k = 0;k < 2 * n;k++)
				{
					a[j][k] = a[j][k] - ratio * a[i][k];
				}
			}
		}
	}
	/* Row Operation to Make Principal Diagonal to 1 */
	for (i = 0; i < n; i++)
	{
		for (j = n; j < 2 * n; j++)
		{
			a[i][j] = a[i][j] / a[i][i];
		}
	}

	///* Displaying Inverse Matrix */
	double** invmatrix = res_matrix;
	if (invmatrix == NULL)
	{
		invmatrix = allocate_matrix(nrow, ncol);
	}

	//cout << endl << "Inverse Matrix is:" << endl;
	for (i = 0; i < n; i++)
	{
		for (j = n;j < 2 * n;j++)
		{
			invmatrix[i][j - n] = a[i][j];
			//cout << a[i][j] << "\t";
		}
		//cout << endl;
	}

	if (__DUMP_LINALG_MESSAGES__)
	{
		fprintf(stderr, "Inverted matrix:\n");
		print_matrix(invmatrix, nrow, ncol);
	}

	return(invmatrix);
}

double** get_diagonal_of_matrix(double** square_matrix, int nrow, int ncol, double** res_matrix)
{
	if (nrow != ncol)
	{
		fprintf(stderr, "Cannot return diagonal of a non-square matrix: %d, %d\n", nrow, ncol);
		exit(1);
	}

	double** res = res_matrix;
	if (res_matrix == NULL)
	{
		res = allocate_matrix(nrow, 1);
	}

	for (int i_diag = 0; i_diag < nrow; i_diag++)
	{
		res[i_diag][0] = square_matrix[i_diag][i_diag];
	}

	return(res);
}

double** process_matrix_elementwise_by_callback(double** matrix, int nrow, int ncol, double(*callback)(double), double** res_matrix)
{
	double** res = res_matrix;
	if (res_matrix == NULL)
	{
		res = allocate_matrix(nrow, ncol);
	}

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			res[i_row][i_col] = callback(matrix[i_row][i_col]);
		} // i_col loop.
	} // i_row loop.

	return(res);
}

double** get_diag_matrix(double** diag_elems, int nrow, int ncol, double** res)
{
	double** diag_matrix = res;

	if (nrow == 1)
	{		
		int l_diag = ncol;
		if (diag_matrix == NULL)
		{
			diag_matrix = allocate_matrix(l_diag, l_diag);
		}

		for (int diag_i = 0; diag_i < l_diag; diag_i++)
		{
			diag_matrix[diag_i][diag_i] = diag_elems[0][diag_i];
		} // diag_i loop.

	}
	else 
	{
		int l_diag = nrow;
		if (diag_matrix == NULL)
		{
			diag_matrix = allocate_matrix(l_diag, l_diag);
		}

		for (int diag_i = 0; diag_i < l_diag; diag_i++)
		{
			diag_matrix[diag_i][diag_i] = diag_elems[diag_i][0];
		} // diag_i loop.

	}

	return(diag_matrix);
}

double** get_diag_matrix(double* diag_elems, int l_diag, double** res_matrix)
{
	double** diag_matrix = res_matrix;
	if (diag_matrix == NULL)
	{
		diag_matrix = allocate_matrix(l_diag, l_diag);
	}

	for (int diag_i = 0; diag_i < l_diag; diag_i++)
	{
		diag_matrix[diag_i][diag_i] = diag_elems[diag_i];
	} // diag_i loop.

	return(diag_matrix);
}

double** matrix_add(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix)
{
	if (nArow != nBrow ||
		nAcol != nBcol)
	{
		fprintf(stderr, "Cannot add matrices, the dimensions non-conformant: %d, %d ; %d, %d\n", nArow, nAcol, nBrow, nBcol);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nArow, nBcol);
	}

	for (int row_i = 0; row_i < nArow; row_i++)
	{
		for (int col_i = 0; col_i < nBcol; col_i++)
		{
			res[row_i][col_i] = A[row_i][col_i] + B[row_i][col_i];
		}
	}

	return(res);
}

double** matrix_subtract(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix)
{
	if (nArow != nBrow ||
		nAcol != nBcol)
	{
		fprintf(stderr, "Cannot subtract matrices, the dimensions non-conformant: %d, %d ; %d, %d\n", nArow, nAcol, nBrow, nBcol);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nArow, nBcol);
	}

	for (int row_i = 0; row_i < nArow; row_i++)
	{
		for (int col_i = 0; col_i < nBcol; col_i++)
		{
			res[row_i][col_i] = A[row_i][col_i] - B[row_i][col_i];
		}
	}

	return(res);
}

double** matrix_multiply_by_scalar(double** A, int nrowA, int ncolA, double scalar, double** res_matrix)
{
	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row_i = 0; row_i < nrowA; row_i++)
	{
		for (int col_i = 0; col_i < ncolA; col_i++)
		{
			res[row_i][col_i] = A[row_i][col_i] * scalar;
		} // col_i loop.
	} // row_i loop.

	return(res);
}

double** matrix_row_by_row_inner_product(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix)
{
	if (nrowA != nrowB ||
		ncolA != ncolB)
	{
		fprintf(stderr, "Dimensions non-conformating for row-by-row multiplication that is not square: %d, %d; %d, %d\n", nrowA, ncolA, nrowB, ncolB);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, 1);
	}

	for (int i_row = 0; i_row < nrowA; i_row++)
	{
		res[i_row][0] = 0;
		for (int i_col = 0; i_col < ncolA; i_col++)
		{
			res[i_row][0] += (A[i_row][i_col] * B[i_row][i_col]);
		}
	} // i_row loop.

	return(res);
}

// Following uses nxn algorithm for fast multiplication with a diagonal matrix.
double** matrix_right_multiply_with_diag_matrix(double** A, int nrowA, int ncolA, double** diag, int nrowdiag, int ncoldiag, double** res_matrix)
{
	if (ncolA != nrowdiag)
	{
		fprintf(stderr, "Dimensions non-conformating for right-diag multiplication that is not square: %d, %d; %d, %d\n", nrowA, ncolA, nrowdiag, ncoldiag);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row = 0; row < nrowA; row++)
	{
		for (int col = 0; col < ncolA; col++)
		{
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
			res[row][col] = A[row][col] * diag[col][col];
		} // col loop.
	} // row loop.

	return(res);
}

double** matrix_right_multiply_with_diag_matrix_as_diagonal_only(double** A, int nrowA, int ncolA, double** diag, int nrowdiag, int ncoldiag, double** res_matrix)
{
	if (ncoldiag != 1 ||
		ncolA != nrowdiag)
	{
		fprintf(stderr, "Dimensions non-conformating for right-diag multiplication that is not square: %d, %d; %d, %d\n", nrowA, ncolA, nrowdiag, ncoldiag);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row = 0; row < nrowA; row++)
	{
		for (int col = 0; col < ncolA; col++)
		{
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
			res[row][col] = A[row][col] * diag[col][0];
		} // col loop.
	} // row loop.

	return(res);
}

// Following uses nxn algorithm for fast multiplication with a diagonal matrix.
double** matrix_left_multiply_with_diag_matrix(double** diag, int nrowdiag, int ncoldiag, double** A, int nrowA, int ncolA, double** res_matrix)
{
	if (ncoldiag != nrowA)
	{
		fprintf(stderr, "Dimensions non-conformating for right-diag multiplication that is not square: %d, %d; %d, %d\n", nrowA, ncolA, nrowdiag, ncoldiag);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row = 0; row < nrowA; row++)
	{
		for (int col = 0; col < ncolA; col++)
		{
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
			res[row][col] = A[row][col] * diag[row][row];
		} // col loop.
	} // row loop.

	return(res);
}

double** matrix_left_multiply_with_diag_matrix_as_diagonal_only(double** diag, int nrowdiag, int ncoldiag, double** A, int nrowA, int ncolA, double** res_matrix)
{
	if (ncoldiag != 1 ||
		ncolA != nrowdiag)
	{
		fprintf(stderr, "Dimensions non-conformating for right-diag multiplication that is not square: %d, %d; %d, %d\n", nrowA, ncolA, nrowdiag, ncoldiag);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row = 0; row < nrowA; row++)
	{
		for (int col = 0; col < ncolA; col++)
		{
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
//#error("MAKE SURE THIS ORDER IS CORRECT FOR RIGHT DIAG MULTIPLICATION")
			res[row][col] = A[row][col] * diag[row][0];
		} // col loop.
	} // row loop.

	return(res);
}

double** matrix_divide_elementwise(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix)
{
	if (nrowA != nrowB ||
		ncolA != ncolB)
	{
		fprintf(stderr, "Dimensions non-conformant for elementwise multiplication: %d, %d; %d, %d\n", nrowA, ncolA, nrowB, ncolB);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row_i = 0; row_i < nrowA; row_i++)
	{
		for (int col_i = 0; col_i < ncolA; col_i++)
		{
			res[row_i][col_i] = A[row_i][col_i] / B[row_i][col_i];
		} // col_i loop.
	} // row_i loop.

	return(res);
}

double** matrix_multiply_elementwise(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix)
{
	if (nrowA != nrowB ||
		ncolA != ncolB)
	{
		fprintf(stderr, "Dimensions non-conformant for elementwise multiplication: %d, %d; %d, %d\n", nrowA, ncolA, nrowB, ncolB);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nrowA, ncolA);
	}

	for (int row_i = 0; row_i < nrowA; row_i++)
	{
		for (int col_i = 0; col_i < ncolA; col_i++)
		{
			res[row_i][col_i] = A[row_i][col_i] * B[row_i][col_i];
		} // col_i loop.
	} // row_i loop.

	return(res);
}

void save_matrix_plain_no_headers(double** matrix, int nrow, int ncol, const char* op_fp)
{
	fprintf(stderr, "Saving %dx%d matrix to %s\n", nrow, ncol, op_fp);
	FILE* f_op = open_f(op_fp, "w");

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			double cur_val = matrix[i_row][i_col];

			if (i_col > 0)
			{
				fprintf(f_op, "\t");
			}

			fprintf(f_op, "%.4f", cur_val);
		}

		fprintf(f_op, "\n");
	} // i_row loop.
	close_f(f_op, op_fp);
}

void save_matrix_plain(double** matrix, int nrow, int ncol, const char* op_fp)
{
	fprintf(stderr, "Saving %dx%d matrix to %s\n", nrow, ncol, op_fp);
	FILE* f_op = open_f(op_fp, "w");

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		fprintf(f_op, "Row %d:", i_row);
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			double cur_val = matrix[i_row][i_col];
			fprintf(f_op, "\t%.4f", cur_val);
		}

		fprintf(f_op, "\n");
	} // i_row loop.
	close_f(f_op, op_fp);
}

// Save matrix in text
void save_matrix_binary(double** matrix, int nrow, int ncol, const char* op_fp)
{
	fprintf(stderr, "Saving %dx%d matrix to %s\n", nrow, ncol, op_fp);
	FILE* f_op = open_f(op_fp, "wb");

	// Write the size of matrix, first.
	fwrite(&nrow, sizeof(int), 1, f_op);
	fwrite(&ncol, sizeof(int), 1, f_op);

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			double cur_val = matrix[i_row][i_col];
			fwrite(&cur_val, sizeof(double), 1, f_op);
		}
	} // i_row loop.
	close_f(f_op, op_fp);
}

double** load_matrix_plain(const char* matrix_fp, int& nrow, int& ncol, bool have_row_ids, bool have_col_ids)
{
	vector<char*>* matrix_lines = buffer_file(matrix_fp);

	int start_l_i = (have_col_ids) ? (1):(0);
	int start_col_i = (have_row_ids) ? (1) : (0);

	nrow = (int)matrix_lines->size();
	if (have_col_ids)
	{
		nrow = (int)matrix_lines->size() - 1;
	}

	double** data_matrix = NULL;
	ncol = 0;
	for (int i_l = start_l_i; i_l < (int)matrix_lines->size(); i_l++)
	{
		t_string_tokens* toks = t_string::tokenize_by_chars(matrix_lines->at(i_l), "\t");

		if (ncol == 0)
		{
			ncol = (int)toks->size();
			if (have_row_ids)
			{
				ncol = (int)toks->size() - 1;
			}

			data_matrix = allocate_matrix(nrow, ncol);
		} // ncol check.

		for (int i_t = start_col_i; i_t < (int)toks->size(); i_t++)
		{
			int i_row = i_l - start_l_i;
			int i_col = i_t - start_col_i;
			data_matrix[i_row][i_col] = atof(toks->at(i_t)->str());
		} // i_t loop.

		t_string::clean_tokens(toks);
	} // i_l loop.

	t_string::clean_string_list(matrix_lines);

	return(data_matrix);
} // load_matrix_plain function.

double** load_matrix_binary(const char* matrix_fp, int& nrow, int& ncol, double** res_matrix)
{
	FILE* f_matrix = open_f(matrix_fp, "rb");

	// Write the size of matrix, first.
	int nrow_f = 0;
	int ncol_f = 0;
	//fread(&nrow_f, sizeof(int), 1, f_matrix);
	nrow_f = read_bin_int(f_matrix);
	//fread(&ncol_f, sizeof(int), 1, f_matrix);
	ncol_f = read_bin_int(f_matrix);

	nrow = nrow_f;
	ncol = ncol_f;

	fprintf(stderr, "Loading %dx%d matrix from %s\n", nrow, ncol, matrix_fp);

	double** matrix = res_matrix;
	if (matrix == NULL)
	{
		matrix = allocate_matrix(nrow, ncol);
	}

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			double cur_val;
			//fread(&cur_val, sizeof(double), 1, f_matrix);
			cur_val = read_bin_double(f_matrix);

			matrix[i_row][i_col] = cur_val;
		}
	} // i_row loop.
	close_f(f_matrix, matrix_fp);

	return(matrix);
} // load_matrix_binary option.

void set_matrix_val(double** matrix, int nrow, int ncol, double val)
{
	fprintf(stderr, "Resetting %dx%d matrix to %.3f elementwise.\n", nrow, ncol, val);
	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			matrix[i_row][i_col] = val;
		}
	}
}

double** load_matrix_binary(const char* matrix_fp, int& nrow, int& ncol)
{
	FILE* f_matrix = open_f(matrix_fp, "rb");

	// Write the size of matrix, first.
	int nrow_f = 0;
	int ncol_f = 0;
	//fread(&nrow_f, sizeof(int), 1, f_matrix);
	nrow_f = read_bin_int(f_matrix);
	//fread(&ncol_f, sizeof(int), 1, f_matrix);
	ncol_f = read_bin_int(f_matrix);

	nrow = nrow_f;
	ncol = ncol_f;

	fprintf(stderr, "Loading %dx%d matrix from %s\n", nrow, ncol, matrix_fp);

	double** matrix = allocate_matrix(nrow, ncol);

	for (int i_row = 0; i_row < nrow; i_row++)
	{
		for (int i_col = 0; i_col < ncol; i_col++)
		{
			double cur_val;
			//fread(&cur_val, sizeof(double), 1, f_matrix);
			cur_val = read_bin_double(f_matrix);

			matrix[i_row][i_col] = cur_val;
		}
	} // i_row loop.
	close_f(f_matrix, matrix_fp);

	return(matrix);
} // load_matrix_binary option.

double get_vec2vec_inner_product(double* vec1, double* vec2, int l_vec)
{
	double prod = 0;
	for (int i = 0; i < l_vec; i++)
	{
		prod += (vec1[i] * vec2[i]);
	} // i loop.

	return(prod);
}

double** matrix_multiply(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix)
{
	if (nAcol != nBrow)
	{
		fprintf(stderr, "Dimensions nonconformant: %d, %d ; %d, %d\n", nArow, nAcol, nBrow, nBcol);
		exit(1);
	}

	double** res = res_matrix;
	if (res == NULL)
	{
		res = allocate_matrix(nArow, nBcol);
	}

	for (int row_i = 0; row_i < nArow; row_i++)
	{
		for (int col_i = 0; col_i < nBcol; col_i++)
		{
			double cur_val = 0;
			for (int i_s = 0; i_s < nAcol; i_s++)
			{
				cur_val += A[row_i][i_s] * B[i_s][col_i];
			} // i_s loop.

			res[row_i][col_i] = cur_val;
		} // col_i loop.
	} // row_i loop.

	return(res);
}