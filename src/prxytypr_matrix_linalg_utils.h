#ifndef __MATRIX_LINALG_UTILS__
#define __MATRIX_LINALG_UTILS__

#include <vector>
using namespace std;

double** allocate_matrix(int nrow, int ncol);
double** allocate_copy_matrix(double** matrix, int nrow, int ncol);
double** allocate_matrix(int nrow, int ncol, double init_val);

void delete_matrix(double** matrix, int n_row, int n_col);

double** process_matrix_elementwise_by_callback(double** matrix, int nrow, int ncol, double(*callback)(double), double** res_matrix);

double** matrix_multiply_elementwise(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix);
double** get_diag_matrix(double** diag_elems, int nrow, int ncol, double** res);

double** matrix_divide_elementwise(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix);
double** matrix_multiply_by_scalar(double** A, int nrowA, int ncolA, double scalar, double** res_matrix);

double** transpose_matrix(double** matrix, int nrow, int ncol, double** res_matrix);

double** get_matrix_per_col_vector(vector<double>* row_vector, double** res_matrix);
double** get_matrix_per_col_vector(double* row_vector, int nrow, double** res_matrix);

double** matrix_row_by_row_inner_product(double** A, int nrowA, int ncolA, double** B, int nrowB, int ncolB, double** res_matrix);
double** matrix_right_multiply_with_diag_matrix(double** A, int nrowA, int ncolA, double** diag, int nrowdiag, int ncoldiag, double** res_matrix);
double** matrix_right_multiply_with_diag_matrix_as_diagonal_only(double** A, int nrowA, int ncolA, double** diag, int nrowdiag, int ncoldiag, double** res_matrix);
double** matrix_left_multiply_with_diag_matrix(double** diag, int nrowdiag, int ncoldiag, double** A, int nrowA, int ncolA, double** res_matrix);
double** matrix_left_multiply_with_diag_matrix_as_diagonal_only(double** diag, int nrowdiag, int ncoldiag, double** A, int nrowA, int ncolA, double** res_matrix);

double** get_matrix_per_row_vector_list(vector<double*>* row_vector_list, int n_col, double** res_matrix);

double get_vec2vec_inner_product(double* vec1, double* vec2, int l_vec);

void print_matrix(double** matrix, int nrow, int ncol);

void copy_matrix(double** src_matrix, int src_nrow, int src_ncol,
	double** dest_matrix, int dest_nrow, int dest_ncol);

double** matrix_add(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix);

void save_matrix_binary(double** matrix, int nrow, int ncol, const char* op_fp);
void save_matrix_plain_no_headers(double** matrix, int nrow, int ncol, const char* op_fp);
void save_matrix_plain(double** matrix, int nrow, int ncol, const char* op_fp);

double** load_matrix_binary(const char* matrix_fp, int& nrow, int& ncol);
double** load_matrix_binary(const char* matrix_fp, int& nrow, int& ncol, double** res_matrix);

void set_matrix_val(double** matrix, int nrow, int ncol, double val);

double** load_matrix_plain(const char* matrix_fp, int& nrow, int& ncol, bool row_ids, bool col_ids);

double** get_submatrix(double** src_matrix, int src_nrow, int src_ncol,
	int nrows_2_copy, int ncols_2_copy, double** res_matrix);

double** invert_matrix_GJ(double** matrix, int nrow, int ncol, double** res_matrix);

double** get_diag_matrix(double* diag_elems, int l_diag, double** res_matrix);

double** get_diagonal_of_matrix(double** square_matrix, int nrow, int ncol, double** res_matrix);

double** matrix_subtract(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix);

double** matrix_multiply(double** A, int nArow, int nAcol,
	double** B, int nBrow, int nBcol, double** res_matrix);

#endif // __MATRIX_LINALG_UTILS__