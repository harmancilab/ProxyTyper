#ifndef __GSL_FFT_FILTER_UTILS__
#define __GSL_FFT_FILTER_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;
class t_rng;
struct t_extrema_node;

double* get_extended_odd_length_gaussian(double sigma, double scale, double n_sigma_per_half_win, int& l_filter);
double* get_scaled_odd_length_rectangular(double scale, double sigma, double n_sigma_per_half_win, int& l_filter);
double* get_scaled_odd_length_gaussian(double sigma, double scale, double n_sigma_per_half_win, int& l_filter);
//vector<double*>* multiscale_filter_data(double* cur_real_track_data, int i_t, int l_track_data, double scale_start, double scale_end, double scale_step,
//	vector<double>* scales_per_i_scale);

vector<double*>* multiscale_gaussian_filter_data(double* cur_real_track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool dump_extrema_regions,
	char* op_file_prefix);

double* median_filter_data(double* track_data,
	int l_track_data,
	int l_averaging_win,
	int skip_value); // This is the value to be skipped from the window values.

double* mean_filter_data(double* signal, int l_sig, int l_bin);

int sort_doubles_descending(const void* p1, const void* p2);

vector<double*>* multiscale_median_filter_multiscale_block_permuted_data(double* original_track_data, 
	int l_original_track_data, 
	t_rng* rng,
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	bool return_filtered_tracks,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix);

vector<double*>* recursive_multiscale_median_filter_data(double* _track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	bool return_filtered_tracks,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix,
	double quantile_fraction);

vector<double*>* multiscale_median_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	bool return_filtered_tracks,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix);

vector<double*>* multiscale_median_filter_data_w_imputation(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	bool return_filtered_tracks,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	double skip_value,
	char* op_file_prefix);

vector<double*>* multiscale_avg_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix);

double* median_filter_data(double* track_data,
	int l_track_data, 
	int l_averaging_win,
	int skip_value); // This is the value to be skipped from the window values.

double* mapability_aware_median_filter(double* signal_profile, int l_profile,
	double* scaled_mapability_profile, int l_mapability_profile,
	double max_mapable_signal_2_use_in_filter,
	int l_mapability_filtering_win);

void get_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale);

void get_mapability_aware_median_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale,
	double* normalized_mapability_signal,
	int l_mapability_profile,
	double max_mapability_signal);

vector<double*>* multiscale_conv_filter_data(double* cur_real_track_data, int i_t, int l_track_data, double scale_start, double scale_end, double log_scale_step, vector<double>* scales_per_i_scale);
double* filter_data_per_odd_length_symmetric_filter(double* cur_real_track_data, int l_track_data, double* cur_filter_array, int l_filter, int n_iters);

double* filter_sub_data(double* data, int l_data, int start, int end, 
	double* filter_params, int l_filter,
	int& l_filtered_signal_length, int& i_filtered_data_start);

bool get_origin_of_symmetry(double* filtered_data, int l_signal, int search_start_i, int& origin, int l_mirror_check);

double* sum_tracks(vector<double*>* multitrack_data, int l_signal);


// Non-linear iterative filters:
double* get_multiscale_max_contribution_weight_per_posn(double* signal_profile, int l_profile, int l_start_win, int l_end_win, double log_step);
double* get_contribution_weight_profile(double* signal_profile, int l_profile, int l_gradient_win);
void iterative_contribution_weighted_mean_filter(double* signal_profile, int l_profile, int n_iterations, int l_gradient_win);
void iterative_contribution_weighted_median_filter(double* signal_profile, int l_profile, int n_iterations, int l_gradient_win);
vector<int>* get_local_gradient_abs_extremas(double* signal_profile, int l_profile, int l_win);
vector<int>* get_local_gradient_abs_extremas(vector<t_extrema_node*>* all_extrema, double* signal_profile, int l_profile, int l_win, bool process_maxima);

void get_next_2_exp(int val, int& larger_exp_val, int& expon);

void iterative_gradient_aided_median_filter(double* signal_profile, int l_profile, int n_iterations, int l_gradient_win);
void iterative_gradient_aided_filter(double* signal_profile, int l_profile, int n_iterations, int l_gradient_win);
double* get_per_posn_mean_gradient_profile_per_win(double* signal_profile, int l_profile, int l_gradient_win);
double* get_per_posn_median_gradient_profile_per_win(double* signal_profile, int l_profile, int l_gradient_win);

void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter(double* signal_profile, int l_profile, 
	int n_iterations, 
	double begin_scale, double end_scale, double scale_step);


void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter_fast(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step);

void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter_log_fast(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step);

void iterative_recursive_multiscale_scale_contribution_weighted_median_filter(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step);


double* get_summit_positioned_gradient_profile(double* gradient_profile, int l_profile);

double get_frac_C(double nabla, double kappa);
double get_exp_C(double nabla, double kappa);
void anisotropic_diffusion_filter(double* signal, int l_signal, 
	int n_iterations,
	double delta_t,
	double kappa,
	double (get_C)(double, double));

#endif // __GSL_FFT_FILTER_UTILS__