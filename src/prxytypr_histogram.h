#ifndef __HISTOGRAM__
#define __HISTOGRAM__

#include <vector>
using namespace std;

#define LOG_SMALL_PROB (-20.0)

const double my_PI = 3.14159265359;

class t_rng;

struct t_prob_info
{
	double* counts;
	double* probs;
};

/* 
Histogram node: Holds the  count/probability information and the following histograms per dimension.
*/
struct t_histogram_node
{
	int cur_prof_min;
	int cur_prof_max;

	// The dimension of the data which this node corresponds to.
	int dim;

	// This is the smoothed total counts value coming from the previous dimension.
	double smoothed_total_counts;

	// These are indexed by their respective min's and maxes.
	t_histogram_node** cur_profile_nodes;

	t_prob_info* cur_profile_prob_info;
};

// Transition probability matrix functions.
struct t_transition_matrix
{
	double** state_transition_probs;

	// Counts from which the probabilities are derived from.
	double** state_transition_counts;
};

void get_mode_online(double* signal, int l_signal, double& mode_posn, int& n_vals_at_mode);

double get_correlation_significance(double* arr1, double* arr2, int n_points, int n_rands, t_rng* rng);

double** resample_data_points_separable_gaussian_kernel(double** data, int n_dims, int n_signal_wins, 
	int n_samples_per_sample, t_rng* rng, double gauss_sigma, 
	int& n_total_new_samples);

void delete_transition_matrix(t_transition_matrix* cur_transition_probs, int n_states);
t_transition_matrix* get_transition_matrix_per_transition_counts(double** state_transition_counts, int n_states);
t_transition_matrix* init_transition_matrix(t_transition_matrix* init_state_transition_counts, int n_states);
void copy_transition_matrix(t_transition_matrix* dest_state_transition_matrix, t_transition_matrix* init_state_transition_matrix, int n_states);
void dump_transition_matrix(t_transition_matrix* transition_matrix, int n_states);
void dump_transition_matrix(t_transition_matrix* transition_matrix, int n_states, char* op_fp);
t_transition_matrix* load_transition_matrix(int n_states, char* op_fp);

// Entropy functions.
double jensen_shannon_per_multi_d_histograms(t_histogram_node* hist1, t_histogram_node* hist2);
double kullback_leibler_per_multi_d_histograms(t_histogram_node* hist_p, t_histogram_node* hist_q);
double kullback_leibler_per_histogram_node(t_histogram_node* hist_p, t_histogram_node* hist_q);
double entropy_per_multi_d_histogram(t_histogram_node* hist);
double entropy_per_histogram_node(t_histogram_node* hist_node);
double multivariate_MI_per_multi_d_histograms(t_histogram_node* hist1, t_histogram_node* hist2);
void get_covariance_matrix_per_multi_d_histogram(t_histogram_node* hist);
double* entropy_per_dimension(t_histogram_node* hist_node, int n_dims);

double get_conditional_entropy_online(double** data, int n_dims, 
	int n_signal_wins, 
	int i_dim_2_constrain, double min_val, double max_val, 
	int& n_constrained_pts);

double get_conditional_entropy_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, int i_dim_2_constrain, double val, int& n_constrained_pts);
double get_conditional_entropy_online(double** data, int n_dims, int n_signal_wins, int i_dim_2_constrain, double val, int& n_constraint_pts);
double get_specific_conditional_entropy_online(double** data1, int n_dims1, double** data2, int n_dims2, int n_signal_wins, int constraining_data2_i, int& n_constrained_pts);
double get_entropy_online(double** data, int n_dims, int n_signal_wins);
double get_entropy_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins);
double get_min_probability_online(double** data, int n_dims, int n_signal_wins, double* min_prob_bin_posn);

int map_0_indexed_continuous_data_2_categorical_data_per_min_n_elements_per_bin(double* data, double* mapped_data, int l_sig, int min_n_elements_per_bin);
void map_0_indexed_continuous_data_2_categorical_data_per_n_bins(double* data, double* mapped_data, int l_sig, int n_bins);

double get_cond_prob_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, 
	double* value_2_constrain, 
	double unconstrained_value,
	double* val_2_probe, 
	double& n_probed_values,
	double& n_constrained_values);

double get_prob_online(double** data, int n_dims, int n_signal_wins, double* val_2_probe);
void dump_prob_distribution_online(char* dist_fp, double** data, int n_dims, int n_signal_wins);
double get_prob_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, double* val_2_probe);
void get_bin_posns_per_max_probability(double** data, int n_dims, int n_signal_wins, 
										double max_prob, 
										vector<double*>* max_prob_bin_posns);

void get_bin_posns_per_min_probability(double** data, int n_dims, int n_signal_wins, 
										double max_prob, 
										vector<double*>* max_prob_bin_posns,
										int& n_total_bins);

void get_bin_posns_per_variant_selector_per_min_probability(double** data, int n_dims, int n_signal_wins, 
	bool* variant_selector,
	double min_prob, 
	vector<double*>* min_prob_bin_posns,
	vector<int>* per_posn_count,
	int& n_total_bins);

// Normalization functions.
void normalize_multi_d_hist_probs_per_hist_counts(t_histogram_node* hist);
void normalize_histogram_node_probs_per_counts(t_histogram_node* node, double normalizing_log_factor);
void amplify_histogram_node_counts(t_histogram_node* hist_node, double log_scaling_val);

// Initialize the histogram per data.
t_histogram_node* init_histogram_limits_per_data(double** profile_data, int n_dims, int n_signal_wins);
t_histogram_node* init_multi_d_histogram(t_histogram_node* init_hist);
t_histogram_node* get_multi_d_hist_per_one_indexed_dat(double** profiles, int n_dims, int n_signal_wins);

// Counting the counts, probability accession.
double get_probability(t_histogram_node* hist, int n_dims, double* data);
double get_total_log_counts_per_multi_d_histogram(t_histogram_node* hist);
double get_total_log_counts_per_node(t_histogram_node* node);
void count_data_per_histogram_node_per_one_indexed_data(t_histogram_node* main_node, double** profiles, int n_dims, int n_signal_wins);

double get_total_log_smoothed_counts_per_node(t_histogram_node* node);
void normalize_multi_d_hist_probs_per_smoothed_hist_counts(t_histogram_node* node);

// Merge histograms.
t_histogram_node* merge_multi_d_hist_per_counts(t_histogram_node* hist1, t_histogram_node* hist2);
t_histogram_node* merge_histogram_nodes_per_counts(t_histogram_node* node1, t_histogram_node* node2);

// Loading histogram from binary file.
t_histogram_node* load_multi_d_histogram(char* multi_d_fp);
t_histogram_node* load_histogram_node(FILE* f_hist_node);
t_histogram_node* load_one_d_histogram_node_per_text_linear_histogram(char* hist_fp);

// Sample from a histogram.
double* sample_multi_d_histogram(t_histogram_node* main_node, int n_dims, t_rng* rng);
void sample_histogram_node(t_histogram_node* main_node, int n_dims, t_rng* rng, double* cur_sampled);

// Sample from a list of probabilities.
int sample_index_per_log_probability_array(double* probability_dist, int n_vals, t_rng* rng);

// Memory cleanup
void delete_histogram_node(t_histogram_node* node);
void delete_multi_d_histogram(t_histogram_node* hist);

// Histogram dumping.
void dump_multi_d_histogram(t_histogram_node* main_node, char* fp);
void dump_histogram_node(t_histogram_node* node, FILE* f_hist_node);
//void dump_single_d_histogram_text(t_histogram_node* main_node, char* fp);
void dump_multi_d_histogram_node_text(t_histogram_node* main_node, FILE* f_hist, int* dist_posn_per_dim, int i_last_dim);
void dump_multi_d_histogram_text(t_histogram_node* main_node, int n_dims, char* fp);

// Memory management functions.
double get_total_memory_per_histogram_node(t_histogram_node* main_node);
void prune_histogram_node_memory(t_histogram_node* main_node, int& n_pruned_nodes);

// Extract per dimension histograms: Marginalize on one dimension.
void get_per_dimension_hist_nodes_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, t_histogram_node** hists_per_dim);
void merge_nodes_per_histogram_node_per_dim(t_histogram_node* main_node, int n_dim_2_merge, t_histogram_node* cur_merged_hist_node);
t_histogram_node* get_2D_hist_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, int i_dim1, int i_dim2);
t_histogram_node* get_multi_d_subset_hist_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, int* sorted_i_sub_dims, int i_sub_dim, int n_sub_dims);

// Histogram smoothing functions.
void smooth_multi_d_histogram_counts(t_histogram_node* main_node, bool* is_log_per_dimension, double log_base, int n_dims, double self_weight, int l_win, int n_iters);
void smooth_histogram_node_counts(t_histogram_node* main_node, int i_prof_node, t_histogram_node** main_nodes_per_dims, bool* is_log_per_dimension, int* log_i_2_consecutive_i, double log_base, double self_weight, int l_win, int n_smoothing_iters);

// Log valued array smoothing functions used to smooth the distributions.
void linear_diff_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters);
void log_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters);
void log_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_iters);
void linear_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters);
void linear_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_smoothing_iters);
void linear_diff_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_smoothing_iters);
void linear_diff_smooth_next_log_indexed_count_array(double* count_array, int* i_2_i_mapped, int max_log_val, double log_base, double self_weight, int l_win, int n_smoothing_iters);

double get_log_gaussian_prob(double lin_mean, double lin_std_dev, double lin_val);
double* generate_odd_length_univariate_zero_mean_log_gaussian_dist(double l_bin, int n_pts);

// Convert log and linear indexed histograms to each other. Note that after one goes to log, cannot recover the same linear histogram since the resolution of linear indexed histograms are larger.
t_histogram_node* get_log_indexed_multi_d_histogram_per_linear_indexed_histogram(t_histogram_node* hist_node, double log_base);
t_histogram_node* get_log_indexed_hist_node_per_linear_indexed_hist_node(t_histogram_node* hist_node, double log_base);
t_histogram_node* get_linear_indexed_multi_d_histogram_per_log_indexed_histogram(t_histogram_node* hist_node);

void get_rank_signal(double* prof, int n_pts, double* rank_prof);
void get_correlation(double* prof1, double* prof2, int l_win, double& cur_win_corr);

// Uniform multi-d histogram building for linear and logarithm indexed histograms.
t_histogram_node* get_uniform_multi_d_histogram_per_limits_per_dim(int* mins_per_dim, int* maxes_per_dim, int n_dims);
void get_uniform_histogram_node_per_limits_per_dim(t_histogram_node* cur_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims);

t_histogram_node* get_uniform_multi_d_log_indexed_histogram_per_limits_per_dim(int* mins_per_dim, int* maxes_per_dim, int n_dims, double log_base);
void get_uniform_log_indexed_histogram_node_per_limits_per_dim(t_histogram_node* cur_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims, double log_base);

void reset_counts_per_hist_node(t_histogram_node* hist_node, double log_count_val);

// Log base processing functions.
bool is_this_a_valid_log(double val, double log_base);
void get_log_mapping(double max_log_val, int* log_i_2_consecutive_i, double log_base);

// Simple data statistics.
void get_stats(double* vals, int n_pts, double& mean, double& var);
void get_stats(vector<double>* energies, double& mean, double& std_dev);
double get_mean_per_log_indexed_histogram(t_histogram_node* one_d_hist_node, double log_base);
double get_std_dev_per_log_indexed_histogram(t_histogram_node* one_d_hist_node, double log_base);

double get_mean_per_linear_indexed_histogram(t_histogram_node* one_d_hist_node);
double get_std_dev_per_linear_indexed_histogram(t_histogram_node* one_d_hist_node);

#endif // __HISTOGRAM__