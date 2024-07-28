#ifndef __FULL_HISTOGRAM__
#define __FULL_HISTOGRAM__

struct t_full_hist_node
{
	t_full_hist_node** cur_profile_nodes;
	double* cur_profile_counts;
};

class t_full_histogram
{
public:
	t_full_histogram(double** data_profiles, int l_profile, int n_profiles);
	t_full_histogram(char* dist_fp);
	t_full_histogram();
	t_full_histogram(t_full_histogram* full_hist);
	~t_full_histogram();

	void update_histogram_counts_per_data_point(double* data_point, int count_2_update);
	t_full_histogram(int min_per_profile, int max_per_profile, int n_profiles);

	double get_entropy();
	double entropy_per_full_histogram_node(t_full_hist_node* hist_node, int i_dim);

	int* mins_per_dim;
	int* maxes_per_dim;
	
	int n_dims;
	double total_counts;

	t_full_hist_node* hist_node;
	void allocate_hist_node(t_full_hist_node* hist_node, int i_dim);

	t_full_histogram* smooth_full_hist_counts(double self_weight, int n_radius_2_process, int n_iterations);
	void normalize_counts_2_probs();
	void normalize_node_counts_2_probs(t_full_hist_node* hist_node, int i_dim, double total_counts);

	bool increment_vector(double* vector);
	bool increment_vector_per_min_max(double* vector_2_inc, double* min_vector, double* max_vector);
	double get_prob(double* data_vector);
	double& get_count(double* data_vector);
	double& get_count_per_hist_node(t_full_hist_node* hist_node, double total_prob, int i_dim, double* data_vector);
	double get_prob_per_hist_node(t_full_hist_node* hist_node, double total_prob, int i_dim, double* data_vector);
	double get_total_counts();
	bool check_limits_per_vector(double* vector);
	void get_total_counts_per_hist_node(double& cur_count, t_full_hist_node* hist_node, int i_dim);
	void delete_hist_node(t_full_hist_node* hist_node, int i_dim);

	bool compare_int_vectors(int* vec1, int* vec2);

	void reset_counts(double val);
	void reset_counts_per_hist_node(t_full_hist_node* hist_node, int i_dim, double val);

	void load_from_file(char* op_fp);
	void dump_to_file(char* op_fp);
	void load_hist_node(t_full_hist_node* cur_hist_node, int i_dim, FILE* f_op);
	void dump_hist_node(t_full_hist_node* cur_hist_node, int i_dim, FILE* f_op);
	void dump_text(char* op_fp);
	void dump_text_full_hist_node(FILE* f_op, int i_dim, t_full_hist_node* cur_hist_node, int* cur_vals);
};

t_full_histogram* merge_full_histograms(t_full_histogram* hist1, t_full_histogram* hist2);

void merge_full_histogram_nodes(t_full_hist_node* hist_node, t_full_hist_node* hist1_node, 
	t_full_hist_node* hist2_node, 
	int* merged_mins_per_dim, int* merged_maxes_per_dim, 
	int* hist1_mins_per_dim, int* hist1_maxes_per_dim,
	int* hist2_mins_per_dim, int* hist2_maxes_per_dim,
	int n_dims, int i_dim);

t_full_hist_node* copy_hist_node(t_full_hist_node* cur_hist_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims);

#endif // __FULL_HISTOGRAM__