#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "full_histogram.h"
#include "../xmath/log/xlog_math.h"
#include "../geometry/multidimensional_geometry.h"
#include "../file/utils.h"

bool __DUMP_FULL_HISTOGRAM_MSGS__ = false;

t_full_histogram::t_full_histogram()
{}

void t_full_histogram::normalize_counts_2_probs()
{
	double total_counts = this->get_total_counts();

	// Recursively normalize the counts per node.
	this->normalize_node_counts_2_probs(this->hist_node, 0, total_counts);
}

void t_full_histogram::normalize_node_counts_2_probs(t_full_hist_node* cur_hist_node, int i_dim, double log_total_counts)
{
	// Go over all the nodes with counts and normalize the counts to probabilities.
	if(cur_hist_node->cur_profile_counts != NULL)
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_hist_node->cur_profile_counts[i] = xlog_div(cur_hist_node->cur_profile_counts[i], log_total_counts);
		} // i loop.
	}
	else
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			this->normalize_node_counts_2_probs(cur_hist_node->cur_profile_nodes[i], i_dim+1, log_total_counts);
		} // i loop.
	}
} // i loop.

t_full_histogram* t_full_histogram::smooth_full_hist_counts(double self_weight, int n_radius_2_process, int n_iterations)
{
	double* cur_vector = new double[this->n_dims];
	double* min_vector = new double[this->n_dims];
	double* cur_smoothing_vector = new double[this->n_dims];
	double* max_vector = new double[this->n_dims];

if(__DUMP_FULL_HISTOGRAM_MSGS__)
	fprintf(stderr, "Smoothing with self weight %lf, radius %d with %d iterations\n", self_weight, n_radius_2_process, n_iterations);

	t_full_histogram* prev_smoothed_hist = new t_full_histogram(this);
	for(int i_iter = 0; i_iter < n_iterations; i_iter++)
	{
		fprintf(stderr, "Smoothing %d. iteration.\n", i_iter);

		// Copy the previously smoothed histogram.
		t_full_histogram* cur_smoothed_hist = new t_full_histogram(prev_smoothed_hist);
		cur_smoothed_hist->reset_counts(xlog(0.0));

		// Set the current vector.
		for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
		{
			cur_vector[i_dim] = this->mins_per_dim[i_dim];
		} // i_sim

		// Go over all the positions in the histogram.
		while(1)
		{
			// Setup the current min and max vectors around the current vector.
			for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
			{
				min_vector[i_dim] = cur_vector[i_dim] - n_radius_2_process;
				max_vector[i_dim] = cur_vector[i_dim] + n_radius_2_process;

				if(this->mins_per_dim[i_dim] > min_vector[i_dim])
				{
					min_vector[i_dim] = this->mins_per_dim[i_dim];
				}

				if(this->maxes_per_dim[i_dim] < max_vector[i_dim])
				{
					max_vector[i_dim] = this->maxes_per_dim[i_dim];
				}

				cur_smoothing_vector[i_dim] = min_vector[i_dim];
			} // i_dim loop.

			// Count the # of vectors around the current vector.
			double n_total_vecs = 0;
			while(1)
			{
				n_total_vecs++;

				if(!this->increment_vector_per_min_max(cur_smoothing_vector, min_vector, max_vector))
				{
					break;
				}
			} // vector counting loop.

			/*
			fprintf(stderr, "Current vector:\n");
			dump_vector(cur_vector, this->n_dims);
			fprintf(stderr, "%lf total vectors in the smoothing cube.\n", n_total_vecs);
			fprintf(stderr, "Min vector:\n");
			dump_vector(min_vector, this->n_dims);
			fprintf(stderr, "Max vector:\n");
			dump_vector(max_vector, this->n_dims);
			getc(stdin);
			*/

			// Setup the current smoothing vector: This is the vector over which the smoothing cube is traced.
			for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
			{
				cur_smoothing_vector[i_dim] = min_vector[i_dim];
			} // i_dim loop.

			fprintf(stderr, "Current vector:\n");
			dump_vector(cur_vector, this->n_dims);
			fprintf(stderr, "%lf total vectors in the smoothing cube.\n", n_total_vecs);
			fprintf(stderr, "Min vector:\n");
			dump_vector(min_vector, this->n_dims);
			fprintf(stderr, "Max vector:\n");
			dump_vector(max_vector, this->n_dims);
			getc(stdin);

			// Diffuse the counts (cur_count) to all the vectors from min_vector to max_vector.
			double cur_count = prev_smoothed_hist->get_count(cur_vector);
			double log_self_weight = xlog(self_weight);
			double log_non_self_weight = xlog((1 - self_weight) / (n_total_vecs-1));
			while(1)
			{
				// Check which weight we should use: If the smoothing_vector is equal to the cur_vector self_weight is used.
				if(compare_vectors(cur_vector, cur_smoothing_vector, this->n_dims))
				{
					// Add the self weight to the current smoothed count.
					cur_smoothed_hist->get_count(cur_smoothing_vector) = xlog_sum(cur_smoothed_hist->get_count(cur_smoothing_vector), xlog_mul(cur_count, log_self_weight));
				}
				else
				{
					// Add the non-self weight to the current smoothed count.
					cur_smoothed_hist->get_count(cur_smoothing_vector) = xlog_sum(cur_smoothed_hist->get_count(cur_smoothing_vector), xlog_mul(cur_count, log_non_self_weight));
				}

				// Is this the end of the smoothing loops?
				if(!this->increment_vector_per_min_max(cur_smoothing_vector, min_vector, max_vector))
				{
					break;
				}		
			} // Go over all the vectors around the current val.
			
			// Update the vector to move to the next vector.
			if(!this->increment_vector(cur_vector))
			{
				break;
			}
		} // Go over all the positions on the smoothed distribution.

		fprintf(stderr, "Smoothed %d. iteration.\n", i_iter);		

		// Replace the smoothed histograms.
		delete prev_smoothed_hist;
		prev_smoothed_hist = cur_smoothed_hist;
	} // i_iter loop.

	prev_smoothed_hist->total_counts = prev_smoothed_hist->get_total_counts();
	return(prev_smoothed_hist);
}

void t_full_histogram::reset_counts(double val)
{
	this->reset_counts_per_hist_node(this->hist_node, 0, val);
}

void t_full_histogram::reset_counts_per_hist_node(t_full_hist_node* hist_node, int i_dim, double val)
{
	if(hist_node->cur_profile_counts == NULL)
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			this->reset_counts_per_hist_node(hist_node->cur_profile_nodes[i], i_dim+1, val);
		} // i loop.
	}
	else
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			hist_node->cur_profile_counts[i] = val;
		} // i loop.
	}
}

bool t_full_histogram::compare_int_vectors(int* vec1, int* vec2)
{
	for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
	{
		if(vec1[i_dim] != vec2[i_dim])
		{
			return(false);
		}
	} // i_dim loop.

	return(true);
}

bool t_full_histogram::increment_vector_per_min_max(double* vector_2_inc, double* min_vector, double* max_vector)
{
	for(int i_dim = this->n_dims-1; i_dim >= 0; i_dim--)
	{
		if(vector_2_inc[i_dim]+1 <= max_vector[i_dim])
		{
			vector_2_inc[i_dim]++;
			return(true);
		}
		else
		{
			if(i_dim == 0)
			{
				return(false);
			}
			else
			{
				vector_2_inc[i_dim] = min_vector[i_dim];
			}
		}
	} // i_dim loop.

	fprintf(stderr, "We are not supposed to be here.\n");
	exit(1);
	return(false);	
}

bool t_full_histogram::increment_vector(double* vector)
{
	for(int i_dim = this->n_dims-1; i_dim >= 0; i_dim--)
	{
		if(vector[i_dim]+1 <= this->maxes_per_dim[i_dim])
		{
			vector[i_dim]++;
			return(true);
		}
		else
		{
			if(i_dim == 0)
			{
				return(false);
			}
			else
			{
				vector[i_dim] = 0;
			}
		}
	} // i_dim loop.

	fprintf(stderr, "We are not supposed to be here.\n");
	exit(1);
	return(false);
}

t_full_histogram::t_full_histogram(char* dist_fp)
{
	this->load_from_file(dist_fp);
}


double t_full_histogram::get_entropy()
{
	return(entropy_per_full_histogram_node(this->hist_node, 0));
}

double t_full_histogram::entropy_per_full_histogram_node(t_full_hist_node* hist_node, int i_dim)
{
	// Make sure that the counts are normalized to probabilities.
	double total_ent = 0.0;

	if(hist_node == NULL)
	{
		return total_ent;
	}

	double cur_profile_min = this->mins_per_dim[i_dim];
	double cur_profile_max = this->maxes_per_dim[i_dim];

	if(hist_node->cur_profile_counts == NULL)
	{
		for(int i = cur_profile_min; i <= cur_profile_max; i++)
		{
			// Recurse to the next dimension.
			total_ent += entropy_per_full_histogram_node(hist_node->cur_profile_nodes[i], i_dim+1);
		} // i loop.
	}
	else
	{
		for(int i = cur_profile_min; i <= cur_profile_max; i++)
		{
			total_ent += (-1 * hist_node->cur_profile_counts[i] * xexp(hist_node->cur_profile_counts[i]));
		} // i loop.
	}

	return(total_ent);
}

t_full_histogram::t_full_histogram(int min_per_profile, int max_per_profile, int n_profiles)
{
	this->n_dims = n_profiles;
	this->mins_per_dim = new int[n_profiles+1];
	this->maxes_per_dim = new int[n_profiles+1];

	for(int i_prof = 0; i_prof < n_profiles; i_prof++)
	{
		this->mins_per_dim[i_prof] = min_per_profile;
		this->maxes_per_dim[i_prof] = max_per_profile;
	} // i_prof loop.

	// Allocate all the nodes.
	fprintf(stderr, "Allocating histogram nodes.\n");
	this->hist_node = new t_full_hist_node();
	this->hist_node->cur_profile_counts = NULL;
	this->hist_node->cur_profile_nodes = NULL;
	this->allocate_hist_node(this->hist_node, 0);
	fprintf(stderr, "Allocated histogram nodes.\n");
}

void t_full_histogram::update_histogram_counts_per_data_point(double* data_point, int count_2_update)
{
	t_full_hist_node* cur_hist_node = this->hist_node;
	int cur_i_dim = 0;
	while(1)
	{
		int cur_i = (int)(floor(data_point[cur_i_dim]));
		if(cur_hist_node->cur_profile_counts == NULL)
		{
			cur_hist_node = cur_hist_node->cur_profile_nodes[cur_i];
			cur_i_dim++;
		}
		else
		{
			cur_hist_node->cur_profile_counts[cur_i] = xlog_sum(xlog(count_2_update), cur_hist_node->cur_profile_counts[cur_i]);
			break;
		}
	} // Go over all the profile points.
}

// Allocate the full histogram.
t_full_histogram::t_full_histogram(double** data_profiles, int l_profile, int n_profiles)
{
	this->n_dims = n_profiles;
	this->mins_per_dim = new int[n_profiles+1];
	this->maxes_per_dim = new int[n_profiles+1];

	// Set the minimums and maximums per dimension.
	for(int i_prof = 0; i_prof < n_profiles; i_prof++)
	{
		this->mins_per_dim[i_prof] = 1000*1000;
		this->maxes_per_dim[i_prof] = -1000*1000;

		// Profile is 1-indexed.
		for(int i_sig = 1; i_sig <= l_profile; i_sig++)
		{
			if(this->mins_per_dim[i_prof] > data_profiles[i_prof][i_sig])
			{
				this->mins_per_dim[i_prof] = (int)(floor(data_profiles[i_prof][i_sig]));
			}

			if(this->maxes_per_dim[i_prof] < data_profiles[i_prof][i_sig])
			{
				this->maxes_per_dim[i_prof] = (int)(floor(data_profiles[i_prof][i_sig]));
			}
		} // i_sig loop.

if(__DUMP_FULL_HISTOGRAM_MSGS__)
		fprintf(stderr, "Profile %d: %d-%d\n", i_prof, this->mins_per_dim[i_prof], this->maxes_per_dim[i_prof]);
	} // i_prof loop.

	//getc(stdin);

	// Allocate all the nodes.
if(__DUMP_FULL_HISTOGRAM_MSGS__)
	fprintf(stderr, "Allocating histogram nodes.\n");

	this->hist_node = new t_full_hist_node();
	this->hist_node->cur_profile_counts = NULL;
	this->hist_node->cur_profile_nodes = NULL;
	this->allocate_hist_node(this->hist_node, 0);

if(__DUMP_FULL_HISTOGRAM_MSGS__)
	fprintf(stderr, "Allocated histogram nodes.\n");
	//getc(stdin);

	// Update the counts: The profiles are 1 based.
	for(int i_sig = 1; i_sig <= l_profile; i_sig++)
	{
		//for(int i_prof = 0; i_prof < n_profiles; i_prof++)
		//{
		t_full_hist_node* cur_hist_node = this->hist_node;
		int cur_i_dim = 0;
		while(1)
		{
			int cur_i = (int)(floor(data_profiles[cur_i_dim][i_sig]));
			if(cur_hist_node->cur_profile_counts == NULL)
			{
				cur_hist_node = cur_hist_node->cur_profile_nodes[cur_i];
				cur_i_dim++;
			}
			else
			{
				cur_hist_node->cur_profile_counts[cur_i] = xlog_increment(cur_hist_node->cur_profile_counts[cur_i]);
				break;
			}
		}
		//} // i_prof loop.
	} // i_sig loop.

if(__DUMP_FULL_HISTOGRAM_MSGS__)
	fprintf(stderr, "Getting total counts.\n");

	this->total_counts = this->get_total_counts();

if(__DUMP_FULL_HISTOGRAM_MSGS__)
	fprintf(stderr, "Built full histogram with %lf total counts with %lf number of signal vectors of %d dimensions.\n", this->total_counts, xlog(l_profile), n_dims); 
}

void t_full_histogram::allocate_hist_node(t_full_hist_node* cur_hist_node, int i_dim)
{
	if(i_dim == this->n_dims-1)
	{
		cur_hist_node->cur_profile_nodes = NULL;

		cur_hist_node->cur_profile_counts = new double[this->maxes_per_dim[i_dim] - this->mins_per_dim[i_dim] + 1];
		cur_hist_node->cur_profile_counts -= this->mins_per_dim[i_dim];
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_hist_node->cur_profile_counts[i] = xlog(0.0);
		} // i loop.
	}
	else
	{
		cur_hist_node->cur_profile_counts = NULL;
		cur_hist_node->cur_profile_nodes = new t_full_hist_node*[this->maxes_per_dim[i_dim] - this->mins_per_dim[i_dim] + 1];
		cur_hist_node->cur_profile_nodes -= this->mins_per_dim[i_dim];

		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_hist_node->cur_profile_nodes[i] = new t_full_hist_node();
			cur_hist_node->cur_profile_nodes[i]->cur_profile_counts = NULL;
			cur_hist_node->cur_profile_nodes[i]->cur_profile_nodes = NULL;
			allocate_hist_node(cur_hist_node->cur_profile_nodes[i], i_dim+1);
		} // i loop.
	}
}

t_full_histogram::t_full_histogram(t_full_histogram* new_hist)
{
	this->n_dims = new_hist->n_dims;
	this->mins_per_dim = new int[this->n_dims + 1];
	this->maxes_per_dim = new int[this->n_dims + 1];
	for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
	{
		this->mins_per_dim[i_dim] = new_hist->mins_per_dim[i_dim];
		this->maxes_per_dim[i_dim] = new_hist->maxes_per_dim[i_dim];
	} // i loop.
	this->hist_node = copy_hist_node(new_hist->hist_node, this->mins_per_dim, this->maxes_per_dim, 0, this->n_dims);

	this->total_counts = new_hist->total_counts;
}

void t_full_histogram::dump_text(char* op_fp)
{
	FILE* f_op = open_f(op_fp, "w");
	int* cur_vals = new int[this->n_dims+1];
	this->dump_text_full_hist_node(f_op, 0, this->hist_node, cur_vals);
	fclose(f_op);
}

void t_full_histogram::dump_text_full_hist_node(FILE* f_op, int i_dim, t_full_hist_node* cur_hist_node, int* cur_vals)
{
	if(cur_hist_node->cur_profile_counts == NULL)
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_vals[i_dim] = i;
			this->dump_text_full_hist_node(f_op, i_dim+1, cur_hist_node->cur_profile_nodes[i], cur_vals);
		} // i loop.
	}
	else
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_vals[i_dim] = i;

			if(cur_hist_node->cur_profile_counts[i] > xlog(0.0))
			{
				for(int j = 0; j <= i_dim; j++)
				{
					fprintf(f_op, "%d\t", cur_vals[j]);
				} // i loop.

				fprintf(f_op, "%lf\n", cur_hist_node->cur_profile_counts[i]);
			}
		} // i loop.
	}
}

t_full_histogram* merge_full_histograms(t_full_histogram* hist1, t_full_histogram* hist2)
{
	if(hist1->n_dims != hist2->n_dims)
	{
		fprintf(stderr, "The histogram dimensions are not same: %d, %d\n", hist1->n_dims, hist2->n_dims);
		exit(1);
	}

	t_full_histogram* merged_hist = new t_full_histogram();
	merged_hist->n_dims = hist2->n_dims;
	merged_hist->mins_per_dim = new int[merged_hist->n_dims+1];
	merged_hist->maxes_per_dim = new int[merged_hist->n_dims+1];
	for(int i_dim = 0; i_dim < merged_hist->n_dims; i_dim++)
	{
		merged_hist->mins_per_dim[i_dim] = MIN(hist1->mins_per_dim[i_dim], hist2->mins_per_dim[i_dim]);
		merged_hist->maxes_per_dim[i_dim] = MAX(hist1->maxes_per_dim[i_dim], hist2->maxes_per_dim[i_dim]);
	} // i_dim loop.

	merged_hist->hist_node = new t_full_hist_node();
	merge_full_histogram_nodes(merged_hist->hist_node, hist1->hist_node, hist2->hist_node, merged_hist->mins_per_dim, merged_hist->maxes_per_dim,
		hist1->mins_per_dim, hist1->maxes_per_dim, hist2->mins_per_dim, hist2->maxes_per_dim, merged_hist->n_dims, 0);

	merged_hist->total_counts = merged_hist->get_total_counts();

	return(merged_hist);
}

t_full_hist_node* copy_hist_node(t_full_hist_node* cur_hist_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims)
{
	t_full_hist_node* copied_hist_node = new t_full_hist_node();
	copied_hist_node->cur_profile_counts = NULL;
	copied_hist_node->cur_profile_nodes = NULL;
	if(cur_hist_node->cur_profile_counts != NULL)
	{
		copied_hist_node->cur_profile_counts = new double[maxes_per_dim[i_dim] - mins_per_dim[i_dim] + 1];
		copied_hist_node->cur_profile_counts -= mins_per_dim[i_dim];
		for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
		{
			copied_hist_node->cur_profile_counts[i] = cur_hist_node->cur_profile_counts[i];
		} // i loop
	}
	else
	{
		copied_hist_node->cur_profile_nodes = new t_full_hist_node*[maxes_per_dim[i_dim] - mins_per_dim[i_dim] + 1];
		copied_hist_node->cur_profile_nodes -= mins_per_dim[i_dim];
		for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
		{
			copied_hist_node->cur_profile_nodes[i] = copy_hist_node(cur_hist_node->cur_profile_nodes[i], mins_per_dim, maxes_per_dim, i_dim+1, n_dims);
		} // i loop
	}

	return(copied_hist_node);
}

void merge_full_histogram_nodes(t_full_hist_node* hist_node, t_full_hist_node* hist1_node, 
	t_full_hist_node* hist2_node, 
	int* merged_mins_per_dim, int* merged_maxes_per_dim, 
	int* hist1_mins_per_dim, int* hist1_maxes_per_dim,
	int* hist2_mins_per_dim, int* hist2_maxes_per_dim,
	int n_dims, int i_dim)
{
	hist_node->cur_profile_counts = NULL;
	hist_node->cur_profile_nodes = NULL;
	//if(i_dim == n_dims - 1)
	if(hist1_node->cur_profile_counts != NULL)
	{
		hist_node->cur_profile_counts = new double[merged_maxes_per_dim[i_dim] - merged_mins_per_dim[i_dim] + 1];
		hist_node->cur_profile_counts -= merged_mins_per_dim[i_dim];

		for(int i = merged_mins_per_dim[i_dim]; i <= merged_maxes_per_dim[i_dim]; i++)
		{
			double merged_count = xlog(0.0);

			if(i >= hist1_mins_per_dim[i_dim] &&
				i <= hist1_maxes_per_dim[i_dim])
			{
				merged_count = xlog_sum(merged_count, hist1_node->cur_profile_counts[i]);
			}

			if(i >= hist2_mins_per_dim[i_dim] &&
				i <= hist2_maxes_per_dim[i_dim])
			{
				merged_count = xlog_sum(merged_count, hist2_node->cur_profile_counts[i]);
			}

			hist_node->cur_profile_counts[i] = merged_count;
		} // i loop.
	} // count exists check.
	else
	{
		hist_node->cur_profile_nodes = new t_full_hist_node*[merged_maxes_per_dim[i_dim] - merged_mins_per_dim[i_dim] + 1];
		hist_node->cur_profile_nodes -= merged_mins_per_dim[i_dim];

		for(int i = merged_mins_per_dim[i_dim]; i <= merged_maxes_per_dim[i_dim]; i++)
		{
			if(i >= hist1_mins_per_dim[i_dim] &&
				i <= hist1_maxes_per_dim[i_dim] &&
				i >= hist2_mins_per_dim[i_dim] &&
				i <= hist2_maxes_per_dim[i_dim])
			{
				hist_node->cur_profile_nodes[i] = new t_full_hist_node();

				merge_full_histogram_nodes(hist_node->cur_profile_nodes[i], 
					hist1_node->cur_profile_nodes[i], hist2_node->cur_profile_nodes[i], 
					merged_mins_per_dim, merged_maxes_per_dim, 
					hist1_mins_per_dim, hist1_maxes_per_dim,
					hist2_mins_per_dim, hist2_maxes_per_dim,
					n_dims, i_dim+1);
			}
			else if(i >= hist1_mins_per_dim[i_dim] &&
				i <= hist1_maxes_per_dim[i_dim])
			{
				hist_node->cur_profile_nodes[i] = copy_hist_node(hist1_node->cur_profile_nodes[i], hist1_mins_per_dim, hist1_maxes_per_dim, i_dim+1, n_dims);
			}
			else if(i >= hist2_mins_per_dim[i_dim] &&
				i <= hist2_maxes_per_dim[i_dim])
			{			
				hist_node->cur_profile_nodes[i] = copy_hist_node(hist2_node->cur_profile_nodes[i], hist2_mins_per_dim, hist2_maxes_per_dim, i_dim+1, n_dims);
			}
		} // i loop.
	} // Limit check.
}

// Free memory.
t_full_histogram::~t_full_histogram()
{
	this->delete_hist_node(this->hist_node, 0);
}

const int NODE_INDICATOR = 1;
const int COUNT_INDICATOR = 2;
void t_full_histogram::dump_to_file(char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	fwrite(&(this->n_dims), sizeof(int), 1, f_op);
	for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
	{
		fwrite(&(this->mins_per_dim[i_dim]), sizeof(int), 1, f_op);
		fwrite(&(this->maxes_per_dim[i_dim]), sizeof(int), 1, f_op);
	} // i loop.

	// Recursively dump the nodes.
	this->dump_hist_node(this->hist_node, 0, f_op);
	fclose(f_op);
}

void t_full_histogram::dump_hist_node(t_full_hist_node* cur_hist_node, int i_dim, FILE* f_op)
{
	if(cur_hist_node->cur_profile_counts == NULL)
	{
		fwrite(&NODE_INDICATOR, sizeof(int), 1, f_op);
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			this->dump_hist_node(cur_hist_node->cur_profile_nodes[i], i_dim+1, f_op);
		} // i loop.
	}
	else
	{
		fwrite(&COUNT_INDICATOR, sizeof(int), 1, f_op);
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			double cur_count = (cur_hist_node->cur_profile_counts[i]);
			fwrite(&cur_count, sizeof(double), 1, f_op);
		} // i loop.
	}
}

void t_full_histogram::load_hist_node(t_full_hist_node* cur_hist_node, int i_dim, FILE* f_op)
{
	// Load the indocator and switch based on its identity.
	int cur_indicator = 0;
	fread(&cur_indicator, sizeof(int), 1, f_op);

	if(cur_indicator == NODE_INDICATOR)
	{
		cur_hist_node->cur_profile_nodes = new t_full_hist_node*[this->maxes_per_dim[i_dim] - this->mins_per_dim[i_dim] + 1];
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			//this->dump_hist_node(cur_hist_node->cur_profile_nodes[i], i_dim+1, f_op);
			cur_hist_node->cur_profile_nodes[i] = new t_full_hist_node();
			cur_hist_node->cur_profile_nodes[i]->cur_profile_counts = NULL;
			cur_hist_node->cur_profile_nodes[i]->cur_profile_nodes = NULL;
			this->load_hist_node(cur_hist_node->cur_profile_nodes[i], i_dim+1, f_op);
		} // i loop.
	}
	else if(cur_indicator == COUNT_INDICATOR)
	{
		cur_hist_node->cur_profile_counts = new double[this->maxes_per_dim[i_dim] - this->mins_per_dim[i_dim] + 1];
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			//double cur_count = (cur_hist_node->cur_profile_counts[i]);
			//fwrite(&cur_count, sizeof(double), 1, f_op);			
			double cur_count = 0;
			fread(&cur_count, sizeof(double), 1, f_op);
			cur_hist_node->cur_profile_counts[i] = cur_count;
		} // i loop.		
	}
}

bool t_full_histogram::check_limits_per_vector(double* vector)
{
	for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
	{
		if(this->mins_per_dim[i_dim] > vector[i_dim] ||
			this->maxes_per_dim[i_dim] < vector[i_dim])
		{
			return(false);
		}
	} // i_dim loop.

	return(true);
}

void t_full_histogram::load_from_file(char* op_fp)
{
	FILE* f_op = open_f(op_fp, "rb");
	fread(&(this->n_dims), sizeof(int), 1, f_op);
	this->mins_per_dim = new int[this->n_dims+1];
	this->maxes_per_dim = new int[this->n_dims+1];
	for(int i_dim = 0; i_dim < this->n_dims; i_dim++)
	{
		fread(&(this->mins_per_dim[i_dim]), sizeof(int), 1, f_op);
		fread(&(this->maxes_per_dim[i_dim]), sizeof(int), 1, f_op);
	} // i loop.

	this->hist_node = new t_full_hist_node();
	this->hist_node->cur_profile_nodes = NULL;
	this->hist_node->cur_profile_counts = NULL;

	// Recursively dump the nodes.
	this->load_hist_node(this->hist_node, 0, f_op);

	this->total_counts = this->get_total_counts();
	fclose(f_op);
}


void t_full_histogram::delete_hist_node(t_full_hist_node* cur_hist_node, int i_dim)
{
	if(cur_hist_node->cur_profile_nodes != NULL)
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			this->delete_hist_node(cur_hist_node->cur_profile_nodes[i], i_dim+1);

			// This is deleted below at the end.
			//delete cur_hist_node->cur_profile_nodes[i];
		} // i loop.

		cur_hist_node->cur_profile_nodes += this->mins_per_dim[i_dim];
		delete [] cur_hist_node->cur_profile_nodes;
	}

	if(cur_hist_node->cur_profile_counts != NULL)
	{
		cur_hist_node->cur_profile_counts += this->mins_per_dim[i_dim];
		delete [] cur_hist_node->cur_profile_counts;
	}

	delete cur_hist_node;
}

double t_full_histogram::get_prob(double* data_vector)
{
	double current_count = this->get_count(data_vector);
	return(xlog_div(current_count, this->total_counts));
}

double& t_full_histogram::get_count(double* data_vector)
{
	return(this->get_count_per_hist_node(this->hist_node, this->total_counts, 0, data_vector));
}

double& t_full_histogram::get_count_per_hist_node(t_full_hist_node* cur_hist_node, double total_prob, int i_dim, double* data_vector)
{
	int cur_hist_i = (int)(floor(data_vector[i_dim]));
	if(cur_hist_node->cur_profile_counts != NULL)
	{
		return(cur_hist_node->cur_profile_counts[cur_hist_i]);
	}
	else
	{
		return(this->get_count_per_hist_node(cur_hist_node->cur_profile_nodes[cur_hist_i], this->total_counts, i_dim+1, data_vector));
	}
}

double t_full_histogram::get_total_counts()
{
	double total_counts = xlog(0.0);

	get_total_counts_per_hist_node(total_counts, this->hist_node, 0);

	return(total_counts);
}

void t_full_histogram::get_total_counts_per_hist_node(double& cur_count, t_full_hist_node* hist_node, int i_dim)
{
	if(hist_node->cur_profile_counts == NULL)
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			this->get_total_counts_per_hist_node(cur_count, hist_node->cur_profile_nodes[i], i_dim+1);
		} // i loop.
	}
	else
	{
		for(int i = this->mins_per_dim[i_dim]; i <= this->maxes_per_dim[i_dim]; i++)
		{
			cur_count = xlog_sum(cur_count, hist_node->cur_profile_counts[i]);
		} // i loop.
	}
}

