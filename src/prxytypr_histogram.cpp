#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prxytypr_histogram.h"
#include "prxytypr_annot_region_tools.h"
#include "prxytypr_xlog_math.h"
#include "prxytypr_utils.h"
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"

#include <vector>
#include <algorithm>

using namespace std;

bool __DUMP_HISTOGRAM_MSGS__ = false;

double get_correlation_significance(double* arr1, double* arr2, int n_points, int n_rands, t_rng* rng)
{
	double real_corr = 0;
	get_correlation(arr1, arr2, n_points, real_corr);

	vector<double>* per_rand_corrs = new vector<double>();
	double* rand_arr1 = new double[n_points + 2];
	for (int i_rand = 0; i_rand < n_rands; i_rand++)
	{
		vector<int>* rand_i = rng->fast_permute_indices(0, n_points);
		for (int i = 0; i < n_points; i++)
		{
			rand_arr1[rand_i->at(i)] = arr1[i];
		} // i loop.

		double cur_corr = 0;
		get_correlation(rand_arr1, arr2, n_points, cur_corr);
		per_rand_corrs->push_back(cur_corr);

		delete rand_i;
	} // i_rand loop.

	double mean_corr, stddev_corr;
	get_stats(per_rand_corrs, mean_corr, stddev_corr);
	delete per_rand_corrs;

	double z_score = 0;
	if (stddev_corr > 0)
	{
		z_score = (real_corr - mean_corr) / stddev_corr;
	}

	return(z_score);
}


void get_stats(double* vals, int n_pts, double& mean, double& var)
{
	// Get the mean.
	double total = 0.0;
	for(int i = 0; i < n_pts; i++)
	{
		total += vals[i];
	} // i loop.

	mean = total / n_pts;

	// Get the variance.
	var = 0.0;
	for(int i = 0; i < n_pts; i++)
	{
		var += (vals[i] - mean) * (vals[i] - mean);
	} // i loop.

	var = var / (n_pts - 1);
}

void get_stats(vector<double>* energies, double& mean, double& std_dev)
{
	mean = 0.0;
	for(int i = 0; i < (int)energies->size(); i++)
	{
		mean += energies->at(i);
	} // i loop.

	mean /= (int)energies->size();

	std_dev = 0.0;
	for(int i = 0; i < (int)energies->size(); i++)
	{
		std_dev += (energies->at(i) - mean) * (energies->at(i) - mean);
	} // i loop.

	// Do unbiased computation for sample variance estimate.
	std_dev /= ((int)energies->size() - 1);

	std_dev = pow(std_dev, 0.5);
}

double entropy_per_multi_d_histogram(t_histogram_node* hist_node)
{
	return(entropy_per_histogram_node(hist_node));
}

double get_min_probability_online(double** data, int n_dims, int n_signal_wins, double* min_prob_bin_posn)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	double min_probability = 1.0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;
			if(min_probability > cur_anchor_prob)
			{
				min_probability = cur_anchor_prob;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					min_prob_bin_posn[i_dim] = data[i_dim][i];
				} // i_dim loop.
			} // min_probability update check.
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(min_probability);
}

void get_bin_posns_per_variant_selector_per_min_probability(double** data, int n_dims, int n_signal_wins, 
	bool* variant_selector,
	double min_prob, 
	vector<double*>* min_prob_bin_posns,
	vector<int>* per_posn_count,
	int& n_total_bins)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	n_total_bins = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			n_total_bins++;

			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_bin_posn_cnt = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(variant_selector[i_dim] &&
						data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_bin_posn_cnt++;
				}
			} // j loop.

			if(cur_bin_posn_cnt == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Check if the probability of the current bin position is smaller than the requested value?
			double cur_anchor_prob = cur_bin_posn_cnt/n_signal_wins;
			if(min_prob < cur_anchor_prob)
			{
				double* cur_bin_posn = new double[n_dims];
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					cur_bin_posn[i_dim] = data[i_dim][i];
				} // i_dim loop.

				min_prob_bin_posns->push_back(cur_bin_posn);
				per_posn_count->push_back(cur_bin_posn_cnt);
			} // min_probability update check.
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
} // 

void get_bin_posns_per_min_probability(double** data, int n_dims, int n_signal_wins, 
	double min_prob, 
	vector<double*>* min_prob_bin_posns,
	int& n_total_bins)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	n_total_bins = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			n_total_bins++;

			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_bin_posn_cnt = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_bin_posn_cnt++;
				}
			} // j loop.

			if(cur_bin_posn_cnt == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Check if the probability of the current bin position is smaller than the requested value?
			double cur_anchor_prob = cur_bin_posn_cnt/n_signal_wins;
			if(min_prob < cur_anchor_prob)
			{
				double* cur_bin_posn = new double[n_dims];
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					cur_bin_posn[i_dim] = data[i_dim][i];
				} // i_dim loop.

				min_prob_bin_posns->push_back(cur_bin_posn);
			} // min_probability update check.
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
} // get_bin_posns_per_min_probability

// Return the bin posns with maximum probability requirement.
void get_bin_posns_per_max_probability(double** data, int n_dims, int n_signal_wins, 
	double max_prob, 
	vector<double*>* max_prob_bin_posns)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_bin_posn_cnt = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_bin_posn_cnt++;
				}
			} // j loop.

			if(cur_bin_posn_cnt == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Check if the probability of the current bin position is smaller than the requested value?
			double cur_anchor_prob = cur_bin_posn_cnt/n_signal_wins;
			if(max_prob > cur_anchor_prob)
			{
				double* cur_bin_posn = new double[n_dims];
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					cur_bin_posn[i_dim] = data[i_dim][i];
				} // i_dim loop.

				max_prob_bin_posns->push_back(cur_bin_posn);
			} // min_probability update check.
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
} // get_bin_posns_per_max_probability

//double** quantize_data(double** data, int n_dims, int l_data, int n_bins)
//{
//	double** quantized_data = new double*[n_dims];
//	double* mins_per_dim = new double[n_dims];
//	double* maxes_per_dim = new double[n_dims];
//	for(int i_dim = 0; i_dim < n_dims; i_dim++)
//	{
//		quantized_data[i_dim] = new double[l_data+2];
//		memset(quantized_data[i_dim], 0, sizeof(double)*(l_data+2));
//		mins_per_dim[i_dim]
//
//		for(int i = 1; i <= l_data; i++)
//		{
//		} // i loop.
//	} // i_dim loop.
//}

double** resample_data_points_separable_gaussian_kernel(double** data, int n_dims, int n_signal_wins, 
	int n_samples_per_sample, t_rng* rng, double gauss_sigma, 
	int& n_total_new_samples)
{
	n_total_new_samples = n_samples_per_sample * n_signal_wins + n_signal_wins + 2;
	fprintf(stderr, "Resampling from a dataset of %d points with %d samples per sample.\n", n_signal_wins, n_samples_per_sample);

	double** resampled_data = new double*[n_dims];
	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		resampled_data[i_dim] = new double[n_total_new_samples];

		// Copy the data.
		for(int i_samp = 1; i_samp <= n_signal_wins; i_samp++)
		{
			resampled_data[i_dim][i_samp] = data[i_dim][i_samp];
		} // i_dim loop.
	} // i_dim loop.

	// Go over all the samples and generate points around them.
	int cur_i = n_signal_wins+1;
	for(int i_samp = 1; i_samp <= n_signal_wins; i_samp++)
	{
		fprintf(stderr, "Resampling %d. sample\r", i_samp);
		for(int i_new_sample = 0; i_new_sample < n_samples_per_sample; i_new_sample++)
		{
			for(int i_dim = 0; i_dim < n_dims; i_dim++)
			{
				// Box-Muller transform to generate Normal distributed data.
				double U1 = rng->random_double_ran3();
				double U2 = rng->random_double_ran3();
				double cur_gaussian_rand = sqrt(-2 * log(U1)) * cos(2*my_PI*U2) * gauss_sigma;

				// Generate the sample around the current sample.
				resampled_data[i_dim][cur_i] = data[i_dim][i_samp] + cur_gaussian_rand;
			} // i_dim loop.

			cur_i++;
		} // i_new_sample loop.
	} // i_samp loop.

	n_total_new_samples = cur_i-1;
	
	return(resampled_data);
}

double get_prob_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, double* val_2_probe)
{
	double n_matches = 0;

	for(int i = 0; i < n_signal_wins; i++)
	{		
		// If the current entry is equal to the anchor value, increment the count.
		bool cur_entry_matches = true;
		for(int i_dim = 0; i_dim < n_dims; i_dim++)
		{
			if(data[i_dim][i] != val_2_probe[i_dim])
			{
				cur_entry_matches = false;
				break;
			}
		} // i_dim comparison loop between ith and jth entries.

		if(cur_entry_matches)
		{
			n_matches++;
		}
	} // i loop.

	return(n_matches / n_signal_wins);
}

double get_prob_per_0_indexed_data_online2(double** data, int n_dims, int n_signal_wins, double* val_2_probe)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	for(int i = 0; i < n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j < n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;

			// Is the current bin position equal to the data?
			bool cur_entry_matches = true;
			for(int i_dim = 0; i_dim < n_dims; i_dim++)
			{
				if(data[i_dim][i] != val_2_probe[i_dim])
				{
					cur_entry_matches = false;
				}
			} // i_dim comparison between the value and bin position.

			if(cur_entry_matches)
			{
				delete [] processed_state;				
				return(cur_anchor_prob);
			}
		} // state processing check.
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(0.0);
}

double get_cond_prob_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, 
	double* value_2_constrain, 
	double unconstrained_value,
	double* val_2_probe, 
	double& n_probed_values,
	double& n_constrained_values)
{
	for(int i = 0; i < n_signal_wins; i++)
	{
		// Check first whether this value satisfies the condition.
		bool cur_value_satisfies = true;
		for(int i_dim = 0; i_dim < n_dims; i_dim++)
		{
			if(value_2_constrain[i_dim] != unconstrained_value &&
				floor(data[i_dim][i]) != floor(value_2_constrain[i_dim]))
			{
				cur_value_satisfies = false;
				break;
			}
		} // i_dim loop.

		if(cur_value_satisfies)
		{
			n_constrained_values++;

			bool cur_val_matches = true;
			for(int i_dim = 0; i_dim < n_dims; i_dim++)
			{
				if(floor(data[i_dim][i]) != floor(val_2_probe[i_dim]))
				{
					cur_val_matches = false;
					break;
				}
			} // i_dim loop.

			if(cur_val_matches)
			{
				n_probed_values++;
			} // value matches check.
		} // value satisfies check.
	} // i loop.

	if(n_constrained_values == 0)
	{
		return(0);
	}
	else
	{
		return(n_probed_values / n_constrained_values);
	}
}

double get_cond_prob_per_0_indexed_data_online2(double** data, int n_dims, int n_signal_wins, 
	double* value_2_constrain, 
	double unconstrained_value,
	double* val_2_probe, 
	double& n_probed_values,
	double& n_constrained_values)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	bool* cond_satistfying_state = new bool[n_signal_wins + 2];
	memset(cond_satistfying_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over all the values and set the condition satisfying state and count the number of condition satisfying values.
	n_constrained_values = 0;
	for(int i = 0; i < n_signal_wins; i++)
	{
		bool cur_value_satisfies = true;
		for(int i_dim = 0; i_dim < n_dims; i_dim++)
		{
			// If there is unconstrained value in the constraining vector or the value matches, we still satisfy the constraint.
			if(value_2_constrain[i_dim] != unconstrained_value &&
				floor(data[i_dim][i]) != floor(value_2_constrain[i_dim]))
			{
				cur_value_satisfies = false;
				break;
			}
		} // i_dim loop.

		if(cur_value_satisfies)
		{
			n_constrained_values++;
		}

		cond_satistfying_state[i] = cur_value_satisfies;
	} // i loop.	

	// Go over each value and process it.
	for(int i = 0; i < n_signal_wins; i++)
	{
		// Make sure that i^th signal satisfies constraint.
		if(!cond_satistfying_state[i])
		{
			continue;
		}

		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			n_probed_values = 0;
			for(int j = i; j < n_signal_wins; j++)
			{
				if(!cond_satistfying_state[j])
				{
					continue;
				}

				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(floor(data[i_dim][i]) != floor(data[i_dim][j]))
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					n_probed_values++;
				}
			} // j loop.

			if(n_probed_values == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = n_probed_values/n_constrained_values;

			// Is the current bin position equal to the data that we are probing?
			bool cur_entry_matches = true;
			for(int i_dim = 0; i_dim < n_dims; i_dim++)
			{
				if(floor(data[i_dim][i]) != floor(val_2_probe[i_dim]))
				{
					cur_entry_matches = false;
				}
			} // i_dim comparison between the value and bin position.

			if(cur_entry_matches)
			{
				delete [] cond_satistfying_state;
				delete [] processed_state;				
				return(cur_anchor_prob);
			}
		} // state processing check.
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] cond_satistfying_state;
	delete [] processed_state;
	return(0.0);
}

void dump_prob_distribution_online(char* dist_fp, double** data, int n_dims, int n_signal_wins)
{
	FILE* f_dist = open_f(dist_fp, "w");

	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		for(int i_samp = 1; i_samp <= n_signal_wins; i_samp++)
		{
			data[i_dim][i_samp] = floor(data[i_dim][i_samp]);
		} // i_samp loop.
	} // i_dim loop.

	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;

			for(int dim_i = 0; dim_i < n_dims; dim_i++)
			{
				fprintf(f_dist, "%.3f\t", data[dim_i][i]);
			} // dim_i loop.

			fprintf(f_dist, "%lf\t%.1f\n", cur_anchor_prob, cur_anchor_count);
		} // state processing check.
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	fclose(f_dist);
}

double get_prob_online(double** data, int n_dims, int n_signal_wins, double* val_2_probe)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;

			// Is the current bin position equal to the data?
			bool cur_entry_matches = true;
			for(int i_dim = 0; i_dim < n_dims; i_dim++)
			{
				if(data[i_dim][i] != val_2_probe[i_dim])
				{
					cur_entry_matches = false;
				}
			} // i_dim comparison between the value and bin position.

			if(cur_entry_matches)
			{
				delete [] processed_state;				
				return(cur_anchor_prob);
			}
		} // state processing check.
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(0.0);
}

// Following computes the specific conditional entropy 
double get_specific_conditional_entropy_online(double** data1, int n_dims1, double** data2, int n_dims2, int n_signal_wins, int constraining_data2_i, int& n_constrained_pts)
{
	bool* data2_matching_state = new bool[n_signal_wins + 2];
	memset(data2_matching_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over all the data2 positions and mark the positions that match.
	n_constrained_pts = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		bool current_posn_matches = true;
		for(int i_dim = 0; i_dim < n_dims2; i_dim++)
		{
			if(data2[i_dim][i] != data2[i_dim][constraining_data2_i])
			{
				current_posn_matches = false;
				break;
			}
		} // i_dim loop.

		data2_matching_state[i] = current_posn_matches;

		if(current_posn_matches)
		{
			n_constrained_pts++;
		}
	} // i loop.

	bool* data1_processed_state = new bool[n_signal_wins + 2];
	memset(data1_processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Go over each value and process it.
	double cond_entropy = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		// Make sure that this value is not processed and that the entry at this position satisfies contraint.
		if(data2_matching_state[i] &&
			!data1_processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims1; i_dim++)
				{
					if(data1[i_dim][i] != data1[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(data2_matching_state[j] && 
					cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					data1_processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_constrained_pts;
			cond_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed/does not match the constraint.
		}
	} // i loop.

	delete [] data2_matching_state;
	delete [] data1_processed_state;
	return(cond_entropy);
}

double get_conditional_entropy_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins, int i_dim_2_constrain, double val, int& n_constrained_pts)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Count the total number of entries that satisfy the contraint.
	n_constrained_pts = 0;
	for(int i = 0; i < n_signal_wins; i++)
	{
		if(data[i_dim_2_constrain][i] == val)
		{
			n_constrained_pts++;
		}
	} // i loop.

	// Go over each value and process it.
	double cond_entropy = 0;
	for(int i = 0; i < n_signal_wins; i++)
	{
		// Make sure that this value is not processed and that the entry at this position satisfies contraint.
		if(!processed_state[i] &&
			data[i_dim_2_constrain][i] == val)
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j < n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_constrained_pts;
			cond_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(cond_entropy);
}


// Compute the conditional entropy, or constrained entropy.
double get_conditional_entropy_online(double** data, int n_dims, int n_signal_wins, int i_dim_2_constrain, double val, int& n_constrained_pts)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Count the total number of entries that satisfy the contraint.
	n_constrained_pts = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(data[i_dim_2_constrain][i] == val)
		{
			n_constrained_pts++;
		}
	} // i loop.

	// Go over each value and process it.
	double cond_entropy = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		// Make sure that this value is not processed and that the entry at this position satisfies contraint.
		if(!processed_state[i] &&
			data[i_dim_2_constrain][i] == val)
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_constrained_pts;
			cond_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(cond_entropy);
}

double get_conditional_entropy_online(double** data, int n_dims, int n_signal_wins, int i_dim_2_constrain, double min_val, double max_val, int& n_constrained_pts)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	// Count the total number of entries that satisfy the contraint.
	n_constrained_pts = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(data[i_dim_2_constrain][i] >= min_val &&
			data[i_dim_2_constrain][i] <= max_val)
		{
			n_constrained_pts++;
		}
	} // i loop.

	// Go over each value and process it.
	double cond_entropy = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		// Make sure that this value is not processed and that the entry at this position satisfies contraint.
		if(!processed_state[i] &&
			data[i_dim_2_constrain][i] >= min_val &&
			data[i_dim_2_constrain][i] <= max_val)
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_constrained_pts;
			cond_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(cond_entropy);
}

double get_entropy_online(double** data, int n_dims, int n_signal_wins)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		for(int i_samp = 1; i_samp <= n_signal_wins; i_samp++)
		{
			data[i_dim][i_samp] = floor(data[i_dim][i_samp]);
		} // i_samp loop.
	} // i_dim loop.

	// Go over each value and process it.
	double total_entropy = 0;
	for(int i = 1; i <= n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j <= n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;
			total_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(total_entropy);
} // get_entropy_online

int map_0_indexed_continuous_data_2_categorical_data_per_min_n_elements_per_bin(double* data, double* mapped_data, int l_sig, int min_n_elements_per_bin)
{
	const int max_n_bins = 5000;
	int n_elements_per_bins[max_n_bins];
	int n_bins = 2;
	while(1)
	{
		memset(n_elements_per_bins, 0, max_n_bins*sizeof(int));

		map_0_indexed_continuous_data_2_categorical_data_per_n_bins(data, mapped_data, l_sig, n_bins);

		// for each bin index, count the # of elements.
		//int n_elements_per_cur_bin = 0;
		for(int i = 0; i < l_sig; i++)
		{
			n_elements_per_bins[(int)(mapped_data[i])]++;
		} // i_bin loop

		bool cur_binning_valid = true;
		for(int i_bin = 0; i_bin < n_bins; i_bin++)
		{
			if(n_elements_per_bins[i_bin] < min_n_elements_per_bin)
			{
				cur_binning_valid = false;
				break;
			}
		} // i_bin loop.

		if(!cur_binning_valid)
		{
			map_0_indexed_continuous_data_2_categorical_data_per_n_bins(data, mapped_data, l_sig, n_bins-1);
			return(n_bins-1);
		}
		else
		{
			n_bins++;
		}
	} // bin size counting loop.
}

void map_0_indexed_continuous_data_2_categorical_data_per_n_bins(double* data, double* mapped_data, int l_sig, int n_bins)
{
	double min_val = 10000;
	double max_val = -10000;

	for(int i = 0; i < l_sig; i++)
	{
		if(data[i] < min_val)
		{
			min_val = data[i];
		}

		if(data[i] > max_val)
		{
			max_val = data[i];
		}
	} // i loop.

	double bin_size = (max_val - min_val) / n_bins;

	for(int i = 0; i < l_sig; i++)
	{
		// Get the bin index for current value.
		double bin_i = (floor((data[i] - min_val) / bin_size));

		if(bin_i == n_bins)
		{
			bin_i = n_bins-1;
		}

		if(mapped_data == NULL)
		{
			data[i] = bin_i;
		}
		else
		{
			mapped_data[i] = bin_i;
		}
	} // i loop
}

double get_entropy_per_0_indexed_data_online(double** data, int n_dims, int n_signal_wins)
{
	bool* processed_state = new bool[n_signal_wins + 2];
	memset(processed_state, 0, sizeof(bool) * (n_signal_wins + 2));

	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		for(int i_samp = 0; i_samp < n_signal_wins; i_samp++)
		{
			data[i_dim][i_samp] = floor(data[i_dim][i_samp]);
		} // i_samp loop.
	} // i_dim loop.

	// Go over each value and process it.
	double total_entropy = 0;
	for(int i = 0; i < n_signal_wins; i++)
	{
		if(!processed_state[i])
		{
			// Set this as the anchor value, go over all the values to the right of the track and count the number of times the entries match this value.
			double cur_anchor_count = 0;
			for(int j = i; j < n_signal_wins; j++)
			{
				// If the current entry is equal to the anchor value, increment the count.
				bool cur_entry_matches = true;
				for(int i_dim = 0; i_dim < n_dims; i_dim++)
				{
					if(data[i_dim][i] != data[i_dim][j])
					{
						cur_entry_matches = false;
						break;
					}
				} // i_dim comparison loop between ith and jth entries.

				if(cur_entry_matches)
				{
					// Set the current entry as processed so it is not re-counted.
					processed_state[j] = true;
					cur_anchor_count++;
				}
			} // j loop.

			if(cur_anchor_count == 0)
			{
				fprintf(stderr, "Sanity check failed: anchor count is 0 for %d. anchor.\n", i);
				exit(1);
			}

			// We have the anchor count. Update the entropy.
			double cur_anchor_prob = cur_anchor_count/n_signal_wins;
			total_entropy += -1 * cur_anchor_prob * xlog(cur_anchor_prob);
		}
		else
		{
			// Skip this entry since it is processed.
		}
	} // i loop.

	delete [] processed_state;
	return(total_entropy);
} // get_entropy_per_0_indexed_data_online


double* entropy_per_dimension(t_histogram_node* hist_node, int n_dims)
{
	t_histogram_node** hist_nodes_per_dim = new t_histogram_node*[n_dims];
	get_per_dimension_hist_nodes_per_multi_d_hist_node(hist_node, n_dims, hist_nodes_per_dim);

	double* entropies_per_dim = new double[n_dims];

	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		normalize_multi_d_hist_probs_per_hist_counts(hist_nodes_per_dim[i_dim]);
		entropies_per_dim[i_dim] = entropy_per_multi_d_histogram(hist_nodes_per_dim[i_dim]);
		delete_multi_d_histogram(hist_nodes_per_dim[i_dim]);
	} // i_dim loop.

	// Return the entropies per dimensions.
	return(entropies_per_dim);
}

double entropy_per_histogram_node(t_histogram_node* hist_node)
{
	double total_ent = 0.0;

	if(hist_node == NULL)
	{
		return total_ent;
	}

	if(hist_node->cur_profile_prob_info == NULL)
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			total_ent += entropy_per_histogram_node(hist_node->cur_profile_nodes[i]);
		} // i loop.
	}
	else
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			total_ent += (-1 * hist_node->cur_profile_prob_info->probs[i] * xexp(hist_node->cur_profile_prob_info->probs[i]));
		} // i loop.
	}

	return(total_ent);
}

// Naming conventions: http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
double kullback_leibler_per_multi_d_histograms(t_histogram_node* hist_p, t_histogram_node* hist_q)
{
	return(kullback_leibler_per_histogram_node(hist_p, hist_q));
}

double kullback_leibler_per_histogram_node(t_histogram_node* hist_p, t_histogram_node* hist_q)
{
	double total_KL = 0.0;

	if((hist_p == NULL && hist_q != NULL) ||
		(hist_p == NULL && hist_q != NULL))
	{
		fprintf(stderr, "One of the nodes is NULL in KL computation.\n");
		exit(1);
	}
	else if(hist_p == NULL)
	{
		return(total_KL);
	}

	if(hist_p->cur_profile_prob_info == NULL)
	{
		for(int i = hist_p->cur_prof_min; i <= hist_p->cur_prof_max; i++)
		{
			total_KL += kullback_leibler_per_histogram_node(hist_p->cur_profile_nodes[i], hist_q->cur_profile_nodes[i]);
		} // i loop.
	}
	else
	{
		for(int i = hist_p->cur_prof_min; i <= hist_p->cur_prof_max; i++)
		{
			if(hist_p->cur_profile_prob_info->probs[i] == xlog(0.0) &&
				hist_q->cur_profile_prob_info->probs[i] == xlog(0.0))
			{
			}
			else if(hist_q->cur_profile_prob_info->probs[i] == xlog(0.0))
			{
				fprintf(stderr, "Q is 0 @ %s(%d).\n", __FILE__, __LINE__);
				exit(1);
			}
			else
			{
				total_KL += xlog_div(hist_p->cur_profile_prob_info->probs[i], hist_q->cur_profile_prob_info->probs[i]) * xexp(hist_p->cur_profile_prob_info->probs[i]);
			}
		} // i loop.
	}

	return(total_KL);
}

// Naming conventions: http://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
double jensen_shannon_per_multi_d_histograms(t_histogram_node* hist_p, t_histogram_node* hist_q)
{
	// Get the total counts for both histograms.
	double log_hist_p_total = get_total_log_counts_per_multi_d_histogram(hist_p);
	double log_hist_q_total = get_total_log_counts_per_multi_d_histogram(hist_q);

	double log_total_count = xlog_sum(log_hist_p_total, log_hist_q_total);

	// Amplify the histograms to match the counts to the total count.
	double q_log_scaling = xlog_div(log_total_count, log_hist_q_total);
	double p_log_scaling = xlog_div(log_total_count, log_hist_p_total);
	amplify_histogram_node_counts(hist_q, q_log_scaling);
	amplify_histogram_node_counts(hist_p, p_log_scaling);

	// Normalize the histograms.
	normalize_multi_d_hist_probs_per_hist_counts(hist_p);
	normalize_multi_d_hist_probs_per_hist_counts(hist_q);

	// Get the merged histogram from the scaled histograms.
	t_histogram_node* merged_pq_dist = merge_multi_d_hist_per_counts(hist_p, hist_q);	

	// Normalize the merged histogram.
	normalize_multi_d_hist_probs_per_hist_counts(merged_pq_dist);

	// Compute the KL divergence for symmetric cases.
	double JS_distance = .5 * kullback_leibler_per_histogram_node(hist_p, merged_pq_dist) + 
						.5 * kullback_leibler_per_histogram_node(hist_q, merged_pq_dist);

	return(JS_distance);
}

void amplify_histogram_node_counts(t_histogram_node* hist_node, double log_scaling_val)
{
	if(hist_node == NULL)
	{
		return;
	}

	if(hist_node->cur_profile_prob_info == NULL)
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			amplify_histogram_node_counts(hist_node->cur_profile_nodes[i], log_scaling_val);
		} // i loop.
	}
	else
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			hist_node->cur_profile_prob_info->counts[i] = xlog_mul(hist_node->cur_profile_prob_info->counts[i], log_scaling_val);
		} // i loop.
	}
}

// Get and set the histogram with limits but do not update the counts of the histogram.
// This makes sure that the histogram has the limits for the data.
t_histogram_node* init_histogram_limits_per_data(double** profiles, int n_dims, int n_signal_wins)
{
	// This is the histogram node corresponding to the first profile.
	t_histogram_node* main_node = new t_histogram_node();
	main_node->cur_prof_min = 1000*1000;
	main_node->cur_prof_max = -1000*1000;
	main_node->dim = 0; // The dimension corresponds to 

	if (__DUMP_HISTOGRAM_MSGS__)
	{
		fprintf(stderr, "Setup the histogram for the main dimension.\n");
	}

/*
	Following loop sets the limits for the histograms for each profile. The main node is the root node for the histogram. 
	For each dimension and for all the data points, first the node corresponding to 
*/
/*
	Following loop sets the next profile limits for each node. For each position in each dimension, it backtracks the
	node to the (next_i_dim-1)st dimension, then sets the limits for that node. next_i_dim indexes the profile for which the 
	mins and maxes are searched for. The (next_i_dim-1)st profile's mins and maxes are set by this search.
*/
	for(int cur_i_dim = 0; cur_i_dim < n_dims; cur_i_dim++)
	{
		if (__DUMP_HISTOGRAM_MSGS__)
		{
			fprintf(stderr, "Setting the next level limits for dimension %d\n", cur_i_dim);
		}

		// Set the limits for the previous dimension histograms.
		for(int i_sig = 1; i_sig <= n_signal_wins; i_sig++)
		{
			// Starting from the nodes for the main histogram, trace this position.
			//t_histogram_node* cur_node = main_node->next_level_nodes[(int)(profiles[0][i_sig])];
			t_histogram_node* cur_node = main_node;
			int traced_i_dim = 0;
			while(1) 
			{
				if(cur_node == NULL)
				{
					fprintf(stderr, "Fatal error 7.\n");
					fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
					exit(1);
				}

				// When the traced dimension reaches i_dim-1, we can exit and set the min and max.
				if(traced_i_dim == cur_i_dim)
				{
					break;
				}

				if(cur_node == NULL || cur_node->cur_profile_nodes == NULL)
				{
					fprintf(stderr, "Fatal error 1: %d, %d, %d, %d, %lf\n", traced_i_dim, i_sig, cur_node->cur_prof_min, cur_node->cur_prof_max,  profiles[traced_i_dim][i_sig]);
					fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
					exit(1);
					break;
				}
				else if(cur_node->cur_prof_max < cur_node->cur_prof_min)
				{
					// The minimum and maximum for the previous dimensions must have been set, if not, there is a problem.
					fprintf(stderr, "Fatal error 5: %d, %d\n", cur_i_dim, traced_i_dim);
					fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
					exit(1);
				}
				else
				{
					// Get the next node: Next value is taken from the next traced dim.
                    if(cur_node->cur_prof_min > (int)profiles[traced_i_dim][i_sig] ||
                            cur_node->cur_prof_max < (int)profiles[traced_i_dim][i_sig])
                    {
						fprintf(stderr, "Fatal error 4: %d, %d, %d, %d, %d, %lf\n", traced_i_dim, cur_i_dim, i_sig, cur_node->cur_prof_min, 
							cur_node->cur_prof_max, profiles[traced_i_dim][i_sig]);	
						fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
						exit(1);
					}

					// Move to the 
					cur_node = cur_node->cur_profile_nodes[(int)(profiles[traced_i_dim][i_sig])];

					if(cur_node == NULL)
					{
						fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
						exit(1);
					}
				}

				traced_i_dim++;
			} // dimension tracing loop.
			
			// Found the node that corresponds to the current dimension.
			if(cur_node->cur_prof_min > (int)(profiles[cur_i_dim][i_sig]))
			{
				cur_node->cur_prof_min = (int)(profiles[cur_i_dim][i_sig]);
			}

			if(cur_node->cur_prof_max < (int)(profiles[cur_i_dim][i_sig]))
			{
				cur_node->cur_prof_max = (int)(profiles[cur_i_dim][i_sig]);
			}
		} // i_sig loop.

		if (__DUMP_HISTOGRAM_MSGS__)
		{
			fprintf(stderr, "Limits are set for profile %d.\n", cur_i_dim);
		}
		
		//getc(stdin);
/*
		for(int i = main_dim_hist->min_val; i <= main_dim_hist->max_val; i++)
		{
			fprintf(stderr, "%d: %d-%d\n", i, main_dim_hist->nodes[i]->cur_prof_min,  main_dim_hist->nodes[i]->cur_prof_max);
		}
		exit(1);
*/
		// All the min-maxes for the current dimension are set, allocate the next node lists for the current dimension.
		for(int i_sig = 1; i_sig <= n_signal_wins; i_sig++)
		{
			// For this position, get the node at the previous dimension.
			t_histogram_node* cur_node = main_node;
			int traced_i_dim = 0;
			while(1)
			{
				// Are we at the node just at the previous dimension compared to the node we are looking for?
				if(traced_i_dim == cur_i_dim)
				{
					break;
				}

				if(cur_node->cur_profile_nodes == NULL)
				{
					fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
					exit(1);
					break;
				}
				else
				{
					// Get the next node.
					if(cur_node->cur_prof_min > (int)profiles[traced_i_dim][i_sig] || 
						cur_node->cur_prof_max < (int)profiles[traced_i_dim][i_sig])	
					{
						fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
						exit(1);
					}
					cur_node = cur_node->cur_profile_nodes[(int)(profiles[traced_i_dim][i_sig])];
				}

				// This is the dimension where the current node is at. Must find the dimension just above the current dimension of data processing.
				traced_i_dim++;
			} // dimension tracing loop.

			// Set up the next level nodes list for the node in the previous dimension.
			if(cur_node->cur_profile_nodes != NULL)
			{
				// This node already has next level nodes setup, no need to allocate it again. This happens for the nodes with multiple entries.
			}
			else
			{
				// If the current node has valid limits, allocate the nodes.
				if(cur_node->cur_prof_max >= cur_node->cur_prof_min)
				{
					//fprintf(stderr, "Setting next level nodes for dimension %d.\n", traced_i_dim);
					cur_node->cur_profile_nodes = new t_histogram_node*[cur_node->cur_prof_max - cur_node->cur_prof_min +1];
					cur_node->cur_profile_nodes -= cur_node->cur_prof_min;

					// These next level nodes will be used to set the limits in the next cur_i_dim iteration.
					for(int i = cur_node->cur_prof_min; i <= cur_node->cur_prof_max; i++)
					{
						cur_node->cur_profile_nodes[i] = new t_histogram_node();
						cur_node->cur_profile_nodes[i]->dim = cur_i_dim+1; // These nodes are at the (cur_i_dim+1)^st dimension.
						cur_node->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
						cur_node->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
					} //i loop.

					if(cur_i_dim != n_dims - 1)
					{
						cur_node->cur_profile_prob_info = NULL;
					}
					else
					{
						cur_node->cur_profile_prob_info = new t_prob_info();
						cur_node->cur_profile_prob_info->counts = new double[cur_node->cur_prof_max - cur_node->cur_prof_min +1];
						cur_node->cur_profile_prob_info->counts -= cur_node->cur_prof_min;
                        cur_node->cur_profile_prob_info->probs = new double[cur_node->cur_prof_max - cur_node->cur_prof_min +1];
                        cur_node->cur_profile_prob_info->probs -= cur_node->cur_prof_min;

						// Initialize the counts.
						for(int i = cur_node->cur_prof_min; i <= cur_node->cur_prof_max; i++)
						{
							cur_node->cur_profile_prob_info->probs[i] = xlog(0.0);
							cur_node->cur_profile_prob_info->counts[i] = xlog(0.0);
						} // i loop.
					}
				} // profile min-max check.
                else
                {
                        fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
                        exit(1);
                }
			} // next level nodes existence check.
		} // i_sig loop.
	} // next_i_dim loop.

	if (__DUMP_HISTOGRAM_MSGS__)
	{
		fprintf(stderr, "Allocation of the node lists are done.\n");
	}
	
	return(main_node);
}

void merge_nodes_per_histogram_node_per_dim_online(t_histogram_node* main_node, t_histogram_node** hist_node_per_dim)
{
	// The current node must exist.
	if(main_node != NULL)
	{
		int i_dim = main_node->dim;

		if(hist_node_per_dim[i_dim] == NULL)
		{
			hist_node_per_dim[i_dim] = new t_histogram_node();
			hist_node_per_dim[i_dim]->cur_prof_min = 1000*1000;
			hist_node_per_dim[i_dim]->cur_prof_max = -1000*1000;
			hist_node_per_dim[i_dim]->cur_profile_nodes = NULL;
			hist_node_per_dim[i_dim]->cur_profile_prob_info = NULL;
			hist_node_per_dim[i_dim]->dim = i_dim;
		}

		// Merge the counts of the current node with the current node's dimension's histogram_node.
		if(main_node->cur_prof_min <= main_node->cur_prof_max)
		{
			t_histogram_node* cur_merged_hist_node = hist_node_per_dim[i_dim];

			// Set the merged limits.
			int merged_min = MIN(main_node->cur_prof_min, cur_merged_hist_node->cur_prof_min);
			int merged_max = MAX(cur_merged_hist_node->cur_prof_max, main_node->cur_prof_max);

			if(merged_max >= merged_min)
			{
				// Update the counts.
				double* merged_counts = new double[merged_max - merged_min + 1];
				merged_counts -= merged_min;
				double* merged_probs = new double[merged_max - merged_min + 1];
				merged_probs -= merged_min;
				t_histogram_node** merged_nodes = new t_histogram_node*[merged_max - merged_min + 1];
				merged_nodes -= merged_min;
				for(int i = merged_min; i <= merged_max; i++)
				{
					merged_nodes[i] = new t_histogram_node();
					merged_nodes[i]->cur_prof_min = 1000*1000;
					merged_nodes[i]->cur_prof_max = -1000*1000;

					merged_counts[i] = xlog(0.0);
					merged_probs[i] = xlog(0.0);
					if(i >= cur_merged_hist_node->cur_prof_min &&
						i <= cur_merged_hist_node->cur_prof_max)
					{
						merged_counts[i] = xlog_sum(merged_counts[i], cur_merged_hist_node->cur_profile_prob_info->counts[i]);
					}

					if(i >= main_node->cur_prof_min &&
						i <= main_node->cur_prof_max)
					{
						if(main_node->cur_profile_prob_info != NULL)
						{
							merged_counts[i] = xlog_sum(merged_counts[i], main_node->cur_profile_prob_info->counts[i]);
						}
						else
						{
							double cur_total = get_total_log_counts_per_multi_d_histogram(main_node->cur_profile_nodes[i]);
							merged_counts[i] = xlog_sum(merged_counts[i], cur_total);
						}
					}
				} // i loop.

				// Delete the memory for the currently merged node arrays.
				if(cur_merged_hist_node->cur_profile_nodes != NULL)
				{
					for(int i = cur_merged_hist_node->cur_prof_min; i <= cur_merged_hist_node->cur_prof_max; i++)
					{
						delete cur_merged_hist_node->cur_profile_nodes[i];
					} // i loop.

					delete [] (cur_merged_hist_node->cur_profile_nodes + cur_merged_hist_node->cur_prof_min);
				}

				if(cur_merged_hist_node->cur_profile_prob_info != NULL)
				{
					delete [] (cur_merged_hist_node->cur_profile_prob_info->counts + cur_merged_hist_node->cur_prof_min);
					delete [] (cur_merged_hist_node->cur_profile_prob_info->probs + cur_merged_hist_node->cur_prof_min);
					delete cur_merged_hist_node->cur_profile_prob_info;
				}		

				// Reallocate the node and update statistics for the merged histogram node.
				cur_merged_hist_node->cur_prof_min = merged_min;
				cur_merged_hist_node->cur_prof_max = merged_max;
				cur_merged_hist_node->cur_profile_prob_info = new t_prob_info();
				cur_merged_hist_node->cur_profile_prob_info->counts = merged_counts;
				cur_merged_hist_node->cur_profile_prob_info->probs = merged_probs;
				cur_merged_hist_node->cur_profile_nodes = merged_nodes;
			} // merged limit check.

			// Recurse on all the nodes.
			for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max;  i++)
			{
				merge_nodes_per_histogram_node_per_dim_online(main_node->cur_profile_nodes[i], hist_node_per_dim);
			} // i loop.
		} // main node limit check.
	} // main_node existence check.
}

void merge_nodes_per_histogram_node_per_dim(t_histogram_node* main_node, int n_dim_2_merge, t_histogram_node* cur_merged_hist_node)
{
	if(main_node->dim == n_dim_2_merge)
	{
		// Merge the current node with cur_merged_hist_node.
        int merged_min = MIN(main_node->cur_prof_min, cur_merged_hist_node->cur_prof_min);
        int merged_max = MAX(cur_merged_hist_node->cur_prof_max, main_node->cur_prof_max);

		if(merged_max >= merged_min)
		{
			// Update the counts.
			double* merged_counts = new double[merged_max - merged_min + 1];
			merged_counts -= merged_min;
			double* merged_probs = new double[merged_max - merged_min + 1];
			merged_probs -= merged_min;
			t_histogram_node** merged_nodes = new t_histogram_node*[merged_max - merged_min + 1];
			merged_nodes -= merged_min;
			for(int i = merged_min; i <= merged_max; i++)
			{
				merged_nodes[i] = new t_histogram_node();
				merged_nodes[i]->cur_prof_min = 1000*1000;
				merged_nodes[i]->cur_prof_max = -1000*1000;

				merged_counts[i] = xlog(0.0);
				merged_probs[i] = xlog(0.0);
				if(i >= cur_merged_hist_node->cur_prof_min &&
					i <= cur_merged_hist_node->cur_prof_max)
				{
					merged_counts[i] = xlog_sum(merged_counts[i], cur_merged_hist_node->cur_profile_prob_info->counts[i]);
				}

				if(i >= main_node->cur_prof_min &&
					i <= main_node->cur_prof_max)
				{
					if(main_node->cur_profile_prob_info != NULL)
					{
						merged_counts[i] = xlog_sum(merged_counts[i], main_node->cur_profile_prob_info->counts[i]);
					}
					else
					{
						double cur_total = get_total_log_counts_per_multi_d_histogram(main_node->cur_profile_nodes[i]);
						merged_counts[i] = xlog_sum(merged_counts[i], cur_total);
					}
				}
			} // i loop.

			// Delete the memory for the currently merged node arrays.
			if(cur_merged_hist_node->cur_profile_nodes != NULL)
			{
				for(int i = cur_merged_hist_node->cur_prof_min; i <= cur_merged_hist_node->cur_prof_max; i++)
				{
					delete cur_merged_hist_node->cur_profile_nodes[i];
				} // i loop.

				delete [] (cur_merged_hist_node->cur_profile_nodes + cur_merged_hist_node->cur_prof_min);
			}

			if(cur_merged_hist_node->cur_profile_prob_info != NULL)
			{
				delete [] (cur_merged_hist_node->cur_profile_prob_info->counts + cur_merged_hist_node->cur_prof_min);
				delete [] (cur_merged_hist_node->cur_profile_prob_info->probs + cur_merged_hist_node->cur_prof_min);
				delete cur_merged_hist_node->cur_profile_prob_info;
			}		

			// Reallocate the node and update statistics for the merged histogram node.
			cur_merged_hist_node->cur_prof_min = merged_min;
			cur_merged_hist_node->cur_prof_max = merged_max;
			cur_merged_hist_node->cur_profile_prob_info = new t_prob_info();
			cur_merged_hist_node->cur_profile_prob_info->counts = merged_counts;
			cur_merged_hist_node->cur_profile_prob_info->probs = merged_probs;
			cur_merged_hist_node->cur_profile_nodes = merged_nodes;
		} // limit check.
	}
	else if(main_node->dim < n_dim_2_merge)
	{
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			merge_nodes_per_histogram_node_per_dim(main_node->cur_profile_nodes[i], n_dim_2_merge, cur_merged_hist_node);
		} // i loop.
	}
}

void get_per_dimension_hist_nodes_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, t_histogram_node** hists_node_per_dim)
{
	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		hists_node_per_dim[i_dim] = NULL;
	} // i_dim loop.

	merge_nodes_per_histogram_node_per_dim_online(main_node, hists_node_per_dim);

	//for(int i_dim = 0; i_dim < n_dims; i_dim++)
	//{
	//	t_histogram_node* cur_merged_hist_node = new t_histogram_node();
	//	cur_merged_hist_node->cur_prof_min = 1000*1000;
	//	cur_merged_hist_node->cur_prof_max = -1000*1000;
	//	cur_merged_hist_node->cur_profile_nodes = NULL;
	//	cur_merged_hist_node->cur_profile_prob_info = NULL;
	//	merge_nodes_per_histogram_node_per_dim(main_node, i_dim, cur_merged_hist_node);
	//	hists_per_dim[i_dim] = cur_merged_hist_node;
	//} // i_dim loop.
}

void get_nodes_per_dimension_per_hist_node(t_histogram_node* main_node, int i_dim, vector<t_histogram_node*>* cur_list_nodes)
{
	if(main_node->dim == i_dim)
	{
		cur_list_nodes->push_back(main_node);
	}
	else if(main_node->dim < i_dim)
	{
		// Call for all the current profile nodes.
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			get_nodes_per_dimension_per_hist_node(main_node->cur_profile_nodes[i], i_dim, cur_list_nodes);
		} // i loop.
	}
	else
	{
		// Do not do anything, we passed the dimension requirement.
	}
}

t_histogram_node* get_multi_d_subset_hist_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, int* sorted_i_sub_dims, int i_sub_dim, int n_sub_dims)
{
	int cur_dim = sorted_i_sub_dims[i_sub_dim];

	if(cur_dim >= n_dims)
	{
		fprintf(stderr, "Cannot extract dimension since it is outside dimension limit for the histogram: %d, %d\n", i_sub_dim, n_dims);
		exit(1);
	}

	t_histogram_node* cur_dimension_node = NULL;
	vector<t_histogram_node*>* next_dimension_nodes = new vector<t_histogram_node*>();
	get_nodes_per_dimension_per_hist_node(main_node, cur_dim, next_dimension_nodes);

	// Merge all the nodes for the currently requested dimensions.
	for(int i_next_node = 0; i_next_node < (int)next_dimension_nodes->size(); i_next_node++)
	{
		// For the current node, get the next-next level nodes, surge them to the next level nodes of the current node, then merge them.
		if(next_dimension_nodes->at(i_next_node)->cur_prof_max >= next_dimension_nodes->at(i_next_node)->cur_prof_min)
		{
			if(cur_dimension_node == NULL)
			{
				cur_dimension_node = init_multi_d_histogram(next_dimension_nodes->at(i_next_node));
			}
			else
			{
				t_histogram_node* new_cur_dimension_node = merge_multi_d_hist_per_counts(cur_dimension_node, next_dimension_nodes->at(i_next_node));
				delete_multi_d_histogram(cur_dimension_node);
				cur_dimension_node = new_cur_dimension_node;
			}
		}
	} // i_next_node loop.

	if(cur_dimension_node != NULL &&
		cur_dimension_node->cur_prof_max >= cur_dimension_node->cur_prof_min)
	{
		// Is this the last dimension to process? If so, add the counts to the currently merged dimensions, delete the next level node information and return.
		if(i_sub_dim == n_sub_dims-1)
		{
			//fprintf(stderr, "Last dimension.\n");

			// Assign the counts to the current profile nodes.
			if(cur_dimension_node->cur_profile_prob_info == NULL)
			{
				cur_dimension_node->cur_profile_prob_info = new t_prob_info();
				cur_dimension_node->cur_profile_prob_info->counts = new double[cur_dimension_node->cur_prof_max - cur_dimension_node->cur_prof_min + 1];
				cur_dimension_node->cur_profile_prob_info->counts -= cur_dimension_node->cur_prof_min;

				cur_dimension_node->cur_profile_prob_info->probs= new double[cur_dimension_node->cur_prof_max - cur_dimension_node->cur_prof_min + 1];
				cur_dimension_node->cur_profile_prob_info->probs -= cur_dimension_node->cur_prof_min;

				for(int i = cur_dimension_node->cur_prof_min; i <= cur_dimension_node->cur_prof_max; i++)
				{
					cur_dimension_node->cur_profile_prob_info->counts[i] = get_total_log_counts_per_multi_d_histogram(cur_dimension_node->cur_profile_nodes[i]);
				} // i loop.
			}
			else
			{
				// There is already probability info, do not do anything, just return.
			}
		}
		else
		{
			// For all the next dimension nodes, get the next-next profile nodes.
			for(int i = cur_dimension_node->cur_prof_min; i <= cur_dimension_node->cur_prof_max; i++)
			{
				t_histogram_node* cur_next_next_profile_node = get_multi_d_subset_hist_per_multi_d_hist_node(cur_dimension_node->cur_profile_nodes[i], n_dims, sorted_i_sub_dims, i_sub_dim+1, n_sub_dims);
				delete_multi_d_histogram(cur_dimension_node->cur_profile_nodes[i]);

				if(cur_next_next_profile_node == NULL)
				{
					cur_dimension_node->cur_profile_nodes[i] = new t_histogram_node();
					cur_dimension_node->cur_profile_nodes[i]->cur_profile_prob_info = NULL;
					cur_dimension_node->cur_profile_nodes[i]->cur_profile_nodes = NULL;
					cur_dimension_node->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
					cur_dimension_node->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
				}
				else
				{
					cur_dimension_node->cur_profile_nodes[i] = cur_next_next_profile_node;
				}
			} // i loop.			
		}
	} // limit check.

	return(cur_dimension_node);
}

// Get the joint probability distribution of the dim1 and dim2 given a multidimensional probability distribution.
t_histogram_node* get_2D_hist_per_multi_d_hist_node(t_histogram_node* main_node, int n_dims, int i_dim1, int i_dim2)
{
	if(i_dim1 == i_dim2)
	{
		fprintf(stderr, "The dimensions must ne different from each other for getting the 2D distribution.\n");
		return(NULL);
	}

	if(i_dim1 > i_dim2)
	{
		int i_temp_dim1 = i_dim1;
		i_dim1 = i_dim2;
		i_dim2 = i_temp_dim1;
	}

	// Start from the bottom node, find the i_dim1 nodes, copy the i_dim1 nodes, then for each 1st dimension nodes, get the i_dim2 nodes, assign these as sub-nodes, 
	// then merge the i_dim1 nodes.
	vector<t_histogram_node*>* dim1_nodes = new vector<t_histogram_node*>();
	get_nodes_per_dimension_per_hist_node(main_node, i_dim1, dim1_nodes);

	fprintf(stderr, "Found %d nodes in dimension %d\n", (int)dim1_nodes->size(), i_dim1);

	t_histogram_node* merged_dim1_node = NULL;

	// For each node, get the i_dim2 nodes.
	for(int i_d1_n = 0; i_d1_n < (int)dim1_nodes->size(); i_d1_n++)
	{
		// Get the i_dim2 nodes per current i_dim1 node.
		t_histogram_node* cur_dim1_node = dim1_nodes->at(i_d1_n);

		// Make sure that the current dim1 node has valid limits.
        if(cur_dim1_node->cur_prof_max >= cur_dim1_node->cur_prof_min)
        {
			t_histogram_node* new_i_dim1_node = new t_histogram_node();
			new_i_dim1_node->cur_prof_min = cur_dim1_node->cur_prof_min;
			new_i_dim1_node->cur_prof_max = cur_dim1_node->cur_prof_max;
			new_i_dim1_node->cur_profile_prob_info = NULL;
			new_i_dim1_node->cur_profile_nodes = new t_histogram_node*[cur_dim1_node->cur_prof_max - cur_dim1_node->cur_prof_min + 1];
			new_i_dim1_node->cur_profile_nodes -= cur_dim1_node->cur_prof_min;
			new_i_dim1_node->dim = 0;
		
			// Go over all the nodes of the cur i_dim1 node, get all the i_dim2 nodes.
			for(int i = cur_dim1_node->cur_prof_min; i <= cur_dim1_node->cur_prof_max; i++)
			{
				if(cur_dim1_node->cur_profile_nodes[i]->cur_prof_max >= cur_dim1_node->cur_profile_nodes[i]->cur_prof_min)
				{
					// Get the i_dim2 nodes for the current i_dim1 nodes.
					vector<t_histogram_node*>* dim2_nodes_per_cur_dim1_node = new vector<t_histogram_node*>();
					get_nodes_per_dimension_per_hist_node(cur_dim1_node->cur_profile_nodes[i], i_dim2, dim2_nodes_per_cur_dim1_node);
                    //fprintf(stderr, "Dimension %d node %d: %d dimension %d nodes.\n", i_dim1, i,  dim2_nodes_per_cur_dim1_node->size(), i_dim2);

					if((int)dim2_nodes_per_cur_dim1_node->size() > 0)
					{
						t_histogram_node* cur_merged_i_dim2_node = NULL;

						// Go over all the i_dim2 nodes and merge them.
						for(int i2 = 0; i2 < (int)dim2_nodes_per_cur_dim1_node->size(); i2++)
						{
							if(cur_merged_i_dim2_node == NULL)
							{
								cur_merged_i_dim2_node = init_multi_d_histogram(dim2_nodes_per_cur_dim1_node->at(i2));
							}
							else
							{
								t_histogram_node* new_merged_i_dim2_node = merge_multi_d_hist_per_counts(cur_merged_i_dim2_node, dim2_nodes_per_cur_dim1_node->at(i2));
								delete_multi_d_histogram(cur_merged_i_dim2_node);
								cur_merged_i_dim2_node = new_merged_i_dim2_node;
							}
						} // i2 loop.

						// Copy the merged node to a new node which will be set as the new i_dim2 node.
						t_histogram_node* new_i_dim2_node = init_multi_d_histogram(cur_merged_i_dim2_node);						

						// At this point, all the i_dim2 nodes are merged for the current i_dim1 node.
						// Now replace the nodes with the counts for the i_dim2 nodes.
						if(new_i_dim2_node->cur_profile_prob_info == NULL)
						{
							new_i_dim2_node->cur_profile_prob_info = new t_prob_info();
							new_i_dim2_node->cur_profile_prob_info->counts = new double[cur_merged_i_dim2_node->cur_prof_max - cur_merged_i_dim2_node->cur_prof_min + 1];
							new_i_dim2_node->cur_profile_prob_info->counts -= cur_merged_i_dim2_node->cur_prof_min;
							new_i_dim2_node->cur_profile_prob_info->probs = new double[cur_merged_i_dim2_node->cur_prof_max - cur_merged_i_dim2_node->cur_prof_min + 1];
							new_i_dim2_node->cur_profile_prob_info->probs -= cur_merged_i_dim2_node->cur_prof_min;
						}

						// Set the current merged i_dim2 nodes: Some node surgery is done to make the i_dim2 node have probability information.
						for(int j = cur_merged_i_dim2_node->cur_prof_min; j <= cur_merged_i_dim2_node->cur_prof_max; j++)
						{
							if(cur_merged_i_dim2_node->cur_profile_prob_info == NULL)
							{
								new_i_dim2_node->cur_profile_prob_info->counts[j] = get_total_log_counts_per_multi_d_histogram(cur_merged_i_dim2_node->cur_profile_nodes[j]);
							}
							else
							{}

							// Free the memory for all the following profiles afterwards.
							delete_multi_d_histogram(new_i_dim2_node->cur_profile_nodes[j]);
							new_i_dim2_node->cur_profile_nodes[j] = new t_histogram_node();
							new_i_dim2_node->cur_profile_nodes[j]->cur_profile_prob_info = NULL;
							new_i_dim2_node->cur_profile_nodes[j]->cur_prof_min = 1000 * 1000;
							new_i_dim2_node->cur_profile_nodes[j]->cur_prof_max = -1000 * 1000;
							new_i_dim2_node->cur_profile_nodes[j]->dim = 2; // Must set the dimension for this node to 2 to make sure it is well defined.
						} // j loop.

						delete_multi_d_histogram(cur_merged_i_dim2_node);
			
						// Set the dimension of the node to 1.
						new_i_dim2_node->dim = 1;

						// Assign the currently merged node to the i_dim1 node. 
						new_i_dim1_node->cur_profile_nodes[i] = new_i_dim2_node;
					} // check if the i_dim2 nodes exist for the current i_dim1 node.
					else
					{
						fprintf(stderr, "The number of i_dim2 nodes is 0.\n");
						exit(1);
					}

					delete dim2_nodes_per_cur_dim1_node;
				} // limits check for the current node.
				else
				{
					t_histogram_node* cur_merged_i_dim2_node = new t_histogram_node();
					cur_merged_i_dim2_node->cur_prof_min = 1000*1000;
					cur_merged_i_dim2_node->cur_prof_max = -1000*1000;
					cur_merged_i_dim2_node->cur_profile_prob_info = NULL;
					cur_merged_i_dim2_node->dim = 1;

					new_i_dim1_node->cur_profile_nodes[i] = cur_merged_i_dim2_node;
				} // limits check for the current node.
			} // i loop.

			// Finally, merge the current i_dim1 node with i_dim2 nodes.
			fprintf(stderr, "Merging/Initing %d. dimension %d node.\n", i_d1_n, i_dim1);
			if(merged_dim1_node == NULL)
			{
				merged_dim1_node = init_multi_d_histogram(new_i_dim1_node);
				delete_multi_d_histogram(new_i_dim1_node);
			}
			else
			{
				t_histogram_node* new_merged_dim1_node = merge_multi_d_hist_per_counts(new_i_dim1_node, merged_dim1_node);
				delete_multi_d_histogram(merged_dim1_node);
				merged_dim1_node = new_merged_dim1_node;
			}
		} // Current i_dim1 limits check.
	} // i_d1_n loop.

	return(merged_dim1_node);
} // get_2D_hist_per_multi_d_hist_node

void linear_diff_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters)
{
	double* cur_node_counts = new double[n_vals + 3];
	double* smoothed_cur_node_counts = new double[n_vals + 2];
	for(int i = 0; i < n_vals; i++)
	{
		cur_node_counts[i] = xexp(count_array[i]);
		smoothed_cur_node_counts[i] = 0.0;
	} // i loop.

	//int l_signal = main_node->cur_prof_max - main_node->cur_prof_min + 1;
	int l_signal = n_vals;
	int l_half_win = (l_win-1)/2;

	double linear_non_self_weight = (1 - self_weight)/(double)(l_win-1);
	//double linear_one_min_self_weight = self_weight - linear_non_self_weight;
	//double log_one_min_self_weight = xlog(linear_one_min_self_weight);
	//double log_non_self_weight = xlog(linear_non_self_weight);

	double* cur_iter_non_self_weighted = new double[l_signal];
	for(int i_iter = 0; i_iter < n_smoothing_iters; i_iter++)
	{
		// Pre-compute the non-self weighted values for all the entries.	
		for(int i_s = 0; i_s < l_signal; i_s++)
		{
			cur_iter_non_self_weighted[i_s] = cur_node_counts[i_s] * linear_non_self_weight;
		} // i_s loop.

		// Compute the first smoothed value.
		// cur_sum keeps the non-weighted sum of all the current values. The centered value is added after cur_sum is updated.
		double cur_sum = 0.0;
		int i = 0;

		int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
		int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);
		for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
		{
			cur_sum += cur_iter_non_self_weighted[i_s];
		} // i_s loop.

		// Fix the self value at 0th value.
		smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

		// Set the previous ends.
		int prev_win_start = cur_win_start;
		int prev_win_end = cur_win_end;
		for(int i = 1; i < l_signal; i++)
		{
			// Update the current ends.
			cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			for(int i_s = prev_win_start; i_s < cur_win_start; i_s++)
			{
				cur_sum -= (cur_iter_non_self_weighted[i_s]); 
			}

			// Add the values to the right of this window.
			for(int i_s = prev_win_end+1; i_s <= cur_win_end; i_s++)
			{
				cur_sum += (cur_iter_non_self_weighted[i_s]);
			} // i_s loop.

			if(cur_sum < -0.00000001)
			{
				fprintf(stderr, "WTF %.20f\n", cur_node_counts[i]);
				exit(1);
			}
			else if(cur_sum < 0)
			{
				cur_sum = 0;
			}

			// Add the self value.
			smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

			// Update the previous window end coordinates.
			prev_win_start = cur_win_start;
			prev_win_end = cur_win_end;
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_counts[i] = smoothed_cur_node_counts[i];
			if(cur_node_counts[i] < -0.000000001)
			{
				fprintf(stderr, "WTF %.20f\n", cur_node_counts[i]);
				exit(1);
			}
			else if(cur_node_counts[i] < 0)
			{
				cur_node_counts[i] = 0;
			}
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_counts[i] = smoothed_cur_node_counts[i];
		} // i loop.
	} // i_iter loop.

	// Copy the counts.
	//for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	for(int i = 0; i < n_vals; i++)
	{
		count_array[i] = xlog(smoothed_cur_node_counts[i]);

		if(count_array[i] == xlog(0.0))
		{
			fprintf(stderr, "Smoothed value is 0.\n");
			exit(1);
		}
	} // i loop.

	delete [] cur_iter_non_self_weighted;
	delete [] cur_node_counts;
	delete [] smoothed_cur_node_counts;
} // linear_diff_smooth_next_linear_indexed_count_array

void linear_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters)
{
	double* cur_node_counts = new double[n_vals + 3];
	double* smoothed_cur_node_counts = new double[n_vals + 2];
	for(int i = 0; i < n_vals; i++)
	{
		cur_node_counts[i] = xexp(count_array[i]);
		smoothed_cur_node_counts[i] = 0.0;
	} // i loop.

	//int l_signal = main_node->cur_prof_max - main_node->cur_prof_min + 1;
	int l_signal = n_vals;
	int l_half_win = (l_win-1)/2;

	double linear_non_self_weight = (1 - self_weight)/(double)(l_win-1);
	double linear_one_min_self_weight = self_weight - linear_non_self_weight;
	//double log_one_min_self_weight = xlog(linear_one_min_self_weight);
	//double log_non_self_weight = xlog(linear_non_self_weight);

	for(int i_iter = 0; i_iter < n_smoothing_iters; i_iter++)
	{
		// Set the previous ends.
		//double cur_log_sum = LOG_OF_ZERO;
		double cur_lin_sum = 0.0;
		for(int i = 0; i < l_signal; i++)
		{
			// Update the current ends.
			int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			cur_lin_sum = 0.0;
			for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
			{
				//cur_log_sum = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i_s], log_non_self_weight)); 
				cur_lin_sum += (cur_node_counts[i_s] * linear_non_self_weight); 
			} // i_s loop.

			// Add the self value: Subtract the non-self part.
			//smoothed_cur_node_log_counts[i] = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i], log_one_min_self_weight));
			smoothed_cur_node_counts[i] = cur_lin_sum + (cur_node_counts[i] * linear_one_min_self_weight);
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_counts[i] = smoothed_cur_node_counts[i];
		} // i loop.
	} // i_iter loop.

	// Copy the counts.
	//for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	for(int i = 0; i < n_vals; i++)
	{
		count_array[i] = xlog(cur_node_counts[i]);

		if(count_array[i] == xlog(0.0))
		{
			fprintf(stderr, "Smoothed value is 0.\n");
			exit(1);
		}
	} // i loop.

	delete [] cur_node_counts;
	delete [] smoothed_cur_node_counts;
} // linear_smooth_next_linear_indexed_count_array 

void log_smooth_next_linear_indexed_count_array(double* count_array, int n_vals, double self_weight, int l_win, int n_smoothing_iters)
{
	double* cur_node_log_counts = new double[n_vals + 3];
	double* smoothed_cur_node_log_counts = new double[n_vals + 2];
	for(int i = 0; i < n_vals; i++)
	{
		cur_node_log_counts[i] = count_array[i];
		smoothed_cur_node_log_counts[i] = xlog(0.0);
	} // i loop.

	//int l_signal = main_node->cur_prof_max - main_node->cur_prof_min + 1;
	int l_signal = n_vals;
	int l_half_win = (l_win-1)/2;

	double linear_non_self_weight = (1 - self_weight)/(double)(l_win-1);
	double linear_one_min_self_weight = self_weight - linear_non_self_weight;
	double log_one_min_self_weight = xlog(linear_one_min_self_weight);
	double log_non_self_weight = xlog(linear_non_self_weight);

	for(int i_iter = 0; i_iter < n_smoothing_iters; i_iter++)
	{
		// Set the previous ends.
		double cur_log_sum = LOG_OF_ZERO;
		for(int i = 0; i < l_signal; i++)
		{
			// Update the current ends.
			int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			cur_log_sum = LOG_OF_ZERO;
			for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
			{
				cur_log_sum = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i_s], log_non_self_weight)); 
			} // i_s loop.

			// Add the self value: Subtract the non-self part.
			smoothed_cur_node_log_counts[i] = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i], log_one_min_self_weight));
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_log_counts[i] = smoothed_cur_node_log_counts[i];
		} // i loop.
	} // i_iter loop.

	// Copy the counts.
	//for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	for(int i = 0; i < n_vals; i++)
	{
		count_array[i] = cur_node_log_counts[i];

		if(count_array[i] == xlog(0.0))
		{
			fprintf(stderr, "Smoothed value is 0.\n");
			exit(1);
		}
	} // i loop.

	delete [] cur_node_log_counts;
	delete [] smoothed_cur_node_log_counts;
} // log_smooth_next_linear_indexed_count_array

void log_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_iters)
{
	if(l_win % 2 != 1)
	{
		fprintf(stderr, "l_win must be an odd number.\n");
		exit(1);
	}

	double* cur_node_log_counts = new double[n_vals + 3];
	double* smoothed_cur_node_log_counts = new double[n_vals + 2];
	int l_half_win = (l_win-1)/2;

	// Map the values to make sure the smoothing is done between existing values.
	int cur_linear_val = 0;
	int last_log_val = -1;

	int* i_2_i_mapped = new int[n_vals + 1];
	for(int i = 0; i < n_vals; i++)
	{
		// Set no non-existing at first.
		i_2_i_mapped[i] = -1;
	} // i loop.
	
	int i_mapped = 0;
	while(1)
	{
		double cur_log_val = floor(xlog(cur_linear_val+1.0) / xlog(log_base));
		if(cur_log_val != last_log_val &&
			cur_log_val >= 0 &&
			cur_log_val < n_vals)
		{	
			last_log_val = cur_log_val;

			// Set the index to the mapped index.
			i_2_i_mapped[(int)cur_log_val] = i_mapped;
			
			// Set the mapped buffers: Copy the value in the 
			cur_node_log_counts[i_mapped] = count_array[(int)cur_log_val];
			smoothed_cur_node_log_counts[i_mapped] = LOG_OF_ZERO;

			//fprintf(stderr, "Mapping: %lf->%d: %lf\n", cur_log_val, i_mapped, cur_node_counts[i_mapped]);

			// Increment the mapped buffer index.
			i_mapped++;
		}
		
		// If we are at the end of the buffer, break.
		if(cur_log_val == n_vals)
		{
			break;
		}

		// Update the current linear value.
		cur_linear_val++;
	} // linear value.	

	// Smooth the i_mapped many values in the histogram.
	int l_signal = i_mapped;

	double log_non_self_weight = xlog((1.0-self_weight) / ((double)(l_win-1)));
	double log_self_weight = xlog(self_weight);

	for(int i_iter = 0; i_iter < n_iters; i_iter++)
	{
		// Set the previous ends.
		double cur_log_sum = LOG_OF_ZERO;
		for(int i = 0; i < l_signal; i++)
		{
			// Update the current ends.
			int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			cur_log_sum = LOG_OF_ZERO;
			for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
			{
				cur_log_sum = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i_s], log_non_self_weight)); 
			} // i_s loop.

			// Add the self value: Subtract the non-self part.
			smoothed_cur_node_log_counts[i] = xlog_sum(cur_log_sum, xlog_mul(cur_node_log_counts[i], xlog_sub(log_self_weight, log_non_self_weight)));
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_log_counts[i] = smoothed_cur_node_log_counts[i];
		} // i loop.
	} // i_iter loop.

	// Copy the counts back using the mapped indices.
	for(int i = 0; i < n_vals; i++)
	{
		if(i_2_i_mapped[i] != -1)
		{
			double current_smoothed_val = cur_node_log_counts[i_2_i_mapped[i]];
		
			if(current_smoothed_val == 0.0)
			{
				fprintf(stderr, "An existing entry is not smoothed @ %d (%d)\n", i_2_i_mapped[i], i);
				exit(0);
			}

			// Copy the value to the soothed counts.
			count_array[i] = current_smoothed_val;
		}
		else
		{
			count_array[i] = LOG_OF_ZERO;
		}
	} // i loop.

	// Free memory.
	delete [] i_2_i_mapped;
	delete [] cur_node_log_counts;
	delete [] smoothed_cur_node_log_counts;
} // log_smooth_next_log_indexed_count_array

void linear_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_iters)
{
	if(l_win % 2 != 1)
	{
		fprintf(stderr, "l_win must be an odd number.\n");
		exit(0);
	}

	double* cur_node_linear_counts = new double[n_vals + 3];
	double* smoothed_cur_node_linear_counts = new double[n_vals + 2];
	int l_half_win = (l_win-1)/2;

	// Map the values to make sure the smoothing is done between existing values.
	int cur_linear_val = 0;
	int last_log_val = -1;

	int* i_2_i_mapped = new int[n_vals + 1];
	for(int i = 0; i < n_vals; i++)
	{
		// Set no non-existing at first.
		i_2_i_mapped[i] = -1;
	} // i loop.
	
	int i_mapped = 0;
	while(1)
	{
		double cur_log_val = floor(xlog(cur_linear_val+1.0) / xlog(log_base));
		if(cur_log_val != last_log_val &&
			cur_log_val >= 0 &&
			cur_log_val < n_vals)
		{	
			last_log_val = cur_log_val;

			// Set the index to the mapped index.
			i_2_i_mapped[(int)cur_log_val] = i_mapped;
			
			// Set the mapped buffers: Copy the value in the 
			cur_node_linear_counts[i_mapped] = xexp(count_array[(int)cur_log_val]);
			smoothed_cur_node_linear_counts[i_mapped] = 0.0;

			// Increment the mapped buffer index.
			i_mapped++;
		}
		
		// If we are at the end of the buffer, break.
		if(cur_log_val == n_vals)
		{
			break;
		}

		// Update the current linear value.
		cur_linear_val++;
	} // linear value.

	// Smooth the i_mapped many values in the histogram.
	int l_signal = i_mapped;

	double non_self_weight = ((1.0-self_weight) / ((double)(l_win-1)));

	for(int i_iter = 0; i_iter < n_iters; i_iter++)
	{
		// Set the previous ends.
		double cur_sum = 0.0;
		for(int i = 0; i < l_signal; i++)
		{
			// Update the current ends.
			int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			cur_sum = 0.0;
			for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
			{
				cur_sum = (cur_sum + (cur_node_linear_counts[i_s] * non_self_weight)); 
			} // i_s loop.

			// Add the self value: Subtract the non-self part.
			smoothed_cur_node_linear_counts[i] = (cur_sum + (cur_node_linear_counts[i] * (self_weight - non_self_weight)));
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_linear_counts[i] = smoothed_cur_node_linear_counts[i];
		} // i loop.
	} // i_iter loop.

	// Copy the counts back using the mapped indices.
	for(int i = 0; i < n_vals; i++)
	{
		if(i_2_i_mapped[i] != -1)
		{
			double current_smoothed_val = cur_node_linear_counts[i_2_i_mapped[i]];
		
			if(current_smoothed_val == 0.0)
			{
				fprintf(stderr, "An existing entry is not smoothed @ %d (%d)\n", i_2_i_mapped[i], i);
				exit(0);
			}

			// Copy the value to the soothed counts.
			count_array[i] = xlog(current_smoothed_val);
		}
		else
		{
			count_array[i] = LOG_OF_ZERO;
		}
	} // i loop.

	// Free memory.
	delete [] i_2_i_mapped;
	delete [] cur_node_linear_counts;
	delete [] smoothed_cur_node_linear_counts;
} // linear_smooth_next_log_indexed_count_array

/*
Count over the linear space to get the mapping of the log array to consecutive array.
The loop counts till max_log_val.
*/
void get_log_mapping(double max_log_val, int* log_i_2_consecutive_i, double log_base)
{
	fprintf(stderr, "Getting log to consecutive array mapping indices for max log value of %lf\n", max_log_val);

	// Initialize the mapped indices.
	for(int cur_log_val = 0; cur_log_val <= max_log_val; cur_log_val++)
	{
		log_i_2_consecutive_i[cur_log_val] = -1;
	} // cur_log_val loop.

	double cur_linear_val = 0.0;
	double last_log_val = -1.0;
	int i_mapped = 0;
	while(1)
	{
		double cur_log_val = floor(xlog(cur_linear_val+1.0) / xlog(log_base));
		if(cur_log_val != last_log_val &&
			cur_log_val >= 0 &&
			cur_log_val <= max_log_val)
		{	
			last_log_val = cur_log_val;

			// Set the index to the mapped index.
			log_i_2_consecutive_i[(int)cur_log_val] = i_mapped;

			// Increment the mapped buffer index.
			i_mapped++;
		}
		
		// If we are at the end of the buffer, break.
		if(cur_log_val == max_log_val)
		{
			break;
		}

		// Update the current linear value.
		cur_linear_val++;
	} // linear value.
}

void linear_diff_smooth_next_log_indexed_count_array(double* count_array, int* log_i_2_consecutive_i, int max_log_val, double log_base, double self_weight, int l_win, int n_smoothing_iters)
{
	if(l_win % 2 != 1)
	{
		fprintf(stderr, "l_win must be an odd number.\n");
		exit(0);
	}

	double* cur_node_counts = new double[max_log_val + 3];
	double* smoothed_cur_node_counts = new double[max_log_val + 2];
	int l_half_win = (l_win-1)/2;

	int l_signal = 0;

	// Copy the values.
	for(int i = 0; i < max_log_val; i++)
	{
		if(log_i_2_consecutive_i[i] != -1)
		{
			smoothed_cur_node_counts[i] = 0.0;
			cur_node_counts[log_i_2_consecutive_i[i]] = xexp(count_array[i]);
			l_signal = log_i_2_consecutive_i[i]+1;
		}
	} // i loop.

	double linear_non_self_weight = (1 - self_weight)/(double)(l_win-1);
	//double linear_one_min_self_weight = self_weight - linear_non_self_weight;
	//double log_one_min_self_weight = xlog(linear_one_min_self_weight);
	//double log_non_self_weight = xlog(linear_non_self_weight);

	double* cur_iter_non_self_weighted = new double[l_signal];
	for(int i_iter = 0; i_iter < n_smoothing_iters; i_iter++)
	{
		// Pre-compute the non-self weighted values for all the entries.	
		for(int i_s = 0; i_s < l_signal; i_s++)
		{
			cur_iter_non_self_weighted[i_s] = cur_node_counts[i_s] * linear_non_self_weight;
		} // i_s loop.

		// Compute the first smoothed value.
		// cur_sum keeps the non-weighted sum of all the current values. The centered value is added after cur_sum is updated.
		double cur_sum = 0.0;
		int i = 0;

		int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
		int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);
		for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
		{
			cur_sum += cur_iter_non_self_weighted[i_s];
		} // i_s loop.

		// Fix the self value at 0th value.
		smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

		// Set the previous ends.
		int prev_win_start = cur_win_start;
		int prev_win_end = cur_win_end;
		for(int i = 1; i < l_signal; i++)
		{
			// Update the current ends.
			cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			for(int i_s = prev_win_start; i_s < cur_win_start; i_s++)
			{
				cur_sum -= (cur_iter_non_self_weighted[i_s]); 
			}

			// Add the values to the right of this window.
			for(int i_s = prev_win_end+1; i_s <= cur_win_end; i_s++)
			{
				cur_sum += (cur_iter_non_self_weighted[i_s]);
			} // i_s loop.

			if(cur_sum < -0.0001)
			{
				fprintf(stderr, "cur_sum is negative!\n");
				exit(0);
			}
			else if(cur_sum < 0.0)
			{
				cur_sum = 0.0;
			}

			// Add the self value.
			smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

			// Update the previous window end coordinates.
			prev_win_start = cur_win_start;
			prev_win_end = cur_win_end;
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_counts[i] = smoothed_cur_node_counts[i];
			if(cur_node_counts[i] < -0.000001)
			{
				fprintf(stderr, "WTF %.20f\n", cur_node_counts[i]);
				exit(0);
			}
			else if(cur_node_counts[i] < 0)
			{
				cur_node_counts[i] = 0.0;
			}
		} // i loop.
	} // i_iter loop.

	// Copy the counts back using the mapped indices.
	for(int i = 0; i < max_log_val; i++)
	{
		if(log_i_2_consecutive_i[i] != -1)
		{
			double current_smoothed_val = cur_node_counts[log_i_2_consecutive_i[i]];
		
			if(current_smoothed_val == 0.0)
			{
				fprintf(stderr, "An existing entry is not smoothed @ %d (%d)\n", log_i_2_consecutive_i[i], i);
				exit(0);
			}

			// Copy the value to the soothed counts.
			count_array[i] = xlog(current_smoothed_val);
		}
		else
		{
			count_array[i] = LOG_OF_ZERO;
		}
	} // i loop.

	// Free memory.
	delete [] cur_node_counts;
	delete [] smoothed_cur_node_counts;
	delete [] cur_iter_non_self_weighted;
}

void linear_diff_smooth_next_log_indexed_count_array(double* count_array, int n_vals, double log_base, double self_weight, int l_win, int n_smoothing_iters)
{
	if(l_win % 2 != 1)
	{
		fprintf(stderr, "l_win must be an odd number.\n");
		exit(0);
	}

	double* cur_node_counts = new double[n_vals + 3];
	double* smoothed_cur_node_counts = new double[n_vals + 2];
	int l_half_win = (l_win-1)/2;

	// Map the values to make sure the smoothing is done between existing values.
	int cur_linear_val = 0;
	int last_log_val = -1;

	int* i_2_i_mapped = new int[n_vals + 1];
	for(int i = 0; i < n_vals; i++)
	{
		// Set no non-existing at first.
		i_2_i_mapped[i] = -1;
	} // i loop.
	
	int i_mapped = 0;
	while(1)
	{
		double cur_log_val = floor(xlog(cur_linear_val+1.0) / xlog(log_base));
		if(cur_log_val != last_log_val &&
			cur_log_val >= 0 &&
			cur_log_val < n_vals)
		{	
			last_log_val = cur_log_val;

			// Set the index to the mapped index.
			i_2_i_mapped[(int)cur_log_val] = i_mapped;
			
			// Set the mapped buffers: Copy the value in the 
			cur_node_counts[i_mapped] = xexp(count_array[(int)cur_log_val]);
			smoothed_cur_node_counts[i_mapped] = 0.0;

			// Increment the mapped buffer index.
			i_mapped++;
		}
		
		// If we are at the end of the buffer, break.
		if(cur_log_val == n_vals)
		{
			break;
		}

		// Update the current linear value.
		cur_linear_val++;
	} // linear value.

	// Smooth the i_mapped many values in the histogram.
	int l_signal = i_mapped;

	//double log_one_min_self_weight = xlog((1.0-self_weight) / ((double)(l_win-1)));
	//double log_self_weight = xlog(self_weight);
	double linear_non_self_weight = (1 - self_weight)/(double)(l_win-1);
	//double linear_one_min_self_weight = self_weight - linear_non_self_weight;
	//double log_one_min_self_weight = xlog(linear_one_min_self_weight);
	//double log_non_self_weight = xlog(linear_non_self_weight);

	double* cur_iter_non_self_weighted = new double[l_signal];
	for(int i_iter = 0; i_iter < n_smoothing_iters; i_iter++)
	{
		// Pre-compute the non-self weighted values for all the entries.	
		for(int i_s = 0; i_s < l_signal; i_s++)
		{
			cur_iter_non_self_weighted[i_s] = cur_node_counts[i_s] * linear_non_self_weight;
		} // i_s loop.

		// Compute the first smoothed value.
		// cur_sum keeps the non-weighted sum of all the current values. The centered value is added after cur_sum is updated.
		double cur_sum = 0.0;
		int i = 0;

		int cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
		int cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);
		for(int i_s = cur_win_start; i_s <= cur_win_end; i_s++)
		{
			cur_sum += cur_iter_non_self_weighted[i_s];
		} // i_s loop.

		// Fix the self value at 0th value.
		smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

		// Set the previous ends.
		int prev_win_start = cur_win_start;
		int prev_win_end = cur_win_end;
		for(int i = 1; i < l_signal; i++)
		{
			// Update the current ends.
			cur_win_start = (i>=l_half_win)?(i-l_half_win):(0);
			cur_win_end = ((i+l_half_win) < l_signal)?(i+l_half_win):(l_signal-1);

			// Subtract the values to the left of this window.
			for(int i_s = prev_win_start; i_s < cur_win_start; i_s++)
			{
				cur_sum -= (cur_iter_non_self_weighted[i_s]); 
			}

			// Add the values to the right of this window.
			for(int i_s = prev_win_end+1; i_s <= cur_win_end; i_s++)
			{
				cur_sum += (cur_iter_non_self_weighted[i_s]);
			} // i_s loop.

			if(cur_sum < -0.0001)
			{
				fprintf(stderr, "cur_sum is negative!\n");
				exit(0);
			}
			else if(cur_sum < 0.0)
			{
				cur_sum = 0.0;
			}

			// Add the self value.
			smoothed_cur_node_counts[i] = cur_sum + (cur_node_counts[i] * self_weight - cur_iter_non_self_weighted[i]);

			// Update the previous window end coordinates.
			prev_win_start = cur_win_start;
			prev_win_end = cur_win_end;
		} // i loop.

		// Copy the smoothed values.
		for(int i = 0; i < l_signal; i++)
		{
			cur_node_counts[i] = smoothed_cur_node_counts[i];
			if(cur_node_counts[i] < -0.000001)
			{
				fprintf(stderr, "WTF %.20f\n", cur_node_counts[i]);
				exit(0);
			}
			else if(cur_node_counts[i] < 0)
			{
				cur_node_counts[i] = 0.0;
			}
		} // i loop.
	} // i_iter loop.

	// Copy the counts back using the mapped indices.
	for(int i = 0; i < n_vals; i++)
	{
		if(i_2_i_mapped[i] != -1)
		{
			double current_smoothed_val = cur_node_counts[i_2_i_mapped[i]];
		
			if(current_smoothed_val == 0.0)
			{
				fprintf(stderr, "An existing entry is not smoothed @ %d (%d)\n", i_2_i_mapped[i], i);
				exit(0);
			}

			// Copy the value to the soothed counts.
			count_array[i] = xlog(current_smoothed_val);
		}
		else
		{
			count_array[i] = LOG_OF_ZERO;
		}
	} // i loop.

	// Free memory.
	delete [] i_2_i_mapped;
	delete [] cur_node_counts;
	delete [] smoothed_cur_node_counts;
	delete [] cur_iter_non_self_weighted;
} // linear_diff_smooth_next_log_indexed_count_array

/*
Smooth each distribution, reassign the counts
*/
void smooth_histogram_node_counts(t_histogram_node* main_node, 
	int i_prof_node, 
	t_histogram_node** main_nodes_per_dims, 
	bool* is_log_per_dimension, 
	int* log_i_2_consecutive_i,
	double log_base, 
	double self_weight, 
	int l_win, 
	int n_smoothing_iters)
{
	// This node is never encountered in the data, skip it.
	if(main_node->cur_prof_max < main_node->cur_prof_min)
	{
		return;
	}

	// Get the smoothed counts from the node.
	double total_smoothed_counts_per_main = main_node->smoothed_total_counts;

	// Assign the smoothed counts to all the nodes.
	double total_counts_per_main = get_total_log_counts_per_multi_d_histogram(main_node);
	
	int i_dim = main_node->dim;

	double* count_array = new double[main_nodes_per_dims[i_dim]->cur_prof_max - main_nodes_per_dims[i_dim]->cur_prof_min + 1];

	if(total_counts_per_main == xlog(0.0))
	{
		if(main_node->cur_prof_max < main_node->cur_prof_min)
		{
			fprintf(stderr, "%s(%d): The limits are not set for a valid node.\n", __FILE__, __LINE__);
			exit(0);
		}

		// The counts will be copied from the per dimension distributio	ns.
		for(int i = main_nodes_per_dims[main_node->dim]->cur_prof_min; i <= main_nodes_per_dims[main_node->dim]->cur_prof_max; i++)
		{
			count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = main_nodes_per_dims[i_dim]->cur_profile_prob_info->counts[i];
		} // i loop.
	}
	else
	{
		// Following loop fills the count values using the counts from the current profile's distribution.
		for(int i = main_nodes_per_dims[i_dim]->cur_prof_min; i <= main_nodes_per_dims[i_dim]->cur_prof_max; i++)
		{
			if(i >= main_node->cur_prof_min && 
				i <= main_node->cur_prof_max)
			{
				if(main_node->cur_profile_nodes[i]->cur_prof_max >= main_node->cur_profile_nodes[i]->cur_prof_min)
				{
					count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = get_total_log_counts_per_multi_d_histogram(main_node->cur_profile_nodes[i]);
				}
				else if(main_node->cur_profile_prob_info != NULL)
				{
					count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = main_node->cur_profile_prob_info->counts[i];
				}
				else
				{
					count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = xlog(0.0);
				}
			}
			else
			{
				count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = xlog(0.0);
			}
		} // i loop.

		// Smooth the counts.
		if(is_log_per_dimension[main_node->dim])
		{
			linear_diff_smooth_next_log_indexed_count_array(count_array, 
				log_i_2_consecutive_i, 
				(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1),
				log_base, 
				self_weight, 
				l_win, 
				n_smoothing_iters);
			//linear_smooth_next_log_indexed_count_array(count_array, 
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	log_base, 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
			//log_smooth_next_log_count_array(count_array, 
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	log_base, 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
		}
		else
		{
			linear_diff_smooth_next_linear_indexed_count_array(count_array,
				(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
				self_weight, 
				l_win, 
				n_smoothing_iters);
			//linear_smooth_next_linear_indexed_count_array(count_array,
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
			//log_smooth_next_linear_count_array(count_array,
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
		}
	} // Total counts check.

	double total_smoothed_counts = xlog(0.0);
	for(int i = main_nodes_per_dims[i_dim]->cur_prof_min; i <= main_nodes_per_dims[i_dim]->cur_prof_max; i++)
	{		
		total_smoothed_counts = xlog_sum(total_smoothed_counts, count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min]);
	} // i loop.

	// This scaling ensures that the smoothed 
	double log_scaling_factor = xlog_div(total_smoothed_counts_per_main, total_smoothed_counts);

	// Scale the counts. When the scaling is done after smoothing, it should be more stable since the scaling may be very small.
	//char cur_smoothed_dist_fp[1000];
	//sprintf(cur_smoothed_dist_fp, "smoothed_dim_%d_i_%d_prof_%d_%d.txt", main_node->dim, i_prof_node, main_node->cur_prof_min, main_node->cur_prof_max);
	//FILE* f_cur_dist = open_f(cur_smoothed_dist_fp, "w");
	for(int i = main_nodes_per_dims[i_dim]->cur_prof_min; i <= main_nodes_per_dims[i_dim]->cur_prof_max; i++)
	{
		if(i >= main_node->cur_prof_min && 
			i <= main_node->cur_prof_max)
		{
			main_node->cur_profile_nodes[i]->smoothed_total_counts = xlog_mul(count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min], log_scaling_factor);
			//fprintf(f_cur_dist, "%d\t%lf\n", i, main_node->cur_profile_nodes[i]->smoothed_total_counts);
		}
	} // i loop.
	//fclose(f_cur_dist);

	// Since we have smoothed values assigned to each node, can call for each recursively.
	for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	{
		if(main_node->cur_profile_prob_info == NULL)
		{
			//fprintf(stderr, "%s(%d): Dimension %d: %d-%d\n", __FILE__, __LINE__, main_node->dim, main_node->cur_prof_min, main_node->cur_prof_max);
			smooth_histogram_node_counts(main_node->cur_profile_nodes[i], i, main_nodes_per_dims, is_log_per_dimension, log_i_2_consecutive_i, log_base, self_weight, l_win, n_smoothing_iters);
		}
		else
		{
			// Assign the probabilities for the current distribution and return: This assigns the smoothed counts to the counts so that the 
			// probabilities can be directly computed.
			//fprintf(stderr, "%s(%d): Dimension %d: %d-%d\n", __FILE__, __LINE__, main_node->dim, main_node->cur_prof_min, main_node->cur_prof_max);
			main_node->cur_profile_prob_info->counts[i] = main_node->cur_profile_nodes[i]->smoothed_total_counts;
		}
	} // i loop.

	delete [] count_array;
} // smooth_log_histogram_node_counts

/*
Recursively call the histogram node smoothing function.
*/
void smooth_multi_d_histogram_counts(t_histogram_node* main_node, bool* is_log_per_dimension, double log_base, int n_dims, double self_weight, int l_win, int n_smoothing_iters)
{
	// Get the per dimension distributions.
	t_histogram_node** main_nodes_per_dims = new t_histogram_node*[n_dims];
	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		main_nodes_per_dims[i_dim] = NULL;
	} // i_dim loop.

	get_per_dimension_hist_nodes_per_multi_d_hist_node(main_node, n_dims, main_nodes_per_dims);

	// Generate mapping indices.
	int* log_i_2_consecutive_i = new int[310];
	get_log_mapping(300.0, log_i_2_consecutive_i, log_base);

	// Smooth all the per dimension distributions so we do not have to smooth them again everytime.
	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		fprintf(stderr, "Smoothing distribution for %d. dimension.\n", i_dim);

		// The counts will be copied from the per dimension distributions.
		double* count_array = new double[main_nodes_per_dims[i_dim]->cur_prof_max - main_nodes_per_dims[i_dim]->cur_prof_min + 1];
		for(int i = main_nodes_per_dims[i_dim]->cur_prof_min; i <= main_nodes_per_dims[i_dim]->cur_prof_max; i++)
		{
			count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min] = main_nodes_per_dims[i_dim]->cur_profile_prob_info->counts[i];
		} // i loop.

		// Smooth the counts.
		if(is_log_per_dimension[i_dim])
		{
			if(main_nodes_per_dims[i_dim]->cur_prof_min != 0)
			{
				fprintf(stderr, "The minimum is not 0.\n");
				exit(0);
			}

			linear_diff_smooth_next_log_indexed_count_array(count_array,
				log_i_2_consecutive_i, 
				(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
				log_base, 
				self_weight, 
				l_win, 
				n_smoothing_iters);
			//log_smooth_next_log_indexed_count_array(count_array,
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	log_base, 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
		}
		else
		{
			linear_diff_smooth_next_linear_indexed_count_array(count_array,
				(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
				self_weight, 
				l_win, 
				n_smoothing_iters);
			//linear_smooth_next_linear_indexed_count_array(count_array,
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
			//log_smooth_next_linear_count_array(count_array,
			//	(main_nodes_per_dims[i_dim]->cur_prof_max-main_nodes_per_dims[i_dim]->cur_prof_min+1), 
			//	self_weight, 
			//	l_win, 
			//	n_smoothing_iters);
		}

		// Copy the smooth counts.
		for(int i = main_nodes_per_dims[i_dim]->cur_prof_min; i <= main_nodes_per_dims[i_dim]->cur_prof_max; i++)
		{
			main_nodes_per_dims[i_dim]->cur_profile_prob_info->counts[i] = count_array[i-main_nodes_per_dims[i_dim]->cur_prof_min];
		} // i loop.

		delete [] count_array;
	} // i_dim loop.

	double total_counts = get_total_log_counts_per_multi_d_histogram(main_node);

	// Recursively smooth the nodes.
	main_node->smoothed_total_counts = total_counts;
	
	smooth_histogram_node_counts(main_node, 0, main_nodes_per_dims, is_log_per_dimension, log_i_2_consecutive_i, log_base, self_weight, l_win, n_smoothing_iters);
	
	// Delete the histogram memory.
	for(int i_dim = 0; i_dim < n_dims; i_dim++)
	{
		delete_multi_d_histogram(main_nodes_per_dims[i_dim]);
	} // i_dim loop.

	delete [] main_nodes_per_dims;

	fprintf(stderr, "Histogram smoothing is finished.\n");
}

t_histogram_node* get_multi_d_hist_per_one_indexed_dat(double** profiles, int n_dims, int n_signal_wins)
{
	//// Floorize all the profiles.
	//for(int i_d = 0; i_d < n_dims; i_d++)
	//{
	//	floorize_profile(profiles[i_d], n_signal_wins);
	//} // i_d loop.

	// Initialize the histogram per data.
	t_histogram_node* main_node = init_histogram_limits_per_data(profiles, n_dims, n_signal_wins);

	// Count the data values.
	count_data_per_histogram_node_per_one_indexed_data(main_node, profiles, n_dims, n_signal_wins);

	return(main_node);
}

void prune_histogram_node_memory(t_histogram_node* main_node, int& n_pruned_nodes)
{
	for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	{
		if(main_node->cur_profile_nodes[i]->cur_prof_max < main_node->cur_profile_nodes[i]->cur_prof_min)
		{
			n_pruned_nodes++;
			delete(main_node->cur_profile_nodes[i]);
			main_node->cur_profile_nodes[i] = NULL;
		}
		else
		{
			prune_histogram_node_memory(main_node->cur_profile_nodes[i], n_pruned_nodes);
		}
	} // i loop.
}

double get_total_memory_per_histogram_node(t_histogram_node* main_node)
{
	double cur_total_mem = 0.0;

	if(main_node == NULL)
	{
		return(cur_total_mem);
	}

	// Add the memory for the current node.
	cur_total_mem += sizeof(t_histogram_node);

	if(main_node->cur_prof_max >= main_node->cur_prof_min)
	{
		cur_total_mem += (main_node->cur_prof_max - main_node->cur_prof_min + 1) * sizeof(t_histogram_node*);
	}

	if(main_node->cur_profile_prob_info != NULL)
	{
		cur_total_mem += sizeof(t_prob_info);
		cur_total_mem += (main_node->cur_prof_max - main_node->cur_prof_min + 1) * sizeof(double);
		cur_total_mem += (main_node->cur_prof_max - main_node->cur_prof_min + 1) * sizeof(double);
	}

	// Get memory for all the next profile nodes.
	for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
	{
		cur_total_mem += get_total_memory_per_histogram_node(main_node->cur_profile_nodes[i]);
	} // i loop.

	return(cur_total_mem);
}

void dump_multi_d_histogram_text(t_histogram_node* main_node, int n_dims, char* fp)
{
	FILE* f_hist = open_f(fp, "w");
	int* dist_posn_per_dim = new int[n_dims];
	dump_multi_d_histogram_node_text(main_node, f_hist, dist_posn_per_dim, 0);
	fclose(f_hist);
}

void dump_multi_d_histogram_node_text(t_histogram_node* main_node, FILE* f_hist, int* dist_posn_per_dim, int i_last_dim)
{
	if(main_node->cur_profile_prob_info != NULL)
	{
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			// dump the values with higher than 0 count.
			if(main_node->cur_profile_prob_info->probs[i] != xlog(0.0))
			{
				// Dump the current list of values, then dump the count.
				for(int i_dim = 0; i_dim < i_last_dim; i_dim++)
				{
					fprintf(f_hist, "%d\t", dist_posn_per_dim[i_dim]);
				} // i_dim loop.
				fprintf(f_hist, "%d\t", i);

				//fprintf(f_hist, "%lf\n", main_node->cur_profile_prob_info->counts[i]);
				fprintf(f_hist, "%lf\n", main_node->cur_profile_prob_info->probs[i]);
			}
		}
	}
	else if(main_node->cur_prof_min <= main_node->cur_prof_max)
	{
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			dist_posn_per_dim[i_last_dim] = i;
			dump_multi_d_histogram_node_text(main_node->cur_profile_nodes[i], f_hist, dist_posn_per_dim, i_last_dim+1);
		} // i loop.
	}
}

//void dump_single_d_histogram_text(t_histogram_node* main_node, char* fp)
//{
//	if(main_node->cur_profile_prob_info == NULL)
//	{
//		fprintf(stderr, "Cannot dump a multi dimensional profile in text. Has multiple dimensions.\n");
//		return;
//	}
//
//	FILE* f_hist_node = open_f(fp, "wb");
//
//	// Recursively dump the nodes.
//	for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
//	{
//		fprintf(f_hist_node, "%d\t%lf\t%lf\n", i, main_node->cur_profile_prob_info->counts[i], main_node->cur_profile_prob_info->probs[i]);
//	} // i loop.
//
//	fclose(f_hist_node);
//}

void dump_multi_d_histogram(t_histogram_node* main_node, char* fp)
{
	FILE* f_hist_node = open_f(fp, "wb");

	// Recursively dump the nodes.
	dump_histogram_node(main_node, f_hist_node);

	fclose(f_hist_node);
}

const int prob_info_marker = 1234;
const int no_prob_info_marker = 4321;
void dump_histogram_node(t_histogram_node* node, FILE* f_hist_node)
{
	// Dump the information for the current node.
	fwrite(&(node->dim), sizeof(int), 1, f_hist_node);
	fwrite(&(node->cur_prof_min), sizeof(int), 1, f_hist_node);
	fwrite(&(node->cur_prof_max), sizeof(int), 1, f_hist_node);	

	// Dump the counts and probabilities if they exist.
	if(node->cur_profile_prob_info != NULL)
	{
		fwrite(&(prob_info_marker), sizeof(int), 1, f_hist_node);
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			fwrite(&(node->cur_profile_prob_info->counts[i]), sizeof(double), 1, f_hist_node);
		} // i loop.
	}
	else
	{
		// Do not load the probability info.
		fwrite(&(no_prob_info_marker), sizeof(int), 1, f_hist_node);
	}

	// Dump the subnodes.
	for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
	{
		dump_histogram_node(node->cur_profile_nodes[i], f_hist_node);
	} // i loop.
}

t_histogram_node* load_multi_d_histogram(char* multi_d_fp)
{
	FILE* f_multi_d = open_f(multi_d_fp, "rb");
	t_histogram_node* main_node = load_histogram_node(f_multi_d);
	fclose(f_multi_d);

	return(main_node);
}

t_histogram_node* load_one_d_histogram_node_per_text_linear_histogram(char* hist_fp)
{
	double hist_min = 1000*1000;
	double hist_max = -1000*1000;

	double cur_val = -1;
	double cur_count;
	double cur_prob;

	// Get the min and max.
	FILE* f_hist = open_f(hist_fp, "r");
	while(1)
	{
		double new_val = 0;
		if(fscanf(f_hist, "%lf %lf %lf", &new_val, &cur_count, &cur_prob) != 3)
		{
			break;
		}

		if(cur_val > 0 && 
			new_val != cur_val+1)
		{
			fprintf(stderr, "The values are not consecutive in the histogram file %s\n", hist_fp);
			exit(0);
		}
		else
		{
			cur_val = new_val;
		}

		if(hist_min > cur_val)
		{
			hist_min = cur_val;
		}

		if(hist_max < cur_val)
		{
			hist_max = cur_val;
		}
	} // file reading loop.
	fclose(f_hist);

	if(hist_max < hist_min)
	{
		fprintf(stderr, "Histogram minimum is larger than histogram max: %lf, %lf\n", hist_min, hist_max);
		exit(0);
	}

	t_histogram_node* cur_hist = new t_histogram_node();
	cur_hist->cur_prof_min = hist_min;
	cur_hist->cur_prof_max = hist_max;
	cur_hist->cur_profile_nodes = new t_histogram_node*[(int)hist_max - (int)hist_min + 1];
	cur_hist->cur_profile_nodes -= (int)hist_min;
	cur_hist->cur_profile_prob_info = new t_prob_info();
	cur_hist->cur_profile_prob_info->counts = new double[(int)hist_max - (int)hist_min + 1];
	cur_hist->cur_profile_prob_info->counts -= (int)hist_min;
	cur_hist->cur_profile_prob_info->probs = new double[(int)hist_max - (int)hist_min + 1];
	cur_hist->cur_profile_prob_info->probs -= (int)hist_min;

	for(int i = hist_min; i <= hist_max; i++)
	{
		cur_hist->cur_profile_nodes[i] = new t_histogram_node();
		cur_hist->cur_profile_nodes[i]->cur_profile_prob_info = NULL;
		cur_hist->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
		cur_hist->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
	} // i loop.

	f_hist = open_f(hist_fp, "r");
	while(1)
	{
		double cur_val;
		double cur_count;
		double cur_prob;
		if(fscanf(f_hist, "%lf %lf %lf", &cur_val, &cur_count, &cur_prob) != 3)
		{
			break;
		}

		// Set the count.
		cur_hist->cur_profile_prob_info->counts[(int)cur_val] = xlog(cur_count);
	} // file reading loop.
	fclose(f_hist);

	return(cur_hist);
}

t_histogram_node* load_histogram_node(FILE* f_hist_node)
{
	t_histogram_node* node = new t_histogram_node();
	node->smoothed_total_counts = 0.0;
	node->cur_profile_nodes = NULL;
	node->cur_profile_prob_info = NULL;

	// Dump the information for the current node.
	//fread(&(node->dim), sizeof(int), 1, f_hist_node);
	node->dim = read_bin_int(f_hist_node);
	//fread(&(node->cur_prof_min), sizeof(int), 1, f_hist_node);
	node->cur_prof_min = read_bin_int(f_hist_node);
	//fread(&(node->cur_prof_max), sizeof(int), 1, f_hist_node);
	node->cur_prof_max = read_bin_int(f_hist_node);

	int cur_prob_info_marker;
	//fread(&(cur_prob_info_marker), sizeof(int), 1, f_hist_node);
	cur_prob_info_marker = read_bin_int(f_hist_node);

	// Dump the counts and probabilities if they exist.
	if(cur_prob_info_marker == prob_info_marker)
	{
		//fwrite(&(prob_info_marker), sizeof(int), 1, f_hist_node);
		node->cur_profile_prob_info = new t_prob_info();
		node->cur_profile_prob_info->counts = new double[node->cur_prof_max - node->cur_prof_min + 1];
		node->cur_profile_prob_info->counts -= node->cur_prof_min;
		node->cur_profile_prob_info->probs = new double[node->cur_prof_max - node->cur_prof_min + 1];
		node->cur_profile_prob_info->probs -= node->cur_prof_min;
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			//fwrite(&(node->cur_profile_prob_info->counts[i]), sizeof(double), 1, f_hist_node);
			//fread(&(node->cur_profile_prob_info->counts[i]), sizeof(double), 1, f_hist_node);
			node->cur_profile_prob_info->counts[i] = read_bin_double(f_hist_node);
			node->cur_profile_prob_info->probs[i] = 0.0;
		} // i loop.
	}

	if(node->cur_prof_max >= node->cur_prof_min)
	{
		node->cur_profile_nodes = new t_histogram_node*[node->cur_prof_max - node->cur_prof_min + 1];
		node->cur_profile_nodes -= node->cur_prof_min;
	}
	else
	{
	}

	// Dump the subnodes.
	for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
	{
		//dump_histogram_node(node->cur_profile_nodes[i], f_hist_node);
		node->cur_profile_nodes[i] = load_histogram_node(f_hist_node);
	} // i loop.

	return(node);
}

// Sample from a histogram.
double* sample_multi_d_histogram(t_histogram_node* main_node, int n_dims, t_rng* rng)
{
	double* sampled_vals = new double[n_dims];

	// Sample the current histogram node.
	sample_histogram_node(main_node, n_dims, rng, sampled_vals);

	return(sampled_vals);
}

void sample_histogram_node(t_histogram_node* main_node, int n_dims, t_rng* rng, double* cur_sampled)
{
	// Get the counts.
	double* cur_node_vals_per_node = new double[main_node->cur_prof_max - main_node->cur_prof_min + 1];
	if(main_node->cur_profile_prob_info != NULL)
	{
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			cur_node_vals_per_node[i - main_node->cur_prof_min] = main_node->cur_profile_prob_info->probs[i];
		} // i loop.

		// Get the sampled index.
		int sampled_index = sample_index_per_log_probability_array(cur_node_vals_per_node, main_node->cur_prof_max - main_node->cur_prof_min + 1, rng);	

		cur_sampled[main_node->dim] = sampled_index + main_node->cur_prof_min;

		delete [] cur_node_vals_per_node;
	}
	else
	{
		for(int i = main_node->cur_prof_min; i <= main_node->cur_prof_max; i++)
		{
			cur_node_vals_per_node[i - main_node->cur_prof_min] = get_total_log_counts_per_multi_d_histogram(main_node->cur_profile_nodes[i]);
		} // i loop.

		// Get the sampled index.
		int sampled_index = sample_index_per_log_probability_array(cur_node_vals_per_node, main_node->cur_prof_max - main_node->cur_prof_min + 1, rng);

		cur_sampled[main_node->dim] = sampled_index + main_node->cur_prof_min;

		delete [] cur_node_vals_per_node;

		// Sample value in the next dimension.
		sample_histogram_node(main_node->cur_profile_nodes[sampled_index+main_node->cur_prof_min], n_dims, rng, cur_sampled);
	}
}

double get_log_gaussian_prob(double lin_mean, double lin_std_dev, double lin_val)
{
	double lin_prob = (1 / (lin_std_dev * pow(2 * my_PI, .5))) * exp(-1 * (lin_val - lin_mean) * (lin_val - lin_mean) / (lin_std_dev * lin_std_dev));
	return(xlog(lin_prob));
}

t_histogram_node* get_log_indexed_multi_d_histogram_per_linear_indexed_histogram(t_histogram_node* linear_indexed_hist_node, double log_base)
{
	t_histogram_node* log_indexed_hist = get_log_indexed_hist_node_per_linear_indexed_hist_node(linear_indexed_hist_node, log_base);

	return(log_indexed_hist);
}

void reset_counts_per_hist_node(t_histogram_node* hist_node, double log_val)
{
	if(hist_node->cur_profile_prob_info != NULL)
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			hist_node->cur_profile_prob_info->counts[i] = log_val;
		} // i loop.
	}
	else
	{
		for(int i = hist_node->cur_prof_min; i <= hist_node->cur_prof_max; i++)
		{
			reset_counts_per_hist_node(hist_node->cur_profile_nodes[i], log_val);
		} // i loop.
	}
}

t_histogram_node* get_log_indexed_hist_node_per_linear_indexed_hist_node(t_histogram_node* linear_indexed_hist_node, double log_base)
{
	if(linear_indexed_hist_node->cur_prof_max > linear_indexed_hist_node->cur_prof_min)
	{
		t_histogram_node* cur_log_indexed_node = new t_histogram_node();
		cur_log_indexed_node->cur_prof_min = xlog(linear_indexed_hist_node->cur_prof_min+1) / xlog(log_base);
		cur_log_indexed_node->cur_prof_max = xlog(linear_indexed_hist_node->cur_prof_max+1) / xlog(log_base);
		cur_log_indexed_node->dim = linear_indexed_hist_node->dim;

		cur_log_indexed_node->cur_profile_nodes = new t_histogram_node*[cur_log_indexed_node->cur_prof_max - cur_log_indexed_node->cur_prof_min + 1];

		// Initialize all the nodes to null first, then build them as needed.
		for(int i = linear_indexed_hist_node->cur_prof_min; i <= linear_indexed_hist_node->cur_prof_max; i++)
		{
			int cur_log_i = xlog(i+1) / xlog(log_base);
			cur_log_indexed_node->cur_profile_nodes[cur_log_i] = NULL;
		} // i loop.

		// Copy the probability info if it exists.
		if(linear_indexed_hist_node->cur_profile_prob_info != NULL)
		{
			cur_log_indexed_node->cur_profile_prob_info = new t_prob_info();
			cur_log_indexed_node->cur_profile_prob_info->counts = new double[cur_log_indexed_node->cur_prof_max - cur_log_indexed_node->cur_prof_min + 1];
			cur_log_indexed_node->cur_profile_prob_info->counts -= cur_log_indexed_node->cur_prof_min;
			cur_log_indexed_node->cur_profile_prob_info->probs = new double[cur_log_indexed_node->cur_prof_max - cur_log_indexed_node->cur_prof_min + 1];
			cur_log_indexed_node->cur_profile_prob_info->probs -= cur_log_indexed_node->cur_prof_min;

			// Set all the counts.
			for(int i = linear_indexed_hist_node->cur_prof_min; i <= linear_indexed_hist_node->cur_prof_max; i++)
			{
				int cur_log_i = xlog(i+1) / xlog(log_base);
				cur_log_indexed_node->cur_profile_prob_info->counts[cur_log_i] = xlog(0.0);
				cur_log_indexed_node->cur_profile_prob_info->probs[cur_log_i] = xlog(0.0);
			} // i loop.

			// Add all the counts. This effectively merges the counts that fall into the same log indexed bin.
			for(int i = linear_indexed_hist_node->cur_prof_min; i <= linear_indexed_hist_node->cur_prof_max; i++)
			{
				int cur_log_i = xlog(i+1);
				cur_log_indexed_node->cur_profile_prob_info->counts[cur_log_i] = xlog_sum(cur_log_indexed_node->cur_profile_prob_info->counts[cur_log_i], 
																							linear_indexed_hist_node->cur_profile_prob_info->counts[i]);
			} // i loop.
		}
		else
		{
			cur_log_indexed_node->cur_profile_prob_info = NULL;
		}
		
		// Add all the profile nodes.
		for(int i = linear_indexed_hist_node->cur_prof_min; i <= linear_indexed_hist_node->cur_prof_max; i++)
		{
			int cur_log_i = xlog(i+1) / xlog(log_base);

			// Make sure that this node is not already allocated yet.
			if(cur_log_indexed_node->cur_profile_nodes[cur_log_i] == NULL)
			{
				cur_log_indexed_node->cur_profile_nodes[cur_log_i] = get_log_indexed_hist_node_per_linear_indexed_hist_node(linear_indexed_hist_node->cur_profile_nodes[i], log_base);
			}
			else
			{
				// Merge the already exiting node.
				t_histogram_node* new_log_node = get_log_indexed_hist_node_per_linear_indexed_hist_node(linear_indexed_hist_node->cur_profile_nodes[i], log_base);
				t_histogram_node* new_merged_log_node = merge_multi_d_hist_per_counts(new_log_node, cur_log_indexed_node->cur_profile_nodes[cur_log_i]);

				delete_multi_d_histogram(new_log_node);
				delete_multi_d_histogram(cur_log_indexed_node->cur_profile_nodes[cur_log_i]);

				cur_log_indexed_node->cur_profile_nodes[cur_log_i] = new_merged_log_node;
			}
		} // i loop.

		return(cur_log_indexed_node);
	}
	else
	{
		t_histogram_node* cur_log_indexed_node = new t_histogram_node();
		cur_log_indexed_node->cur_prof_min = linear_indexed_hist_node->cur_prof_min;
		cur_log_indexed_node->cur_prof_max = linear_indexed_hist_node->cur_prof_max;
		cur_log_indexed_node->cur_profile_nodes = NULL;
		cur_log_indexed_node->cur_profile_prob_info = NULL;
		cur_log_indexed_node->dim = linear_indexed_hist_node->dim;

		return(cur_log_indexed_node);
	}
}

double get_mean_per_linear_indexed_histogram(t_histogram_node* one_d_hist_node)
{
	normalize_multi_d_hist_probs_per_hist_counts(one_d_hist_node);

	//double* probs = new double[one_d_hist_node->cur_prof_max - one_d_hist_node->cur_prof_min + 1];
	double cur_hist_mean = 0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		cur_hist_mean += (i * one_d_hist_node->cur_profile_prob_info->probs[i]);
	} // i loop.

	return(cur_hist_mean);
}

double get_std_dev_per_linear_indexed_histogram(t_histogram_node* one_d_hist_node)
{
	normalize_multi_d_hist_probs_per_hist_counts(one_d_hist_node);

	//double* probs = new double[one_d_hist_node->cur_prof_max - one_d_hist_node->cur_prof_min + 1];
	double cur_hist_mean = 0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		cur_hist_mean += (i * one_d_hist_node->cur_profile_prob_info->probs[i]);
	} // i loop.

	double cur_hist_var = 0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		cur_hist_var += ((i - cur_hist_mean) * (i - cur_hist_mean) * one_d_hist_node->cur_profile_prob_info->probs[i]);
	} // i loop.

	return(pow(cur_hist_var, .5));
}

/*
In log indexing, make sure that only valid log indices are used in computationa of mean/variance.
*/
double get_mean_per_log_indexed_histogram(t_histogram_node* one_d_hist_node, double log_base)
{
	normalize_multi_d_hist_probs_per_hist_counts(one_d_hist_node);

	double* mapped_probs = new double[one_d_hist_node->cur_prof_max - one_d_hist_node->cur_prof_min + 1];
	int i_mapped = 0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		if(is_this_a_valid_log(i, log_base))
		{
			mapped_probs[i_mapped] = one_d_hist_node->cur_profile_prob_info->probs[i];
			i_mapped++;
		}
	} // i loop.

	// Get the mean of the counts.
	double mapped_mean = 0;
	for(int i = 0; i < i_mapped; i++)
	{
		mapped_mean += (i * xexp(mapped_probs[i]));
	} // i loop.

	delete [] mapped_probs;

	mapped_mean = floor(mapped_mean);

	i_mapped = 0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		if(is_this_a_valid_log(i, log_base))
		{
			if(mapped_mean == i_mapped)
			{
				// Return the histogram index.
				return(i);
			}

			i_mapped++;
		}
	} // i loop.

	fprintf(stderr, "Could not find the mean position.\n");
	exit(0);
}

double get_std_dev_per_log_indexed_histogram(t_histogram_node* one_d_hist_node, double log_base)
{
	normalize_multi_d_hist_probs_per_hist_counts(one_d_hist_node);

	double hist_mean = get_mean_per_log_indexed_histogram(one_d_hist_node, log_base);

	double variance = 0.0;
	for(int i = one_d_hist_node->cur_prof_min; i <= one_d_hist_node->cur_prof_max; i++)
	{
		if(is_this_a_valid_log(i, log_base))
		{
			variance += (i - hist_mean) * (i - hist_mean) * xexp(one_d_hist_node->cur_profile_prob_info->probs[i]);
		}
	} // i loop.

	//delete [] counts;
	return(pow(variance, .5));
}
/*
Check if 
*/
bool is_this_a_valid_log(double val, double log_base)
{
	// First floorize the value.
	val = floor(val);

	double pow_val = pow(log_base, val);
	double next_pow_val = pow(log_base, val+1);

	if(floor(pow_val) != floor(next_pow_val))
	{
		return(true);
	}
	else
	{
		return(false);
	}

	////xlog(val+1)/xlog(log_base)

	////	pow(log_base, val)-1

	//double pow_val = pow(log_base, val)-1;
	//if(floor(pow_val) == pow_val)
	//{
	//	return(true);
	//}
	//else
	//{
	//	return(false);
	//}
}

t_histogram_node* get_uniform_multi_d_log_indexed_histogram_per_limits_per_dim(int* mins_per_dim, int* maxes_per_dim, int n_dims, double log_base)
{
	t_histogram_node* per_limit_uniform_log_indexed_hist_node = new t_histogram_node();
	per_limit_uniform_log_indexed_hist_node->cur_prof_min = 1000*1000;
	per_limit_uniform_log_indexed_hist_node->cur_prof_max = -1000*1000;
	per_limit_uniform_log_indexed_hist_node->dim = 0;
	get_uniform_log_indexed_histogram_node_per_limits_per_dim(per_limit_uniform_log_indexed_hist_node, mins_per_dim, maxes_per_dim, 0, n_dims, log_base);

	return(per_limit_uniform_log_indexed_hist_node);
}

void get_uniform_log_indexed_histogram_node_per_limits_per_dim(t_histogram_node* cur_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims, double log_base)
{
	if(i_dim == n_dims)
	{
		return;
	}

	// Set the limits for the current node.
	cur_node->cur_prof_min = mins_per_dim[i_dim];
	cur_node->cur_prof_max = maxes_per_dim[i_dim];
	cur_node->cur_profile_nodes = new t_histogram_node*[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
	cur_node->cur_profile_nodes -= cur_node->cur_prof_min;
	cur_node->cur_profile_prob_info = NULL;

	// Setup the probability information if necessary.
	if(i_dim == n_dims - 1)
	{
		cur_node->cur_profile_prob_info = new t_prob_info();
		cur_node->cur_profile_prob_info->counts = new double[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
		cur_node->cur_profile_prob_info->counts -= cur_node->cur_prof_min;

		cur_node->cur_profile_prob_info->probs = new double[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
		cur_node->cur_profile_prob_info->probs -= cur_node->cur_prof_min;

		for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
		{
			if(is_this_a_valid_log((double)i, log_base))
			{
				cur_node->cur_profile_prob_info->counts[i] = xlog(10.0);
				cur_node->cur_profile_prob_info->probs[i] = xlog(0.0);
			}
			else
			{
				cur_node->cur_profile_prob_info->counts[i] = xlog(0.0);
			}
		} // i loop.
	}
	
	// Go over all the next dimension nodes.
	for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
	{
		cur_node->cur_profile_nodes[i] = new t_histogram_node();
		cur_node->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
		cur_node->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
		cur_node->cur_profile_nodes[i]->dim = i_dim+1;

		// Get the next node.
		if(is_this_a_valid_log((double)i, log_base))
		{
			get_uniform_log_indexed_histogram_node_per_limits_per_dim(cur_node->cur_profile_nodes[i], mins_per_dim, maxes_per_dim, i_dim+1, n_dims, log_base);
		}
	} // i loop.
}

struct t_r_s_n
{
	double sig;
	int i;
};

bool sort_rsn_nodes_per_decreasing_signal(t_r_s_n* n1, t_r_s_n* n2)
{
	return(n1->sig < n2->sig);
}

void get_rank_signal(double* prof, int n_pts, double* rank_prof)
{
	vector<t_r_s_n*>* nodes = new vector<t_r_s_n*>();
	for (int i_pt = 0; i_pt < n_pts; i_pt++)
	{
		t_r_s_n* cur_node = new t_r_s_n();
		cur_node->i = i_pt;
		cur_node->sig = prof[i_pt];
		nodes->push_back(cur_node);
	} // i_pt loop.

	sort(nodes->begin(), nodes->end(), sort_rsn_nodes_per_decreasing_signal);

	for (int rank_i = 0; rank_i < n_pts; rank_i++)
	{
		rank_prof[nodes->at(rank_i)->i] = rank_i;
		delete nodes->at(rank_i);
	} // i_pt loop
	delete nodes;
}

void get_correlation(double* prof1, double* prof2, int l_win, double& cur_win_corr)
{
	double prof1_mean;
	double prof1_var;
	get_stats(prof1, l_win, prof1_mean, prof1_var);

	double prof2_mean;
	double prof2_var;
	get_stats(prof2, l_win, prof2_mean, prof2_var);

	double cur_cross_cov = 0.0;
	for(int i = 0; i < l_win; i++)
	{
		cur_cross_cov += (prof1[i]-prof1_mean) * (prof2[i]-prof2_mean);
	} // i loop.

	if(prof1_var > 0 && prof2_var > 0)
	{
		cur_win_corr = cur_cross_cov / ((l_win-1)*pow(prof1_var * prof2_var, .5));
	}
	else
	{
		cur_win_corr = 0.0;
	}
	//double curr_inner_prod = inner_product(prof1, prof2, l_win);
	//double l_vec1 = inner_product(prof1, prof1, l_win);
	//double l_vec2 = inner_product(prof2, prof2, l_win);

	//cur_win_corr = (curr_inner_prod / pow(l_vec1 * l_vec2, .5));
}

t_histogram_node* get_uniform_multi_d_histogram_per_limits_per_dim(int* mins_per_dim, int* maxes_per_dim, int n_dims)
{
	t_histogram_node* per_limit_uniform_hist_node = new t_histogram_node();
	per_limit_uniform_hist_node->cur_prof_min = 1000*1000;
	per_limit_uniform_hist_node->cur_prof_max = -1000*1000;
	per_limit_uniform_hist_node->dim = 0;
	get_uniform_histogram_node_per_limits_per_dim(per_limit_uniform_hist_node, mins_per_dim, maxes_per_dim, 0, n_dims);

	return(per_limit_uniform_hist_node);
}

void get_uniform_histogram_node_per_limits_per_dim(t_histogram_node* cur_node, int* mins_per_dim, int* maxes_per_dim, int i_dim, int n_dims)
{
	if(i_dim == n_dims)
	{
		return;
	}

	// Set the limits for the current node.
	cur_node->cur_prof_min = mins_per_dim[i_dim];
	cur_node->cur_prof_max = maxes_per_dim[i_dim];
	cur_node->cur_profile_nodes = new t_histogram_node*[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
	cur_node->cur_profile_nodes -= cur_node->cur_prof_min;
	cur_node->cur_profile_prob_info = NULL;

	// Setup the probability information if necessary.
	if(i_dim == n_dims - 1)
	{
		cur_node->cur_profile_prob_info = new t_prob_info();
		cur_node->cur_profile_prob_info->counts = new double[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
		cur_node->cur_profile_prob_info->counts -= cur_node->cur_prof_min;

		cur_node->cur_profile_prob_info->probs = new double[cur_node->cur_prof_max - cur_node->cur_prof_min + 1];
		cur_node->cur_profile_prob_info->probs -= cur_node->cur_prof_min;

		for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
		{
			cur_node->cur_profile_prob_info->counts[i] = xlog(10.0);
			cur_node->cur_profile_prob_info->probs[i] = xlog(0.0);
		} // i loop.
	}
	
	// Go over all the next dimension nodes.
	for(int i = mins_per_dim[i_dim]; i <= maxes_per_dim[i_dim]; i++)
	{
		cur_node->cur_profile_nodes[i] = new t_histogram_node();
		cur_node->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
		cur_node->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
		cur_node->cur_profile_nodes[i]->dim = i_dim+1;

		// Get the next node 
		get_uniform_histogram_node_per_limits_per_dim(cur_node->cur_profile_nodes[i], mins_per_dim, maxes_per_dim, i_dim+1, n_dims);
	} // i loop.
}

double* generate_odd_length_univariate_zero_mean_log_gaussian_dist(double l_bin, int n_pts)
{
	if(n_pts % 2 != 1)
	{
		fprintf(stderr, "Make sure the gaussian length is odd.\n");
		exit(0);
	}

	double* log_probs = new double[n_pts+2];

	int n_half_pts = (n_pts - 1) / 2;

	// Move the pointer to handle the indices.
	log_probs -= (-1*n_half_pts);

	// Generate the half.
	double cur_val = 0;
	log_probs[0] = get_log_gaussian_prob(0.0, 1.0, cur_val);
	cur_val += l_bin;
	for(int i_pt = 1; i_pt < n_half_pts; i_pt++)
	{
		log_probs[i_pt] = get_log_gaussian_prob(0, 1.0, cur_val);

		// Copy the negative value.
		log_probs[-1 * i_pt] = log_probs[i_pt];
		cur_val += l_bin;
	} // i_pt loop.
	
	// Move te pointer back.
	log_probs += (-1*n_half_pts);
	return(log_probs);
}

// Sample from a list of probabilities.
int sample_index_per_log_probability_array(double* probability_dist, int n_vals, t_rng* rng)
{
	// Get the total val.
	double total_log_counts = xlog(0.0);
	for(int i = 0; i < n_vals; i++)
	{
		total_log_counts = xlog_sum(total_log_counts, probability_dist[i]);
	} // i loop.

	//fprintf(stderr, "Probability sum: %lf\n", total_log_counts);

	double rand_normalized_val = xlog(rng->random_double_ran3()) + total_log_counts;

	double* cdf = new double[n_vals+1];
	cdf[0] = probability_dist[0];
	//fprintf(stderr, "CDFs:\n");
	for(int i = 1; i < n_vals; i++)
	{	
		cdf[i] = xlog_sum(cdf[i-1], probability_dist[i]);
		//fprintf(stderr, "%d: %lf\n", i, cdf[i]);
	} // i loop.

	for(int i = 0; i < n_vals; i++)
	{
		if(rand_normalized_val <= cdf[i])
		{
			delete [] cdf;
			return(i);
		}
	} // i loop.

	fprintf(stderr, "We are not supposed to be here: %lf, %lf\n", rand_normalized_val, cdf[n_vals-1]);
	exit(0);
	//return(n_vals-1);
}

// For an initialized histogram, count the data entries.
void count_data_per_histogram_node_per_one_indexed_data(t_histogram_node* main_node, double** profiles, int n_dims, int n_signal_wins)
{
	// At this point, all the next level nodes are allocated, we can do the counting.
	for(int i_sig = 1; i_sig <= n_signal_wins; i_sig++)
	{
		t_histogram_node* cur_node = main_node;

		// For 2 dimensions, this loop must be executed 1ce. For 3, 2.
		int cur_i_dim = 0;
		while(cur_node->cur_profile_prob_info == NULL)
		{
			if(cur_node->cur_prof_min > (int)(profiles[cur_i_dim][i_sig]) || 
				cur_node->cur_prof_max < (int)(profiles[cur_i_dim][i_sig]))
			{
				fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}

			if(cur_node->dim != cur_i_dim)
			{
				fprintf(stderr, "Fatal error: %s(%d)\n", __FILE__, __LINE__);
				exit(0);
			}

			cur_node = cur_node->cur_profile_nodes[(int)(profiles[cur_i_dim][i_sig])];
			cur_i_dim++;
		}

		// Sanity check on the dimension of the current node.
		if(cur_node->dim != n_dims-1)
		{
			fprintf(stderr, "The node to be updated is not at the lowest level: %d, %d\n", cur_node->dim, n_dims-1);
			exit(0);
		}

		//fprintf(stderr, "Found probs @ %d\n", i_dim);
		//getc(stdin);
/*
		for(int i_dim = 1; i_dim < n_dims; i_dim++)
		{
			if(cur_node == NULL)
			{
				fprintf(stderr, "WTF??\n");
				exit(0);
			}

			if((int)(profiles[i_dim][i_sig]) > cur_node->cur_prof_max ||
				(int)(profiles[i_dim][i_sig]) < cur_node->cur_prof_min)
			{
				fprintf(stderr, "WTF8??\n");
				exit(0);
			}
			cur_node = cur_node->next_level_nodes[(int)(profiles[i_dim][i_sig])];
		} // i_dim loop.
*/
		//fprintf(stderr, "Setting i_sig: %d (%lf)\n", i_sig, profiles[i_dim][i_sig]);

		// At this point, we have the node at the lowest profile.
		cur_node->cur_profile_prob_info->counts[(int)profiles[n_dims-1][i_sig]] = xlog_increment(cur_node->cur_profile_prob_info->counts[(int)profiles[n_dims-1][i_sig]]);

		//if(cur_node->cur_profile_prob_info->counts[(int)profiles[n_dims-1][i_sig]] > 0.0)
		//{
		//	bool profile_non_zero = false;
  //          for(int i_dim = 0; i_dim < n_dims; i_dim++)
  //          {
  //              if(profiles[i_dim][i_sig] > 4)
		//		{
		//			profile_non_zero = true;
		//		}
		//	}

		//	if(profile_non_zero)
		//	{
		//		fprintf(stderr, "%d. position has higher than 1:\n", i_sig);
		//		for(int i_dim = 0; i_dim < n_dims; i_dim++)
		//		{
		//			fprintf(stderr, "%d: %lf\n", i_dim, profiles[i_dim][i_sig]);
		//		}	
		//		getc(stdin);
		//	}
		//}
	} // i_sig loop.

	double total_counts = get_total_log_counts_per_multi_d_histogram(main_node);

	if (__DUMP_HISTOGRAM_MSGS__)
	{
		fprintf(stderr, "Total number of points: %lf (%lf), %d data points\n", total_counts, exp(total_counts), n_signal_wins);
	}
}

double get_total_log_counts_per_multi_d_histogram(t_histogram_node* hist)
{
	return(get_total_log_counts_per_node(hist));
}

double get_total_log_smoothed_counts_per_node(t_histogram_node* node)
{
	double total_log_counts = xlog(0.0);
	if(node == NULL)
	{
		return(total_log_counts);
	}

	if(node->cur_profile_prob_info != NULL && node->cur_profile_prob_info->counts != NULL)
	{
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			total_log_counts = xlog_sum(total_log_counts, node->cur_profile_prob_info->counts[i]);
		} // i loop.
	}
	else
	{
		if(node->cur_prof_min <= node->cur_prof_max)
		{
			for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
			{
				total_log_counts = xlog_sum(total_log_counts, get_total_log_smoothed_counts_per_node(node->cur_profile_nodes[i]));
			} // i loop.
		}
		else
		{
			total_log_counts = xlog_sum(total_log_counts, node->smoothed_total_counts);
		}
	} 

	return(total_log_counts);
}

double get_total_log_counts_per_node(t_histogram_node* node)
{
	double total_log_counts = xlog(0.0);
	if(node == NULL)
	{
		return(total_log_counts);
	}

	if(node->cur_profile_prob_info != NULL && node->cur_profile_prob_info->counts != NULL)
	{
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			total_log_counts = xlog_sum(total_log_counts, node->cur_profile_prob_info->counts[i]);
		} // i loop.
	}
	else
	{
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			total_log_counts = xlog_sum(total_log_counts, get_total_log_counts_per_node(node->cur_profile_nodes[i]));
		} // i loop.
	} 

	return(total_log_counts);
}

t_histogram_node* copy_histogram_node(t_histogram_node* node)
{
	if(node == NULL)
	{
		return(NULL);
	}

	t_histogram_node* new_node = new t_histogram_node();
	new_node->dim = node->dim;
	new_node->cur_prof_min = node->cur_prof_min;
	new_node->cur_prof_max = node->cur_prof_max;
	new_node->cur_profile_prob_info = NULL;
	
	if(node->cur_prof_max >= node->cur_prof_min)
	{
		new_node->cur_profile_nodes = new t_histogram_node*[new_node->cur_prof_max - new_node->cur_prof_min + 1];
		new_node->cur_profile_nodes -= new_node->cur_prof_min;

		// Copy all the next level nodes recursively.
		for(int i = new_node->cur_prof_min; i <= new_node->cur_prof_max; i++)
		{
			// Copy all the next level histogram nodes.
			new_node->cur_profile_nodes[i] = copy_histogram_node(node->cur_profile_nodes[i]);
		} // i loop.

		// Copy all the probabilities.
		if(node->cur_profile_prob_info != NULL)
		{
			new_node->cur_profile_prob_info = new t_prob_info(); 
			new_node->cur_profile_prob_info->counts = new double[new_node->cur_prof_max - new_node->cur_prof_min + 1];
			new_node->cur_profile_prob_info->counts -= new_node->cur_prof_min;
			new_node->cur_profile_prob_info->probs = new double[new_node->cur_prof_max - new_node->cur_prof_min + 1];
			new_node->cur_profile_prob_info->probs -= new_node->cur_prof_min;
			for(int i = new_node->cur_prof_min; i <= new_node->cur_prof_max; i++)
			{
				new_node->cur_profile_prob_info->probs[i] = node->cur_profile_prob_info->probs[i];
				new_node->cur_profile_prob_info->counts[i] = node->cur_profile_prob_info->counts[i];
			} // i loop.
		}
		else
		{
			new_node->cur_profile_prob_info = NULL;
		}
	}

	return(new_node);
}

t_histogram_node* init_multi_d_histogram(t_histogram_node* init_hist)
{
	return(copy_histogram_node(init_hist));
}

void normalize_multi_d_hist_probs_per_hist_counts(t_histogram_node* node)
{
	double normalizing_log_factor = get_total_log_counts_per_node(node);
	normalize_histogram_node_probs_per_counts(node, normalizing_log_factor);
}

void normalize_multi_d_hist_probs_per_smoothed_hist_counts(t_histogram_node* node)
{
	double normalizing_log_factor = get_total_log_counts_per_node(node);
	normalize_histogram_node_probs_per_counts(node, normalizing_log_factor);
}

//void smoothing_check(t_histogram_node* node)
//{
//	double total_smoothed_per_counting = get_total_log_smoothed_counts_per_node(node);
//	double total_smoothed_count = node->smoothed_total_counts;
//
//	fprintf(stderr, "%lf, %lf\n", total_smoothed_per_counting, total_smoothed_count);
//	exit(0);
//}

void normalize_histogram_node_probs_per_counts(t_histogram_node* node, double normalizing_log_factor)
{
	if(node == NULL)
	{
		return;
	}

	for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
	{
		normalize_histogram_node_probs_per_counts(node->cur_profile_nodes[i], normalizing_log_factor);
	} //

	// If the prob. info exists, normalize it.
	if(node->cur_profile_prob_info != NULL)
	{
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			node->cur_profile_prob_info->probs[i] = xlog_div(node->cur_profile_prob_info->counts[i], normalizing_log_factor);
		}
	}
}

double get_probability(t_histogram_node* hist, int n_dims, double* data)
{
	int cur_dim = 0;
	t_histogram_node* cur_node = hist;
	while(cur_node->cur_profile_prob_info == NULL)
	{
		if(cur_node == NULL)
		{
			fprintf(stderr, "Traced node is null, dimension node valid.\n");
			exit(0);
		}
		else if((int)(data[cur_dim]) >= cur_node->cur_prof_min &&
			(int)(data[cur_dim]) <= cur_node->cur_prof_max)
		{
			cur_node = cur_node->cur_profile_nodes[(int)(data[cur_dim])];
			cur_dim++;
		}
		else
		{
			fprintf(stderr, "Dimension not valid.\n");
			exit(0);
		}
	} // cur_node loop.

	if(cur_dim != n_dims-1)
	{
		fprintf(stderr, "Dimension after node tracing is not consistent: %d, %d\n", cur_dim, n_dims);
		exit(0);
	}

	if((int)(data[cur_dim]) >= cur_node->cur_prof_min &&
			(int)(data[cur_dim]) <= cur_node->cur_prof_max)
	{
		double cur_prob = cur_node->cur_profile_prob_info->probs[(int)(data[cur_dim])];
		if(cur_prob > xlog(0.0))
		{
			return(cur_prob);
		}
		else
		{
			return(xlog(0.0));
			//fprintf(stderr, "Zero probability!\n");
			//exit(0);
		}
	}
	else
	{
		fprintf(stderr, "Non existing dimension!\n");
		exit(0);
	}
}

void delete_histogram_node(t_histogram_node* node)
{
	// Delete probability info if it exists.
	if(node->cur_profile_prob_info != NULL)
	{
		delete [] (node->cur_profile_prob_info->counts + node->cur_prof_min);
		delete [] (node->cur_profile_prob_info->probs + node->cur_prof_min);

		delete node->cur_profile_prob_info;
	}

	// Delete the profile nodes.
	if(node->cur_prof_max >= node->cur_prof_min)
	{
		for(int i = node->cur_prof_min; i <= node->cur_prof_max; i++)
		{
			delete_histogram_node(node->cur_profile_nodes[i]);
		} // i loop.

		// Delete the list.
		delete [] (node->cur_profile_nodes + node->cur_prof_min);
	}

	delete node;
}

void delete_multi_d_histogram(t_histogram_node* hist)
{
	delete_histogram_node(hist);
}

// Merge two histograms.
t_histogram_node* merge_multi_d_hist_per_counts(t_histogram_node* hist1, t_histogram_node* hist2)
{
	return(merge_histogram_nodes_per_counts(hist1, hist2));
}

t_histogram_node* merge_histogram_nodes_per_counts(t_histogram_node* node1, t_histogram_node* node2)
{
	if(node1->dim != node2->dim)
	{
		fprintf(stderr, "Merging two nodes at different dimensions: %d, %d\n", node1->dim, node2->dim);
		exit(0);
	}

	int merged_node_min = MIN(node1->cur_prof_min, node2->cur_prof_min);
	int merged_node_max = MAX(node1->cur_prof_max, node2->cur_prof_max);
	t_histogram_node* merged_hist = new t_histogram_node();
	merged_hist->cur_prof_max = merged_node_max;
	merged_hist->cur_prof_min = merged_node_min;
	merged_hist->cur_profile_prob_info = NULL;
	merged_hist->dim = node1->dim;

	if(merged_node_max >= merged_node_min)
	{	
		merged_hist->cur_profile_nodes = new t_histogram_node*[merged_node_max - merged_node_min + 1];
		merged_hist->cur_profile_nodes -= merged_node_min;
	}

	// Copy the nodes.
	for(int i = merged_node_min; i <= merged_node_max; i++)
	{
		if(i >= node1->cur_prof_min && i <= node1->cur_prof_max)
		{
			if(i >= node2->cur_prof_min && i <= node2->cur_prof_max)
			{
				// Both matches: Merge the nodes, recursively.
				t_histogram_node* cur_i_merged_hist = merge_multi_d_hist_per_counts(node1->cur_profile_nodes[i], node2->cur_profile_nodes[i]);
				merged_hist->cur_profile_nodes[i] = cur_i_merged_hist;
			}
			else
			{
				// The index is valid for only histogram 1, copy the histogram from histogram 1.
				merged_hist->cur_profile_nodes[i] = init_multi_d_histogram(node1->cur_profile_nodes[i]);
			}
		} // i matches node1 coordinates.
		else
		{
			if(i >= node2->cur_prof_min && i <= node2->cur_prof_max)
			{
				merged_hist->cur_profile_nodes[i] = init_multi_d_histogram(node2->cur_profile_nodes[i]);
			}
			else
			{
				// No matches, this happens when the limits are not overlapping. This merges the limits. In this case, we need to setup an empty node. 
				// Empty node adding is ok since this point should never be encountered.
				merged_hist->cur_profile_nodes[i] = new t_histogram_node();
				merged_hist->cur_profile_nodes[i]->cur_prof_min = 1000*1000;
				merged_hist->cur_profile_nodes[i]->cur_prof_max = -1000*1000;
				merged_hist->cur_profile_nodes[i]->cur_profile_prob_info = NULL;
				merged_hist->cur_profile_nodes[i]->dim = merged_hist->dim+1;
				//fprintf(stderr, "The index does not match with both node1 and node2.\n");
				//exit(0);
			}
		} // i does not match node1 coordinates.
	} // i loop.

	// Update the probability information: Copy all the values.
	if(node1->cur_profile_prob_info != NULL ||
		node2->cur_profile_prob_info != NULL)
	{
		merged_hist->cur_profile_prob_info = new t_prob_info();
		merged_hist->cur_profile_prob_info->counts = new double[merged_node_max - merged_node_min + 1];
		merged_hist->cur_profile_prob_info->counts -= merged_node_min;
		merged_hist->cur_profile_prob_info->probs = new double[merged_node_max - merged_node_min + 1];
		merged_hist->cur_profile_prob_info->probs -= merged_node_min;

		// Initialize the counts/probs.
		for(int i = merged_node_min; i <= merged_node_max; i++)
		{
			merged_hist->cur_profile_prob_info->counts[i] = xlog(0.0);
			merged_hist->cur_profile_prob_info->probs[i] = xlog(0.0);

			if(node1->cur_profile_prob_info != NULL &&
				i >= node1->cur_prof_min && i <= node1->cur_prof_max)
			{
				merged_hist->cur_profile_prob_info->counts[i] = xlog_sum(merged_hist->cur_profile_prob_info->counts[i], node1->cur_profile_prob_info->counts[i]);
			}
			
			if(node2->cur_profile_prob_info != NULL &&
				i >= node2->cur_prof_min && i <= node2->cur_prof_max)
			{
				merged_hist->cur_profile_prob_info->counts[i] = xlog_sum(merged_hist->cur_profile_prob_info->counts[i], node2->cur_profile_prob_info->counts[i]);
			}
		} // i loop.
	}

	return(merged_hist);
}


void dump_transition_matrix(t_transition_matrix* transition_matrix, int n_states, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	fwrite(&n_states, sizeof(int), 1, f_op);

	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			//fprintf(stderr, "%lf ", transition_matrix->state_transition_counts[i_st][j_st]);
			fwrite(&(transition_matrix->state_transition_counts[i_st][j_st]), sizeof(double), 1, f_op);
		} // j_st loop.

		//fprintf(stderr, "\n");
	} // i_st loop.

	//fprintf(stderr, "Transition Probabilities:\n");
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			//fprintf(stderr, "%lf ", transition_matrix->state_transition_probs[i_st][j_st]);
			fwrite(&(transition_matrix->state_transition_probs[i_st][j_st]), sizeof(double), 1, f_op);
		} // j_st loop.

		//fprintf(stderr, "\n");
	} // i_st loop.

	fclose(f_op);
}

t_transition_matrix* load_transition_matrix(int n_states, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "rb");
	//int _n_states;
	//fread(&_n_states, sizeof(int), 1, f_op);
	//_n_states = read_bin_double(f_op);
	read_bin_double(f_op);
	t_transition_matrix* transition_matrix = init_transition_matrix(NULL, n_states);

	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			//fprintf(stderr, "%lf ", transition_matrix->state_transition_counts[i_st][j_st]);
			//fwrite(&(transition_matrix->state_transition_counts[i_st][j_st]), sizeof(double), 1, f_op);
			//fread(&(transition_matrix->state_transition_counts[i_st][j_st]), sizeof(double), 1, f_op);
			transition_matrix->state_transition_counts[i_st][j_st] = read_bin_double(f_op);
		} // j_st loop.

		//fprintf(stderr, "\n");
	} // i_st loop.

	//fprintf(stderr, "Transition Probabilities:\n");
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			//fprintf(stderr, "%lf ", transition_matrix->state_transition_probs[i_st][j_st]);
			//fwrite(&(transition_matrix->state_transition_probs[i_st][j_st]), sizeof(double), 1, f_op);
			//fread(&(transition_matrix->state_transition_probs[i_st][j_st]), sizeof(double), 1, f_op);
			transition_matrix->state_transition_probs[i_st][j_st] = read_bin_double(f_op);
		} // j_st loop.

		//fprintf(stderr, "\n");
	} // i_st loop.

	fclose(f_op);

	return(transition_matrix);
}

void dump_transition_matrix(t_transition_matrix* transition_matrix, int n_states)
{
	fprintf(stderr, "Transition Counts:\n");
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			fprintf(stderr, "%lf ", transition_matrix->state_transition_counts[i_st][j_st]);
		} // j_st loop.

		fprintf(stderr, "\n");
	} // i_st loop.

	fprintf(stderr, "Transition Probabilities:\n");
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			fprintf(stderr, "%lf ", transition_matrix->state_transition_probs[i_st][j_st]);
		} // j_st loop.

		fprintf(stderr, "\n");
	} // i_st loop.
}

// Get the mode of the distribution.
void get_mode_online(double* signal, int l_signal, double& mode_posn, int& n_vals_at_mode)
{
	bool* processed = new bool[l_signal+2];
	memset(processed, 0, sizeof(bool) * (l_signal+2));

	double cur_max_count = 0.0;
	mode_posn = 0.0;
	n_vals_at_mode = 0;

	for(int i = 1; i <= l_signal; i++)
	{
		if(!processed[i])
		{
			//double cur_bin_val = signal[i];

			int cur_bin_posn_cnt = 0;
			for(int j = i; j <= l_signal; j++)
			{
				if(!processed[j] && 
					signal[j] == signal[i])
				{
					cur_bin_posn_cnt++;
					processed[j] = true;
				}
			} // j loop.

			if(cur_bin_posn_cnt > cur_max_count)
			{
				cur_max_count = cur_bin_posn_cnt;
				mode_posn = signal[i];
				n_vals_at_mode = cur_bin_posn_cnt;
			} // max_prob check.

			processed[i] = true;
		} // processed check for current bin posn.
		else
		{
			// This value is processed.
		}
	} // i loop.

	fprintf(stderr, "Mode: %lf, Count: %d\n", mode_posn, n_vals_at_mode);
} // get_mode_online function.

void copy_transition_matrix(t_transition_matrix* dest_state_transition_matrix, t_transition_matrix* init_state_transition_matrix, int n_states)
{
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			dest_state_transition_matrix->state_transition_probs[i_st][j_st] = init_state_transition_matrix->state_transition_probs[i_st][j_st];
			dest_state_transition_matrix->state_transition_counts[i_st][j_st] = init_state_transition_matrix->state_transition_counts[i_st][j_st];
		} // j_st loop.
	} // i_st loop.
}

t_transition_matrix* init_transition_matrix(t_transition_matrix* init_state_transition_counts, int n_states)
{
	t_transition_matrix* trans_matrix = new t_transition_matrix();
	trans_matrix->state_transition_counts = new double*[n_states+1];
	trans_matrix->state_transition_probs = new double*[n_states+1];

	for(int i_st = 0; i_st < n_states; i_st++)
	{
		trans_matrix->state_transition_counts[i_st] = new double[n_states];
		trans_matrix->state_transition_probs[i_st] = new double[n_states];

		for(int j_st = 0; j_st < n_states; j_st++)
		{
			if(init_state_transition_counts != NULL)
			{
				trans_matrix->state_transition_probs[i_st][j_st] = init_state_transition_counts->state_transition_probs[i_st][j_st];
				trans_matrix->state_transition_counts[i_st][j_st] = init_state_transition_counts->state_transition_counts[i_st][j_st];
			}
			else
			{
				trans_matrix->state_transition_probs[i_st][j_st] = 1.0 / ((double)n_states);
				trans_matrix->state_transition_counts[i_st][j_st] = 1;
			}			
		} // i_st2 loop.
	} // i_st loop.

	return(trans_matrix);
}

void delete_transition_matrix(t_transition_matrix* cur_transition_probs, int n_states)
{
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		delete [] cur_transition_probs->state_transition_counts[i_st];
		delete [] cur_transition_probs->state_transition_probs[i_st];
	} // i_st loop.

	delete [] cur_transition_probs->state_transition_counts;
	delete [] cur_transition_probs->state_transition_probs;

	delete cur_transition_probs;
}

t_transition_matrix* get_transition_matrix_per_transition_counts(double** state_transition_counts, int n_states)
{
	t_transition_matrix* trans_matrix = new t_transition_matrix();
	trans_matrix->state_transition_counts = new double*[n_states+1];
	trans_matrix->state_transition_probs = new double*[n_states+1];

	for(int i_st = 0; i_st < n_states; i_st++)
	{
		trans_matrix->state_transition_counts[i_st] = new double[n_states];
		trans_matrix->state_transition_probs[i_st] = new double[n_states];

		for(int j_st = 0; j_st < n_states; j_st++)
		{
			trans_matrix->state_transition_probs[i_st][j_st] = 0.0;
			trans_matrix->state_transition_counts[i_st][j_st] = 0;
		} // i_st2 loop.
	} // i_st loop.

	// Update the emission and transition stats.
	fprintf(stderr, "Transition statistics:\n");
	for(int i_st = 0; i_st < n_states; i_st++)
	{
		// Get the total counts and update the counts.
		double cur_state_cnts = 0;
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			trans_matrix->state_transition_counts[i_st][j_st] = state_transition_counts[i_st][j_st];
			cur_state_cnts += state_transition_counts[i_st][j_st];
		} // j loop.

		// Update the transition matrix based on the counts.
		for(int j_st = 0; j_st < n_states; j_st++)
		{
			trans_matrix->state_transition_probs[i_st][j_st] = state_transition_counts[i_st][j_st] / cur_state_cnts;
		} // j loop.
	} // i loop.

	return(trans_matrix);
}

