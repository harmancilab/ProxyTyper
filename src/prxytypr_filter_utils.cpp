#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "filter_utils.h"
#include "../../../lib/utils/file/utils.h"
#include "../../../lib/utils/rng/rng.h"
#include "../../../lib/utils/ansi_string/ansi_string.h"
#include "../../../lib/genomics_coords/genomics_coords.h"
#include "../../../lib/genomics_utils/annotation/annot_region_tools.h"
#include "../../../lib/genomics_utils/signal_track/signal_track_tools.h"
#include "../../../lib/genomics_utils/signal_processing/min_max_utils.h"
#include <algorithm>
#include "../../../lib/utils/xmath/log/xlog_math.h"
#include "string.h"

using namespace std;

bool __DUMP_FILTER_MSGS__ = false;

#ifdef __unix__
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_sort.h>
	//#include <gsl/gsl_wavelet.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_fft_complex.h>

	#define REAL(z,i) ((z)[2*(i)])
	#define IMAG(z,i) ((z)[2*(i)+1])
#endif

struct t_ext_double_array
{
	double* buffer;	
	int n_vals;
	int buffer_length;
};

t_ext_double_array* alloc_ext_array()
{
	t_ext_double_array* ext_array = new t_ext_double_array();
	ext_array->buffer = new double[1000];
	ext_array->buffer_length = 1000;
	ext_array->n_vals = 0;

	return(ext_array);
}

void anisotropic_diffusion_filter(double* signal, int l_signal, 
	int n_iterations,
	double delta_t,
	double kappa,
	double (get_C)(double, double))
{
	double nablaW;
	double nablaE;
	double cW;
	double cE;
	for(int i_iter = 0; i_iter < n_iterations; i_iter++)
	{
		if(i_iter % 1000 == 0)
		{
			fprintf(stderr, "%d. iteration.           \r", i_iter);
		}

		// Compute current nabla arrays.
		for(int i = 0; i < l_signal; i++)
		{
			if(i > 0)
			{
				nablaW = signal[i-1] - signal[i];
			}
			else
			{
				nablaW = -1 * signal[i];
			}

			if(i < l_signal-1)
			{
				nablaE = signal[i+1] - signal[i];
			}
			else
			{
				nablaE = -1 * signal[i];
			}

			// Compute conductivities.
			cW = get_C(nablaW, kappa);
			cE = get_C(nablaE, kappa);

			signal[i] += delta_t * (nablaW * cW + nablaE * cE);
		} // i loop.

		// Dump the signal.
if(__DUMP_FILTER_MSGS__)
{
		char cur_iter_data_fp[1000];
		sprintf(cur_iter_data_fp, "diff_data_%d.txt", i_iter);
		FILE* f_cur_iter = open_f(cur_iter_data_fp, "w");
		for(int i = 0; i < l_signal; i++)
		{
			fprintf(f_cur_iter, "%.3f ", signal[i]);
		} // i loop.
		fclose(f_cur_iter);
}
	} // i_iter loop.
}

double get_exp_C(double nabla, double kappa)
{
	return(exp(-(nabla/kappa) * (nabla/kappa)));
}

double get_frac_C(double nabla, double kappa)
{
	return(1 / (1 + (nabla/kappa) * (nabla/kappa)));
}

double* get_per_posn_median_gradient_profile_per_win(double* signal_profile, int l_profile, int l_gradient_win)
{
	if(l_gradient_win % 2 != 1)
	{
		l_gradient_win++;
	}

	int half_l_gradient_win = l_gradient_win / 2;
	int min_val_int = -10000;
	int max_val_int = 10000;
	int* left_pdf = new int[max_val_int - min_val_int + 2];
	memset(left_pdf, 0, sizeof(int) * (max_val_int - min_val_int + 2));
	left_pdf -= min_val_int;
	int* right_pdf = new int[max_val_int - min_val_int + 2];
	right_pdf -= min_val_int;
	memset(right_pdf, 0, sizeof(int) * (max_val_int - min_val_int + 2));

	int cur_win_mid = 1;
	int cur_win_start = 1;
	int cur_win_end = 1 + half_l_gradient_win;
	double cur_gradient = 0;

	// Fill the left pdf.
	int cur_win_min = 10000;
	int cur_win_max = -10000;
	for(int i = cur_win_start; i <= cur_win_mid; i++)
	{
		if(cur_win_min > (int)signal_profile[i])
		{
			cur_win_min = (int)signal_profile[i];
		}

		if(cur_win_max < (int)signal_profile[i])
		{
			cur_win_max = (int)signal_profile[i];
		}

		int cur_val = (int)floor(signal_profile[i]);
		left_pdf[cur_val]++;
	} // i loop.

	// Fill the right pdf.
	for(int i = cur_win_mid; i <= cur_win_end; i++)
	{
		if(cur_win_min > (int)signal_profile[i])
		{
			cur_win_min = (int)signal_profile[i];
		}

		if(cur_win_max < (int)signal_profile[i])
		{
			cur_win_max = (int)signal_profile[i];
		}

		int cur_val = (int)floor(signal_profile[i]);
		right_pdf[cur_val]++;
	} // i loop.

	int prev_win_start = cur_win_start;
	int prev_win_end = cur_win_end;
	int prev_win_mid = cur_win_mid;

	// Do a simple right to left gradient estimation per position.
	double* gradient_profile = new double[l_profile + 3];
	memset(gradient_profile, 0, sizeof(double) * (l_profile+2));
	for(cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
	{
		// Current gradient is the right half signal minus left half signal.
		cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
		cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);
		
		// Remove the values that left the window.
		for(int i = prev_win_start; i < cur_win_start; i++)
		{
			left_pdf[(int)(signal_profile[i])]--;
		} // i loop.

		for(int i = prev_win_mid+1; i <= cur_win_mid; i++)
		{
			left_pdf[(int)(signal_profile[i])]++;
		} // i loop.

		// Update the right pdf values.
		for(int i = prev_win_mid; i < cur_win_mid; i++)
		{
			right_pdf[(int)(signal_profile[i])]--;
		} // i loop.

		// Add the values that entered the window on the right.
		for(int i = prev_win_end+1; i <= cur_win_end; i++)
		{
			if(cur_win_min > (int)signal_profile[i])
			{
				cur_win_min = (int)signal_profile[i];
			}

			if(cur_win_max < (int)signal_profile[i])
			{
				cur_win_max = (int)signal_profile[i];
			}

			right_pdf[(int)(signal_profile[i])]++;
		} // i loop.

		int left_pdf_median = 0;
		int cur_left_pdf_total = 0;
		for(int i = cur_win_min; i <= cur_win_max; i++)
		{
			cur_left_pdf_total += left_pdf[i];
			if(cur_left_pdf_total >= l_gradient_win/4)
			{
				left_pdf_median = i;
				break;
			}
		} // i loop.

		int right_pdf_median = 0;
		int cur_right_pdf_total = 0;
		for(int i = cur_win_min; i <= cur_win_max; i++)
		{
			cur_right_pdf_total += right_pdf[i];
			if(cur_right_pdf_total >= l_gradient_win/4)
			{
				right_pdf_median = i;
				break;
			}
		} // i loop.

		// Update the profile
		gradient_profile[cur_win_mid] = right_pdf_median - left_pdf_median;

		// Update the previous window coordinates.
		prev_win_start = cur_win_start;
		prev_win_end = cur_win_end;

		// Update the previous window middle point
		prev_win_mid = cur_win_mid;
	} // cur_win_mid loop.

	return(gradient_profile);
}

double* get_per_posn_mean_gradient_profile_per_win(double* signal_profile, int l_profile, int l_gradient_win)
{
	if(l_gradient_win % 2 != 1)
	{
		l_gradient_win++;
	}

	int half_l_gradient_win = l_gradient_win / 2;

	int cur_win_mid = 1;
	int cur_win_start = 1;
	int cur_win_end = 1 + half_l_gradient_win;
	double cur_gradient = 0;

	for(int i = cur_win_start; i <= cur_win_end; i++)
	{
		if(i < cur_win_mid)
		{
			cur_gradient -= signal_profile[i];
		}
		else
		{
			cur_gradient += signal_profile[i];
		}	
	} // i loop.

	int prev_win_start = cur_win_start;
	int prev_win_end = cur_win_end;
	int prev_win_mid = cur_win_mid;

	// Do a simple right to left gradient estimation per position.
	double* gradient_profile = new double[l_profile + 2];
	memset(gradient_profile, 0, sizeof(double) * l_profile);
	for(cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
	{
		// Current gradient is the right half signal minus left half signal.
		cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
		cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);
		
		// Remove the values that left the window.
		for(int i = prev_win_start; i < cur_win_start; i++)
		{
			cur_gradient += signal_profile[i];
		} // i loop.

		// The middle value in the current window also moved.
		cur_gradient -= signal_profile[prev_win_mid];
		cur_gradient -= signal_profile[prev_win_mid];

		// Add the values that entered the window on the right.
		for(int i = prev_win_end+1; i <= cur_win_end; i++)
		{
			cur_gradient += signal_profile[i];
		} // i loop.

		// Update the profile
		gradient_profile[cur_win_mid] = cur_gradient / (half_l_gradient_win);

		//// Following is a sanity check.
		//double check_gradient = 0;
		//for(int i = cur_win_start; i <= cur_win_end; i++)
		//{
		//	if(i < cur_win_mid)
		//	{
		//		check_gradient -= signal_profile[i];
		//	}
		//	else
		//	{
		//		check_gradient += signal_profile[i];
		//	}	
		//} // i loop.

		//if(fabs(check_gradient - cur_gradient) > 0.1)
		//{
		//	fprintf(stderr, "Gradients do not match @ %d\n", cur_win_mid);
		//	exit(0);
		//}

		// Update the previous window coordinates.
		prev_win_start = cur_win_start;
		prev_win_end = cur_win_end;
		prev_win_mid = cur_win_mid;
	} // i loop.

	return(gradient_profile);
}

vector<int>* get_local_gradient_abs_extremas(double* signal_profile, int l_profile, int l_win)
{
	// Compute the extrema.
	// Get the extrema regions for the current filtered regions.
	vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
	vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
	int* derivative_map = new int[l_profile+2];
	memset(derivative_map, 0, sizeof(int) * (l_profile+2));
	get_extrema_per_plateaus(signal_profile, l_profile, maxima, minima, derivative_map, 0);

	sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);
	sort(minima->begin(), minima->end(), sort_extremas_per_posn);

	vector<int>* extrema_extremas = new vector<int>();

	vector<int>* maxima_extremas = get_local_gradient_abs_extremas(maxima, signal_profile, l_profile, l_win, true);
	for(int i = 0; i < maxima_extremas->size(); i++)
	{
		extrema_extremas->push_back(maxima->at(maxima_extremas->at(i))->extrema_posn);
	}
	vector<int>* minima_extremas = get_local_gradient_abs_extremas(minima, signal_profile, l_profile, l_win, false);
	for(int i = 0; i < minima_extremas->size(); i++)
	{
		extrema_extremas->push_back(minima->at(minima_extremas->at(i))->extrema_posn);
	}

	delete_extrema_nodes(maxima);
	delete_extrema_nodes(minima);

	return(extrema_extremas);
}

vector<int>* get_local_gradient_abs_extremas(vector<t_extrema_node*>* all_extrema, double* signal_profile, int l_profile, int l_win, bool process_maxima)
{
	bool* posn_is_maximum = new bool[all_extrema->size() + 2];
	memset(posn_is_maximum, 0, all_extrema->size());

	int cur_extrema_start_i = 0;
	int cur_extrema_end_i = 0;
	int ext_i = 0;

	vector<int>* gradient_abs_maxima_posns = new vector<int>();
	while(cur_extrema_start_i < all_extrema->size())
	{
		// Find the end extrema position.
		int cur_extrema_end_i = cur_extrema_start_i;
		while(cur_extrema_end_i < all_extrema->size() && 
			all_extrema->at(cur_extrema_end_i)->extrema_posn <= (all_extrema->at(cur_extrema_start_i)->extrema_posn+l_win))
		{
			cur_extrema_end_i++;
		} //  ext_i loop.

		if(cur_extrema_end_i >= all_extrema->size())
		{
			break;
		}

		// Find the maximum position and set its flag.
		int climax_posn_i = -1;
		double climax_posn_abs_height = (process_maxima)?(-1000):(1000);
		for(int ext_i = cur_extrema_start_i; ext_i <= cur_extrema_end_i; ext_i++)
		{
			if(process_maxima && all_extrema->at(ext_i)->height_at_extrema > climax_posn_abs_height ||
				!process_maxima && all_extrema->at(ext_i)->height_at_extrema < climax_posn_abs_height)
			{
				climax_posn_abs_height = (all_extrema->at(ext_i)->height_at_extrema);

				// Set this position as the local absolute maximum.
				climax_posn_i = ext_i;
				posn_is_maximum[ext_i] = true;
				climax_posn_abs_height = (all_extrema->at(ext_i)->height_at_extrema);
			}
		} // ext_i loop.

		cur_extrema_start_i++;
	} // cur_extrema_start_i loop.

	cur_extrema_start_i = 0;
	while(cur_extrema_start_i < all_extrema->size())
	{
		// Find the end extrema position.
		int cur_extrema_end_i = cur_extrema_start_i;
		while(cur_extrema_end_i < all_extrema->size() && 
			all_extrema->at(cur_extrema_end_i)->extrema_posn <= (all_extrema->at(cur_extrema_start_i)->extrema_posn+l_win))
		{
			cur_extrema_end_i++;
		} //  ext_i loop.

		if(cur_extrema_end_i >= all_extrema->size())
		{
			break;
		}

		// Find the maximum position and set it's flag.
		double climax_posn_abs_height = (process_maxima)?(-1000):(1000);
		for(int ext_i = cur_extrema_start_i; ext_i <= cur_extrema_end_i; ext_i++)
		{
			if(posn_is_maximum[ext_i] && 
				(process_maxima && all_extrema->at(ext_i)->height_at_extrema > climax_posn_abs_height ||
				!process_maxima && all_extrema->at(ext_i)->height_at_extrema < climax_posn_abs_height))
			{
				// Set this position as the local absolute maximum.
				climax_posn_abs_height = (all_extrema->at(ext_i)->height_at_extrema);
			}
		} // ext_i loop.

		// Following loop unsets the positions that have are not locally optimal.
		for(int ext_i = cur_extrema_start_i; ext_i <= cur_extrema_end_i; ext_i++)
		{
			if(posn_is_maximum[ext_i] && 
				(process_maxima && all_extrema->at(ext_i)->height_at_extrema < climax_posn_abs_height ||
				!process_maxima && all_extrema->at(ext_i)->height_at_extrema > climax_posn_abs_height))
			{
				posn_is_maximum[ext_i] = false;
			}
		} // ext_i loop.

		cur_extrema_start_i++;
	} // cur_extrema_start_i loop.

	for(int ext_i = 0; ext_i < all_extrema->size(); ext_i++)
	{
		if(posn_is_maximum[ext_i] &&
			(process_maxima && all_extrema->at(ext_i)->height_at_extrema > 0 ||
			!process_maxima && all_extrema->at(ext_i)->height_at_extrema < 0))
		{
			gradient_abs_maxima_posns->push_back(ext_i);
		}
	} // ext_i loop.

	delete [] posn_is_maximum;

	return(gradient_abs_maxima_posns);
}

double* get_contribution_weight_profile(double* signal_profile, int l_profile, int l_gradient_win)
{
	double* gradient_profile = get_per_posn_mean_gradient_profile_per_win(signal_profile, l_profile, l_gradient_win);

	double max_abs_gradient = 0;
	int max_grad_posn = 0;
	for(int i = 1; i <= l_profile; i++)
	{
		if(fabs(gradient_profile[i]) > max_abs_gradient)
		{
			max_abs_gradient = fabs(gradient_profile[i]);
			max_grad_posn = i;
		}
	} // i loop.

	fprintf(stderr, "Maximum absolute gradient is %lf (%d, %lf)\n", max_abs_gradient, max_grad_posn, signal_profile[max_grad_posn]);

	// Transform the per position gradient profile.
	double* contribution_weight_profile = new double[l_profile + 2];
	
	for(int i = 1; i <= l_profile; i++)
	{
		contribution_weight_profile[i] = get_exp_C(gradient_profile[i], max_abs_gradient/2);
		//contribution_weight_profile[i] = get_frac_C(gradient_profile[i], max_abs_gradient);
	} // i loop.

	delete [] gradient_profile;
	return(contribution_weight_profile);
}

double* get_multiscale_max_contribution_weight_per_posn(double* signal_profile, int l_profile, int l_start_win, int l_end_win, double log_step)
{
	double cur_scale_l_win = (double)l_start_win;

	double* max_contribution_per_posn = new double[l_profile + 10];
	memset(max_contribution_per_posn, 0, sizeof(double) * (l_profile+2));

	int n_scales = 0;

	while(cur_scale_l_win <= l_end_win)
	{
		fprintf(stderr, "Processing gradient @ %lf\n", cur_scale_l_win);

		double* cur_scale_gradient = get_per_posn_mean_gradient_profile_per_win(signal_profile, l_profile, (int)cur_scale_l_win);

		double max_abs_gradient = 0;
		for(int i = 1; i <= l_profile; i++)
		{
			if(fabs(cur_scale_gradient[i]) > max_abs_gradient)
			{
				max_abs_gradient = fabs(cur_scale_gradient[i]);
			}
		} // i loop.
		fprintf(stderr, "Maximum absolute gradient for current scale is %lf\n", max_abs_gradient);

		for(int i = 1; i <= l_profile; i++)
		{
			double cur_posn_contrib = get_exp_C(cur_scale_gradient[i], max_abs_gradient/2);

			if(cur_posn_contrib < 0)
			{
				cur_posn_contrib = 0;
			}

			max_contribution_per_posn[i] = MAX(max_contribution_per_posn[i], cur_posn_contrib);
			//max_contribution_per_posn[i] = MIN(max_contribution_per_posn[i], cur_posn_contrib);
			//max_contribution_per_posn[i] *= pow(cur_posn_contrib, .5);
		} // i loop.

if(__DUMP_FILTER_MSGS__)
{
		char cur_gradient_profile_fp[1000];
		sprintf(cur_gradient_profile_fp, "gradient_l_win_%lf.bin", cur_scale_l_win);
		dump_per_nucleotide_binary_profile(cur_scale_gradient, l_profile, cur_gradient_profile_fp);
}

		delete [] cur_scale_gradient;

		cur_scale_l_win *= log_step;
	} // cur_scale loop.

	return(max_contribution_per_posn);
}

void iterative_contribution_weighted_median_filter(double* signal_profile, int l_profile, int n_iterations, int l_gradient_win)
{
	fprintf(stderr, "THIS FUNCTION IS OBSOLETE, USE iterative_recursive_multiscale_scale_contribution_weighted_mean_filter(...)\n");
	exit(1);

	int half_l_gradient_win = l_gradient_win / 2;

	for(int cur_win_mid = 1; cur_win_mid <= 115229903-10000; cur_win_mid++)
	{
		signal_profile[cur_win_mid] = 0;
	} // cur_win_mid loop.

	for(int cur_win_mid = 115332220+10000; cur_win_mid <= l_profile; cur_win_mid++)
	{
		signal_profile[cur_win_mid] = 0;
	} // cur_win_mid loop.

	double* cur_contrib_pdf = new double[20000+2];
	memset(cur_contrib_pdf, 0, 20000 * sizeof(double));
	cur_contrib_pdf -= -10000;

	double* cur_iter_smoothed_signal = new double[l_profile + 2];

	for(int i_iter = 0; i_iter < n_iterations; i_iter++)
	{
		fprintf(stderr, "Iteration %d\n", i_iter);

		double* cur_contrib_weight_profile = get_contribution_weight_profile(signal_profile, l_profile, l_gradient_win);		

		// For each position, sum up the contribution weights.
		for(int cur_win_mid = 115229903; cur_win_mid <= 115332220; cur_win_mid++)
		{
			int cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
			int cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);
			double cur_cumulative_weight = 1;

			double total_contribution = 0;
			int min_val_int = 10000;
			int max_val_int = -10000;

			// Go from middle to ends to simulate the end effects.
			for(int i = cur_win_mid; i >= cur_win_start; i--)
			{
				///cur_cumulative_weight *= cur_contrib_weight_profile[i] * cur_contrib_weight_profile[i];
				cur_cumulative_weight *= cur_contrib_weight_profile[i];

				cur_contrib_pdf[int(signal_profile[i])] += cur_cumulative_weight;

				if(min_val_int > (int)(signal_profile[i]))
				{
					min_val_int = (int)(signal_profile[i]);
				}

				if(max_val_int < (int)(signal_profile[i]))
				{
					max_val_int = (int)(signal_profile[i]);
				}

				total_contribution += cur_cumulative_weight;

				if(cur_cumulative_weight < 0.000001)
				{
					break;
				}
			} // i loop.

			// Update the cumulative weight: This is necessary since we are setting a new start from the middle.
			cur_cumulative_weight = cur_contrib_weight_profile[cur_win_mid];
			for(int i = cur_win_mid+1; i <= cur_win_end; i++)
			{
				//cur_cumulative_weight *= cur_contrib_weight_profile[i] * cur_contrib_weight_profile[i];
				cur_cumulative_weight *= cur_contrib_weight_profile[i];

				// Update the contribution pdf 
				cur_contrib_pdf[int(signal_profile[i])] += cur_cumulative_weight;

				if(min_val_int > (int)(signal_profile[i]))
				{
					min_val_int = (int)(signal_profile[i]);
				}

				if(max_val_int < (int)(signal_profile[i]))
				{
					max_val_int = (int)(signal_profile[i]);
				}

				total_contribution += cur_cumulative_weight;

				if(cur_cumulative_weight < 0.000001)
				{
					break;
				}
			} // i loop.

			// Get the median of the contribution pdf.
			double cur_total = 0;
			for(int i = min_val_int; i <= max_val_int; i++)
			{
				cur_total += cur_contrib_pdf[i];
				if(cur_total >= total_contribution / 2)
				{
					cur_iter_smoothed_signal[cur_win_mid] = i;
					break;
				}
			} // i loop.

			// Reset the pdf.
			for(int i = MAX(min_val_int-10, -10000); i <= MIN(max_val_int+10, 10000); i++)
			{
				cur_contrib_pdf[i] = 0;
			} // i loop.
		} // cur_win_mid loop.

		// Copy the smoothed signal.
		for(int i = 115229903; i <= 115332220; i++)
		{
			signal_profile[i] = cur_iter_smoothed_signal[i];
		} // i loop.

		if(i_iter % 5 == 0)
		{
			fprintf(stderr, "Dumping data\n");

			char cur_iter_signal_profile_fp[1000];
			sprintf(cur_iter_signal_profile_fp, "iter_%d.bin", i_iter);
			dump_per_nucleotide_binary_profile(signal_profile, l_profile, cur_iter_signal_profile_fp);

			char cur_iter_contribution_weight_profile_fp[1000];
			sprintf(cur_iter_contribution_weight_profile_fp, "contribution_weight_profile_%d.bin", i_iter);

			// Dump the contribution weight profile.
			dump_per_nucleotide_binary_profile(cur_contrib_weight_profile, l_profile, cur_iter_contribution_weight_profile_fp);

			fprintf(stderr, "Done dumping data\n");
		}

		// Free memory.
		delete [] cur_contrib_weight_profile;
	} // i_iter loop.
}

void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step)
{
	for(int cur_win_mid = 1; cur_win_mid <= 113886644-10000; cur_win_mid++)
	{
		signal_profile[cur_win_mid] = 0;
	} // cur_win_mid loop.

	for(int cur_win_mid = 115523741+10000; cur_win_mid <= l_profile; cur_win_mid++)
	{
		signal_profile[cur_win_mid] = 0;
	} // cur_win_mid loop.

	double* cur_iter_smoothed_signal = new double[l_profile + 2];

	int i_iter = 0;
	double cur_scale = begin_scale;
	while(cur_scale <= end_scale)
	{
		fprintf(stderr, "Processing scale %d\n", (int)cur_scale);

		// Process the current scale's iterations.
		for(int cur_scale_i_iter = 0; cur_scale_i_iter < n_iters_per_scale; cur_scale_i_iter++)
		{
			fprintf(stderr, "Processing iteration %d\n", cur_scale_i_iter);

			int l_gradient_win = (int)cur_scale;
			int half_l_gradient_win = l_gradient_win/ 2;

			double* cur_contrib_weight_profile = get_contribution_weight_profile(signal_profile, l_profile, l_gradient_win);

			// For each position, sum up the contribution weights.
			//for(int cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
			for(int cur_win_mid = 113886644; cur_win_mid <= 115523741; cur_win_mid++)
			{
				int cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				int cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);
				double cur_weighted_sum = 0;
				double cur_total_weight = 0;
				double cur_cumulative_weight = 1;

				// Go from middle to ends to simulate the end effects.
				for(int i = cur_win_mid; i >= cur_win_start; i--)
				{
					cur_cumulative_weight *= cur_contrib_weight_profile[i];

					cur_weighted_sum += signal_profile[i] * cur_cumulative_weight;
					cur_total_weight += cur_cumulative_weight;

					if(cur_cumulative_weight < 0.001)
					{
						break;
					}
				} // i loop.

				// Update the cumulative weight: This is necessary since we are setting a new start from the middle.
				cur_cumulative_weight = cur_contrib_weight_profile[cur_win_mid];
				for(int i = cur_win_mid+1; i <= cur_win_end; i++)
				{
					cur_cumulative_weight *= cur_contrib_weight_profile[i];

					cur_weighted_sum += signal_profile[i] * cur_cumulative_weight;
					cur_total_weight += cur_cumulative_weight;

					if(cur_cumulative_weight < 0.001)
					{
						break;
					}
				} // i loop.

				cur_iter_smoothed_signal[cur_win_mid] = cur_weighted_sum / cur_total_weight;
			} // i loop.

			// Copy the smoothed signal.
			for(int i = 1; i <= l_profile; i++)
			{
				signal_profile[i] = cur_iter_smoothed_signal[i];
			} // i loop.

			// Free memory.
			delete [] cur_contrib_weight_profile;
		} // cur_scale_i_iter loop.

		// Dump the current values.
		fprintf(stderr, "Dumping data for scale %lf\n", cur_scale);

		char cur_iter_signal_profile_fp[1000];
		sprintf(cur_iter_signal_profile_fp, "scale_%d.bin", (int)i_iter);
		dump_per_nucleotide_binary_profile(signal_profile, l_profile, cur_iter_signal_profile_fp);

		// Update the iteration number.
		i_iter++;

		// Updte the scale.
		cur_scale *= scale_step;
	} // cur_scale loop.
}

// An attempt to make the weighted filter fast.
void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter_fast(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step)
{
	fprintf(stderr, "THIS FUNCTION DOES NOT WORK CURRENTLY, NEED TO FIX THE UPDATES.\n");
	double* cur_iter_smoothed_signal = new double[l_profile + 2];

	int i_iter = 0;
	double cur_scale = begin_scale;
	while(cur_scale <= end_scale)
	{
		fprintf(stderr, "Processing scale %d\n", (int)cur_scale);

		// Process the current scale's iterations.
		for(int cur_scale_i_iter = 0; cur_scale_i_iter < n_iters_per_scale; cur_scale_i_iter++)
		{
			fprintf(stderr, "Processing iteration %d\n", cur_scale_i_iter);

			int l_gradient_win = (int)cur_scale;
			int half_l_gradient_win = l_gradient_win/ 2;

			double* cur_contrib_weight_profile = get_contribution_weight_profile(signal_profile, l_profile, l_gradient_win);

			for(int i = 0; i <= l_profile; i++)
			{
				if(cur_contrib_weight_profile[i] < 0.05)
				{
					cur_contrib_weight_profile[i] = 0.05;
				}
			}

			int cur_win_mid = 0;
			int cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
			int cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

			double cur_left_cumulative_weight = 1;
			double cur_left_sum = 0;
			double cur_left_total_weight = 0;

			int prev_win_start = cur_win_start;
			int prev_win_end = cur_win_end;
			int prev_win_mid = cur_win_mid;

			// Do the left computes: Slide right.
            double* cur_iter_left_total_weights = new double[l_profile + 2];
            memset(cur_iter_left_total_weights, 0, sizeof(double) * (l_profile+1));
			double* cur_iter_left_sums = new double[l_profile + 2];
            memset(cur_iter_left_sums, 0, sizeof(double) * (l_profile+1));
			for(cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
			{
				if(cur_win_mid % 1000000 == 0)
				{
					fprintf(stderr, "Processing left half window %d                \r", cur_win_mid);
				}

				cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

				// Process the full window unless the window is out of the bounds.
				bool window_out_of_bounds = (cur_win_mid <= l_profile - half_l_gradient_win - 10 && cur_win_mid >= half_l_gradient_win + 10);
				if(!window_out_of_bounds)
				{
					cur_left_cumulative_weight = 1;
					cur_left_sum = 0;
					cur_left_total_weight = 0;
					for(int i = cur_win_mid; i >= cur_win_start; i--)
					{
						cur_left_cumulative_weight *= cur_contrib_weight_profile[i];
						cur_left_sum += signal_profile[i] * cur_left_cumulative_weight;
						cur_left_total_weight += cur_left_cumulative_weight;
					} // i loop.
				}
				else
				{
					// Update the current left sum: 
					cur_left_total_weight -= cur_left_cumulative_weight;
					cur_left_total_weight *= cur_contrib_weight_profile[cur_win_mid];
					cur_left_total_weight += cur_contrib_weight_profile[cur_win_mid];

					cur_left_sum -= cur_left_cumulative_weight * signal_profile[prev_win_start];
					cur_left_sum *= cur_contrib_weight_profile[cur_win_mid];
					cur_left_sum += cur_contrib_weight_profile[cur_win_mid] * signal_profile[cur_win_mid];

					cur_left_cumulative_weight /= cur_contrib_weight_profile[prev_win_start];
					cur_left_cumulative_weight *= cur_contrib_weight_profile[cur_win_mid];
				}

				cur_iter_left_sums[cur_win_mid] = cur_left_sum;
				cur_iter_left_total_weights[cur_win_mid] = cur_left_total_weight;

				// Update the previous window coordinates.
				prev_win_mid = cur_win_mid;
				prev_win_start = cur_win_start;
				prev_win_end = cur_win_end;

#undef __LEFT_CHECK__
#ifdef __LEFT_CHECK__
				if(cur_win_mid % 100000 == 0)
				{
					// Go from middle to ends to simulate the end effects.
					double temp_cur_left_cumulative_weight = 1;
					double temp_cur_left_sum = 0;
					double temp_cur_left_total_weight = 0;
					for(int i = cur_win_mid; i >= cur_win_start; i--)
					{
						temp_cur_left_cumulative_weight *= cur_contrib_weight_profile[i];
						temp_cur_left_sum += signal_profile[i] * temp_cur_left_cumulative_weight;
						temp_cur_left_total_weight += temp_cur_left_cumulative_weight;
					} // i loop.

					 double weight_ratio = fabs(cur_left_total_weight - temp_cur_left_total_weight) / temp_cur_left_total_weight;
					 double sum_ratio = fabs(cur_left_sum - temp_cur_left_sum) / temp_cur_left_sum;
					fprintf(stderr, "%.15f, %.15f                   \r", weight_ratio, sum_ratio);

					if(weight_ratio > 0.1 ||
						sum_ratio > 0.1)
					{
						fprintf(stderr, "Check failed.\n");
						exit(0);
					}
				} // position check.
#endif // __LEFT_CHECK__
			} // cur_win_mid loop.

			double cur_right_cumulative_weight = 1;
			double cur_right_sum = 0;
			double cur_right_total_weight = 0;

			// Do the right computes: Do left sliding window to compute right weights.
            double* cur_iter_right_total_weights = new double[l_profile + 2];
            memset(cur_iter_right_total_weights, 0, sizeof(double) * (l_profile+1));
			double* cur_iter_right_sums = new double[l_profile + 2];
            memset(cur_iter_right_sums, 0, sizeof(double) * (l_profile+1));
			for(cur_win_mid = l_profile; cur_win_mid >= 1; cur_win_mid--)
			{
				if(cur_win_mid % 1000000 == 0)
				{
					fprintf(stderr, "Processing right half window %d                \r", cur_win_mid);
				}

				cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

				// Process the full window unless the window is out of the bounds.
				bool window_out_of_bounds = (cur_win_mid <= l_profile - half_l_gradient_win - 10 && cur_win_mid >= half_l_gradient_win + 10);
				if(!window_out_of_bounds)
				{
					cur_right_cumulative_weight = 1;
					cur_right_sum = 0;
					cur_right_total_weight = 0;
					for(int i = cur_win_mid+1; i <= cur_win_end; i++)
					{
						cur_right_cumulative_weight *= cur_contrib_weight_profile[i];
						cur_right_sum += signal_profile[i] * cur_right_cumulative_weight;
						cur_right_total_weight += cur_right_cumulative_weight;
					} // i loop.
				}
				else
				{
					cur_right_total_weight -= cur_right_cumulative_weight;
					cur_right_total_weight *= cur_contrib_weight_profile[cur_win_mid+1];
					cur_right_total_weight += cur_contrib_weight_profile[cur_win_mid+1];

					cur_right_sum -= cur_right_cumulative_weight * signal_profile[prev_win_end];
					cur_right_sum *= cur_contrib_weight_profile[cur_win_mid+1];
					cur_right_sum += cur_contrib_weight_profile[cur_win_mid+1] * signal_profile[cur_win_mid+1];

					cur_right_cumulative_weight /= cur_contrib_weight_profile[prev_win_end];
					cur_right_cumulative_weight *= cur_contrib_weight_profile[cur_win_mid+1];
				}

				cur_iter_right_total_weights[cur_win_mid] = cur_right_total_weight;
				cur_iter_right_sums[cur_win_mid] = cur_right_sum;

				prev_win_end = cur_win_end;
				prev_win_start = cur_win_start;
				prev_win_mid = cur_win_mid;

#undef __RIGHT_CHECK__
#ifdef __RIGHT_CHECK__
				if(cur_win_mid % 100000 == 0)
				{
					double temp_cur_right_cumulative_weight = 1;
					double temp_cur_right_sum = 0;
					double temp_cur_right_total_weight = 0;
					for(int i = cur_win_mid+1; i <= cur_win_end; i++)
					{
						temp_cur_right_cumulative_weight *= cur_contrib_weight_profile[i];
						temp_cur_right_sum += signal_profile[i] * temp_cur_right_cumulative_weight;
						temp_cur_right_total_weight += temp_cur_right_cumulative_weight;
					} // i loop.

					 double weight_ratio = fabs(cur_right_total_weight - temp_cur_right_total_weight) / temp_cur_right_total_weight;
					 double sum_ratio = fabs(cur_right_sum - temp_cur_right_sum) / temp_cur_right_sum;
					fprintf(stderr, "%.15f, %.15f                   \r", weight_ratio, sum_ratio);

					if(weight_ratio > 0.1 ||
							sum_ratio > 0.1)
					{
						fprintf(stderr, "Check failed.\n");
						exit(0);
					}
				} // position check.
#endif // __RIGHT_CHECK__
			} // cur_win_mid loop.

			// Compute the smoothed signal.
			for(int i = 1; i <= l_profile; i++)
			{
				signal_profile[i] = (cur_iter_right_sums[i] + cur_iter_left_sums[i]) / (cur_iter_left_total_weights[i] + cur_iter_right_total_weights[i]);
			} // i loop.

			// Free memory.
			delete [] cur_iter_right_sums;
			delete [] cur_iter_left_sums;
			delete [] cur_iter_left_total_weights;
			delete [] cur_iter_right_total_weights;

			// Free memory.
			delete [] cur_contrib_weight_profile;
		} // cur_scale_i_iter loop.

		// Dump the current values.
		fprintf(stderr, "Dumping data for scale %lf\n", cur_scale);

		char cur_iter_signal_profile_fp[1000];
		sprintf(cur_iter_signal_profile_fp, "scale_%d.bin", (int)i_iter);
		dump_per_nucleotide_binary_profile(signal_profile, l_profile, cur_iter_signal_profile_fp);

		// Update the iteration number.
		i_iter++;

		// Updte the scale.
		cur_scale *= scale_step;
	} // cur_scale loop.
}

/*
Following tests the log versus linear accuracy for reconstructing values.
	double step = 0.99;
	double* cumulative_step_vals = new double[5000];
	double* rand_vals = new double[5000];
	srand(0);
	double weighted_total = 0;
	double cur_cumulative_weight = 1;

	double log_weighted_total = xlog(0);
	double cur_log_cumulative_weight = 0;
	for(int i = 0; i < 5000; i++)
	{
		rand_vals[i] = (double)rand();

		cur_cumulative_weight *= step;
		weighted_total += cur_cumulative_weight * rand_vals[i];

		cur_log_cumulative_weight += xlog(step);
		log_weighted_total = xlog_sum(log_weighted_total, cur_log_cumulative_weight + log(rand_vals[i]));
	} // i loop.

	// Recover the last value from the log and linear values.
	cur_cumulative_weight = 1;
	cur_log_cumulative_weight = 0;
	for(int i = 0; i < 2000-1; i++)
	{
		fprintf(stderr, "Removing %d. value.\n", i);
		weighted_total /= step;
		weighted_total -= rand_vals[i];

		log_weighted_total -= xlog(step);
		log_weighted_total = xlog_sub(log_weighted_total, log(rand_vals[i]));
	} // i loop.

	fprintf(stderr, "Actual: %.15f\nLinear: %.15f\nLog: %.15f\n",
	rand_vals[1999],
	weighted_total / step,
	exp(log_weighted_total - xlog(step)));

	exit(0);
*/
void iterative_recursive_multiscale_scale_contribution_weighted_mean_filter_log_fast(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step)
{
	fprintf(stderr, "THIS FUNCTION IS NOT TESTED EXTENSIVELY, MAY NEED MORE TESTING, MAY HAVE UNSTABLE RESULTS.\n");
	double* cur_iter_smoothed_signal = new double[l_profile + 2];
	double* log_signal_profile = new double[l_profile + 2];

	for(int i = 1; i <= l_profile; i++)
	{
		log_signal_profile[i] = xlog(signal_profile[i]);
	} // i loop.

	int i_iter = 0;
	double cur_scale = begin_scale;
	while(cur_scale <= end_scale)
	{
		fprintf(stderr, "Processing scale %d\n", (int)cur_scale);

		// Process the current scale's iterations.
		for(int cur_scale_i_iter = 0; cur_scale_i_iter < n_iters_per_scale; cur_scale_i_iter++)
		{
			fprintf(stderr, "Processing iteration %d\n", cur_scale_i_iter);

			int l_gradient_win = (int)cur_scale;
			int half_l_gradient_win = l_gradient_win/ 2;

			double* cur_contrib_weight_profile = get_contribution_weight_profile(signal_profile, l_profile, l_gradient_win);

			for(int i = 0; i <= l_profile; i++)
			{
				if(cur_contrib_weight_profile[i] < 0.05)
				{
					cur_contrib_weight_profile[i] = 0.05;
				}

				cur_contrib_weight_profile[i] = xlog(cur_contrib_weight_profile[i]);
			}

			int cur_win_mid = 0;
			int cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
			int cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

			double cur_left_cumulative_weight = 1;
			double cur_left_sum = 0;
			double cur_left_total_weight = 0;

			int prev_win_start = cur_win_start;
			int prev_win_end = cur_win_end;
			int prev_win_mid = cur_win_mid;

			// Do the left computes: Slide right.
            double* cur_iter_left_total_weights = new double[l_profile + 2];
            memset(cur_iter_left_total_weights, 0, sizeof(double) * (l_profile+1));
			double* cur_iter_left_sums = new double[l_profile + 2];
            memset(cur_iter_left_sums, 0, sizeof(double) * (l_profile+1));
			for(cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
			{
				if(cur_win_mid % 1000000 == 0)
				{
					fprintf(stderr, "Processing left half window %d                \r", cur_win_mid);
				}

				cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

				// Process the full window unless the window is out of the bounds.
				bool window_out_of_bounds = (cur_win_mid <= l_profile - half_l_gradient_win - 10 && cur_win_mid >= half_l_gradient_win + 10);
				if(!window_out_of_bounds)
				{
					cur_left_cumulative_weight = 0;
					cur_left_sum = xlog(0);
					cur_left_total_weight = xlog(0);
					for(int i = cur_win_mid; i >= cur_win_start; i--)
					{
						cur_left_cumulative_weight += cur_contrib_weight_profile[i];
						cur_left_sum = xlog_sum(cur_left_sum, log_signal_profile[i] + cur_left_cumulative_weight);
						cur_left_total_weight = xlog_sum(cur_left_total_weight, cur_left_cumulative_weight);
					} // i loop.
				}
				else
				{
					// Update the current left sum: 
					if(cur_left_total_weight > cur_left_cumulative_weight)
						cur_left_total_weight = xlog_sub(cur_left_total_weight, cur_left_cumulative_weight);
					else
						cur_left_cumulative_weight = xlog(0);

					cur_left_total_weight += cur_contrib_weight_profile[cur_win_mid];
					cur_left_total_weight = xlog_sum(cur_left_total_weight, cur_contrib_weight_profile[cur_win_mid]);

					if(cur_left_sum > cur_left_cumulative_weight + log_signal_profile[prev_win_start])
						cur_left_sum = xlog_sub(cur_left_sum, cur_left_cumulative_weight + log_signal_profile[prev_win_start]);
					else
						cur_left_sum = xlog(0);
					cur_left_sum += cur_contrib_weight_profile[cur_win_mid];
					cur_left_sum = xlog_sum(cur_left_sum, cur_contrib_weight_profile[cur_win_mid] + log_signal_profile[cur_win_mid]);

					cur_left_cumulative_weight -= cur_contrib_weight_profile[prev_win_start];
					cur_left_cumulative_weight += cur_contrib_weight_profile[cur_win_mid];
				}

				cur_iter_left_sums[cur_win_mid] = cur_left_sum;
				cur_iter_left_total_weights[cur_win_mid] = cur_left_total_weight;

				// Update the previous window coordinates.
				prev_win_mid = cur_win_mid;
				prev_win_start = cur_win_start;
				prev_win_end = cur_win_end;

#define __LEFT_CHECK__
#ifdef __LEFT_CHECK__
				if(cur_win_mid % 10000 == 0)
				{
					// Go from middle to ends to simulate the end effects.
					double temp_cur_left_cumulative_weight = 1;
					double temp_cur_left_sum = 0;
					double temp_cur_left_total_weight = 0;
					for(int i = cur_win_mid; i >= cur_win_start; i--)
					{
						temp_cur_left_cumulative_weight *= exp(cur_contrib_weight_profile[i]);
						temp_cur_left_sum += signal_profile[i] * temp_cur_left_cumulative_weight;
						temp_cur_left_total_weight += temp_cur_left_cumulative_weight;
					} // i loop.

					if(temp_cur_left_sum > 1.0 && 
						temp_cur_left_sum / cur_left_total_weight > 1)
					{
						double weight_ratio = fabs(exp(cur_left_total_weight) - temp_cur_left_total_weight) / temp_cur_left_total_weight;
						double sum_ratio = fabs(exp(cur_left_sum) - temp_cur_left_sum) / temp_cur_left_sum;

						if(weight_ratio > 0.00001 ||
							sum_ratio > 0.00001)
						{
							fprintf(stderr, "Check failed @ %d:\n");
							fprintf(stderr, "%.15f, %.15f\n", weight_ratio, sum_ratio);
							fprintf(stderr, "Weights: %.15f, %.15f\nSums: %.15f, %.15f\n", exp(cur_left_total_weight), temp_cur_left_total_weight, 
								exp(cur_left_sum), temp_cur_left_sum);
							exit(0);
						}
					}
				} // position check.
#endif // __LEFT_CHECK__
			} // cur_win_mid loop.

			double cur_right_cumulative_weight = 0;
			double cur_right_sum = xlog(0);
			double cur_right_total_weight = xlog(0);

			// Do the right computes: Do left sliding window to compute right weights.
            double* cur_iter_right_total_weights = new double[l_profile + 2];
            memset(cur_iter_right_total_weights, 0, sizeof(double) * (l_profile+1));
			double* cur_iter_right_sums = new double[l_profile + 2];
            memset(cur_iter_right_sums, 0, sizeof(double) * (l_profile+1));
			for(cur_win_mid = l_profile; cur_win_mid >= 1; cur_win_mid--)
			{
				if(cur_win_mid % 1000000 == 0)
				{
					fprintf(stderr, "Processing right half window %d                \r", cur_win_mid);
				}

				cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);

				// Process the full window unless the window is out of the bounds.
				bool window_out_of_bounds = (cur_win_mid <= l_profile - half_l_gradient_win - 10 && cur_win_mid >= half_l_gradient_win + 10);
				if(!window_out_of_bounds)
				{
					cur_right_cumulative_weight = 0;
					cur_right_sum = xlog(0);
					cur_right_total_weight = xlog(0);
					for(int i = cur_win_mid+1; i <= cur_win_end; i++)
					{
						cur_right_cumulative_weight += cur_contrib_weight_profile[i];
						cur_right_sum = xlog_sum(cur_right_sum, log_signal_profile[i] + cur_right_cumulative_weight);
						cur_right_total_weight = xlog_sum(cur_right_total_weight, cur_right_cumulative_weight);
					} // i loop.
				}
				else
				{
					if(cur_right_total_weight > cur_right_cumulative_weight)
						cur_right_total_weight = xlog_sub(cur_right_total_weight, cur_right_cumulative_weight);
					else
						cur_right_total_weight = xlog(0);
					cur_right_total_weight += cur_contrib_weight_profile[cur_win_mid+1];
					cur_right_total_weight = xlog_sum(cur_right_total_weight, cur_contrib_weight_profile[cur_win_mid+1]);

					if(cur_right_sum >  cur_right_cumulative_weight + log_signal_profile[prev_win_end])
					{
						cur_right_sum = xlog_sub(cur_right_sum, cur_right_cumulative_weight + log_signal_profile[prev_win_end]);
					}
					else
					{
						cur_right_sum = xlog(0);
					}
					cur_right_sum += cur_contrib_weight_profile[cur_win_mid+1];
					cur_right_sum = xlog_sum(cur_right_sum, cur_contrib_weight_profile[cur_win_mid+1] + log_signal_profile[cur_win_mid+1]);

					cur_right_cumulative_weight -= cur_contrib_weight_profile[prev_win_end];
					cur_right_cumulative_weight += cur_contrib_weight_profile[cur_win_mid+1];
				}

				cur_iter_right_total_weights[cur_win_mid] = cur_right_total_weight;
				cur_iter_right_sums[cur_win_mid] = cur_right_sum;

				prev_win_end = cur_win_end;
				prev_win_start = cur_win_start;
				prev_win_mid = cur_win_mid;

#define __RIGHT_CHECK__
#ifdef __RIGHT_CHECK__
				if(cur_win_mid % 10000 == 0)
				{
					double temp_cur_right_cumulative_weight = 1;
					double temp_cur_right_sum = 0;
					double temp_cur_right_total_weight = 0;
					for(int i = cur_win_mid+1; i <= cur_win_end; i++)
					{
						temp_cur_right_cumulative_weight *= exp(cur_contrib_weight_profile[i]);
						temp_cur_right_sum += signal_profile[i] * temp_cur_right_cumulative_weight;
						temp_cur_right_total_weight += temp_cur_right_cumulative_weight;
					} // i loop.

					if(temp_cur_right_sum > 1.0 && 
						temp_cur_right_sum / cur_right_total_weight > 1)
					{
						double weight_ratio = fabs(exp(cur_right_total_weight) - temp_cur_right_total_weight) / temp_cur_right_total_weight;
						double sum_ratio = fabs(exp(cur_right_sum) - temp_cur_right_sum) / temp_cur_right_sum;

						if(weight_ratio > 0.00001 ||
							sum_ratio > 0.00001)
						{
							fprintf(stderr, "Check failed @ %d:\n");
							fprintf(stderr, "%.15f, %.15f\n", weight_ratio, sum_ratio);
							fprintf(stderr, "Weights: %.15f, %.15f\nSums: %.15f, %.15f\n", exp(cur_right_total_weight), temp_cur_right_total_weight, 
								exp(cur_right_sum), temp_cur_right_sum);
							exit(0);
						}
					}
				} // position check.
#endif // __RIGHT_CHECK__
			} // cur_win_mid loop.

			// Update the signal profile, and crecompute the log signal.
			fprintf(stderr, "Computing the signal from left and right components\n");
			for(int i = 1; i <= l_profile; i++)
			{
				signal_profile[i] = (exp(cur_iter_right_sums[i]) + exp(cur_iter_left_sums[i])) / (exp(cur_iter_left_total_weights[i]) + exp(cur_iter_right_total_weights[i]));
				log_signal_profile[i] = xlog(signal_profile[i]);
			} // i loop.
			
			// Free memory.
			delete [] cur_iter_right_sums;
			delete [] cur_iter_left_sums;
			delete [] cur_iter_left_total_weights;
			delete [] cur_iter_right_total_weights;

			// Free memory.
			delete [] cur_contrib_weight_profile;
		} // cur_scale_i_iter loop.

		// Dump the current values.
		fprintf(stderr, "Dumping data for scale %lf\n", cur_scale);

		char cur_iter_signal_profile_fp[1000];
		sprintf(cur_iter_signal_profile_fp, "scale_%d.bin", (int)i_iter);
		dump_per_nucleotide_binary_profile(signal_profile, l_profile, cur_iter_signal_profile_fp);

		// Update the iteration number.
		i_iter++;

		// Updte the scale.
		cur_scale *= scale_step;
	} // cur_scale loop.

	delete [] log_signal_profile;
}

void iterative_recursive_multiscale_scale_contribution_weighted_median_filter(double* signal_profile, int l_profile, 
	int n_iters_per_scale, 
	double begin_scale, double end_scale, double scale_step)
{
	fprintf(stderr, "THIS FUNCTION IS VERY SLOW FOR PRACTICAL PURPOSES.\n");

	for(int i = 10*1000*1000+1; i <= l_profile; i++)
	{
		signal_profile[i] = 0;
	} // i loop.

	double* cur_iter_smoothed_signal = new double[l_profile + 2];

	double* cur_win_pdf = new double[20000 + 2];
	memset(cur_win_pdf, 0, sizeof(double) * (20000 + 2));
	cur_win_pdf -= -10000;

	int i_iter = 0;
	double cur_scale = begin_scale;
	while(cur_scale <= end_scale)
	{
		fprintf(stderr, "Processing scale %d\n", (int)cur_scale);

		// Process the current scale's iterations.
		for(int cur_scale_i_iter = 0; cur_scale_i_iter < n_iters_per_scale; cur_scale_i_iter++)
		{
			fprintf(stderr, "Processing iteration %d\n", cur_scale_i_iter);

			int l_gradient_win = (int)cur_scale;
			int half_l_gradient_win = l_gradient_win/ 2;

			double* cur_contrib_weight_profile = get_contribution_weight_profile(signal_profile, l_profile, l_gradient_win);

			// For each position, sum up the contribution weights.
			//for(int cur_win_mid = 1; cur_win_mid <= l_profile; cur_win_mid++)
			for(int cur_win_mid = 1; cur_win_mid <= 10*1000*1000; cur_win_mid++)
			{
				if(cur_win_mid % 1000000 == 0)
				{
					fprintf(stderr, "Processing %d. value            \r", cur_win_mid);
				}

				int cur_win_start = MAX(1, cur_win_mid - half_l_gradient_win);
				int cur_win_end = MIN(l_profile, cur_win_mid + half_l_gradient_win);
				double cur_total_weight = 0;
				double cur_cumulative_weight = 1;

				int cur_win_min_int = 10000;
				int cur_win_max_int = -10000;

				// Go from middle to ends to simulate the end effects.
				for(int i = cur_win_mid; i >= cur_win_start; i--)
				{
					cur_cumulative_weight *= cur_contrib_weight_profile[i];

					cur_total_weight += cur_cumulative_weight;

					cur_win_pdf[(int)(signal_profile[i])] += cur_cumulative_weight;

					if(cur_win_min_int > (int)(signal_profile[i]))
					{
						cur_win_min_int = (int)(signal_profile[i]);
					}

					if(cur_win_max_int < (int)(signal_profile[i]))
					{
						cur_win_max_int = (int)(signal_profile[i]);
					}

					if(cur_cumulative_weight < 0.000001)
					{
						break;
					}
				} // i loop.

				// Update the cumulative weight: This is necessary since we are setting a new start from the middle.
				cur_cumulative_weight = 1;
				for(int i = cur_win_mid+1; i <= cur_win_end; i++)
				{
					cur_cumulative_weight *= cur_contrib_weight_profile[i];

					cur_total_weight += cur_cumulative_weight;

					cur_win_pdf[(int)(signal_profile[i])] += cur_cumulative_weight;

					if(cur_win_min_int > (int)(signal_profile[i]))
					{
						cur_win_min_int = (int)(signal_profile[i]);
					}

					if(cur_win_max_int < (int)(signal_profile[i]))
					{
						cur_win_max_int = (int)(signal_profile[i]);
					}

					if(cur_cumulative_weight < 0.000001)
					{
						break;
					}
				} // i loop.

				// Find the weighted median in the current window.
				double cur_sum = 0;
				for(int sig = cur_win_min_int; sig <= cur_win_max_int; sig++)
				{
					cur_sum += cur_win_pdf[sig];
					if(cur_sum >= cur_total_weight / 2)
					{
						cur_iter_smoothed_signal[cur_win_mid] = sig;
						break;
					}
				} // sig loop.

				cur_win_pdf += -10000;
				memset(cur_win_pdf, 0, sizeof(double) * (20000 + 2));
				cur_win_pdf -= -10000;
			} // i loop.

			// Copy the smoothed signal.
			for(int i = 1; i <= l_profile; i++)
			{
				signal_profile[i] = cur_iter_smoothed_signal[i];
			} // i loop.

			// Free memory.
			delete [] cur_contrib_weight_profile;
		} // cur_scale_i_iter loop.

		// Dump the current values.
		fprintf(stderr, "Dumping data for scale %lf\n", cur_scale);

		char cur_iter_signal_profile_fp[1000];
		sprintf(cur_iter_signal_profile_fp, "scale_%d.bin", (int)i_iter);
		dump_per_nucleotide_binary_profile(signal_profile, l_profile, cur_iter_signal_profile_fp);

		// Update the iteration number.
		i_iter++;

		// Updte the scale.
		cur_scale *= scale_step;
	} // cur_scale loop.
}

void add_val(t_ext_double_array* ext_array, double val)
{
	// Check if the buffer memory will be overwhelmed by adding the new value.
	if(ext_array->n_vals+1 >= ext_array->buffer_length)
	{
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Extending the array: %d -> %d\n", ext_array->buffer_length, ext_array->buffer_length * 2);
		int l_ext_ext_buffer = ext_array->buffer_length * 2;
		double* new_ext_buffer = new double [l_ext_ext_buffer];

		// Copy the values.
		for(int i = 0; i < ext_array->n_vals; i++)
		{
			new_ext_buffer[i] = ext_array->buffer[i];
		} // i loop.

		// Copy the new buffer the extended array.		
		delete [] ext_array->buffer;
		ext_array->buffer_length = l_ext_ext_buffer;
		ext_array->buffer = new_ext_buffer;
	}

	// Add the value and return.
	ext_array->buffer[ext_array->n_vals] = val;
	ext_array->n_vals++;
}

double* get_scaled_odd_length_rectangular(double scale, double sigma, double n_sigma_per_half_win, int& l_filter)
{
	t_ext_double_array* ext_array = alloc_ext_array();

	double win = n_sigma_per_half_win * scale * sigma;

	double delta_t = 1.0;
	int n_vals = 0;
	for(double t = -1.0*win; t <= win; t += delta_t)
	{
		n_vals++;
	} // t loop.

	if(n_vals % 2 == 0)
	{
		n_vals++;
	}

	for(double t = -1.0*win; t <= win; t += delta_t)
	{
		double cur_g_val = 1/(double)n_vals;

		// Add the value.
		add_val(ext_array, cur_g_val);
	} // t loop.

	double* rect_buffer = ext_array->buffer;
	delete(ext_array);

	return(rect_buffer);
}

const double _PI = 3.141592653589;
double* get_extended_odd_length_gaussian(double sigma, double scale, double n_sigma_per_half_win, int& l_filter)
{
	t_ext_double_array* ext_array = alloc_ext_array();
	double win = n_sigma_per_half_win * scale * sigma;

	double denom = scale * pow(2*_PI, .5) * sigma;
	int l_side = 0;
	int l_skipped = 0;
	//double delta_t = 0.006125;
	double delta_t = 0.0125;

	double total_filter = 0.0;
	for(double t = -1.0*win; t <= 0; t += delta_t)
	{
		double cur_sigma_squared = (sigma * sigma);
		double cur_g_val = (1/denom) * delta_t * exp(-1 * t * t / (2 * cur_sigma_squared));

		// Add the value.
		if(cur_g_val > pow(10.0, -11.0))
		{
			//fprintf(f_half_gaus, "%lf\n", cur_g_val);
			add_val(ext_array, cur_g_val);
			l_side++;
			total_filter += cur_g_val;
		}
		else
		{
			l_skipped++;
		}
	} // t loop.

	if((l_skipped + l_side + l_side) % 2 == 0)
	{
		l_skipped++;
	}
	
	l_filter = l_skipped + l_side + l_side;
	double max_val = ext_array->buffer[l_side-1];

	total_filter = total_filter + total_filter + l_skipped * max_val;

	double* extended_gauss = new double[l_filter];
	int i_arr = 0;
	for(int i = 0; i < l_side; i++)
	{
		extended_gauss[i_arr] = ext_array->buffer[i] / total_filter;
		i_arr++;
	}

	for(int i = 0; i < l_skipped; i++)
	{
		extended_gauss[i_arr] = max_val / total_filter;
		i_arr++;
	}

	for(int i = 0; i < l_side; i++)
	{
		extended_gauss[i_arr] = ext_array->buffer[l_side-i-1] / total_filter;
		i_arr++;
	}

	delete [] ext_array->buffer;
	delete ext_array;	

	return(extended_gauss);
}

double* get_scaled_odd_length_gaussian(double sigma, double scale, double n_sigma_per_half_win, int& l_filter)
{
	t_ext_double_array* ext_array = alloc_ext_array();

	double win = n_sigma_per_half_win * scale * sigma;

	double delta_t = 1.0;
	double denom = scale * pow(2*_PI, .5) * sigma;

	// Generate the mesh that corresponds to the half of the gaussian. This is used to generate the symmetric halves of the discrete gaussian filter.
	vector<double>* gauss_mesh = new vector<double>();
	for(double t = delta_t; t <= 1.0*win; t += delta_t)
	{
		gauss_mesh->push_back(t);
	} // t loop.

	// This adds the negative side: Go back over the mesh.
	for(int i_mesh = (int)gauss_mesh->size()-1; i_mesh >= 0; i_mesh--)
	{
		double t = gauss_mesh->at(i_mesh) * -1.0;
		double cur_sigma_squared = ((sigma * scale) * (sigma * scale));
		double cur_g_val = (1/denom) * delta_t * exp(-1 * t * t / (2 * cur_sigma_squared));

		// Add the value.
		add_val(ext_array, cur_g_val);
	} // t loop.

	// Add middle entry.
	double cur_sigma_squared = ((sigma * scale) * (sigma * scale));
	double cur_g_val = (1/denom) * delta_t * exp(-1 * 0 * 0 / (2 * cur_sigma_squared));

	// Add the value.
	add_val(ext_array, cur_g_val);

	// This adds the positive side: Go forward over the mesh.
	double total_val = 0;
	for(int i_mesh = 0; i_mesh < (int)gauss_mesh->size(); i_mesh++)
	{
		double t = gauss_mesh->at(i_mesh);
		double cur_sigma_squared = ((sigma * scale) * (sigma * scale));
		double cur_g_val = (1/denom) * delta_t * exp(-1 * t * t / (2 * cur_sigma_squared));

		// Add the value.
		add_val(ext_array, cur_g_val);

		total_val += cur_g_val;
	} // t loop.

	// Normalize the gaussian.
	for(int i_mesh = 0; i_mesh < (int)gauss_mesh->size(); i_mesh++)
	{
		ext_array->buffer[i_mesh] = ext_array->buffer[i_mesh] / total_val;
	}

	// The length is the number of values.
	l_filter = ext_array->n_vals;

	if(l_filter % 2 != 1)
	{
		fprintf(stderr, "The filter is not of odd length.\n");
		exit(0);
	}

	// Delete the mesh.
	delete(gauss_mesh);

	// This is a little tricky: Free the extended array buffer without freeing the buffer memory.
	double* gauss_buffer = ext_array->buffer;
	delete(ext_array);

	return(gauss_buffer);
}

#ifdef __unix__
gsl_vector* get_gsl_vector_per_array(double* data, int l_data)
{
	// Allocate the block.
	gsl_block* vec_block = new gsl_block();
	vec_block->size = l_data;
	vec_block->data = data;

	// Allocate the gsl_vector, then copy the data.
	gsl_vector* vec = new gsl_vector();
	vec->size = l_data;
	vec->data = data;
	vec->block = vec_block;
	vec->owner = 0;
	vec->stride = 1;

	return(vec);
}
#endif

void get_next_2_exp(int val, int& larger_exp_val, int& expon)
{
	int cur_val = 1;
	int cur_exp = 0;
	while(cur_val < val)
	{
		cur_exp++;
		cur_val *= 2;
	}

	expon = cur_exp;
	larger_exp_val = cur_val;
}

vector<double*>* multiscale_conv_filter_data(double* cur_real_track_data, int i_t, int l_track_data, double scale_start, double scale_end, double log_scale_step, vector<double>* scales_per_i_scale)
{	
	//for(int i = 0; i < 1000; i++)
	//{
	//	fprintf(stderr, "%d: %lf\n", i, cur_real_track_data[i]);
	//} // i loop.

	//getc(stdin);

	vector<double*>* decomps = new vector<double*>();

	// Take the FFT of track outside the loop. Copy the track everytime.
	int l_pre_ext = l_track_data/2;
	int l_ext_buf = l_track_data+l_pre_ext+l_pre_ext;
	int l_ext = (l_ext_buf - l_track_data) / 2;

	if(l_ext > l_track_data)
	{
		fprintf(stderr, "Must increase the track data length, it is shorter than the extension length.\n");
		exit(0);
	}

	fprintf(stderr, "In multiscale_conv_filter_data: %d, %d, %d\n", l_track_data, l_ext_buf, l_ext);

	double* buffered_track_data = new double[l_ext_buf+2];
	if(buffered_track_data == NULL)
	{
		fprintf(stderr, "Could not allocate data.\n");
		exit(0);
	}

	//memset(buffered_track_data, 0, sizeof(double) * l_ext_buf);

	// Do mirror image extension at the ends with the preset extension length.
    for(int i = 0; i < l_ext; i++)
    {
            buffered_track_data[i] = cur_real_track_data[l_ext-1-i];
    }
    for(int i = l_ext; i < l_ext + l_track_data; i++)
    {
            buffered_track_data[i] = cur_real_track_data[i-l_ext];
    }
    for(int i = l_ext+l_track_data; i < l_ext+l_track_data+l_ext; i++)
    {
            buffered_track_data[i] = cur_real_track_data[l_track_data-(i-l_ext-l_track_data)-1];
    }

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		i_scale++;
		scale *= log_scale_step;

		//if(i_scale >= 15)
		//{
		//	break;
		//}

		if(scale_end > 0 && 
			scale >= scale_end)
		{
			break;
		}

		scales_per_i_scale->push_back(scale);

		fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);
		// At this point, we have two track data: One is the loaded real data, other is the complex based data.

		// Get the filter for the current scale: The gaussian filter.
		int l_filter = 0;
		double* cur_filter_array = get_scaled_odd_length_gaussian(scale, 1.0, 8.0, l_filter);
		//fprintf(stderr, "%d filter weights.\n", l_filter);
		//getc(stdin);
		//double* cur_filter_array = get_extended_odd_length_gaussian(1.0, scale, 8.0, l_filter);

if(__DUMP_FILTER_MSGS__)
{
		char cur_filt_op_fp[1000];
		sprintf(cur_filt_op_fp, "filter_%lf.txt", scale);
		FILE* f_filt = open_f(cur_filt_op_fp, "w");
		for(int i = 0; i < l_filter; i++)
		{
			fprintf(f_filt, "%lf\n", cur_filter_array[i]);
		} // i loop.
		fclose(f_filt);
}
		if(l_filter > l_ext)
		{
			delete [] cur_filter_array;
			break;
		}

		double* cur_filter_params = new double[l_filter];
		memset(cur_filter_params, 0, sizeof(double) * l_filter);
		for(int i = 0; i < l_filter; i++)
		{
			cur_filter_params[i] = cur_filter_array[i];
		} // i loop.

		double* filtered_track_data = new double[l_ext_buf];
		memset(filtered_track_data, 0, l_ext_buf * sizeof(double));
		for(int i_d = 0; i_d < l_ext_buf; i_d++)
		{
			double cur_filtered_val = 0;

			// Compute the filtered value at i^{th} val.
			for(int i_c = 0; i_c < l_filter; i_c++)
			{
				if(i_d > i_c)
				{
					cur_filtered_val += (cur_filter_array[i_c] * buffered_track_data[i_d - i_c]);
				}
			} // i_c loop.

			filtered_track_data[i_d] = cur_filtered_val;
		} // i loop.


if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Done!\n");	

		// Copy the valid region.
		//t_ext_double_array* filtered_data = alloc_ext_array();
		int i_start = l_filter/2 + l_ext;
		//int i_start = l_ext;
		int i_end = i_start + l_track_data;

		// Do copying from start to end.
        double* cur_filtered_track = new double[i_end - i_start + 1];

		//for(int i = i_start-30; i <= i_start+30; i++)
		//{
		//	fprintf(stderr, "%d: %lf\n", i, filtered_track_data[i]);
		//}
		//getc(stdin);

		for(int i = i_start; i < i_end; i++)
		{
			//add_val(filtered_data, filtered_signal_fft[2*i]);
            cur_filtered_track[i-i_start] = filtered_track_data[i];
		} // i loop.

if(__DUMP_FILTER_MSGS__)
{
		char cur_decomp_fp[1000];
		sprintf(cur_decomp_fp, "decomp_%ld.txt", decomps->size());

		FILE* f_cur_decomp = open_f(cur_decomp_fp, "w");

		// Dump the decomposition and the extrema.
		for(int i = 0; i < l_track_data; i++)
		{
			//if(i < 30)
			//{
			//	fprintf(stderr, "%lf: %d: %lf\n", scale, i, cur_filtered_track[i]);
			//}
			fprintf(f_cur_decomp, "%.15f ", cur_filtered_track[i]);
		} // i loop.

		fclose(f_cur_decomp);
		//fprintf(f_decomp, "\n");
} // __DUMP_FILTER_MSGS__

		// Copy the data.
		decomps->push_back(cur_filtered_track);		

		delete [] cur_filter_params;
		delete [] cur_filter_array;
	} // scale loop.

	delete [] buffered_track_data;

	return(decomps);
}

vector<double*>* multiscale_gaussian_filter_data(double* cur_real_track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool dump_extrema_regions,
	char* op_file_prefix)
{
#ifdef WIN32
	return(NULL);
#elif defined __unix__
	vector<double*>* decomps = new vector<double*>();

	// Take the FFT of track outside the loop. Copy the track everytime.
	//int l_pre_ext = 1000*1000;
	int l_pre_ext = l_track_data/20;
	int l_ext_buf = l_track_data+l_pre_ext+l_pre_ext;
	int larger_exp_val = 0;
	int expon = 0;
	get_next_2_exp(l_ext_buf, larger_exp_val, expon);
	l_ext_buf = larger_exp_val;
	int l_ext = (l_ext_buf - l_track_data) / 2;

	if(l_ext > l_track_data)
	{
		fprintf(stderr, "Must increase the track data length, it is shorter than the extension length.\n");
		exit(0);
	}

	// Set up the workspace and wavetable.
	gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(l_ext_buf);
	gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(l_ext_buf);

	double* buffered_track_data = new double[2 * l_ext_buf];
	memset(buffered_track_data, 0, sizeof(double) * 2 * l_ext_buf);

	// Do mirror image extension at the ends with the preset extension length.
    for(int i = 0; i < l_ext; i++)
    {
            buffered_track_data[2*i] = cur_real_track_data[l_ext-i-1];
    }
    for(int i = l_ext; i < l_ext + l_track_data; i++)
    {
            buffered_track_data[2*i] = cur_real_track_data[i-l_ext];
    }
    for(int i = l_ext+l_track_data; i < l_ext+l_track_data+l_ext; i++)
    {
            buffered_track_data[2*i] = cur_real_track_data[l_track_data-(i-l_ext-l_track_data)-1];
    }

if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "FFT'ing buffered track data.\n");
	gsl_fft_complex_forward (buffered_track_data, 1, l_ext_buf,
						wavetable, workspace);
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Finished FFT'ing buffered track data.\n");

	// Start from the first scale, go over all the scales and filter the data.
	//double scale = 1.0 / log_scale_step;
	double scale = scale_start / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		if(scale_end > 0 && 
			scale >= scale_end)
		{
			break;
		}

		scales_per_i_scale->push_back(scale);

		fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);
		// At this point, we have two track data: One is the loaded real data, other is the complex based data.

		// Get the filter for the current scale: The gaussian filter.
		int l_filter = 0;
		double* cur_filter_array = get_scaled_odd_length_gaussian(scale, 1.0, 8.0, l_filter);
		//double* cur_filter_array = get_extended_odd_length_gaussian(1.0, scale, 8.0, l_filter);

if(__DUMP_FILTER_MSGS__)
{
		char cur_filt_op_fp[1000];
		sprintf(cur_filt_op_fp, "filter_%lf.txt", scale);
		FILE* f_filt = open_f(cur_filt_op_fp, "w");
		for(int i = 0; i < l_filter; i++)
		{
			fprintf(f_filt, "%lf\n", cur_filter_array[i]);
		} // i loop.
		fclose(f_filt);
}

		if(l_filter > l_ext)
		{
			delete [] cur_filter_array;
			break;
		}

		double* cur_filter_params = new double[2 * l_ext_buf];
		memset(cur_filter_params, 0, sizeof(double) * 2 * l_ext_buf);
		for(int i = 0; i < l_filter; i++)
		{
			cur_filter_params[2 * i] = cur_filter_array[i];
		} // i loop.

if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Starting FFT'ing filter parameters at scale %lf.\n", scale);
		gsl_fft_complex_forward(cur_filter_params, 1, l_ext_buf,
							wavetable, workspace);

if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Finished FFT'ing filter parameters.\n");

		// Multiply the filter and the signal track fft's.
		fprintf(stderr, "Retrieving real and imaginary parts of the track and filter data.\n");
		double* filtered_signal_fft = new double[2 * l_ext_buf];
		for(int i = 0; i < l_ext_buf; i+=1)
		{
			double cur_real_val1 = buffered_track_data[2*i] * cur_filter_params[2*i] - (buffered_track_data[2*i+1] * cur_filter_params[2*i+1]);
			double cur_imag_val1 = buffered_track_data[2*i+1] * cur_filter_params[2*i] + (buffered_track_data[2*i] * cur_filter_params[2*i+1]);
			filtered_signal_fft[2*i] = cur_real_val1;
			filtered_signal_fft[2*i+1] = cur_imag_val1;
		} // i loop.

		// Take the inverse fft, and we are done.
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Taking inverse fft of multiplication.\n");

		gsl_fft_complex_inverse(filtered_signal_fft, 1, l_ext_buf, wavetable, workspace);

if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Done!\n");	

		// Copy the valid region.
		//t_ext_double_array* filtered_data = alloc_ext_array();
		int i_start = l_filter/2 + l_ext;
		int i_end = i_start + l_track_data;

if(__DUMP_FILTER_MSGS__)
{
		double* data_copy = new double[l_ext_buf];
		for(int i = 0; i < l_ext_buf; i++)
		{
			data_copy[i] = filtered_signal_fft[2*i];
		}

		fprintf(stderr, "Start:\n");
		for(int i = i_start-30; i <= i_start+30; i++)
		{
			fprintf(stderr, "%d: %.5f\n", i-i_start, data_copy[i]);
		}

		fprintf(stderr, "%d-%d\n", i_start, i_end);
		FILE* f_data_copy = open_f("data_copy.txt", "w");
		for(int i = 0; i < l_ext_buf; i++)
		{
			fprintf(f_data_copy, "%.5f\n",  data_copy[i]);
		}
		fclose(f_data_copy);

		fprintf(stderr, "Dumped the whole filtered data to data_copy.txt");
		getc(stdin);

		delete [] data_copy;
} // whole data dump check.

		// Do copying from start to end.
        double* cur_filtered_track = new double[i_end - i_start + 1];

		for(int i = i_start; i < i_end; i++)
		{
			//add_val(filtered_data, filtered_signal_fft[2*i]);
            cur_filtered_track[i-i_start] = filtered_signal_fft[2*i];

			//// Make sure that the imaginary value is 0.
			//if(fabs(filtered_signal_fft[2*i+1] > 0.000001))
			//{
			//	fprintf(stderr, "IFFT problem.\n");
			//	exit(0);
			//}
		} // i loop.

if(__DUMP_FILTER_MSGS__)
{
		char cur_decomp_fp[1000];
		sprintf(cur_decomp_fp, "decomp_%ld.txt", decomps->size());

		FILE* f_cur_decomp = open_f(cur_decomp_fp, "w");

		// Dump the decomposition and the extrema.
		for(int i = 0; i < l_track_data; i++)
		{
			fprintf(f_cur_decomp, "%.5f ", cur_filtered_track[i]);
		} // i loop.

		fclose(f_cur_decomp);
		//fprintf(f_decomp, "\n");
} // __DUMP_FILTER_MSGS__

		// Copy the data.
		if(dump_decomposition)
		{
			// Dump and free memory.
			char cur_decomp_fp[1000];
			//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
			sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
			fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
			double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
			dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
			fprintf(stderr, "Closing file\n");

			// Clean track memory.
			delete [] _cur_filtered_track;
			delete [] cur_filtered_track;
		}
		else if(dump_extrema_regions)
		{
			// Get the minima and maxima, then dump.
			fprintf(stderr, "Getting the extrema regions.\n");

			// Get the minima and maxima, then dump.
			char cur_minima_regs_fp[1000];
			char cur_maxima_regs_fp[1000];
			sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
			sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

			// Get the extrema regions for the current filtered regions.
			vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
			vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
			int* derivative_map = new int[l_track_data+2];
			memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
			get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

			// Sort the extremas.
			sort(minima->begin(), minima->end(), sort_extremas_per_posn);
			sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

			fprintf(stderr, "Dumping the extrema regions.\n");
			FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
			for(int i_m = 0; i_m < (int)maxima->size(); i_m++)
			{
				fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
			} // i_m loop.
			fclose(f_maxes);

			FILE* f_mins = open_f(cur_minima_regs_fp, "w");
			for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
			{
				fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
			} // i_m loop.
			fclose(f_mins);

			// Clean extrema nodes.
			delete_extrema_nodes(minima);
			delete_extrema_nodes(maxima);

			// Free the current decomposition.
			delete [] derivative_map;
			delete [] cur_filtered_track;
		}
		else
		{
			decomps->push_back(cur_filtered_track);
		}

		i_scale++;

		delete [] cur_filter_params;
		delete [] cur_filter_array;
		delete [] filtered_signal_fft;
	} // scale loop.

	// Free workspace memory.
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);

	delete [] buffered_track_data;

	return(decomps);
#endif // __unix__ check.	
}

int sort_doubles_descending(const void* p1, const void* p2)
{
	double val1 = *((double*)p1);
	double val2 = *((double*)p2);

	if(val1 < val2)
	{
		return(-1);
	}
	else if(val1 == val2)
	{
		return(0);
	}
	else
	{
		return(1);
	}
}

double* median_filter_data(double* track_data,
	int l_track_data,
	int l_averaging_win,
	int skip_value) // This is the value to be skipped from the window values.
{
	// Do copying from start to end.
	double* cur_filtered_track = new double[l_track_data + 2];

	int n_signal_wins = 0;
	int half_l_averaging = l_averaging_win / 2;

	// Set the maximum value for the histogram.
	int MAX_VAL = -1000 * 1000;
	int MIN_VAL = 1000 * 1000;
	for (int i_sig = 1; i_sig <= l_track_data; i_sig++)
	{
		if (MAX_VAL < (int)(track_data[i_sig]))
		{
			MAX_VAL = (int)(track_data[i_sig]);
		}

		if (MIN_VAL >(int)(track_data[i_sig]))
		{
			MIN_VAL = (int)(track_data[i_sig]);
		}
	} // i_sig loop.

	  // Make sure that this is the maximum.
	MAX_VAL += 1000;
	MIN_VAL -= 1000;

	// Initialize the current pdf.
	int* cur_win_pdf = NULL;
	int cur_win_max = 0;
	int cur_win_min = 1000 * 1000;
	double* cur_win_vals = new double[l_averaging_win + 2];
	int prev_avg_start = 0;
	int prev_avg_end = 0;

	// Go over all the positions as the middle of the filtering window.
	for (int cur_win_mid = 1; cur_win_mid <= l_track_data; cur_win_mid++)
	{
		int cur_avg_start = (cur_win_mid > half_l_averaging) ? (cur_win_mid - half_l_averaging) : (1);
		int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data) ? (cur_win_mid + half_l_averaging) : (l_track_data);
		if (cur_win_pdf == NULL)
		{
			cur_win_pdf = new int[MAX_VAL - MIN_VAL + 2];
			memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
			//t_string::set_byte_buffer(cur_win_pdf, sizeof(int) * (MAX_VAL - MIN_VAL + 1), 0);
			cur_win_pdf -= MIN_VAL;

			// Generate the pdf, get the minimum and maximum in the current window.
			for (int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if (cur_win_max < track_data[i])
					{
						cur_win_max = track_data[i];
					}

					if (cur_win_min > track_data[i])
					{
						cur_win_min = track_data[i];
					}
				} // skip_Value check.
			} // i loop.
		} // cur_win_pdf is NULL check.
		else
		{
			// Remove the old values from the pdf, add the new values.
			for (int i = prev_avg_start; i < cur_avg_start; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]--;
				}
			} // i loop.

			  // Update the min and max only for the values that are new in this window.
			for (int i = prev_avg_end + 1; i <= cur_avg_end; i++)
			{
				if (track_data[i] != skip_value)
				{
					cur_win_pdf[(int)(track_data[i])]++;

					if (cur_win_max < track_data[i])
					{
						cur_win_max = (int)track_data[i];
					}

					if (cur_win_min > track_data[i])
					{
						cur_win_min = (int)track_data[i];
					}
				} // skip_value check.
			}

			// Sanity Check: The total # of points must be equal to the window length.
			int n_total_pts = 0;
			for (int i = cur_win_min; i <= cur_win_max; i++)
			{
				if (i != skip_value)
				{
					n_total_pts += cur_win_pdf[i];
				}
			} // i loop.
		} // cur_win_pdf is NULL check.

		  // Count the total number of points without the skip value.
		int n_total = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				n_total += cur_win_pdf[i];
			}
		} // i loop.

		  // Generate the window median.
		int cur_win_median = 0;
		int n_cur_total = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				n_cur_total += cur_win_pdf[i];
				if (n_cur_total > n_total / 2)
				{
					// We found the median, can break out of the loop.
					cur_win_median = i;
					break;
				}
			}
		} // i loop.

		  // Track the minimum and maximum from the histogram to update them for the current window.
		int updated_win_min = 0;
		for (int i = cur_win_min; i <= cur_win_max; i++)
		{
			if (i != skip_value)
			{
				if (cur_win_pdf[i] > 0)
				{
					updated_win_min = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		int updated_win_max = 0;
		for (int i = cur_win_max; i >= cur_win_min; i--)
		{
			if (i != skip_value)
			{
				if (cur_win_pdf[i] > 0)
				{
					updated_win_max = i;
					break;
				}
			} // skip_value check.
		} // i loop.

		  // Set the median.
		  //int median_per_pdf = cur_win_median;
		cur_filtered_track[cur_win_mid] = cur_win_median;

		// Update the previous averaging window limits.
		prev_avg_start = cur_avg_start;
		prev_avg_end = cur_avg_end;
		cur_win_min = updated_win_min;
		cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
		// Get the median via qsort and compare as a sanity check.
		int n_valid_vals = 0;
		for (int i = cur_avg_start; i <= cur_avg_end; i++)
		{
			if (track_data[i] != skip_value)
			{
				cur_win_vals[i - cur_avg_start] = track_data[i];
				n_valid_vals++;
			}
		} // i loop.
		qsort(cur_win_vals, n_valid_vals, sizeof(double), sort_doubles_descending);
		//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

		if (cur_win_vals[n_valid_vals / 2] != median_per_pdf)
		{
			fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[n_valid_vals / 2], median_per_pdf);
			for (int i = cur_avg_start; i <= cur_avg_end; i++)
			{
				fprintf(stderr, "%lf ", track_data[i]);
			} // i loop.
			fprintf(stderr, "\n");
			getc(stdin);
		}
#endif // _QSORT_CHECK_

		n_signal_wins++;
	} // main signal filtering loop.

	  // Free pdf memory.
	delete[] cur_win_vals;
	delete[](cur_win_pdf + MIN_VAL);

	return(cur_filtered_track);
}

double* mean_filter_data(double* signal, int l_sig, int l_bin)
{
	int prev_win_start = 0;
	int prev_win_end = 0;
	double cur_total_win_signal = -1;
	if(l_bin % 2 == 0)
	{
		fprintf(stderr, "Resetting bin length to %d\n", l_bin);
		l_bin++;
	}

	double* filtered_dat = new double[l_sig + 2];
	memset(filtered_dat, 0, sizeof(double) * (l_sig + 1));

	int half_l_win = (l_bin - 1) / 2;

	for(int i_mid = 1; i_mid <= l_sig; i_mid++)
	{
		int cur_win_start = MAX(1, i_mid - half_l_win);
		int cur_win_end = MIN(l_sig, i_mid + half_l_win);

		if(cur_total_win_signal == -1)
		{
			for(int i = cur_win_start; i <= cur_win_end; i++)
			{
				cur_total_win_signal += signal[i];
			} // i loop.
		}
		else
		{
			// Remove the exiting signals.
			for(int i = prev_win_start; i < cur_win_start; i++)
			{
				cur_total_win_signal -= signal[i];
			} // i loop.

			// Add the newly entered signals.
			for(int i = prev_win_end+1; i <= cur_win_end; i++)
			{
				cur_total_win_signal += signal[i];
			} // i loop.
		}

		// Update the previous window coordinates.
		prev_win_start = cur_win_start;
		prev_win_end = cur_win_end;

		// Set the filtered data.
		filtered_dat[i_mid] = cur_total_win_signal / l_bin;
	} // i_mid loop.

	return(filtered_dat);
}

vector<double*>* multiscale_avg_filter_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	bool dump_decomposition,
	bool compute_extrema_regions,
	bool dump_extrema_regions,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions, 
	vector<vector<t_annot_region*>*>* per_scale_maxima_regions,
	char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
{
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);
}

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Initialize the current pdf.
			//int* cur_win_pdf = NULL;
			//int cur_win_max = 0;
			//int cur_win_min = 1000*1000;
			//double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			double cur_win_total = -123123;
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_total == -123123)
				{
					// Get the total by a loop for the current window.
					cur_win_total = 0.0;
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_total += track_data[i];
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_total -= track_data[i];
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_total += track_data[i];
					}
				} // cur_win_pdf is NULL check.

				// Set the median.
				cur_filtered_track[cur_win_mid] = cur_win_total / l_averaging_win;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;

#undef _MANUAL_TOTAL_CHECK_
#ifdef _MANUAL_TOTAL_CHECK_
				// Get the median via qsort and compare as a sanity check.
				double manual_total = 0.0;
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					manual_total += track_data[i];
				} // i loop.

				if(manual_total != cur_win_total)
				{
					fprintf(stderr, "Sanity check failed at manual total check, %s(%d)\n", __FILE__, __LINE__);
					exit(0);
				}
#endif // _MANUAL_TOTAL_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				fprintf(stderr, "Closing file\n");

				// Clean track memory.
				delete [] _cur_filtered_track;
				delete [] cur_filtered_track;
			}
			else if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size(); i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size(); i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
				delete [] cur_filtered_track;
			}
			else
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
} // multiscale_avg_filter_data


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
	char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			double prev_non_skip_median_value = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						if(track_data[i] != skip_value)
						{
							cur_win_pdf[(int)(track_data[i])]++;

							if(cur_win_max < track_data[i])
							{
								cur_win_max = track_data[i];
							}
					
							if(cur_win_min > track_data[i])
							{
								cur_win_min = track_data[i];
							}
						} // skip_Value check.
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						if(track_data[i] != skip_value)
						{
							cur_win_pdf[(int)(track_data[i])]--;
						}
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						if(track_data[i] != skip_value)
						{
							cur_win_pdf[(int)(track_data[i])]++;

							if(cur_win_max < track_data[i])
							{
								cur_win_max = (int)track_data[i];
							}
					
							if(cur_win_min > track_data[i])
							{
								cur_win_min = (int)track_data[i];
							}
						} // skip_value check.
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						if(i != skip_value)
						{
							n_total_pts += cur_win_pdf[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.

				// Count the total number of points without the skip value.
				int n_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(i != skip_value)
					{
						n_total += cur_win_pdf[i];
					}
				} // i loop.
		
				// Generate the window median.
				int cur_win_median = 0;
				int n_cur_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(i != skip_value)
					{
						n_cur_total += cur_win_pdf[i];
						if(n_cur_total > n_total/2)
						{
							// We found the median, can break out of the loop.
							cur_win_median = i;
							break;
						}
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(i != skip_value)
					{
						if(cur_win_pdf[i] > 0)
						{
							updated_win_min = i;
							break;
						}
					} // skip_value check.
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(i != skip_value)
					{
						if(cur_win_pdf[i] > 0)
						{
							updated_win_max = i;
							break;
						}
					} // skip_value check.
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				if(n_cur_total > 0)
				{
					cur_filtered_track[cur_win_mid] = cur_win_median;
				}
				else
				{
					cur_filtered_track[cur_win_mid] = prev_non_skip_median_value;
				}

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

				if(n_cur_total > 0)
				{
					prev_non_skip_median_value = cur_win_median;
				}

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				int n_valid_vals = 0;
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					if(track_data[i] != skip_value)
					{
						cur_win_vals[i-cur_avg_start] = track_data[i];
						n_valid_vals++;
					}
				} // i loop.
				qsort(cur_win_vals, n_valid_vals, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[n_valid_vals/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[n_valid_vals/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

			// Sanity check on the filtered track: Make sure that there is no skip value in the filtered track data.
			for(int i = 1; i <= l_track_data; i++)
			{
				if(cur_filtered_track[i] == skip_value)
				{
					fprintf(stderr, "Imputed MS filtering failed.\n");
					exit(0);
				}
			} // i loop.

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				delete [] _cur_filtered_track;
				fprintf(stderr, "Closing file\n");
			}
			
			if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
			}
			
			// Do we need the filtered track back?
			if(return_filtered_tracks)
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}
			else
			{
				// Clean track memory.
				delete [] cur_filtered_track;
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
}

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
	char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];
				} // i loop.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int cur_n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					cur_n_total += cur_win_pdf[i];
					if(cur_n_total > n_total/2)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				delete [] _cur_filtered_track;
				fprintf(stderr, "Closing file\n");
			}
			
			if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
			}
			
			// Do we need the filtered track back?
			if(return_filtered_tracks)
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}
			else
			{
				// Clean track memory.
				delete [] cur_filtered_track;
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
}

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
	double quantile_fraction)
{
	double* track_data = new double[l_track_data + 3];
	for(int i = 1; i <= l_track_data; i++)
	{
		track_data[i] = _track_data[i];
	} // i loop.

	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				int n_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];
				} // i loop.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int cur_n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					cur_n_total += cur_win_pdf[i];
					if(cur_n_total > n_total * quantile_fraction)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				delete [] _cur_filtered_track;
				fprintf(stderr, "Closing file\n");
			}
			
			if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
			}
		
			// Copy the current filtered track data to the track data.
			for(int i = 1; i <= l_track_data; i++)
			{
				track_data[i] = cur_filtered_track[i];
			}

			// Do we need the filtered track back?
			if(return_filtered_tracks)
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}
			else
			{
				// Clean track memory.
				delete [] cur_filtered_track;
			}

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	delete [] track_data;

	return(decomps);
}



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
	char* op_file_prefix)
{
	vector<double*>* decomps = new vector<double*>();

	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;
	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			decomps->push_back(NULL);

			if(compute_extrema_regions)
			{
				if(per_scale_minima_regions != NULL)
				{
					per_scale_minima_regions->push_back(new vector<t_annot_region*>());
				}

				if(per_scale_maxima_regions != NULL)
				{
					per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
				}
			}

			i_scale++;
			continue;
		}
		else
		{

if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Permute the current track data.
			int l_track_data;
			double* track_data = get_block_permute_profile(original_track_data, rng, l_original_track_data, l_track_data, int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < (int)track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > (int)track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];
					if(n_total > l_averaging_win/2)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

			// Copy the data.
			if(dump_decomposition)
			{
				// Dump and free memory.
				char cur_decomp_fp[1000];
				//sprintf(cur_decomp_fp, "decomp_%d_%d_%d.bin", i_scale, cur_win_start, cur_win_start + cur_l_win);
				sprintf(cur_decomp_fp, "%s_scale_%d.bin", op_file_prefix, i_scale);
				fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
				double* _cur_filtered_track = get_one_indexed_per_zero_indexed_data(cur_filtered_track, l_track_data);
				dump_per_nucleotide_binary_profile(_cur_filtered_track, l_track_data, cur_decomp_fp);
				fprintf(stderr, "Closing file\n");

				// Clean track memory.
				delete [] _cur_filtered_track;
			}
			
			if(compute_extrema_regions)
			{
if(__DUMP_FILTER_MSGS__)
				fprintf(stderr, "Getting the extrema regions.\n");

				// Get the extrema regions for the current filtered regions.
				vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
				vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
				int* derivative_map = new int[l_track_data+2];
				memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
				get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

				// Sort the extremas.
				sort(minima->begin(), minima->end(), sort_extremas_per_posn);
				sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

				if(dump_extrema_regions)
				{
					// Get the minima and maxima, then dump.
					char cur_minima_regs_fp[1000];
					char cur_maxima_regs_fp[1000];
					sprintf(cur_minima_regs_fp, "%s_scale_%d_mins.bed", op_file_prefix, i_scale);
					sprintf(cur_maxima_regs_fp, "%s_scale_%d_maxes.bed", op_file_prefix, i_scale);

					fprintf(stderr, "Dumping the extrema regions.\n");
					FILE* f_maxes = open_f(cur_maxima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						fprintf(f_maxes, "XX\t%d\t%d\n", translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_maxes);

					FILE* f_mins = open_f(cur_minima_regs_fp, "w");
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						fprintf(f_mins, "XX\t%d\t%d\n", translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
							translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
					} // i_m loop.
					fclose(f_mins);
				}

				if(per_scale_minima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
					{
						t_annot_region* new_minima = get_empty_region();
						new_minima->chrom = t_string::copy_me_str("XX");
						new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_minima->strand = '+';

						cur_scale_minima_regions->push_back(new_minima);
					} // i_min loop.

					per_scale_minima_regions->push_back(cur_scale_minima_regions);
				}
			
				if(per_scale_maxima_regions != NULL)
				{
					vector<t_annot_region*>* cur_scale_maxima_regions = new vector<t_annot_region*>();
					for(int i_m = 0; i_m < (int)maxima->size()-1; i_m++)
					{
						t_annot_region* new_maxima = get_empty_region();
						new_maxima->chrom = t_string::copy_me_str("XX");
						new_maxima->start = translate_coord(maxima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
						new_maxima->end = translate_coord(maxima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
						new_maxima->strand = '+';

						cur_scale_maxima_regions->push_back(new_maxima);
					} // i_m loop.

					per_scale_maxima_regions->push_back(cur_scale_maxima_regions);
				}

				delete_extrema_nodes(minima);
				delete_extrema_nodes(maxima);
			
				// Free the current decomposition.
				delete [] derivative_map;
			}
						
			// Do we need the filtered track back?
			if(return_filtered_tracks)
			{
				// Store the decompositions.
				decomps->push_back(cur_filtered_track);
			}
			else
			{
				// Clean track memory.
				delete [] cur_filtered_track;
			}			

			// Delete the block shuffled track data.
			delete [] track_data;

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

	return(decomps);
}


double* mapability_aware_median_filter(double* signal_profile, int l_profile,
	double* scaled_mapability_profile, int l_mapability_profile,
	double max_mapable_signal_2_use_in_filter,
	int l_mapability_filtering_win)
{
	// Generate the signal profile to input to the filtering.
	double* input_signal_profile = new double[l_profile + 2];
	for(int i = 1; i <= l_profile; i++)
	{
		if(l_mapability_profile > i)
		{
			if(scaled_mapability_profile[i] > max_mapable_signal_2_use_in_filter)
			{
				input_signal_profile[i] = -1;
			}
			else
			{
				input_signal_profile[i] = signal_profile[i];
			}
		}
		else
		{
			input_signal_profile[i] = signal_profile[i];
		}
	} // i loop.

	fprintf(stderr, "Mapability aware filtering.\n");
	double* mapable_median_profile = median_filter_data(input_signal_profile,
													l_profile, 
													l_mapability_filtering_win,
													-1);

	double* filtered_signal = new double[l_profile + 5];

	// Compare the actual signal profile with the 
	for(int i = 1; i <= l_profile; i++)
	{
		filtered_signal[i] = MAX(mapable_median_profile[i], signal_profile[i]);
	} // i loop.

	delete [] mapable_median_profile;
	delete [] input_signal_profile;

	return(filtered_signal);
}
	







void get_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale)
{
	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;

	FILE* f_feature_signal_vals = NULL;
if(__DUMP_FILTER_MSGS__)
	f_feature_signal_vals = open_f("features_signal_values.bed", "w");

	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);


		//if(scale < scale_start)
		//{
		//	decomps->push_back(NULL);

		//	if(compute_extrema_regions)
		//	{
		//		if(per_scale_minima_regions != NULL)
		//		{
		//			
		//		}

		//		if(per_scale_maxima_regions != NULL)
		//		{
		//			per_scale_maxima_regions->push_back(new vector<t_annot_region*>());
		//		}
		//	}

		//	i_scale++;
		//	continue;
		//}


		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			i_scale++;
			per_scale_minima_regions->push_back(new vector<t_annot_region*>());
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						cur_win_pdf[(int)(track_data[i])]--;
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						cur_win_pdf[(int)(track_data[i])]++;

						if(cur_win_max < track_data[i])
						{
							cur_win_max = (int)track_data[i];
						}
					
						if(cur_win_min > track_data[i])
						{
							cur_win_min = (int)track_data[i];
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						n_total_pts += cur_win_pdf[i];
					} // i loop.

					// Sanity check on the number of points to be processed in the current window.
					if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					{
						fprintf(stderr, "Sanity check failed at %d-%d (%d-%d): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

						for(int i = cur_win_min; i <= cur_win_max; i++)
						{
							fprintf(stderr, "%d ", cur_win_pdf[i]);
						} // i loop.

						FILE* f_op = open_f("op.txt", "w");
						for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
						{
							fprintf(f_op, "%lf ", track_data[i_sig]);
						} // i_sig loop.
						fclose(f_op);

						exit(0);
					}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				int cur_win_median = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					n_total += cur_win_pdf[i];

					//if(n_total > l_averaging_win/2)
					if(n_total > ((cur_avg_end - cur_avg_start + 1))/2)
					{
						// We found the median, can break out of the loop.
						cur_win_median = i;
						break;
					}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_min = i;
						break;
					}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					if(cur_win_pdf[i] > 0)
					{
						updated_win_max = i;
						break;
					}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);

if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Getting the extrema regions.\n");

			// Get the extrema regions for the current filtered regions.
			vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
			vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
			int* derivative_map = new int[l_track_data+2];
			memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
			get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

			// Sort the extremas.
			sort(minima->begin(), minima->end(), sort_extremas_per_posn);
			sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

			// Get the minima regions.
			vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();

			fprintf(stderr, "Generating filtered minima regions from %d minima.\n", (int)minima->size());

			for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
			{
				t_annot_region* new_minima = get_empty_region();
				new_minima->chrom = t_string::copy_me_str("XX");
				new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
				new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
				new_minima->strand = '+';

				// The maximum in the original signal is found within the effective data that contributes to the filtering process.
				int i_eff_data_start = (1 > new_minima->start-(int)scale/2)?(1):(new_minima->start-(int)scale/2);
				int i_eff_data_end = (l_track_data < new_minima->end+scale/2)?(l_track_data):(new_minima->end+(int)scale/2);
				double original_max = 0.0;
				for(int i = i_eff_data_start; i <= i_eff_data_end; i++)
				{
					if(original_max < track_data[i])
					{
						original_max = track_data[i];
					}
				}

				double filtered_max = 0.0;
				for(int i = new_minima->start; i <= new_minima->end; i++)
				{
					if(filtered_max < cur_filtered_track[i])
					{
						filtered_max = cur_filtered_track[i];
					}
				} // i loop.

if(__DUMP_FILTER_MSGS__)
				fprintf(f_feature_signal_vals, "XX\t%d\t%d\t%lf\t%lf\t%lf\n", new_minima->start, new_minima->end, filtered_max, original_max, scale);

				// This is where the filtering is done.
				if(filtered_max >= original_max / min_allowed_filtered_value_scale)
				{
					cur_scale_minima_regions->push_back(new_minima);
				}
				else
				{					
					delete new_minima;
				}
			} // i_min loop.

			fprintf(stderr, "%ld minima regions are generated.\n", cur_scale_minima_regions->size());
			per_scale_minima_regions->push_back(cur_scale_minima_regions);

			delete_extrema_nodes(minima);
			delete_extrema_nodes(maxima);

if(__DUMP_FILTER_MSGS__)
{
			// Dump the current filtered track.
			char cur_filtered_data_op_fp[1000];
			sprintf(cur_filtered_data_op_fp, "filtered_%lf.bin", scale);
			dump_per_nucleotide_binary_profile(cur_filtered_track, l_track_data, cur_filtered_data_op_fp);
}

			// Free the current decomposition.
			delete [] derivative_map;
			delete [] cur_filtered_track;

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

if(__DUMP_FILTER_MSGS__)
		fclose(f_feature_signal_vals);
} // -get_filtered_maxima_regions_multiscale_filtered_data option.














void get_mapability_aware_median_filtered_maxima_regions_multiscale_filtered_data(double* track_data, 
	int l_track_data, 
	double scale_start, double scale_end, double log_scale_step, 
	vector<double>* scales_per_i_scale,
	vector<vector<t_annot_region*>*>* per_scale_minima_regions,
	double min_allowed_filtered_value_scale,
	double* normalized_mapability_signal,
	int l_mapability_profile,
	double max_mapability_signal)
{
	// Start from the first scale, go over all the scales and filter the data.
	double scale = 1.0 / log_scale_step;
	int i_scale = 0;

	FILE* f_feature_signal_vals = NULL;
if(__DUMP_FILTER_MSGS__)
	f_feature_signal_vals = open_f("features_signal_values.bed", "w");

	while(1)
	{
		scale *= log_scale_step;

		scales_per_i_scale->push_back(scale);

		// Are we going to process a scale within the requested limits? If not, we will skip this computation to save time.
		if(scale < scale_start)
		{
			i_scale++;
			per_scale_minima_regions->push_back(new vector<t_annot_region*>());
			continue;
		}
		else
		{
if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Processing scale %lf (%lf)\n", scale, scale_end);

			// Get the filter for the current scale: The gaussian filter.
			int int_scale = (int)(scale);
			int l_averaging_win = ((int_scale % 2) == 1)?(int_scale-1):(int_scale);

			// Do copying from start to end.
			double* cur_filtered_track = new double[l_track_data + 2];

			int n_signal_wins = 0;
			int half_l_averaging = l_averaging_win / 2;

			// Set the maximum value for the histogram.
			int MAX_VAL = -1000*1000;
			int MIN_VAL = 1000*1000;
			for(int i_sig = 1; i_sig <= l_track_data; i_sig++)
			{
				if(MAX_VAL < (int)(track_data[i_sig]))
				{
					MAX_VAL = (int)(track_data[i_sig]);
				}

				if(MIN_VAL > (int)(track_data[i_sig]))
				{
					MIN_VAL = (int)(track_data[i_sig]);
				}
			} // i_sig loop.

			// Make sure that this is the maximum.
			MAX_VAL += 1000;
			MIN_VAL -= 1000;
		
			// Initialize the current pdf.
			int* cur_win_pdf = NULL;
			int cur_win_max = 0;
			int cur_win_min = 1000*1000;
			double* cur_win_vals = new double[l_averaging_win + 2];
			int prev_avg_start = 0;
			int prev_avg_end = 0;

			// Go over all the positions as the middle of the filtering window.
			for(int cur_win_mid = 0; cur_win_mid < l_track_data; cur_win_mid++)
			{
				int cur_avg_start = (cur_win_mid > half_l_averaging)?(cur_win_mid - half_l_averaging):(1);
				int cur_avg_end = (cur_win_mid + half_l_averaging <= l_track_data)?(cur_win_mid + half_l_averaging):(l_track_data);
				if(cur_win_pdf == NULL)
				{
					cur_win_pdf = new int[MAX_VAL - MIN_VAL+ 2];
					memset(cur_win_pdf, 0, sizeof(int) * (MAX_VAL - MIN_VAL + 1));
					cur_win_pdf -= MIN_VAL;

					// Generate the pdf, get the minimum and maximum in the current window.
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						//if(track_data[i] != skip_value)
						if(i >= l_mapability_profile ||
							normalized_mapability_signal[i] < max_mapability_signal)
						{
							cur_win_pdf[(int)(track_data[i])]++;

							if(cur_win_max < track_data[i])
							{
								cur_win_max = (int)track_data[i];
							}
					
							if(cur_win_min > track_data[i])
							{
								cur_win_min = (int)track_data[i];
							}
						}
					} // i loop.
				} // cur_win_pdf is NULL check.
				else
				{
					// Remove the old values from the pdf, add the new values.
					for(int i = prev_avg_start; i < cur_avg_start; i++)
					{
						//if(track_data[i] != skip_value)
						if(i >= l_mapability_profile ||
							normalized_mapability_signal[i] < max_mapability_signal)
						{
							cur_win_pdf[(int)(track_data[i])]--;
						}
					} // i loop.

					// Update the min and max only for the values that are new in this window.
					for(int i = prev_avg_end+1; i <= cur_avg_end; i++)
					{
						//if(track_data[i] != skip_value)
						if(i >= l_mapability_profile ||
							normalized_mapability_signal[i] < max_mapability_signal)
						{
							cur_win_pdf[(int)(track_data[i])]++;

							if(cur_win_max < track_data[i])
							{
								cur_win_max = (int)track_data[i];
							}
					
							if(cur_win_min > track_data[i])
							{
								cur_win_min = (int)track_data[i];
							}
						}
					}

					// Sanity Check: The total # of points must be equal to the window length.
					int n_total_pts = 0;
					for(int i = cur_win_min; i <= cur_win_max; i++)
					{
						/*if(i != skip_value)
						{*/
							n_total_pts += cur_win_pdf[i];
						//}
					} // i loop.

					//// Sanity check on the number of points to be processed in the current window.
					//if(n_total_pts != (cur_avg_end - cur_avg_start + 1))
					//{
					//	fprintf(stderr, "Sanity check failed at %d-%d (%lf-%lf): %d, %d\n", cur_avg_start, cur_avg_end, cur_win_min, cur_win_max, n_total_pts, (cur_avg_end - cur_avg_start + 1));

					//	for(int i = cur_win_min; i <= cur_win_max; i++)
					//	{
					//		fprintf(stderr, "%d ", cur_win_pdf[i]);
					//	} // i loop.

					//	FILE* f_op = open_f("op.txt", "w");
					//	for(int i_sig = cur_avg_start; i_sig <= cur_avg_end; i_sig++)
					//	{
					//		fprintf(f_op, "%lf ", track_data[i_sig]);
					//	} // i_sig loop.
					//	fclose(f_op);

					//	exit(0);
					//}
				} // cur_win_pdf is NULL check.

				// At this point, the histogram is updated, now must find the median for the histogram.
				int n_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					/*if(i != skip_value)
					{*/
						n_total += cur_win_pdf[i];
					//}
				} // i loop.

				// Generate the window median.
				int cur_win_median = 0;
				int n_cur_total = 0;
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					/*if(i != skip_value)
					{*/
						n_cur_total += cur_win_pdf[i];
						if(n_cur_total > n_total/2)
						{
							// We found the median, can break out of the loop.
							cur_win_median = i;
							break;
						}
					//}
				} // i loop.

				// Track the minimum and maximum from the histogram to update them for the current window.
				int updated_win_min = 0;			
				for(int i = cur_win_min; i <= cur_win_max; i++)
				{
					/*if(i != skip_value)
					{*/
						if(cur_win_pdf[i] > 0)
						{
							updated_win_min = i;
							break;
						}
					//}
				} // i loop.

				int updated_win_max = 0;
				for(int i = cur_win_max; i >= cur_win_min; i--)
				{
					/*if(i != skip_value)
					{*/
						if(cur_win_pdf[i] > 0)
						{
							updated_win_max = i;
							break;
						}
					//}
				} // i loop.

				// Set the median.
				//int median_per_pdf = cur_win_median;
				//cur_filtered_track[cur_win_mid] = (cur_win_median > track_data[cur_win_mid])?(cur_win_median):(track_data[cur_win_mid]);
				cur_filtered_track[cur_win_mid] = cur_win_median;

				// Update the previous averaging window limits.
				prev_avg_start = cur_avg_start;
				prev_avg_end = cur_avg_end;
				cur_win_min = updated_win_min;
				cur_win_max = updated_win_max;

#undef _QSORT_CHECK_
#ifdef _QSORT_CHECK_
				// Get the median via qsort and compare as a sanity check.
				for(int i = cur_avg_start; i <= cur_avg_end; i++)
				{
					cur_win_vals[i-cur_avg_start] = track_data[i];
				} // i loop.
				qsort(cur_win_vals, l_averaging_win, sizeof(double), sort_doubles_descending);
				//per_win_profile[n_signal_wins] = cur_win_vals[l_averaging_win/2];

				if(cur_win_vals[l_averaging_win/2] != median_per_pdf)
				{
					fprintf(stderr, "Medians do not match: %d, %d:\n", (int)cur_win_vals[l_averaging_win/2], median_per_pdf);
					for(int i = cur_avg_start; i <= cur_avg_end; i++)
					{
						fprintf(stderr, "%lf ", track_data[i]);
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
				}
#endif // _QSORT_CHECK_

				n_signal_wins++;
			} // main signal filtering loop.

if(1)
{
			// Dump and free memory.
			char cur_decomp_fp[1000];
			sprintf(cur_decomp_fp, "filtered_data_scale_%d.bin", i_scale);
			fprintf(stderr, "Dumping %s\n", cur_decomp_fp);
			dump_per_nucleotide_binary_profile(cur_filtered_track, l_track_data, cur_decomp_fp);
			fprintf(stderr, "Closing file\n");
}

			// Free pdf memory.
			delete [] cur_win_vals;
			delete [] (cur_win_pdf+MIN_VAL);


if(__DUMP_FILTER_MSGS__)
			fprintf(stderr, "Getting the extrema regions.\n");

			// Get the extrema regions for the current filtered regions.
			vector<t_extrema_node*>* maxima = new vector<t_extrema_node*>();
			vector<t_extrema_node*>* minima = new vector<t_extrema_node*>();
			int* derivative_map = new int[l_track_data+2];
			memset(derivative_map, 0, sizeof(int) * (l_track_data+2));
			get_extrema_per_plateaus(cur_filtered_track, l_track_data, maxima, minima, derivative_map, 0);

			// Sort the extremas.
			sort(minima->begin(), minima->end(), sort_extremas_per_posn);
			sort(maxima->begin(), maxima->end(), sort_extremas_per_posn);

			// Get the minima regions.
			vector<t_annot_region*>* cur_scale_minima_regions = new vector<t_annot_region*>();

			fprintf(stderr, "Generating filtered minima regions from %ld minima.\n", minima->size());

			for(int i_m = 0; i_m < (int)minima->size()-1; i_m++)
			{
				t_annot_region* new_minima = get_empty_region();
				new_minima->chrom = t_string::copy_me_str("XX");
				new_minima->start = translate_coord(minima->at(i_m)->extrema_posn+1, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
				new_minima->end = translate_coord(minima->at(i_m+1)->extrema_posn+1, CODEBASE_COORDS::end_base, BED_COORDS::end_base);
				new_minima->strand = '+';

				// The maximum in the original signal is found within the effective data that contributes to the filtering process.
				int i_eff_data_start = (1 > new_minima->start-(int)scale/2)?(1):(new_minima->start-(int)scale/2);
				int i_eff_data_end = (l_track_data < new_minima->end+scale/2)?(l_track_data):(new_minima->end+(int)scale/2);
				double original_max = 0.0;
				for(int i = i_eff_data_start; i <= i_eff_data_end; i++)
				{
					if(original_max < track_data[i])
					{
						original_max = track_data[i];
					}
				}

				double filtered_max = 0.0;
				for(int i = new_minima->start; i <= new_minima->end; i++)
				{
					if(filtered_max < cur_filtered_track[i])
					{
						filtered_max = cur_filtered_track[i];
					}
				} // i loop.

if(__DUMP_FILTER_MSGS__)
				fprintf(f_feature_signal_vals, "XX\t%d\t%d\t%lf\t%lf\t%lf\n", new_minima->start, new_minima->end, filtered_max, original_max, scale);

				// This is where the filtering is done.
				if(filtered_max >= original_max / min_allowed_filtered_value_scale)
				{
					cur_scale_minima_regions->push_back(new_minima);
				}
				else
				{					
					delete new_minima;
				}
			} // i_min loop.

			fprintf(stderr, "%ld minima regions are generated.\n", cur_scale_minima_regions->size());
			per_scale_minima_regions->push_back(cur_scale_minima_regions);

			delete_extrema_nodes(minima);
			delete_extrema_nodes(maxima);
			
			// Free the current decomposition.
			delete [] derivative_map;
			delete [] cur_filtered_track;

			// Was the last processed value equal to or larger than what was requested?
			if(scale_end > 0 && 
				scale >= scale_end)
			{
				break;
			}

			i_scale++;
		} // scale computation check.
	} // scale loop.

if(__DUMP_FILTER_MSGS__)
		fclose(f_feature_signal_vals);
} // -get_filtered_maxima_regions_multiscale_filtered_data option.




















double* filter_data_per_odd_length_symmetric_filter(double* cur_real_track_data, int l_track_data, double* cur_filter_array, int l_filter, int n_iters)
{
#ifdef WIN32
	return(NULL);
#elif defined __unix__
	// Take the FFT of track outside the loop. Copy the track everytime.
	//int l_pre_ext = 1000*1000;
	int l_pre_ext = l_track_data/20;
	int l_ext_buf = l_track_data+l_pre_ext+l_pre_ext;
	int larger_exp_val = 0;
	int expon = 0;
	get_next_2_exp(l_ext_buf, larger_exp_val, expon);
	l_ext_buf = larger_exp_val;
	int l_ext = (l_ext_buf - l_track_data) / 2;

	if(l_ext > l_track_data)
	{
		fprintf(stderr, "Must increase the track data length, it is shorter than the extension length.\n");
		exit(0);
	}

	// Set up the workspace and wavetable.
	gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(l_ext_buf);
	gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(l_ext_buf);

//if(__DUMP_FILTER_MSGS__)
//	fprintf(stderr, "Finished FFT'ing buffered track data.\n");

	double* cur_filter_params = new double[2 * l_ext_buf];
	memset(cur_filter_params, 0, sizeof(double) * 2 * l_ext_buf);
	for(int i = 0; i < l_filter; i++)
	{
		cur_filter_params[2 * i] = cur_filter_array[i];
	} // i loop.

	// Get the FFT of the data.
	gsl_fft_complex_forward(cur_filter_params, 1, l_ext_buf,
						wavetable, workspace);

	double* buffered_track_data = new double[2 * l_ext_buf];
	memset(buffered_track_data, 0, sizeof(double) * 2 * l_ext_buf);

	// Do mirror image extension at the ends with the preset extension length.
	for(int i = 0; i < l_ext; i++)
	{
			buffered_track_data[2*i] = cur_real_track_data[l_ext-i-1];
	}
	for(int i = l_ext; i < l_ext + l_track_data; i++)
	{
			buffered_track_data[2*i] = cur_real_track_data[i-l_ext];
	}
	for(int i = l_ext+l_track_data; i < l_ext+l_track_data+l_ext; i++)
	{
			buffered_track_data[2*i] = cur_real_track_data[l_track_data-(i-l_ext-l_track_data)-1];
	}

	// Do FFT of the mirror image extended signal track.
	gsl_fft_complex_forward (buffered_track_data, 1, l_ext_buf,
						wavetable, workspace);

	// Filter the data n_iters times.
	double* filtered_signal_fft = new double[2 * l_ext_buf];
	double* cur_filtered_track_data = new double[l_track_data+1];
	for(int i_iter = 0; i_iter < n_iters; i_iter++)
	{
		// Multiply the filter and the signal track fft's.
		//fprintf(stderr, "Retrieving real and imaginary parts of the track and filter data.\n");
		fprintf(stderr, "%d. filtering iteration.           \r", i_iter);
		
		for(int i = 0; i < l_ext_buf; i+=1)
		{
			double cur_real_val1 = buffered_track_data[2*i] * cur_filter_params[2*i] - (buffered_track_data[2*i+1] * cur_filter_params[2*i+1]);
			double cur_imag_val1 = buffered_track_data[2*i+1] * cur_filter_params[2*i] + (buffered_track_data[2*i] * cur_filter_params[2*i+1]);
			filtered_signal_fft[2*i] = cur_real_val1;
			filtered_signal_fft[2*i+1] = cur_imag_val1;
		} // i loop.

		// Take the inverse FFT.
		gsl_fft_complex_inverse(filtered_signal_fft, 1, l_ext_buf, wavetable, workspace);

		// Copy the valid region.
		int half_l_filter = (l_filter-1)/2; 
		int i_start = half_l_filter + l_ext;
		int i_end = i_start + l_track_data;

		// Copy the real part of the filtered data back.
		for(int i = i_start; i < i_end; i++)
		{
			cur_filtered_track_data[i-i_start] = filtered_signal_fft[2*i];

			if(filtered_signal_fft[2*i+1] > 0.01)
			{
				fprintf(stderr, "FFT problem.\n");
				exit(0);
			}
		} // i loop.

		// Do mirror image extension at the ends with the preset extension length.
		memset(buffered_track_data, 0, sizeof(double) * 2 * l_ext_buf);
		for(int i = 0; i < l_ext; i++)
		{
			buffered_track_data[2*i] = cur_filtered_track_data[l_ext-i-1];
		}
		for(int i = l_ext; i < l_ext + l_track_data; i++)
		{
			buffered_track_data[2*i] = cur_filtered_track_data[i-l_ext];
		}
		for(int i = l_ext+l_track_data; i < l_ext+l_track_data+l_ext; i++)
		{
			buffered_track_data[2*i] = cur_filtered_track_data[l_track_data-(i-l_ext-l_track_data)-1];
		}

		// Do FFT of the mirror image extended signal track.
		gsl_fft_complex_forward (buffered_track_data, 1, l_ext_buf,
							wavetable, workspace);
	} // i_iter loop.

	delete [] cur_filter_params;
	delete [] filtered_signal_fft;

	// Free workspace memory.
	gsl_fft_complex_wavetable_free (wavetable);
	gsl_fft_complex_workspace_free (workspace);

	delete [] buffered_track_data;
	
	return(cur_filtered_track_data);
#endif // __unix__ check.	
}
























bool get_origin_of_symmetry(double* filtered_data, int l_signal, int search_start_i, int& origin, int l_mirror_check)
{
	if(search_start_i < l_mirror_check)
	{
		search_start_i = l_mirror_check;
	}

	for(int i_c = l_mirror_check; i_c < l_signal-l_mirror_check; i_c++)
	{
		bool current_origin_check = true;
		for(int i = 0; 
			current_origin_check && i < l_mirror_check; 
			i++)
		{
			if(fabs(filtered_data[i_c+i] - filtered_data[i_c-i]) > 0.00001)
			{
				current_origin_check = false;
			}
		} // i loop.

		if(current_origin_check)
		{
			fprintf(stderr, "Found origin of symmetry @ %d\n", i_c);
			origin = i_c;
			return(true);
		}
	} // i_c loop.

	return(false);
}

double* sum_tracks(vector<double*>* multitrack_data, int l_signal)
{
	int n_tracks = (int)multitrack_data->size();
	fprintf(stderr, "Summing %d tracks of %d data points.\n", n_tracks, l_signal);
	double* aggregate_track = new double[l_signal];

	for(int i = 0; i < l_signal; i++)
	{
		double cur_aggregated_val = 0.0;

		// Go over all the tracks and sum their values at this position.
		for(int i_t = 0; i_t < n_tracks; i_t++)
		{
			cur_aggregated_val += (multitrack_data->at(i_t))[i];
		} // i_t loop.

		aggregate_track[i] = cur_aggregated_val;
	} // i loop.

	return(aggregate_track);
}

double* conv_filter_sub_data(double* data, int l_data, int start, int end, 
	double* filter_params, int l_filter,
	int& l_filtered_signal_length, int& i_filtered_data_start)
{
	int l_sub_data = (end-start+1);
	int l_ext = l_filter;
	int l_ext_buf = l_ext+l_ext+l_sub_data;

	double* buffered_track_data = new double[l_ext_buf];

    for(int i = 0; i < l_ext; i++)
    {
		buffered_track_data[i] = 0;
	}

	for(int i = l_ext; i < l_ext + l_sub_data; i++)
	{
		buffered_track_data[i] = data[i+start-l_ext];
	}

	for(int i = l_ext+l_sub_data; i < l_ext+l_sub_data+l_ext; i++)
	{
		buffered_track_data[i] = 0;
	}

	l_filtered_signal_length = l_ext_buf;
	i_filtered_data_start = l_ext + l_filter/2;

	double* filtered_track_data = new double[l_ext_buf];
	memset(filtered_track_data, 0, sizeof(double) * l_ext_buf);
	for(int i_d = 0; i_d < l_ext_buf; i_d++)
	{
		double cur_filtered_val = 0;

		// Compute the filtered value at i^{th} val.
		for(int i_c = 0; i_c < l_filter; i_c++)
		{
			if(i_d > i_c)
			{
				cur_filtered_val += (filter_params[i_c] * buffered_track_data[i_d - i_c]);
			}
		} // i_c loop.

		filtered_track_data[i_d] = cur_filtered_val;
	} // i loop.

	delete [] buffered_track_data;

	return(filtered_track_data);
}

double* filter_sub_data(double* data, int l_data, int start, int end, 
	double* filter_params, int l_filter,
	int& l_filtered_signal_length, int& i_filtered_data_start)
{
if(__DUMP_FILTER_MSGS__)
	fprintf(stderr, "Filtering data of length %d subdata from %d to %d.\n", l_data, start, end);

	// Determine the extension length.
	int l_sub_data = (end-start+1);
	//int l_pre_ext = l_filter;
	int l_pre_buf_ext = l_filter+l_filter+l_sub_data;

	int larger_exp_val = 0;
	int expon = 0;
	get_next_2_exp(l_pre_buf_ext, larger_exp_val, expon);
	int l_ext_buf = larger_exp_val;
	int l_ext = (l_ext_buf - l_sub_data) / 2;

#ifdef WIN32
	double* sub_filtered_data = conv_filter_sub_data(data, l_data, start, end, filter_params, l_filter, l_filtered_signal_length, i_filtered_data_start);
	return(sub_filtered_data);
#elif defined __unix__
	// Set up the workspace and wavetable.
	gsl_fft_complex_wavetable* wavetable = gsl_fft_complex_wavetable_alloc(l_ext_buf);
	gsl_fft_complex_workspace* workspace = gsl_fft_complex_workspace_alloc(l_ext_buf);

	// Generate the buffer with the extension length: Append 0's.
	double* ext_data_buffer = new double[2 * l_ext_buf];
	memset(ext_data_buffer, 0, sizeof(double) * 2 * l_ext_buf);
	for(int i = l_ext; i < l_ext+l_sub_data; i++)
	{
		ext_data_buffer[2*i] = data[start+i-l_ext];
		ext_data_buffer[2*i+1] = 0;
	} // i loop.

	double* ext_filter_buffer = new double[2*l_ext_buf];
	memset(ext_filter_buffer, 0, sizeof(double) * 2 * l_ext_buf);
	for(int i = 0; i < l_filter; i++)
	{
		ext_filter_buffer[2*i] = filter_params[i];
		ext_filter_buffer[2*i+1] = 0;
	} // i loop.

	// Take the FFT's.
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "FFT'ing extended data.\n");
		gsl_fft_complex_forward(ext_data_buffer, 1, l_ext_buf,
							wavetable, workspace);
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Done FFT'ing extended.\n");

if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "FFT'ing filter parameters.\n");
		gsl_fft_complex_forward(ext_filter_buffer, 1, l_ext_buf,
							wavetable, workspace);
if(__DUMP_FILTER_MSGS__)
		fprintf(stderr, "Done FFT'ing filter parameters.\n");

	// Multiply.
	double* filtered_signal_fft = new double[2 * l_ext_buf];
	for(int i = 0; i < l_ext_buf; i+=1)
	{
		double cur_real_val1 = ext_data_buffer[2*i] * ext_filter_buffer[2*i] - (ext_data_buffer[2*i+1] * ext_filter_buffer[2*i+1]);
		double cur_imag_val1 = ext_data_buffer[2*i+1] * ext_filter_buffer[2*i] + (ext_data_buffer[2*i] * ext_filter_buffer[2*i+1]);
		filtered_signal_fft[2*i] = cur_real_val1;
		filtered_signal_fft[2*i+1] = cur_imag_val1;
	} // i loop.

	// Take IFFT's.
if(__DUMP_FILTER_MSGS__)
	fprintf(stderr, "IFFT'ing of multiplication.\n");

	gsl_fft_complex_inverse(filtered_signal_fft, 1, l_ext_buf, wavetable, workspace);

if(__DUMP_FILTER_MSGS__)
	fprintf(stderr, "Done IFFT'ing the multiplication.\n");	

	// Set the correct indices and return the result.
	//t_ext_double_array* filtered_data = alloc_ext_array();
	l_filtered_signal_length = l_ext_buf;
	i_filtered_data_start = l_filter/2 + l_ext;
	int i_start = 0;
	int i_end = l_ext_buf;
	double* cur_filtered_track = new double[l_ext_buf];

	FILE* f_signal = NULL;
if(__DUMP_FILTER_MSGS__)
	f_signal = fopen("filtered_signal.txt", "w");

	for(int i = i_start; i < i_end; i++)
	{
		//add_val(filtered_data, filtered_signal_fft[2*i]);
		cur_filtered_track[i] = filtered_signal_fft[2*i];

		// Make sure that the imaginary value is 0.
		if(fabs(filtered_signal_fft[2*i+1] > 0.000001))
		{
			fprintf(stderr, "IFFT problem.\n");
			exit(0);
		}

if(__DUMP_FILTER_MSGS__)
		fprintf(f_signal, "%.10f %10f\n", REAL(filtered_signal_fft, i), IMAG(filtered_signal_fft, i));

	} // i loop.

if(__DUMP_FILTER_MSGS__)
	fclose(f_signal);

	// Copy the data.
	//double* cur_filtered_track = filtered_data->buffer;
	//delete(filtered_data);

	// Free workspace memory.
	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
	delete [] ext_filter_buffer;
	delete [] ext_data_buffer;
	delete [] filtered_signal_fft;

	return(cur_filtered_track);
#endif // __unix__
}
