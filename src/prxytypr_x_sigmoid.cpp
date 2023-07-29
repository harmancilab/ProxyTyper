#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_file_utils.h"
#include "prxytypr_histogram.h"
#include <math.h>
#include <string.h>

#include "prxytypr_x_sigmoid.h"

bool __DUMP_X_SIGMOID_MESSAGES__ = false;

double get_log_val_protected(double param)
{
	if (param > 0)
	{
		return(log(param));
	}
	else
	{
		fprintf(stderr, "Encountered negatives for log, returning 0!\n");
		return(0);
	}
}

double get_sigmoid_val_per_feat_comb(double feat_comb)
{
	double sig_val = exp(-1 * feat_comb);

	return(1.0 / (1.0 + sig_val));
}

double get_logit(double p)
{
	return(log(p / (1 - p)));
}

double get_sigmoid_val_per_feats_weights(double* feats, int n_feats, double* per_feat_weights)
{
	double feat_comb = 0;
	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		feat_comb += feats[feat_i] * per_feat_weights[feat_i];
	} // feat_i loop.

	double sig_val = exp(-1 * feat_comb);

	return(1.0 / (1.0 + sig_val));
}

double get_TenSeal_poly_approx_sigmoid_per_feat_comb(double feat_comb)
{
	if (feat_comb > 5)
	{
		return(1.0);
	}
	else if (feat_comb < -5)
	{
		return(0);
	}

	// https://github.com/OpenMined/TenSEAL/blob/main/tutorials/Tutorial%201%20-%20Training%20and%20Evaluation%20of%20Logistic%20Regression%20on%20Encrypted%20Data.ipynb
	double sig_approx = 0.5 + 0.197 * feat_comb - 0.004 * (feat_comb * feat_comb * feat_comb);

	return(sig_approx);
}

double get_TenSeal_poly_approx_sigmoid_per_feats_weights(double* feats, int n_feats, double* per_feat_weights)
{
	double feat_comb = 0;
	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		feat_comb += feats[feat_i] * per_feat_weights[feat_i];
	} // feat_i loop.

	double sig_approx = get_TenSeal_poly_approx_sigmoid_per_feat_comb(feat_comb);

	if (__DUMP_X_SIGMOID_MESSAGES__)
	{
		fprintf(stderr, "TenSeal approx(%lf)=%lf\n", feat_comb, sig_approx);
	}
	
	return(sig_approx);
}

double get_Kim_etal_poly_approx_sigmoid_per_feats_weights(double* feats, int n_feats, double* per_feat_weights)
{
	double feat_comb = 0;
	for (int feat_i = 0; feat_i < n_feats; feat_i++)
	{
		feat_comb += feats[feat_i] * per_feat_weights[feat_i];
	} // feat_i loop.

	double sig_approx = get_Kim_etal_poly_approx_sigmoid_per_feat_comb(feat_comb);
	return(sig_approx);
}

double get_Kim_etal_poly_approx_sigmoid_per_feat_comb(double feat_comb)
{
	// https://eprint.iacr.org/2018/074.pdf
	//double a1 = 0.25;
	//double a3 = -1.0 * (double)1.0 / 48;
	//double a5 = (double)1.0 / 480;
	//double a7 = -1.0 * (double)17.0 / 80640;
	//double a9 = (double)31.0 / 1451520;

	//double sig_approx = 0.5 + a1 * feat_comb +
	//					a3 * pow(feat_comb, 3) +
	//					a5 * pow(feat_comb, 5) + 
	//					a7 * pow(feat_comb, 7) + 
	//					a9 * pow(feat_comb, 9);

	double a0 = 0.5;
	double a1 = 1.73496;
	double a3 = -1.0 * 4.19407;
	double a5 = 5.43402;
	double a7 = -1.0 * 2.50739;

	feat_comb /= 8;

	double sig_approx = a0 +
		a1 * feat_comb +
		a3 * pow(feat_comb, 3) +
		a5 * pow(feat_comb, 5) +
		a7 * pow(feat_comb, 7);

	return(sig_approx);
}
