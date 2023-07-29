#ifndef __X_SIGMOID__
#define __X_SIGMOID__

double get_logit(double p);

double get_log_val_protected(double param);

double get_sigmoid_val_per_feats_weights(double* feats, int n_feats, double* per_feat_weights);
double get_sigmoid_val_per_feat_comb(double feat_comb);

double get_Kim_etal_poly_approx_sigmoid_per_feats_weights(double* feats, int n_feats, double* per_feat_weights);
double get_Kim_etal_poly_approx_sigmoid_per_feat_comb(double feat_comb);

double get_TenSeal_poly_approx_sigmoid_per_feats_weights(double* feats, int n_feats, double* per_feat_weights);
double get_TenSeal_poly_approx_sigmoid_per_feat_comb(double feat_comb);

#endif // __X_SIGMOID__
