#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "prxytypr_rng.h"
#include "prxytypr_seed_manager.h"
#include "prxytypr_file_utils.h"
#include "prxytypr_histogram.h"
#include <math.h>
#include <string.h>
#include "prxytypr_x_chisqr.h"

double approx_gamma(double Z)
{
	const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
	const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

	double D = 1.0 / (10.0 * Z);
	D = 1.0 / ((12 * Z) - D);
	D = (D + Z) * RECIP_E;
	D = pow(D, Z);
	D *= sqrt(TWOPI / Z);

	return D;
}

double x_gamma(double N)
{
	const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;

	long double A = 5;
	long double Z = (long double)N;
	long double Sc = powl((Z + A), (Z + 0.5));
	Sc *= expl(-1.0 * (Z + A));
	Sc /= Z;

	long double F = 1.0;
	long double Ck;
	long double Sum = SQRT2PI;


	for (int K = 1; K < A; K++)
	{
		Z++;
		Ck = powl(A - K, K - 0.5);
		Ck *= expl(A - K);
		Ck /= F;

		Sum += (Ck / Z);

		F *= (-1.0 * K);
	}

	return (double)(Sum * Sc);
}

static double igf(double S, double Z)
{
	if (Z < 0.0)
	{
		return 0.0;
	}
	double Sc = (1.0 / S);
	Sc *= pow(Z, S);
	Sc *= exp(-Z);

	double Sum = 1.0;
	double Nom = 1.0;
	double Denom = 1.0;

	for (int I = 0; I < 200; I++)
	{
		Nom *= Z;
		S++;
		Denom *= S;
		Sum += (Nom / Denom);
	}

	return Sum * Sc;
}

double chisqr(int Dof, double Cv)
{
	if (Cv < 0 || Dof < 1)
	{
		return 0.0;
	}
	double K = ((double)Dof) * 0.5;
	double X = Cv * 0.5;
	if (Dof == 2)
	{
		return exp(-1.0 * X);
	}

	double PValue = igf(K, X);
	if (isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
	{
		return 1e-14;
	}

	PValue /= x_gamma(K);
	//PValue /= tgamma(K); 

	return (1.0 - PValue);
}