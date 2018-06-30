#ifndef PARAM_HEAD
#define PARAM_HEAD

#include <R.h>
#include <Rmath.h>

#include "global.h"

class Param{
public:
    //param
    double m_bart;
    double k_bart;
    double mu_prec;
    double alpha_bart;
    double beta_bart;
    double lambda_bart;
    double log_lambda_bart;

    double nn_prior_term;
    double half_log_2pi;
};

void InitParam(Param* param, double ntree, double kfac, double sigdf, double sigquant);

double compute_gamma_param(double min_val, double alpha, double q, double init_val = -1);

extern Param param;

#endif
