#include "Param.h"
#include "Lib.h"

void InitParam(Param* param, double ntree, double kfac, double sigdf, double sigquant){
    param->m_bart = ntree;
    param->k_bart = kfac;
    param->mu_prec = ntree * 4 * kfac * kfac;
    param->alpha_bart = sigdf;

    double mean = 0, d;
    double var = 0;
    for(int i = 1; i <= NumObs; i++)
        mean += YDat1[i];
    mean /= NumObs;
    for(int i = 1; i <= NumObs; i++){
        d = YDat1[i] - mean;
        var += d * d;
    }
    var /= NumObs;
    double prec = 1 / var;
    param->lambda_bart = 2 * prec;
    param->log_lambda_bart = log(param->lambda_bart);
    param->beta_bart = compute_gamma_param(prec, sigdf, sigquant);

    param->nn_prior_term = 0.5 * log(param->mu_prec);
    double pi = 3.14159265359;
    param->half_log_2pi = 0.5 * log(2 * pi);

    //Rprintf("prec = %f\n", prec);
    //Rprintf("beta_bart = %f\n", param->beta_bart);
    //Rprintf("lambda_bart = %f\n", param->lambda_bart);

}

double compute_gamma_param(double min_val, double alpha, double q, double init_val){
    if (init_val < 0) {
      init_val = alpha / 3.0 / min_val;
    }
    double f1, f2 = 0;
    double x1 = init_val;
    double x2 = init_val;
    do {
      if (x2 < 0){
        x1 = 0;
      }
      else {
        x1 = x2;
      }
      f1 = pgamma(min_val, alpha, 1 / x1, 0, 0) - q;
      f2 = (pgamma(min_val, alpha, 1 / (x1 + 0.0005), 0, 0) - pgamma(min_val, alpha, 1 / x1, 0, 0)) / 0.0005;
      x2 = x1 - f1 / f2;

    } while ((x1 - x2 > 1.49012e-8) || (x2 - x1 > 1.49012e-8));
    try{
      if (abs(pgamma(min_val, alpha, 1 / x1, 0, 0) - q) > 1e-3) {
        throw x1;
      }
    }
    catch (double) {
      double new_init = 0.001 > init_val * 0.9 ? 0.001 : init_val * 0.9;
      x1 = compute_gamma_param(min_val, alpha, q, new_init);
    }
    return x1;
}
