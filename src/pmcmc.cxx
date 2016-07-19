#include "pmcmc.h"

void sample_smc_approximation( gsl_rng *r, 
                               int K,
                               int N,
                               int T,
                               Tensor<double, 4> &x, 
                               Tensor<double, 2> &w,
                               Tensor<double, 3> &X )
{

    double p[N];
    
    for(int i=0; i<N; ++i)
    {
        p[i] = w(T-1, i);
    }
    
    unsigned int particle_index;
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, p);
    particle_index = gsl_ran_discrete(r, g);
    
    for(int k=0; k<K; ++k)
    for(int t=0; t<T; ++t)
    {
        X(k, t,0) = x(k, t, particle_index, 0);
        X(k, t,1) = x(k, t, particle_index, 1);
    }

}

double calculate_log_acceptance_probability( double log_marginal_likelihood_c,
                                             double log_marginal_likelihood_c_star,
                                             Tensor<double, 1> &C, 
                                             Tensor<double, 1> &C_star,
                                             double proposal_c_variance )
{
    return  log_marginal_likelihood_c_star-log_marginal_likelihood_c\
            +(0.5/proposal_c_variance)*( C(0)*C(0) - C_star(0)*C_star(0) );
}
