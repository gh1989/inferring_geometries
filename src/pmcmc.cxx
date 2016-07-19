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
        p[i] = exp( w(T-1, i) );
    
    unsigned int particle_index;
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, p);
    particle_index = gsl_ran_discrete(r, g);
    
    for(int k=0; k<K; ++k)
    for(int t=0; t<T; ++t)
    {
        X(k, t, 0) = x(k, t, particle_index, 0);
        X(k, t, 1) = x(k, t, particle_index, 1);
    }

}

double calculate_log_acceptance_probability( double log_marginal_likelihood_c,
                                             double log_marginal_likelihood_c_star,
                                             Tensor<double, 1> &C, 
                                             Tensor<double, 1> &C_star,
                                             double proposal_c_variance, 
                                             int parameters )
{

    double log_prob = 0;
    double proposal_constant = 0.5/proposal_c_variance;
    
    log_prob += log_marginal_likelihood_c_star;
    log_prob -= log_marginal_likelihood_c;
    
    for( int i=0; i<parameters; ++i )
    {
       log_prob += proposal_constant*( C(i)*C(i) );
       log_prob -= proposal_constant*( C_star(i)*C_star(i) );
    }
    return log_prob;
}
