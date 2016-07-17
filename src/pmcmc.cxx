#include "pmcmc.h"

void sample_smc_approximation( gsl_rng *r, 
                               int N,
                               int T,
                               Tensor<double, 3> &x, 
                               Tensor<double, 2> &w,
                               Tensor<double, 2> &X )
{

    double p[N];
    
    for(int i=0; i<N; ++i)
    {
        p[i] = w(T-1, i);
    }
    
    unsigned int particle_index;
    gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, p);
    particle_index = gsl_ran_discrete(r, g);
    
    for(int t=0; t<T; ++t)
    {
        X(t,0) = x(t, particle_index, 0);
        X(t,1) = x(t, particle_index, 1);
    }

}

double calculate_log_acceptance_probability( double marginal_likelihood_c,
                                             double marginal_likelihood_c_star,
                                             Tensor<double, 1> &C, 
                                             Tensor<double, 1> &C_star,
                                             double proposal_c_variance )
{
    return (log(marginal_likelihood_c_star) - log(marginal_likelihood_c)) +\
            (log_prior_c(C_star, proposal_c_variance) - log_prior_c(C, proposal_c_variance)) +\
            (log_q( C, C_star, proposal_c_variance ) - log_q( C_star, C, proposal_c_variance));
}

double log_q(   Tensor<double, 1> &C, 
                Tensor<double, 1> &C_star,
                double proposal_c_variance)
{
    double diff = C(0)-C_star(0);
    return -(0.5/proposal_c_variance)*diff*diff;
}

double log_prior_c( Tensor<double, 1> &C,
                    double proposal_c_variance )
{
    return -(0.5/proposal_c_variance)*C(0)*C(0);
}