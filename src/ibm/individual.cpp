#include <vector>
#include <random>
#include <cassert>
#include "individual.hpp"

Individual::Individual(Parameters const &par) :
    pars(par),
    W(par.L, std::vector< double >(par.L) ),
    S(par.L, std::vector< double >(par.L) )
{} // end Individual() initialization constructor

Individual::Individual(Individual const &other) :
    pars(other.pars),
    W(other.W),
    S(other.S)
{} // end copy constructor

// the birth constructor: this function is used to
// build an individual out of its two parents
Individual::Individual(Individual const &mom,
        Individual const &dad,
        Parameters const &par,
        std::mt19937 &rng_r) :
    pars(par)
{
    // set up distribution functions
    std::uniform_real_distribution uniform{0.0,1.0};
    std::normal_distribution<double> standard_normal{0.0,1.0};

    // inherit GRN gene loci
    for (unsigned row_idx{0}; row_idx < pars.L; ++row_idx)
    {
        for (unsigned col_idx{0}; col_idx < pars.L; ++row_idx)
        {
            // Mendelian transmission of a haploid locus
            // while Odorico et al 2018 study both diploidy and haploidy
            // let's start simple. 
            W[row_idx][col_idx] = uniform(rng_r) < 0.5 ? 
                dad.W[row_idx][col_idx] 
                :
                mom.W[row_idx][col_idx];

            // mutate the allele after inheritance
            if (uniform(rng_r) < pars.mu_w)
            {
                W[row_idx][col_idx] += standard_normal(rng_r) * pars.sdmu_w;
            }
        }
    } // end for unsigned
    
    // now biparental transmission of the gene expression
    // products, contained in the phenotype vector S
    for (unsigned int row_idx{0}; row_idx < pars.L; ++row_idx)
    {
        S[0][row_idx] = (1.0 - pars.p_nongenetic) * pars.a
            + pars.p_nongenetic * (
                    pars.p_maternal * mom.S[pars.max_dev_time_step - 1][row_idx]
                    +
                    (1.0 - pars.p_maternal) * dad.S[pars.max_dev_time_step - 1][row_idx]
                );
    } // end transmission of S
} // end birth constructor

// assignment operator overloaded
void Individual::operator=(Individual const &other)
{
    W = other.W;
    S = other.S;
} // end operator=()

void Individual::average_phenotype()
{
    // empty vector for the averages
    std::fill(Sbar.begin(),Sbar.end(),

    for (unsigned t_prior{pars.max_dev_time_step - 1 - pars.max_dev_time_step_stats}; 
            t_prior < pars.max_dev_time_step; ++t_prior)
    {
        assert(t_prior >= 0);
        assert(t_prior < S.size());
        assert(S[0].size() == pars.L);

        for (unsigned s_idx{0}; s_idx < pars.L; ++s_idx)
        {
            Sbar[s_idx] += S[t_prior][s_idx];
        } // end for s_idx
    } // end for t_prior
    

    for (unsigned s_idx{0}; s_idx < pars.L; ++s_idx)
    {
        Sbar[s_idx] += S[t_prior][s_idx] / 
            pars.max_dev_time_step_stats;
    } // end for s_idx
} // end function average phenotype

void Individual::update_phenotype(
        unsigned const dev_time_step)
{
    // TODO: sort out how to quickly loop over matrices in C++
    // TODO: reduce size of S by only focusing on current time step 
    // and then time steps for stats
    //
    // S(t+1) vector, all set to 0
    std::vector<double> Stplus1(pars.L, 0.0);

    for (unsigned int row_idx{0}; row_idx < pars.L; ++row_idx)
    {
        for (unsigned int col_idx{0}; col_idx < pars.L; ++col_idx)
        {
            S[dev_time_step + 1][row_idx] += gene_expression_sigmoid(
                    W[row_idx][col_idx] * S[dev_time_step][col_idx],
                    pars.a
                    );
        } // end for col_idx
    } // end for row_idx
} // update_phenotype()

// the sigmoid function to relate raw gene expression values (x)
// to a value of 0 and 1. See Odorico et al 2018 J Evol Biol eq. (2) and Runneburger
// and Le Rouzic (2016) BMC Evol Biol eq. (2)
// important to realise that x here is a scalar, but the result of 
// sum_{j = 0}^{L} w_ij S_j (i.e., the summed product of weights and current
// phenotype, see Runneburger & Le Rouzic (2016) eq. (1)
double Individual::gene_expression_sigmoid(double const x, double const a)
{
    double xtplus1 = a / (a + (1.0 - a) * std::exp(-x/(a * (1.0 - a))));

    return(xtplus1);
}// end gene expression sigmoid
