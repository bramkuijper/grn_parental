#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <vector>
#include <random>
#include "parameters.hpp"

class Individual
{
    private:
        // reference to parameter object. This assumes that the parameter
        // object that is referenced lives longer than the individual
        // otherwise undefined behaviour. The assumption is valid tho,
        // as the original parameter object lives for the duration of the
        // simulation.
        //
        // We don't want a non-reference parameter object in each and every
        // individual because it would mean we have to copy all its data members
        // upon the creation of each and every individual
        Parameters const &pars;
    public:
        // the GRN represented by a 2-dimensional matrix
        // each element is a haploid allele that is inherited
        // according to classical Mendelian inheritance in 
        // a sexual population
        //
        // element W_12 can be thought of as the influence
        // of current transcription (at time step t) of Gene_2 (measured by
        // S_2(t)) on the transcription of Gene_1 in the next
        // time step t+1 , measured by S_1(t+1)
        std::vector < std::vector < double > > W{};

        // the eventual phenotype,
        // each row here is a phenotype vector at time t = 0,
        // 1, 2, ..., max_tS.
        // where the vector S[0] is the initial
        // vector of gene expression levels at birth.
        // By contrast, the vector S[max_tS - 1] 
        // is the most recent
        // phenotype vector and the one used for fitness 
        // calculations.  
        std::vector < std::vector < double > > S{};

        // average phenotype where Sbar[i] is the time-average
        // phenotype of trait i
        std::vector <double> Sbar{};
        // variance in phenotype where Vbar[i] is the variance
        // of phenotype of trait i, measured over time
        std::vector <double> V{};
       
        // initialisation constructor, 
        // used at the start of the simulation
        Individual(Parameters const &par);

        // copy constructor, which is used when you assign
        // individuals to vectors and for other vector-related operations
        Individual(Individual const &other);

        // birth constructor 
        Individual(Individual const &mother, 
                Individual const &father,
                Parameters const &params, 
                std::mt19937 &rng_r);

        void operator=(Individual const &other);

        // the sigmoid function to relate a raw gene expression value (x)
        // to a value of 0 and 1. See Odorico et al 2018 J Evol Biol eq. (2) and 
        // Runneburger and Le Rouzic (2016) BMC Evol Biol eq. (2)
        double gene_expression_sigmoid(double const x, double const a);

        // update the individual's S vector, which contains the expression
        // levels of L genes
        void update_phenotype(unsigned const dev_time_step);

        void average_phenotype();

}; // end class individual

#endif
