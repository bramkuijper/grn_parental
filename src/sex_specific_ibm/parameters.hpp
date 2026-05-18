#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>
#include <vector>

class Parameters
{
    public:
        // data is written to file each nth generation
        unsigned data_output_interval{1};

        // population size
        unsigned N{5000};

        // population size when sampling individuals 
        // to assess genetic canalization
        unsigned N_genetic_canalization{1000};
        
        // population size when sampling individuals 
        // to assess environmental canalization
        unsigned N_environmental_canalization{1000};

        // number of genes in network
        unsigned L{5};

        // index of the sex-specific locus
        unsigned sex_specific_locus_idx;

        double sk[2]{0.0,0.0};

        // the locus of the sex-determining 
        // master gene, 
        // ideally this is not one of the loci under selection
        unsigned master_gene;

        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        // [Wij drawn on N(0,0.1)]
        double sd_init_strength_w{0.1};

        // initial value of each element of the phenotype
        double s_init{0.0};

        // constitutive gene expression value
        // see Odorico et al 2018 p 689 1st column, 1st para
        // in the absence of any GRN, this is the value
        // attained by each element of S
        double a{0.2};

        // proportion inheritance that is nongenetic
        double p_nongenetic{0.0};

        // if proportion inheritance > 0
        // then p_maternal determines proportion inherited
        // from mother vs inherited from father
        //
        // alternative ideas would be to designate specific
        // genes to maternal vs paternal transmission
        double p_maternal{0.0};

        // maximum number of time steps that simulation lasts
        unsigned long max_time_step{10};
        
        // number of developmental time steps
        // has to be larger than 0
        unsigned max_dev_time_step{10};
       
        // number of developmental time steps in case
        // of environmental canalization
        unsigned max_dev_time_step_envt_canalize{100};

        // maximum developmental time step when we test for
        // stability, we want this to be as long as possible
        unsigned max_dev_time_step_stability{10000};

        // number of developmental time steps
        // recorded for stats from time T backwards. 
        // These are the final
        // n time steps before development has been
        // completed.
        // this parameter is used to calculate
        // the values of Sbar_i (mean gene expression for locus i) 
        // and V_i (the variance of gene expression)
        //
        // a value of 1 is basically the current value
        unsigned max_dev_time_step_stats{2};

        // mutation rate of W matrix entries
        double mu_w{0.01};
       
        // standard deviation of W matrix entries
        double sdmu_w{0.5};

        // optimum for each locus for females
        std::vector < std::vector < double > > theta;

        // strength of selection on matching the optimum
        // see Odorico et al. eq (4)
        // we are going to make this trait-specific, so that
        // we can turn off selection on specific traits
        std::vector <double> s;

        double baseline_fitness{1.0};

        // threshold below which gene expression at a 
        // particular locus is considered the same or different
        // relative to that of a reference individual
        // see Odorico p690 1st col 3rd para
        double canalization_threshold{0.01};

        // stability threshold, i.e., a difference S_i(t+1) - S_i(t)
        // smaller than this threshold means that we have approximately
        // S_i(t+1) = S_i(t)
        double stability_threshold{0.001};

        // strength of selection on variance 
        // of gene expression
        // see Odorico et al. eq (4)
        // apparently this is locus specific
        std::vector <double> sprime;

        std::string file_name{"grn_simulation.csv"};
        std::string file_name_individuals{"grn_simulation_individuals.csv"};

        Parameters(unsigned const Lval) : 
            L(Lval),
            sex_specific_locus_idx(Lval - 1),
            theta(2, std::vector < double > (Lval,0.0)), // initialize all the values for theta
            s(L, 0.0), // initialize all the values for s, the strength of selection on each trait (0 = no selection)
            sprime(L, 0.0) // initialize all the values for s, the strength of selection on each trait (0 = no selection)
        {}
};


#endif 
