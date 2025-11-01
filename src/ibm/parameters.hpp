#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>
#include <vector>

class Parameters
{
    public:
        // data is written to file each nth generation
        unsigned data_output_interval{10};

        // population size
        unsigned N{1000};

        unsigned N_canalize{1000};

        // number of genes in network
        unsigned L{5};

        // initial value of each node in the grn
        double w_init{0.0};

        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        // [Wij drawn on N(0,0.1)]
        double sd_init_strength_w{0.1};

        // initial value of each element of the phenotype
        double s_init{0.0};

        // constitutive gene expression value
        // in the absence of any GRN, this is the value
        // attained by each element of S
        double a{0.0};

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

        // number of developmental time steps
        // recorded for stats. These are the final
        // n time steps before development has been
        // completed.
        unsigned max_dev_time_step_stats{4};

        // mutation rate of W matrix entries
        double mu_w{0.02};
        
        // standard deviation of W matrix entries
        double sdmu_w{0.02};

        // optimum for each locus
        std::vector <double> theta;

        // strength of selection on matching the optimum
        // see Odorico et al. eq (4)
        // we are going to make this trait-specific, so that
        // we can turn off selection on specific traits
        std::vector <double> s;


        // strength of selection on variance 
        // of gene expression
        // see Odorico et al. eq (4)
        double sprime{46000.0};

        std::string file_name{"grn_simulation.csv"};
        std::string file_name_individuals{"grn_simulation_individuals.csv"};

        Parameters() :
            theta(L, 0.0), // initialize all the values for theta
            s(L, 0.0) // initialize all the values for s
        {}
};


#endif 
