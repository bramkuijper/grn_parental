#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

class Parameters
{
    public:
        // data is written to file each nth generation
        unsigned data_output_interval{10};

        // population size
        unsigned N{1000};

        // number of genes in network
        unsigned L{5};

        std::string file_name{"grn_simulation.csv"};
};


#endif _PARAMETERS_HPP_
