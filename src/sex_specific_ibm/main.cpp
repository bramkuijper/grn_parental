//#include <iostream>
#include "parameters.hpp"
#include "gene_network_matpat.hpp"

int main(int argc, char **argv)
{
    size_t argv_idx{1};
    // we always first need to retrieve the number of loci
    // with which we can instantiate the parameter object
    unsigned Lval = static_cast<unsigned>(std::stoi(argv[argv_idx++]));
    Parameters pars(Lval);

    for (size_t theta_idx{0}; theta_idx < pars.L; ++theta_idx)
    {
        pars.theta[0][theta_idx] = std::stod(argv[argv_idx++]);
        pars.theta[1][theta_idx] = std::stod(argv[argv_idx++]);
    }

    pars.max_time_step = static_cast<unsigned>(std::stoi(argv[argv_idx++]));
    pars.p_nongenetic = std::stod(argv[argv_idx++]);
    pars.p_maternal = std::stod(argv[argv_idx++]);
    pars.max_dev_time_step = static_cast<unsigned>(std::stoi(argv[argv_idx++]));

    for (size_t s_idx{0}; s_idx < pars.L; ++s_idx)
    {
        pars.s[s_idx] = std::stod(argv[argv_idx++]);
    }
    
    for (size_t s_idx{0}; s_idx < pars.L; ++s_idx)
    {
        pars.sprime[s_idx] = std::stod(argv[argv_idx++]);
    }

    pars.data_output_interval = static_cast<unsigned>(std::stoi(argv[argv_idx++]));
    pars.file_name = argv[argv_idx++];
    pars.file_name_individuals = argv[argv_idx++];

    GRN_MatPat Simulation(pars);

    return 0;
}
