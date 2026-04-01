//#include <iostream>
#include "parameters.hpp"
#include "gene_network_matpat.hpp"

int main(int argc, char **argv)
{
    // we always first need to retrieve the number of loci
    // with which we can instantiate the parameter object
    unsigned Lval = static_cast<unsigned>(std::stoi(argv[1]));
    double theta = std::stod(argv[2]);
    Parameters pars(Lval, theta);

    pars.max_time_step = static_cast<unsigned>(std::stoi(argv[3]));
    pars.p_nongenetic = std::stod(argv[4]);
    pars.p_maternal = std::stod(argv[5]);
    pars.max_dev_time_step = static_cast<unsigned>(std::stoi(argv[6]));

    size_t argv_last{7};

    for (size_t s_idx{0}; s_idx < pars.L; ++s_idx)
    {
        pars.s[s_idx] = std::stod(argv[argv_last + s_idx]);
    }

    pars.sprime = std::stod(argv[argv_last + pars.L]);
    pars.data_output_interval = static_cast<unsigned>(std::stoi(argv[argv_last + pars.L + 1]));
    pars.file_name = argv[argv_last + pars.L + 2];
    pars.file_name_individuals = argv[argv_last + pars.L + 3];

    GRN_MatPat Simulation(pars);

    return 0;
}
