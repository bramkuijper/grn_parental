#include <iostream>
#include "parameters.hpp"
#include "gene_network_matpat.hpp"

int main(int argc, char **argv)
{
    int Lval = std::stoi(argv[2]);
    Parameters pars(Lval);

    pars.max_time_step = std::stoi(argv[1]);
    pars.p_nongenetic = std::stod(argv[3]);
    pars.p_maternal = std::stod(argv[4]);
    pars.max_dev_time_step = std::stoi(argv[5]);

    int argv_last = 6;

    for (int s_idx{0}; s_idx < pars.L; ++s_idx)
    {
        pars.s[s_idx] = std::stod(argv[argv_last + s_idx]);
    }

    pars.file_name = argv[argv_last + pars.L];
    pars.file_name_individuals = argv[argv_last + pars.L + 1];

    GRN_MatPat Simulation(pars);

    return 0;
}
