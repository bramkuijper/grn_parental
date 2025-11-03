#include "parameters.hpp"
#include "gene_network_matpat.hpp"

int main(int argc, char **argv)
{
    Parameters pars;

    pars.max_time_step = std::stoi(argv[1]);
    pars.L = std::stoi(argv[2]);
    pars.p_nongenetic = std::stod(argv[3]);
    pars.p_maternal = std::stod(argv[4]);
    pars.max_dev_time_step = std::stoi(argv[5]);
    pars.file_name = argv[6];
    pars.file_name_individuals = argv[7];


    GRN_MatPat Simulation(pars);

    return 0;
}
