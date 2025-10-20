#include "gene_network_matpat.hpp"

GRN_MatPat::GRN_MatPat(Parameters const &par) :
    rd{},
    seed{rd()},
    rng_r{seed},
    data_file{par.file_name},
    males(par.N/2, Individual(par)),
    females(par.N/2, Individual(par))
{}

// run the actual simulation
void GRN_MatPat::run()
{
    for (time_step = 0; 
            time_step < par.max_time_step; ++time_step)
    {
        survive();
        reproduce();

        // write out the data every nth generation
        if (time_step % par.data_output_interval == 0)
        {
            write_data();
        }
    } // end for
} // end run()

void GRN_MatPat::survive()
{


} // end survive()
