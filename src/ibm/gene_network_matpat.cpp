#include <vector>
#include <cassert>
#include "individual.hpp"
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
        reproduce();

        // write out the data every nth generation
        if (time_step % par.data_output_interval == 0)
        {
            write_data();
        }
    } // end for
} // end run()



void GRN_MatPat::reproduce()
{
    // clear out any old juveniles
    juveniles.clear();

    std::vector <double> male_fitnesses{};
    // make fitness distribution of males
    for (auto ind_iter{males.begin()}; ind_iter != males.end(); )
    {
        male_fitnesses.push_back(ind_iter->fitness());
    }

    // set up fitness distributions for males
    // so that those individuals with larger 
    // fitness values are more
    // likely to be sampled
    std::discrete_distribution<unsigned> male_sampler(
            male_fitnesses.begin(),
            male_fitnesses.end());

    
    std::vector <double> female_fitnesses{};
    // make fitness distribution of males
    for (auto ind_iter{females.begin()}; ind_iter != females.end(); )
    {
        // if individual dies swap
        female_fitnesses.push_back(ind_iter->fitness());
    }

    // set up fitness distributions for females
    // so that those individuals with larger 
    // fitness values are more
    // likely to be sampled
    std::discrete_distribution<unsigned> female_sampler(
            female_fitnesses.begin(),
            female_fitnesses.end());

    // aux variables for sampling parents
    unsigned father_idx, mother_idx;

    for (unsigned offspring_idx{0}; 
            offspring_idx < par.N; ++offspring_idx)
    {
        father_idx = male_sampler(rng_r);
        mother_idx = female_sampler(rng_r);

        assert(father_idx < males.size());
        assert(mother_idx < females.size());

        juveniles.push_back(
                Individual(             // use birth constructor 
                    males[father_idx],
                    females[mother_idx],
                    par,
                    rng_r
                    )
                );
    }

    males.clear();
    females.clear();

    for (unsigned adult_idx{0};
            adult_idx < par.N; ++adult_idx)
    {
        if (uniform(rng_r) < 0.5)
        {
            males.push_back(juveniles[adult_idx]);
        }
        else
        {
            females.push_back(juveniles[adult_idx]);
        }
    }
} // end reproduce()

void GRN_MatPat::write_data()
{

} // end write_data()
