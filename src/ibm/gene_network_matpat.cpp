#include <vector>
#include <cassert>
#include "individual.hpp"
#include "gene_network_matpat.hpp"

GRN_MatPat::GRN_MatPat(Parameters const &par) :
    rd{},
    seed{rd()},
    rng_r{seed},
    data_file{par.file_name},
    data_file_individuals{par.file_name_individuals},
    males(par.N/2, Individual(par)),
    females(par.N/2, Individual(par)),
    meanW(par.L, std::vector < double >(par.L))
{}

// give each individual a random W matrix
void GRN_MatPat::initialize_population()
{
    for (auto male_itaterator{males.begin()}; 
            male_iterator != males.end();
            ++males)
    {
        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                male_iterator->w[row_idx][col_idx] = 
                    par.sd_init_strength_w * normal(rng_r);
            }
        }
    }
    
    for (auto female_itaterator{females.begin()}; 
            female_iterator != females.end();
            ++females)
    {
        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                female_iterator->w[row_idx][col_idx] = 
                    par.sd_init_strength_w * normal(rng_r);
            }
        }
    }
} // end initialize_population()

// run the actual simulation
void GRN_MatPat::run()
{
    initialize_population();

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

    juveniles.clear();
} // end reproduce()

// give the column headers to the data file
void GRN_MatPat::write_data_headers()
{
    data_file << "time;";

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file << "sbar" << (row_idx + 1) << ";"
            << "varsbar" << (row_idx + 1) << ";"
            << "v" << (row_idx + 1) << ";"
            << "varv" << (row_idx + 1) << ";";

        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            data_file << "w" << (row_idx + 1) << (col_idx + 1) << ";"
                << "varw" << (row_idx + 1) << (col_idx + 1) << ";";

        }
    }
}// end write_data_headers()

void GRN_MatPat::write_data()
{
    // allocate a 1d matrix to store the population stats for S
    std::vector < double > meanSbar(par.L, 0.0); 
    std::vector < double > varSbar(par.L, 0.0); 
    
    // allocate a 1d matrix to store the population stats for V
    // the within-individual variance over time in S
    std::vector < double > meanV(par.L, 0.0); 
    std::vector < double > varV(par.L, 0.0); 

    // set all elements of the 2d matrix to store 
    // the population stats for W to 0
    // we keep meanW global as we also use it to calculate 
    // canalization
    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            meanW[row_idx][col_idx] = 0.0;
        }
    }

    // initialize a vector for the variance in W with 0s
    std::vector < std::vector < double > > 
        varW(par.L, std::vector< double >(par.L, 0.0));

    double x;  // aux variable to temporarily store each trait value
    
    // first go over all males
    for (auto male_iterator{males.begin()}; male_iterator != males.end();
            ++male_iterator)
    {
        // calculate W and S statistics
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            x = male_iterator->Sbar[row_idx];

            meanSbar[row_idx] += x;
            varSbar[row_idx] += x*x;

            x = male_iterator->V[row_idx];

            meanV[row_idx] += x;
            varV[row_idx] += x*x;

            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                x = male_iterator->W[row_idx][col_idx];

                meanW[row_idx][col_idx] += x;
                varW[row_idx][col_idx] += x*x;
            }
        }
    } // end for male_iterator
    
    for (auto female_iterator{females.begin()}; female_iterator != females.end();
            ++female_iterator)
    {
        // calculate W and S statistics
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            x = female_iterator->Sbar[row_idx];

            meanSbar[row_idx] += x;
            varSbar[row_idx] += x*x;

            x = female_iterator->V[row_idx];

            meanV[row_idx] += x;
            varV[row_idx] += x*x;

            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                x = female_iterator->W[row_idx][col_idx];

                meanW[row_idx][col_idx] += x;
                varW[row_idx][col_idx] += x*x;
            }
        }
    } // end for female_iterator

    // make variables for numbers of individuals
    // rather than having to call .size() umpteen times
    unsigned nf{static_cast<unsigned>(females.size())};
    unsigned nm{static_cast<unsigned>(males.size())};

    // begin the actual output
    data_file << time_step << ";" << std::endl;

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        meanSbar[row_idx] /= nf + nm;
        varSbar[row_idx] = varSbar[row_idx] / (nf + nm) - 
            meanSbar[row_idx] * meanSbar[row_idx];
        
        meanV[row_idx] /= nf + nm;
        varV[row_idx] = varV[row_idx] / (nf + nm) - 
            meanV[row_idx] * meanV[row_idx];

        data_file << meanSbar[row_idx] << ";"
            << varSbar[row_idx] << ";"
            << meanV[row_idx] << ";"
            << varV[row_idx] << ";";

        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            meanW[row_idx][col_idx] /= nf + nm;
            varW[row_idx][col_idx] = varW[row_idx][col_idx] / (nf + nm) -
                meanW[row_idx][col_idx] * meanW[row_idx][col_idx];

            data_file << meanW[row_idx][col_idx] << ";" 
                << varW[row_idx][col_idx] << ";";
        }
    }

    // finish  by starting a new line
    data_file << std::endl;

} // end write_data()

// headers to the individual data file
void GRN_MatPat::write_out_all_individuals_headers()
{
    data_file_individuals << "time;id;sex;";
        
    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            data_file_individuals << "W" << (row_idx + 1) << (col_idx + 1) << ";";
        }
    }

    for (unsigned t_idx{0}; t_idx < par.max_dev_time_step; ++t_idx)
    {
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            data_file_individuals << "S" << (t_idx + 1) << (row_idx + 1) << ";";
        }
    }

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file_individuals << "Sbar" << (row_idx + 1) << ";";
    }

    data_file_individuals << std::endl;
}

// for inspection purposes, write out all individuals
void GRN_MatPat::write_out_all_individuals()
{
    unsigned individual_idx{0};

    for (auto male_iterator{males.begin()}; male_iterator != males.end();
            ++male_iterator)
    {
        data_file_individuals << time_step << ";" << individual_idx << ";male;";
        ++individual_idx;

        // output W
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                data_file_individuals << male_iterator->W[row_idx][col_idx] << ";";
            }
        }

        // now output S
        for (unsigned t_idx{0}; t_idx < par.max_dev_time_step; ++t_idx)
        {
            for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
            {
                data_file_individuals << male_iterator->S[t_idx][row_idx] << ";";
            }
        }

        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            data_file_individuals << male_iterator->Sbar[row_idx] << ";";
        }

        data_file_individuals << std::endl;
    } // end for male_iterator
    
    for (auto female_iterator{females.begin()}; female_iterator != females.end();
            ++female_iterator)
    {
        data_file_individuals << time_step << ";" << individual_idx << ";female;";
        ++individual_idx;

        // output W
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                data_file_individuals << female_iterator->W[row_idx][col_idx] << ";";
            }
        }

        // now output S
        for (unsigned t_idx{0}; t_idx < par.max_dev_time_step; ++t_idx)
        {
            for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
            {
                data_file_individuals << female_iterator->S[t_idx][row_idx] << ";";
            }
        }
        
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            data_file_individuals << female_iterator->Sbar[row_idx] << ";";
        }

        data_file_individuals << std::endl;
    } // end for female_iterator
}// end write_out_all_individuals

void GRN_MatPat::write_parameters()
{
    data_file << std::endl
        << std::endl
        << "N;" << par.N << std::endl
        << "L;" << par.L << std::endl
        << "w_init;" << par.w_init << std::endl
        << "s_init;" << par.s_init << std::endl
        << "sd_init_strength_w;" << par.sd_init_strength_w << std::endl
        << "a;" << par.a << std::endl
        << "p_nongenetic;" << par.p_nongenetic << std::endl
        << "p_maternal;" << par.p_maternal << std::endl
        << "max_time_step;" << par.max_time_step << std::endl
        << "max_dev_time_step;" << par.max_dev_time_step << std::endl
        << "max_dev_time_step_nstats;" << par.max_dev_time_step_nstats << std::endl
        << "mu_w;" << par.mu_w << std::endl
        << "mu_w;" << par.sdmu_w << std::endl;
        << "sprime;" << par.sprime << std::endl;

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file << "theta" << (row_idx + 1) << ";" << theta[row_idx] << std::endl;
        data_file << "s" << (row_idx + 1) << ";" << s[row_idx] << std::endl;
    }

} // end write_parameters

// calculate canalization as 
void GRN_MatPat::mean_genetic_canalization()
{
    // ok randomly select one of the par.L^2 values of w_ij for 
    // mutation
    std::uniform_int_distribution<unsigned> node_sampler{0, par.L * par.L - 1};

    // generate 1000 clones with 1 mutation. I take it mutation works 
    // by simply changing one of the wij's by mutation
    for (unsigned int individual_idx{0}; 
            individual_idx < par.N_canalize; ++individual_idx)
    {
        Individual clone(par);

        // randomly sample a number between 0 and par.L^2, which is
        // sampling the id of the w_ij to mutate
        unsigned node_sequential_id{node_sampler(rng_r)};

        // find the column and the row that belongs to this value
        // use modulo operator to get column, as: 
        // number               0 1 2 3 4 5 6 7 ... par.L - 1 
        // number % par.L       0 1 2 3 4 0 1 2 ... par.L - 1
        //
        // we use floor of number / L to get the row
        // number               0 1 2 3 4 5 6 7 ... par.L - 1 
        // floor(number/L)      0 0 0 0 0 1 1 1 ... par.L - 1
        unsigned col_idx_to_mutate{node_sequential_id % par.L};
        unsigned row_idx_to_mutate{std::floor(
                static_cast<double>(node_sequential_id) / 6)};

        // first assign this individual the average W
        for (unsigned int row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned int col_idx{0}; col_idx < par.L; ++col_idx)
            {
                clone.W[row_idx][col_idx] = meanW[row_idx][col_idx];

            }
        }
       
        // then change one element by mutation
        clone.W[row_idx_to_mutate][col_idx_to_mutate] += sdmu_w * normal(rng_r);

        // TODO: now remainder of canalization calculations
    }
} // end mean_canalization 
