#include <vector>
#include <cassert>
#include <iostream>
#include "individual.hpp"
#include "gene_network_matpat.hpp"

//TODO:offspring size, have size affect the selective optimum
//then we have three channels: genes, RNA, size
//remaining channel: epigenetics

GRN_MatPat::GRN_MatPat(Parameters const &par) :
    rd{}, // random device to get a random seed
    seed{rd()}, // get the seed
    rng_r{seed}, // initialise the random-number generator
    par{par}, // initialize parameters
    data_file{par.file_name}, // data for averages
    data_file_individuals{par.file_name_individuals}, // data for individuals
    males(par.N/2, Individual(par, false)), // initialize the males
    females(par.N/2, Individual(par, true)), // initialize females
    meanW(par.L, std::vector < double >(par.L)), // the matrix with mean values for each matrix element (for canalization tests)
    meanS(2, std::vector<double>(par.L,0.0)), // vector with mean values of gene expression at the end of development
    C(2,std::vector < double >(par.L)), // the vector with the percentage of mutations in network without effect on gene expression
    Ce(2,std::vector < double >(par.L)) // the vector with the percentage of networks that are unchanged despite envtal perturbation 
{
    run();
}


// run the actual simulation
void GRN_MatPat::run()
{
    initialize_population();
    write_data_headers();

    write_out_all_individuals_headers();
    for (time_step = 0; 
            time_step <= par.max_time_step; ++time_step)
    {
        develop();

        // final timestep: run the statistics
        if (time_step == par.max_time_step)
        {
            write_out_all_individuals();
            genetic_canalization();
            environmental_canalization();
            time_to_stability();
        }

        // write out the data every nth generation
        if (time_step % par.data_output_interval == 0)
        {
            write_data();
        }

        reproduce();

    } // end for

    write_parameters();
} // end run()

// give each individual a random W matrix
void GRN_MatPat::initialize_population()
{
    for (auto male_iterator{males.begin()}; 
            male_iterator != males.end();
            ++male_iterator)
    {
        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                
                if (row_idx == par.sex_specific_locus_idx
                        && col_idx != row_idx)
                {
                    male_iterator->W[row_idx][col_idx] = 0.0;
                    continue;
                }
                
                male_iterator->W[row_idx][col_idx] = 
                    par.sd_init_strength_w * normal(rng_r);
            }
        }
    }
    
    for (auto female_iterator{females.begin()}; 
            female_iterator != females.end();
            ++female_iterator)
    {
        // initialize for each individual the elements of w
        // as per p689 2nd column 1st paragraph in Odorico et al
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                if (row_idx == par.sex_specific_locus_idx
                        && col_idx != row_idx)
                {
                    female_iterator->W[row_idx][col_idx] = 0.0;
                    continue;
                }

                female_iterator->W[row_idx][col_idx] = 
                    par.sd_init_strength_w * normal(rng_r);
            }
        }
    }
} // end initialize_population()


// develop all the individuals from birth
// until their gene expression levels equilibriate 
// or until a maximum number of time steps is reached
void GRN_MatPat::develop()
{
    for (auto male_iterator{males.begin()};
            male_iterator != males.end();
            ++male_iterator)
    {
        // update phenotypes during each developmental time step
        for (unsigned time_idx{0}; 
                time_idx < par.max_dev_time_step; 
                ++time_idx)
        {
            male_iterator->update_phenotype(time_idx);
        }

	male_iterator->average_phenotype();
    } 

    
    for (auto female_iterator{females.begin()};
            female_iterator != females.end();
            ++female_iterator)
    {
        for (unsigned time_idx{0}; time_idx < par.max_dev_time_step; ++time_idx)
        {
            female_iterator->update_phenotype(time_idx);
        }
	female_iterator->average_phenotype();
    } 
} // end GRN_MatPat::develop


void GRN_MatPat::reproduce()
{
    assert(females.size() > 0);
    assert(males.size() > 0);

    // clear out any old juveniles
    juveniles.clear();

    std::vector <double> male_fitnesses{};

    // make fitness distribution of males
    for (auto ind_iter{males.begin()}; 
            ind_iter != males.end(); 
            ++ind_iter)
    {
        male_fitnesses.push_back(ind_iter->fitness());
    }

    assert(male_fitnesses.size() >= 1);

    // set up fitness distributions for males
    // so that those individuals with larger 
    // fitness values are more
    // likely to be sampled
    std::discrete_distribution<unsigned> male_sampler(
            male_fitnesses.begin(),
            male_fitnesses.end());

    
    std::vector <double> female_fitnesses{};

    // make fitness distribution of males
    for (auto ind_iter{females.begin()}; 
            ind_iter != females.end();
            ++ind_iter)
    {
        // if individual dies swap
        female_fitnesses.push_back(ind_iter->fitness());
    }

    assert(females.size() == female_fitnesses.size());
    assert(males.size() == male_fitnesses.size());

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
                    females[mother_idx],
                    males[father_idx],
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
        if (juveniles[adult_idx].is_female)
        {
            females.push_back(juveniles[adult_idx]);
        }
        else
        {
            males.push_back(juveniles[adult_idx]);
        }
    }

    juveniles.clear();
} // end reproduce()

// give the column headers to the data file
void GRN_MatPat::write_data_headers()
{
    data_file << "time;"
            << "time_stability" << ";"
            << "fitness" << ";"
            << "varfitness" << ";"
            << "distance_optimum" << ";"
            << "var_distance_optimum" << ";";

    std::string sex_specifier;

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        for (unsigned sex_idx{0}; sex_idx < 2; ++sex_idx)
        {
            sex_specifier = sex_idx == 0 ? "m" : "f";

            data_file << "sbar" << "_" << sex_specifier << (row_idx + 1) << ";"
                << "varsbar" << "_" << sex_specifier << (row_idx + 1) << ";"
                << "v" << "_" << sex_specifier << (row_idx + 1) << ";"
                << "varv" << "_" << sex_specifier << (row_idx + 1) << ";"
                << "C" << "_" << sex_specifier << (row_idx + 1) << ";"
                << "Ce" << "_" << sex_specifier << (row_idx + 1) << ";";
        }

        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            data_file << "w" << (row_idx + 1) << (col_idx + 1) << ";"
                << "varw" << (row_idx + 1) << (col_idx + 1) << ";";

        }
    }

    data_file << std::endl;
}// end write_data_headers()

void GRN_MatPat::write_data()
{
    // allocate a 1d matrix to store the population stats for S
    // for each of the two sexes and then for each phenotypic trait
    std::vector < std::vector < double > > meanSbar(2, std::vector < double > (par.L, 0.0)); 
    std::vector < std::vector < double > > varSbar(2, std::vector < double > (par.L, 0.0)); 
    
    // allocate a 1d matrix to store the population stats for V
    // the within-individual variance over time in S
    std::vector < std::vector < double > > meanV(2, std::vector < double > (par.L, 0.0)); 
    std::vector < std::vector < double > > varV(2, std::vector < double > (par.L, 0.0)); 

    // parameters for mean fitness and variance in fitness
    double mean_fitness{0.0};
    double var_fitness{0.0};
    
    // parameters for distance and variance in distance from optimum 
    double mean_distance{0.0};
    double var_distance{0.0};

    // set all elements of the 2d matrix to store 
    // the population stats for W to 0
    // we need to keep meanW global 
    // as we also use it to calculate 
    // canalization of the average W matrix
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

    // make variables for numbers of individuals
    // rather than having to call .size() umpteen times
    unsigned nf{static_cast<unsigned>(females.size())};
    unsigned nm{static_cast<unsigned>(males.size())};
    
    // first go over all males
    for (auto male_iterator{males.begin()}; male_iterator != males.end();
            ++male_iterator)
    {
        // calculate W and S statistics
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            x = male_iterator->Sbar[row_idx];
            meanSbar[false][row_idx] += x;
            varSbar[false][row_idx] += x*x;

            x = male_iterator->V[row_idx];
            meanV[false][row_idx] += x;
            varV[false][row_idx] += x*x;

            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                x = male_iterator->W[row_idx][col_idx];

                meanW[row_idx][col_idx] += x;
                varW[row_idx][col_idx] += x*x;
            }
        }

        x = male_iterator->fitness();
        mean_fitness += x;
        var_fitness += x*x;
        
        x = male_iterator->mean_distance_to_optimum();
        mean_distance += x;
        var_distance += x*x;
    } // end for male_iterator

    for (auto female_iterator{females.begin()}; female_iterator != females.end();
            ++female_iterator)
    {
        // calculate W and S statistics
        for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
        {
            x = female_iterator->Sbar[row_idx];

            meanSbar[true][row_idx] += x;
            varSbar[true][row_idx] += x*x;

            x = female_iterator->V[row_idx];

            meanV[true][row_idx] += x;
            varV[true][row_idx] += x*x;

            for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
            {
                x = female_iterator->W[row_idx][col_idx];

                meanW[row_idx][col_idx] += x;
                varW[row_idx][col_idx] += x*x;
            }
        }
        x = female_iterator->fitness();
        mean_fitness += x;
        var_fitness += x*x;

        x = female_iterator->mean_distance_to_optimum();
        mean_distance += x;
        var_distance += x * x;
    } // end for female_iterator

    // begin the actual output
    data_file << time_step << ";";
    mean_fitness /= nf + nm;
    var_fitness = var_fitness / (nf + nm) - mean_fitness * mean_fitness;
    
    mean_distance /= (nf + nm);
    var_distance = var_distance / (nf + nm) - mean_distance * mean_distance;

    data_file 
        << t_stability << ";"
        << mean_fitness << ";" 
        << var_fitness << ";"
        << mean_distance << ";"
        << var_distance << ";"; 

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        for (unsigned sex_idx{0}; sex_idx < 2; ++sex_idx)
        {
            meanSbar[sex_idx][row_idx] /= nf + nm;
            varSbar[sex_idx][row_idx] = varSbar[sex_idx][row_idx] / (nf + nm) - 
                meanSbar[sex_idx][row_idx] * meanSbar[sex_idx][row_idx];
            
            meanV[sex_idx][row_idx] /= nf + nm;
            varV[sex_idx][row_idx] = varV[sex_idx][row_idx] / (nf + nm) - 
                meanV[sex_idx][row_idx] * meanV[sex_idx][row_idx];

            data_file << meanSbar[sex_idx][row_idx] << ";"
                << varSbar[sex_idx][row_idx] << ";"
                << meanV[sex_idx][row_idx] << ";"
                << varV[sex_idx][row_idx] << ";"
                << C[sex_idx][row_idx] << ";"
                << Ce[sex_idx][row_idx] << ";";
        }

        for (unsigned col_idx{0}; col_idx < par.L; ++col_idx)
        {
            meanW[row_idx][col_idx] /= nf + nm;
            varW[row_idx][col_idx] = varW[row_idx][col_idx] / (nf + nm) -
                meanW[row_idx][col_idx] * meanW[row_idx][col_idx];

            data_file << meanW[row_idx][col_idx] << ";" 
                << varW[row_idx][col_idx] << ";";
        }
    }

    meanS = meanSbar;

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
// if one wants to obtain information about the development of phenotypes
// then make sure to call this after development(), rather than
// after reproduce()
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
        << "N_genetic_canalization;" << par.N_genetic_canalization << std::endl
        << "N_environmental_canalization;" << par.N_environmental_canalization << std::endl
        << "L;" << par.L << std::endl
        << "sd_init_strength_w;" << par.sd_init_strength_w << std::endl
        << "a;" << par.a << std::endl
        << "max_time_step;" << par.max_time_step << std::endl
        << "max_dev_time_step;" << par.max_dev_time_step << std::endl
        << "max_dev_time_step_nstats;" << par.max_dev_time_step_stats << std::endl
        << "max_dev_time_step_envt_canalize;" << par.max_dev_time_step_envt_canalize << std::endl
        << "mu_w;" << par.mu_w << std::endl
        << "mu_w;" << par.sdmu_w << std::endl
        << "init_sk_f;" << par.init_sk[1] << std::endl
        << "init_sk_m;" << par.init_sk[0] << std::endl;

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file << "theta_f" << (row_idx + 1) << ";" << par.theta[true][row_idx] << std::endl
            << "theta_m" << (row_idx + 1) << ";" << par.theta[false][row_idx] << std::endl
            << "s" << (row_idx + 1) << ";" << par.s[row_idx] << std::endl
            << "sprime" << (row_idx + 1) << ";" << par.sprime[row_idx] << std::endl
            << "p_nongenetic" << (row_idx + 1) << ";" << par.p_nongenetic[row_idx] << std::endl
            << "p_maternal" << (row_idx + 1) << ";" << par.p_maternal[row_idx] << std::endl;
    }

} // end write_parameters

void GRN_MatPat::make_reference_individual(Individual &ref)
{
    // first assign this individual the average W
    ref.W = meanW;

    double initial_expression;

    bool is_female = ref.is_female;
    
    // then initilize gene loci
    for (unsigned int row_idx{0}; row_idx < par.L; ++row_idx)
    {
        // if the locus is the sex-specific one, then set
        // the initial level of expression accordingly,
        // otherwise just set it to be equal to the constitutive 
        // gene expression
        initial_expression = 
            row_idx == par.sex_specific_locus_idx ? par.init_sk[is_female] : par.a;

        // now perform development given an average S0
        //
        // because the mean phenotype is the same for males and females
        // no need to focus on paternal vs maternal effects. However, 
        // once we study sexual dimorphism, we can fully implement this
        ref.S[0][row_idx] = (1.0 - par.p_nongenetic[row_idx]) * initial_expression
            + par.p_nongenetic[row_idx] * par.p_maternal[row_idx] * meanS[true][row_idx]
            + par.p_nongenetic[row_idx] * (1.0 - par.p_maternal[row_idx]) * meanS[false][row_idx];
    } // end for unsigned row_idx
        
    // development, yay
    for (unsigned dev_time_step_idx{0}; 
            dev_time_step_idx < par.max_dev_time_step;
            ++dev_time_step_idx)
    {
        ref.update_phenotype(dev_time_step_idx);
    }
} // end make_reference_individual()

// calculate genetic canalization by mutating individuals
// and look at their phenotype development, see p690 1st col, 3rd para
// in Odorico et al
void GRN_MatPat::genetic_canalization()
{
    // each network has par.L^2 number of nodes of which we 
    // need only one to mutate. Here we build a uniform sample distribution
    // between 0 and (par.L^2 - 1) to sample the node that will be mutated
    // later on we will then calculate a row and a column index for this number
    // (see below)
    std::uniform_int_distribution<unsigned> node_sampler{0, par.L * par.L - 1};

    // vector for male and female reference individuals
    // against which clones can be compared
    std::vector < Individual > reference_individuals{};

    // male on position 1
    reference_individuals.push_back(Individual(par, false));
    make_reference_individual(reference_individuals[0]);
    // female on position 2
    reference_individuals.push_back(Individual(par, true));
    make_reference_individual(reference_individuals[1]);

    // some error checking, to ensure the mean values of W
    // are indeed properly copied over to the reference individual
    assert(reference_individuals[0].W[par.L - 1][1] == meanW[par.L - 1][1]);
    assert(reference_individuals[0].W[0][par.L - 1] == meanW[0][par.L - 1]);
    assert(reference_individuals[1].W[par.L - 1][1] == meanW[par.L - 1][1]);
    assert(reference_individuals[1].W[0][par.L - 1] == meanW[0][par.L - 1]);

    // reset C to 0.0
    C = std::vector < std::vector < double > >(2, std::vector <double> (par.L,0.0));

    double expression_difference_i, initial_expression;

    unsigned node_sequential_id, 
             col_idx_to_mutate, 
             row_idx_to_mutate;

    for (unsigned sex_idx{0}; sex_idx < 2; ++sex_idx)
    {
        // now generate 1000 clones with 1 mutation. I take it mutation works 
        // by simply changing one of the wij's by mutation
        for (unsigned int individual_idx{0}; 
                individual_idx < par.N_genetic_canalization; ++individual_idx)
        {
            // make a clone
            Individual clone(par, static_cast<bool>(sex_idx));

            do 
            {
                // randomly sample a number between 0 and par.L^2, which is
                // sampling the id of the w_ij to mutate
                node_sequential_id = node_sampler(rng_r);

                // find the column and the row that belongs to this value
                // use modulo operator to get column, as: 
                // number               0 1 2 3 4 5 6 7 ... par.L^2 - 1 
                // number % par.L       0 1 2 3 4 5 0 1 ... par.L - 1
                //
                // we use floor of number / L to get the row
                // number               0 1 2 3 4 5 6 7 ... par.L - 1 
                // floor(number/L)      0 0 0 0 0 1 1 1 ... par.L - 1
                col_idx_to_mutate = node_sequential_id % par.L;

                row_idx_to_mutate = 
                    static_cast<unsigned>(
                                std::floor(
                                    static_cast<double>(
                                        node_sequential_id) / par.L
                                    )
                                );

                // sample another locus if this is part of the 
                // sex specific hierarchy
            } while(row_idx_to_mutate == par.sex_specific_locus_idx
                    && row_idx_to_mutate != col_idx_to_mutate);

            // do bounds checking on the col and row
            // of the element that is selected for mutation
            assert(col_idx_to_mutate >= 0);
            assert(col_idx_to_mutate < par.L);

            assert(row_idx_to_mutate >= 0);
            assert(row_idx_to_mutate < par.L);

            // assign all the weights
            clone.W = meanW;

            // first assign this individual the average W
            for (unsigned int row_idx{0}; row_idx < par.L; ++row_idx)
            {
                initial_expression = 
                    row_idx == par.sex_specific_locus_idx ? 
                        par.init_sk[sex_idx] 
                        : 
                        par.a;
                // now perform development given an average S0
                //
                // because the mean phenotype is the same for males and females
                // no need to focus on paternal vs maternal effects. However, 
                // once we study sexual dimorphism, we can fully implement this
                clone.S[0][row_idx] = 
                    (1.0 - par.p_nongenetic[row_idx]) * initial_expression
                    + par.p_nongenetic[row_idx] * par.p_maternal[row_idx] * meanS[true][row_idx]
                    + par.p_nongenetic[row_idx] * (1.0 - par.p_maternal[row_idx]) * meanS[false][row_idx];
            } // end unsigned int row_idx
           
            // then change one element by mutation
            clone.W[row_idx_to_mutate][col_idx_to_mutate] += normal(rng_r) * par.sdmu_w;

            // have the clone develop over the developmental time steps
            for (unsigned dev_time_step_idx{0}; 
                    dev_time_step_idx < par.max_dev_time_step;
                    ++dev_time_step_idx)
            {
                clone.update_phenotype(dev_time_step_idx);
            }

            // now look at gene loci at the final dev time step
            // and see how it differs from the reference individual,
            // that is the value of C_i as per page 690 first col
            // 3rd paragraph in Odorico et al
            for (unsigned int locus_idx{0}; locus_idx < par.L; ++locus_idx)
            {
                // calculate expression difference between current and reference
                // individual
                expression_difference_i = std::fabs(clone.S[par.max_dev_time_step - 1][locus_idx] - 
                    reference_individuals[sex_idx].S[par.max_dev_time_step - 1][locus_idx]);

                // only calculate those individuals whose gene expression is
                // unaltered despite mutation
                if (expression_difference_i < 
                        par.canalization_threshold)
                {
                    C[sex_idx][locus_idx] += 1.0;
                }
            } // end for unsigned locus_idx
        } // end for unsigned individual_idx

        // now finallly transform C (which is counts up to now)
        // to contain percentages
        for (unsigned int locus_idx{0}; 
                locus_idx < par.L; ++locus_idx)
        {
            C[sex_idx][locus_idx] /= par.N_genetic_canalization;
        }
    }// end unsigned sex_idx
} // end mean_canalization 

// obtain statistics on the amount of environmental
// canalization, as per Odorico et al p.690, 1st col,
// 4th paragraph
void GRN_MatPat::environmental_canalization()
{
    // reset Ce vector to 0
    Ce = std::vector < std::vector < double > >(
            2, 
            std::vector < double >(par.L,0.0));

    for (unsigned int sex_idx{0}; sex_idx < 2; ++sex_idx)
    {
        // sample individuals and give
        // then give them random initial expression of gene loci
        for (unsigned individual_idx{0};
                individual_idx < par.N_environmental_canalization;
                ++individual_idx)
        {
            Individual clone(par, static_cast<bool>(sex_idx));
            clone.W = meanW;

            // assign individuals a different S vector so that it can
            // deal with different number of time steps
            clone.S = std::vector < 
                std::vector < double > >
                        (par.max_dev_time_step_envt_canalize, std::vector < double > (par.L, 0.0));

            // sample loci from a uniform distribution
            for (unsigned int locus_idx{0};
                    locus_idx < par.L;
                    ++locus_idx)
            {
                if (locus_idx == par.sex_specific_locus_idx)
                {
                    clone.S[0][locus_idx] = par.init_sk[sex_idx];
                } 
                else
                {
                    clone.S[0][locus_idx] = uniform(rng_r);
                }
            } // end for locus_idx

            assert(clone.S.size() == par.max_dev_time_step_envt_canalize);

            for (unsigned int dev_time_step_idx{0}; 
                    dev_time_step_idx < par.max_dev_time_step_envt_canalize;
                    ++dev_time_step_idx)
            {
                clone.update_phenotype(dev_time_step_idx);
    //            std::cout << "clone " << individual_idx << " ";

    //            for (unsigned int locus_idx{0};
    //                    locus_idx < par.L;
    //                    ++locus_idx)
    //            {
    //                std::cout << locus_idx << " " << clone.S[dev_time_step_idx][locus_idx] << " ";
    //            }
    //
    //            std::cout << std::endl;
            } // end for unsigned int dev_time_step_idx

            // now calculate difference with reference individual
            for (unsigned int locus_idx{0};
                    locus_idx < par.L;
                    ++locus_idx)
            {
    //            std::cout << locus_idx << ": " 
    //                << clone.S[0][locus_idx] << " " 
    //                << clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] << " " 
    //                << meanS[true][locus_idx] << " "
    //                << std::fabs(clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] -
    //                meanS[locus_idx]) << " ";

                if (std::fabs(clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] -
                    meanS[clone.is_female][locus_idx]) < par.stability_threshold)
                {
                    Ce[clone.is_female][locus_idx] += 1.0;
                }
            } // end for locus_idx

            std::cout << std::endl;
        } // end for individual_idx
        
        for (unsigned int locus_idx{0};
                locus_idx < par.L;
                ++locus_idx)
        {
            Ce[sex_idx][locus_idx] /= par.N_environmental_canalization;
        }
    }// end for sex_idx
}// end environmental_canalization()

// assess time to equilibrium for phenotypes of
// the average individual in the population
void GRN_MatPat::time_to_stability()
{
    Individual reference_individual(par, true);
    reference_individual.W = meanW;

    // reinitialize S, the per-locus gene expression 
    // for each developmental time step
    // as it now needs to deal with a different number of
    // time steps to assess stability
    reference_individual.S = std::vector < std::vector < double > >(
            par.max_dev_time_step_stability, std::vector < double > (par.L, 0.0)
            );

    for (unsigned int locus_idx{0}; 
            locus_idx < par.L; ++locus_idx)
    {
        reference_individual.S[0][locus_idx] = par.a;
    }

    bool converged;

    reference_individual.update_phenotype(0);

    unsigned int time_idx;
    for (time_idx = 1;
            time_idx < par.max_dev_time_step_stability;
            ++time_idx)
    {
        converged = true;
        assert(time_idx < reference_individual.S.size());
        reference_individual.update_phenotype(time_idx);

        for (unsigned int locus_idx{0};
                locus_idx < par.L;
                ++locus_idx)
        {
            if (std::fabs(reference_individual.S[time_idx][locus_idx] - 
                        reference_individual.S[time_idx - 1][locus_idx]) > 
                            par.stability_threshold)
            {
                converged = false;
                break;
            }
        }

        if (converged)
        {
            break;
        }
    }

    t_stability = time_idx;
} // end time_to_stability()
