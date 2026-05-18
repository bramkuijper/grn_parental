#include <vector>
#include <cassert>
#include <iostream>
#include "individual.hpp"
#include "gene_network_matpat.hpp"

// TODO: should parental effects also influence
// sex-specific regulatory locus, at the moment they do but 
// perhaps it should not.
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
    meanS(par.L, 0.0), // vector with mean values of gene expression at the end of development
    C(par.L, 0.0), // the vector with the percentage of mutations in network without effect on gene expression
    Ce(par.L, 0.0) // the vector with the percentage of networks that are unchanged despite envtal perturbation 
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
//        write_out_all_individuals();

        // final timestep: run the statistics
        if (time_step == par.max_time_step)
        {
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

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file << "sbar" << (row_idx + 1) << ";"
            << "varsbar" << (row_idx + 1) << ";"
            << "v" << (row_idx + 1) << ";"
            << "varv" << (row_idx + 1) << ";"
            << "C" << (row_idx + 1) << ";"
            << "Ce" << (row_idx + 1) << ";";

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
    std::vector < double > meanSbar(par.L, 0.0); 
    std::vector < double > varSbar(par.L, 0.0); 
    
    // allocate a 1d matrix to store the population stats for V
    // the within-individual variance over time in S
    std::vector < double > meanV(par.L, 0.0); 
    std::vector < double > varV(par.L, 0.0); 

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
        meanSbar[row_idx] /= nf + nm;
        varSbar[row_idx] = varSbar[row_idx] / (nf + nm) - 
            meanSbar[row_idx] * meanSbar[row_idx];
        
        meanV[row_idx] /= nf + nm;
        varV[row_idx] = varV[row_idx] / (nf + nm) - 
            meanV[row_idx] * meanV[row_idx];

        data_file << meanSbar[row_idx] << ";"
            << varSbar[row_idx] << ";"
            << meanV[row_idx] << ";"
            << varV[row_idx] << ";"
            << C[row_idx] << ";"
            << Ce[row_idx] << ";";

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
        << "s_init;" << par.s_init << std::endl
        << "sd_init_strength_w;" << par.sd_init_strength_w << std::endl
        << "a;" << par.a << std::endl
        << "p_nongenetic;" << par.p_nongenetic << std::endl
        << "p_maternal;" << par.p_maternal << std::endl
        << "max_time_step;" << par.max_time_step << std::endl
        << "max_dev_time_step;" << par.max_dev_time_step << std::endl
        << "max_dev_time_step_nstats;" << par.max_dev_time_step_stats << std::endl
        << "max_dev_time_step_envt_canalize;" << par.max_dev_time_step_envt_canalize << std::endl
        << "mu_w;" << par.mu_w << std::endl
        << "mu_w;" << par.sdmu_w << std::endl;

    for (unsigned row_idx{0}; row_idx < par.L; ++row_idx)
    {
        data_file << "theta_f" << (row_idx + 1) << ";" << par.theta_f[row_idx] << std::endl;
        data_file << "theta_m" << (row_idx + 1) << ";" << par.theta_m[row_idx] << std::endl;
        data_file << "s" << (row_idx + 1) << ";" << par.s[row_idx] << std::endl;
        data_file << "sprime" << (row_idx + 1) << ";" << par.sprime[row_idx] << std::endl;
    }

} // end write_parameters

void GRN_MatPat::make_reference_individual(Individual &ref)
{
    // first assign this individual the average W
    ref.W = meanW;
    
    // then initilize gene loci
    for (unsigned int row_idx{0}; row_idx < par.L; ++row_idx)
    {
        // now perform development given an average S0
        //
        // because the mean phenotype is the same for males and females
        // no need to focus on paternal vs maternal effects. However, 
        // once we study sexual dimorphism, we can fully implement this
        ref.S[0][row_idx] = (1.0 - par.p_nongenetic) * par.a
            + par.p_nongenetic * meanS[row_idx];
    }
        
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

    // first make a reference individual
    // against which all the clones will be compared
    Individual reference_individual(par);
    make_reference_individual(reference_individual);
   
    // some error checking, to ensure the mean values of W
    // are indeed properly copied over to the reference individual
    assert(reference_individual.W[par.L - 1][1] == meanW[par.L - 1][1]);
    assert(reference_individual.W[0][par.L - 1] == meanW[0][par.L - 1]);

    // reset C to 0.0
    fill(C.begin(), C.end(), 0.0);

    double expression_difference_i;

    // now generate 1000 clones with 1 mutation. I take it mutation works 
    // by simply changing one of the wij's by mutation
    for (unsigned int individual_idx{0}; 
            individual_idx < par.N_genetic_canalization; ++individual_idx)
    {
        // make a clone
        Individual clone(par);

        // randomly sample a number between 0 and par.L^2, which is
        // sampling the id of the w_ij to mutate
        unsigned node_sequential_id{node_sampler(rng_r)};

        // find the column and the row that belongs to this value
        // use modulo operator to get column, as: 
        // number               0 1 2 3 4 5 6 7 ... par.L^2 - 1 
        // number % par.L       0 1 2 3 4 5 0 1 ... par.L - 1
        //
        // we use floor of number / L to get the row
        // number               0 1 2 3 4 5 6 7 ... par.L - 1 
        // floor(number/L)      0 0 0 0 0 1 1 1 ... par.L - 1
        unsigned col_idx_to_mutate{node_sequential_id % par.L};

        unsigned row_idx_to_mutate{
            static_cast<unsigned>(
                        std::floor(
                            static_cast<double>(
                                node_sequential_id) / par.L
                            )
                        )
                    };

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
            // now perform development given an average S0
            //
            // because the mean phenotype is the same for males and females
            // no need to focus on paternal vs maternal effects. However, 
            // once we study sexual dimorphism, we can fully implement this
            clone.S[0][row_idx] = (1.0 - par.p_nongenetic) * par.a
                + par.p_nongenetic * meanS[row_idx];
        }
       
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
                reference_individual.S[par.max_dev_time_step - 1][locus_idx]);

            // only calculate those individuals whose gene expression is
            // unaltered despite mutation
            if (expression_difference_i < 
                    par.canalization_threshold)
            {
                C[locus_idx] += 1.0;
            }
        } // end for unsigned locus_idx
    } // end for unsigned individual_idx

    // now finallly transform C (which is counts up to now)
    // to contain percentages
    for (unsigned int locus_idx{0}; 
            locus_idx < par.L; ++locus_idx)
    {
        C[locus_idx] /= par.N_genetic_canalization;
    }
} // end mean_canalization 

// obtain statistics on the amount of environmental
// canalization, as per Odorico et al p.690, 1st col,
// 4th paragraph
void GRN_MatPat::environmental_canalization()
{
    // reset the Ce vector
    fill(Ce.begin(), Ce.end(), 0.0);
    // sample individuals and give
    // then give them random initial expression of gene loci
    for (unsigned individual_idx{0};
            individual_idx < par.N_environmental_canalization;
            ++individual_idx)
    {
        Individual clone(par);
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
            clone.S[0][locus_idx] = uniform(rng_r);
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
            std::cout << locus_idx << ": " 
                << clone.S[0][locus_idx] << " " 
                << clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] << " " 
                << meanS[locus_idx] << " "
                << std::fabs(clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] -
                meanS[locus_idx]) << " ";

            if (std::fabs(clone.S[par.max_dev_time_step_envt_canalize - 1][locus_idx] -
                meanS[locus_idx]) < par.stability_threshold)
            {
                Ce[locus_idx] += 1.0;
            }
        } // end for locus_idx

        std::cout << std::endl;
    } // end for individual_idx
    
    for (unsigned int locus_idx{0};
            locus_idx < par.L;
            ++locus_idx)
    {
        Ce[locus_idx] /= par.N_environmental_canalization;
    }
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
