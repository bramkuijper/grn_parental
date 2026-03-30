#!/usr/bin/env Rscript

# number of replicates
nrep <- 20

# maximum amount of generations a simulation runs
maxgen <-10000

# the number of phenotypes, L
# note that the GRN is size L*L
L <- c(6)

a <- 0.2

# fraction of phenotype that is nongenetically inherited
p_nongenetic <- c(0.0)

p_maternal <- c(0.0)

# number of time steps before development is completed
dev_time <- c(4)

# generation of output files
# generate a date_time stamp as a character
date_str <- format(Sys.time(), "%d%m%Y_%H%M%S")

# generate names for the output file based on date and time 
# to ensure some sort of uniqueness, to avoiding that output
# files from different runs overwrite eachother
output_file_prefix <- paste0("sim_grn_",date_str)

counter <- 0

exe = "./gene_network_matpat.exe"


# 2 out of L values are having stabilizing selection
# for which s > 0
# see Odorico et al p 690, 1st col, final para
s_values <- numeric(length=L[[1]])

# set the first and the last one to nonzero values
s_values[[1]] <- 0.5
s_values[[L[[1]]]] <- 0.5

batch_file_contents <- ""

for (rep_i in 1:nrep)
{
    for (L_i in L)
    {
        for (p_nongenetic_i in p_nongenetic)
        {
            for (p_maternal_i in p_maternal)
            {
                for (dev_time_i in dev_time)
                {
                    counter <- counter + 1
                    file_name_i <- paste0(
                            output_file_prefix,"_",counter)

                    file_name_individuals_i <- paste0(
                            output_file_prefix,"_individuals_",counter)

                    # add a line that spits out a number to the screen
                    # to inform the user how many simulations have been
                    # processed
                    echo_str <- paste("echo",counter)

                    command_str <- paste(exe,
                                    maxgen,
                                    L_i,
                                    a,
                                    p_nongenetic_i,
                                    p_maternal_i,
                                    dev_time_i,
                                    paste(s_values, collapse=" "),
                                    file_name_i,
                                    file_name_individuals_i)

                    # append to batch file contents
                    batch_file_contents <- paste0(batch_file_contents
                            ,"\n"
                            ,echo_str
                            ,"\n"
                            ,command_str)
                }
            } # end dev_time
      } # end p_nongenetic_i
    } # end L_i 
} # end for loop

writeLines(text=batch_file_contents)
