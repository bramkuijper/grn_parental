#!/usr/bin/env Rscript

# number of replicates
nrep <- 3

# maximum amount of generations a simulation runs
maxgen <- 10000

# output interval in the data file, i.e., 
# how many timesteps until the next statistics are written out
output_interval <- 10

# the number of phenotypes, L
# note that the GRN is size L*L
L <- 6

# expression in the absence of positive
# or negative regulation on gene expression
a <- 0.2

# see eq. (4)  in Odorico et al for values of theta
# the first two are given nonnegative values, the other 
# loci stabilize around 0 (if they are under selection at all,
# see the s vector) 
theta <- numeric(length = L)

# see Odorico et al Fig S3 for parameter values
theta[1] <- 0.33
theta[2] <- 2*theta[1]

# fraction of phenotype that is nongenetically inherited
p_nongenetic <- seq(0,1,0.1)

p_maternal <- seq(0,1,0.1)

# number of time steps before development is completed
dev_time <- c(24)

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

# first allocate 0-filled vector of s values
s_values <- numeric(length=L)

# set the first and the second one to nonzero values
s_values[1:6] <- 10

# selection on developmental stability
# of each of the loci 
sprime <- numeric(length = L)
sprime[1:6] <- 1000

batch_file_contents <- ""

command_str <- paste("#",
                     exe,
                     "L",
                     "theta",
                     "\t",
                     "maxgen",
                     "\t",
                     "pn,pm",
                     "dev_time",
                     "s",
                     "\t",
                     "sprime",
                     "\t",
                     "file",
                     "file_individuals")

batch_file_contents <- command_str

# now vary all the possible combinations
# and see what happens
for (rep_i in 1:nrep)
{
      for (p_nongenetic_i in p_nongenetic)
      {
	      p_maternal_x <- p_maternal
	      # if p_nongenetic is 0 no need
	      # to vary over all p_maternal combinations
	      # as there is no nongenetic effect anywayz
	      if (p_nongenetic_i == 0.0)
	      {
		      p_maternal_x <- c(0.0)
	      }

          for (p_maternal_i in p_maternal_x)
          {
              for (dev_time_i in dev_time)
              {
                  counter <- counter + 1
                  
                  # dynamically generate a filename
                  # for this particular parameter combination
                  file_name_i <- paste0(
                          output_file_prefix,"_",counter)

                  file_name_individuals_i <- paste0(
                          output_file_prefix,"_individuals_",counter)

                  # add a line that spits out a number to the screen
                  # to inform the user how many simulations have been
                  # processed
                  echo_str <- paste("echo",counter)
                  
                  command_str <- paste(exe,
                                  L,
                                  "\t",
                                  paste(theta, collapse=" "),
                                  "\t",
                                  maxgen,
                                  "\t",
                                  p_nongenetic_i,
                                  p_maternal_i,
                                  "\t",
                                  dev_time_i,
                                  paste(s_values, collapse=" "),
                                  "\t",
                                  paste(sprime,collapse=" "),
                                  "\t",
                                  output_interval,
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
} # end for loop

writeLines(text=batch_file_contents)


