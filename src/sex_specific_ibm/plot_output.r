#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}
# use tidyverse's sym() to transform strings to symbols 
transform.sym  <- function(x) 
{
    if (!is.null(x))
    {
        sym(x)
    }
}

distribution_file_prefix = "distribution_"

if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)


# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line
        ,col_names=c("name","value")
        )


data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

if (nrow(data.tibble) > 5000)
{
    data.tibble <- data.tibble[round(seq(1,nrow(data.tibble),length.out=1000)),]
}

row_params <- data.tibble.params$value
names(row_params) <- data.tibble.params$name

n_loci <- as.numeric(row_params["L"])

names_sbar <- c()
names_vbar <- c()
names_diagonal <- c()
names_lower_diagonal <- c()
names_upper_diagonal <- c()

name <- ""
for (locus_idx in 1:n_loci)
{
    names_sbar <- c(names_sbar, paste0("sbar",locus_idx))
    names_vbar <- c(names_vbar, paste0("v",locus_idx))
    for (locus_idx2 in 1:n_loci)
    {
        name <- paste0("w",locus_idx,locus_idx2)

        if (locus_idx2 == locus_idx)
        {
            names_diagonal <- c(names_diagonal, name)
        } else if (locus_idx2 < locus_idx)
        {
            names_upper_diagonal <- c(names_upper_diagonal, name)
        } else
        {
            names_lower_diagonal <- c(names_lower_diagonal, name)
        }
    }
}

jsonstuff <- paste0('[
    {"xvar" : "time",
        "yvar" : ["',paste(names_diagonal,collapse="\",\""),
        '"]
    },
    {
        "xvar" : "time",
        "yvar" : ["',paste(names_lower_diagonal, collapse="\",\""),
            '"]
    },
    {
        "xvar" : "time",
        "yvar" : ["',paste(names_upper_diagonal, collapse="\",\""),
            '"]
    },
    {
        "xvar" : "time",
        "yvar" : ["',paste(names_sbar, collapse="\",\""),
            '"]
    },
    {
        "xvar" : "time",
        "yvar" : ["',paste(names_vbar, collapse="\",\""),
            '"]
    }
]')

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list()

plot.list.idx <- 1

single.plot <- function(xvar, yvar, sub.data.tibble) {
    
    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=sub.data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        # get rid of aes_string like this: 
        # https://stackoverflow.com/questions/74414272/how-to-replace-the-deprecated-ggplot2-function-aes-string-accepting-an-arbitrar/74414389#74414389 

        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=xvar,
                        y=yvar_values),
                        FUN=transform.sym)

        # aes arguments for the geom_line() call
        line_args <- lapply(X=list(
                        colour=yvar_name),
                FUN=transform.sym)

        return(ggplot(data=sub.data
                ,mapping=aes(!!!plot_args)) + 
                    geom_line(mapping=aes(!!!line_args)))
    } else {
        
        # aes arguments for the ggplot() call
        plot_args <- lapply(X=list(
                        x=xvar,
                        y=yvar
                        ),
                FUN=transform.sym)

        return(ggplot(data=sub.data.tibble
                ,mapping=aes(!!!plot_args)) + geom_line())
    }
} # end single.plot

# first the 'normal' plots over the whole of
# evolutionary time
# then later on during the more 'dense' data sampling
# period around the change point
for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    xvar <- unlist(plot.structure[[plot_struct_idx]]$xvar)
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    plot.list[[plot.list.idx]] <- single.plot(xvar, yvar, data.tibble)

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

base <- basename(file.name)
dir <- dirname(file.name)

distribution_file_name <- paste0(dir,"/",distribution_file_prefix,base)

title <- ""

wrap_plots(plot.list,ncol=2, byrow=F)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height=20, width = 20)


