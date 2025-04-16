########## Script to create input for Power ###################
# Combine MAF, effect size and N

######################
#### Packages ########
#######################
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(stringr)) install.packages("stringr")

shh(library(data.table))
shh(library(dplyr))
shh(library(stringr))

###########################
### Inputs and outpus ####
############################
# Inputs
n_traits <- snakemake@input[["n_traits"]]

# Params
maf <- snakemake@params[["maf"]]
effect <- snakemake@params[["effect"]]
alpha <- snakemake@params[["alpha"]]
n_linear <- snakemake@params[["n_linear"]]
n_logistic <- snakemake@params[["n_logistic"]]
n_ordinal <- snakemake@params[["n_ordinal"]]

# Output
linear <- snakemake@output[["linear"]]
logistic <- snakemake@output[["logistic"]]
ordinal <- snakemake@output[["ordinal"]]

###########################
### Import inputs ####
############################
n_traits <- as.data.frame(fread(n_traits))

maf <- stringr::str_split(maf, pattern = ";")[[1]]
effect <- stringr::str_split(effect, pattern = ";")[[1]]
alpha <- stringr::str_split(alpha, pattern = ";")[[1]]

##################################################
### Create all combinations of parameters ####
##################################################

# df with N, method, maf, effect, alpha - get all combinations
df.combined <- data.frame(N = numeric(), MAF = numeric(), alpha = numeric(), effect = numeric())
for (i in 1:nrow(n_traits)) {  
    for (m in maf) {       
        for (e in effect) {
            for (a in alpha) {
                df <- n_traits[i, ]
                df$MAF <- m
                df$effect <- e
                df$alpha <- a 
                df.combined <- rbind(df.combined, df)              
            }
        }
    }    
}

#################
### Split df ####
#################
# Split data by method
df.linear <- df.combined[df.combined$method %in% "LinearRegression", ]
df.logistic <- df.combined[df.combined$method %in% "LogisticRegression", ]
df.ordinal <- df.combined[df.combined$method %in% "OrdinalRegression", ]

# Function to divide dataframe per row
split_dataframe <- function(data, n_chunk, output_paths) {
  for (i in 1:n_chunk) {
    chunk_data <- data[i, ]
    fwrite(chunk_data, file = output_paths[i], na = NA, sep = "\t", quote = FALSE)
  }
}

#################
### Save dfs ####
#################
split_dataframe(df.linear, n_linear, linear)
split_dataframe(df.logistic, n_logistic, logistic)
split_dataframe(df.ordinal, n_ordinal, ordinal)