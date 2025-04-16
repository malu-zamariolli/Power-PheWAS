#### Filter Power Results and add trait Names


######################
#### Packages ########
#######################
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")

shh(library(data.table))
shh(library(dplyr))

###########################
### Inputs and outpus ####
############################
res_all <- snakemake@input[["res_all"]]
trait_n <- snakemake@input[["trait_n"]]

# Params
power <- snakemake@params[["power"]]

# Output
res_all_traits  <- snakemake@output[["res_all_traits"]]
res_traits_filter <- snakemake@output[["res_traits_filter"]]

###########################
### Load inputs ####
############################
results <- as.data.frame(fread(res_all))
trait_n <- as.data.frame(fread(trait_n))

###########################
### Combine results with trait_names ####
############################
df <- inner_join(results, trait_n)


###########################
### Filter by power ####
############################
df.filtered <- df[df$PowerPercentage >= power, ]

#############
### Save ####
##############
fwrite(df, res_all_traits, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t", na = NA)
fwrite(df.filtered, res_traits_filter, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t", na = NA)