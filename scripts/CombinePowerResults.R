#### Combine results from Power Calculation


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
linear <- snakemake@input[["linear"]]
logistic <- snakemake@input[["logistic"]]
ordinal <- snakemake@input[["ordinal"]]

# Output
res_all  <- snakemake@output[["res_all"]]

########################
### Load inputs 
###########################
input_files <- c(linear, logistic, ordinal)
print(input_files)
print(linear)
print(logistic)
print(ordinal)

# Initialize an empty list to store data frames
data_frames <- list()

# Iterate over the input files and read them into data frames
for (file in input_files) {
  # Read the file into a data frame
  df <- as.data.frame(fread((file)))
  # Append the data frame to the list
  data_frames <- append(data_frames, list(df))
}


########################
### Combine dataframes
###########################
df.all <- do.call(rbind, data_frames) 


########################
### Save output
###########################
fwrite(df.all, res_all, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t", na = NA)