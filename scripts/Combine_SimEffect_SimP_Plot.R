#### Combine results and plot
# Plot EffectSize vs p-values for all simulations


######################
#### Packages ########
#######################
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(cowplot)) install.packages("cowplot")
if(!require(gridExtra)) install.packages("gridExtra")

shh(library(data.table))
shh(library(dplyr))
shh(library(ggplot2))
shh(library(cowplot))
shh(library(gridExtra))

###########################
### Inputs and outpus ####
############################
linear <- snakemake@input[["linear"]]
logistic <- snakemake@input[["logistic"]]
ordinal <- snakemake@input[["ordinal"]]

# Output
output  <- snakemake@output[["plot"]]

########################
### Load inputs 
###########################
input_files <- c(linear, logistic, ordinal)

# Initialize an empty list to store data frames
data_frames <- list()

# Iterate over the input files and read them into data frames
for (file in input_files) {
  # Read the file into a data frame
  df <- as.data.frame(fread((file)))
  # Append the data frame to the list
  data_frames <- append(data_frames, list(df))
}

################################
#### Plots ####################
#############################

plot_list <- list()

for(i in 1:length(data_frames)) {
    df_plot <- data_frames[[i]]
    main_title <- paste0("N:", unique(df_plot$N), "-", "MAF", 
        unique(df_plot$MAF), "-", "P:", unique(df_plot$Power), "Eff:", unique(df_plot$TestedEffectSize)) 
    p1 <- ggplot(df_plot, aes(x = EffectSize)) +
        geom_histogram() +
        labs(title = main_title) +
        theme_minimal()
    p2 <- ggplot(df_plot, aes(x=pvalue, y = EffectSize)) +
        geom_point() +
        labs(title = main_title) +
        scale_x_continuous(breaks = c(0,0.05, 0.25, 0.5, 0.75, 1)) +
        theme_minimal()
    plot <- plot_grid(p1, p2, nrow = 1, ncol = 2)
    # Save list
    plot_list[[i]] <- plot
}

########################
### Save output
###########################

ggsave(filename = output, 
       plot = marrangeGrob(plot_list, ncol = 1, nrow = 1), 
       width = 15, height = 10)

