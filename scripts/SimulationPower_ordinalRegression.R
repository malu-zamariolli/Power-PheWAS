##### Simulation to calculate Power for Logistic Regression

######################
#### Packages ########
#######################
shh <- suppressPackageStartupMessages # https://stackoverflow.com/questions/18931006/how-to-suppress-warning-messages-when-loading-a-library

if(!require(data.table)) install.packages("data.table")
if(!require(dplyr)) install.packages("dplyr")
if(!require(stringr)) install.packages("stringr")
if(!require(MASS)) install.packages("MASS")
if(!require(car)) install.packages("car")

shh(library(data.table))
shh(library(dplyr))
shh(library(stringr))
shh(library(MASS))
shh(library(car))

###########################
### Inputs and outpus ####
############################
n <- snakemake@input[["n"]]

# Params
n_sim <- snakemake@params[["n_sim"]]

# Output
res <- snakemake@output[["res"]]
plot_input <- snakemake@output[["plot_input"]]


###########################
### Load input ####
############################
n <- as.data.frame(fread(n))

###########################
### Params for power ####
############################
# MAF
MAF <- n$MAF 
# effect size (beta)
b_snp <- n$effect
# alpha
alpha <- n$alpha

# Sample size
N <- n$N
group_sizes <- as.numeric(stringr::str_split(N, pattern = ";")[[1]])
group_labels <- 1:length(group_sizes)

N <- sum(group_sizes)

###########################
### Simulation ####
############################
# Simulation function
run_simulation <- function(N, MAF, b_snp, alpha, group_sizes) { 
  # Calculate cumulative proportions from group sizes
  cumulative_proportions <- cumsum(group_sizes) / sum(group_sizes)  
  # Calculate thresholds based on cumulative proportions
  thresholds <- qnorm(cumulative_proportions[-length(cumulative_proportions)]) # define k-1 thr, k is number of groups  
  # Create a temporary data frame for storing SNP genotypes and phenotype
  df <- as.data.frame(matrix(nrow = N, ncol = 4, dimnames = list(NULL, c("SNP", "error", "Phenotype", "OrdinalPhenotype"))))  
  # Simulate SNP genotypes: 0 for AA, 1 for AB, 2 for BB based on MAF 
  df[, "SNP"] <- runif(N, min = 0, max = 1)
  # Frequencies of each genotype
  freq_AA <- (1 - MAF)^2
  freq_AB <- 2 * (1 - MAF) * MAF
  freq_BB <- MAF^2
  # Atribute SNPs to uniform distribution
  df[which(df$SNP >= (1-freq_BB)), "SNP"] <- 2
  df[which(df$SNP <= freq_AA), "SNP"] <- 0
  df[which(df$SNP < (1-freq_BB) & df$SNP > freq_AA), "SNP"] <- 1
  
  # Count genotypes
  # Calculate genotype counts
  genotype_counts <- table(df$SNP)
  ma_count <- as.numeric(genotype_counts[3])
  
  # Define minimum count threshold
  min_count_threshold <- 5
  
  # Check if any genotype group has fewer than the threshold
  if (is.na(ma_count) | ma_count < min_count_threshold) {
    return(NULL)
  } else {
    # Simulate the error: phenotype variance is assumed to be 1 and is composed of the variance explained by SNPs and the variance due to error
    var_error <- 1 - var(df$SNP) * b_snp^2
    df[, "error"] <- rnorm(N, mean = 0, sd = sqrt(var_error)) # standard deviation is the square root of the variance
    
    # Simulate phenotypes (Phenotype): Phenotype = SNP * b_snp + error
    df[, "Phenotype"] <- df[, "SNP"] * b_snp + df[, "error"]
    
    # Stratify phenotype into specific groups based on thresholds
    df$OrdinalPhenotype <- cut(df$Phenotype, breaks = c(-Inf, thresholds, Inf), labels = group_labels)
    
    # Convert to factor
    df$OrdinalPhenotype <- factor(df$OrdinalPhenotype, levels = group_labels)
    
    # Fit ordinal logistic regression
    model <- MASS::polr(OrdinalPhenotype ~ SNP, data = df, Hess = TRUE)
    
    # Proportional Odds assumption with error handling
    prop.odds <- tryCatch({
      car::poTest(model)
    }, error = function(e) {
      return(NULL)  # Return NULL if poTest fails
    })
    
    # If poTest fails, assign p_value_propOdds = 0
    p_value_propOdds <- if (is.null(prop.odds)) 0 else {
      chi_square_value <- prop.odds$chisq
      df_overall <- prop.odds$df
      1 - pchisq(chi_square_value, df_overall)
    }
    
    if (p_value_propOdds <= 0.05) {
      return(NULL)
    } else {
      # Store Odds Ratio from model
      effect_size <- summary(model)$coefficients[1, 1]  
      
      # p-value
      t_value <- summary(model)$coefficients[1, 3]
      p_val <- pnorm(abs(t_value), lower.tail = FALSE) * 2
      return(list(p_val = p_val, effect_size = effect_size))
    }
  }
}

## Run simulations
results <- list()
valid_simulations <- 0

# Set seed 
set.seed(123)

print("----- Starting simulation  -----------")  

while (valid_simulations < n_sim) {
  sim_result <- run_simulation(N, MAF, b_snp, alpha, group_sizes)
  if (!is.null(sim_result)) {
    valid_simulations <- valid_simulations + 1
    results[[valid_simulations]] <- sim_result
  }
}

print("----- End of simulation -----------")

###########################
### Obtain results ####
############################
# Extract p-values and Odds Ratios
p_val <- sapply(results, function(x) x$p_val)
effect_size <- sapply(results, function(x) x$effect_size)

# Calculate power from simulations
sim_power <- mean(p_val < alpha, na.rm = TRUE)
sim_power <- round(sim_power * 100, digits = 2)

# Calculate percentage of simulations with higher effect size
high_effect <- mean(effect_size >= b_snp, na.rm = TRUE)
      
# Obtain mean, min and max effect size
mean_eff <- mean(effect_size, na.rm = TRUE)
min_eff <- min(effect_size, na.rm = TRUE)
max_eff <- max(effect_size, na.rm = TRUE)

results <- data.frame(N = paste0(group_sizes, collapse = ";"), 
    MAF = MAF, EffectSize = b_snp, alpha = alpha,  
    PowerPercentage = sim_power, MeanSimEffect = mean_eff, MinSimEffect = min_eff,
    MaxSimEffect = max_eff, method = "OrdinalRegression")

# OR and p-values for plotting
df_plot <- data.frame(EffectSize = effect_size, pvalue = p_val, N = N, 
    MAF = MAF, TestedEffectSize = b_snp, alpha = alpha, Power = sim_power, method = "OrdinalRegression")

#######################
### Save Results ###
######################
fwrite(results, res, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")

fwrite(df_plot, plot_input, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")




