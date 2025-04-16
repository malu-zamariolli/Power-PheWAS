##### Simulation to calculate Power for Linear Regression

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
# Sample size
N <- n$N
# MAF
MAF <- n$MAF 
# effect size (beta)
b_snp <- n$effect
# alpha
alpha <- n$alpha

###########################
### Simulation ####
############################
# Set seed 
set.seed(123)

# Initialize vector to store p-values for simulations
p_val <- numeric(n_sim)
effect_size <- numeric(n_sim)
      
# Simulation loop
print("----- Starting simulation loop -----------")
for (i in 1:n_sim) {
    # Create a temporary data frame for SNP genotypes and phenotype
    df <- data.frame(SNP = numeric(N), error = numeric(N), Phenotype = numeric(N))
        
    # Simulate SNP genotypes based on MAF
    df$SNP <- runif(N, min = 0, max = 1)
    freq_AA <- (1 - MAF)^2
    freq_AB <- 2 * (1 - MAF) * MAF
    freq_BB <- MAF^2
    df$SNP[df$SNP >= (1 - freq_BB)] <- 2
    df$SNP[df$SNP < (1 - freq_BB) & df$SNP > freq_AA] <- 1
    df$SNP[df$SNP <= freq_AA] <- 0
        
    # Simulate error and phenotype
    var_error <- 1 - var(df$SNP) * b_snp^2
    df$error <- rnorm(N, mean = 0, sd = sqrt(var_error))
    df$Phenotype <- df$SNP * b_snp + df$error
        
    # Apply Inverse Rank Normal Transformation
    df$Phenotype <- qnorm((rank(df$Phenotype, na.last = "keep") - 0.5) / N)
        
    # Perform linear regression and store p-value
    p_val[i] <- summary(lm(Phenotype ~ SNP, data = df))$coefficients[2, 4]
    effect_size[i] <- summary(lm(Phenotype ~ SNP, data = df))$coefficients[2, 1]
}

print("----- End of simulation loop -----------")
      
###########################
### Obtain results ####
############################
# Calculate power from simulations
sim_power <- mean(p_val < alpha, na.rm = TRUE)
sim_power <- round(sim_power * 100, digits = 2)

# Calculate percentage of simulations with higher effect size
high_effect <- mean(effect_size >= b_snp, na.rm = TRUE)
      
# Obtain mean, min and max effect size
mean_eff <- mean(effect_size, na.rm = TRUE)
min_eff <- min(effect_size, na.rm = TRUE)
max_eff <- max(effect_size, na.rm = TRUE)

results <- data.frame(N = N, 
    MAF = MAF, EffectSize = b_snp, alpha = alpha,  
    PowerPercentage = sim_power, MeanSimEffect = mean_eff, MinSimEffect = min_eff,
    MaxSimEffect = max_eff, method = "LinearRegression")

# OR and p-values for plotting
df_plot <- data.frame(EffectSize = effect_size, pvalue = p_val, N = N, 
    MAF = MAF, TestedEffectSize = b_snp, alpha = alpha, Power = sim_power, method = "LinearRegression")

#######################
### Save Results ###
######################
fwrite(results, res, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")

fwrite(df_plot, plot_input, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")