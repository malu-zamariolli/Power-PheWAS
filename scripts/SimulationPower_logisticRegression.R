##### Simulation to calculate Power for Logistic Regression

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
N <- stringr::str_split(N, pattern = ";")[[1]]
N_controls <- as.numeric(N[[1]])
N_cases <- as.numeric(N[[2]])

###########################
### Simulation ####
############################
# Set seed 
set.seed(123)

# Initialize vector to store p-values for simulations
p_val <- numeric(n_sim)
effect_size <- numeric(n_sim)

        
# Frequencies of each genotype based on MAF
freq_AA <- (1 - MAF)^2
freq_AB <- 2 * MAF * (1 - MAF)
freq_BB <- MAF^2
        

# Simulation loop
print("----- Starting simulation loop -----------")        
for (i in 1:n_sim) {
        # Simulate SNP genotypes for cases and controls separately
        df_cases <- data.frame(SNP = runif(N_cases))
        df_controls <- data.frame(SNP = runif(N_controls))
          
        df_cases$SNP[df_cases$SNP >= (1 - freq_BB)] <- 2
        df_cases$SNP[df_cases$SNP <= freq_AA] <- 0
        df_cases$SNP[df_cases$SNP < (1 - freq_BB) & df_cases$SNP > freq_AA] <- 1
          
        df_controls$SNP[df_controls$SNP >= (1 - freq_BB)] <- 2
        df_controls$SNP[df_controls$SNP <= freq_AA] <- 0
        df_controls$SNP[df_controls$SNP < (1 - freq_BB) & df_controls$SNP > freq_AA] <- 1
          
        # Combine cases and controls
        df <- rbind(df_cases, df_controls)
          
        # Simulate probability of disease based on genotype based on ESPRESSO
        pheno.prev <- N_cases/(N_cases + N_controls)
        intercept <- log(pheno.prev / (1 - pheno.prev))
        lp <- intercept + b_snp * df$SNP  # Linear predictor
        mu <- exp(lp) / (1 + exp(lp))  # Logistic transformation
          
        # Generate binary phenotype (case/control)
        df$Phenotype <- rbinom(N_cases + N_controls, 1, mu)
          
        # Fit logistic regression model and store p-value
        mod <- glm(Phenotype ~ SNP, family = binomial(link = "logit"), data = df)
        p_val[i] <- summary(mod)$coefficients[2, 4]
        effect_size[i] <- summary(mod)$coefficients[2, 1]
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

results <- data.frame(N = paste0(N_controls, ";", N_cases), 
    MAF = MAF, EffectSize = b_snp, alpha = alpha,  
    PowerPercentage = sim_power, MeanSimEffect = mean_eff, MinSimEffect = min_eff,
    MaxSimEffect = max_eff, method = "LogisticRegression")

# OR and p-values for plotting
df_plot <- data.frame(EffectSize = effect_size, pvalue = p_val, N = N, 
    MAF = MAF, TestedEffectSize = b_snp, alpha = alpha, Power = sim_power, method = "LogisticRegression")

#######################
### Save Results ###
######################
fwrite(results, res, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")

fwrite(df_plot, plot_input, quote = FALSE, col.names = TRUE, row.names = FALSE, na = NA, sep = "\t")
        

