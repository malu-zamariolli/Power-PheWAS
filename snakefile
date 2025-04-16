############# Power in PheWAS ########################
# A. Config file
configfile: "config_power.yaml"

# B. Define wildcard
wildcard_constraints:
     snp = config["snp"]

# C. non-wildcard definitions
trait_n = config["trait_n"]
unique_n = config["unique_n"]

# a per SNP approach


#------------------------ Rule all -------------------------------------------#
rule all:
    input:
        expand("{snp}/{snp}_all_PowerCalc_TraitName.txt", snp = config["snp"]),
        expand("{snp}/{snp}_all_PowerCalc.txt", snp = config["snp"]),
        expand("{snp}/simulation_p/plots/SimEffect_SimP_plots.pdf", snp = config["snp"])


#------------------------ Rules -------------------------------------------#

# Splits all possible combinations of MAF, effect size and N into single dfs for parallelization
rule split:
    input:
        n_traits = f'{unique_n}'
    output:
        linear = temp(expand("{snp}/temp/{chunk}_linear.txt", snp = config["snp"], chunk = range(1, int(config["n_linear"]) + 1))),
        logistic = temp(expand("{snp}/temp/{chunk}_logistic.txt", snp = config["snp"], chunk = range(1, int(config["n_logistic"]) + 1))),
        ordinal = temp(expand("{snp}/temp/{chunk}_ordinal.txt", snp = config["snp"], chunk = range(1, int(config["n_ordinal"]) + 1)))
    params:
        maf = config["maf"],
        effect = config["effect"],
        alpha = config["alpha"],
        n_linear = config["n_linear"],
        n_logistic = config["n_logistic"], 
        n_ordinal = config["n_ordinal"]
    script:
        "scripts/Split_runs.R"

# Performs Power calculation in linear regression
rule linear:
    input:
        n = "{snp}/temp/{chunk}_linear.txt"
    output:
        res = temp("{snp}/temp/{chunk}_linear_power_calc.txt"),
        plot_input = temp("{snp}/simulation_p/{chunk}_linear_simulation.txt")
    params:
        n_sim = config["n_sim"]
    benchmark:
        "{snp}/benchmarks/chunks/{chunk}.linear.benchmark.txt"
    resources:
        mem_mb=9000
    script:
        "scripts/SimulationPower_linearRegression.R"

# Performs power calculation in ordinal regression
rule ordinal:
    input:
        n = "{snp}/temp/{chunk}_ordinal.txt"
    output:
        res = temp("{snp}/temp/{chunk}_ordinal_power_calc.txt"),
         plot_input = temp("{snp}/simulation_p/{chunk}_ordinal_simulation.txt")
    params:
        n_sim = config["n_sim"]
    benchmark:
        "{snp}/benchmarks/chunks/{chunk}.ordinal.benchmark.txt"
    resources:
        mem_mb=10000
    script:
        "scripts/SimulationPower_ordinalRegression.R"

# Performs power calculation in binary logistic regression
rule logistic:
    input:
        n = "{snp}/temp/{chunk}_logistic.txt"
    output:
        res = temp("{snp}/temp/{chunk}_logistic_power_calc.txt"),
        plot_input = temp("{snp}/simulation_p/{chunk}_logistic_simulation.txt")
    params:
        n_sim = config["n_sim"]
    benchmark:
        "{snp}/benchmarks/chunks/{chunk}.logistic.benchmark.txt"
    resources:
        mem_mb=9000
    script:
        "scripts/SimulationPower_logisticRegression.R"

# Combines all results from power calculations
rule combine_power_calc:
    input:
        linear = expand("{snp}/temp/{chunk}_linear_power_calc.txt", snp = config["snp"], chunk = range(1, int(config["n_linear"] + 1))),
        logistic = expand("{snp}/temp/{chunk}_logistic_power_calc.txt", snp = config["snp"], chunk = range(1, int(config["n_logistic"] + 1))),
        ordinal = expand("{snp}/temp/{chunk}_ordinal_power_calc.txt", snp = config["snp"], chunk = range(1, int(config["n_ordinal"] + 1)))
    output:
        res_all = "{snp}/{snp}_all_PowerCalc.txt"
    script:
        "scripts/CombinePowerResults.R"

# Adds trait names to corresponding sample sizes and filter power
rule filter_power:
    input:
        res_all = "{snp}/{snp}_all_PowerCalc.txt", 
        trait_n = f'{trait_n}'
    output:
        res_all_traits = "{snp}/{snp}_all_PowerCalc_TraitName.txt",
        res_traits_filter = "{snp}/{snp}_all_PowerCalc_TraitName_filtered.txt"
    params:
        power = config["power"]
    script:
        "scripts/FilterPower.R"

# Plot effect size x p-value of all simulations
rule combine_plot_input:
    input:
        linear = expand("{snp}/simulation_p/{chunk}_linear_simulation.txt", snp = config["snp"], chunk = range(1, int(config["n_linear"] + 1))),
        ordinal = expand("{snp}/simulation_p/{chunk}_ordinal_simulation.txt", snp = config["snp"], chunk = range(1, int(config["n_ordinal"] + 1))),
        logistic = expand("{snp}/simulation_p/{chunk}_logistic_simulation.txt", snp = config["snp"], chunk = range(1, int(config["n_logistic"] + 1)))
    output:
        plot = "{snp}/simulation_p/plots/SimEffect_SimP_plots.pdf"
    script:
        "scripts/Combine_SimEffect_SimP_Plot.R"


#----------------------------- Save config --------------------------------------------------#
#Copy config file to corresponding folder, so config informations are saved to each run
rule config:
    input:
        "config_power.yaml"
    output:
        "{snp}/config/config_power.yaml"
    params:
        path = "{snp}/config/"
    shell:
        """
        cp {input} {params.path}
        """
