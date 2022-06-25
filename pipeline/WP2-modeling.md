# Temperature-specific model generation and analysis

Files in this repository represent the genome-scale model of the deep sea
psychrophilic bacterium, *Shewanella psychrophila* WP2. A detailed account on
the reconstruction and modeling of the WP2 model is summarized in the following
 manuscript:

Dufault-Thompson K§, Nie C, Jian H, Wang F, and Zhang Y. (2022)
Reconstruction and analysis of thermodynamically-constrained models reveal 
metabolic responses of a deep-sea bacterium to temperature perturbations. 

Please note that the data provided to run all the Rscripts is included but you
will need to set the working directory to be the directory as the directory that
the script is located in.

## List of model files
* model.yaml
  >Model definition file for the overall WP2 model

* model-??C.yaml
  >Model definition files for the temperature specific models.

* WP2-reactions.yaml
  >Reaction definition YAML file.

* compounds.tsv
  >Compound feature table for the model in tsv format.

* WP2-biomass.yaml
  >Biomass and macromolecule synthesis reactions.

* WP2-gene-associations.tsv
  >Gene associations for reactions in the model.

* WP2-Expression-TPM-table.tsv
  >Tab separated table of TPM values for all genes in the WP2 genome based
on WP2 transcriptomes from growth at 4C, 15C, and 20C in late log phase.

* model-def.tsv
  >Model definition files defining the reactions present in each of the
temperature based models.

* exchange/*
  >Exchange files representing the LMO 812 media, with temperature specific
   representations of the O2 availability.

* limits/*
  >Flux limits settings for each temperature based model.

* tmfa-inputs/*
  >Input files to run TMFA simulations on temperature specific models.

* scripts/*
  >scripts used for the simulation of temperature-specific models.

## Outline of the procedure used to simulate temperature specific models

A [PSAMM virtual environment](https://psamm.readthedocs.io/en/stable/tutorial/psamm-install.html#psamm-installation) should be set up before running the following procedure.

### Determining ATP maintenance fluxes

The temperature specific WP2 model uses parameterization that fits the
experimentally measured biomass production rates. Specifically the fluxes of
an ATP maintenance reaction was fitted to the temperature-specific growth rate
measured from experimental studies to model the non-growth-related maintenance
cost under each temperature. The fitted ATP maintenance flux for each
temperature condition is summarized in the corresponding YAML files under the
`limits/` directory, and a brief reference of the fitting procedure is described
below:

Under each temperature condition, a robustness analysis was performed using
the TMFA protocol, where fluxes of the ATP maintenance reaction `ATPM`
was increased at increments of 0.05 while the biomass flux was maximized at
each step. This analysis was implemented in the script `TMFA-robust.py`,
provided in the `scripts/` directory. The script takes a number of
input parameters, including the path to the model.yaml file, the varying
reaction (in robustness analysis), the biomass or objective reaction that is
optimized in each step of the robustness analysis, the path to the TMFA
configuration file (in the tmfa-inputs/ directory), temperature (in Celsius),
the lower and upper bounds of intracellular pH, and the lower and upper bounds
of extracellular pH.

***NOTE: When running this script make sure that no flux limits are placed
on ATPM***

  `python ./scripts/TMFA-Robust.py \
    --model model-04C.yaml --varying ATPM --biomass Biomass_WP2 \
    --config tmfa-inputs/config-04C.yaml --temp 4 --phin 6,8 --phout 6,8 \
    > scripts/ATPM/04C-ATPM-Robustness.tsv`

  `python ./scripts/TMFA-Robust.py \
    --model model-15C.yaml --varying ATPM --biomass Biomass_WP2 \
    --config tmfa-inputs/config-04C.yaml --temp 15  --phin 6,8 --phout 6,8 \
    > scripts/ATPM/15C-ATPM-Robustness.tsv`

  `python ./scripts/TMFA-Robust.py \
    --model model-20C.yaml --varying ATPM --biomass Biomass_WP2 \
    --config tmfa-inputs/config-04C.yaml --temp 20  --phin 6,8 --phout 6,8 \
    > scripts/ATPM/20C-ATPM-Robustness.tsv`

The commands above will generate a two column table, with the first column
representing the varying flux through the ATP maintenance reaction and the
second column showing the maximized biomass flux in each step. (Note: the
solver may generate slightly variable results in solving the TMFA problems.
We observed a fluctuation of less than 1E-12 across different computers and
less than 1E-6 across different Cplex versions in the solution.)

A linear fitting
of the ATP maintenance flux vs the biomass fluxes was then done to determine
the ATP maintenance cost for each temperature that would result in the model
matching the experimental biomass production. The fitting and estimation of
the ATP maintenance flux can be done using the ATPM-fit.R script by providing
the script with the directory that contains the output from TMFA-Robust.py.

The Rscript and data used to perform the fitting of the growth curve and ATPM robustness
simulations is in the ./scripts/ATPM directory. The script plot_Fig1.r will read
in the associated tables from the directory and plot out a summarized growth curve
and a plot of the ATPM robustness simulations for each temperature. The script
will produce summarized data tables and a PDF image file Figure1.pdf.

The range of ATPM values for each temperature was determined to be:

   | Temp   | Min   | Max   |
   |--------|-------|-------|
   | 4C     | 2.71  | 2.82  |
   | 15C    | 1.32  | 1.36  |
   | 20C    | 7.41  | 11.55 |

Flux limits files with the
final values for the ATP maintenance reaction are provided in the
limits/ directory of the WP2 model repository.


## Simulating Growth of WP2 using random simulations

Simulations of growth for each temperature were performed using the TMFA
approach as implemented in PSAMM. The input files for running these
simulations are provided in the tmfa-inputs/ directory of the WP2 model
repository. Running the TMFA simulations can be done through the following
commands:

  `psamm-model --model model-04C.yaml tmfa --single-solution random --solver threads=1 \
    --config tmfa-inputs/config-04C.yaml --threshold {value} \
    --temp 4 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/04C-TMFA.tsv`

  `psamm-model --model model-15C.yaml tmfa --single-solution random --solver threads=1 \
    --config tmfa-inputs/config-15C.yaml --threshold {value} \
    --temp 15 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/15C-TMFA.tsv`

  `psamm-model --model model-20C.yaml tmfa --single-solution random  --solver threads=1 \
    --config tmfa-inputs/config-20C.yaml --threshold {value} \
    --temp 20 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/20C-TMFA.tsv`

*NOTE: It is suggested that these simulations be run with the --solver threads=1 option
if you are using the CPLEX linear programming solver.*

The output from these commands will be the variability analysis results for
each of the variables in the TMFA problem. Results for each variable will be
printed on separate lines in the following order: variable type, variable ID,
lower bound, and upper bound. Variables without any solved value will have 'NA'
printed in the lower and upper bound columns (for example on the DeltaG line
for reactions that are not thermodynamically constrained). Pre-generated output
is provided in the scripts/TMFA-results/ directory of the WP2 model repository.

The outputs from each of the three temperatures can be condensed into one
table using the 'TMFA-Condense.py' script provided in the script/ directory of
the WP2 model repository. This script will take all three simulation results as
inputs and then produce a set of tables that summarize the TMFA results:

  `python ./scripts/TMFA-Condense.py \
    --four ./scripts/TMFA-results/04C-TMFA.tsv \
    --fifteen ./scripts/TMFA-results/15C-TMFA.tsv \
    --twenty ./scripts/TMFA-results/20C-TMFA.tsv \
    --out ./scripts/TMFA-results/TMFA-condensed`

The script will produce six files that summarize the TMFA results across
all three temperatures. The 'Flux_Variability.tsv' file summarizes the flux
comparisons across the temperatures. The 'Compound_Variability.tsv' file
summarizes the concentration differences across the temperature simulations.
The 'deltaG_Variability.tsv' summarizes the Gibbs Free Energy Change
differences between the temperatures. The 'carbNorm_Variability.tsv' file
provides a summary of the fluxes with all reaction fluxes normalized by
the carbon uptake flux. the 'absMin_Variability.tsv' file provides a summary
of the fluxes normalized by biomass flux for each temperature.
The last 'SummaryCounts.tsv' file is an overall
summary of the changes between temperatures for the fluxes, compound
concentrations, and Gibbs Free Energy changes.

## Calculating Metabolic Efficiency Metrics

The metabolic efficiency metrics are calculated based on the flux simulation
results for each temperature. A script, 'efficiency-plotting.R', to calculate
these metrics is provided in the scripts/efficiency directory of the WP2 model
repository. This script will calculate the efficiency metrics and plot out
the metrics for each temperature based on the flux results from the
TMFA simulations.her

  `Rscript scripts/efficiency/efficiency-plotting.R ./scripts/efficiency/`

The out from running these commands is a single random solution from within the
solution space of the TMFA problem with the biomass fixed to the value specified
with the '--threshold' option. Each variable in the problem (flux, deltaG, and
compound concentration) is reported with its value for this single solution in
the output of these commands. For the simulations performed in the study the
biomass was fixed at a random value between the min and max biomass of each model:

   | Temp   | Min    | Max   |
   |--------|--------|-------|
   | 4C     | 0.110  | 0.112 |
   | 15C    | 0.254  | 0.255 |
   | 20C    | 0.042  | 0.135 |

For each temperature 1000 random simulations were performed and then the results
were summarized into summary tables provided in the ./scripts/TMFA-results/
directory. This directory contains the files needed to run the subsequent analyses

## Calculating Metabolic Efficiency Metrics

The metabolic efficiency metrics are calculated based on the flux simulation
results for each temperature. The scripts and data associated with the generating
the efficiency metrics is included in the ./script/efficiency/ directory. The
calculate-effic.py script takes the 'atp-rxns', 'carbon-exchange', temperature,
and path to a directory containing random simulation results as inputs and produces
a summarized file containing the CUE, ATP per biomass, and ATP per carbon values
for each random simulation in the directory supplied. This scrip produces a table
like the one in the random-effic-values.tsv file.

The 'plot_Fig2.R' script can then be used to plot out boxplots showing the distributions
of each efficiency metric based on the results in the random-effic-values.tsv file.
This script produces a PDF image file called Figure2.pdf.

## Analysis of Flux, DeltaG, and Metabolite distributions
R scripts detailing the workflow for the analysis of the simulation results are
included in the ./scripts/analysis/ directory. The data read in by these analysis
scripts is all supplied in the ./scripts/TMFA-results/ directory. This directory
contains tables that summarize the compound concentration, deltaG, and flux
values from the random simulations of each temperature. They also include the
gene count table associated with the WP2 transcriptome.

The first script '1_flux_gibbs_expression.R', include the analysis workflow related to the distribution
of fluxes and deltaG values, along with the comparisons of the modeling and differential
expression results. The second script '2_cpd_analysis.R' includes the analysis of the
metabolite concentrations across the random simulations. These scripts will both produce
summary tables that provide information on the summarized values for each simulation variable
alongside statistical test results that were used to identify the significant differences between
temperatures. The scripts also produce multiple figures that summarize the comparison between
fluxes and deltaG values, look at gene expression in different conditions, and compare the
ratios of different metabolite concentrations.