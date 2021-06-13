# Temperature-specific model generation and analysis

This is an outline of the procedure used perform simulations on the model
of *Shewanella psychrophila* WP2.

All functions can be run as detailed in this document from the root folder of
the WP2 model git repository.

## Determining ATP maintenance fluxes

The WP2 model is adjusted to match the
experimental biomass production rates through the addition of an ATP
maintenance reaction. Flux is forced through this ATP maintenance reaction
by a flux limits setting until the biomass production in the model matches
the experimental values for each temperature.

To determine the amount of flux through the ATP maintenance reaction 'ATPM',
was increased at increments of 0.05 flux, calculating the biomass flux in the
at each step. A script named 'TMFA-robust.py' is provided in the scripts/
directory of the WP2 model repository. The script can be run by providing the
path to the model.yaml file, path to the configuration file (in the
tmfa-inputs/ directory), temperature, and the IDs of the biomass reaction and
ATP maintenance reaction). When running this script make sure that no flux
limits are currently placed on the ATP maintenance reaction.

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

This script will generate a two column table, with the first column having
the flux through the ATP maintenance reaction and the second column having
the corresponding biomass flux. A linear fitting of the ATP maintenance
flux vs the biomass fluxes was then done to determine the value of ATP
maintenance flux for each temperature that would result in the model
matching the experimental biomass production. The fitting and estimation of
the ATP maintenance flux can be done using theATPM-fit.R script by providing
the script with the directory that contains the output from TMFA-Robust.py.

  `Rscript scripts/ATPM/ATPM-fit.R ./scripts/ATPM/`

The output of this script will be a PDF showing the plots of the ATP maintenance
fluxes ('scripts/ATPM/ATPM-fittings.pdf') and a table will be printed to the
standard output that gives the estimated ATP maintenance fluxes for each
temperature:

  ` [1,] "4C"  "2.77793659458319"
    [2,] "15C" "1.33139444782041"
    [3,] "20C" "9.40642097520033"`

Flux limits files with the
final values for the ATP maintenance reaction are provided in the
limits/ directory of the WP2 model repository.


## Simulating Growth of WP2

Simulations of growth for each temperature were performed using the TMFA
approach as implemented in PSAMM. The input files for running these
simulations are provided in the tmfa-inputs/ directory of the WP2 model
repository. Running the TMFA simulations can be done through the following
commands:

  `psamm-model --model model-04C.yaml tmfa --solver threads=1 \
    --config tmfa-inputs/config-04C.yaml \
    --temp 4 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/04C-TMFA.tsv`
  `

  `psamm-model --model model-15C.yaml tmfa --solver threads=1 \
    --config tmfa-inputs/config-15C.yaml \
    --temp 15 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/15C-TMFA.tsv`

  `psamm-model --model model-20C.yaml tmfa --solver threads=1 \
    --config tmfa-inputs/config-20C.yaml \
    --temp 20 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/20C-TMFA.tsv`

*NOTE: It is suggested that these simulations be run with the --solver threads=1 option
if you are using the CPLEX linear programming solver. 

The output from these commands will be the variability analysis results for
each of the variables in the TMFA problem. Results for each variable will be
printed on separate lines in the following order: variable type, variable ID,
lower bound, and upper bound. Variables without any solved value will have 'NA'
printed in the lower rand upper bound columns (for example on the DeltaG line
for reactions that are not thermodynamically constrained). Pre-generated output
is provided in the scripts/TMFA-results/ directory of the WP2 model repository.

The outputs from each of the three temperatures can be condensed into one
table using the 'TMFA-Condense.py' script provided in the script/ directory of
the WP2 model repository. This script will take all three simulation results as
inputs and then produce a set of tables that summarize the TMFA results:

  `python ./scripts/TMFA-Condense.py \
    --four ./scripts/TMFA-results/04C-TMFA.tsv --fifteen ./scripts/TMFA-results/15C-TMFA.tsv \
    --twenty ./scripts/TMFA-results/20C-TMFA.tsv --out ./scripts/TMFA-results/TMFA-condensed`

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


# Comparing Flux Variability

Flux variability between the models can be plotted and compared using a Pairwise Wilcoxan Rank Sum test using the 'scripts/variability/flux-variability.R' script. This will produce the plot used as Figure S4.

  `Rscript scripts/variability/flux-variability.R ./scripts/variability/`


# Plotting Concentration, Flux, and DeltaG Ranges

The variability of selected metabolite concentrations, fluxes, and Gibbs free energy values shown in Supplemental Figures 1-3 can be plotted using the 'scripts/variability/plot-fluxes-barbell.R'

  `Rscript scripts/variability/plot-fluxes-barbell.R ./scripts/variability/`
