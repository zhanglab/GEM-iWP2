# Temperature-specific model generation and analysis

This is an outline of the procedure used perform simulations on the model
of *Shewanella psychrophila* WP2.

Please note that the data provided to run all the Rscripts is included but you
will need to set the working directory to be the directory as the directory that
the script is located in.

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
ATP maintenance reaction). *When running this script make sure that no flux
limits are currently placed on the ATP maintenance reaction.*

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
  `

  `psamm-model --model model-15C.yaml tmfa --single-solution random --solver threads=1 \
    --config tmfa-inputs/config-15C.yaml --threshold {value} \
    --temp 15 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/15C-TMFA.tsv`

  `psamm-model --model model-20C.yaml tmfa --single-solution random  --solver threads=1 \
    --config tmfa-inputs/config-20C.yaml --threshold {value} \
    --temp 20 --phin 6,8 --phout 6,8  simulation > ./scripts/TMFA-results/20C-TMFA.tsv`

*NOTE: It is suggested that these simulations be run with the --solver threads=1 option
if you are using the CPLEX linear programming solver.

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
