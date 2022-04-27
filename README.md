# 2EvolutionModes

This project contains:

- code to simulate evolution under directional selection and Gamma-distributed mutation effects ("Amicone_DirectionalSelection_Gamma.R").

- metadata file to run the R script ("Metadata.txt") under the specified conditions.

- outcome of simulations with Gamma(Shape=1), N=10^7, Ub=10^(-9)->10^(-5), mean(s)=0.05, generations=10000 ("DirSel_Simulations.rar").

- scripts to compute and plot the Sojourn times from simulation outcome ("Sojour_Analysis.Rmd")

- scripts to compute and plot the frequency propagator from simulation outcome ("Frequency_Proppagator.Rmd")

- scripts to compute and plot the p-tau from simulation outcome ("p-tau_Stat.Rmd")

- summary plots of sojourn times, frequency propagator and p-tau statistics ("*.PNG")

- tables of raw data points for the figures (".xlsx")

 

Instructions to reproduce the final figures:

1) run simulations from bash command line as:

bash Metadata

(this should produce the same outcome that is stored here in the .rar archive)

2) Run "Sojour_Analysis.Rmd" and "Frequency_Proppagator.Rmd" to compute, plot and export raw data points of the figures.

3) Run "p-tau_Stat.Rmd" to compute, plot and export raw data points of the figures.

 

NOTE: you can start reproducing the figures at any step (1-3)!

Example you can jump on step 3 and plot directly the p-tau statistics because you already have the excel table with the pre-computed data points.
