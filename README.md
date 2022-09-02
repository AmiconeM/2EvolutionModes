# 2EvolutionModes
[![DOI](https://zenodo.org/badge/486260829.svg)](https://zenodo.org/badge/latestdoi/486260829)

This project contains:

- code to simulate evolution under directional selection and Gamma-distributed mutation effects ("DirectionalSelection_Gamma.R").

- code to simulate evolution under competition for resources, adapted from Amicone and Gordo 2021, https://doi.org/10.1111/evo.14315 ("EcologicalSelection_Amicone2021.R").

- metadata file to run the R scripts ("Metadata_*.txt") under the two different models and the conditions reported in Extended Data Figure 10 of Frazao et al. 2022.

- outcome of simulations produced when running "DirectionalSelection_Gamma.R" with conditions specified in "Metadata_DirSel.txt".

Instructions to reproduce the row data presented in Extended Data Figure 10 of Frazao et al. 2022:
1) Install R (We used R version 3.6.1 (2019-07-05))

2) Access via terminal the folder containing the scripts and metadata files.
3) Run simulations from bash command line as:

bash Metadata_DirSel.txt
(this should take less than 10 minutes and produce the outcome of simulations under directional selection)

bash Metadata_EcoSel.txt
(this should take less than 1 hour and produce the outcome of simulations under competition for resources)

4) Access the output folders to see the tables of raw data of mutation trajectories across population.

In order to produce the final statistics reported in Frazao et al. 2022, these raw data were processed by the code described in the supplementary materials of Frazao et al. 2022.

