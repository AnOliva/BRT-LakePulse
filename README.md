# BRT-LakePulse

Author: Anaïs Oliva

Correspondence: anais.oliva@usherbrooke.ca

## Overview - Usage
This repository contains scripts used in Oliva *et al.* (*Water Research*, submitted and under review).
The study objectives were (1) to explore the diversity of Potentially Pathogenic Bacteria (PPB); (2) to build a fecal multi-indicator from a cluster of co-occurring PPB; and (3) to predict the fecal multi-indicator over thousands of lakes. 

Boosted Regression Tree (BRT) models were applied over 1000 bootstrap samples through a boostrap aggregating (or bagging) ensemble method. This helped to determine the most influent environmental variables related to the abundance of the bacterial clusters and to make predictions.

The folder provides all the script used from the raw data recovered and pre-processed to the BRT model outputs.

## Folder contents

1. *1_microbes_preprocessing*
   - *1_rarefying.R*: rarefaction of sampled 16S rRNA sequences.
   - *2_pathogen_extraction.R* : extraction of the potential pathogens based on the ePathogen datase and a partial match algorithm.
   - *3_clustering_analysis.R*: clustering analysis to form groups of similar pathogens based on Sorensen dissimilarity index.
   
3. *2_BRT_preprocessing*
   - *1_combined_variables.R*: load all the independent variables, verify the flags and normalize the data.
   - *2_collinearity_CDOM-TSS.R*: check Pearson correlation and the generalized variance inflation factor (GVIF) for TSS and CDOM datasets.
   - *2_collinearity_QC.R*: check Pearson correlation and the GVIF for the quality controlled dataset .
   - *2_collinearity_upscaled.R*: check Pearson correlation and the GVIF for the upscaled dataset.
    
5. *BRT* (either *BRT_aCDOM* or *BRT_species*)
   - *1_normalization.R*: set/edit the destination directory for BRT results, loads required datasets and make preparation for the BRT scripts.
   - *2_BRT_preliminary.R*: calculate 5 consecutive BRT for each combination of BRT parameters defined in an hypergrid. Select the best combination.
   - *3_BRT_bootstrap.R*: using the best tuned parameters, calculate 1000 bootstrapped (with replacement) BRTs for each dependent variable.

## Notes
For the BRT, *bigmemory* and *foreach* packages were used to share the memory and work between cores and speed up the process.
Those might be memory consuming (32 GB of RAM were used for the script *3_BRT_bootstrap.R*).

## Citations
Oliva, A., Onana, E.V., Garner, R.E., Kraemer, S.A., Fradette, M., Walsh, D.A., Huot, Y., 2023. Geospatial analysis reveals a hotspot of fecal bacteria in Canadian prairie lakes linked to agricultural non-point sources. Water Res. *under review*.

Oliva, A., Garner, R.E., Walsh, D., Huot, Y., 2022. The occurrence of potentially pathogenic fungi and protists in Canadian lakes predicted using geomatics, in situ and satellite-derived variables: Towards a tele-epidemiological approach. Water Res. 209, 117935. https://doi.org/10.1016/j.watres.2021.117935
