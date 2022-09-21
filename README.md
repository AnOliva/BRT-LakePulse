# BRT-LakePulse

Author: Anaïs Oliva

Correspondence: anais.oliva@usherbrooke.ca

## Overview - Usage
This repository contains scripts used in Oliva et al. (submitted and under review) to: (1) explore the diversity of potentially pathogenic bacteria (PPB); (2) build a multi-indicator model of anthropogenic fecal contamination from a cluster of PPB; and (3) predict the bacterial multi-indicator over thousands of lakes.

This study is based on the abundance of 16S rRNA amplicon sequences of PPB sampled in 413 lakes within 8 southern Canadian ecozones and representing a wide diversity of lakes and watershed land use

The folders are devided as follow:

1. *1_microbes_preprocessing*
   - *1_rarefaction.R*: rarefaction of sampled 16S rRNA sequences.
   - *2_pathogen_extraction.R* : extraction of the potential pathogens based on the ePathogen datase and a partial match algorithm.
   - *3_clustering_analysis.R*: clustering analysis to form groups of similar pathogens based on Sorensen dissimilarity index.
   - 
3. *2_BRT_preprocessing*
   - *1_combined_variables.R*: load all the independent variables, verify the flags and normalize the data.
   - *2_collinearity_CDOM-TSS.R*: check Pearson correlation and the generalized variance inflation factor (GVIF) for TSS and CDOM datasets.
   - *2_collinearity_QC.R*: check Pearson correlation and the GVIF for the quality controlled dataset .
   - *2_collinearity_upscaled.R*: check Pearson correlation and the GVIF for the upscaled dataset.
   - 
5. *BRT*
   - *1_normalization.R*: set/edit the destination directory for BRT results, loads required datasets and make preparation for the BRT scripts.
   - *2_BRT_preliminary.R*: calculate 5 consecutive BRT for each combination of BRT parameters defined in an hypergrid. Select the best combination.
   - *3_BRT_bootstrap.R*: using the best tuned parameters, calculate 1000 bootstrapped (with replacement) BRTs for each dependent variable.

## Citations
Oliva, A., Onana, E.V., Garner, R.E., Kraemer, S.A., Fradette, M., Walsh, D.A., Huot, Y., 2022. A multi-indicator mapping analysis reveals a spatial hotspot of a putative bacterial pathogens assemblage in Canadian prairie lakes that is linked to anthropogenically-altered landscapes. Water Res. *submitted*.

Oliva, A., Garner, R.E., Walsh, D., Huot, Y., 2022. The occurrence of potentially pathogenic fungi and protists in Canadian lakes predicted using geomatics, in situ and satellite-derived variables: Towards a tele-epidemiological approach. Water Res. 209, 117935. https://doi.org/10.1016/j.watres.2021.117935
