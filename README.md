# Structural_coexistence_multifunctionality

In this folder you may find the data and code for the analysis of "Stably coexisting communities deliver higher ecosystem multifunctionality".

## What can you find in this project?

### a) data

Please refer to the METADATA text file that regroups all information on the data folder, and contains a list of abbreviations for each variable that can be found in this folder.

#### *coexistence and experimental design*

In the data folder, you can find a subfolder called "species alphas and lambdas" containing all competition matrices used during the analysis, as text files which names start with "biommatrix", and their associated standard errors, which are named "SE_biommatrix". The species intrinsic growth rates are also present in this folder under the name starting with "biomintrinsic".These matrices are derived over the whole year 2020 (combining June and August sampling for this year), or only for June 2020 sampling, as this month was used to compute predictions for 2022 species evenness.

Two text files, "plot_info" and "structural_coexistence_experimental_design", describe the organisation of the plots in our study, which species compose them and what are their predicted coexistence strength and mechanisms. The "structural_coexistence_experimental_design" file itself is derived from another text file "all_triplets_experiment" which lists structural coexistence values for all possible triplets from our pool of species. The environment image "structural_coexistence_image.Rdata" also contains the same information but pre-loaded for convenient sourcing.

#### *ecosystem functions*

A set of ecosystem functions and data sampled in our triplet experiment may be found in this data folder. Each text file name starts with the month and date of its sampling. The functions measured are:

-   biomass in May and August 2022

-   betaglucosidase, in May and August 2022

-   phosphatase in May and August 2022

-   herbivory damage in May 2022

-   pathogen damage in May 2022

-   root biomass in August 2022

-   decomposition in April 2023

Aside from these ecosystem function, the percentage cover of each species present in the plot was also recorded in June and August 2022.

### b) R functions

This repository contains two folders that regroups functions from previous studies:

-   the folder "anisoFun functions" contains functions from the paper: "Structural asymmetry in biotic interactions as a tool to understand and predict ecological persistence.", by Allen‐Perkins, Alfonso, et al., *Ecology Letters* (2023). This code was used to compute the minimum distance to exclusion variable, a measure of coexistence strength.

-   the folder "structural coexistence functions" contains a set of functions from the paper: "A structural approach for understanding multispecies coexistence." by Saavedra, Serguei, et al., *Ecological Monographs* (2017). This code was used to compute all measure of coexistence mechanisms: niche and fitness differences as well as indirect interactions.

### c) code

#### *main analyses*

This project contains a set of original code. These code files are numbered from (1) to (5) and they contain the main analyses of the study. Each file can be run independently, as they all source all files that are needed to run the analysis. However, the code was build in this specific order and the logic of the analysis might be easier to follow if each file is read in this order.

#### *extra code*

Finally, the repository contains some files that are satellites to the main analyses. These files don't need to be read for the main analysis to function, as their output is already saved in the "data" folder.

-   The first one, "structural_coexistence_control_nitrogen.qmd" is a quarto document that I used to compute all coexistence mechanisms and minimum distance to exclusion for the set of 48 triplets in this experiment. The format of the code file is different because this was a file created before the start of the experiment, in 2021.

-   The second one, "time_to_extinction_tool_no_SE" was authored by Oscar Godoy and modified by me, in order to compute the predicted evenness and abundances for our communities. Its output are the files "bh_predicted_evenness_15y_2022" and "abundances_predictedBH_all_plots" and they are analysed in the main code file number (5).

## Contact information

If you wish to get more information on the code or if you have questions, feel free to contact me using this email address:

*caroline.daniel\@unibe.ch*
