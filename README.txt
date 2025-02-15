#################
DATA PACKAGE FROM
#################

"Plant and animal functional diversity drive mutualistic network assembly across an elevational gradient"

by Albrecht et al. Nature Communications 9:3177 (2018)

https://doi.org/10.1038/s41467-018-05610-w

###########
DESCRIPTION
###########

This data package contains all data files and computer code (in R language) that is needed to reproduce the analyses presented in the paper.
In particular, the data package contains the following objects (Note: asterisks [*] denote existing files, arrows [->] denote output files generated during analysis):

* 'runScript.R' (file) - This is the central file that needs to be executed to reproduce the analysis from the paper. The file will execute the scripts in the folder 'include', use the data in the folder 'data' and store results of the analysis in the folder 'output'.


* 'data' (folder) - A folder containing the following data objects:
  *  'dataSet1.RData' (file) - RData-file that contains the raw data for the analysis in R.
  *  'rawData' (folder) - Contains the raw data files in .csv-format.
  *  'metaData' (folder) - Contains the meta data for each file, containing descriptions of all variables and how they were measured.
  -> 'jagsData1.RData' (file) - RData-file containing the prepared data (produced by script '01_prepareData.R').
  -> 'jagsModelFit1.RData' (file) - RData-file containing the results of the Bayesian analysis (produced by script '03_fitModel.R')


* 'include' (folder) - A folder containing the following R-scripts to reproduce the analysis.
  *  '00_sourceFunctions.R' (file) - R-script containing all functions written for the analysis.
  *  '01_prepareData.R' (file) - R-script to prepare & aggregate the data for the analysis (NOTE: this script also executes the RLQ analysis).
  *  '02_setupModel.R' (file) - R-script that contains the JAGS code for the Bayesian hierarchical structural equation model.
  *  '03_fitModel.R' (file) - R-script to fit the JAGS model.
  *  '04_inferenceModel.R' (file) - R-script to extract, plot and summarize the results from the RLQ analysis and the structural equation model.
  -> 'jagsModel1.txt' (file) - File in .txt-format containing the JAGS model code (produced by script '02_setupModel.R').


* 'output' (folder) - A folder in which all results of the analysis after execution of the file 'runScript.R' are stored (see above).
  -> '01_descriptive_summary.txt' (file) - A file in .txt-format that contains all main results that are reported in the paper.
  -> '02_rlqTable.csv' (file) - A table (in .csv-format) containing the correlations of individual plant and animal traits with the first and second RLQ-axes with fourth-corner statistic (see paper).
  -> '03_rlqMeanTable.csv' (file) - A table (in .csv-format) containing the correlations of plant and animal traits aggregated by trait type RLQ-axis with Moran's Test statistic (see paper).
  -> '04_pathTable1.csv' (file) - A table (in .csv-format) containing information for path coefficients of the structural equation model for niche breadth of plants and animals
  -> '05_pathTable2.csv' (file) - A table (in .csv-format) containing information for path coefficients of the structural equation model for niche partitioning of plants and animals
  -> 'Figures 2 & 3' & 'Supplementary Figures 1-4' from the paper (produced by script '04_inferenceModel.R').

	
#######################
HOW TO RUN THE ANALYSIS
#######################

To run the analysis you simply need to set your working directory in R to ".../data_package/" and run the script file named "runScript.R" in the main folder of the package.

############
SESSION INFO
############

R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS  10.13.5

locale:
[1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] qgraph_1.4.3         rjags_4-6            coda_0.19-1          bipartite_2.08       sna_2.4              network_1.13.0      
 [7] statnet.common_4.0.0 FD_1.0-12            geometry_0.3-6       magic_1.5-6          abind_1.4-5          ape_4.1             
[13] ade4_1.7-8           vegan_2.4-4          lattice_0.20-35      permute_0.9-4       

loaded via a namespace (and not attached):
 [1] maps_3.2.0          splines_3.3.2       ellipse_0.3-8       gtools_3.5.0        Formula_1.2-2       stats4_3.3.2       
 [7] latticeExtra_0.6-28 d3Network_0.5.2.1   pbivnorm_0.6.0      backports_1.1.0     quadprog_1.5-5      digest_0.6.12      
[13] RColorBrewer_1.1-2  checkmate_1.8.3     ggm_2.3             minqa_1.2.4         colorspace_1.3-2    htmltools_0.3.6    
[19] Matrix_1.2-11       plyr_1.8.4          psych_1.7.3.21      pkgconfig_2.0.1     corpcor_1.6.9       scales_0.5.0       
[25] whisker_0.3-2       glasso_1.8          jpeg_0.1-8          fdrtool_1.2.15      lme4_1.1-15         huge_1.2.7         
[31] arm_1.9-3           tibble_1.3.4        htmlTable_1.9       mgcv_1.8-19         ggplot2_2.2.1       nnet_7.3-12        
[37] lazyeval_0.2.0      mnormt_1.5-5        survival_2.41-3     magrittr_1.5        nlme_3.1-131        MASS_7.3-47        
[43] foreign_0.8-69      data.table_1.10.4   tools_3.3.2         stringr_1.2.0       munsell_0.4.3       cluster_2.0.6      
[49] sem_3.1-9           rlang_0.1.2         grid_3.3.2          nloptr_1.0.4        rjson_0.2.15        htmlwidgets_0.9    
[55] spam_1.4-0          igraph_1.1.2        lavaan_0.5-23.1097  base64enc_0.1-3     boot_1.3-20         mi_1.0             
[61] gtable_0.2.0        reshape2_1.4.2      gridExtra_2.2.1     knitr_1.17          Hmisc_4.0-3         stringi_1.1.5      
[67] matrixcalc_1.0-3    Rcpp_0.12.12        fields_9.0          rpart_4.1-11        acepack_1.4.1       png_0.1-7 