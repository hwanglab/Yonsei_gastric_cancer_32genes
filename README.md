**Overview**

This zip file contains MATLAB and R scripts that reproduce the main results of the paper. 

**System requirements**

The codes require only a standard computer with enough RAM to support the in-memory operations.

**OS requirements**

The codes are tested on Windows 10 and Linux (Ubuntu 20.04)

**Installation guide** 

All the scripts are tested on MATLAB (R2018-b, R2021a) and R (4.0.4). It might take few minutes to complete to install all the required packages.  

The MATLAB scripts require few functions, nmf and LIBSVM packages. The additional functions and the nmf package which can be used without compilation have been included in “utils” directory and included in the scripts. You need to install LIBSVM package in MATLAB. Please follow the instructions below. 
Installation of LIBSVM
-	Windows: pre-compiled mex files are included in "util" directory. The path where the package is located has been include in the scripts  
-	Linux: The scripts interact with LIBSVM in MATLAB interface. Please visit the Github page for LIBSVM (https://github.com/cjlin1/libsvm) and 
		   check out the installation instructions under “MATLAB/OCTAVE Interface” (After downloading the package, please go to "matlab" directory and 
		   run "make.m" in MATLAB. Please make sure that C/C++ compilers are properly installed in your MATLAB. 
		   Please add "matlab" directory as the path for the package in MATLAB after "make.m" runs successfully.)  

The R script requires survminer, survival, ggplot2 and gdata packages. You can use install.packages("package name") to install the packages. 

**Demo**

Please find each script file in "codes" directory and run it in MATLAB or R. The datasets which are required to run each script have been included in "data" directory. 

For each script, the expected output and run time are as follows

script_consensus_clustering.m 
- Kaplan-meier (KM) curves for overall survival stratified by the consensus clustering described in the main text (please see Figure 2-BCD)
- It might take few minutes due to the bootstrapping steps (where NMF runs 1000 times)

script_riskscore_prediction.m
- KM curves for overall survival stratified by risk group in ACRG + MD Anderson + TCGA combined dataset (please see Figure 3-B)
- It might take less than a minutes 

script_consensus_clustering.m 
- KM curves for overall survival stratified by the multiclass classifier trained using the Yonsei cohort 
- It might take less than a minutes 

script_5Fu_comparison.R
- Adjust KM curves for overall survival in the Yonsei cohort with the adjuvant chemotherapy information (please see Figure 4)


 
