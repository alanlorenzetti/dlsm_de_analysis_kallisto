# alorenzetti 202101

# description ####
# this script is the scaffold
# running all the scripts
# required to perform
# the differential expression analysis

# loading required variables ####
# setting up number of threads
threads=8

# should we run quality control?
reportFlag="no"

# DE analysis thresholds
# DESeq2 adjusted pval
padjthreshold = 0.05

# DeSeq2 log2FoldChange
log2fcthreshold = 0

# creating data directory
if(!dir.exists("data")){dir.create("data")}

# creating results directory
if(!dir.exists("results")){dir.create("results")}

# sourcing ####
# loading libs
source("scripts/01_loadingLibs.R")
if(reportFlag == "yes"){source("scripts/02_qualityControl.R")}
source("scripts/03_deAnalysis.R")
source("scripts/04_results.R")
