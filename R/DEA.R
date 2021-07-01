##########################################################################################
#                                                                                        #
##              Benchmarking the cost-effectiveness of commercial fisheries             ##
##                       and scientific survey sampling programmes                      ##
##                     Script No. 3: Data Envelopment Analysis (DEA)                     #
#                                                                                        #
##########################################################################################


# Last update: June, 2021

# Code written and mantained by Marie-Christine Rufener
# Contact < macrufener@gmail.com > for any query or to report code issues.


# This script performs the actual cost-effectiveness analysis.
# We will use a DEA model from the Benchmarking R-package.
# The following information will be needed:
# (i) Costs of each monitoring program
# (ii) No. of days that each monitoring program spent at sea collecting samples 
# (iii) Simulated abundances 


# For exemplification, we will use again use information from the Danish fisheries sampling programs
# higlighted in the paper "Benchmarking the cost-effectiveness of commercial fisheries and scientific survey sampling programmes".
# This includes data from (i) IBTS survey (Q1 and Q4) and (ii) on-board observers program. 


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


rm(list=ls())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [1] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [1] #Choose data from which results should be loaded; default is both


#~~~~~~~~~~~~~~~~~~~~~~
# 2) Load R packages
#~~~~~~~~~~~~~~~~~~~~~~
library(Benchmarking)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Load results from simulated abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are the results that were processed in the previous script (AbuProcessing.R)


# 3.1) Define WD 
#~~~~~~~~~~~~~~~~~
# Based on default inputs
setwd(paste("~","FishCost", "Results", "SimuAbundance", SPECIES , YEAR, "Processed",sep="/"))


# 3.2) Load results
#~~~~~~~~~~~~~~~~~~~~
results <- readRDS(paste0("Processed_SimuResults_", DATA, ".rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Result post-processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results[,c("Data","Scenario")] <- lapply(results[,c("Data","Scenario")],factor) #Set correct columns to factor


# 4.1) Remove outliers
#~~~~~~~~~~~~~~~~~~~~~~~~
# In this section we remove outstanding outliers following Tukey's fences method.
# This is the same method used in the original paper.
# However, note that this step is NOT MANDATORY.
# The outlier removal will differ from case study to case study.

## Caclulate Tukey's fences by sampling strategy 
results2 <- as.data.frame(results %>%
                        # Compute per group mean values of columns
                        group_by(Scenario) %>%
                        mutate(iqr=IQR(TotAbu), Q1=summary(TotAbu)[2],Q3=summary(TotAbu)[5],
                               lowerTukeyFence=summary(TotAbu)[2]-(1.5*IQR(TotAbu)),
                               upperTukeyFence=summary(TotAbu)[5]+(1.5*IQR(TotAbu))) %>%
                        ungroup())

## Define which rows will be removed based on Tukey's fence criteria
results2$remove <- ifelse(results2$TotAbu <= results2$lowerTukeyFence | results2$TotAbu >= results2$upperTukeyFence,"remove","keep")  

## Remove outliers
results2 <- subset(results2,remove=="keep")



# 4.2) Calculate CV and variance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Info will be used later for the DEA

results2 <- as.data.frame(results2 %>%
                           # Compute per group mean values of columns
                           group_by(Scenario) %>%
                           mutate(CV=sd(TotAbu)/mean(TotAbu), Variance=var(TotAbu)) %>%
                           ungroup())
