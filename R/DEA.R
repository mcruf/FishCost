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

#~~~~~~~~~~~~~~~~~~~~~~~~~~

# For exemplification, we will use again use information from the Danish fisheries sampling programs
# higlighted in the paper "Benchmarking the cost-effectiveness of commercial fisheries and scientific survey 
# sampling programmes". 
# This includes data from (i) IBTS survey Q1, (ii) IBTS survey Q4, and (iii) on-board observers program. 
# Remember: only variable costs were considered!
# Also: Overall costs of each sampling strategy were converted to daily average costs
# (i.e., total average costs divided by the amount of days spent at sea).
# On-board observers spent 55 days at sea, while the research vessel in the 1st and 4th quarter surveys
# spent 19 and 18 days at sea, respectively. Lastly, to relate the average daily costs to the actual
# sampling size, it is assumed that the number of days spent at sea would also be proportionally affected 
# by the sampling size. 


#### Costs and No. of sampling days per monitoring program

#========================================================
# Scenarios        Days at Sea        Total costs (DKK)
#========================================================
#       On-board observers (daily cost = 4,922 DKK)     
#--------------------------------------------------------
#   100%              55                 270,710.00
#    75%              41                 203,032.50
#    50%              28                 135,355.00
#    25%              14                  67,677.50
#--------------------------------------------------------
#          Survey-Q1 (daily cost = 38,522 DKK)
#--------------------------------------------------------
#   100%              19                 731,909.67
#    75%              14                 539,301.86
#    50%              10                 385,215.62
#    25%               5                 192,607.81
#--------------------------------------------------------
#          Survey-Q4 (daily cost = 42,028 DKK)
#--------------------------------------------------------
#   100%              18                 756,497.94
#    75%              14                 588,387.29
#    50%               9                 378,248.97
#    25%               5                 210,138.32
#========================================================



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


rm(list=ls())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [1] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015



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
resfiles <- list.files(pattern="*.rds") 
reslist <- lapply(resfiles, readRDS)
results <- do.call("rbind",reslist); rm(reslist)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Results post-processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results[,c("Data","SamplingStrategy")] <- lapply(results[,c("Data","SamplingStrategy")],factor) #Set correct columns to factor


# 4.1) Remove outliers
#~~~~~~~~~~~~~~~~~~~~~~~~
# In this section we remove outstanding outliers following Tukey's fences method.
# This is the same method used in the original paper.
# However, note that this step is NOT MANDATORY.
# The outlier removal will differ from case study to case study.


########################    IGNORE THIS SECTION FOR THE PRESENT TOY EXAMPLE    ######################
# (Otherwise it filters out the sampling strategy that will be used later as a baseline in section 4.6)


# ## Calculate Tukey's fences by sampling strategy 
# results2 <- as.data.frame(results %>%
#                         # Compute per group mean values of columns
#                         group_by(SamplingStrategy) %>%
#                         mutate(iqr=IQR(TotAbu), Q1=summary(TotAbu)[2],Q3=summary(TotAbu)[5],
#                                lowerTukeyFence=summary(TotAbu)[2]-(1.5*IQR(TotAbu)),
#                                upperTukeyFence=summary(TotAbu)[5]+(1.5*IQR(TotAbu))) %>%
#                         ungroup())
# 
# ## Define which rows will be removed based on Tukey's fence criteria
# results2$Outlier <- ifelse(results2$TotAbu <= results2$lowerTukeyFence | results2$TotAbu >= results2$upperTukeyFence,"yes","no")  
# 
# ## Remove outliers
# results2 <- subset(results2,Outlier=="no")

results2 <- results #Ignore this line if not running script on toy example AND applying an outlier filtering method similar as the lines above



# 4.2) Calculate CV and variance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Info will be used later for the DEA
results2 <- as.data.frame(results2 %>%
                           # Compute per group mean values of columns
                           group_by(SamplingStrategy) %>%
                           mutate(CV=sd(TotAbu)/mean(TotAbu), Variance=var(TotAbu)) %>%
                           ungroup())


# 4.3) Get abundance performance metrics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This includes: mean abundance, median abundance, maximum abundance, CV and variance.
# We need to aggregate the data to get average values per sampling strategy
results2$SamplingStrategy <- factor(results2$SamplingStrategy)

dfsimu <- ddply(results2, .(SamplingStrategy), summarize,  
            meanTotAbu=mean(TotAbu), 
            medianTotAbu=median(TotAbu),
            maxAbu=max(TotAbu),
            CV=mean(CV),
            Variance=mean(Variance))


# 4.4) Include scenario-specific columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SCENARIO <- strsplit(as.character(dfsimu$SamplingStrategy),"_")
scenario_cols <- as.data.frame(do.call("rbind",SCENARIO))
colnames(scenario_cols) <- c("Scom","Ssurq1","Ssurq4")

dfsimu <- cbind(dfsimu, scenario_cols)
dfsimu[,c("Scom","Ssurq1","Ssurq4")] <- apply(dfsimu[,c("Scom","Ssurq1","Ssurq4")], 2, function(x) as.numeric(as.character(x)))



# 4.4) Include data input identifier column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dfsimu$DataInput <- ifelse(dfsimu$Scom == 0 & dfsimu$Ssurq1 == 0,"SQ4",
                     ifelse(dfsimu$Scom == 0 & dfsimu$Ssurq4 == 0,"SQ1",
                            ifelse(dfsimu$Ssurq1 > 0 & dfsimu$Ssurq4 > 0 & dfsimu$Scom == 0 ,"SQ1+SQ4",
                                   ifelse(dfsimu$Ssurq1 == 0 & dfsimu$Ssurq4 == 0,"COM",
                                          ifelse(dfsimu$Scom > 0 & dfsimu$Ssurq1 == 0 & dfsimu$Ssurq4 > 0,"COM+SQ4",
                                                 ifelse(dfsimu$Scom > 0 & dfsimu$Ssurq1 > 0 & dfsimu$Ssurq4 == 0,"COM+SQ1",
                                                        "COM+SQ1+SQ4"))))))

# Similar to the column above, but on the level of COM/SUR/BOTH
dfsimu$DataInput2 <- ifelse(dfsimu$DataInput == "SQ1","Sur",
                     ifelse(dfsimu$DataInput == "SQ4","Sur",
                           ifelse(dfsimu$DataInput == "SQ1+SQ4","Sur",
                                  ifelse(dfsimu$DataInput == "COM","Com",
                                         ifelse(dfsimu$DataInput == "COM+SQ1","Both",
                                                ifelse(dfsimu$DataInput == "COM+SQ4","Both", 
                                                       "Both"))))))

# Set scenario labels back to original labels (for convenience)
dfsimu$Scom <- scenario_cols$Scom
dfsimu$Ssurq1 <- scenario_cols$Ssurq1
dfsimu$Ssurq4 <- scenario_cols$Ssurq4



# 4.5) Include the costs realted to each sampling strategy 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For details, see costs reported at the start of this script
## (alternatively, costs reported in the paper)

dfsimu$cost_COM <- NA
dfsimu$cost_SQ1<- NA
dfsimu$cost_SQ4<- NA

# Costs of the commercial scenarios
dfsimu$cost_COM <- ifelse(dfsimu$Scom=="0",0,
                      ifelse(dfsimu$Scom=="025",67677.5,
                             ifelse(dfsimu$Scom=="05",135355,
                                    ifelse(dfsimu$Scom=="075",203032.5, 270710))))

# Costs of the survey-Q1 scenarios
dfsimu$cost_SQ1 <- ifelse(dfsimu$Ssurq1=="0",0,
                      ifelse(dfsimu$Ssurq1=="025",192607.81,
                             ifelse(dfsimu$Ssurq1=="05",385215.62,
                                    ifelse(dfsimu$Ssurq1=="075",539301.86, 731909.67))))

# Costs of the survey-Q4 scenarios
dfsimu$cost_SQ4 <- ifelse(dfsimu$Ssurq4=="0",0,
                      ifelse(dfsimu$Ssurq4=="025",210138.32,
                             ifelse(dfsimu$Ssurq4=="05",378248.97,
                                    ifelse(dfsimu$Ssurq4=="075",588387.29, 756497.94))))

# Calculate total costs for a given sampling strategy
dfsimu$totCost <- rowSums(dfsimu[12:14])



# 4.6) Calculate median and maximum bias
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The bias will be calculated by taking the mean abundance of a given sampling strategy
# substracted from the total abundance of the baseline sampling strategy.
# Note: for the paper the assumed baseline was assumed as the sampling stragey containing
# the full amount of data from both surveys (i.e, 0_1_1).
# However, the baseline is likely to change from case study to case study, depending on the research question!

baseline <- subset(dfsimu,SamplingStrategy == "0_1_1")$meanTotAbu #Baseline abundance


dfsimu$biasMedian <- abs(dfsimu$medianTotAbu-baseline) #bias in relation to the median of the total abundance of sampling strategy
dfsimu$biasMax <- abs(dfsimu$maxAbu-baseline) #bias in relation to the maximum value of the total abundance of each sampling strategy


# 4.7) Remove non-simuluated sampling strategies 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are all sampling strategies containing scenarios with 100% of data selection.
# Necessary step, otherwise it will create problems for the DEA
# (no variances and CVs can be calculated from these sampling strategies)

#     - (Com=0, SurQ1=1, SurQ4=1)
#     - (Com=0, SurQ1=0, SurQ4=1)
#     - (Com=0, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=0)
#     - (Com=1, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=1)

out <- c("0_1_1", "0_0_1", "0_1_0", "1_0_0", "1_1_0", "1_0_1", "1_1_1")

dfsimu <- subset(dfsimu, !(SamplingStrategy %in% out))



#~~~~~~~~~~~~~~~~~~~~
# 5) Go for the DEA 
#~~~~~~~~~~~~~~~~~~~~

# First we need to remove all sampling strategies containing scenarios with 100% of 



# 5.1) The most cost-effective sampling strategies (DEA case-1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identifies the most cost-effective sampling strategies
# Input: bias (mean) and variance 
# output: total costs


## Prepare input and output matrices
input1 <- as.matrix(dfsimu[,c("biasMedian","Variance")]) 
rownames(input1) <- dfsimu$SamplingStrategy


output1 <- as.matrix(dfsimu[,"totCost"])
rownames(output1) <- dfsimu$SamplingStrategy
class(output1)  <- "numeric"


## Run DEA
DEA1 <- dea(input1,output1, RTS="vrs", ORIENTATION="in") #VRS method


## Get raw efficiecny scores 
EFF1 <- eff(DEA1) 
sort(round(EFF1,2))


## Bootstrap efficiency scores to get uncertainty
EFF1boot <- dea.boot(input1,output1, RTS="vrs", ORIENTATION="in", NREP = 500) #VRS method

DEA1df <- data.frame(EFF_score=EFF1boot$eff.bc, SamplingStrategy=names(EFF1boot$eff),
                   CI_lower=EFF1boot$conf.int[,2], CI_upper=EFF1boot$conf.int[,1])




# 5.2) The most risk-averse sampling strategies (DEA case-2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identifies the most risk-averse sampling strategies
# Input: bias (maximum) and variance 
# output: total cost


## Prepare input and output matrices
input2 <- as.matrix(dfsimu[,c("biasMax","Variance")]) 
rownames(input2) <- dfsimu$SamplingStrategy


output2 <- as.matrix(dfsimu[,"totCost"])
rownames(output2) <- dfsimu$SamplingStrategy
class(output2)  <- "numeric"


## Run DEA
DEA2 <- dea(input2,output2, RTS="vrs", ORIENTATION="in") #VRS method


## Get raw efficiecny scores 
EFF2 <- eff(DEA2) 
sort(round(EFF2,2))


## Bootstrap efficiency scores to get uncertainty
EFF2boot <- dea.boot(input2,output2, RTS="vrs", ORIENTATION="in", NREP = 500) #VRS method

DEA2df <- data.frame(EFF_score=EFF2boot$eff.bc, SamplingStrategy=names(EFF2boot$eff),
                     CI_lower=EFF2boot$conf.int[,2], CI_upper=EFF2boot$conf.int[,1])




# 5.3) Sampling strategies with the best overall trade-off (DEA case-2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This approach considers as an input measure the cost, variance and bias (either max or median)
# and constant of 1 as an output. This creats a 3D frontier, where everyone has the same starting point
# at 1 and aims to work-out what´s the best combination across that frontier, and how much of each 
# sampling strategy could be improved. I.e., among the 3D surface, bias, variance and cost are minimized,
# returining the best set of combination.


## Prepare input and output matrices
input3 <- as.matrix(dfsimu[,c("biasMedian","Variance","totCost")])
rownames(input3) <- dfsimu$SamplingStrategy


output3 <- as.matrix(rep(1, nrow(dfsimu)))
rownames(output3) <- dfsimu$SamplingStrategy


## Run DEA
DEA3 <- dea(input3,output3, RTS="vrs", ORIENTATION="in") #VRS method


## Get raw efficiecny scores 
EFF3 <- eff(DEA3) 
sort(round(EFF3,2))


## Bootstrap efficiency scores to get uncertainty
EFF3boot <- dea.boot(input3,output3, RTS="vrs", ORIENTATION="in", NREP = 500) #VRS method

DEA3df <- data.frame(EFF_score=EFF3boot$eff.bc, SamplingStrategy=names(EFF3boot$eff),
                     CI_lower=EFF3boot$conf.int[,2], CI_upper=EFF3boot$conf.int[,1])



#~~~~~~~~~~~~
# 6) Plots
#~~~~~~~~~~~~

# Plot efficiency scores + confidence intervals 
# For illustration, the plot below is only made for DEA case-1.

## Copy info related to input data
DEA1df$DataInput <- dfsimu$DataInput

## Reorder levels of Sampling strategies (for visualization purpose)
## A tedious step, though....
DEA1df$SamplingStrategy <- as.factor(DEA1df$SamplingStrategy )
DEA1df$SamplingStrategy <- ordered(DEA1df$SamplingStrategy, levels=c(
                                    "025_0_0", "05_0_0", "075_0_0", #commercial
  
                                    "0_025_0","0_05_0","0_075_0", #survey SQ1
                                    "0_0_025", "0_0_05", "0_0_075", #survey SQ4
                                    "0_025_025","0_025_05","0_025_075","0_025_1", #survey full
                                    "0_05_025","0_05_05","0_05_075","0_05_1", #survey full
                                    "0_075_025","0_075_05","0_075_075","0_075_1", #survey full
                                    "0_1_025","0_1_05","0_1_075", #survey full
                                    
                                    "025_025_0","025_05_0","025_075_0","025_1_0", #both
                                    "05_025_0","05_05_0","05_075_0","05_1_0", #both
                                    "075_025_0","075_05_0","075_075_0","075_1_0", #both
                                    "1_025_0","1_05_0","1_075_0", #both
                                    "025_0_025", "025_0_05", "025_0_075","025_0_1", #both
                                    "05_0_025", "05_0_05", "05_0_075","05_0_1", #both
                                    "075_0_025", "075_0_05", "075_0_075","075_0_1", #both
                                    "1_0_025", "1_0_05", "1_0_075", #both
                                    "025_025_025", "025_025_05", "025_025_075","025_025_1", #both
                                    "025_05_025", "025_05_05", "025_05_075","025_05_1", #both
                                    "025_075_025", "025_075_05", "025_075_075","025_075_1", #both
                                    "025_1_025", "025_1_05", "025_1_075","025_1_1", #both
                                    "05_025_025", "05_025_05", "05_025_075","05_025_1", #both
                                    "05_05_025", "05_05_05", "05_05_075","05_05_1", #both
                                    "05_075_025", "05_075_05", "05_075_075","05_075_1", #both
                                    "05_1_025", "05_1_05", "05_1_075","05_1_1", #both
                                    "075_025_025", "075_025_05", "075_025_075","075_025_1", #both
                                    "075_05_025", "075_05_05", "075_05_075","075_05_1", #both
                                    "075_075_025", "075_075_05", "075_075_075","075_075_1", #both
                                    "075_1_025", "075_1_05", "075_1_075","075_1_1", #both
                                    "1_025_025", "1_025_05", "1_025_075","1_025_1", #both
                                    "1_05_025", "1_05_05", "1_05_075","1_05_1", #both
                                    "1_075_025", "1_075_05", "1_075_075","1_075_1",#both
                                    "1_1_025", "1_1_05", "1_1_075"))
  

## Replace underlines by colon (only for visualization purposes)
DEA1df$SamplingStrategy <- as.factor(gsub("_",":",DEA1df$SamplingStrategy))


## Reorder also levels of input data ID
DEA1df$DataInput <- as.factor(DEA1df$DataInput)
DEA1df$DataInput <- ordered(DEA1df$DataInput , levels=c("COM","SQ1","SQ4","SQ1+SQ4","COM+SQ1","COM+SQ4","COM+SQ1+SQ4"))



jco <- c("#A73030FF","#4A6990FF","#7AA6DCFF","#B09C85FF","#79AF97FF","#DF8F44FF","#868686FF")
ggplot(DEA1df,aes(x=SamplingStrategy, y=EFF_score,col=DataInput)) + 
  geom_pointrange(aes(ymin=CI_lower, ymax=CI_upper),size=.75,shape=16,position=position_dodge(0.7)) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.5, position=position_dodge(0.7)) +
  scale_color_manual(values=jco) +
  theme_pubclean() +
  ylab("Eff scores") + xlab("Sampling strategies") +
  theme(#legend.position="none",
    plot.margin = unit(c(1,1,1,1),"cm"),
    axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1.2, size=11, margin=margin(15,0,0,0)),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14,margin = margin(t = 0, r = 20, b = 0, l = 15)),
    axis.title.x = element_text(size=14,vjust=-4,margin = margin(t = 0, r = 20, b = 20, l = 0)),
    legend.position = "right",
    legend.title = element_blank())

