##################################################################################
#                                                                                #
#                                                                                #
#                       Cost-effectiveness analysis - Plaice                     #
#                                 (PhD project II)                               #
#                                                                                #
#                                                                                #
##################################################################################


# In this script we perform the cost-effectiveness analysis.
# This is basically a linear model, whose response variable is variance/cost, and
# predictor variables are the costs associated to each scenario and response variable.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load data & additional data processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.1) Read data
#~~~~~~~~~~~~~~~~~~
#setwd("D:/PhD/Project II/2016/Plaice/Finalized/")
setwd("D:/PhD/Project II/2015/Plaice/Finalized/")



temp <- list.files(pattern=".rds")
temp2 <- lapply(temp, readRDS)
results <- do.call("rbind",temp2)



# 1.2) Remove outliers
#~~~~~~~~~~~~~~~~~~~~~~~~
# Same procedure as done in the "plotting_result" script;
# Too big/small outliers are driving a lot the final results; hence, needs to be removed.


## For 2016 data!!!
###################
# idxout <- which(results$cv>0.3)
# idxout2 <- which(results$totAbund >= 8000)
# results <- results[-c(idxout,idxout2),]



## For 2015 data!!!
###################
idxout <- which(results$totAbund >= 7000)
results <- results[-idxout,]



#The CV and VARIANCE needs to be calculated again, as outliers were removed
results <- as.data.frame(results %>%
                           # Compute per group mean values of columns
                           group_by(scenario) %>%
                           mutate(CV2=sd(totAbund)/mean(totAbund), variance2=var(totAbund)) %>%
                           ungroup())



# 1.3) Data processing
#~~~~~~~~~~~~~~~~~~~~~~~

## Rename the dataID levels
results$dataID <- as.factor(results$dataID)
#levels(results$dataID) <- c("both","both","com","sur") #Extra care here!!!! Double check all levels
levels(results$dataID) <- c("both","both","com","com","sur","sur") #2015



## Include scenario-specific column
SCENARIO <- strsplit(results$scenario,"_")
scenario_cols <- as.data.frame(do.call("rbind",SCENARIO))
colnames(scenario_cols) <- c("Scom","Ssurq1","Ssurq4")

results <- cbind(results, scenario_cols)

## Include total numbers of hauls for each scenario
results$NH_COM <- 0 # No. of hauls of commercial data
results$NH_SQ1 <- 0 # No. of hauls of survey Q1 data
results$NH_SQ4 <- 0 # No. of hauls of survey Q4 data



## FOR 2016 DATA!!!
#####################
# results$NH_COM <- ifelse(results$Scom=="0",0,
#                          ifelse(results$Scom=="025",22,
#                                 ifelse(results$Scom=="05",45,
#                                        ifelse(results$Scom=="075",68, 90))))
# 
# results$NH_SQ1 <- ifelse(results$Ssurq1=="0",0,
#                          ifelse(results$Ssurq1=="025",27,
#                                 ifelse(results$Ssurq1=="05",54,
#                                        ifelse(results$Ssurq1=="075",80,
#                                               107))))
# 
# results$NH_SQ4 <- ifelse(results$Ssurq4=="0",0,
#                          ifelse(results$Ssurq4=="025",28,
#                                 ifelse(results$Ssurq4=="05",55,
#                                        ifelse(results$Ssurq4=="075",83,
#                                               110))))



## FOR 2015 DATA!!!
#####################
results$NH_COM <- ifelse(results$Scom=="0",0,
                         ifelse(results$Scom=="025",9,
                                ifelse(results$Scom=="05", 18,
                                       ifelse(results$Scom=="075", 26,
                                              35))))

results$NH_SQ1 <- ifelse(results$Ssurq1=="0",0,
                         ifelse(results$Ssurq1=="025",12,
                                ifelse(results$Ssurq1=="05",24,
                                       ifelse(results$Ssurq1=="075",35,
                                              47))))

results$NH_SQ4 <- ifelse(results$Ssurq4=="0",0,
                         ifelse(results$Ssurq4=="025",13,
                                ifelse(results$Ssurq4=="05",26,
                                       ifelse(results$Ssurq4=="075",40,
                                              53))))




## Incude detailed scenario information 
results$dataID2 <- ifelse(results$NH_COM == 0 & results$NH_SQ1 == 0,"SQ4",
                          ifelse(results$NH_COM == 0 & results$NH_SQ4 == 0,"SQ1",
                                 ifelse(results$NH_SQ1 > 0 & results$NH_SQ4 > 0 & results$NH_COM == 0 ,"SQ1+SQ4",
                                        ifelse(results$NH_SQ1 == 0 & results$NH_SQ4 == 0,"COM",
                                               ifelse(results$NH_COM > 0 & results$NH_SQ1 == 0 & results$NH_SQ4 > 0,"COM+SQ4",
                                                      ifelse(results$NH_COM > 0 & results$NH_SQ1 > 0 & results$NH_SQ4 == 0,"COM+SQ1",
                                                             "COM+SQ1+SQ4"))))))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Include the scenario-specific costs 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Commercial scenarios
#~~~~~~~~~~~~~~~~~~~~~~
# COM_1 = 270710.00
# COM_075 = 203032.50
# COM_05 = 135355.00
# COM_025 = 67677.50
# COM_0 = 0

# Survey Q1 scenarios
#~~~~~~~~~~~~~~~~~~~~~~
# SURQ1_1 = 731909.67
# SURQ1_075 = 539301.86
# SURQ1_05 = 385215.62
# SURQ1_025 = 192607.81
# SURQ1_0 = 0


# Survey Q4 scenarios
#~~~~~~~~~~~~~~~~~~~~~~
# SURQ4_1 = 756497.94
# SURQ4_075 = 588387.29
# SURQ4_05 = 378248.97
# SURQ4_025 = 210138.32
# SURQ4_0 = 0


results$cost_COM <- NA
results$cost_SQ1<- NA
results$cost_SQ4<- NA

results$cost_COM <- ifelse(results$Scom=="0",0,
                           ifelse(results$Scom=="025",67677.5,
                                  ifelse(results$Scom=="05",135355,
                                         ifelse(results$Scom=="075",203032.5, 270710))))

results$cost_SQ1 <- ifelse(results$Ssurq1=="0",0,
                           ifelse(results$Ssurq1=="025",192607.81,
                                  ifelse(results$Ssurq1=="05",385215.62,
                                         ifelse(results$Ssurq1=="075",539301.86, 731909.67))))

results$cost_SQ4 <- ifelse(results$Ssurq4=="0",0,
                           ifelse(results$Ssurq4=="025",210138.32,
                                  ifelse(results$Ssurq4=="05",378248.97,
                                         ifelse(results$Ssurq4=="075",588387.29, 756497.94))))

results$totCost <- rowSums(results[19:21])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Calculate the response variable (variance/cost)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#results$response <- results$totCost/results$variance 
results$response <- results$totCost/results$variance2 


df <- results[!duplicated(results$scenario),]

df$response <- round(df$response,2)


library(ggplot2)
library(ggpubr)

ggplot(df, aes(variance2, totCost)) +
  geom_point() + geom_smooth(method = "loess") +
  theme_pubr() +
  scale_y_continuous(labels = scales::scientific) +
  labs(x="Variance",y="Cost (DKK) ") +
  theme(
    axis.text.x = element_text(face="bold",size=11),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.y = element_text(margin=margin(t=0,r=20,b=0,l=0),size=12,face="bold"),
    axis.title.x = element_text(margin=margin(t=20,r=0,b=0,l=0),size=12,face="bold"),
    plot.margin = unit(c(1,1,1,1),"cm"))



#~~~~~~~~~~~~~~~~~~~~~~~
# 4) Go for the LM
#~~~~~~~~~~~~~~~~~~~~~~~

model <- lm(response ~ Scom + Ssurq1 + Ssurq4, data=df)
summary(model)
par(mfrow=c(2,2))
plot(model)


# response on log-scale
model2 <- lm(log1p(response) ~ Scom + Ssurq1 + Ssurq4, data=df)
summary(model2)
par(mfrow=c(2,2))
plot(model2)


library(ggiraphExtra)
ggPredict(model,interactive = TRUE)
