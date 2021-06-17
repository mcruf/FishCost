##########################################################################################
#                                                                                        #
##              Benchmarking the cost-effectiveness of commercial fisheries             ##
##                       and scientific survey sampling programmes                      ##
##                    Script No. 2: Processing of simulated abundances                   #
#                                                                                        #
##########################################################################################

# Last update: June, 2021


# Code written and mantained by Marie-Christine Rufener
# Contact < macrufener@gmail.com > for any query or to report code issues.


# This script is to be performed after the abundance simulations (Script No. 1 - SimuAbu)
# The LGNB model provides the abundance estimates for each grid point and input time-period.
# As the model was set on a quarterly level, this means that we will have an abundance estimate
# per quarter and grid point. We need to convert this into a "total" abundance with its standard deviation.
# This is what the present script does.



rm(list=ls())





setwd("~/FishCost/Results/") # Select appropriate directory


temp <- list.files(pattern=".RData")


DATAID <- c("commercial","survey","both")[3]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load data and set into list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pat   <- gregexpr(pattern ='_scenario',temp)
d <- temp[which(pat==10)] #1-9
dd <- temp[which(pat==11)] #10-99
ddd <- temp[which(pat==12)] #100-500


# Load datafiles with numbers 1-9
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:length(d)){
  nsimu <- substr(d[i],9,9)
  
  if(DATAID=="commercial"){
  filename <- paste("com",nsimu, sep="_")
  } else if(DATAID=="survey"){
  filename <- paste("sur",nsimu, sep="_")
  } else if(DATAID=="both"){
  filename <- paste("both",nsimu, sep="_")
  }
  
  assign(paste(filename), readRDS(d[i]))
}


# Load datafiles with numbers 10-99
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:length(dd)){
  nsimu <- substr(dd[i],9,10)
  if(DATAID=="commercial"){
  filename <- paste("com",nsimu, sep="_")
  } else if(DATAID=="survey"){
  filename <- paste("sur",nsimu, sep="_")
  } else if(DATAID=="both"){
  filename <- paste("both",nsimu, sep="_")
  }
  
  assign(paste(filename), readRDS(dd[i]))
}

# Load datafiles with numbers 100-500
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:length(ddd)){
  nsimu <- substr(ddd[i],9,11)
  if(DATAID=="commercial"){
  filename <- paste("com",nsimu, sep="_")
  } else if(DATAID=="survey"){
  filename <- paste("sur",nsimu, sep="_")
  } else if(DATAID=="both"){
  filename <- paste("both",nsimu, sep="_")
  }
  
  assign(paste(filename), readRDS(ddd[i]))
}


# Set everyhting into a list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reslist  <- list()


if(DATAID=="commercial"){
reslist <- mget(ls(pattern =  "*com_"))
} else if (DATAID=="survey"){
reslist <- mget(ls(pattern =  "*sur_"))
} else if (DATAID=="both"){
reslist <- mget(ls(pattern =  "*both_"))
}

if(DATAID=="commercial"){
rm(list=ls(pattern=c("com_","d","dd","ddd")))
} else if(DATAID=="survey"){
rm(list=ls(pattern=c("sur_","d","dd","ddd")))
} else if(DATAID=="both"){
rm(list=ls(pattern=c("both_","d","dd","ddd")))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Get mean abundance and std. error for each grid point
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# First we need to convert the abundances to natural scale.
# Not sure about the SE, though. I will apply for it as well.

# for(i in seq_along(reslist)){
# reslist[[i]]$meanAbu <- rowMeans(reslist[[i]][3:6])
# reslist[[i]]$meanSE <- rowMeans(reslist[[i]][7:10])
# }

# 2.1) Convert to natural scale (exp values)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reslistexp <- list()
for(i in seq_along(reslist)){
  reslistexp[[i]] = apply(reslist[[i]][3:10],2,exp)
}
reslistexp <- lapply(reslistexp ,as.data.frame) #convert to a df, otherwise we will have trouble in the next step



# 2.2) Go for the mean abundance and SE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(reslistexp)){
  reslistexp[[i]]$meanAbu <- rowMeans(reslistexp[[i]][1:4])
  reslistexp[[i]]$meanSE <- rowMeans(reslistexp[[i]][5:8])
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Get total mean abundance and std. error for the whole year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = data.frame(totAbund = rep(0, 500), meanSEmodel = rep(0,500), medianSEmodel = rep(0,500))

for(i in seq_along(reslistexp)){
  totAbund <- colSums(reslistexp[[i]][9],na.rm=T)
  meanSEmodel <- colMeans(reslistexp[[i]][10], na.rm=T)
  medianSEmodel <- median(reslistexp[[i]]$meanSE,na.rm=T)
  df[i,] <- c(totAbund, meanSEmodel, medianSEmodel)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Include additional information on the dataframe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wd  <- getwd() 


if(DATAID=="commercial"){
#pat2 <- substr(wd,82,105) #commercial 2016
pat2 <- substr(wd,91,120) #commercial 2015
} else if(DATAID=="survey"){
#pat2 <- substr(wd,82,105) #survey 2016
pat2 <- substr(wd,91,120) #survey 2015
} else if(DATAID=="both"){
#pat2B <- substr(wd,87,105) #both 2016
pat2B <- substr(wd,96,120) #both 2015
}


if(DATAID=="commercial"){
scenario <- substr(pat2,16,25) # commercial
} else if(DATAID=="survey"){
scenario <- substr(pat2,12,25) #survey
} else if(DATAID=="both"){
scenario <- substr(pat2B,6,25) # both
}


df$scenario <- paste(scenario)



# Get overall SD, SE, CV and VARIANCE from the abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
se <- function(x) sd(x)/sqrt(length(x)) #function to calculate standard error


df$sd <- sd(df$totAbund)
df$se <- se(df$totAbund)
df$cv <- sd(df$totAbund)/mean(df$totAbund)
df$variance <- var(df$totAbund) 


if(DATAID=="commercial"){
datid <- substr(pat2,1,3)
} else if(DATAID=="survey"){
datid <- substr(pat2,1,3)
} else if (DATAID=="both"){
datid <- substr(pat2B,1,4)
}


df$dataID <- datid


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Save the final results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_02/Results/Finalized/")
setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_02/Results/Cod/2015/Finalized/")


#saveRDS(df, paste(datid, scenario, sep="_",".rds"))
saveRDS(df, paste(scenario,".rds",sep=""))



#################################################################################



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Extra: Reference models (do not have 500 runs)
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_02/Results/Reference_models/") # Select appropriate directory
# 
# tmp  <- readRDS("BOTH_1_0_1.RData")
# 
# 
# 
# # Get mean abundance and std. error for each grid point
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # First convert to natural scale (exp values)
# tmp2 <- as.data.frame(apply(tmp[3:10],2,exp))
# 
# tmp2$meanAbu <- rowMeans(tmp2[,1:4])
# tmp2$meanSE <- rowMeans(tmp2[,5:8])
# 
# 
# # Get total mean abundance and std. error for the whole year
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df = data.frame(totAbund = rep(0, 1), meanSEmodel = rep(0,1), medianSEmodel = rep(0,1))
# 
# df$totAbund <- sum(tmp2$meanAbu,na.rm=T)
# df$meanSEmodel <- mean(tmp2$meanSE, na.rm=T)
# df$medianSEmodel <- median(tmp2$meanSE,na.rm=T)
# 
# 
# 
# # Include additional information on the dataframe
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df$scenario <- paste("1_0_1")
# 
# 
# # Get overall SD and SE from the abundances
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# se <- function(x) sd(x)/sqrt(length(x)) #function to calculate standard error
# 
# 
# df$sd <- sd(df$totAbund)
# df$se <- se(df$totAbund)
# df$cv <- sd(df$totAbund)/mean(df$totAbund)
# df$variance <- var(df$totAbund) 
# 
# df$dataID <- "both"
# 
# 
# 
# # Save the final results
# #~~~~~~~~~~~~~~~~~~~~~~~~~
# setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/PhD/Manuscript_02/Results/Finalized/")
# 
# saveRDS(df, "1_0_1.rds")
