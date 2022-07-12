##########################################################################################
#                                                                                        #
##              Benchmarking the cost-effectiveness of commercial fisheries             ##
##                       and scientific survey sampling programmes                      ##
##                    Script No. 2: Processing of simulated abundances                   #
#                                                                                        #
##########################################################################################

# Last update: July, 2022


# Code written and mantained by Marie-Christine Rufener
# Contact < macrufener@gmail.com > for any query or to report code issues.


# This script is to be performed after the abundance simulations (Script No. 1 - SimuAbu)
# The LGNB model provides abundance estimates for each grid point and input time-period.
# As the default model configuration was set on a quarterly level, this means that we will 
# have an abundance estimate per quarter and grid point.
# We need to convert this into a "total" abundance with its standard deviation.
# This is what the present script does.


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


rm(list=ls())



#~~~~~~~~~~~~~~~~~~
# 1) Load results
#~~~~~~~~~~~~~~~~~~

#Code lines below are not the pretties, but they work :)



# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [1] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [1] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("~","FishCost", "Results", "SimuAbundance", SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) List directories within selected working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmpdir <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
tmpdir2 <- tmpdir[2:length(tmpdir)] #Keep out result from fixed scenario (those that were not simulated)
tmpdir2 <- sub('..', '', tmpdir2) #Remove first two characters of directories


# 1.4) List result files
#~~~~~~~~~~~~~~~~~~~~~~~~
tmpfile <- list()
pat <- list()


for(i in seq_along(tmpdir2)){ 
  
  setwd(paste(wd,tmpdir2[i], sep="/"))
  
  tmpfile[[i]] <- list.files(pattern=".RData")
  pat[[i]]   <- gregexpr(pattern ='_scenario',tmpfile[[i]])
  
}



# 1.5) Load scenario-specific results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#We need to indicate which group of simulation indices we want to load.
#There are 3 groups of simulation: 1<i<9, 10<i<99, and i<100<500 (where i=simulation number)
#For the toy example, only 5 simulations were run (i=1,...,5). 
#However, in original paper i=1,...,500. The code lines below thus expands to cases where 1<i<500 


d <- list()
dd <- list()
ddd <- list()

for(i in seq_along(tmpfile)){

d[[i]] <- tmpfile[[i]][which(pat[[i]]==10)] #Get simulation indices from 1-9
dd[[i]] <- tmpfile[[i]][which(pat[[i]]==11)] #...from 10-99
ddd[[i]] <- tmpfile[[i]][which(pat[[i]]==12)] #...from 100-500
  
}



## Load results from simulation indices 1-9
for(j in seq_along(d)){
  for(s in seq_along(d[[j]])){
    setwd(paste(wd,tmpdir2[j], sep="/"))
    nsimu <- substr(d[[j]][s],9,9)
    
    filename <- paste(tmpdir2[[j]],"_","S",nsimu,sep="")
    assign(paste(filename), readRDS(d[[j]][s]))
    
  }
}



## Load results from simulation indices 10-99
for(j in seq_along(dd)){
  for(s in seq_along(dd[[j]])){
    setwd(paste(wd,tmpdir2[j], sep="/"))
    nsimu <- substr(dd[[j]][s],9,10)
    
    filename <- paste(tmpdir2[[j]],"_","S",nsimu,sep="")
    assign(paste(filename), readRDS(dd[[j]][s]))
    
  }
}


## Load datafiles with numbers 100-500
for(j in seq_along(ddd)){
  for(s in seq_along(ddd[[j]])){
    setwd(paste(wd,tmpdir2[j], sep="/"))
    nsimu <- substr(ddd[[j]][s],9,11)
    
    filename <- paste(tmpdir2[[j]],"_","S",nsimu,sep="")
    assign(paste(filename), readRDS(ddd[[j]][s]))
    
  }
}


# 1.6) Set scenario-specific results into a list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reslist  <- list()
p <- list()

for(i in seq_along(tmpdir2)){
  
  if(DATA == "both"){
  p[[i]] <- substr(tmpdir2[[i]],6,20)

  } else if(DATA=="commercial"){
  p[[i]] <- substr(tmpdir2[[i]],12,20)

  } else if(DATA == "survey"){
  p[[i]] <- substr(tmpdir2[[i]],8,20)
  }

 reslist[[i]] <- mget(ls(pattern =  p[[i]]))
 #reslist[[i]] <- mget(ls(pattern =  grep("^",p[[i]],"$",value=T))) #FIXME 


}
### FIXME   
## Fix bug above - for now the bug is fixed manually later below in the script".
# Problem happens when sampling strategies ends with "0" (e.g., 025_025_0, 025_025_025, ..._05,..._075). 
# For those cases, the mget function above inlcudes all sampling strategies in the list where the last scenario
# includes any level of haul selections (except when 100%) - for the example above: 025_025_0, 025_025_025, 025_025_05, 025_025_075, 025_025_1
# see: lapply(reslist, print(length)) #in some cases length>nsimu. length==nsimu


names(reslist) <- tmpdir2





# 1.7) Remove unused objects from environment
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#To make space
rm(list=setdiff(ls(), c("reslist","tmpdir2","DATA","SPECIES","YEAR")))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Get mean abundance for each grid point
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2.1) Convert to natural scale (exp values)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First we need to convert the abundances to natural scale.
# LGNB output is given on a log-scale!

reslistexp <- vector('list', length(reslist))
for(i in seq_along(reslist)){
  for(j in seq_along(reslist[[i]])){
  reslistexp[[i]][[j]] = exp(reslist[[i]][[j]][3:6]) #Select only columns with abundance estimates
  }
}



# 2.2) Get mean abundance 
#~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(reslistexp)){
  for(j in seq_along(reslistexp[[i]])){
  reslistexp[[i]][[j]]$meanAbu <- rowMeans(reslistexp[[i]][[j]][1:4])
  }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Get total abundance for the whole year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsimu <- 5 #Define here the total amount of simulation; toy example nsim=5, but in paper nsimu=500


df <- vector('list', length(reslistexp))
for(i in seq_along(reslistexp)){
  for(j in seq_along(reslistexp[[i]])){
    df[[i]][[j]]<- colSums(reslistexp[[i]][[j]][5],na.rm=T)#take column with "mean Abundance"
  }
}


df2 <- list()
for(i in seq_along(reslistexp)){
  df2[[i]] <- as.data.frame(df[[i]]) 
  colnames(df2[[i]]) <- "TotAbu"
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Some additional manipulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# 4.1) Include additional information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(df2)){
  
  if(DATA == "both"){
    df2[[i]]$Data <- substr(tmpdir2[[i]],1,4) 
    df2[[i]]$SamplingStrategy <- substr(tmpdir2[[i]],6,25) 
    df2[[i]]$nsimu <- seq(1:nsimu)
    
  } else if(DATA == "commercial"){
    df2[[i]]$Data <- substr(tmpdir2[[i]],1,10) 
    df2[[i]]$SamplingStrategy <- substr(tmpdir2[[i]],12,25) 
    df2[[i]]$nsimu <- seq(1:nsimu)

  } else if(DATA == "survey"){
    df2[[i]]$Data <- substr(tmpdir2[[i]],1,6) 
    df2[[i]]$SamplingStrategy <- substr(tmpdir2[[i]],8,25)  
    df2[[i]]$nsimu <- seq(1:nsimu)
  }
}




# 4.2) Fix bug reported in section 1.6
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df3 <- list()
for(i in seq_along(df2)){
  df3[[i]] <- df2[[i]][1:nsimu,] #Keep only the first rows with length=nsimu in each data frame
}

lapply(df3, dim) #Check



# 4.3) Bind lists into a dingle df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- do.call("rbind", df3)
stopifnot(isTRUE(dim(dat)[1]==(length(tmpdir2)*nsimu)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Include not simulated results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Results related to those scenarios where no simulation were needed.
#Includes the following scenarios:
#     - (Com=0, SurQ1=1, SurQ4=1)
#     - (Com=0, SurQ1=0, SurQ4=1)
#     - (Com=0, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=0)
#     - (Com=1, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=1)
#     - (Com=1, SurQ1=1, SurQ4=1)



# 5.1) Get file names
#~~~~~~~~~~~~~~~~~~~~~~
setwd("../")
extrafiles <- list.files(pattern="*.RData") #Replace with rds if necessary


# 5.2) Load results
#~~~~~~~~~~~~~~~~~~~~
resextra <- list()
for(i in seq_along(extrafiles)){
  filename <- sub('\\.RData$', '', extrafiles[i]) #Replace with rds if necessary
  resextra[[i]]<-assign(paste(filename), readRDS(extrafiles[i]))
}



# 5.3) Convert to natural scale (exp values)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
resextraexp <- vector('list', length(resextra))
for(i in seq_along(resextra)){
    resextraexp[[i]] = exp(resextra[[i]][3:6]) #Select only columns with abundance estimates
}


# 5.4) Get mean abundance 
#~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(resextraexp)){
  resextraexp[[i]]$meanAbu <- rowMeans(resextraexp[[i]][1:4])
}



# 5.5) Get total abundance for the whole year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dfextra <- list()
for(i in seq_along(resextraexp)){
    dfextra[[i]]<- as.data.frame(colSums(resextraexp[[i]][5],na.rm=T))#take column with "mean Abundance"
    colnames(dfextra[[i]]) <- "TotAbu"
    rownames(dfextra[[i]])<- NULL
}



# 5.6) Include additional information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(dfextra)){
  filename[[i]] <- sub('\\.RData$', '', extrafiles[i]) 
  
  if(DATA == "both"){
    dfextra[[i]]$Data <- substr(filename[i],1,4) 
    dfextra[[i]]$SamplingStrategy <- substr(filename[[i]],6,25) 
    dfextra[[i]]$nsimu <- NA
  
  } else if(DATA == "commercial"){
    dfextra[[i]]$Data <- substr(filename[i],1,10) 
    dfextra[[i]]$SamplingStrategy <- substr(filename[[i]],12,25) 
    dfextra[[i]]$nsimu <- NA
    
  } else if(DATA == "survey"){
    dfextra[[i]]$Data <- substr(filename[i],1,6) 
    dfextra[[i]]$SamplingStrategy <- substr(filename[[i]],8,25) 
    dfextra[[i]]$nsimu <- NA
  }
 
}




# 5.7) Bind lists into a dingle df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datext <- do.call("rbind", dfextra)




#~~~~~~~~~~~~~~~~~~~~~~
# 6) Bind the two dfs
#~~~~~~~~~~~~~~~~~~~~~~
dfall <- rbind(dat, datext) 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7) Finally, save the results...
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(.Platform$OS.type == "windows") setwd("~/FishCost/Results/SimuAbundance")

wd <- getwd()
dir.create(paste0(wd,"/",SPECIES,"/",YEAR,"/","Processed"),recursive=T, showWarnings = FALSE)

setwd(paste0(wd,"/",SPECIES,"/",YEAR,"/","Processed"))
OUTFILE  <- paste0("Processed_SimuResults_", DATA, ".rds")
saveRDS(dfall,file=OUTFILE)


