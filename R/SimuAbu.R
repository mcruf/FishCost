##########################################################################################
#                                                                                        #
##              Benchmarking the cost-effectiveness of commercial fisheries             ##
##                       and scientific survey sampling programmes                      ##
##                          Script No. 1: Simulating abundances                          #
#                                                                                        #
##########################################################################################

# Last update: June, 2021


# Code written and mantained by Marie-Christine Rufener
# Contact < macrufener@gmail.com > for any query or to report code issues.



# This script selects randomly samples (i.e. hauls) given the input scenario* to
# generate abundance indices. Abundances are estimated according to the LGNB spatio-temporal 
# species distibution model (for details see Rufener et al., 2021; doi: ..... ).
# For exemplification, data from the Danish sampling programmes are used along this GitHub repository (as in the paper).
# This includes data from (i) IBTS survey and (ii) on-board observers program**. 



## *Note:
##  A scenario could be, for example, defined as selecting 75% of the commercial data,
##  25% of the survey Q1 data and 50% of the survey Q4 data. The simulation is then made upon 
##  selecting randomly the data percentages within each scenario.

## **Note: 
##   DTU Aqua has a data agreement with the Danish Ministry for Food, Agriculture and Fisheries  
##   where DTU receive commercial fisheries data as part of an agreement on science based advice, 
##   to be used for obligations under the EU Data Collection Framework (EU 2017/1004), advice  and research.
##   DTU Aqua does not have permission to forward these data un-aggregated to third part due to data
##   sensitivity under the GDPR regulation. As a toy dataset for the current repository, a toy dataset is provided
##   which is based on a smaller fraction of the original data, where critical information were omitted and 
##   haul locations fully simulated.


#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


#~~~~~~~~~~~~~~~~~~~~
# 1) Default inputs
#~~~~~~~~~~~~~~~~~~~~

nsimu <- 10 #No. of simulations; in original paper, nsimu=500



# For scripting
#~~~~~~~~~~~~~~~~
# Especially useful when running script on hpc through Makefile
YEAR <- c("2015","2016")[2] #Choose year from data for which abundance index should be estimated; default is 2016
SPECIES <- c("cod","plaice","herring") #Choose species of interest
DATA  <- c("commercial", "survey", "both") [3] #

input <- parse(text=Sys.getenv("SCRIPT_INPUT")) #See Makefile in folder
print(input)
eval(input)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Loading basic functions and libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(.Platform$OS.type == "windows") setwd("~/FishCost/src/")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utilities.R")

#devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct")
mLoad(raster,rgeos,maptools,maps,data.table,dplyr,TMB,sp,DATRAS,gridConstruct,rgdal,geosphere,devtools,plyr,fields,forcats,gtools)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Setting-up the scenarios
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We have in total 125 scenarios (5*5*5); However, 6 scenarios will be removed, namely:
# 1) scenario where no data are selected at all from the three cases (Com=0, SurQ1=0, SurQ4=0)
# 2) scenario where the full dataset are selected for all three cases (Com=1, SurQ1=1, SurQ4=1); This needs to be run only once, and not 500 times
# 3) scenario where the full dataset is selected for either commercial or survey when the oter is null
# (Com=0, SurQ1=1, SurQ4=1), (Com=0, SurQ1=0, SurQ4=1), (Com=0, SurQ1=1, SurQ4=0), (Com=1, SurQ1=0, SurQ4=0), 



# Generate all possible combinations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos <- c(0, 0.25, 0.50, 0.75, 1)
scenarios <- as.data.frame(permutations(n=5,r=3,v=pos,repeats.allowed=T)) # Create all possible combinations (125 scenarios)
colnames(scenarios) <- c("COM","SURQ1","SURQ4")

# Remove the some unnecessary scenarios 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx0 <- which(rowSums(scenarios)==0) #(Com=0, SurQ1=0, SurQ4=0) 
idx1 <- which(rowSums(scenarios)==3) #(Com=1, SurQ1=1, SurQ4=1)

scenarios <- scenarios[-c(idx0,idx1),]
rownames(scenarios) <- NULL


# Include a data-specific identifier column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idxSur <- which(scenarios[,1]==0) #Identifies scenarios where no commercial data is used (thus, only survey)
idxCom <- which(scenarios[,2]==0 & scenarios[,3]==0) #Identifies scenarios where no survey data is used (thus, only commercial)

scenarios$Data <- "both"
scenarios[idxSur,"Data"] <- "survey"
scenarios[idxCom,"Data"] <- "commercial"

scenarios$structure <- paste(scenarios[,1],scenarios[,2],scenarios[,3],scenarios[,4],sep="_")


# Remove additional unnecessary scenarios
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx2 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==1) #(Com=0, SurQ1=1, SurQ4=1) 
idx3 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==0) #(Com=1, SurQ1=0, SurQ4=0) 
idx4 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==0) #(Com=0, SurQ1=1, SurQ4=0) 
idx5 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==1) #(Com=0, SurQ1=0, SurQ4=1) 
idx6 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==0) #(Com=1, SurQ1=1, SurQ4=0) 
idx7 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==1) #(Com=1, SurQ1=0, SurQ4=1) 

scenarios <- scenarios[-c(idx2,idx3,idx4,idx5,idx6,idx7),]
rownames(scenarios) <- NULL


scenarios <- subset(scenarios, Data==DATA)
rownames(scenarios) <- NULL



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Load fishery-depedent and independent datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load fishery-depedent dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("commercial_NEW.RData"); comFULL <- tester2; rm(tester2)  #obo-dfad dataset
#comFULL <- readRDS("commercial.rds") #Full commercial dataset (obo-dfad-vms-logbook)


#Subset dataset according to pre-defined inputs
com <- subset(comFULL,Year == YEAR) # if only WBS cod should be modelled



## Create 6+ group
com[,"age_6"] <- rowSums(com[,c("age_6","age_7","age_8")])
com$age_7 <- NULL
com$age_8 <- NULL


## Create a Ntot column (sum over all ages, total abundance)
com$Ntot <- rowSums(com[,c("age_1","age_2","age_3","age_4","age_5","age_6")])

## Include an age-0 column (to be equal as the survey)
# NOTE: I HAVE TO CHECK WHETER THIS HAS AN EFFECT WHEN MODELLING age_0 
com$age_0 <- rep(0,nrow(com))




# Load fishery-indepedent dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sur <- readRDS("survey_NEW.rds") 
sur$Data  <- as.factor(rep("survey",nrow(sur)))
sur$timeyear  <- as.factor(paste(sur$Year,sur$Data,sep=":"))
sur$timeyear2  <- as.factor(paste(sur$Year,sur$Quarter,sur$Data,sep=":"))






## Create a Ntot column (sum over all ages, total abundance)
survey$Ntot <- rowSums(survey[,c("age_0","age_1","age_2","age_3","age_4","age_5","age_6")])




stopifnot(levels(factor(com$Year))==levels(factor(survey$Year))) # Security check





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Build grid for the study area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Creating a dataframe with mean values of long and lat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comFULL$lon_mean <- rowMeans(comFULL[,c("lonStart", "lonEnd")])
comFULL$lat_mean <- rowMeans(comFULL[,c("latStart", "latEnd")])
df <- data.frame(lon=comFULL$lon_mean, lat=comFULL$lat_mean)


# Building the grid
#~~~~~~~~~~~~~~~~~~~~
grid <- gridConstruct3(df,km=5,scale=1.2)
gr <- gridFilter(grid,df,icesSquare = T,connected=T)
# plot(gr)
# plot(DK_map,add=T,fill=T,col="grey70")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Start with the simulations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(s in 1:nrow(scenarios)){
  for(i in 1:nsimu){
    
    if(.Platform$OS.type == "windows") setwd("C:/Users/mruf/Documents/PHD_projects/Proj_2/CostEffect/Cod")
    
    setwd("/zhome/25/f/124809/PHD_projects/Proj_2/CostEffect/Cod")
    
    
    
    # 6.1) Define the scenario to be simulated
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCENARIO_COM <- strsplit(scenarios$structure, "_")[[s]][1]
    SCENARIO_SURQ1 <- strsplit(scenarios$structure, "_")[[s]][2]
    SCENARIO_SURQ4 <- strsplit(scenarios$structure, "_")[[s]][3]
    #DATA <- strsplit(scenarios$structure, "_")[[s]][4]
    
    ## Validate script input
    stopifnot(DATA %in% c("commercial", "survey", "both"))
    
    
    # 6.2) Define the support area
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUPPORT_AREA <- c("One","Several")[1] #Choose whether to use one or several support areas to describe the commercial fisheries data; Survey data is described by only single support area
    ALPHA <- c("No","Single I","Single II", "Multi")[2] #Choose whether anhd how model should consider preferrential sampling parameter (alpha) or not.
    # Check for the NBCP script to see the alpha & support_area speicifcations
    
    # We go from the assumption that survey data has ALWAYS only one support area; one can choose whether to estimate the alpha-parameter or not. 
    # if(DATA=="survey"){
    #   SUPPORT_AREA <- "One" 
    #   ALPHA <- c("No","Single I")[2] # Default is Multi alphas for commercial; this automatically implies in several support areas
    # }
    
    
    
    
    # 6.3) Subset dataframes according to the simulated scenario 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    # 6.3.1) Fishery dependent data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Sample only a fraction of the commercial data (fractions are the SCENARIOS_COM stated in the beginning)
    commercial <- as.data.frame(com %>% sample_frac(as.numeric(SCENARIO_COM)))
    
    
    ## Remove unused columns
    commercial[,c("logBldNr","sampleId","cruise","catch_kg_observer","lon_mean","lat_mean")] <- NULL
    commercial[,c("Area","efid","Quarter","Day","Month","Year","HLID","metiers_g1","metiers_g2","Data","timeyear","timeyear2")] <- lapply(commercial[,c("Area","efid","Quarter","Day","Month","Year","HLID","metiers_g1","metiers_g2","Data","timeyear","timeyear2")],factor)
    # OBS.: NOte that we can have problems when a scenario includes only one level either in VE_LENcat or metiers; thats why I don not drop those levels in the above line
    
    
    # 6.3.2) Fishery independent data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ## Subset data for first and fourth quarter, sample a fraction of each data, and then bind them together again
    surveyQ1 <- subset(survey, Quarter =="1")
    surveyQ4 <- subset(survey, Quarter =="4")
    
    surveyQ1 <- as.data.frame(surveyQ1 %>% sample_frac(as.numeric(SCENARIO_SURQ1)))
    surveyQ4 <- as.data.frame(surveyQ4 %>% sample_frac(as.numeric(SCENARIO_SURQ4)))
    
    survey2  <- rbind(surveyQ1,surveyQ4)
    
    
    ## Subset the most relevant columns from the survey data / reorder columns
    colsel_sur <- c("Quarter","Year","Month","Day","Area","HaulDur","Depth",
                    "Data","haul.id","lon","lat","timeyear","timeyear2","Ntot",paste("age",0:6,sep="_")) #column selection for survey data
    
    survey2 <- subset(survey2, select=colsel_sur)
    
    survey2[,c("Month","Year","Quarter","Area","Data","timeyear","timeyear2")] <- lapply(survey2[,c("Month","Year","Quarter","Area","Data","timeyear","timeyear2")], factor)
    
    
    
    
    # 6.4) Bind survey and commercial datasets into a single dataset
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    age_com   <- transform(commercial, HaulDur=haulduration_hours, numYear = as.numeric( as.character(Year)))
    age_sur   <- transform(survey2, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,
                           HLID=haul.id, numYear = as.numeric(as.character(Year)))
    
    if (DATA == "both"){
      datatot <- mybind(age_com, age_sur)
    } else if (DATA == "survey") {
      datatot <- age_sur
    } else if (DATA == "commercial")
      datatot <- age_com
    
    
    
    ## Create equally time spaced intervals - VERY important for the AR1 process
    timeLevels <- as.vector(t(outer(min(datatot$numYear):max(datatot$numYear), 1:4, paste))) # Decide: quarterly or monthly basis
    datatot$YearQuarter <- factor(paste(datatot$Year, datatot$Quarter), levels=timeLevels)
    
    
    
    
    # 6.5) Define response variable 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if(AGE == "A0"){
      datatot$Response <- as.numeric(paste(datatot$age_0))
    } else if(AGE == "A1"){
      datatot$Response <- as.numeric(paste(datatot$age_1))
    } else if (AGE == "A2") {
      datatot$Response <- as.numeric(paste(datatot$age_2))
    } else if (AGE == "A3"){
      datatot$Response <- as.numeric(paste(datatot$age_3))
    } else if (AGE == "A4"){
      datatot$Response <- as.numeric(paste(datatot$age_4))
    } else if (AGE == "A5"){
      datatot$Response <- as.numeric(paste(datatot$age_5))  
    } else if (AGE == "A6"){
      datatot$Response <- as.numeric(paste(datatot$age_6))
    } else if(AGE == "Ntot"){
      datatot$Response <- as.numeric(paste(datatot$Ntot))
    }
    
    
    
    # 6.6) Discretize and associate hauls along grid cells 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    # 6.6.1) Create a data frame containing the trawl ID, starting long/lat and ending long/lat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dat <- data.frame(sampleID=datatot$HLID, start_long=datatot$lonStart,
                      start_lat=datatot$latStart, end_long=datatot$lonEnd, end_lat=datatot$latEnd)
    #segments(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat, col="red",lwd=1.2)
    
    
    ## Define matrix for start/end long & lat
    p1 <- matrix(c(dat$start_long,dat$start_lat), ncol=2)
    p2 <- matrix(c(dat$end_long,dat$end_lat), ncol=2)
    
    
    # 6.6.2) Interpolate intermediate points at regular distance (default is 1km)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nbpts <- floor(distance(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat) / 1) # one point every 1 km ~ reasonable for a 5x5 km grid
    
    ## Get intermediate points 
    inter_pts <- gcIntermediate(p1,p2, n=nbpts, addStartEnd=FALSE) #inter_pts returns a list within several list. We need to convert this list to a single dataframe
    
    
    # 6.6.3) Create a dataframe to specify the frequency of match between each haul and grid ID. 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datatot$HLID <- factor(datatot$HLID ) #MUST be a factor
    stopifnot(nlevels(datatot$HLID)==dim(datatot)[1])
    
    tmp <- lapply(1:length(inter_pts), function(i) {
      print(i)
      x <- inter_pts[[i]]
      colnames(x) <- c("lon", "lat")
      x <- as.data.frame(x)
      haul.id <- datatot$HLID[i]
      ind <- gridLocate(gr, x)
      data.frame(haul.id=haul.id, ind=ind, rowID = i) 
    })
    tmp2 <- do.call("rbind", tmp)
    tmp2$haulid <- factor(tmp2$haul.id)
    tmp2$gf <- factor(tmp2$ind, 1:nrow(gr))
    tmp2 <- tmp2[c("haulid","gf","rowID")] 
    
    
    
    
    # 6.7) Define support areas (to be used later in association to the alpha-parameter)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 6.7.1) For single support area
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Applied when using only one alpha parameter for each dataset)
    datatot$split_area <- ifelse(as.character(datatot$Data)=="commercial", "commercial", "survey") #Same configuration as model with a single alpha parameter, where we have only one single support area describes the commercial data by aggregating all hauls of the time series.
    datatot$split_area <- as.factor(datatot$split_area)
    
    kk <- tmp2; kk$split <- datatot$split_area[tmp2$rowID]
    SupportAreaMatrix    <- table(kk$gf, kk$split)
    SupportAreaMatrix[]  <- SupportAreaMatrix>0
    SupportAreaMatrix    <- ifelse(SupportAreaMatrix==0,FALSE,TRUE)
    
    
    # 6.7.2) For multiple support areas
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Applied when using either one alpha parameter (collapsing all alphas into one), 
    # or multiple alphas for commercial data, and one alpha parameter for survey)
    
    if(SUPPORT_AREA == "Several"){
      datatot$split_area2 <- ifelse(as.character(datatot$Data)=="commercial", as.character(datatot$Quarter), "survey") #Same configuration as model with a single alpha parameter, where we have only one single support area describing the commercial data.
      datatot$split_area2 <- as.factor(datatot$split_area2)
      
      kk2 <- tmp2; kk2$split <- datatot$split_area2[tmp2$rowID]
      SupportAreaMatrix2 <- table(kk2$gf, kk2$split)
      SupportAreaMatrix2[] <- SupportAreaMatrix2>0
      SupportAreaMatrix2 <- ifelse(SupportAreaMatrix2==0,FALSE,TRUE)
    }
    
    
    # 6.7.3) Setting support areas based on chosen input
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(DATA == "commercial" & SUPPORT_AREA == "Several"){
      SupportAreaMatrix2[] <- SupportAreaMatrix[,1]
      SupportAreaMatrix <- SupportAreaMatrix2 
    } else if(DATA == "both" & SUPPORT_AREA == "Several"){
      SupportAreaMatrix2[,1:(ncol(SupportAreaMatrix2)-1)] <- SupportAreaMatrix[,1]
      SupportAreaMatrix <- SupportAreaMatrix2 
    } else if(DATA == "commercial" & SUPPORT_AREA == "One"){
      SupportAreaMatrix <- SupportAreaMatrix
    } else if (DATA == "both" & SUPPORT_AREA =="One"){
      SupportAreaMatrix <- SupportAreaMatrix
    }
    
    
    ## Map haulid to data.frame rows (VERY IMPORTANT: haulid must match with the dataframe's row number (in increasing order))
    rowid <- match(as.character(tmp2$haulid),as.character(datatot$HLID))
    stopifnot(all(is.finite(rowid)))
    rowid <- factor(rowid)
    
    
    
    # 6.8) Specifying spatio-temporal model 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 6.8.1) Sparse matrices for GMRF: Q = Q0+delta*I
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Q0 <- -attr(gr,"pattern")
    diag(Q0) <- 0
    diag(Q0) <- -rowSums(Q0)
    I <- .symDiagonal(nrow(Q0))
    
    
    # 6.8.2) Complile TMB
    #~~~~~~~~~~~~~~~~~~~~~
    compile("model.cpp")
    dyn.load(dynlib("model"))
    
    
    
    # 6.8.3) TMB Data
    #~~~~~~~~~~~~~~~~
    # TMB data are set in such way that it recognizes automatically wheter one is using FD or FID data. 
    data <- list(
      time = datatot$YearQuarter[rowid],
      gf = tmp2$gf  ,              
      rowid = rowid  ,
      response = datatot$Response,
      Q0 = Q0,
      I = I,
      Xpredict = matrix(0,0,0),
      Apredict = factor(numeric(0)),
      SupportAreaMatrix = SupportAreaMatrix,
      SupportAreaGroup = if(SUPPORT_AREA == "Several"){
        as.factor(datatot$split_area2)
      } else if(SUPPORT_AREA == "One"){
        as.factor(datatot$split_area)
      },
      Data=datatot$Data,
      h = mean(summary(as.polygons(gr))$side.length) #To be used later to plot the spatial decorrelation as a function of distance
    )
    
    
    # Guess time effects (needed by index calculation)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    guess_time_effects_in_beta <- function(data, time_levels) {
      cn <- colnames(data$X)
      cn <- gsub(" ",":",time_levels)
      tl <- gsub(" ",":",time_levels)
      li <- lapply(tl, function(nm) grep(nm, cn) )
      if (any( sapply(li, length) > 1 )) {
        print(lapply(li, function(i)cn[i]))
        warning("Time effect is not unique. Index calc will be omitted")
        return (integer(0))
      }
      found <- sapply(li, function(x)x[1])
      found[is.na(found)] <- 0
      ## NOTE: Output length = length(time_levels)
      ## NOTE: NAs coded as "-1" ===> SKIP index calc or crash !!!
      as.integer(found) - 1L
    }
    
    
    
    # 6.8.4) Fitting the model
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ## Specify internal aspects of the model fitting
    fit_model <- function(data, model_struct=NULL, with_static_field=FALSE) {
      time_levels <- levels(data$time) # There must be at least one time levels; They MUST match between the two data types!
      grid_nlevels <- nlevels(data$gf)
      data$doPredict <- 0
      data$offset <- model_struct$offset 
      
      if(FALSE) { ## FIXME: prediction disabled
        ## Stuff for prediction: Dummy dataset that matches the space time grid:
        ## Xpredict: design matrix for prediction
        ## Apredict: Area factor for prediction
        ##DFpredict <- expand.grid(gf=levels(data$gf), time=levels(data$time))
        DFpredict <- expand.grid(gf=levels(data$gf), time=levels(data$time))
        ## FIXME: We should include depth and covariates here !
        ##        But that requires depth on the entire grid...
        Xpredict <- model.matrix(~time, data=DFpredict) ## <-- gear removed
        stopifnot( all( colnames(Xpredict) %in% colnames(data$X) ) ) ## Validity check
        tmp <- matrix(0, nrow(Xpredict), ncol(data$X))
        tmp[,match(colnames(Xpredict), colnames(data$X))] <- Xpredict
        data$Xpredict <- tmp
        data$Apredict <- factor(icesSquare(gr))
      }
      ## Perhaps we want measure all indices relative to a fixed reference square:
      ##   plot(cod, plot.response = FALSE)
      ## "42G1" seems appropriate
      ## data$refindex <- which(levels(data$Apredict) == "42G1")
      data$refindex <- 0 ## <-- Disable
      
      parameters <- list(
        eta_density = matrix(0,nrow(Q0),length(time_levels)),
        eta_nugget = numeric(0),
        logdelta = -4,       # Check values
        logscale = 0,         # Check values
        logsd_nugget = 0,    # Check values
        time_corr = 2,       # Check values
        beta = rep(0, ncol(data$X)),
        logphi = rep(0, nlevels(data$Data) ),
        alpha = rep(0, nlevels(data$SupportAreaGroup))
      )
      parameters$eta_static <- rep(0, grid_nlevels * with_static_field )
      parameters$logdelta_static <- rep(0, 1 * with_static_field )
      parameters$logscale_static <- rep(0, 1 * with_static_field )
      
      ## Plugin model specification
      if(!is.null(model_struct)) {
        data$X                  <- model_struct$Xf
        data$beta_r_fac         <- model_struct$beta_r_fac
        parameters$beta         <- model_struct$beta
        parameters$beta_r       <- model_struct$beta_r
        parameters$beta_r_logsd <- model_struct$beta_r_logsd
        data$which_beta_time    <- guess_time_effects_in_beta(data, time_levels)
      }
      
      ## Prior std dev on fixed effects (for robustness only)
      data$huge_sd <- 100
      
      map <- list()
      if(TRUE) map$logsd_nugget <- factor(NA)
      
      if(ALPHA == "No" & DATA == "both"){
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(ALPHA == "No" & DATA == "commercial") {
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(ALPHA == "No" & DATA == "survey"){
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(ALPHA == "Single I"){
        # map$alpha = factor(map$alpha)
      } else if(ALPHA == "Single II"){
        map$alpha <- factor(levels(data$SupportAreaGroup) != "survey") 
        map$alpha = factor(map$alpha)
      } else if(ALPHA == "Multi"){
        # map$alpha = factor(map$alpha)
      } 
      
      
      if(length(parameters$beta)>0 || length(parameters$beta)>0) {
        profile <- c("beta")
      } else {
        profile <- NULL
      }
      
      obj <- MakeADFun(data, parameters, random=c("eta_density","eta_nugget","eta_static","beta_r"),
                       profile = profile,
                       map=map, DLL="model")
      
      obj$env$tracepar <-TRUE
      
      .GlobalEnv$obj <- obj ## Copy for debugging
      runSymbolicAnalysis(obj)
      fit <- nlminb(obj$par,obj$fn,obj$gr)
      if(FALSE) { ## FIXME: disabled
        rep <- obj$report(obj$env$last.par.best)
        rownames(rep$logindex) <- levels(data$Apredict)
        colnames(rep$logindex) <- levels(data$time)
      }
      
      sdr <- sdreport(obj)
      s <- summary(sdr)
      est <- s[rownames(s) != "eta_density" & rownames(s) != "eta_nugget",]
      s <- summary(sdr,p.value = TRUE)
      s1 <- s[rownames(s) == "beta",]
      rownames(s1) <- head(colnames(data$X), nrow(s1)) # extracting names of the fixed effects
      s.fixed <- s1
      s1 <- s[rownames(s) == "beta_r",]
      rownames(s1) <- tail(colnames(data$X), nrow(s1)) # extracting names of random effects
      s.random <- s1
      
      s1
      s2 <- s[rownames(s) == "eta_density",]
      return(environment())
    }
    
    
    
    
    ## Model configurations
    buildModelMatrices <- function(fixed, random=NULL, ..., offset=NULL, data) {
      mm <- function(formula, data) {
        ##myna<-function(object,...){object[is.na(object)]<-0; object}
        mf <- model.frame(formula, data, na.action=na.pass)
        ans <- model.matrix(formula, mf)
        ans[is.na(ans)] <- 0
        ans
      }
      Xf <- mm(fixed, data=data)
      if(!is.null(random))
        Xr <- lapply(list(random, ...), mm, data=data)
      else
        Xr <- list(matrix(NA, nrow(Xf), 0))
      nr <- sapply(Xr, ncol)
      nf <- ncol(Xf)
      beta <- rep(0, nf)
      beta_r <- rep(0, sum(nr))
      beta_r_fac <- factor(rep(1:length(nr), nr))
      beta_r_logsd <- rep(0, nlevels(beta_r_fac))
      offset <- eval(offset, data)
      if(!is.null(offset)) stopifnot(is.numeric(offset)) else offset <- numeric(0)
      list(Xf=cbind(Xf, do.call("cbind", Xr)), nr=nr, nf=nf, beta=beta, beta_r=beta_r, beta_r_fac=beta_r_fac, beta_r_logsd=beta_r_logsd, offset=offset)
    }
    
    
    
    ## Fit model
    if (MODEL_FORMULA == "m2") {
      if(DATA == "commercial"){
        m2 <- buildModelMatrices(~  -1 + YearQuarter + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
      } else if(DATA == "survey"){
        m2 <- buildModelMatrices(~ -1 + YearQuarter + offset(log(HaulDur)), data=datatot)
      } else {
        m2 <- buildModelMatrices(~ -1 + Data + YearQuarter + offset(log(HaulDur)), ~metiers-1 + VE_LENcat-1, data=datatot)
      }
      env2 <- tryCatch(fit_model(data, m2, with_static_field = F), error=function(e) {})
    }
    
    
    
    if(!is.null(env2)){
      #~~~~~~~~~~~~~~~~~~~
      # 6.9) Save results
      #~~~~~~~~~~~~~~~~~~~
      #rm(list=setdiff(ls(), ls(pattern="datatot|gr"))) #CHECK IF IT WORKS AFTER MODELS HAVE BEEN RUN.
      
      pl <- as.list(env2$sdr,"Estimate"); pl <- as.data.frame(pl$eta_density) #Estimates
      colnames(pl) <- paste(levels(env2$data$time),"est",sep=".")
      
      pl.sd <- as.list(env2$sdr,"Std. Error"); pl.sd <- as.data.frame(pl.sd$eta_density) #Std. Error
      colnames(pl.sd) <- paste(levels(env2$data$time),"se",sep=".")
      
      
      abuindex <- as.numeric(summary(env2$sdr,"report")[,1])
      logindex <- c(abuindex,rep("NA",dim(gr)[1]-length(abuindex)))
      
      expindex <- exp(dfsimu$logIndex)
      totabu <- sum(expindex,na.rm=T)/length(abuindex)
      
      
      dfsimu <- cbind(gr,pl,pl.sd)
      dfsimu$logIndex <- as.numeric(as.character(logindex))
      dfsimu$totabu <- as.numeric(paste(totabu))
      
      
      
      SCENARIO <- gsub(",", ".", gsub("\\.", "", scenarios[s,5])) # s is the index for the scenario stated in the beginning of the script
      
      
      
      #rm(list=setdiff(ls(), "dfsimu"))
      
      if(.Platform$OS.type == "windows") setwd("C:/Users/mruf/Desktop/TMP")
      
      wd <- getwd()
      setwd(paste0(wd,"/","results_SIMU",sep=""))
      
      OUTFILE  <- paste0("results_", paste(i,"scenario",SCENARIO,sep="_"), ".RData")
      saveRDS(dfsimu,file=OUTFILE)
      #rm(list=ls())
    }
  }
  
}

beep("mario")
