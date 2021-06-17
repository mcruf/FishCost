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

nsimu <- 5 #No. of simulations;  nsimu=500 in original paper


# For scripting
#~~~~~~~~~~~~~~~~
# Especially useful when running script on hpc through Makefile
YEAR <- c("2015","2016")[1] #Choose year from data for which abundance index should be estimated; default is 2015
SPECIES <- c("cod","plaice","herring") [1] #Choose species of interest; default is cod
DATA  <- c("both", "commercial", "survey") [2] #Choose input data for which scenarios will be run; default is both (commercial + survey)
PS <- c("No","One","Two")[1] #Define how the sampling nature should be accounted for; default is that no preferential sampling is accounted for in the commercial
TIME <- c("YearQuarter","YearMonth")[1] #Define the temporal resolution for the model; default is set on a quarterly basis as in paper



input <- parse(text=Sys.getenv("SCRIPT_INPUT")) #See Makefile in folder
print(input)
eval(input)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Loading basic functions and libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(.Platform$OS.type == "windows") setwd("~/FishCost/src/")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utilities.R")

#devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct") #Install gridConstruct package to build grid later below
mLoad(raster,rgeos,maptools,maps,data.table,dplyr,TMB,sp,DATRAS,gridConstruct,rgdal,geosphere,devtools,plyr,fields,forcats,gtools)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Setting-up the scenarios
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We have in total 125 scenarios (5^3) based on the level of the amount of data (in %) that will be selected in
# the dataset (Commercial, Survey-Q1, Survey-Q4). From these scenarios, 8 will be removed, namely:
# 1) Case where no data are selected in any of the three datasets (when Com=0%, SurQ1=0%, SurQ4=0%)
# 2) Cases where the full amount of data are selected in all three datasets (Com=100%, SurQ1=100%, SurQ4=100%)
#     - This is the baseline adopted in the paper and needs to be run only once, and not nsim times!
# 3) Cases where the full amount of data is selected in one of the datasets and the other is not selected at all.
#     - (Com=0, SurQ1=1, SurQ4=1)
#     - (Com=0, SurQ1=0, SurQ4=1)
#     - (Com=0, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=0)
#     - (Com=1, SurQ1=1, SurQ4=0)
#     - (Com=1, SurQ1=0, SurQ4=1)
#     *Note: Any of the above options (3) could be chosen as baseline instead of the case where Com=1, SurQ1=1, SurQ4=1



# 3.1) Generate all possible data selection combination (i.e., scenrio)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos <- c(0, 0.25, 0.50, 0.75, 1) #From 0 to 100%, with 25% increase of data selection
scenarios <- as.data.frame(permutations(n=5,r=3,v=pos,repeats.allowed=T)) # Create all possible combinations (125 scenarios)
colnames(scenarios) <- c("COM","SURQ1","SURQ4")


# 3.2) Include a data-specific identifier column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idxSur <- which(scenarios[,1]==0) #Identifies scenarios where no commercial data is used (thus, only survey)
idxCom <- which(scenarios[,2]==0 & scenarios[,3]==0) #Identifies scenarios where no survey data is used (thus, only commercial)

scenarios$Data <- "both"
scenarios[idxSur,"Data"] <- "survey"
scenarios[idxCom,"Data"] <- "commercial"

scenarios$structure <- paste(scenarios[,1],scenarios[,2],scenarios[,3],scenarios[,4],sep="_")


# 3.3) Remove unused scenarios
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idx1 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==0) #(Com=0, SurQ1=0, SurQ4=0) 
idx2 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==1) #(Com=0, SurQ1=1, SurQ4=1) 
idx3 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==0) #(Com=1, SurQ1=0, SurQ4=0) 
idx4 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==0) #(Com=0, SurQ1=1, SurQ4=0) 
idx5 <- which(scenarios[,"COM"]==0 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==1) #(Com=0, SurQ1=0, SurQ4=1) 
idx6 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==0) #(Com=1, SurQ1=1, SurQ4=0) 
idx7 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==0 & scenarios[,"SURQ4"]==1) #(Com=1, SurQ1=0, SurQ4=1) 
idx8 <- which(scenarios[,"COM"]==1 & scenarios[,"SURQ1"]==1 & scenarios[,"SURQ4"]==1) #(Com=1, SurQ1=0, SurQ4=1) 

scenarios <- scenarios[-c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8),]
rownames(scenarios) <- NULL


# 3.4) Subset scenarios according to input data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#To be used later in the simulation loop
scenarios <- subset(scenarios, Data==DATA)
rownames(scenarios) <- NULL



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Load fishery-depedent and independent datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("~/FishCost/Data/") #Set path to Data folder


# 4.1) Load fishery-depedent data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
com <- readRDS(paste("com","_", SPECIES,".rds",sep=""))

com <- subset(com,Year == YEAR); com$Year <- factor(com$Year) #Subset dataset according to pre-defined inputs



# 4.2) Load fishery-indepedent data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
surFULL <- readRDS(paste("sur","_", SPECIES,".rds",sep=""))

sur <- subset(surFULL,Year == YEAR); sur$Year <- factor(sur$Year) #Subset dataset according to pre-defined inputs




### Checkpoint ###
stopifnot(levels(com$Year)==levels(sur$Year)) #Check that both datasets were subsetted for the same year




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Build grid for the study area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note: Grid for both commercial and survey data have to be the same.
# It doesn't matter wheter the grid is constructed upon the commercial or survey data, 
# as long as it is consistent across data types (commercial, survey or both)

grid <- GridConstruct(surFULL[,c("lon","lat")],km=5,scale=1.2) #Modified function from original gridConstruct. See "utilities.R" to see changes
gr <- GridFilter(grid,df,icesSquare = T,connected=T) # filter out unnecessary spatial extensions; #Modified function from original gridConstruct. See "utilities.R" to see changes
# plot(gr)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) Go for the simulations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(s in 1:nrow(scenarios)){
  for(i in 1:nsimu){
    
    
    #FIXME: CHECK WHETHER THE TWO LINES BELOW ARE REALLY NECESSARY
    #if(.Platform$OS.type == "windows") setwd("~/FishCost/")
    
    #setwd("PATH_IN_HPC") #Turn-on if running on hpc
    
    
    
    # 6.1) Define the scenario to be simulated
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCENARIO_COM   <- strsplit(scenarios$structure, "_")[[s]][1]
    SCENARIO_SURQ1 <- strsplit(scenarios$structure, "_")[[s]][2]
    SCENARIO_SURQ4 <- strsplit(scenarios$structure, "_")[[s]][3]
    #DATA <- strsplit(scenarios$structure, "_")[[s]][4]
    
    
    
    ### Checkpoint ###
    stopifnot(DATA %in% c("commercial", "survey", "both")) #Validate script input
    
    
    
    # 6.2) Subset dataframes according to input scenario 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ### 6.2.1) Fishery dependent data
    ##Sample only a fraction of the commercial data (fractions are the SCENARIOS_COM stated in the beginning)
    commercial <- as.data.frame(com %>% sample_frac(as.numeric(SCENARIO_COM)))
  
    
    ### 6.2.2) Fishery independent data
    ##Subset survey data first by quarters, and then sample a fraction of each data, and then bind them together again
    surQ1 <- subset(sur, Quarter =="1")
    surQ4 <- subset(sur, Quarter =="4")
    
    surQ1 <- as.data.frame(surQ1 %>% sample_frac(as.numeric(SCENARIO_SURQ1)))
    surQ4 <- as.data.frame(surQ4 %>% sample_frac(as.numeric(SCENARIO_SURQ4)))
    
    survey  <- rbind(surQ1,surQ4)
    
    #FIXME: Optimize above lines with something like the lines below
    #test <- split(sur, factor(sur$Quarter))
    #test2 <- lapply(teste, sample_frac, as.numeric(c(SCENARIO_SURQ1,SCENARIO_SURQ4)))
    
    
    
    # 6.3) Bind survey and commercial datasets into a single dataset
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ### 6.3.1) Standardize essential columns that are shared across the two dataset
    df_com   <- transform(commercial, HaulDur=haulduration_hours, numYear = as.numeric( as.character(Year)))
    df_sur   <- transform(survey, latStart=lat, lonStart=lon, latEnd=lat, lonEnd=lon,
                           HLID=haul.id, numYear = as.numeric(as.character(Year)))
    
    
    ### 6.3.2) Bind the two datasets
    if (DATA == "both"){
      datatot <- mybind(df_com, df_sur)
    } else if (DATA == "survey") {
      datatot <- df_sur
    } else if (DATA == "commercial")
      datatot <- df_com
    
    
    
    # 6.4) Create equally time spaced intervals 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##VERY important for the AR1 process
    if(TIME=="YearMonth"){
      timeLevels <- as.vector(t(outer(min(datatot$numYear):max(datatot$numYear), 1:12, paste)))
      datatot$YearMonth <- factor(paste(datatot$Year, datatot$Month), levels=timeLevels)
    } else if(TIME=="YearQuarter"){
      timeLevels <- as.vector(t(outer(min(datatot$numYear):max(datatot$numYear), 1:4, paste)))
      datatot$YearQuarter <- factor(paste(datatot$Year, datatot$Quarter), levels=timeLevels)
    } 
    
   
    
    # 6.5) Discretize and associate hauls along grid cells 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ### 6.5.1) Setting a data frame containing the haul ID, and start and end long/lat of the haul
    dat <- data.frame(sampleID=datatot$HLID, start_long=datatot$lonStart,
                      start_lat=datatot$latStart, end_long=datatot$lonEnd, end_lat=datatot$latEnd)
    #segments(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat, col="red",lwd=1.2)
    
    
    ### 6.5.2) Define a matrix for the haul´s start and end position
    p1 <- matrix(c(dat$start_long,dat$start_lat), ncol=2)
    p2 <- matrix(c(dat$end_long,dat$end_lat), ncol=2)
    
    
    ### 6.5.3) Interpolate points at regular distance (default is 1km)
    nbpts <- floor(distance(dat$start_long,dat$start_lat, dat$end_long, dat$end_lat) / 1) # one point every 1 km ~ reasonable for a 5x5 km grid
    inter_pts <- gcIntermediate(p1,p2, n=nbpts, addStartEnd=FALSE)
    
    
    
    # 6.6) Associate the discretized hauls to the grid ID 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##Note that the haul Id MUST be a factor, where each level is the frequency 
    ##of a particular haul crossing a specific grid ID
    
    
    ### Checkpoint ###
    stopifnot(is.factor(datatot$HLID))


    tmp <- lapply(1:length(inter_pts), function(i) {
      print(i)
      x <- inter_pts[[i]]
      colnames(x) <- c("lon", "lat") #Needs to be the same names as in inter_pts
      x <- as.data.frame(x)
      haul.id <- datatot$HLID[i] #Pick the specific haul id
      ind <- gridLocate(gr, x) #Locate the grid ID
      data.frame(haul.id=haul.id, ind=ind, rowID = i) 
    })
    tmp2 <- do.call("rbind", tmp)
    tmp2$haulid <- factor(tmp2$haul.id)
    tmp2$gf <- factor(tmp2$ind, 1:nrow(gr))
    tmp2 <- tmp2[c("haulid","gf","rowID")] 
    
    
    
    # 6.7) Defining the preferential sampling (PS)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##If a preferential sampling behaviour in the commercial and/or survey data
    ##is accounted for, we have to define the so-called sampling support area.
    ##See Rufener et al. (2021) and github.com/mcruf/LGNB for more details.
    
  
    ### 6.7.1) For single support area
    # Applied when using only one alpha parameter for each dataset)
    datatot$split_area <- ifelse(as.character(datatot$Data)=="commercial", "commercial", "survey") #Takes the haul positions of the aggregated time-series 
    datatot$split_area <- as.factor(datatot$split_area) #IMPORTANT - needs to be a factor!!
    tmpOne <- tmp2; tmpOne$split <- datatot$split_area[tmp2$rowID]
    SupportAreaMatrix    <- table(tmpOne$gf, tmpOne$split)
    SupportAreaMatrix[]  <- SupportAreaMatrix>0
    SupportAreaMatrix    <- ifelse(SupportAreaMatrix==0,FALSE,TRUE)
    #levels(datatot$split_area); levels(datatot$Data) #Check
    #image(gr, SupportAreaMatrix[,1]) #To see the progress..
    
    
    # 6.8) Map haulid to data.frame rows 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## VERY IMPORTANT: haulid must match with the dataframe's row number (in increasing order)
    rowid <- match(as.character(tmp2$haulid),as.character(datatot$HLID))
    stopifnot(all(is.finite(rowid)))
    rowid <- factor(rowid)
    
    
    
    # 6.9) TMB processing
    #~~~~~~~~~~~~~~~~~~~~~~
    ### 6.9.1) Sparse matrices for GMRF: Q = Q0+delta*I
    Q0 <- -attr(gr,"pattern")
    diag(Q0) <- 0
    diag(Q0) <- -rowSums(Q0)
    I <- .symDiagonal(nrow(Q0))
    
    
    ### 6.9.2) Complile LGNB model
    Sys.setenv(PATH="%PATH%;C:/Rtools/gcc-4.6.3/bin;c:/Rtools/bin") #Run only when on windows
    setwd("~/FishCost/src/")
    compile("LGNB.cpp")
    dyn.load(dynlib("LGNB"))
    
    
    
    ### 6.9.3) Prepare TMB Data
    # TMB data are set in such way that it automatically
    # recognizes the data-specific inputs. 
    
    
    #### Linking R data to TMB data
    data <- list(
      time = if(TIME=="YearQuarter") {datatot$YearQuarter[rowid]} #Quarterly resolution
                 else if(TIME=="YearMonth"){datatot$YearMonth[rowid]}, #Monthly resolution
      gf = tmp2$gf,              
      rowid = rowid,
      response = datatot$Ntot,
      Q0 = Q0,
      I = I,
      Xpredict = matrix(0,0,0),
      Apredict = factor(numeric(0)),
      SupportAreaMatrix = SupportAreaMatrix,
      SupportAreaGroup = as.factor(datatot$split_area),
      Data=datatot$Data,
      h = mean(summary(as.polygons(gr))$side.length) #To be used later to plot the spatial decorrelation as a function of distance
    )
    
    
    
    
    # 6.10) Optimize and minimize the objective function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fit_model <- function(data, model_struct=NULL, with_static_field=FALSE, profile=TRUE) {
      time_levels <- levels(data$time) # There must be at least one time levels; They MUST match between the two data types!
      grid_nlevels <- nlevels(data$gf)
      data$doPredict <- 0
      data$offset <- model_struct$offset 
      
      if(FALSE) { ## FIXME: prediction disabled
        ## Stuff for prediction: Dummy dataset that matches the space time grid:
        ## Xpredict: design matrix for prediction
        ## Apredict: Area factor for prediction
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
        data$Xs                 <- model_struct$Xs
        data$X                  <- model_struct$Xf
        data$beta_r_fac         <- model_struct$beta_r_fac
        parameters$beta         <- model_struct$beta
        parameters$beta_r       <- model_struct$beta_r
        parameters$beta_s       <- model_struct$beta_s
        parameters$beta_r_logsd <- model_struct$beta_r_logsd
        data$which_beta_time    <- guess_time_effects_in_beta(data, time_levels)
      }
      
      ## Prior std dev on fixed effects (for robustness only)
      data$huge_sd <- 100
      
      map <- list()
      if(TRUE) map$logsd_nugget <- factor(NA)
      
      if(PS == "No" & DATA == "both"){
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(PS == "No" & DATA == "commercial") {
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(PS == "No" & DATA == "survey"){
        map$alpha <- factor(rep(NA,nlevels(data$SupportAreaGroup)))
      } else if(PS == "One"){
        lev <- levels(data$SupportAreaGroup)
        lev[lev=="survey"] <- NA
        map$alpha <- factor(lev)
      } else if(PS == "Two"){
        lev <- levels(data$SupportAreaGroup)
        #lev[lev=="survey"] <- NA
        map$alpha <- factor(lev)
        # map$alpha = factor(map$alpha)
      } 
      
      
      if(profile && (length(parameters$beta)>0 || length(parameters$beta)>0) ) {
        profile <- c("beta")
      } else {
        profile <- NULL
      }
      
      obj <- MakeADFun(data, parameters, random=c("eta_density","eta_nugget","eta_static","beta_r"),
                       profile = profile,
                       map=map, DLL="LGNB")
      
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
    
    
    
    # 6.11) Fit model
    #~~~~~~~~~~~~~~~~~~
    ## Reverse levels (to re-parameterize)
    data$Data <- factor(data$Data, rev(levels(data$Data))) #Make sure that survey data comes always first!
    datatot$Data <- factor(datatot$Data, rev(levels(datatot$Data))) #Make sure that survey data comes always first!
    
 
      if(DATA == "commercial"){
        
               if(TIME == "YearMonth"){
          m1 <- buildModelMatrices(fixed = ~  -1 + YearMonth, random = ~metiers-1 + VE_LENcat-1, offset = quote(log(HaulDur)), data=datatot)
        } else if(TIME == "YearQuarter"){
          m1 <- buildModelMatrices(fixed = ~  -1 + YearQuarter, random = ~ metiers-1 + VE_LENcat-1, offset = quote(log(HaulDur)), data=datatot)
        } 
        
      } else if(DATA == "survey"){
        
               if(TIME == "YearMonth"){
          m1 <- buildModelMatrices(fixed = ~ -1 + YearMonth, offset = quote(log(HaulDur)), data=datatot)
        } else if(TIME == "YearQuarter"){
          m1 <- buildModelMatrices(fixed = ~ -1 + YearQuarter, offset = quote(log(HaulDur)), data=datatot)
        } 
        
      } else {
        
                if(TIME == "YearMonth"){
          m1 <- buildModelMatrices(fixed = ~ -1 + YearMonth + Data, random = ~metiers-1 + VE_LENcat-1, offset= quote(log(HaulDur)), data=datatot)
        } else if (TIME == "YearQuarter"){
          m1 <- buildModelMatrices(fixed = ~ -1 + YearQuarter + Data, random = ~metiers-1 + VE_LENcat-1, offset= quote(log(HaulDur)), data=datatot)
        } 
      }
    
      env1 <- tryCatch(fit_model(data, m1, with_static_field = F), error=function(e) {})
    
    
    
    if(!is.null(env1)){
      
      
  
      # 6.12) Save results
      #~~~~~~~~~~~~~~~~~~~
      pl <- as.list(env1$sdr,"Estimate"); pl <- as.data.frame(pl$eta_density) #Estimates
      colnames(pl) <- paste(levels(env1$data$time),"est",sep=".")
      
      #pl.sd <- as.list(env1$sdr,"Std. Error"); pl.sd <- as.data.frame(pl.sd$eta_density) #Std. Error
      #colnames(pl.sd) <- paste(levels(env1$data$time),"se",sep=".")
      
      #dfsimu <- cbind(gr,pl,pl.sd)
      dfsimu <- cbind(gr,pl)
      
      
      SCENARIO_FULL <- gsub(",", ".", gsub("\\.", "", scenarios[s,5])) #s is the index for the scenario stated in the beginning of the script
      SCENARIO <- gsub(",",".", gsub("\\.","",paste(scenarios[s,1],scenarios[s,2],scenarios[s,3],sep="_")))

      #rm(list=setdiff(ls(), "dfsimu"))
      
      if(.Platform$OS.type == "windows") setwd("~/FishCost/Results/SimuAbundance")
      
      wd <- getwd()
      dir.create(paste0(wd,"/",DATA,"/",paste(DATA,SCENARIO,sep="_")),recursive=T)
      
      setwd(paste0(wd,"/",DATA,"/",paste(DATA,SCENARIO,sep="_")))
      OUTFILE  <- paste0("results_", paste(i,"scenario",SCENARIO_FULL,sep="_"), ".rds")
      saveRDS(dfsimu,file=OUTFILE)
      #rm(list=ls())
    }
  }
  
}

beep("mario")
