rm(list=ls())


library(ggplot2)
library(ggpubr)
library(purrr)
library(rgdal)
library(geosphere)


#~~~~~~~~~~~~~~~~~~
# 1) Load results
#~~~~~~~~~~~~~~~~~~


##############
##############
##    Cod   ##
##############
##############


############
#   2015   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [1] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results", SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_cod <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)




############
#   2016   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YEAR <- c("2015","2016")[2] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results",SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df16_cod <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Additional data manipulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.1) Define year
#~~~~~~~~~~~~~~~~~~
df15_cod$Year <- as.factor(paste("2015"))
df16_cod$Year <- as.factor(paste("2016"))


# 2.2) Define area of ICES polygon excluding WBS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("E:/OneDrive - DTU/PhD/GIS/ICES_area/")
ICES <- readOGR(dsn = ".", layer = "ICESareas_WB_and_Kattegat")
proj4string(ICES) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

ICES$MA_AREA
ICES$AR_WAT_KM2 # WB (ICES area 24 ((WB) has the largest area)

study.area16 <- sum(areaPolygon(ICES)[1:4]/1000000) #Include WBS area
study.area15 <- sum(areaPolygon(ICES)[2:4]/1000000) #Exclude WBS area as there are no survey sampling points in this area for the other species and years


df16_cod$StudyArea <- study.area16
df15_cod$StudyArea <- study.area15


# 2.3) Calculate area ratio between convexhull and study area 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_cod$AreaRatio <- df15_cod$ChullArea/df15_cod$StudyArea
df16_cod$AreaRatio <- df16_cod$ChullArea/df16_cod$StudyArea



# 2.4) Set scenario column
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_cod$COM <- gsub("//.*","",df15_cod$COM)
df15_cod$SURQ1 <- gsub("//.*","",df15_cod$SURQ1)
df15_cod$SURQ4 <- gsub("//.*","",df15_cod$SURQ4)

df15_cod$Scenario <- paste(df15_cod$COM, df15_cod$SURQ1, df15_cod$SURQ4,sep=":")


df16_cod$COM <- gsub("//.*","",df16_cod$COM)
df16_cod$SURQ1 <- gsub("//.*","",df16_cod$SURQ1)
df16_cod$SURQ4 <- gsub("//.*","",df16_cod$SURQ4)

df16_cod$Scenario <- paste(df16_cod$COM, df16_cod$SURQ1, df16_cod$SURQ4,sep=":")



# 2.5) Bind the two dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_cod <- rbind(df15_cod,df16_cod)



# 2.6) Reorder scenario levels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_cod$Scenario <- as.factor(df_full_cod$Scenario)
df_full_cod$Scenario <- ordered(df_full_cod$Scenario, levels=c("0:025:0","0:05:0","0:075:0","0:1:0", #survey SQ1
                                                               "0:0:025", "0:0:05", "0:0:075","0:0:1", #survey SQ4
                                                               
                                                               "0:025:025","0:025:05","0:025:075","0:025:1", #survey full
                                                               "0:05:025","0:05:05","0:05:075","0:05:1", #survey full
                                                               "0:075:025","0:075:05","0:075:075","0:075:1", #survey full
                                                               "0:1:025","0:1:05","0:1:075","0:1:1" #survey full
))






df_full_cod$dataID <- as.factor(ifelse(df_full_cod$COM==0 & df_full_cod$SURQ1==0,"SQ4",
                                       ifelse(df_full_cod$COM==0 & df_full_cod$SURQ4==0,"SQ1", "SQ1+SQ4")))



###############################################################################

#################
#################
##    Plaice   ##
#################
#################


############
#   2015   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [2] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results", SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_plaice <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)




############
#   2016   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YEAR <- c("2015","2016")[2] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results",SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df16_plaice <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Additional data manipulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.1) Define year
#~~~~~~~~~~~~~~~~~~
df15_plaice$Year <- as.factor(paste("2015"))
df16_plaice$Year <- as.factor(paste("2016"))


# 2.2) Define area of ICES polygon excluding WBS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

study.area <- sum(areaPolygon(ICES)[2:4]/1000000) #Exclude WBS area as there are no survey sampling points in this area for the other species and years

df15_plaice$StudyArea <- study.area
df16_plaice$StudyArea <- study.area


# 2.3) Calculate area ratio between convexhull and study area 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_plaice$AreaRatio <- df15_plaice$ChullArea/df15_plaice$StudyArea
df16_plaice$AreaRatio <- df16_plaice$ChullArea/df16_plaice$StudyArea



# 2.4) Set scenario column
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_plaice$COM <- gsub("//.*","",df15_plaice$COM)
df15_plaice$SURQ1 <- gsub("//.*","",df15_plaice$SURQ1)
df15_plaice$SURQ4 <- gsub("//.*","",df15_plaice$SURQ4)

df15_plaice$Scenario <- paste(df15_plaice$COM, df15_plaice$SURQ1, df15_plaice$SURQ4,sep=":")


df16_plaice$COM <- gsub("//.*","",df16_plaice$COM)
df16_plaice$SURQ1 <- gsub("//.*","",df16_plaice$SURQ1)
df16_plaice$SURQ4 <- gsub("//.*","",df16_plaice$SURQ4)

df16_plaice$Scenario <- paste(df16_plaice$COM, df16_plaice$SURQ1, df16_plaice$SURQ4,sep=":")



# 2.5) Bind the two dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_plaice <- rbind(df15_plaice,df16_plaice)



# 2.6) Reorder scenario levels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_plaice$Scenario <- as.factor(df_full_plaice$Scenario)
df_full_plaice$Scenario <- ordered(df_full_plaice$Scenario, levels=c("0:025:0","0:05:0","0:075:0","0:1:0", #survey SQ1
                                                                     "0:0:025", "0:0:05", "0:0:075","0:0:1", #survey SQ4
                                                                     
                                                                     "0:025:025","0:025:05","0:025:075","0:025:1", #survey full
                                                                     "0:05:025","0:05:05","0:05:075","0:05:1", #survey full
                                                                     "0:075:025","0:075:05","0:075:075","0:075:1", #survey full
                                                                     "0:1:025","0:1:05","0:1:075", "0:1:1" #survey full
))






df_full_plaice$dataID <- as.factor(ifelse(df_full_plaice$COM==0 & df_full_plaice$SURQ1==0,"SQ4",
                                          ifelse(df_full_plaice$COM==0 & df_full_plaice$SURQ4==0,"SQ1", "SQ1+SQ4")))





###############################################################################

##################
##################
##    Herring   ##
##################
##################


############
#   2015   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SPECIES <- c("cod","plaice","herring") [3] #Choose species from which results should be loaded; default is cod
YEAR <- c("2015","2016")[1] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results", SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_herring <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)




############
#   2016   #
############


# 1.1) Define default inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YEAR <- c("2015","2016")[2] #Choose year from wich simulation results should be loaded; default is 2015
DATA  <- c("both", "commercial", "survey") [3] #Choose data from which results should be loaded; default is both



# 1.2) Define working directory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Based on default inputs
setwd(paste("E:","OneDrive - DTU", "PhD", "Manuscript_02","CJFAS submission_02",
            "Review_n1","Results",SPECIES , YEAR, DATA, sep="/"))
wd <- getwd() #To be used later in the loop



# 1.3) Load all results into single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df16_herring <- list.files( path = wd, pattern = "*.rds", full.names = TRUE ) %>%
  map_dfr(readRDS)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Additional data manipulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.1) Define year
#~~~~~~~~~~~~~~~~~~
df15_herring$Year <- as.factor(paste("2015"))
df16_herring$Year <- as.factor(paste("2016"))


# 2.2) Define area of ICES polygon excluding WBS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

study.area <- sum(areaPolygon(ICES)[2:4]/1000000) #Exclude WBS area as there are no survey sampling points in this area for the other species and years

df15_herring$StudyArea <- study.area
df16_herring$StudyArea <- study.area


# 2.3) Calculate area ratio between convexhull and study area 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_herring$AreaRatio <- df15_herring$ChullArea/df15_herring$StudyArea
df16_herring$AreaRatio <- df16_herring$ChullArea/df16_herring$StudyArea



# 2.4) Set scenario column
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df15_herring$COM <- gsub("//.*","",df15_herring$COM)
df15_herring$SURQ1 <- gsub("//.*","",df15_herring$SURQ1)
df15_herring$SURQ4 <- gsub("//.*","",df15_herring$SURQ4)

df15_herring$Scenario <- paste(df15_herring$COM, df15_herring$SURQ1, df15_herring$SURQ4,sep=":")


df16_herring$COM <- gsub("//.*","",df16_herring$COM)
df16_herring$SURQ1 <- gsub("//.*","",df16_herring$SURQ1)
df16_herring$SURQ4 <- gsub("//.*","",df16_herring$SURQ4)

df16_herring$Scenario <- paste(df16_herring$COM, df16_herring$SURQ1, df16_herring$SURQ4,sep=":")



# 2.5) Bind the two dataset
#~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_herring <- rbind(df15_herring,df16_herring)



# 2.6) Reorder scenario levels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_herring$Scenario <- as.factor(df_full_herring$Scenario)
df_full_herring$Scenario <- ordered(df_full_herring$Scenario, levels=c("0:025:0","0:05:0","0:075:0","0:1:0", #survey SQ1
                                                                       "0:0:025", "0:0:05", "0:0:075","0:0:1", #survey SQ4
                                                                       
                                                                       "0:025:025","0:025:05","0:025:075","0:025:1", #survey full
                                                                       "0:05:025","0:05:05","0:05:075","0:05:1", #survey full
                                                                       "0:075:025","0:075:05","0:075:075","0:075:1", #survey full
                                                                       "0:1:025","0:1:05","0:1:075", "0:1:1" #survey full
))






df_full_herring$dataID <- as.factor(ifelse(df_full_herring$COM==0 & df_full_herring$SURQ1==0,"SQ4",
                                           ifelse(df_full_herring$COM==0 & df_full_herring$SURQ4==0,"SQ1", "SQ1+SQ4")))



###############################################################################



###############################
###############################
##                           ##
##          Plotting         ##
##                           ##
###############################
###############################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Some extra data manipulation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1.1) Include species-specific column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_full_cod$Species <- as.factor(paste("Cod"))
df_full_plaice$Species <- as.factor(paste("Plaice"))
df_full_herring$Species <- as.factor(paste("Herring"))


# 1.2) Bind all dfs into a single df
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_all <- rbind(df_full_cod, df_full_plaice, df_full_herring)

df_all$dataID <- factor(df_all$dataID, levels = c("SQ1", "SQ4", "SQ1+SQ4"))



#~~~~~~~~~~~~~~~~~~~~
# 2) Go for the plot
#~~~~~~~~~~~~~~~~~~~~


# jco_#B09C85FF: baige
# jco_#CD534CFF: red; 
# jco_#4A6990FF: dark blue; 
# jco_#DF8F44FF: orange
# jco_#79AF97FF: green
# jco_#7AA6DCFF: blue
# jco_#868686FF: gray

jco <- c("#4A6990FF","#7AA6DCFF","#B09C85FF")


baseline <- c("0:0:1","0:1:0","0:1:1")
df_all_base <- subset(df_all, Scenario %in% baseline)


df_all_base$Scenario <- factor(df_all_base$Scenario)





## All togehter
ggplot(df_all, aes(x=Scenario, y=AreaRatio)) + 
      geom_boxplot(aes(fill=dataID)) +
      geom_point(df_all_base, mapping=aes(x=Scenario, y=AreaRatio,fill=dataID), size=3, pch=21, col="black") +
  
      facet_grid(Species~Year)+
      scale_y_continuous(breaks=seq(0,1,0.1))+

      theme_pubclean()+
      scale_fill_manual(values=jco) + #could be a nice option, but reverting color order
      labs(x = "Scenario", y="Area ratio", fill="Data") +
      
      theme(#legend.position="none",  
        plot.margin = unit(c(1,2,1,1),"cm"),
        axis.text.x = element_text(size=13,angle = 90,vjust=0.5, hjust = 1,face="bold", margin=margin(6,0,0,0)),
        axis.text.y = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=17,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 15)),
        axis.title.x = element_text(size=17,face="bold", vjust=-4,margin = margin(t = 0, r = 20, b = 20, l = 0)),
        legend.position = "right",
        
        strip.background.y = element_rect(
          fill="gray96", size=1.5),
        
        strip.background.x = element_blank(),
        
        strip.text.y = element_text(size=17,face="bold"),
        strip.text.x = element_text(size=19,face="bold", margin = margin(0.5,0,0.5,0, "cm")),
        
        legend.title = element_text(size=15,face="bold"),
        legend.title.align=0.5,
        legend.text = element_text(size=13),
        legend.key = element_rect(color = NA, fill = NA),
        #legend.key.size = unit(1, "cm"),
        panel.border = element_rect(colour = "gray40", fill=NA, size=1.5))

  




setwd("E:/OneDrive - DTU/PhD/Manuscript_02/CJFAS submission_02/Review_n1/Results")
ggsave("Area_ratio_simulations.png", width = 30, height = 25, units = "cm", dpi=300)




## On a species level
#~~~~~~~~~~~~~~~~~~~~~~~~~~

# ggplot(df_full_cod, aes(x=Scenario, y=AreaRatio)) + 
#   geom_boxplot(aes(fill=dataID),notch=T) +
#   facet_grid(Year~.)+
#   scale_y_continuous(breaks=seq(0,1,0.1))+
#   #scale_y_continuous(breaks=seq(0,1.2,0.1))+
#   theme_pubclean() +
#   theme(#legend.position="none",  
#     plot.margin = unit(c(1,2,1,1),"cm"),
#     axis.text.x = element_text(size=30,angle = 90,vjust=0.5, hjust = 1,face="bold", margin=margin(15,0,0,0)),
#     axis.text.y = element_text(size=30,face="bold"),
#     axis.title.y = element_text(size=34,face="bold", margin = margin(t = 0, r = 20, b = 0, l = 15)),
#     axis.title.x = element_text(size=34,face="bold", vjust=-4,margin = margin(t = 0, r = 20, b = 20, l = 0)),
#     #legend.position = "none",
#     
#     # legend.position = "right",
#     # legend.title = element_blank(),
#     # legend.text = element_text(size=28),
#     # legend.key = element_rect(color = NA, fill = NA),
#     # legend.key.size = unit(2, "cm"),
#     panel.border = element_rect(colour = "black", fill=NA, size=1.5))



