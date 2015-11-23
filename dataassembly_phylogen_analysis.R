#phylogenetic distance project
#assembling the seedling transplant experiment dataset (longformat of all variables)

#set up work environment
library(ggplot2)
library(lubridate)
library(reshape2)
library(plyr)
library(dplyr)
library(rjags)
library(multcompView)
library(readr)

##########
  rm(list = ls())
  

#seedling data  
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/seedling data")
    seedling <- read.csv("allseedlings150528.csv")
    seedling <- seedling[-13428,]  #remove the last row (sort test row)    
 
#excluding transplant shock
  seedling <- subset(seedling, timealive > 14)  #not including seedlings that may have died within the first 2 weeks

#size of database
  N <- nrow(seedling)

##########
#data cleaning
##########

#standardizing censuses (removing 2010 and 2011 extra censuses from certain locations so all years/regions have the same number of censuses)
    #survival
      #2010 C3
      seedling$y.2010.C3.alive[seedling$region== "ann arbor"]<-seedling$y.2010.C5.alive[seedling$region== "ann arbor"]
      
      #2011 C1
      seedling$y.2011.C1.alive[seedling$region== "ann arbor"]<-seedling$y.2011.C2.alive[seedling$region== "ann arbor"]
      
      #2011 C2
      seedling$y.2011.C2.alive[seedling$region== "ann arbor"]<-seedling$y.2011.C3.alive[seedling$region== "ann arbor"]
      seedling$y.2011.C2.alive[seedling$region== "umbs"]<-seedling$y.2011.C3.alive[seedling$region== "umbs"]
      
      #2011 C3
      seedling$y.2011.C3.alive[seedling$region== "ann arbor"]<-seedling$y.2011.C5.alive[seedling$region== "ann arbor"]
      seedling$y.2011.C3.alive[seedling$region== "umbs"]<-seedling$y.2011.C5.alive[seedling$region== "umbs"]
      
    #dates
      seedling[,96:119] <- lapply(seedling[,96:119], as.character)   #making sure that the date columns aren't factors
    
      seedling$y.2010.C3.date[seedling$region== "ann arbor"] <- seedling$y.2010.C5.date[seedling$region== "ann arbor"]
      
      #2011 C1
      seedling$y.2011.C1.date[seedling$region== "ann arbor"]<-seedling$y.2011.C2.date[seedling$region== "ann arbor"]
      
      #2011 C2
      seedling$y.2011.C2.date[seedling$region== "ann arbor"]<-seedling$y.2011.C3.date[seedling$region== "ann arbor"]
      seedling$y.2011.C2.date[seedling$region== "umbs"]<-seedling$y.2011.C3.date[seedling$region== "umbs"]
      
      #2011 C3
      seedling$y.2011.C3.date[seedling$region== "ann arbor"]<-seedling$y.2011.C5.date[seedling$region== "ann arbor"]
      seedling$y.2011.C3.date[seedling$region== "umbs"]<-seedling$y.2011.C5.date[seedling$region== "umbs"]

    #heights
      #2010
      seedling$y.2010.C1.height[seedling$region== "ann arbor"] <- seedling$y.2010.C2.height[seedling$region== "ann arbor"]
      seedling$y.2010.height_fall <- seedling$y.2010.C5.height
      seedling$y.2010.height_fall[seedling$region== "umbs"] <- seedling$y.2010.C3.height[seedling$region== "umbs"]
      
      #2011
      seedling$y.2011.height_fall <- seedling$y.2011.C6.height
      seedling$y.2011.height_fall[seedling$region== "umbs"] <- seedling$y.2011.C4.height[seedling$region== "umbs"]

    #diameters      
      #2011
      seedling$y.2011.diameter_fall <- seedling$y.2011.C6.diameter
      seedling$y.2011.diameter_fall[seedling$region== "umbs"] <- seedling$y.2011.C4.diameter[seedling$region== "umbs"]


    #herbiv
      #2010 C3
      seedling$y.2010.C3.h.avg[seedling$region== "ann arbor"]<-seedling$y.2010.C5.h.avg[seedling$region== "ann arbor"]
      
      #2011 C1
      seedling$y.2011.C1.h.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C2.h.avg[seedling$region== "ann arbor"]
      
      #2011 C2
      seedling$y.2011.C2.h.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C3.h.avg[seedling$region== "ann arbor"]
      seedling$y.2011.C2.h.avg[seedling$region== "umbs"]<-seedling$y.2011.C3.h.avg[seedling$region== "umbs"]
      
      #2011 C3
      seedling$y.2011.C3.h.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C5.h.avg[seedling$region== "ann arbor"]
      seedling$y.2011.C3.h.avg[seedling$region== "umbs"]<-seedling$y.2011.C5.h.avg[seedling$region== "umbs"]
      
    #pathogens/disease
      #2010 C3
      seedling$y.2010.C3.p.avg[seedling$region== "ann arbor"]<-seedling$y.2010.C5.p.avg[seedling$region== "ann arbor"]
      
      #2011 C1
      seedling$y.2011.C1.p.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C2.p.avg[seedling$region== "ann arbor"]
      
      #2011 C2
      seedling$y.2011.C2.p.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C3.p.avg[seedling$region== "ann arbor"]
      seedling$y.2011.C2.p.avg[seedling$region== "umbs"]<-seedling$y.2011.C3.p.avg[seedling$region== "umbs"]
      
      #2011 C3
      seedling$y.2011.C3.p.avg[seedling$region== "ann arbor"]<-seedling$y.2011.C5.p.avg[seedling$region== "ann arbor"]
      seedling$y.2011.C3.p.avg[seedling$region== "umbs"]<-seedling$y.2011.C5.h.avg[seedling$region== "umbs"]
  

    #number of leaves
      #2010 C3
      seedling$y.2010.C3.nleaves[seedling$region== "ann arbor"]<-seedling$y.2010.C5.nleaves[seedling$region== "ann arbor"]
      
      #2011 C1
      seedling$y.2011.C1.nleaves[seedling$region== "ann arbor"]<-seedling$y.2011.C2.nleaves[seedling$region== "ann arbor"]
      
      #2011 C2
      seedling$y.2011.C2.nleaves[seedling$region== "ann arbor"]<-seedling$y.2011.C3.nleaves[seedling$region== "ann arbor"]
      seedling$y.2011.C2.nleaves[seedling$region== "umbs"]<-seedling$y.2011.C3.nleaves[seedling$region== "umbs"]
      
      #2011 C3
      seedling$y.2011.C3.nleaves[seedling$region== "ann arbor"]<-seedling$y.2011.C5.nleaves[seedling$region== "ann arbor"]
      seedling$y.2011.C3.nleaves[seedling$region== "umbs"]<-seedling$y.2011.C5.nleaves[seedling$region== "umbs"]


    #mammal browse  #standardizing census dates
      seedling$mam_clean <- as.character(seedling$mambrowse_yr_cen)

        #2010 C3
          seedling$mam_clean[seedling$region == "ann arbor" & seedling$mam_clean[] =="y.2010.C5.date"] <- "y.2010.C3.date"

        #2011 C1
            seedling$mam_clean[seedling$region == "ann arbor" & seedling$mam_clean[] =="y.2011.C2.date"] <- "y.2011.C1.date"

        #2011 C2
            seedling$mam_clean[seedling$region == "ann arbor" & seedling$mam_clean[] =="y.2011.C3.date"] <- "y.2011.C2.date"
            seedling$mam_clean[seedling$region == "umbs" & seedling$mam_clean[] =="y.2011.C3.date"] <- "y.2011.C2.date"
        #2011 C3
            seedling$mam_clean[seedling$region == "ann arbor" & seedling$mam_clean[] =="y.2011.C5.date"] <- "y.2011.C3.date"
            seedling$mam_clean[seedling$region == "umbs" & seedling$mam_clean[] =="y.2011.C5.date"] <- "y.2011.C3.date"

        #other missing bits
            seedling$mam_clean[seedling$mam_clean[] =="y.2011.C4.date"] <- "y.2011.C2.date"
            seedling$mam_clean[seedling$mam_clean[] =="y.2011.C6.date"] <- "y.2011.C3.date"
            seedling$mam_clean[seedling$mam_clean[] =="y.2011.C7.date"] <- "y.2011.C3.date"

#############
#survival data
#############

   #assembling dataframes for time delimited variables
       survmat<-seedling[,c("y.2010.C1.alive","y.2010.C2.alive","y.2010.C3.alive",
                           "y.2011.C1.alive","y.2011.C2.alive","y.2011.C3.alive",
                           "y.2012.C1.alive","y.2012.C2.alive","y.2012.C3.alive",
                           "y.2013.C1.alive","y.2013.C2.alive","y.2013.C3.alive",
                           "y.2014.C1.alive","y.2014.C2.alive","y.2014.C3.alive")]

#########
#dates and ages
#########
    datemat<-seedling[,c(
                  "y.2010.C1.date","y.2010.C2.date","y.2010.C3.date",
                  "y.2011.C1.date","y.2011.C2.date","y.2011.C3.date",
                  "y.2012.C1.date","y.2012.C2.date","y.2012.C3.date",
                  "y.2013.C1.date","y.2013.C2.date","y.2013.C3.date",
                  "y.2014.C1.date","y.2014.C2.date","y.2014.C3.date")]

    #year matrix
        yearmat <- matrix( rep(2010:2014,each = nrow(datemat)*3),
          nrow= nrow(datemat), ncol = (ncol(datemat)))

    #season of census 0 = spring, 1 = summer 2 = fall
        smat <- data.frame(matrix(NA, nrow=N, ncol=15))
        smat[,c(1,4,7,10,13)] <- 0 
        smat[,c(2,5,8,11,14)] <- 1
        smat[,c(3,6,9,12,15)] <- 2

    #planting season
         seedling$approx_plantingdate_ch <- as.character(seedling$approx_plantingdate)
                    seedling$approx_plantingdate2 <- mdy(seedling$approx_plantingdate_ch)
                    seedling$planting_season <- month(seedling$approx_plantingdate2)
                    seedling$planting_season[seedling$planting_season[] < 8] <- 1
                    seedling$planting_season[seedling$planting_season[] > 7] <- 2
  
    #agedaymat
        agedaymat <- datemat
        agedaymat[] <- lapply(datemat,as.character)
        agedaymat[] <- lapply(agedaymat,mdy)
        
        agedaymat2 <- matrix(NA, nrow=N, ncol=15)
        for(i in 1:ncol(agedaymat2)){
          agedaymat2[,i] <- difftime(agedaymat[,i],seedling$approx_plantingdate2, units="days")
        }

    #agemat seedlings that are first vs. 2+ years 
        seedling$yearplanted <- as.numeric(year(seedling$approx_plantingdate2))
        agematyear <- yearmat[,1:15] - seedling$yearplanted +1
        agematyear[agematyear[] < 1] <- NA
        agematyear1vs2 <- agematyear
        agematyear1vs2[agematyear1vs2[]>1] <- 2

#########
#height
########
        
  #height matrix
    seedling$offset <- NA
    heightmat <- seedling[,c("y.2010.C1.height","offset","y.2010.height_fall",
                          "y.2011.C1.height","offset","y.2011.height_fall",
                          "y.2012.C1.height","offset","y.2012.C3.height",
                          "y.2013.C1.height","offset","y.2013.C3.height",
                          "y.2014.C1.height","offset","y.2014.C3.height")]
    heightmat[ ,2] <- (heightmat[ ,1] + heightmat[ ,3] )/2  #mid season census as mean of spring and fall census
    heightmat[ ,5] <- (heightmat[ ,4] + heightmat[ ,6] )/2
    heightmat[ ,8] <- (heightmat[ ,7] + heightmat[ ,9] )/2
    heightmat[ ,11] <- (heightmat[ ,10] + heightmat[ ,11] )/2
    heightmat[ ,14] <- (heightmat[ ,13] + heightmat[ ,15] )/2

  #difference in height
    difheight <- data.frame(matrix(NA,nrow=nrow(seedling),ncol=4))
    names(difheight) <- c("dh2011","dh2012", "dh2013", "dh2014")
    
    difheight$dh2011 <- seedling$y.2011.height_fall - seedling$y.2010.height_fall
    difheight$dh2012 <- seedling$y.2012.C3.height - seedling$y.2011.height_fall
    difheight$dh2013 <- seedling$y.2013.C3.height - seedling$y.2012.C3.height
    difheight$dh2014 <- seedling$y.2014.C3.height - seedling$y.2013.C3.height   

    #adding in marked heights
    #compiling height marked
      difheight2013_marked <- seedling$y.2013.C3.height_marked - seedling$Y.2012.C1.height_marked
      difheight2014_marked <- seedling$y.2014.C3.height_marked - seedling$Y.2013.C3.height_marked
      
      for(i in 1:nrow(seedling)){
        if(!is.na(difheight2013_marked[i])){ difheight$dh2013[i] <- difheight2013_marked[i] }
        if(!is.na(difheight2014_marked[i])){ difheight$dh2014[i] <- difheight2014_marked[i] }
      }
    
  #adding in heights for fall plantings
    difheight2012_fallplant <- seedling$y.2012.C3.height - seedling$y.2012.C1.height
    difheight2013_fallplant <- seedling$y.2013.C3.height - seedling$y.2013.C1.height
    difheight2014_fallplant <- seedling$y.2014.C3.height - seedling$y.2014.C1.height
    
    difheight$dh2012[seedling$y.planted == 2011 & seedling$planting_season == 2]  <- 
      difheight2012_fallplant[seedling$y.planted == 2011 & seedling$planting_season == 2] 
    
    difheight$dh2013[seedling$y.planted == 2012 & seedling$planting_season == 2]  <- 
      difheight2013_fallplant[seedling$y.planted == 2012 & seedling$planting_season == 2] 
    
    difheight$dh2014[seedling$y.planted == 2013 & seedling$planting_season == 2]  <- 
      difheight2014_fallplant[seedling$y.planted == 2013 & seedling$planting_season == 2] 
    
  #making a version that will be compatible with the censuses
    difheightb <- data.frame(matrix(NA,nrow = nrow(seedling),ncol=15))
    difheightb[,4:6] <- difheight$dh2011  
    difheightb[,7:9] <- difheight$dh2012  
    difheightb[,10:12] <- difheight$dh2013  
    difheightb[,13:15] <- difheight$dh2014  

#########
#diameter
########
    
 #adding in marked diameters for seedlings where there wasn't an unmarked version
      seedling$y.2013.C3.diameter[is.na(seedling$y.2013.C3.diameter) & !is.na(seedling$y.2013.C3.diameter_marked)] <- 
        seedling$y.2013.C3.diameter_marked[is.na(seedling$y.2013.C3.diameter) & !is.na(seedling$y.2013.C3.diameter_marked)] 

      seedling$y.2014.C1.diameter <- as.numeric(as.character(seedling$y.2014.C1.diameter))
      seedling$y.2014.C1.diameter[is.na(seedling$y.2014.C1.diameter) & !is.na(seedling$y.2014.C1.diameter_marked)] <- 
        seedling$y.2014.C1.diameter_marked[is.na(seedling$y.2014.C1.diameter) & !is.na(seedling$y.2014.C1.diameter_marked)] 

  #diameter matrix
  diammat <- seedling[,c("y.2010.C2.diameter","offset","y.2010.C7.diameter",
                           "y.2011.C1.diameter","offset","y.2011.diameter_fall",
                           "y.2012.C1.diameter","offset","y.2012.C3.diameter",
                           "y.2013.C1.diameter","offset","y.2013.C3.diameter",
                           "y.2014.C1.diameter","offset","y.2014.C3.diameter")]
  diammat[ ,2] <- (diammat[ ,1] + diammat[ ,3] )/2  #mid season census as mean of spring and fall census
  diammat[ ,5] <- (diammat[ ,4] + diammat[ ,6] )/2
  diammat[ ,8] <- (diammat[ ,7] + diammat[ ,9] )/2
  diammat[ ,11] <- (diammat[ ,10] + diammat[ ,11] )/2
  diammat[ ,14] <- (diammat[ ,13] + diammat[ ,15] )/2
  
  #difference in diam
  difdiam <- data.frame(matrix(NA,nrow=nrow(seedling),ncol=4))
  names(difdiam) <- c("dh2011","dh2012", "dh2013", "dh2014")
  
  difdiam$dh2011 <- seedling$y.2011.diameter_fall - seedling$y.2010.C7.diameter
  difdiam$dh2012 <- seedling$y.2012.C3.diameter - seedling$y.2011.diameter_fall
  difdiam$dh2013 <- seedling$y.2013.C3.diameter - seedling$y.2012.C3.diameter
  difdiam$dh2014 <- seedling$y.2014.C3.diameter - seedling$y.2013.C3.diameter
    
  #adding in diams for fall plantings
  difdiam2012_fallplant <- seedling$y.2012.C3.diameter - seedling$y.2012.C1.diameter
  difdiam2013_fallplant <- seedling$y.2013.C3.diameter - seedling$y.2013.C1.diameter
  difdiam2014_fallplant <- seedling$y.2014.C3.diameter - seedling$y.2014.C1.diameter
  
  difdiam$dh2012[seedling$y.planted == 2011 & seedling$planting_season == 2]  <- 
    difdiam2012_fallplant[seedling$y.planted == 2011 & seedling$planting_season == 2] 
  
  difdiam$dh2013[seedling$y.planted == 2012 & seedling$planting_season == 2]  <- 
    difdiam2013_fallplant[seedling$y.planted == 2012 & seedling$planting_season == 2] 
  
  difdiam$dh2014[seedling$y.planted == 2013 & seedling$planting_season == 2]  <- 
    difdiam2014_fallplant[seedling$y.planted == 2013 & seedling$planting_season == 2] 
  
  
  #making a version that will be compatible with the censuses
  difdiamb <- data.frame(matrix(NA,nrow = nrow(seedling),ncol=15))
  difdiamb[,4:6] <- difdiam$dh2011  
  difdiamb[,7:9] <- difdiam$dh2012  
  difdiamb[,10:12] <- difdiam$dh2013  
  difdiamb[,13:15] <- difdiam$dh2014  


#########
#treatment
#########

#treatment; making a matrix
  tmtmat <- data.frame(matrix(0, nrow= N, ncol = 15))
  tmtmat[,1:3] <- seedling$tmt2010
  tmtmat[,4:6] <- seedling$tmt2011
  tmtmat[,7:9] <- seedling$tmt2012
  tmtmat[,10:12] <- seedling$tmt2013
  tmtmat[,13:15] <- seedling$tmt2014

    tmtmat[tmtmat[]=="control"] <- 1
    tmtmat[tmtmat[]=="partial exclosure"] <- 2
    tmtmat[tmtmat[]=="full exclosure"] <- 3
    tmtmat[tmtmat[]=="pesticide"] <- 4

#########
#herbivory
#########
    seedling$y.2014.C1.h.avg<-NA  #not measured, but putting in a placeholder
    seedling$y.2014.C3.h.avg<-NA  #not measured, but putting in a placeholder
    hmat<-seedling[,c(
                  "y.2010.C1.h.avg","y.2010.C2.h.avg","y.2010.C3.h.avg",
                  "y.2011.C1.h.avg","y.2011.C2.h.avg","y.2011.C3.h.avg",
                  "y.2012.C1.h.avg","y.2012.C2.h.avg","y.2012.C3.h.avg",
                  "y.2013.C1.h.avg","y.2013.C2.h.avg","y.2013.C3.h.avg",
                  "y.2014.C1.h.avg","y.2014.C2.h.avg","y.2014.C3.h.avg")]

      #hmat # hmat <- as.matrix(hmat)   hmat <- hmat*.01  #proportion instead of percent
          hmat2<-unlist(hmat); hist(hmat2)
          hmat3 <- hmat2 * 0.01
          hmat4 <- hmat
          hmat4[]<-hmat3
          hmat4[hmat4[] == 0.001] <- 0



#########
#pathogens
#########
    seedling$y.2014.C1.p.avg <- NA  #not measured, but putting in a placeholder
    seedling$y.2014.C3.p.avg <- NA  #not measured, but putting in a placeholder
    seedling$y.2013.C1.p.avg <- NA  #was too early of a census to catch p
    
    pmat<-seedling[,c(
                      "y.2010.C1.p.avg","y.2010.C2.p.avg","y.2010.C3.p.avg",
                      "y.2011.C1.p.avg","y.2011.C2.p.avg","y.2011.C3.p.avg",
                      "y.2012.C1.p.avg","y.2012.C2.p.avg","y.2012.C3.p.avg",
                      "y.2013.C1.p.avg","y.2013.C2.p.avg","y.2013.C3.p.avg",
                      "y.2014.C1.p.avg","y.2014.C2.p.avg","y.2014.C3.p.avg")]
    
    #pmat # pmat <- as.matrix(pmat)   pmat <- pmat*.01  #proportion instead of percent
    pmat2<-unlist(pmat); hist(pmat2)
    pmat3 <- pmat2 * 0.01
    pmat4 <- pmat
    pmat4[]<-pmat3
    pmat4[pmat4[] == 0.001] <- 0



############
#mammal browse
############
      seedling$mam_cnumber <-  seedling$mam_clean
      datenames <- names(datemat)

      for(i in 1:length(datenames)){
      seedling$mam_cnumber[seedling$mam_cnumber==datenames[i]] <- which(datenames == datenames[i])
      }
      
      unique(seedling$mam_cnumber)
      seedling$mam_cnumber <- as.numeric(seedling$mam_cnumber)

      #creating matrix to house results
          browsemat <- data.frame(matrix(0,nrow(datemat),ncol(datemat)))
          names(browsemat) <- names(datemat)

      #populating browsemat
          for(i in 1:length(seedling$mam_cnumber)){
            if(!is.na(seedling$mam_cnumber[i])){browsemat[i,seedling$mam_cnumber[i]] <- 1}           
          }

      #Get the NA's into browsemat, using blank cells from datemat
                #getting NAs from datemat
                datemat2_char <- datemat
                datemat2_char[] <- lapply(datemat, as.character)  #switching it to character
                datemat2_char[datemat2_char[]==""] <- NA

            #the more elegant solution didn't work, so just going through it with for loops to install NAs
                for(i in 1:nrow(datemat)){
                  for(j in 1:ncol(datemat)){
                    if(is.na(datemat2_char[i,j])){browsemat[i,j] <- NA}             
                  }
                }



############
#environmental covariates
############
#soil nutrients
  #adding in soil data
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/environmental data/resin capsules and nutrient analysis")
    resin<-read.csv("resindata140501.csv")
    resin<-subset(resin,PI=="dk")
  
  #average per plot
      nutr<-ddply(resin,c("plot"),summarise,
                  Ntotal=mean(na.omit(TotalN)),NO3_N=mean(na.omit(NO3_N)), NH4_N=mean(na.omit(NH4_N)),
                  P=mean(na.omit(P)),Ca=mean(na.omit(Ca)),Mg=mean(na.omit(Mg)),Mn=mean(na.omit(Mn)))
      
      seedling <- join(seedling, nutr, by = "plot", type ="left")

  #light
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/environmental data/canopy photos")
    light<-read.csv("plotlight.csv")

    seedling$plot.letter <- tolower(as.character(seedling$plotletter))
    light$plot.letter <- tolower(as.character(light$plot.letter))

    seedling <- join(seedling,light, by = c("plot","plot.letter"), type = "left")

  #soil moisture
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/environmental data/soil moisture")
    soilm <- read.csv("soilmoisture_census_plotsection150325.csv", stringsAsFactors = FALSE)  #this was produced by soilmoistureintegration150325.R    
    soilm$y.2012.C3.date[soilm$y.2012.C3.date < 0] <- 0  #setting all negative values as 0
    soilm$plot.letter <- tolower(as.character(soilm$plotletter))

    plotsm <- subset(soilm, select =-c(X, region ,site,plot, plotsection, plotletter, plot.letter))  #getting rid of non-date columns
    plotsm2<-unlist(plotsm); hist(plotsm2)
    plotsm3<-(plotsm2-mean(na.omit(plotsm2)))/(2*sd(na.omit(plotsm2))); hist(plotsm3)
    plotsm4<-plotsm
    plotsm4[]<-plotsm3
          #plotsm4[,16] <- NA
          #plotsm5 <- plotsm4
          #plotsm5[is.na(plotsm5)]<-0
  
    #required later for sm integration
      seedling$plotsection <- NA
      seedling$plotsection[seedling$row <6 ] <- "front"
      seedling$plotsection[seedling$row >5 & seedling$row < 11] <- "middle"
      seedling$plotsection[seedling$row >10 ] <- "back"

    #getting a version that lines up with seedling
      soilm_seedling <- seedling[, c("plot", "plot.letter", "plotsection","col","row")]
      soilm_seedling <- join(soilm_seedling, soilm, by = c("plot", "plot.letter", "plotsection"),type= "left")

      #soil moisture sd
        soilm_sd <- read.csv("soilmoisture_census_plotsection_sd150325.csv", stringsAsFactors = FALSE)  #this was produced by soilmoistureintegration150325.R    
        soilm_sd$plot.letter <- tolower(as.character(soilm_sd$plotletter))

      plotsm_sd <- subset(soilm_sd, select =-c(X, region ,site,plot, plotsection, plotletter, plot.letter))  #getting rid of non-date columns
      
      #getting a version that lines up with seedling
      soilm_sd_seedling <- seedling[, c("plot", "plot.letter", "plotsection","col","row")]
      soilm_sd_seedling <- join(soilm_sd_seedling, soilm_sd, by = c("plot", "plot.letter", "plotsection"),type= "left")
      
      

############
# seedling: nleaves
############

  nleaves <- seedling[,c("y.2010.C1.nleaves","y.2010.C2.nleaves","y.2010.C3.nleaves",
                    "y.2011.C1.nleaves","y.2011.C2.nleaves","y.2011.C3.nleaves",
                    "y.2012.C1.nleaves","y.2012.C2.nleaves","y.2012.C3.nleaves",
                    "y.2013.C1.nleaves","y.2013.C2.nleaves","y.2013.C3.nleaves",
                    "offset","y.2014.C2.nleaves","offset")]
  
      nleaves[survmat[,] == 0] <- NA    #turning dead seedling nleaves to NA
      nleaves[is.na(survmat[,])] <- NA  #turning any missing data to NA



###########
#stem mapping data #neigh
##########

  #adding in seedling coordinates, neighborhood index, plot basal area
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/chapter 2 neighborhood analysis")
    ni <- read.csv("seedling_ni_ba_d_10m150319.csv")
        ni$plot.letter <- tolower(as.character(ni$plotletter))
        ni$col_order <- ni$col
    seedling <- join(seedling, ni, by = c("plot","plot.letter","col_order","row"), type ="left" )

  #conspecific present in site?
    setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/stem maps")
    stems <- read.csv("stemmaps141014.csv")
    
    stems$treeba <- (stems$tree.DBH^2)*0.00007854  #calculating ba from dbh
    
    stemssite <- ddply(stems,c("site","tree.sp."),summarise,  #getting the number and ba of each species of tree at each site
                          ba = sum(treeba),
                          nobs = length(na.omit(distance)))

    stemssitetot_ba <- ddply(stemssite, "site", summarise,  #getting total ba at each site
                             tot_ba = sum(ba))
    
    stemssite <- join(stemssite, stemssitetot_ba, by = "site", type="left")  #merging in total ba
    stemssite$relba <- stemssite$ba / stemssite$tot_ba  #relative ba
    
       stemssite <- subset(stemssite, tree.sp. == "Acru" | tree.sp. == "Cagl" | tree.sp. =="Litu" | tree.sp. =="Qual" | tree.sp. =="Quru" |
                                   tree.sp. =="Quve" |  tree.sp. =="Rops")  #only including study species

    stemssite$sp <- tolower(stemssite$tree.sp.)  #converting tree.sp. to lower case
    stemssite$conspecific_relba_site <- stemssite$relba

    stemssite <- stemssite[,c("site", "sp", "conspecific_relba_site")]  #cutting down dataframe to just the parts required to merge it in to seedling
    
    seedling <- join(seedling, stemssite, by = c("site","sp"), type = "left")
    seedling$conspecific_relba_site[is.na(seedling$conspecific_relba_site)] <- 0  
    seedling$conspecific_at_site  <- "not present"
    seedling$conspecific_at_site[seedling$conspecific_relba_site > 0.0001] <- "present"

    #conspecific ni
      seedling$conspecific_ni <- NA
      seedling$conspecific_ni[seedling$sp == "acru"] <- seedling$ni_acru[seedling$sp == "acru"]
      seedling$conspecific_ni[seedling$sp == "cagl"] <- seedling$ni_casp[seedling$sp == "cagl"]
      seedling$conspecific_ni[seedling$sp == "litu"] <- seedling$ni_litu[seedling$sp == "litu"]
      seedling$conspecific_ni[seedling$sp == "qual"] <- seedling$ni_qual[seedling$sp == "qual"]
      seedling$conspecific_ni[seedling$sp == "quru"] <- seedling$ni_quru[seedling$sp == "quru"]
      seedling$conspecific_ni[seedling$sp == "quve"] <- seedling$ni_quve[seedling$sp == "quve"]
      seedling$conspecific_ni[seedling$sp == "rops"] <- seedling$ni_rops[seedling$sp == "rops"]
      seedling$conspecific_ni[seedling$sp == "beth"] <- 0
      seedling$conspecific_ni[seedling$sp == "ceor"] <- 0
      seedling$conspecific_ni[seedling$sp == "elum"] <- 0
      
    #conspecific d
      seedling$conspecific_d <- NA
      seedling$conspecific_d[seedling$sp == "acru"] <- seedling$d_acru[seedling$sp == "acru"]
      seedling$conspecific_d[seedling$sp == "cagl"] <- seedling$d_casp[seedling$sp == "cagl"]
      seedling$conspecific_d[seedling$sp == "litu"] <- seedling$d_litu[seedling$sp == "litu"]
      seedling$conspecific_d[seedling$sp == "qual"] <- seedling$d_qual[seedling$sp == "qual"]
      seedling$conspecific_d[seedling$sp == "quru"] <- seedling$d_quru[seedling$sp == "quru"]
      seedling$conspecific_d[seedling$sp == "quve"] <- seedling$d_quve[seedling$sp == "quve"]
      seedling$conspecific_d[seedling$sp == "rops"] <- seedling$d_rops[seedling$sp == "rops"]
      seedling$conspecific_d[seedling$sp == "beth"] <- 0
      seedling$conspecific_d[seedling$sp == "ceor"] <- 0
      seedling$conspecific_d[seedling$sp == "elum"] <- 0
    
    #conspecific ba
      seedling$conspecific_ba <- NA
      seedling$conspecific_ba[seedling$sp == "acru"] <- seedling$ba_acru[seedling$sp == "acru"]
      seedling$conspecific_ba[seedling$sp == "cagl"] <- seedling$ba_casp[seedling$sp == "cagl"]
      seedling$conspecific_ba[seedling$sp == "litu"] <- seedling$ba_litu[seedling$sp == "litu"]
      seedling$conspecific_ba[seedling$sp == "qual"] <- seedling$ba_qual[seedling$sp == "qual"]
      seedling$conspecific_ba[seedling$sp == "quru"] <- seedling$ba_quru[seedling$sp == "quru"]
      seedling$conspecific_ba[seedling$sp == "quve"] <- seedling$ba_quve[seedling$sp == "quve"]
      seedling$conspecific_ba[seedling$sp == "rops"] <- seedling$ba_rops[seedling$sp == "rops"]
      seedling$conspecific_ba[seedling$sp == "beth"] <- 0
      seedling$conspecific_ba[seedling$sp == "ceor"] <- 0
      seedling$conspecific_ba[seedling$sp == "elum"] <- 0
  
    #congeneric ba
      seedling$congeneric_ba <- NA
      seedling$congeneric_ba[seedling$sp == "acru"] <-  seedling$ba_acsp[seedling$sp == "acru"] - seedling$ba_acru[seedling$sp == "acru"]
      seedling$congeneric_ba[seedling$sp == "cagl"] <-  0  #I don't trust MR, SW, MD to distinguish between sp.  
      seedling$congeneric_ba[seedling$sp == "litu"] <-  0 #no congenerics
      seedling$congeneric_ba[seedling$sp == "rops"] <-  0 #no congenerics
      seedling$congeneric_ba[seedling$sp == "qual"] <-  seedling$ba_qusp[seedling$sp == "qual"] - seedling$ba_qual[seedling$sp == "qual"]
      seedling$congeneric_ba[seedling$sp == "quru"] <-  seedling$ba_qusp[seedling$sp == "quru"] - seedling$ba_quru[seedling$sp == "quru"]
      seedling$congeneric_ba[seedling$sp == "quve"] <-  seedling$ba_qusp[seedling$sp == "quve"] - seedling$ba_quve[seedling$sp == "quve"]
      seedling$congeneric_ba[seedling$sp == "elum"] <-  0
      seedling$congeneric_ba[seedling$sp == "beth"] <-  0
      seedling$congeneric_ba[seedling$sp == "ceor"] <-  0


#######
#conspecific density at multiple scales
#######
  seedling$multiscaled <- NA
  seedling$multiscaled[seedling$migstat == "exotic"] <- "exotic"
  seedling$multiscaled[seedling$migstat == "migrant"] <- "conspecific not in region"
  seedling$multiscaled[seedling$migstat == "native"] <- "conspecific in region"
  seedling$multiscaled[seedling$conspecific_at_site == "present"] <- "conspecific present at site"
  seedling$multiscaled[seedling$conspecific_ba > 0] <- "conspecific within 10 m"
  
  seedling$multiscaled <- factor(seedling$multiscaled, levels = c("conspecific within 10 m",
                                                                "conspecific present at site","conspecific in region", 
                                                                "conspecific not in region", "exotic"))
  


###########
#assembling a long format dataset 
############
seedling$seedid <- as.numeric(row.names(seedling))

#survival
    surv2 <- survmat
    surv2$seedid <- seedling$seedid
    surv2_long <- melt (surv2, id.vars= "seedid")
    names(surv2_long) <- c("seedid", "census", "surv")

#height
    heightmat3 <- heightmat 
    heightmat3$seedid <- seedling$seedid
    heightmat3_long <- melt (heightmat3, id.vars= "seedid")  #heightmat3_long$height
    names(heightmat3_long) <- c("seedid", "census", "height")

#difference in height
    difheightmatb3 <- difheightb
    difheightmatb3$seedid <- seedling$seedid
    difheightmatb3_long <- melt (difheightmatb3, id.vars= "seedid")
    names(difheightmatb3_long) <- c("seedid", "census", "difheight")

#diameter  
    diammat3 <- diammat  #
    diammat3$seedid <- seedling$seedid
    diammat3_long <- melt (diammat3, id.vars= "seedid")  #diammat3_long$diam
    names(diammat3_long) <- c("seedid", "census", "diam")

#difference in diam
    difdiammatb3 <- difdiamb
    difdiammatb3$seedid <- seedling$seedid
    difdiammatb3_long <- melt (difdiammatb3, id.vars= "seedid")
    names(difdiammatb3_long) <- c("seedid", "census", "difdiam")

#herbivory
    hmat9 <- hmat4  #hmat 5 has herbivory smeared, hmat7 does not
    hmat9$seedid <- seedling$seedid
    hmat9_long <- melt(hmat9,id.vars="seedid")
    names(hmat9_long) <- c("seedid", "census", "h")

#pathogens
    pmat9 <- pmat4
    pmat9$seedid <- seedling$seedid
    pmat9_long <- melt(pmat9,id.vars="seedid")
    names(pmat9_long) <- c("seedid", "census", "p")

#age
    agemat10 <- as.data.frame(agematyear1vs2)  
    agemat10$seedid <- seedling$seedid
    agemat10_long <- melt(agemat10, id.vars = "seedid")
    names(agemat10_long) <- c("seedid", "census", "age")

#agedaymat
    agedaymat4 <- as.data.frame(agedaymat2)
    agedaymat4$seedid <- seedling$seedid
    agedaymat4_long <- melt(agedaymat4, id.vars = "seedid")
    names(agedaymat4_long) <- c("seedid", "census", "age_days")

#census season
    smat3 <- smat
    smat3$seedid <- seedling$seedid
    smat3_long <- melt(smat3, id.vars = "seedid")
    names(smat3_long) <- c("seedid", "census", "season")

#datemat
    datemat2 <- datemat
    datemat2$seedid <- seedling$seedid
    datemat2_long <- melt(datemat2, id.vars = "seedid")
    datemat2_long$value <- mdy(as.character(datemat2_long$value))
    names(datemat2_long) <- c("seedid", "census", "date")

#julian
    datemat_julian_long <- datemat2_long
    datemat_julian_long$date <- yday(datemat_julian_long$date)
    names(datemat_julian_long) <- c("seedid", "census", "julianday")

#nleaves
    nleaves4 <- nleaves
    nleaves4$seedid <- seedling$seedid
    nleaves4_long <- melt(nleaves4, id.vars = "seedid")
    names(nleaves4_long) <- c("seedid", "census", "nleaves")

#tmt
    tmtmat3 <- tmtmat
    tmtmat3$seedid <- seedling$seedid
    tmtmat3_long <- melt(tmtmat3, id.vars = "seedid")
    names(tmtmat3_long) <- c("seedid", "census", "tmt")

#browse
    browsemat3 <- browsemat
    browsemat3$seedid <- seedling$seedid
    browsemat3_long <- melt(browsemat3, id.vars = "seedid")
    names(browsemat3_long) <- c("seedid", "census", "browse")

#soil moisture #mean
    soilm_seedling$seedid <- seedling$seedid
    soilm_seedling2 <- subset(soilm_seedling, select =-c(plot, plot.letter, plotletter, plotsection,region,site, X, col, row))
    soilm_seedling2_long <- melt(soilm_seedling2 , id.vars = "seedid")
    names(soilm_seedling2_long) <- c("seedid","census","sm_mean")

#soil moisture #sd
    soilm_sd_seedling$seedid <- seedling$seedid
    soilm_sd_seedling2 <- subset(soilm_sd_seedling, select =-c(plot, plot.letter, plotletter, plotsection,region,site, X, col, row))
    soilm_sd_seedling2_long <- melt(soilm_sd_seedling2 , id.vars = "seedid")
    names(soilm_sd_seedling2_long) <- c("seedid","census","sm_sd")

#gsf
    #NOT INCLUDING DIFFERENCES IN TIME YET

#year
    yearmat3 <- as.data.frame(yearmat)
    yearmat3$seedid <- seedling$seedid
    yearmat3_long <- melt(yearmat3, id.vars = "seedid")
    names(yearmat3_long) <- c("seedid", "census", "year")

#######
#assembling them into one dataframe
#######

#single column variables from seedling database
    singlevars <- data.frame(seedling$seedid, seedling$planting_season, seedling$prepheight, seedling$region, seedling$site, 
                             seedling$region_label, seedling$site_label, seedling$plot, seedling$gsf_mean, seedling$gsf_sd,
                        seedling$Ntotal, seedling$NO3_N, seedling$NH4_N, seedling$P, seedling$Ca, seedling$Mg, seedling$Mn,
                        seedling$sp, seedling$migstat,
                        seedling$x, seedling$y, seedling$plot.letter,
                        seedling$ni_ntrees10, seedling$ni_allsp, seedling$ni_quru, seedling$ni_quve, seedling$ni_qual, seedling$ni_qusp,
                        seedling$ni_casp, seedling$ni_litu, seedling$ni_rops, seedling$ni_acru, seedling$ni_acsp, seedling$ni_ntrees,
                        seedling$ba_ntrees10, seedling$ba_allsp, seedling$ba_quru, seedling$ba_quve, seedling$ba_qual, seedling$ba_qusp, 
                        seedling$ba_casp, seedling$ba_litu, seedling$ba_rops, seedling$ba_acru, seedling$ba_acsp, 
                        seedling$d_allsp, seedling$d_quru, seedling$d_quve, seedling$d_qual, seedling$d_qusp, seedling$d_casp, 
                        seedling$d_litu, seedling$d_rops, seedling$d_acru, seedling$d_acsp,
                        seedling$conspecific_ba, seedling$conspecific_ni, seedling$conspecific_d, seedling$congeneric_ba,
                        seedling$conspecific_relba_site, seedling$conspecific_at_site,
                        seedling$multiscaled
                        )
    names(singlevars) <- gsub("seedling.","",names(singlevars))

survexp <-cbind(surv2_long, hmat9_long, pmat9_long, agemat10_long, 
                tmtmat3_long, browsemat3_long, soilm_seedling2_long, soilm_sd_seedling2_long, 
                heightmat3_long, difheightmatb3_long, diammat3_long, difdiammatb3_long,
                yearmat3_long,datemat2_long, datemat_julian_long, nleaves4_long, smat3_long, agedaymat4_long)
survexp <- survexp[, c(1,2,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48, 51, 54)]
survexp <- join(survexp, singlevars, by = "seedid", type = "left")

setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/phylogenetic distance")
write.csv(survexp, "diss_long151123.csv")
write_csv(survexp, "diss_long151123.csv")



