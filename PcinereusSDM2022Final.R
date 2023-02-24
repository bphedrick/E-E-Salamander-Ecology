#Plethodon cinereus SDM Paper
#Ecology and Evolution, Hedrick et al.
#Revised Jan 17, 2023 (Final)

#Download presence data from https://doi.org/10.15468/dl.d3j5ys

####CONTENTS####
#Lines 68-138
    #Load in data
    #Make plots of data bins (Fig S1)
#Lines 141 - 207
    #Define study area and split data into new (2001-2020) and old (1961-1980) presence data
    #Make plot showing the buffered area (Fig S2)
#Lines 210 - 401
    #Using all the variables from WorldClim included in the analysis, run GLM, GAM, and Maxent analyses
    #Generate comparison plots (Fig S3) between probabilities for GLM and GAM (Maxent compared with favorability later)
    #Calculate correlations between the three models (Table S1)
#Lines 404 - 472
    #Using the full set of variables (rather than the reduced model from multGLM), calculate favorability for GLM and GAM to compare with Maxent
    #Uses spc_GLM_P rather than Plci_P (which is from the multGLM model)
    #Make Fig S4 showing these comparisons
    #Calculate comparisons between these three models (GLM_F, GAM_F, Maxent) (Table S1)
#Lines 475 - 537
    #Build the multGLM trimming correlated/irrelevant variables (Table 2)
    #add the multGLM predictions to the data table (both the probabilities (Plci_P) and favorability (Plci_F))
#Lines 540 - 706
    #Calculate model evaluation metrics using different thresholds with the probabilities from the multGLM model (Plci_P)
    #Evaluation metrics figure (Fig S5)
#Lines 709 - 797
    #Compare new and old data (Fig 1)
    #Generate fuzzy range change metrics between 1961-1980 and 2001-2020
    #Plot expansion and contraction (Fig 2)
    #Calculate overlap between new and old data (Table S3)
#Lines 800 - 945
    #Download future layers from worldclim (CCSM4 and MIROC from 2070)
    #Crop those layers to only include the buffered range defined in our study
#Lines 948 - 1040
    #Predict probabilities of presence on future layers using the multGLM present model for CCSM4 and MIROC
    #Plot to visually compare probabilities
#Lines 1043 - 1225
    #Get favorability data for each future map using the multGLM model to predict
    #Save all these data into a single spatial dataframe (futureDatSpatial)
    #Plot CCSM4 under different RCPs (Fig S6) and MIROC (Fig S7)
    #Plot a subset of these RCP26 and RCP85 for both models for Fig 3
#Lines 1228 - 1371
    #Do fuzzy overlap operations on the maps (expansion, contraction) and plot (Fig S8, Table S3)
    #Calculate model overlap metrics (Table S3)
    #Plot fuzzy range change metrics (Fig 4)
#Lines 1374 - end
    #Run spatial autocorrelation


setwd("~/Desktop/Current Projects/RBSDistributionModeling/RBSDistributionModeling-RCode")

library(modEvA)
library(fuzzySim)
library(sp)
library(tidyr)
library(raster)
library(sdmpredictors)
library(gam)
library(rpart)
library(maxnet)
library(scrubr)
library(rworldmap)
library(sf)
library(rgdal)
library(spdep)

####Load in data####

salCoords <- read.csv("GBIFPCinereusSept22021.csv") #Load in data from gbif downloaded on Sept 2, 2021
      species <- c("Plethodon cinereus") #Define species name

      salCoords <- salCoords[!duplicated(salCoords$catalogNumber), ] #Remove specimens with identical catalog numbers

      #Make a simplified dataframe with only relevant columns
      salCoords <- data.frame(species = salCoords$species, 
                        decimalLongitude = salCoords$decimalLongitude, 
                        decimalLatitude = salCoords$decimalLatitude, 
                        occurrenceStatus = salCoords$occurrenceStatus, 
                        coordinateUncertaintyInMeters = salCoords$coordinateUncertaintyInMeters,
                        year = salCoords$year,
                        catalogNumber = salCoords$catalogNumber)

#load in map of countries with state boundaries      
statesMap <- readOGR("~/Desktop/Current Projects/RBSDistributionModeling/RBSDistributionModeling-RCode/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.shp")
    plot(statesMap, 
         xlim = range(salCoords$decimalLongitude), 
         ylim = range(salCoords$decimalLatitude))
    points(salCoords[ , c("decimalLongitude", "decimalLatitude")], col = "blue")

#Check for errors
nrow(salCoords)
      sort(unique(salCoords$occurrenceStatus))  # check for different indications of "absent", which could be in different languages
      absence_rows <- which(salCoords$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
            if (length(absence_rows) > 0)  salCoords <- salCoords[-absence_rows, ]
      salCoords <- dframe(salCoords) %>% coord_imprecise(which = "has_dec") #has decimals
      salCoords <- dframe(salCoords) %>% coord_imprecise(which = "no_zeros") #all zeros
      salCoords <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(salCoords))))
      salCoords <- coord_uncertain(salCoords, coorduncertainityLimit = 7000) #remove highly uncertain data
nrow(salCoords)

points(salCoords[ , c("decimalLongitude", "decimalLatitude")], col = "turquoise")

#Remove additional outliers from data that are outside of northeastern North America
salCoords <- salCoords[-(which(salCoords$decimalLongitude <= -94)), ]
salCoords <- salCoords[-(which(salCoords$decimalLatitude <= 32)), ]

#Check final map
plot(statesMap, xlim = range(salCoords$decimalLongitude), ylim = range(salCoords$decimalLatitude))
    points(salCoords[ , c("decimalLongitude", "decimalLatitude")], col = "green")


####Make plots of data bins (Fig S1)####
hist(salCoords$year) #Subset data to look at years

salCoords <- subset(salCoords, year >= 1800) #This removes data from pre-1800 and year == NA
      min(salCoords$year) #oldest record

nrow(subset(salCoords, year <= 1980)) / nrow(salCoords) #Ratio of before 1980 to total
nrow(subset(salCoords, year >= 1980 & year <= 2010)) / nrow(salCoords) #Ratio between 1980 and 2010

#Subset data into two categories
salCoordsOld <- subset(salCoords, year >= 1961 & year <= 1980)
    nrow(salCoordsOld)

salCoordsNew <- subset(salCoords, year >= 2001 & year <= 2020)
    nrow(salCoordsNew)

#Figure S1, make a comparison plot across the year bins, histogram, 1961-1980, and 2001-2020
par(mfrow = c(3, 1), mar = c(2, 2, 2, 1))

hist(salCoords$year, main = "")
plot(statesMap, xlim = range(salCoords$decimalLongitude), ylim = range(salCoords$decimalLatitude),
     main = "1961 - 1980")
  points(salCoordsOld[ , c("decimalLongitude", "decimalLatitude")], col = "red")
plot(statesMap, xlim = range(salCoords$decimalLongitude), ylim = range(salCoords$decimalLatitude),
     main = "2001 - 2020")
  points(salCoordsNew[ , c("decimalLongitude", "decimalLatitude")], col = "blue")  
  
  
####Define Study Area####
#Used same buffer for old (1961-1980) and new (2001-2020) datasets based on the all presences buffer

salPoints <- salCoords #Change to a spatial points dataframe
    coordinates(salPoints) <- salPoints[ , c("decimalLongitude", "decimalLatitude")] 
    crs(salPoints) <- "+proj=longlat" #Set a coordinate reference system
    plot(salPoints, col = "red")
          plot(statesMap, lwd = 1, add = TRUE)
  
#Use buffer to get study area
    presBuff <- raster::buffer(salPoints, width = 350 * 1000, dissolve = TRUE) #Use a buffer around the datapoints
          plot(presBuff, lwd = 2)
          plot(salPoints, col = "blue", add = TRUE)
          plot(statesMap, border = "black", add = TRUE)
    studyArea <- presBuff 

#now import and cut (crop + mask) the variable maps to the extent of the study area defined above

#Import the worldclim data, exclude bio3, bio14, bio15 (low correlation between current and future variables)    
layerNames <- c("WC_alt", "WC_bio1", "WC_bio2", "WC_bio4", "WC_bio5", "WC_bio6", "WC_bio7",
                "WC_bio8", "WC_bio9", "WC_bio10", "WC_bio11", "WC_bio12", "WC_bio13", "WC_bio16",
                "WC_bio17", "WC_bio18", "WC_bio19")

layers <- load_layers(layerNames)
    plot(layers) #Check that layers loaded in properly
layersCut <- mask(crop(layers, studyArea), studyArea) #Crop them to the defined study area
    plot(layersCut) #Replot as just the study area

#Figure S2, example buffer figure using altitude map, plotted at 1200 width rather than as a pdf
par(mfrow = c(2, 1), mar = c(2, 2, 2, 1))

plot(layersCut[[1]])
    plot(statesMap, border = "grey", add = TRUE)
    plot(studyArea, add = TRUE)

plot(layersCut[[1]])
    plot(statesMap, border = "grey", add = TRUE)
    plot(studyArea, add = TRUE)
    plot(salPoints, add = TRUE, pch = 1)

#set the appropriate spatial resolution
    
#Set species names as abbreviations (Plci)
salCoordsNew$speciesCode <- spCodes(species = salCoordsNew$species, nchar.gen = 2, nchar.sp = 2)
salCoordsOld$speciesCode <- spCodes(species = salCoordsOld$species, nchar.gen = 2, nchar.sp = 2)
    salSpeciesCodes <- spCodes(species = species, nchar.gen = 2, nchar.sp = 2)
          salSpeciesCodes 
          
#this function (soon to be included in the fuzzySim package) is a wrapper for fuzzySim::gridRecords to 
#grid more than one species at a time.
source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/gridRecords_mult")  

#For new (2001-2020)
datNew <- gridRecords_mult(rst = layersCut, 
                        sp.data = as.data.frame(salCoordsNew), 
                        sp.col = "speciesCode", 
                        coord.cols = c("decimalLongitude", "decimalLatitude"))

#For old (1961-1980)
datOld <- gridRecords_mult(rst = layersCut, 
                        sp.data = as.data.frame(salCoordsOld), 
                        sp.col = "speciesCode", 
                        coord.cols = c("decimalLongitude", "decimalLatitude"))

    nrow(datOld)  # should be the same number as:
    nrow(datNew)  # should be the same number as:
    sum(!is.na(getValues(layersCut[[1]])))

    
####Build Distribution Models (GLM, GAM, Maxent) and Compare Them####

names(datNew)
spcCols <- 2  # species codes are in these columns IN THIS DATASET (change as necessary)
varCols <- 5:21  # variables are in these columns IN THIS DATASET (change as necessary)
names(datNew)[spcCols]  # check OK
names(datNew)[varCols]  # check OK

names(datOld)[spcCols]  # check OK
names(datOld)[varCols]  # check OK

spc <- salSpeciesCodes[1] 

####GLM for both the new (2001-2020) and old (1961-1980) data
#Make a formula that includes all of the variance columns with the given species as the independent variable
formGLMNew <- as.formula(paste(spc, "~", paste(names(datNew)[varCols], collapse = " + ")))
    formGLMNew
    
    # GLM with all specified variables, this is binomial because it is based on presence/absence data
    spcGLMNew <- glm(formula = formGLMNew, family = binomial, data = datNew)
        summary(spcGLMNew)  

#Do this again for old data (1961-1980)
formGLMOld <- as.formula(paste(spc, "~", paste(names(datOld)[varCols], collapse = " + ")))
    formGLMOld
    
    spcGLMOld <- glm(formula = formGLMOld, family = binomial, data = datOld)
        summary(spcGLMOld)  

#apply the GLM predictions to the data table (variables must have exact same names as in the model):
datNew$spc_GLM_P <- predict(spcGLMNew, newdata = datNew, type = "response")  #"response" applies the link function, i.e. yields results as probability values
    head(datNew)

datOld$spc_GLM_P <- predict(spcGLMOld, newdata = datOld, type = "response") 
    head(datOld)
    
#map GLM predictions for this species:
datSpatialNew <- datNew #convert 'dat' to a spatial object
        coordinates(datSpatialNew) <- datSpatialNew[ , c("x", "y")]
    spplot(datSpatialNew, zcol = "spc_GLM_P", cex = 0.5) 

datSpatialOld <- datOld
        coordinates(datSpatialOld) <- datSpatialOld[ , c("x", "y")]
    spplot(datSpatialOld, zcol = "spc_GLM_P", cex = 0.5) 

####GAM for both the new (2001-2020) and old (1961-1980) data
formGAMNew <- as.formula(paste(spc, "~", paste0("s(", names(datNew)[varCols], ")", collapse = "+")))  # GAM with smoothing splines ('s')
        formGAMNew
    spcGAMNew <- gam(formula = formGAMNew, family = binomial, data = datNew)
        summary(spcGAMNew)

formGAMOld <- as.formula(paste(spc, "~", paste0("s(", names(datOld)[varCols], ")", collapse = "+")))  # GAM with smoothing splines ('s')
        formGAMOld
    spcGAMOld <- gam(formula = formGAMOld, family = binomial, data = datOld)
        summary(spcGAMOld)

#apply predictions to the data table and map them
datNew$spc_GAM_P <- predict(spcGAMNew, datNew, type = "response")
          head(datNew)
    datSpatialNew$spc_GAM_P <- datNew$spc_GAM_P
          spplot(datSpatialNew, zcol = "spc_GAM_P", cex = 0.5)

datOld$spc_GAM_P <- predict(spcGAMOld, datOld, type = "response")
          head(datOld)
    datSpatialOld$spc_GAM_P <- datOld$spc_GAM_P
          spplot(datSpatialOld, zcol = "spc_GAM_P", cex = 0.5)
    
####Maxent (presence-background based model) for both the new (2001-2020) and old (1961-1980) data
          
#argument "p" should be "a vector of 1 (for presence) or 0 (for background)", not 0 for absence as we have in 'dat'
#so let's add to the table a replicate of the presence rows and turn them to zeros, to include them in the background:
nrow(datNew)
nPresNew <- sum(datNew[ , spc])
    nPresNew
spcPres0New <- subset(datNew, get(spc) == 1)
    spcPres0New[ , spc] <- 0 # convert the replicated presences to "0" for background:
    spcPbNew <- rbind(datNew, spcPres0New) # join the original table to the newly created background:

    #check that everything went as expected:
    nrow(spcPbNew) == nrow(datNew) + nPresNew  # should be TRUE
    sum(datNew[ , spc]) == sum(spcPbNew[ , spc])  # should also be TRUE

#Do this again for the old data (1961-1980)
nrow(datOld)
nPresOld <- sum(datOld[ , spc])
    nPresOld
spcPres0Old <- subset(datOld, get(spc) == 1)
    spcPres0Old[ , spc] <- 0
    spcPbOld <- rbind(datOld, spcPres0Old)
    nrow(spcPbOld) == nrow(datOld) + nPresOld  # should be TRUE
    sum(datOld[ , spc]) == sum(spcPbOld[ , spc])  # should also be TRUE

#Run a Maxent model with the presence-background data for 'spc' for the new data (2001-2020)
spcMXNew <- maxnet(p = spcPbNew[ , spc], 
                 data = spcPbNew[ , varCols], 
                 f = maxnet.formula(spcPbNew[ , spc], 
                                    spcPbNew[ , varCols], 
                                    classes = "lq"))  # using linear ('l') and quadratic ('q') features
spcMXNew

#get Maxent predictions on the 'dat' table and map
datNew$spc_MX_P <- as.vector(predict(spcMXNew, newdata = datNew, type = "cloglog"))
        head(datNew)
    datSpatialNew$spc_MX_P <- datNew$spc_MX_P
        spplot(datSpatialNew, zcol = "spc_MX_P", cex = 0.5)

#Run a Maxent model with the presence-background data for 'spc' for the old data (1961-1980)
spcMXOld <- maxnet(p = spcPbOld[ , spc], 
                 data = spcPbOld[ , varCols], 
                 f = maxnet.formula(spcPbOld[ , spc], 
                                    spcPbOld[ , varCols], 
                                    classes = "lq"))  # using linear ('l') and quadratic ('q') features
spcMXOld

datOld$spc_MX_P <- as.vector(predict(spcMXOld, newdata = datOld, type = "cloglog"))
        head(datOld)
    datSpatialOld$spc_MX_P <- datOld$spc_MX_P
        spplot(datSpatialOld, zcol = "spc_MX_P", cex = 0.5)
    
        
##Compare the predictions from the three modeling methods (GLM, GAM, Maxent)
predCols <- c("spc_GLM_P", "spc_GAM_P", "spc_MX_P")

#map all model predictions together with the same color scale:
spplot(datSpatialNew, zcol = predCols) #New data (2001-2020)
spplot(datSpatialOld, zcol = predCols) #Old data (1961-1980)

#predict to the raster maps of the variables for the new data (2001-2020)
GLM_P_rastNew <- predict(layersCut, spcGLMNew, type = "response")
GAM_P_rastNew <- predict(layersCut, spcGAMNew, type = "response")
MX_P_rastNew <- predict(layersCut, spcMXNew, type = "cloglog", clamp = FALSE)

#plot the raster predictions (note different scales)
par(mfrow = c(3, 1), mar = c(3, 3, 2, 2))
    plot(GLM_P_rastNew, main = "GLM (2001-2020)")
    plot(GAM_P_rastNew, main = "GAM (2001-2020)")
    plot(MX_P_rastNew, main = "Maxent (2001-2020)")

#predict to the raster maps of the variables for the old data (1961-1980)
GLM_P_rastOld <- predict(layersCut, spcGLMOld, type = "response")
GAM_P_rastOld <- predict(layersCut, spcGAMOld, type = "response")
MX_P_rastOld <- predict(layersCut, spcMXOld, type = "cloglog", clamp = FALSE)

par(mfrow = c(3, 1), mar = c(3, 3, 2, 2)) #Different scales
    plot(GLM_P_rastOld, main = "GLM (1961-1980)")
    plot(GAM_P_rastOld, main = "GAM (1961-1980)")
    plot(MX_P_rastOld, main = "Maxent (1961-1980)")

#Plot for Figure S3 to show the GLM and GAM comparisons for both time slices. Maxent is compared later with favorability maps
par(mfrow = c(2, 2), mar = c(3, 3, 2, 2))

plot(GLM_P_rastOld,
     main = "GLM (1961-1980)", 
     breaks = seq(0, 0.75, by = 0.001), 
     col = terrain.colors(700, rev = TRUE),
     legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
plot(GLM_P_rastNew,
     main = "GLM (2001-2020)", 
     breaks = seq(0, 0.75, by = 0.001), 
     col = terrain.colors(700, rev = TRUE),
     legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
plot(GAM_P_rastOld,
     main = "GAM (1961-1980)", 
     breaks = seq(0, 0.75, by = 0.001), 
     col = terrain.colors(700, rev = TRUE),
     legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
plot(GAM_P_rastNew,
     main = "GAM (2001-2020)", 
     breaks = seq(0, 0.75, by = 0.001), 
     col = terrain.colors(700, rev = TRUE),
     legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
  
#check correlations among model predictions, Table S1
pairs(datNew[ , predCols], pch = 20, cex = 0.5) #Plot correlations between three models
    predCorrsNew <- cor(datNew[ , predCols]) #Check correlations for all three models
    predCorrsNew 
    min(predCorrsNew, na.rm = TRUE) #What is the smallest correlation?

pairs(datOld[ , predCols], pch = 20, cex = 0.5) 
    predCorrsOld <- cor(datOld[ , predCols]) 
    predCorrsOld 
    min(predCorrsOld, na.rm = TRUE)

#given that correlations among predictions are high, choosing a particular modeling method does not make a difference


####Probability and Favorability####
#Convert presence probability values to favorability values. Can't use maxent since it is suitability rather than favorability

#Generate favorability data from the GLM model and the GAM model. 
#Note that this does not use the decreased model from the multGLM function
#This is merely to compare favorability from the GLM model to that of the GAM model and to the Maxent model

#First for the new data (2001-2020)
datNew$spc_GLM_F <- Fav(spcGLMNew)
datNew$spc_GAM_F <- Fav(spcGAMNew)

head(datNew)

# plot sorted predictions for each model:
par(mfrow = c(2, 1), mar = c(2, 2, 1, 1)) 
plot(sort(datNew$spc_GLM_P))
plot(sort(datNew$spc_GLM_F))
plot(sort(datNew$spc_GAM_P))
plot(sort(datNew$spc_GAM_F))

#Do the same thing with the old data (1961-1980)
datOld$spc_GLM_F <- Fav(spcGLMOld)
datOld$spc_GAM_F <- Fav(spcGAMOld)

head(datOld)

#plot sorted predictions for each model:
par(mfrow = c(2, 1), mar = c(2, 2, 1, 1)) 
plot(sort(datOld$spc_GLM_P))
plot(sort(datOld$spc_GLM_F))
plot(sort(datOld$spc_GAM_P))
plot(sort(datOld$spc_GAM_F))

#map the predictions comparing the data from the probabilities to the favorability
datSpatialNew <- datNew
names(datSpatialNew)
coordinates(datSpatialNew) <- datSpatialNew[ , c("x", "y")]
spplot(datSpatialNew, zcol = c("spc_GLM_P", "spc_GAM_P", "spc_GLM_F", "spc_GAM_F"), cex = 0.3)

datSpatialOld<- datOld
names(datSpatialOld)
coordinates(datSpatialOld) <- datSpatialOld[ , c("x", "y")]
spplot(datSpatialOld, zcol = c("spc_GLM_P", "spc_GAM_P", "spc_GLM_F", "spc_GAM_F"), cex = 0.3)

#Make Figure S4 comparing favorability for GLM and GAM to Maxent
spplot(datSpatialNew, 
       zcol = c("spc_GLM_F", "spc_GAM_F", "spc_MX_P"), 
       cex = 0.3,
       col.regions =  terrain.colors(1000, rev = TRUE))

spplot(datSpatialOld, 
       zcol = c("spc_GLM_F", "spc_GAM_F", "spc_MX_P"), 
       cex = 0.3,
       col.regions =  terrain.colors(1000, rev = TRUE))

#check correlations among model predictions for new and old data, Table S1
predColsFav <- c("spc_GLM_F", "spc_GAM_F", "spc_MX_P")

pairs(datNew[ , predColsFav], pch = 20, cex = 0.5) #Plot correlations between three models
predCorrsNewFav <- cor(datNew[ , predColsFav]) #Check correlations for all three models
predCorrsNewFav 
min(predCorrsNewFav, na.rm = TRUE) #What is the smallest correlation?

pairs(datOld[ , predColsFav], pch = 20, cex = 0.5) #Plot correlations between three models
predCorrsOldFav <- cor(datOld[ , predColsFav]) #Check correlations for all three models
predCorrsOldFav 
min(predCorrsOldFav, na.rm = TRUE) #What is the smallest correlation?

#This demonstrates that model correlations are very high across the favorability models and maxent


####Build the multGLM model for the old and the new data####
set.seed(123)  #set a particular seed so the next command consistently chooses the same test data sample and produces the same results

modsMultGLMNew <- multGLM(data = datNew[complete.cases(datNew[ , varCols]), ],
                        sp.cols = salSpeciesCodes,  # or choose your species, e.g. sp.cols = c(27, 41, 50)
                        var.cols = varCols,
                        id.col = 1,
                        test.sample = 0.2,  # reserve 20% data for model testing (optional)
                        FDR = TRUE,  # remove variables following false discovery rate
                        corSelect = TRUE,  # filter correlated variables
                        cor.thresh = 0.8,  # correlation threshold
                        step = TRUE,  # do stepwise AIC variable selection
                        start = "null.model",  # start with no variables
                        direction = "both",  # forward-backward stepwise selection
                        trim = TRUE)  # remove non-significant variables

set.seed(123)

modsMultGLMOld <- multGLM(data = datOld[complete.cases(datOld[ , varCols]), ],
                        sp.cols = salSpeciesCodes,  # or choose your species, e.g. sp.cols = c(27, 41, 50)
                        var.cols = varCols,
                        id.col = 1,
                        test.sample = 0.2,  # reserve 20% data for model testing (optional)
                        FDR = TRUE,  # don't remove variables following false discovery rate
                        corSelect = TRUE,  # filter correlated variables
                        cor.thresh = 0.8,  # correlation threshold
                        step = TRUE,  # do stepwise AIC variable selection
                        start = "null.model",  # start with no variables
                        direction = "both",  # forward-backward stepwise selection
                        trim = TRUE)

# look at the multGLM results (Table 2)
names(modsMultGLMNew)  # 3 resulting objects
modsMultGLMNew$variables #Different variables used for different species based on the stepwise AICs
names(modsMultGLMNew$predictions)
head(modsMultGLMNew$predictions) #Created a training and test dataset, as well as probability of presence and favorability using all three models. 20% will be test and 80% will be train, test should be between 15-30%

lapply(modsMultGLMNew$models, summary) #Model summary for new data (2001-2020)

names(modsMultGLMOld)  
modsMultGLMOld$variables #Different variables used for different species based on the stepwise AICs
names(modsMultGLMOld$predictions)
head(modsMultGLMOld$predictions) #Created a training and test dataset, as well as probability of presence and favorability using all three models. 20% will be test and 80% will be train, test should be between 15-30%

lapply(modsMultGLMOld$models, summary) #Model summary for old data (1961-1980)

# add the multGLM predictions to the data table (both the predictions (Plci_P) and favorability (Plci_F))
datNew <- merge(datNew, modsMultGLMNew$predictions, by = "cells", sort = FALSE, all.x = TRUE)
head(datNew)

datOld <- merge(datOld, modsMultGLMOld$predictions, by = "cells", sort = FALSE, all.x = TRUE)
head(datOld)

#map the multGLM probability ("_P") predictions and the species occurrences:
datSpatialNew <- datNew
    names(datSpatialNew)
        coordinates(datSpatialNew) <- datSpatialNew[ , c("x", "y")]
    spplot(datSpatialNew[ , c("Plci_P")], cex = 0.3)

datSpatialOld <- datOld
    names(datSpatialOld)
        coordinates(datSpatialOld) <- datSpatialOld[ , c("x", "y")]
    spplot(datSpatialOld[ , c("Plci_P")], cex = 0.3)


####Model Evaluation####
#Data go into Table S3
    
#separate the model training and testing data for evaluation:
head(datNew)
    datTrainNew <- subset(datNew, sample == "train")
    datTestNew <- subset(datNew, sample == "test")
        nrow(datTrainNew)
        nrow(datTestNew) #In this case, test is 20% of train

head(datOld)
    datTrainOld <- subset(datOld, sample == "train")
    datTestOld <- subset(datOld, sample == "test")
        nrow(datTrainOld)
        nrow(datTestOld) #In this case, test is 20% of train

        
#Run model evaluation for the new data (2001-2020)        
#compute some threshold-based classification measures of model performance using prevalence as the threshold
with(datTestNew, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = "preval", 
                             ylim = c(0, 1),
                             standardize = FALSE, #Get regular kappa and TSS rather than skappa and sTSS when = FALSE
                             main = "Prevalence (2001-2020)", las = 0))  # here we are using species prevalence as the classification threshold

#Try optimizing the threshold based on a particular metric or criterion:
with(datTestNew, optiThresh(obs = Plci, 
                         pred = Plci_P))
optimizeTSSNew <- with(datTestNew, 
                     optiThresh(obs = Plci, 
                                pred = Plci_P, 
                                measures = "TSS", 
                                optimize = "each", 
                                interval = 0.001))
    optimizeTSSNew

threshMaxTSSNew <- optimizeTSSNew$optimals.each$threshold
    threshMaxTSSNew

#Repeat the 'threshMeasures' command using, as the 'thresh' argument, thresh_maxTSS instead of "preval"
with(datTestNew, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = threshMaxTSSNew, 
                             ylim = c(0, 1),
                             standardize = FALSE,
                             main = "threshMaxTSS (2001-2020)", las = 0)) 

# optimizing TSS is the same as optimizing sensitivity and specificity:
with(datTestNew, optiPair(obs = Plci, 
                        pred = Plci_P, 
                        measures = c("Sensitivity", "Specificity"), 
                        interval = 0.001))

# compute overall discrimination capacity with the area under the ROC curve:
with(datTestNew, AUC(obs = Plci, pred = Plci_P))

# compute explained deviance and pseudo-R-squared metrics:
with(datTestNew, Dsquared(obs = Plci, pred = Plci_P, family = "binomial"))
with(datTestNew, RsqGLM(obs = Plci, pred = Plci_P))

# compute correlation between observations and predictions:
with(datTestNew, cor.test(Plci, Plci_P)) #Check that these are significantly correlated.

# compute the Miller calibration line and parameters:
par(mfrow = c(3, 1))

with(datNew, MillerCalib(obs = Plci, pred = Plci_P))
with(datTrainNew, MillerCalib(obs = Plci, pred = Plci_P))  # with GLMs, calibration is always perfect on the training data
with(datTestNew, MillerCalib(obs = Plci, pred = Plci_P)) #This is the important one. Well-calibrated models have slope close to 1


#Run model evaluation for the old data (1961-1980)    
with(datTestOld, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = "preval", 
                             ylim = c(0, 1),
                             standardize = FALSE,
                             main = "Prevalence (1961-1980)", las = 0))  # here we are using species prevalence as the classification threshold

with(datTestOld, optiThresh(obs = Plci, 
                         pred = Plci_P))
optimizeTSSOld <- with(datTestOld, 
                     optiThresh(obs = Plci, 
                                pred = Plci_P, 
                                measures = "TSS", 
                                optimize = "each", 
                                interval = 0.001))
    optimizeTSSOld

threshMaxTSSOld <- optimizeTSSOld$optimals.each$threshold
    threshMaxTSSOld

# you could now repeat the 'threshMeasures' command using, as the 'thresh' argument, thresh_maxTSS instead of "preval"
with(datTestOld, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = threshMaxTSSOld, 
                             ylim = c(0, 1),
                             standardize = FALSE,
                             main = "threshMaxTSS (1961-1980)", las = 0)) 

# optimizing TSS is the same as optimizing sensitivity and specificity:
with(datTestOld, optiPair(obs = Plci, 
                        pred = Plci_P, 
                        measures = c("Sensitivity", "Specificity"), 
                        interval = 0.001))

# compute overall discrimination capacity with the area under the ROC curve:
with(datTestOld, AUC(obs = Plci, pred = Plci_P))

# compute explained deviance and pseudo-R-squared metrics:
with(datTestOld, Dsquared(obs = Plci, pred = Plci_P, family = "binomial"))
with(datTestOld, RsqGLM(obs = Plci, pred = Plci_P))

# compute correlation between observations and predictions:
with(datTestOld, cor.test(Plci, Plci_P)) #Check that these are significantly correlated. They are. Awesome.

# compute the Miller calibration line and parameters:
par(mfrow = c(3, 1))

with(datOld, MillerCalib(obs = Plci, pred = Plci_P))
with(datTrainOld, MillerCalib(obs = Plci, pred = Plci_P))  # with GLMs, calibration is always perfect on the training data
with(datTestOld, MillerCalib(obs = Plci, pred = Plci_P)) #This is the important one
# well-calibrated models have slope close to 1

#Figure S5, a plot showing all of this together
par(mfrow = c(3, 2))

with(datTestOld, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = "preval", 
                             standardize = FALSE,
                             ylim = c(0, 1),
                             main = "Prevalence (1961-1980)", las = 0))

with(datTestNew, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = "preval", 
                             standardize = FALSE,
                             ylim = c(0, 1),
                             main = "Prevalence (2001-2020)", las = 0))

with(datTestOld, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = threshMaxTSSOld, 
                             standardize = FALSE,
                             ylim = c(0, 1),
                             main = "threshMaxTSS (1961-1980)", las = 0)) 

with(datTestNew, threshMeasures(obs = Plci, 
                             pred = Plci_P, 
                             measures = c("CCR", "Sensitivity", "Specificity", "Precision", "kappa", "TSS"), 
                             thresh = threshMaxTSSNew, 
                             standardize = FALSE,
                             ylim = c(0, 1),
                             main = "threshMaxTSS (2001-2020)", las = 0)) 

with(datTestOld, AUC(obs = Plci, pred = Plci_P))

with(datTestNew, AUC(obs = Plci, pred = Plci_P))


####Compare between 1961-1980 and 2001-2020####
predNewGLM <- predict(layersCut, modsMultGLMNew$models[["Plci"]], type = "response")
predOldGLM <- predict(layersCut, modsMultGLMOld$models[["Plci"]], type = "response")

#Plot of GLMs based on presences from multGLM 
plot(predOldGLM, 
         main = "GLM 1961-1980", 
         breaks = seq(0, 0.5, by = 0.001), 
         col = terrain.colors(500, rev = TRUE),
         legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
plot(predNewGLM, 
         main = "GLM 2001-2020", 
         breaks = seq(0, 0.5, by = 0.001), 
         col = terrain.colors(500, rev = TRUE),
         legend = TRUE) 
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)

#Make one spatial dataset for old and new using favorability generated by the multGLM function with variables removed
datSpatialOldNew <- datSpatialNew[, 1:4]
    datSpatialOldNew$Old_F <- datOld$Plci_F 
    datSpatialOldNew$New_F <- datNew$Plci_F

#Make one dataset of non-spatial values using favorability generated by the multGLM function with variables removed
datOldNew <- datNew[, 1:4]
    datOldNew$Old_F <- datOld$Plci_F 
    datOldNew$New_F <- datNew$Plci_F    
    
##Figure 1##
#Plot of favorability from the multGLM function, use to make Figure 1
#New (2001-2020)
spplot(datSpatialOldNew, 
            zcol = c("New_F"), 
            cex = 0.3,
            col.regions = terrain.colors(1000, rev = TRUE))

#Old (1961-1980)
spplot(datSpatialOldNew, 
       zcol = c("Old_F"), 
       cex = 0.3,
       col.regions = terrain.colors(1000, rev = TRUE))

dev.off() #Print this as well and put on spplot maps in photoshop
  plot(presBuff, lwd = 2) 
  plot(statesMap, border = "grey", add = TRUE)


#Fuzzy range change metrics using the favorability data
datOldNew$expansion <- fuzzyOverlay(data = datOldNew, 
                                            overlay.cols = c("Old_F", "New_F"), 
                                            op = "expansion")
datOldNew$contraction <- fuzzyOverlay(data = datOldNew, 
                                            overlay.cols = c("Old_F", "New_F"), 
                                            op = "contraction")

#Make a spatial dataframe for looking at fuzzy range change metrics and then plot
datFuzzyMetricsOldNew <- datOldNew
    coordinates(datFuzzyMetricsOldNew) <- datFuzzyMetricsOldNew[ , c("x", "y")]
    names(datFuzzyMetricsOldNew)

##Figure 2##
spplot(datFuzzyMetricsOldNew, 
       c("expansion"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "Expansion from 1961-1980 and 2001-2020",
       cex = 0.3)
spplot(datFuzzyMetricsOldNew, 
       c("contraction"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "Contraction from 1961-1980 and 2001-2020",
       cex = 0.3)


dev.off() #Reset sizing

#Table S3
#Calculate fuzzy range change metrics between old (1961-1980) and new (2001-2020) models
fuzzyRangeChange(datOldNew$Old_F, datOldNew$New_F,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.25, 0.5),
                 main = "Range Change between 1961-1980 and 2001-2020",
                 las = 2)

#Check overlap between the old and new data
modOverlap(datOldNew$Old_F, datOldNew$New_F) #Schoener's D, Warren's I, Hellinger Distance
with(datOldNew, fuzSim(Old_F, New_F, method = "Jaccard"))
with(datOldNew, fuzSim(Old_F, New_F, method = "Baroni"))


####Download Future layers####

#Train the future climate model using the past data that were generated
#I am using the rcp model here for 2070 under 4 GHG emission rates, this eliminates bio3, 14, 15 (see methods)

futureLayersRCP26 <- get_future_layers(c("WC_bio1", "WC_bio2", "WC_bio4", "WC_bio5", "WC_bio6",
                                         "WC_bio7", "WC_bio8", "WC_bio9", "WC_bio10","WC_bio11",
                                         "WC_bio12", "WC_bio13", "WC_bio16", "WC_bio17", "WC_bio18",
                                         "WC_bio19"), scenario = "rcp26", year = 2070)

futureLayersRCP45 <- get_future_layers(c("WC_bio1", "WC_bio2", "WC_bio4", "WC_bio5", "WC_bio6",
                                         "WC_bio7", "WC_bio8", "WC_bio9", "WC_bio10", "WC_bio11",
                                         "WC_bio12", "WC_bio13", "WC_bio16", "WC_bio17", "WC_bio18",
                                         "WC_bio19"), scenario = "rcp45", year = 2070)

futureLayersRCP60 <- get_future_layers(c("WC_bio1", "WC_bio2", "WC_bio4", "WC_bio5", "WC_bio6",
                                         "WC_bio7", "WC_bio8", "WC_bio9", "WC_bio10", "WC_bio11",
                                         "WC_bio12", "WC_bio13", "WC_bio16", "WC_bio17", "WC_bio18",
                                         "WC_bio19"), scenario = "rcp60", year = 2070)

futureLayersRCP85 <- get_future_layers(c("WC_bio1", "WC_bio2", "WC_bio4", "WC_bio5", "WC_bio6",
                                         "WC_bio7", "WC_bio8", "WC_bio9", "WC_bio10", "WC_bio11",
                                         "WC_bio12", "WC_bio13", "WC_bio16", "WC_bio17", "WC_bio18",
                                         "WC_bio19"), scenario = "rcp85", year = 2070)


#CCSM4
futureLayersCCSM4rcp26 <- futureLayersRCP26[which(futureLayersRCP26$model == "CCSM4"), ]
futureLayersCCSM4rcp45 <- futureLayersRCP45[which(futureLayersRCP45$model == "CCSM4"), ]
futureLayersCCSM4rcp60 <- futureLayersRCP60[which(futureLayersRCP60$model == "CCSM4"), ]
futureLayersCCSM4rcp85 <- futureLayersRCP85[which(futureLayersRCP85$model == "CCSM4"), ]

#MIROC-ESM
futureLayersMIROCrcp26 <- futureLayersRCP26[which(futureLayersRCP26$model == "MIROC-ESM"), ]
futureLayersMIROCrcp45 <- futureLayersRCP45[which(futureLayersRCP45$model == "MIROC-ESM"), ]
futureLayersMIROCrcp60 <- futureLayersRCP60[which(futureLayersRCP60$model == "MIROC-ESM"), ]
futureLayersMIROCrcp85 <- futureLayersRCP85[which(futureLayersRCP85$model == "MIROC-ESM"), ]

#This will load the layers and overwrite the variables outlined above replacing them as raster stacks
futureLayersCCSM4rcp26 <- load_layers(futureLayersCCSM4rcp26)
futureLayersCCSM4rcp45 <- load_layers(futureLayersCCSM4rcp45)
futureLayersCCSM4rcp60 <- load_layers(futureLayersCCSM4rcp60)
futureLayersCCSM4rcp85 <- load_layers(futureLayersCCSM4rcp85)

futureLayersMIROCrcp26 <- load_layers(futureLayersMIROCrcp26)
futureLayersMIROCrcp45 <- load_layers(futureLayersMIROCrcp45)
futureLayersMIROCrcp60 <- load_layers(futureLayersMIROCrcp60)
futureLayersMIROCrcp85 <- load_layers(futureLayersMIROCrcp85)

##Crop CCSM4 and MIROC layers to fit the study areas
#Cut and fix cc26
futureLayersCCSM4rcp26Cut <- mask(crop(futureLayersCCSM4rcp26, studyArea), studyArea)
    plot(futureLayersCCSM4rcp26Cut)
        names(futureLayersCCSM4rcp26Cut) <- gsub("_cc26_2070", "", names(futureLayersCCSM4rcp26Cut)) #Change names to match other data
            names(futureLayersCCSM4rcp26Cut)

futureLayersCCSM4rcp26Cut <- stack(layersCut[["WC_alt"]], futureLayersCCSM4rcp26Cut) #Do alt separately
    names(futureLayersCCSM4rcp26Cut)

#Reorder the layers so that they match
futureLayersCCSM4rcp26Cut <- subset(futureLayersCCSM4rcp26Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))    
      names(layersCut) == names(futureLayersCCSM4rcp26Cut) #Check that they match the original layer order

#Cut and fix cc45
futureLayersCCSM4rcp45Cut <- mask(crop(futureLayersCCSM4rcp45, studyArea), studyArea)
    plot(futureLayersCCSM4rcp45Cut)
        names(futureLayersCCSM4rcp45Cut) <- gsub("_cc45_2070", "", names(futureLayersCCSM4rcp45Cut))
            names(futureLayersCCSM4rcp45Cut)

futureLayersCCSM4rcp45Cut <- stack(layersCut[["WC_alt"]], futureLayersCCSM4rcp45Cut)
    names(futureLayersCCSM4rcp45Cut)

futureLayersCCSM4rcp45Cut <- subset(futureLayersCCSM4rcp45Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))    
    names(layersCut) == names(futureLayersCCSM4rcp45Cut)

#Cut and fix cc60
futureLayersCCSM4rcp60Cut <- mask(crop(futureLayersCCSM4rcp60, studyArea), studyArea)
    plot(futureLayersCCSM4rcp60Cut)
        names(futureLayersCCSM4rcp60Cut) <- gsub("_cc60_2070", "", names(futureLayersCCSM4rcp60Cut))
            names(futureLayersCCSM4rcp60Cut)

futureLayersCCSM4rcp60Cut <- stack(layersCut[["WC_alt"]], futureLayersCCSM4rcp60Cut)
    names(futureLayersCCSM4rcp60Cut)
    
futureLayersCCSM4rcp60Cut <- subset(futureLayersCCSM4rcp60Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))    
    names(layersCut) == names(futureLayersCCSM4rcp60Cut)

#Cut and fix cc85
futureLayersCCSM4rcp85Cut <- mask(crop(futureLayersCCSM4rcp85, studyArea), studyArea)
    plot(futureLayersCCSM4rcp85Cut)
        names(futureLayersCCSM4rcp85Cut) <- gsub("_cc85_2070", "", names(futureLayersCCSM4rcp85Cut))
            names(futureLayersCCSM4rcp85Cut)

futureLayersCCSM4rcp85Cut <- stack(layersCut[["WC_alt"]], futureLayersCCSM4rcp85Cut)
    names(futureLayersCCSM4rcp85Cut)

futureLayersCCSM4rcp85Cut <- subset(futureLayersCCSM4rcp85Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))    
    names(layersCut) == names(futureLayersCCSM4rcp85Cut)

#Cut and fix MIROC26
futureLayersMIROCrcp26Cut <- mask(crop(futureLayersMIROCrcp26, studyArea), studyArea)
    plot(futureLayersMIROCrcp26Cut)
        names(futureLayersMIROCrcp26Cut) <- gsub("_mr26_2070", "", names(futureLayersMIROCrcp26Cut))
            names(futureLayersMIROCrcp26Cut)

futureLayersMIROCrcp26Cut <- stack(layersCut[["WC_alt"]], futureLayersMIROCrcp26Cut)
    names(futureLayersMIROCrcp26Cut)

futureLayersMIROCrcp26Cut <- subset(futureLayersMIROCrcp26Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))        
    names(layersCut) == names(futureLayersMIROCrcp26Cut)

#Cut and fix MIROC45
futureLayersMIROCrcp45Cut <- mask(crop(futureLayersMIROCrcp45, studyArea), studyArea)
    plot(futureLayersMIROCrcp45Cut)
        names(futureLayersMIROCrcp45Cut) <- gsub("_mr45_2070", "", names(futureLayersMIROCrcp45Cut))
            names(futureLayersMIROCrcp45Cut)

futureLayersMIROCrcp45Cut <- stack(layersCut[["WC_alt"]], futureLayersMIROCrcp45Cut)
    names(futureLayersMIROCrcp45Cut)

futureLayersMIROCrcp45Cut <- subset(futureLayersMIROCrcp45Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))  
    names(layersCut) == names(futureLayersMIROCrcp45Cut)

#Cut and fix MIROC60
futureLayersMIROCrcp60Cut <- mask(crop(futureLayersMIROCrcp60, studyArea), studyArea)
    plot(futureLayersMIROCrcp60Cut)
        names(futureLayersMIROCrcp60Cut) <- gsub("_mr60_2070", "", names(futureLayersMIROCrcp60Cut))
            names(futureLayersMIROCrcp60Cut)

futureLayersMIROCrcp60Cut <- stack(layersCut[["WC_alt"]], futureLayersMIROCrcp60Cut)
    names(futureLayersMIROCrcp60Cut)

futureLayersMIROCrcp60Cut <- subset(futureLayersMIROCrcp60Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))  
    names(layersCut) == names(futureLayersMIROCrcp60Cut)

#Cut and fix MIROC85
futureLayersMIROCrcp85Cut <- mask(crop(futureLayersMIROCrcp85, studyArea), studyArea)
    plot(futureLayersMIROCrcp85Cut)
        names(futureLayersMIROCrcp85Cut) <- gsub("_mr85_2070", "", names(futureLayersMIROCrcp85Cut))
            names(futureLayersMIROCrcp85Cut)

futureLayersMIROCrcp85Cut <- stack(layersCut[["WC_alt"]], futureLayersMIROCrcp85Cut)
    names(futureLayersMIROCrcp85Cut)

futureLayersMIROCrcp85Cut <- subset(futureLayersMIROCrcp85Cut, c(1, 2, 11, 12, 13, 14, 15, 16, 17, 3, 4, 5, 6, 7, 8, 9, 10))  
    names(layersCut) == names(futureLayersMIROCrcp85Cut)

    
####Look at future changes in probability for the CCSM4 and MIROC models####
#These plots were not used in the paper and are just for visualization here
    
#This will plot all RCP probability values for the CCSM4 (2070) 
predFutureGLMCCSM426 <- predict(futureLayersCCSM4rcp26Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMCCSM445 <- predict(futureLayersCCSM4rcp45Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMCCSM460 <- predict(futureLayersCCSM4rcp60Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMCCSM485 <- predict(futureLayersCCSM4rcp85Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predPresentGLM <- predict(layersCut, modsMultGLMNew$models[["Plci"]], type = "response")
    
#note that I turned the legends off here, but they are all set to have the same scale
par(mfrow = c(3, 2), mar = c(2, 2, 1, 1))
    plot(predPresentGLM, 
         main = "GLM Present (2001-2020)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Present
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMCCSM426, 
         main = "GLM Future (2070CE) (CCSM4 RCP26 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMCCSM445, 
         main = "GLM Future (2070CE) (CCSM4 RCP45 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMCCSM460, 
         main = "GLM Future (2070CE) (CCSM4 RCP60 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend =  FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMCCSM485, 
         main = "GLM Future (2070CE) (CCSM4 RCP85 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = TRUE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
          

#This will plot all RCP probability values for the MIROC (2070) 
predFutureGLMMIROC26 <- predict(futureLayersMIROCrcp26Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMMIROC45 <- predict(futureLayersMIROCrcp45Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMMIROC60 <- predict(futureLayersMIROCrcp60Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predFutureGLMMIROC85 <- predict(futureLayersMIROCrcp85Cut, modsMultGLMNew$models[["Plci"]], type = "response")
predPresentGLM <- predict(layersCut, modsMultGLMNew$models[["Plci"]], type = "response")
    
#note that I turned the legends off here, but they are all set to have the same scale
par(mfrow = c(3, 2), mar = c(2, 2, 1, 1))
    plot(predPresentGLM, 
         main = "GLM Present (2001-2020)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Present
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMMIROC26, 
         main = "GLM Future (2070CE) (MIROC-ESM RCP26 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMMIROC45, 
         main = "GLM Future (2070CE) (MIROC-ESM RCP45 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMMIROC60, 
         main = "GLM Future (2070CE) (MIROC-ESM RCP60 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = FALSE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
    plot(predFutureGLMMIROC85, 
         main = "GLM Future (2070CE) (MIROC-ESM RCP85 model)", 
         breaks = seq(0, 0.65, by = 0.001), 
         col = terrain.colors(650, rev = TRUE),
         legend = TRUE) #Future
          plot(statesMap, border = "grey", add = TRUE)
          plot(studyArea, add = TRUE)
  
          
####Future Favorability####
#Incorporate favorability from the multGLM present model for each GHG model (CCSM4 RCP 26)
futureDatCCSM4rcp26 <- gridRecords_mult(rst = futureLayersCCSM4rcp26Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      #predict using the multGLM model
      futureDatCCSM4rcp26$spc_GLM_P_CCSM4rcp26 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatCCSM4rcp26, 
                                                          type = "response")

      #the prevalence still needs to be the one in the model training data to get future favorability
      futureDatCCSM4rcp26$spc_GLM_F_CCSM4rcp26 <- with(futureDatCCSM4rcp26, 
                                                      Fav(pred = spc_GLM_P_CCSM4rcp26, sample.preval = prevalence(datTrainNew[ , spc])))  

      #then convert to spatial dataset and map both favorability and probability
      futureDatSpatialCCSM4rcp26 <- futureDatCCSM4rcp26
            coordinates(futureDatSpatialCCSM4rcp26) <- futureDatSpatialCCSM4rcp26[ , c("x", "y")]
      spplot(futureDatSpatialCCSM4rcp26, zcol = c("spc_GLM_P_CCSM4rcp26", "spc_GLM_F_CCSM4rcp26"), cex = 0.3)          
          
#Incorporate favorability from the multGLM present model for each GHG model (CCSM4 RCP 45)
futureDatCCSM4rcp45 <- gridRecords_mult(rst = futureLayersCCSM4rcp45Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
      futureDatCCSM4rcp45$spc_GLM_P_CCSM4rcp45 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatCCSM4rcp45, type = "response")
      futureDatCCSM4rcp45$spc_GLM_F_CCSM4rcp45 <- with(futureDatCCSM4rcp45, 
                                                       Fav(pred = spc_GLM_P_CCSM4rcp45, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialCCSM4rcp45 <- futureDatCCSM4rcp45
            coordinates(futureDatSpatialCCSM4rcp45) <- futureDatSpatialCCSM4rcp45[ , c("x", "y")]
      spplot(futureDatSpatialCCSM4rcp45, zcol = c("spc_GLM_P_CCSM4rcp45", "spc_GLM_F_CCSM4rcp45"), cex = 0.3)         

#Incorporate favorability from the multGLM present model for each GHG model (CCSM4 RCP 60)
futureDatCCSM4rcp60 <- gridRecords_mult(rst = futureLayersCCSM4rcp60Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatCCSM4rcp60$spc_GLM_P_CCSM4rcp60 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatCCSM4rcp60, type = "response")
      futureDatCCSM4rcp60$spc_GLM_F_CCSM4rcp60 <- with(futureDatCCSM4rcp60, 
                                                       Fav(pred = spc_GLM_P_CCSM4rcp60, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialCCSM4rcp60 <- futureDatCCSM4rcp60
            coordinates(futureDatSpatialCCSM4rcp60) <- futureDatSpatialCCSM4rcp60[ , c("x", "y")]
      spplot(futureDatSpatialCCSM4rcp60, zcol = c("spc_GLM_P_CCSM4rcp60", "spc_GLM_F_CCSM4rcp60"), cex = 0.3)   

#Incorporate favorability from the multGLM present model for each GHG model (CCSM4 RCP 85)
futureDatCCSM4rcp85 <- gridRecords_mult(rst = futureLayersCCSM4rcp85Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatCCSM4rcp85$spc_GLM_P_CCSM4rcp85 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatCCSM4rcp85, type = "response")
      futureDatCCSM4rcp85$spc_GLM_F_CCSM4rcp85 <- with(futureDatCCSM4rcp85, 
                                                       Fav(pred = spc_GLM_P_CCSM4rcp85, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialCCSM4rcp85 <- futureDatCCSM4rcp85
            coordinates(futureDatSpatialCCSM4rcp85) <- futureDatSpatialCCSM4rcp85[ , c("x", "y")]
      spplot(futureDatSpatialCCSM4rcp85, zcol = c("spc_GLM_P_CCSM4rcp85", "spc_GLM_F_CCSM4rcp85"), cex = 0.3)


#Now do it for MIROC
#Incorporate favorability for each GHG model (MIROC RCP 26)
futureDatMIROCrcp26 <- gridRecords_mult(rst = futureLayersMIROCrcp26Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatMIROCrcp26$spc_GLM_P_MIROCrcp26 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatMIROCrcp26, type = "response")
      futureDatMIROCrcp26$spc_GLM_F_MIROCrcp26 <- with(futureDatMIROCrcp26, 
                                                       Fav(pred = spc_GLM_P_MIROCrcp26, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialMIROCrcp26 <- futureDatMIROCrcp26
            coordinates(futureDatSpatialMIROCrcp26) <- futureDatSpatialMIROCrcp26[ , c("x", "y")]
      spplot(futureDatSpatialMIROCrcp26, zcol = c("spc_GLM_P_MIROCrcp26", "spc_GLM_F_MIROCrcp26"), cex = 0.3)          
          
#Incorporate favorability for each GHG model (MIROC RCP 45)
futureDatMIROCrcp45 <- gridRecords_mult(rst = futureLayersMIROCrcp45Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatMIROCrcp45$spc_GLM_P_MIROCrcp45 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatMIROCrcp45, type = "response")
      futureDatMIROCrcp45$spc_GLM_F_MIROCrcp45 <- with(futureDatMIROCrcp45, 
                                                       Fav(pred = spc_GLM_P_MIROCrcp45, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialMIROCrcp45 <- futureDatMIROCrcp45
            coordinates(futureDatSpatialMIROCrcp45) <- futureDatSpatialMIROCrcp45[ , c("x", "y")]
      spplot(futureDatSpatialMIROCrcp45, zcol = c("spc_GLM_P_MIROCrcp45", "spc_GLM_F_MIROCrcp45"), cex = 0.3)         

#Incorporate favorability for each GHG model (MIROC RCP 60)
futureDatMIROCrcp60 <- gridRecords_mult(rst = futureLayersMIROCrcp60Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatMIROCrcp60$spc_GLM_P_MIROCrcp60 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatMIROCrcp60, type = "response")
      futureDatMIROCrcp60$spc_GLM_F_MIROCrcp60 <- with(futureDatMIROCrcp60, 
                                                       Fav(pred = spc_GLM_P_MIROCrcp60, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialMIROCrcp60 <- futureDatMIROCrcp60
            coordinates(futureDatSpatialMIROCrcp60) <- futureDatSpatialMIROCrcp60[ , c("x", "y")]
      spplot(futureDatSpatialMIROCrcp60, zcol = c("spc_GLM_P_MIROCrcp60", "spc_GLM_F_MIROCrcp60"), cex = 0.3)   

#Incorporate favorability for each GHG model (MIROC RCP 85)
futureDatMIROCrcp85 <- gridRecords_mult(rst = futureLayersMIROCrcp85Cut, 
                              sp.data = as.data.frame(salCoordsNew), 
                              sp.col = "speciesCode", 
                              coord.cols = c("decimalLongitude", "decimalLatitude"))
          
      futureDatMIROCrcp85$spc_GLM_P_MIROCrcp85 <- predict(modsMultGLMNew$models[["Plci"]], 
                                                          newdata = futureDatMIROCrcp85, type = "response")
      futureDatMIROCrcp85$spc_GLM_F_MIROCrcp85 <- with(futureDatMIROCrcp85, 
                                                       Fav(pred = spc_GLM_P_MIROCrcp85, sample.preval = prevalence(datTrainNew[ , spc])))  
      futureDatSpatialMIROCrcp85 <- futureDatMIROCrcp85
            coordinates(futureDatSpatialMIROCrcp85) <- futureDatSpatialMIROCrcp85[ , c("x", "y")]
      spplot(futureDatSpatialMIROCrcp85, zcol = c("spc_GLM_P_MIROCrcp85", "spc_GLM_F_MIROCrcp85"), cex = 0.3)

#Make one dataset with all of the future data
futureDatSpatial <- futureDatSpatialCCSM4rcp26[, 1:4]
    futureDatSpatial$present_F <- datSpatialNew$Plci_F
    futureDatSpatial$spc_GLM_F_CCSM4rcp26 <- futureDatSpatialCCSM4rcp26$spc_GLM_F_CCSM4rcp26
    futureDatSpatial$spc_GLM_F_CCSM4rcp45 <- futureDatSpatialCCSM4rcp45$spc_GLM_F_CCSM4rcp45
    futureDatSpatial$spc_GLM_F_CCSM4rcp60 <- futureDatSpatialCCSM4rcp60$spc_GLM_F_CCSM4rcp60
    futureDatSpatial$spc_GLM_F_CCSM4rcp85 <- futureDatSpatialCCSM4rcp85$spc_GLM_F_CCSM4rcp85
    futureDatSpatial$spc_GLM_F_MIROCrcp26 <- futureDatSpatialMIROCrcp26$spc_GLM_F_MIROCrcp26
    futureDatSpatial$spc_GLM_F_MIROCrcp45 <- futureDatSpatialMIROCrcp45$spc_GLM_F_MIROCrcp45
    futureDatSpatial$spc_GLM_F_MIROCrcp60 <- futureDatSpatialMIROCrcp60$spc_GLM_F_MIROCrcp60
    futureDatSpatial$spc_GLM_F_MIROCrcp85 <- futureDatSpatialMIROCrcp85$spc_GLM_F_MIROCrcp85

#Figure S6, Plot CCSM4 dataset (10 x 7 PDFs)
spplot(futureDatSpatial, 
            zcol = c("spc_GLM_F_CCSM4rcp60",
                     "spc_GLM_F_CCSM4rcp85",
                     "spc_GLM_F_CCSM4rcp26",
                     "spc_GLM_F_CCSM4rcp45",
                     "present_F"), 
            cex = 0.3,
            col.regions = terrain.colors(1000, rev = TRUE))

#Figure S7, Plot MIROC dataset (10 x 7 PDFs)
spplot(futureDatSpatial, 
            zcol = c("spc_GLM_F_MIROCrcp60",
                     "spc_GLM_F_MIROCrcp85",
                     "spc_GLM_F_MIROCrcp26",
                     "spc_GLM_F_MIROCrcp45",
                     "present_F"), 
            cex = 0.3,
            col.regions = terrain.colors(1000, rev = TRUE))    

##Figure 3      
spplot(futureDatSpatial, 
       zcol = c("spc_GLM_F_MIROCrcp26",
                "spc_GLM_F_MIROCrcp85",
                "spc_GLM_F_CCSM4rcp26",
                "spc_GLM_F_CCSM4rcp85"), 
       cex = 0.3,
       col.regions = terrain.colors(1000, rev = TRUE))

#Make blank map of US/Canada
dev.off()

plot(predFutureGLMMIROC85, 
     main = "GLM Future (2070CE) (MIROC85 model)", 
     breaks = seq(0, 0.001, by = 0.001), 
     col = terrain.colors(0.001, rev = TRUE),
     legend = TRUE) #Future
    plot(statesMap, border = "black", add = TRUE, lwd = 1)
    plot(studyArea, add = TRUE, lwd = 2)

#Generate non-spatial dataset
futureDat <- datNew[, 1:4]
    futureDat$present_F <- datNew$Plci_F
    futureDat$spc_GLM_F_CCSM4rcp26 <- futureDatCCSM4rcp26$spc_GLM_F_CCSM4rcp26
    futureDat$spc_GLM_F_CCSM4rcp45 <- futureDatCCSM4rcp45$spc_GLM_F_CCSM4rcp45
    futureDat$spc_GLM_F_CCSM4rcp60 <- futureDatCCSM4rcp60$spc_GLM_F_CCSM4rcp60
    futureDat$spc_GLM_F_CCSM4rcp85 <- futureDatCCSM4rcp85$spc_GLM_F_CCSM4rcp85
    futureDat$spc_GLM_F_MIROCrcp26 <- futureDatMIROCrcp26$spc_GLM_F_MIROCrcp26
    futureDat$spc_GLM_F_MIROCrcp45 <- futureDatMIROCrcp45$spc_GLM_F_MIROCrcp45
    futureDat$spc_GLM_F_MIROCrcp60 <- futureDatMIROCrcp60$spc_GLM_F_MIROCrcp60
    futureDat$spc_GLM_F_MIROCrcp85 <- futureDatMIROCrcp85$spc_GLM_F_MIROCrcp85

    
####Fuzzy Overlap Operations####

#Compare present favorability with CCSM4 RCP 85 and MIROC RCP 85
futureDat$CCSM4expansion <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_CCSM4rcp85"), 
                                            op = "expansion")
futureDat$CCSM4contraction <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_CCSM4rcp85"), 
                                            op = "contraction")
futureDat$CCSM4maintenance <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_CCSM4rcp85"), 
                                            op = "maintenance")    
futureDat$CCSM4change <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_CCSM4rcp85"), 
                                            op = "change")
futureDat$MIROCexpansion <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_MIROCrcp85"), 
                                            op = "expansion")
futureDat$MIROCcontraction <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_MIROCrcp85"), 
                                            op = "contraction")
futureDat$MIROCmaintenance <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_MIROCrcp85"), 
                                            op = "maintenance")    
futureDat$MIROCchange <- fuzzyOverlay(data = futureDat, 
                                            overlay.cols = c("present_F", "spc_GLM_F_MIROCrcp85"), 
                                            op = "change")

#Change to a spatial dataframe so that we can make maps
futureDatFuzzyMetrics <- futureDat
    coordinates(futureDatFuzzyMetrics) <- futureDatFuzzyMetrics[ , c("x", "y")]
    names(futureDatFuzzyMetrics)

#Figure S8 Maps    
spplot(futureDatFuzzyMetrics, 
       c("CCSM4expansion"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "CCSM4 Expansion")
spplot(futureDatFuzzyMetrics, 
       c("CCSM4contraction"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "CCSM4 Contraction")
spplot(futureDatFuzzyMetrics, 
       c("MIROCexpansion"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "MIROC Expansion")
spplot(futureDatFuzzyMetrics, 
       c("MIROCcontraction"),
       col.regions = terrain.colors(1000, rev = TRUE),
       main = "MIROC Contraction")

#Calculate model overlap (Schoener's D, Warren's I, Hellinger Distance), Table S3
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp26)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp45)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp60)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp85)

modOverlap(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp26)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp45)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp60)
modOverlap(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp85)

#Compare overlap of CCSM4 and MIROC models at RCP 85
modOverlap(futureDat$spc_GLM_F_CCSM4rcp26, futureDat$spc_GLM_F_MIROCrcp26)
modOverlap(futureDat$spc_GLM_F_CCSM4rcp45, futureDat$spc_GLM_F_MIROCrcp45)
modOverlap(futureDat$spc_GLM_F_CCSM4rcp60, futureDat$spc_GLM_F_MIROCrcp60)
modOverlap(futureDat$spc_GLM_F_CCSM4rcp85, futureDat$spc_GLM_F_MIROCrcp85) #Table S3

#Additional overlap metrics (Jaccard and Baroni), Table S3
#CCSM4 Models
with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp26, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp26, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp45, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp45, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp60, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp60, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp85, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_CCSM4rcp85, method = "Baroni"))

#MIROC Models
with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp26, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp26, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp45, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp45, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp60, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp60, method = "Baroni"))

with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp85, method = "Jaccard"))
with(futureDat, fuzSim(present_F, spc_GLM_F_MIROCrcp85, method = "Baroni"))


#Calculate fuzzy range change measures, Gain is gained presences, Loss is lost presences
#Presences that remained as presences, Absences that remained as absences, Balance is balance of gained and lost presences

#Figure 4, Plotted as 6 x 10
par(mfrow = c(2, 4), mar = c(10, 3, 2, 3))

fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp26,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "CCSM4 RCP2.6",
                 xlab = "",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp45,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "CCSM4 RCP4.5",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp60,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "CCSM4 RCP6.0",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_CCSM4rcp85,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "CCSM4 RCP8.5",
                 las = 2)

fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp26,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "MIROC RCP2.6",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp45,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "MIROC RCP4.5",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp60,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "MIROC RCP6.0",
                 las = 2)
fuzzyRangeChange(futureDat$present_F, futureDat$spc_GLM_F_MIROCrcp85,
                 measures = c("Gain", "Loss", "Stable presence", "Stable absence", "Balance"),
                 ylim = c(-0.4, 0.4),
                 main = "MIROC RCP8.5",
                 las = 2)


####Check for Spatial Autocorrelation####

#Using the multGLM model for the period 2001-2020, add in the coords for the train data and the residuals
#Note that this will take a very long time to run and needs to be run on a computer with 128GB of RAM

datResids <- cbind(modsMultGLMNew$models$Plci$data, datTrainNew$x, datTrainNew$y, modsMultGLMNew$models$Plci$residuals)
    head(datResids)
    names(datResids)[9] <- "x"
    names(datResids)[10] <- "y"
    names(datResids)[11] <- "resids"

library(ape) # Moran's Index
library(ncf) # Correlogram

matriz.dist <- as.matrix(dist(cbind(datResids$x, datResids$y))) #Make a distance matrix
    head(matriz.dist)
    inv_matriz.dist <- 1/matriz.dist

# Define the  DIAGONAL as 0
diag(inv_matriz.dist) <- 0

mPlci <- Moran.I(datResids$Plci, inv_matriz.dist)
mRes <- Moran.I(datResids$res, inv_matriz.dist) 







