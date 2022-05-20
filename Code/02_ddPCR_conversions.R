########################################
### convert ddPCR output from copies per
### reactions to copies per mL sample,
### and remove bad assays
### Lou LaMartina, May 17, 2022
########################################


setwd("~/Desktop/Lab/Projects/Resistance")


#################
### load data ###
#################

# ddpcr results (see 01_ddPCR_dataPrep.R)
time_ddpcr <- read.csv("./RData/01_timeSeries_ddPCR_output.csv")
cities_ddpcr <- read.csv("./RData/01_cities_ddPCR_output.csv")
neigh_ddpcr <- read.csv("./RData/01_neighborhoods_ddPCR_output.csv")


# sample set info
time_info <- read.csv("./Data/Metadata/TimeSeries_sample_info.csv")
cities_info <- read.csv("./Data/Metadata/Cities_sample_info.csv")
neigh_info <- read.csv("./Data/Metadata/Neighborhoods_sample_info.csv")


# additional assay information (eg, dilution in sample name)
cities_assays <- read.csv("./Data/Metadata/Cities_ddPCR_assays.csv")




#################
### prep data ###
#################

# combine
cities_ddpcr <- merge(cities_assays, cities_ddpcr, by = "Sample")
time_ddpcr <- merge(time_info[c("Sample", "Sample_name")], time_ddpcr, by = "Sample")
neigh_ddpcr <- merge(neigh_info[c("Sample", "Sample_name")], neigh_ddpcr, by = "Sample")


# make assay variable
time_ddpcr <- data.frame(Assay = paste(time_ddpcr$Sample_name, time_ddpcr$Target), time_ddpcr)
cities_ddpcr <- data.frame(Assay = paste(cities_ddpcr$Sample_name, cities_ddpcr$Target), cities_ddpcr)
neigh_ddpcr <- data.frame(Assay = paste(neigh_ddpcr$Sample_name, neigh_ddpcr$Target), neigh_ddpcr)


# add elutions -
# add volume eluted (ul) during extraction -
# always 150 except when katie extracted south shores for me
cities_ddpcr$elute_uL <- 150
neigh_ddpcr$elute_uL <- 150
time_ddpcr$elute_uL <- 150
time_ddpcr$elute_uL[grep("South", time_ddpcr$Sample_name)] <- 100


# add dilution - always 100 except for cities, specified already in "assays" df
time_ddpcr$dilute_x <- 100
neigh_ddpcr$dilute_x <- 100


# keep relevant data
ddpcr_stats <- rbind(time_ddpcr[colnames(time_ddpcr) %in% c("Assay", "Sample_name", "Target", "dilute_x", 
                                                            "elute_uL", "CopiesPer20uLWell", "Plate")],
                 cities_ddpcr[colnames(cities_ddpcr) %in% c("Assay", "Sample_name", "Target", "dilute_x", 
                                                            "elute_uL", "CopiesPer20uLWell", "Plate")],
                   neigh_ddpcr[colnames(neigh_ddpcr) %in% c("Assay", "Sample_name", "Target", "dilute_x", 
                                                            "elute_uL", "CopiesPer20uLWell", "Plate")])




#################
### mean & CV ###
#################

# calculate mean & cv
means.ls <- list()
cv.ls <- list()

for(i in ddpcr_stats$Assay){
  means.ls[[i]] <- mean(subset(ddpcr_stats, Assay == i)$CopiesPer20uLWell, na.rm = TRUE)
  cv.ls[[i]] <- sd(subset(ddpcr_stats, Assay == i)$CopiesPer20uLWell, na.rm = TRUE) /
    mean(subset(ddpcr_stats, Assay == i)$CopiesPer20uLWell, na.rm = TRUE)
}


# combine
mean_cv <- data.frame(Assay = names(means.ls), meanPerWell = unlist(means.ls), CV = unlist(cv.ls))
ddpcr_stats <- merge(ddpcr_stats, mean_cv, by = "Assay")




#####################
### copies per mL ###
#####################

# calculate
ddpcr_stats$copiespermL <- 
  
  # mean of duplicate assays, times
  ddpcr_stats$meanPerWell *
  
  # dilution factor, divided by vol added to pcr (always 2ul), times
  ddpcr_stats$dilute_x / 2  *
  
  # DNA extract elution volume, divided by sample vol filtered (always 25 ml)
  ddpcr_stats$elute_uL / 25




###############
### cleanup ###
###############

# remove failed assays
nrow(ddpcr_stats) # [1] 1184
ddpcr_stats <- ddpcr_stats[is.na(ddpcr_stats$CopiesPer20uLWell) == F,]
nrow(ddpcr_stats) # [1] 1123


# subset assays without successful duplicates
fails <- ddpcr_stats[ddpcr_stats$Assay %in% names(which(table(ddpcr_stats$Assay) == 1)),]
fails$Redo <- "no duplicate"
ddpcr_stats <- ddpcr_stats[ddpcr_stats$Assay %in% names(which(table(ddpcr_stats$Assay) > 1)),]
nrow(ddpcr_stats) # [1] 1110


# remove assays with CV > 0.3
hiCV <- ddpcr_stats[ddpcr_stats$Assay %in% mean_cv$Assay[mean_cv$CV > 0.3],]
hiCV$Redo <- "high CV"
fails <- rbind(fails, hiCV)
ddpcr_stats <- ddpcr_stats[! ddpcr_stats$Assay %in% mean_cv$Assay[mean_cv$CV > 0.3],]
nrow(ddpcr_stats) # [1] 1011




############
### save ###
############

write.csv(ddpcr_stats, "./RData/02_ddPCR_conversions_all.csv", row.names = F, na = "")

write.csv(unique(ddpcr_stats[! colnames(ddpcr_stats) %in% c("Plate", "CopiesPer20uLWell")]), 
          "./RData/02_ddPCR_conversions_uniq.csv", row.names = F, na = "")

# save this to plan redos later
# adding "x" to names so excel stops converting names to dates
write.csv(fails, "./RData/failed_ddPCR.csv", row.names = F, na = "")
length(unique(fails$Sample_name)) # [1] 36
length(unique(fails$Target)) # [1] 7
length(unique(fails$Assay)) # [1] 46

