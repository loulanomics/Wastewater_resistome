####################################
### read & organize ARG ddPCR output
### Lou LaMartina, May 17, 2022
####################################


setwd("~/Desktop/Lab/Projects/Resistance")


#################
### prep data ###
#################

# find files
output_files <- list.files("./Data/ddPCR_output", pattern = ".csv", full.names = TRUE, recursive = T)


# read in all ddPCR output files
output_files.ls <- lapply(output_files, function(i) read.csv(i, header = F))


# name each dataset
names(output_files.ls) <- gsub(" ", "", gsub("-", "", gsub("_", "", gsub(".csv", "", basename(output_files)))))


# there is a "ghost" last column (BL if opened in excel)
cols <- c(as.character(output_files.ls[[1]][1,])[1:63], "ghost")


# remove it
for (i in names(output_files.ls)) {
  j <- output_files.ls[[i]]
  colnames(j) <- cols[1:ncol(j)]
  j <- j[! grepl("Well", j$Well),]
  j <- j[colnames(j) != "ghost"]
  output_files.ls[[i]] <- j
  rm(j)
}


# add plate name
for (i in names(output_files.ls)) {
  output_files.ls[[i]]$Plate <- i
}


# combine into one data frame
output_all <- do.call(rbind, output_files.ls)




#######################
### subset by study ###
#######################

# remove irrelevant assays
output <- subset(output_all, ! Plate %in% c("Lou21819", "Lou2519","Lou91520", "Lou91620", "Lou122017"))
output <- output[output$Target %in% c("sul1", "tetO", "int1", "qnrS", "tetA", "tetW", "intI1", "GES-1", "KPC"),]


# fix intI1
output$Target[output$Target == "int1"] <- "intI1"


# remove plasmidome
plasmid_grep <- paste(c("plasmid", "total", "gblock"), collapse = "|")
output <- output[grepl(plasmid_grep, output$Sample) == F,]


# pattern matching sample names
ts_grep <- paste(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 
                       collapse = "|")
cities_grep <- paste(state.abb, collapse = "|")
neigh_grep <- paste(c("EG", "MI", "NB", "INF", "EFF"), collapse = "|")


# time series assays
ts_assays <- output$Sample[grepl(ts_grep, output$Sample, ignore.case = T)]
ts_assays <- ts_assays[grepl(" ", ts_assays) == F]
sort(unique(ts_assays))


# cities assays
cities_assays <- output$Sample[grepl(cities_grep, output$Sample, ignore.case = T)]
cities_assays <- cities_assays[! cities_assays %in% ts_assays]
cities_assays <- cities_assays[grepl(neigh_grep, cities_assays, ignore.case = T) == F]
sort(unique(cities_assays))


# neighborhood assays
neigh_assays <- output$Sample[grepl(neigh_grep, output$Sample, ignore.case = T)]
sort(unique(neigh_assays))


# subset
ts_output <- subset(output, Sample %in% ts_assays)
cities_output <- subset(output, Sample %in% cities_assays)
neigh_output <- subset(output, Sample %in% neigh_assays)


# leftovers should be just experimental serial dilutions, negative controls
leftover_output <- output[! output$Sample %in% c(ts_assays, cities_assays, neigh_assays),]
unique(leftover_output$Sample)


# fix neighborhood names
neigh_output$Sample[neigh_output$Sample == "MI"] <- "MI1A"
neigh_output$Sample[neigh_output$Sample == "EG"] <- "EG1A"
neigh_output$Sample[neigh_output$Sample == "NB"] <- "NB1A"
neigh_output <- neigh_output[grepl(":", neigh_output$Sample) == F,] # remove dilutions, they failed anyway




############
### save ###
############

write.csv(ts_output, "./RData/01_timeSeries_ddPCR_output.csv", row.names = F, na = "")
write.csv(cities_output, "./RData/01_cities_ddPCR_output.csv", row.names = F, na = "")
write.csv(neigh_output, "./RData/01_neighborhoods_ddPCR_output.csv", row.names = F, na = "")

