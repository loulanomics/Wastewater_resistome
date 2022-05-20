###########################################
### Cleaning up output files from programs:
### - RGI: resistance gene identifier
### - mobileOG: mobile element database
### - Kaiju: taxonomy classifier
### - PlasmidFinder: plasmid identifier
### Lou LaMartina, May 18, 2022
###########################################


setwd("~/Desktop/Lab/Projects/Resistance")

# sample info
info <- read.csv("./Data/Metadata/Metagenomes_sample_info.csv")


###########
### RGI ###
###########
# https://card.mcmaster.ca/analyze/rgi

# read files
arg_files <- list.files("./Data/output_RGI", pattern = ".txt", full.names = TRUE)


# read all files from above into data frames
arg_files.ls <- list()
for (i in arg_files) {
  arg_files.ls[[i]] <- read.delim(i, sep = "\t")
}


# combine into single data frame
arg_results <- dplyr::bind_rows(arg_files.ls, .id = "column_label")


# change to NA
arg_results[arg_results == "n/a"] <- NA
arg_results[arg_results == ""] <- NA


# change to underscore
colnames(arg_results) <- gsub("\\.", "_", colnames(arg_results))


# add sample, dataset, biome, location variables
arg_results$File_name <- paste0(sapply(strsplit(gsub("RGI_", "", 
                                                     basename(arg_results$column_label)), "_"), '[', 1), "_",
                                sapply(strsplit(gsub("RGI_", "", 
                                                     basename(arg_results$column_label)), "_"), '[', 2))
arg_results <- merge(info, arg_results, by = "File_name")


# remove sequences to make smaller file
arg_results <- arg_results[! colnames(arg_results) %in%
                             c("Predicted_DNA", "Predicted_Protein", "CARD_Protein_Sequence")]


# split contig from contig_ORF
arg_results$Contig_ORF <- arg_results$Contig
arg_results$Contig <- paste0(sapply(strsplit(arg_results$Contig_ORF, "_"), '[', 1), "_",
                             sapply(strsplit(arg_results$Contig_ORF, "_"), '[', 2), "_",
                             sapply(strsplit(arg_results$Contig_ORF, "_"), '[', 3))


# remove empties
arg_results <- arg_results[is.na(arg_results$Drug_Class) == F,]


# clean up
arg_results <- arg_results[! colnames(arg_results) %in% c("column_label", "Note")]
rm(arg_files.ls, arg_files)




###################
### mobileOG-db ###
###################
# https://mobileogdb.flsi.cloud.vt.edu

# read files
mge_files <- list.files("./Data/output_mobileOG", pattern = "mobileOG.Alignment.Out.csv", 
                        full.names = TRUE, recursive = TRUE)


# read all files from above into data frames
mge_files.ls <- list()
for (i in mge_files) {
  mge_files.ls[[i]] <- read.csv(i)
}


# combine into single data frame
mge_results <- do.call(rbind, mge_files.ls)


# remove empty column
mge_results <- mge_results[-1]


# change to underscore
colnames(mge_results) <- gsub("\\.", "_", colnames(mge_results))


# add sample, dataset, biome, location variables
mge_results$File_name <- paste0(sapply(strsplit(mge_results$Specific_Contig, "_"), '[', 1), "_",
                                sapply(strsplit(mge_results$Specific_Contig, "_"), '[', 2))
mge_results <- merge(info, mge_results, by = "File_name")


# split by gene category - there are several but this one seems most relevant to ARGs
int_results <- mge_results[mge_results$Major_mobileOG_Category == "integration/excision",]  


# clean up
rm(mge_results, mge_files.ls, mge_files)


# replace NA minor categories with major categories -
# note that "integration/excision" was the major category,
# and any minor categories that were blank are now this
int_results$Minor_mobileOG_Category[is.na(int_results$Minor_mobileOG_Category)] <- 
  int_results$Major_mobileOG_Category[is.na(int_results$Minor_mobileOG_Category)]




#############
### kaiju ###
#############
# https://kaiju.binf.ku.dk

###################
### alignment stats

# read files
kaiju_files <- list.files("./Data/output_Kaiju", pattern = "align.txt", 
                                full.names = TRUE, recursive = TRUE)


# read all files from above into data frames - ignore warnings
kaiju_files.ls <- list()
for (i in kaiju_files) {
  kaiju_files.ls[[i]] <- try(read.delim(i, fill = T, header = F))
}


# combine into single data frame
kaiju_results <- do.call(rbind, kaiju_files.ls)


# add column names
colnames(kaiju_results) <- c("classed_or_unclassed", "Contig", "NCBI_ID", "Score",
                                   "taxon_ID", "accesion_num", "prot_seq_match")
rownames(kaiju_results)<- 1:nrow(kaiju_results)


# get file names
kaiju_results <- data.frame(File_name = 
                                    paste0(sapply(strsplit(kaiju_results$Contig, "_"), '[', 1), "_", 
                                           sapply(strsplit(kaiju_results$Contig, "_"), '[', 2)),
                                  kaiju_results)


# change empties to NA
kaiju_results[kaiju_results == ""] <- NA


# remove unclassified
kaiju_results <- kaiju_results[kaiju_results$classed_or_unclassed == "C",]


# # need to associate contigs with taxonomy - see next section & 04_ncbi_IDs.sh
# # why? because the command line version has an add taxonomy option,
# # but i didn't use the command line version (requires downloading
# # enormous ncbi nr database), and the online tool does not have this option.
# # it does give you a file with taxonomy info, but not per-contig
# ncbi <- unique(as.numeric(kaiju_results$NCBI_ID))
# write(ncbi, "./RData/04_output/04_ncbi_data.txt", sep = "\n")


# clean up
rm(kaiju_files, kaiju_files.ls)



#################
### ncbi taxonomy

# read files
ncbi_files <- list.files("./Data/output_NCBI", pattern = ".txt", full.names = TRUE)


# you don't want subgroups, suborders, unclassifieds, etc
remove <- c("group", "Group", "IncertaeSedis", "incertaesedis", 
            "unclassified", "complex", "division", "cluster",
            "Sphaerobacteridae", "Sorangiineae", "Roseiflexineae",
            "Rhizobiaceae", "Rickettsieae", "Cystobacterineae",
            "Nannocystineae", "Sphaerobacterineae", 
            "Wolbachieae", "Ureaplasma", "Chloroflexineae")


# # read all files from above into data frames - takes a minute
# ncbi_files.ls <- list()
# for (i in ncbi_files) {
#   ncbi_files.ls[[i]] <- scan(i, character(), sep = ",", quiet = T, quote = "")
#   if (TRUE %in% grepl(paste(remove, collapse = "|"), fixed = F, ncbi_files.ls[[i]])) {
#     ncbi_files.ls[[i]] <- ncbi_files.ls[[i]][-grep(paste(remove, collapse = "|"), fixed = F, ncbi_files.ls[[i]])]
#   }
#   length(ncbi_files.ls[[i]]) <- 20
# }
# saveRDS(ncbi_files.ls, "./RData/04_output/04_ncbi_files_list.RData")
ncbi_files.ls <- readRDS("./RData/04_output/04_ncbi_files_list.RData")


# combine into single data frame
ncbi_results <- data.frame(do.call(rbind, ncbi_files.ls))


# remove column specifying cells or viruses
ncbi_results <- ncbi_results[-1]


# adjust kindgom, keep bacteria
colnames(ncbi_results)[1:8] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Genus__Species", "Subgroup")
ncbi_results$Kingdom[ncbi_results$Kingdom == "Bacteria(eubacteria)"] <- "Bacteria"
ncbi_results$Kingdom[ncbi_results$Kingdom == "Eukaryota(eucaryotes)"] <- "Eukaryota"
ncbi_results <- ncbi_results[ncbi_results$Kingdom == "Bacteria",]


# remove empty columns
ncbi_results <- ncbi_results[1:8]


# remove "sub" species
ncbi_results$Subgroup[grepl("sub", ncbi_results$Subgroup)] <- NA
ncbi_results$Subgroup[grepl("pv", ncbi_results$Subgroup)] <- NA
ncbi_results$Genus__Species[grepl("sub", ncbi_results$Genus__Species)] <- NA
ncbi_results$Subgroup[grepl("etal", ncbi_results$Subgroup)] <- NA
ncbi_results$Genus__Species[grepl("__", ncbi_results$Genus__Species)] <- 
  sapply(strsplit(ncbi_results$Genus__Species[grepl("__", ncbi_results$Genus__Species)], "__"), '[', 1)
ncbi_results$Genus__Species[is.na(ncbi_results$Genus__Species)] <- ncbi_results$Subgroup[is.na(ncbi_results$Genus__Species)]


# remove host
ncbi_results$Subgroup[ncbi_results$Genus == "ant"] <- NA
ncbi_results$Genus__Species[ncbi_results$Genus == "ant"] <- NA
ncbi_results$Genus[ncbi_results$Genus == "ant"] <- NA
ncbi_results[ncbi_results == "environmentalsamples"] <- NA


# some have [genusName]speciesName or [familyName]genusName,
# so extract info from within the brackets
ncbi_results$Genus[grep("]", ncbi_results$Genus__Species)] <- 
  ncbi_results$Genus__Species[grep("]", ncbi_results$Genus__Species)]

ncbi_results$Genus[grep("]", ncbi_results$Family)] <- 
  ncbi_results$Family[grep("]", ncbi_results$Family)]

ncbi_results$Family[grep("]", ncbi_results$Family)] <- NA

ncbi_results$Genus[grep("]", ncbi_results$Genus)] <-
  gsub("\\[", "", sapply(strsplit(ncbi_results$Genus[grep("]", ncbi_results$Genus)], "]"), '[', 1))

ncbi_results$Genus__Species[grep("]", ncbi_results$Genus__Species)] <-
  gsub("]", "", gsub("\\[", "", ncbi_results$Genus__Species[grep("]", ncbi_results$Genus__Species)]))


# separate genus/species from subgroups
for (i in 1:nrow(ncbi_results)) {
  ncbi_results$Subgroup[i] <- gsub(ncbi_results$Genus__Species[i], "", ncbi_results$Subgroup[i])
}


# separate genus from subgroups
for (i in 1:nrow(ncbi_results)) {
  ncbi_results$Subgroup[i] <- gsub(ncbi_results$Genus[i], "", ncbi_results$Subgroup[i])
}


# remove genus from species
for (i in 1:nrow(ncbi_results)) {
  ncbi_results$Species[i] <- gsub(ncbi_results$Genus[i], "", ncbi_results$Genus__Species[i])
}


# remove redundant columns and reorder
ncbi_results <- ncbi_results[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subgroup")]


# remove "Candidatus"
for (i in colnames(ncbi_results)) {
  ncbi_results[[i]] <- gsub("Candidatus", "", ncbi_results[[i]])
}


# remove paranthesis and anything inside (ususally authors or host info)
for (i in colnames(ncbi_results)) {
  ncbi_results[[i]] <- sapply(strsplit(ncbi_results[[i]], "\\("), '[', 1)
}


# redo removing genus from species now that there are more matches
for (i in 1:nrow(ncbi_results)) {
  ncbi_results$Species[i] <- gsub(ncbi_results$Genus[i], "", ncbi_results$Species[i])
}


# yay! now add NCBI ID
ncbi_results <- data.frame(NCBI_ID = gsub(".txt", "", 
                                          gsub("./Data/output_NCBI/sed-lineage-tax_", "", 
                                               rownames(ncbi_results))), ncbi_results)
kaiju_results <- merge(ncbi_results, kaiju_results, by = "NCBI_ID")


# remove anything unclassified at kingdom or phylum level
kaiju_results <- kaiju_results[is.na(kaiju_results$Kingdom) == F,]
kaiju_results <- kaiju_results[kaiju_results$Kingdom != "",]
kaiju_results <- kaiju_results[is.na(kaiju_results$Phylum) == F,]
kaiju_results <- kaiju_results[kaiju_results$Phylum != "Bacteriacandidatephyla",]


# clean up
kaiju_results <- kaiju_results[colnames(kaiju_results) != "prot_seq_match"]
rm(ncbi, ncbi_files, ncbi_files.ls, ncbi_results, remove)




#####################
### PlasmidFinder ###
#####################
# https://cge.cbs.dtu.dk/services/PlasmidFinder/
# see 04_plasmidFinder.sh

# read files
plas_files <- list.files("./Data/output_plasmidFinder", pattern = ".txt", 
                        full.names = TRUE, recursive = TRUE)


# read all files from above into data frames
# these are huge so just grab the the best alignments
plas_files.ls <- list()
for (i in plas_files) {
  file <- read.delim(i, sep = "\t")
  plas_files.ls[[i]] <- file[file$pident == 100 & file$evalue < 0.05,]
  colnames(plas_files.ls[[i]])[1] <- "qseqid"
  rm(file)
}


# combine into single data frame
plas_results <- do.call(rbind, plas_files.ls)


# add sample, dataset, biome, location variables
plas_results$File_name <- paste0(sapply(strsplit(plas_results$qseqid, "_"), '[', 1), "_",
                                sapply(strsplit(plas_results$qseqid, "_"), '[', 2))
plas_results <- merge(info, plas_results, by = "File_name")


# extract info from name
plasIDs <- data.frame(sseqid = unique(plas_results$sseqid),
                      Gene_name = unique(plas_results$sseqid))


# the last item is the NCBI accession #
plasIDs$Accession <- sapply(strsplit(plasIDs$Gene_name, "__"), '[', 2)
plasIDs2 <- plasIDs[is.na(plasIDs$Accession),]
plasIDs <- plasIDs[is.na(plasIDs$Accession) == F,]

for(i in 1:nrow(plasIDs2)) {
  chars <- unlist(strsplit(plasIDs2$Gene_name[i], "_"))
  j = length(chars)
  plasIDs2$Accession[i] <- chars[j]
}
plasIDs <- rbind(plasIDs, plasIDs2); rm(plasIDs2)


# remove accession # from gene name
for (i in 1:nrow(plasIDs)) {
  plasIDs$Gene_name[i] <- gsub(paste0("__", plasIDs$Accession[i]), "", plasIDs$Gene_name[i])
  plasIDs$Gene_name[i] <- gsub(paste0("_", plasIDs$Accession[i]), "", plasIDs$Gene_name[i])
}


# NC accessions have an underscore, annoyingly
plasIDs$Accession[grep("_NC$", plasIDs$Gene_name)] <- 
  paste0("NC_", plasIDs$Accession[grep("_NC$", plasIDs$Gene_name)])
for (i in 1:nrow(plasIDs)) {
  plasIDs$Gene_name[i] <- gsub(paste0("__", plasIDs$Accession[i]), "", plasIDs$Gene_name[i])
  plasIDs$Gene_name[i] <- gsub(paste0("_", plasIDs$Accession[i]), "", plasIDs$Gene_name[i])
}


# the _number_ is think is an occurrence thing, just remove it
plasIDs$Subtype <- sapply(strsplit(plasIDs$Gene_name, "_[[:digit:]]_"), '[', 2)
plasIDs$Gene_name <- gsub("_[[:digit:]]$", "", plasIDs$Gene_name)


# remove subtype from gene name
for (i in 1:nrow(plasIDs)) {
  plasIDs$Gene_name[i] <- gsub(paste0("_[[:digit:]]_", plasIDs$Subtype[i]), "", plasIDs$Gene_name[i])
}


# remove subtype # from gene name, in paranthesis
for (i in 1:nrow(plasIDs)) {
  plasIDs$Gene_name2[i] <- gsub(plasIDs$Subtype[i], "", plasIDs$Gene_name[i])
}
plasIDs$Gene_name2 <- gsub("\\()", "", plasIDs$Gene_name2)
plasIDs$Gene_name[is.na(plasIDs$Gene_name2) == F] <- plasIDs$Gene_name2[is.na(plasIDs$Gene_name2) == F]
plasIDs <- plasIDs[-5]


# more NC
plasIDs$Subtype[grep("_NC$", plasIDs$Subtype)] <- gsub("_NC$", "", plasIDs$Subtype[grep("_NC$", plasIDs$Subtype)])


# one left with accession
plasIDs$Gene_name[plasIDs$Gene_name == "IncFIB(AP001918)"] <- "IncFIB"


# make categories
plasIDs$Category <- NA
plasIDs$Category[grep("Inc", plasIDs$Gene_name)] <- substr(plasIDs$Gene_name[grep("Inc", plasIDs$Gene_name)], 1, 4)
plasIDs$Category[grep("Col", plasIDs$Gene_name)] <- "Col"
plasIDs$Category[grep("(^p|^P)", plasIDs$Gene_name)] <- "p"
plasIDs$Category[grep("Trf", plasIDs$Gene_name)] <- "Trf"
plasIDs$Category[grep("(^rep|^Rep)", plasIDs$Gene_name)] <- "rep"


# merge
plas_results <- merge(plasIDs, plas_results, by = "sseqid")


# clean up
rm(plas_files.ls, plas_files, plasIDs)




############
### save ###
############

write.csv(arg_results, "./RData/04_output/04_RGI_compiled_data.csv", row.names = F, na = "")
write.csv(int_results, "./RData/04_output/04_integron_compiled_data.csv", row.names = F, na = "")
write.csv(kaiju_results, "./RData/04_output/04_kaiju_compiled_data.csv", row.names = F, na = "")
write.csv(plas_results, "./RData/04_output/04_plasmid_compiled_data.csv", row.names = F, na = "")
