###########################################
### Correlate ARGs + MGEs + plasmids + taxa
### mined from metagenomes
### Lou LaMartina, May 18, 2022
###########################################


setwd("~/Desktop/Lab/Projects/Resistance")

library(vegan)
library(reshape2)


#################
### load data ###
#################
# see 05_metagenome_geneCounts.R
# identities = 100%
# taxa bitscore > 99%-tile

# gene counts normalized to #16S genes in samples
arg_counts <- read.csv("./RData/05_output/05_ARGs_relabun_good.csv")
rownames(arg_counts) <- arg_counts$File_name
arg_counts <- arg_counts[-1]

int_counts <- read.csv("./RData/05_output/05_int_relabun_good.csv")
rownames(int_counts) <- int_counts$File_name
int_counts <- int_counts[-1]

plas_counts <- read.csv("./RData/05_output/05_plasmid_relabun_good.csv")
rownames(plas_counts) <- plas_counts$File_name
plas_counts <- plas_counts[-1]

tax_counts <- read.csv("./RData/05_output/05_taxa_relabun_good.csv")
rownames(tax_counts) <- tax_counts$File_name
tax_counts <- tax_counts[-1]


# gene info, from outputs of respective programs
arg_info <- read.csv("./RData/05_output/05_ARGs_geneInfo.csv")
int_info <- read.csv("./RData/05_output/05_int_geneInfo.csv")
plas_info <- read.csv("./RData/05_output/05_plasmid_geneInfo.csv")
tax_info <- read.csv("./RData/05_output/05_taxa_geneInfo.csv")


# sample info
info <- read.csv("./Data/Metadata/Metagenomes_sample_info.csv")




############################
### combine same samples ###
############################

# make list to loop through
counts.ls <- list(plas_counts, int_counts, arg_counts, tax_counts)
names(counts.ls) <- c("plas_counts", "int_counts", "arg_counts", "tax_counts")
datasets <- names(counts.ls)


# SandIsland__HI__2018_10_16 was split in filter + filtrate
for (i in datasets) {
  smps <- info$File_name[grep("SandIsland__HI__2018_10_16", info$Sample_name)]
  df <- counts.ls[[i]]
  print(nrow(df))
  S251_S2__S262_S13 <- colSums(df[rownames(df) %in% smps,])
  df <- rbind(df, S251_S2__S262_S13)
  rownames(df)[nrow(df)] <- "S251_S2__S262_S13"
  counts.ls[[i]] <- df
}


# SandIsland__HI__2018_10_16 was split in filter + filtrate
# AND two contigs files
for (i in datasets) {
  smps <- info$File_name[grep("SandIsland__HI__2018_11_28", info$Sample_name)]
  df <- counts.ls[[i]]
  print(nrow(df))
  S10B_S4__S253_S4__S263_S14 <- colSums(df[rownames(df) %in% smps,])
  df <- rbind(df, S10B_S4__S253_S4__S263_S14)
  rownames(df)[nrow(df)] <- "S10B_S4__S253_S4__S263_S14"
  counts.ls[[i]] <- df
}


# two contigs files
for (i in datasets) {
  smps <- info$File_name[info$FT_number == "25457"]
  df <- counts.ls[[i]]
  print(nrow(df))
  S257_S8__S199B_S8 <- colSums(df[rownames(df) %in% smps,])
  df <- rbind(df, S257_S8__S199B_S8)
  rownames(df)[nrow(df)] <- "S257_S8__S199B_S8"
  counts.ls[[i]] <- df
}


# keep just those samples
for (i in datasets) {
  dups <- c("S199B_S8", "S10B_S4", "S262_S13", "S251_S2", "S253_S4", "S263_S14", "S257_S8")
  df <- counts.ls[[i]]
  print(nrow(df))
  df <- df[! rownames(df) %in% dups,]
  rownames(df)[rownames(df) %in%
                 c("S251_S2__S262_S13", "S10B_S4__S253_S4__S263_S14", "S257_S8__S199B_S8")] <-
    c("S251_S2", "S253_S4", "S257_S8")
  df$File_name <- rownames(df)
  df <- merge(info[c("File_name", "Sample_name")], df, by = "File_name")
  rownames(df) <- df$Sample_name
  df <- df[-c(1:2)]
  counts.ls[[i]] <- df
}


# extract
plas_counts <- counts.ls[["plas_counts"]]; print(nrow(plas_counts))
int_counts <- counts.ls[["int_counts"]]; print(nrow(int_counts))
arg_counts <- counts.ls[["arg_counts"]]; print(nrow(arg_counts))
tax_counts <- counts.ls[["tax_counts"]]; print(nrow(tax_counts))


# cleanup
rm(df, smps, i, dups, datasets, counts.ls,
   S10B_S4__S253_S4__S263_S14, S251_S2__S262_S13, S257_S8__S199B_S8)


### all possible combinations?
combos <- paste0(t(combn(c("TAX", "INT", "PLAS", "ARG"), 2))[,1], 
                 "_", t(combn(c("TAX", "INT", "PLAS", "ARG"), 2))[,2])
# "TAX_INT" "TAX_PLAS"  "TAX_ARG" "INT_PLAS"  "INT_ARG" "PLAS_ARG" 




#################
### TAX ~ INT ###
#################

# combine & convert to binary data
tax_int <- cbind(tax_counts, int_counts)
tax_int[tax_int > 0] <- 1


# jaccard distances
tax_int.cor <- vegdist(t(tax_int), method = "jaccard")
tax_int <- data.frame(scores(tax_int.cor))
tax_int$Gene_ID <- rownames(tax_int)
tax_int <- melt(tax_int, variable.name = "INT", value.name = "Dist")


# subset int ~ tax
tax_int <- tax_int[tax_int$Gene_ID %in% colnames(tax_counts),]
tax_int <- tax_int[tax_int$INT %in% colnames(int_counts),]


# focus on 100% similarity
tax_int <- tax_int[tax_int$Dist == 1,]


# subset by phylum
tax_int.ls <- list()
for(i in unique(tax_info$Phylum)) {
  tax_int.ls[[i]] <- tax_info$Gene_ID[tax_info$Phylum == i]
  tax_int.ls[[i]] <- tax_int[tax_int$Gene_ID %in% tax_int.ls[[i]],]
  if (nrow(tax_int.ls[[i]]) == 0) {
    tax_int.ls[[i]] <- NA
  }
}
tax_int.ls <- tax_int.ls[is.na(tax_int.ls) == F]


# add tax info
for(i in names(tax_int.ls)) {
  temp <- tax_info[tax_info$Phylum == i,  c(1,6:13)]
  temp <- unique(temp[temp$Gene_ID %in% tax_int.ls[[i]]$Gene_ID,])
  tax_int.ls[[i]] <- merge(tax_int.ls[[i]], temp, by = "Gene_ID")
  rm(temp)
}


# add int info
for(i in names(tax_int.ls)) {
  colnames(tax_int.ls[[i]])[1:2] <- c("TAX", "Gene_ID")
  temp <- int_info[c(1,6,7)]
  temp <- unique(temp[temp$Gene_ID %in% tax_int.ls[[i]]$Gene_ID,])
  tax_int.ls[[i]] <- merge(tax_int.ls[[i]], temp, by = "Gene_ID")
  colnames(tax_int.ls[[i]])[1:2] <- c("INT", "TAX")
  rm(temp)
}


# save - now & not the end, these are big files
path = "./RData/06_output/TAX_INT/"
dir.create(path)
for(i in names(tax_int.ls)) {
  write.csv(tax_int.ls[[i]], paste0(path, "TAX_INT_", i, ".csv"), row.names = F, na = "")
}


# cleanup
rm(tax_int.cor, tax_int.ls, tax_int)




##################
### TAX ~ PLAS ###
##################

# combine & convert to binary data
tax_plas <- cbind(tax_counts, plas_counts)
tax_plas[tax_plas > 0] <- 1


# jaccard distances
tax_plas.cor <- vegdist(t(tax_plas), method = "jaccard")
tax_plas <- data.frame(scores(tax_plas.cor))
tax_plas$Gene_ID <- rownames(tax_plas)
tax_plas <- melt(tax_plas, variable.name = "PLAS", value.name = "Dist")


# subset pl ~ tax
tax_plas <- tax_plas[tax_plas$Gene_ID %in% colnames(tax_counts),]
tax_plas <- tax_plas[tax_plas$PLAS %in% colnames(plas_counts),]


# focus on 100% similarity
tax_plas <- tax_plas[tax_plas$Dist == 1,]


# subset phylum
tax_plas.ls <- list()
for(i in unique(tax_info$Phylum)) {
  tax_plas.ls[[i]] <- tax_info$Gene_ID[tax_info$Phylum == i]
  tax_plas.ls[[i]] <- tax_plas[tax_plas$Gene_ID %in% tax_plas.ls[[i]],]
  if (nrow(tax_plas.ls[[i]]) == 0) {
    tax_plas.ls[[i]] <- NA
  }
}
tax_plas.ls <- tax_plas.ls[is.na(tax_plas.ls) == F]


# add tax info
for(i in names(tax_plas.ls)) {
  temp <- tax_info[tax_info$Phylum == i,  c(1,6:13)]
  temp <- unique(temp[temp$Gene_ID %in% tax_plas.ls[[i]]$Gene_ID,])
  tax_plas.ls[[i]] <- merge(tax_plas.ls[[i]], temp, by = "Gene_ID")
  rm(temp)
}


# add plasmid info
for(i in names(tax_plas.ls)) {
  colnames(tax_plas.ls[[i]])[1:2] <- c("TAX", "Gene_ID")
  temp <- plas_info[c(1,6,7)]
  temp <- unique(temp[temp$Gene_ID %in% tax_plas.ls[[i]]$Gene_ID,])
  tax_plas.ls[[i]] <- merge(tax_plas.ls[[i]], temp, by = "Gene_ID")
  colnames(tax_plas.ls[[i]])[1:2] <- c("PLAS", "TAX")
  rm(temp)
}


# save - now & not the end, these are big files
path = "./RData/06_output/TAX_PLAS/"
dir.create(path)
for(i in names(tax_plas.ls)) {
  write.csv(tax_plas.ls[[i]], paste0(path, "TAX_PLAS_", i, ".csv"), row.names = F, na = "")
}


# cleanup
rm(tax_plas.cor, tax_plas.ls, tax_plas)




#################
### TAX ~ ARG ###
#################

# combine & convert to binary data
tax_arg <- cbind(tax_counts, arg_counts)
tax_arg[tax_arg > 0] <- 1


# jaccard distances
tax_arg.cor <- vegdist(t(tax_arg), method = "jaccard")
tax_arg <- data.frame(scores(tax_arg.cor))
tax_arg$Gene_ID <- rownames(tax_arg)
tax_arg <- melt(tax_arg, variable.name = "ARG", value.name = "Dist")


# subset arg ~ tax
tax_arg <- tax_arg[tax_arg$Gene_ID %in% colnames(tax_counts),]
tax_arg <- tax_arg[tax_arg$ARG %in% colnames(arg_counts),]


# focus on 100% similarity
tax_arg <- tax_arg[tax_arg$Dist == 1,]


# subset phylum
tax_arg.ls <- list()
for(i in unique(tax_info$Phylum)) {
  tax_arg.ls[[i]] <- tax_info$Gene_ID[tax_info$Phylum == i]
  tax_arg.ls[[i]] <- tax_arg[tax_arg$Gene_ID %in% tax_arg.ls[[i]],]
  if (nrow(tax_arg.ls[[i]]) == 0) {
    tax_arg.ls[[i]] <- NA
  }
}
tax_arg.ls <- tax_arg.ls[is.na(tax_arg.ls) == F]


# add tax info
for(i in names(tax_arg.ls)) {
  temp <- tax_info[tax_info$Phylum == i,  c(1,6:13)]
  temp <- unique(temp[temp$Gene_ID %in% tax_arg.ls[[i]]$Gene_ID,])
  tax_arg.ls[[i]] <- merge(tax_arg.ls[[i]], temp, by = "Gene_ID")
  rm(temp)
}


# add arg info
for(i in names(tax_arg.ls)) {
  colnames(tax_arg.ls[[i]])[1:2] <- c("TAX", "Gene_ID")
  temp <- arg_info[c(1,6,7)]
  temp <- unique(temp[temp$Gene_ID %in% tax_arg.ls[[i]]$Gene_ID,])
  tax_arg.ls[[i]] <- merge(tax_arg.ls[[i]], temp, by = "Gene_ID")
  colnames(tax_arg.ls[[i]])[1:2] <- c("ARG", "TAX")
  rm(temp)
}


# save - now & not the end, these are big files
path = "./RData/06_output/TAX_ARG/"
dir.create(path)
for(i in names(tax_arg.ls)) {
  write.csv(tax_arg.ls[[i]], paste0(path, "TAX_ARG_", i, ".csv"), row.names = F, na = "")
}


# cleanup
rm(tax_arg.cor, tax_arg.ls, tax_arg)




#################
### INT ~ ARG ###
#################

# combine & convert to binary data
int_arg <- cbind(arg_counts, int_counts)
int_arg[int_arg > 0] <- 1


# jaccard distances
int_arg.cor <- vegdist(t(int_arg), method = "jaccard")
int_arg <- data.frame(scores(int_arg.cor))
int_arg$Gene_ID <- rownames(int_arg)
int_arg <- melt(int_arg, variable.name = "INT", value.name = "Dist")


# subset int ~ arg
int_arg <- int_arg[int_arg$Gene_ID %in% colnames(arg_counts),]
int_arg <- int_arg[int_arg$INT %in% colnames(int_counts),]


# focus on 100% similarity
int_arg <- int_arg[int_arg$Dist == 1,]


# new category names
arg_info$altARG <- gsub(" antibiotic", "", arg_info$Category)
arg_info$altARG <- gsub("; ", "_", arg_info$altARG)
arg_info$altARG <- gsub("-", "", gsub(" ", "", arg_info$altARG))


# subset by arg category
int_arg.ls <- list()
for(i in unique(arg_info$altARG)) {
  int_arg.ls[[i]] <- int_arg[arg_info$altARG == i,]
}


# add arg info
for(i in names(int_arg.ls)) {
  temp <- arg_info[c(1,6,7)]
  temp <- unique(temp[temp$Gene_ID %in% int_arg.ls[[i]]$Gene_ID,])
  int_arg.ls[[i]] <- merge(int_arg.ls[[i]], temp, by = "Gene_ID")
  rm(temp)
}



# add int info
for(i in names(int_arg.ls)) {
  colnames(int_arg.ls[[i]]) <- c("ARG", "Gene_ID", "Dist", "ARG_name", "ARG_category")
  temp <- int_info[c(1,6,7)]
  temp <- unique(temp[temp$Gene_ID %in% int_arg.ls[[i]]$Gene_ID,])
  int_arg.ls[[i]] <- merge(int_arg.ls[[i]], temp, by = "Gene_ID")
  colnames(int_arg.ls[[i]])[c(1,6,7)] <- c("INT","INT_name", "INT_category")
}


# save - now & not the end, these are big files
path = "./RData/06_output/INT_ARG/"
dir.create(path)
for(i in names(int_arg.ls)) {
  write.csv(int_arg.ls[[i]], paste0(path, "INT_ARG_", i, ".csv"), row.names = F, na = "")
}


# cleanup
rm(int_arg.ls, int_arg, int_arg.cor)



##################
### INT ~ PLAS ###
##################

# combine & convert to binary data
int_plas <- cbind(int_counts, plas_counts)
int_plas[int_plas > 0] <- 1


# jaccard distances
int_plas.cor <- vegdist(t(int_plas), method = "jaccard")
int_plas <- data.frame(scores(int_plas.cor))
int_plas$Gene_ID <- rownames(int_plas)
int_plas <- melt(int_plas, variable.name = "PLAS", value.name = "Dist")


# subset plasmid ~ int
int_plas <- int_plas[int_plas$Gene_ID %in% colnames(int_counts),]
int_plas <- int_plas[int_plas$PLAS %in% colnames(plas_counts),]


# focus on 100% similarity
int_plas <- int_plas[int_plas$Dist == 1,]


# add INT info
temp <- unique(int_info[int_info$Gene_ID %in% int_plas$Gene_ID, c(1,6,7)])
int_plas <- merge(int_plas, temp, by = "Gene_ID")
colnames(int_plas) <- c("INT", "Gene_ID", "Dist", "INT_name", "INT_category")


# add PLAS info
temp <- unique(plas_info[plas_info$Gene_ID %in% int_plas$Gene_ID, c(1,6,7)])
int_plas <- merge(int_plas, temp, by = "Gene_ID")
colnames(int_plas)[c(1,6,7)] <- c("PLAS", "PLAS_name", "PLAS_category")


# save
write.csv(int_plas, "./RData/06_output/INT_PLAS.csv", row.names = F, na = "")

# cleanup
rm(int_plas.cor, int_plas, temp)




##################
### PLAS ~ ARG ###
##################

# combine & convert to binary data
plas_arg <- cbind(arg_counts, plas_counts)
plas_arg[plas_arg > 0] <- 1


# jaccard distances
plas_arg.cor <- vegdist(t(plas_arg), method = "jaccard")
plas_arg <- data.frame(scores(plas_arg.cor))
plas_arg$Gene_ID <- rownames(plas_arg)
plas_arg <- melt(plas_arg, variable.name = "PLAS", value.name = "Dist")


# subset plasmid ~ arg
plas_arg <- plas_arg[plas_arg$Gene_ID %in% colnames(arg_counts),]
plas_arg <- plas_arg[plas_arg$PLAS %in% colnames(plas_counts),]


# focus on 100% similarity
plas_arg <- plas_arg[plas_arg$Dist == 1,]


# add ARG info
temp <- unique(arg_info[arg_info$Gene_ID %in% plas_arg$Gene_ID, c(1,6,7)])
plas_arg <- merge(plas_arg, temp, by = "Gene_ID")
colnames(plas_arg) <- c("ARG", "Gene_ID", "Dist", "ARG_name", "ARG_category")


# add PLAS info
temp <- unique(plas_info[plas_info$Gene_ID %in% plas_arg$Gene_ID, c(1,6,7)])
plas_arg <- merge(plas_arg, temp, by = "Gene_ID")
colnames(plas_arg)[c(1,6,7)] <- c("PLAS", "PLAS_name", "PLAS_category")


# save
write.csv(plas_arg, "./RData/06_output/PLAS_ARG.csv", row.names = F, na = "")


# cleanup
rm(plas_arg.cor, plas_arg, temp)

