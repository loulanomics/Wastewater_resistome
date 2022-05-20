#######################################
### organize metagenome program outputs
### into gene counts per sample
### Lou LaMartina, May 18, 2022
########################################


setwd("~/Desktop/Lab/Projects/Resistance")

library(reshape2)


############
### load ###
############
# see 04_metagenome_infoPrep.R

# sample info
info <- read.csv("./Data/Metadata/Metagenomes_sample_info.csv")


# RGI results - mined ARGs
arg_info <- read.csv("./RData/04_output/04_RGI_compiled_data.csv")
arg_info[arg_info == ""] <- NA


# mobileOG-db results - mined integrons
int_info <- read.csv("./RData/04_output/04_integron_compiled_data.csv")
int_info[int_info == ""] <- NA
int_info <- int_info[! int_info$Gene_Name %in% c("", "NA:Keyword"),]


# kaiju results - taxonomy from mined 16s genes
tax_info <- read.csv("./RData/04_output/04_kaiju_compiled_data.csv")
tax_info[tax_info == ""] <- NA


# plasmid finder results - mined plasmids
plas_info <- read.csv("./RData/04_output/04_plasmid_compiled_data.csv")
plas_info[plas_info == ""] <- NA




###################
### unique IDs ####
###################

# args - use ARO ID as unique identifier
arg_info$Gene_ID <- paste0("ARG_", arg_info$ARO)


# mges - use mobileOG ID as unique identifier
int_info$Gene_ID <- gsub("mobileOG_", "INT_", int_info$mobileOG_ID)


# taxa - NCBI ID as unique identifier
temp <- data.frame(NCBI_ID = unique(tax_info$NCBI_ID))
temp$Gene_ID <- paste0("TAX_", temp$NCBI_ID)
tax_info <- merge(temp, tax_info, by = "NCBI_ID")
rm(temp)


# plasmids - make own
temp <- unique(plas_info[c("sseqid", "Category")])
for (i in unique(temp$Category)) {
  j = temp$sseqid[temp$Category == i]
  k = length(j)
  temp$Gene_ID[temp$Category == i] <- paste0("PLAS_", i, sprintf("%02d", 1:length(j)))
}
plas_info <- merge(temp[-2], plas_info, by = "sseqid")




##################
### all counts ###
##################

# ARGs
arg_counts <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = arg_info)
rownames(arg_counts) <- arg_counts$File_name
arg_counts <- arg_counts[-1]
ncol(arg_counts) # 1029


# MGEs
int_counts <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = int_info)
rownames(int_counts) <- int_counts$File_name
int_counts <- int_counts[-1]
ncol(int_counts) # 16907


# plasmids
plas_counts <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = plas_info)
rownames(plas_counts) <- plas_counts$File_name
plas_counts <- plas_counts[-1]
ncol(plas_counts) # 121


# taxa
tax_counts <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID", length, data = tax_info)
rownames(tax_counts) <- tax_counts$File_name
tax_counts <- tax_counts[-1]
ncol(tax_counts) # 24341


# 16s genes per sample
tax_sums <- data.frame(File_name = rownames(tax_counts), Sum = rowSums(tax_counts))
mean(tax_sums$Sum) # 100018.8




#######################
### good alignments ###
#######################

# ARGs
arg_good <- arg_info[arg_info$Best_Identities == 100,]
arg_good <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = arg_good)
rownames(arg_good) <- arg_good$File_name
arg_good <- arg_good[-1]
ncol(arg_good) # 942


# MGEs - 
int_good <- int_info[int_info$Pident == 100,]
int_good <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = int_good)
rownames(int_good) <- int_good$File_name
int_good <- int_good[-1]
ncol(int_good) # 5503


# plasmids (already 100%)
plas_good <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID",  length, data = plas_info)
rownames(plas_good) <- plas_good$File_name
plas_good <- plas_good[-1]
ncol(plas_good) # 121


# taxa
tax_good <- tax_info[tax_info$Score > quantile(tax_info$Score, 0.99),]
tax_good <- dcast(File_name ~ Gene_ID, value.var = "Gene_ID", length, data = tax_good)
rownames(tax_good) <- tax_good$File_name
tax_good <- tax_good[-1]
ncol(tax_good) # 946


# 16s genes per sample
tax_good_sums <- data.frame(File_name = rownames(tax_good), Sum = rowSums(tax_good))
mean(tax_good_sums$Sum) # 999.0556




#################
### normalize ###
#################

#######
### all

# args
relabuns.ls <- list()
for (i in tax_sums$File_name) {
  smp <- arg_counts[rownames(arg_counts) == i,]
  relabuns.ls[[i]] <- smp / tax_sums$Sum[tax_sums$File_name == i]
}
arg_relabun <- do.call(rbind, relabuns.ls)


# ints
relabuns.ls <- list()
for (i in tax_sums$File_name) {
  smp <- int_counts[rownames(int_counts) == i,]
  relabuns.ls[[i]] <- smp / tax_sums$Sum[tax_sums$File_name == i]
}
int_relabun <- do.call(rbind, relabuns.ls)


# plasmids
relabuns.ls <- list()
for (i in rownames(plas_counts)) {
  smp <- plas_counts[rownames(plas_counts) == i,]
  relabuns.ls[[i]] <- smp / tax_sums$Sum[tax_sums$File_name == i]
}
plas_relabun <- do.call(rbind, relabuns.ls)
rm(relabuns.ls, smp)


# relative abundance of taxa
tax_relabun <- tax_counts / rowSums(tax_counts)



########
### good

# args
relabuns.ls <- list()
for (i in tax_good_sums$File_name) {
  smp <- arg_good[rownames(arg_good) == i,]
  relabuns.ls[[i]] <- smp / tax_good_sums$Sum[tax_good_sums$File_name == i]
}
arg_good_relabun <- do.call(rbind, relabuns.ls)


# ints
relabuns.ls <- list()
for (i in tax_good_sums$File_name) {
  smp <- int_good[rownames(int_good) == i,]
  relabuns.ls[[i]] <- smp / tax_good_sums$Sum[tax_good_sums$File_name == i]
}
int_good_relabun <- do.call(rbind, relabuns.ls)


# plasmids
relabuns.ls <- list()
for (i in tax_good_sums$File_name) {
  smp <- plas_good[rownames(plas_good) == i,]
  relabuns.ls[[i]] <- smp / tax_good_sums$Sum[tax_good_sums$File_name == i]
}
plas_good_relabun <- do.call(rbind, relabuns.ls)
rm(relabuns.ls, smp)


# relative abundance of taxa
tax_good_relabun <- tax_good / rowSums(tax_good)




#####################
### simplify data ###
#####################

# args
arg_sub <- arg_info[c("Gene_ID", "File_name", "Sample_name", "Contig", 
                          "Contig_ORF", "Best_Hit_ARO", "Drug_Class", "Best_Identities")]
colnames(arg_sub) <- c("Gene_ID", "File_name", "Sample_name", "Contig", 
                           "ORF", "Gene_name", "Category", "Identity")
arg_sub$ORF <- as.numeric(sapply(strsplit(arg_sub$ORF, "_"), '[', 4))


# mges
int_sub <- int_info[c("Gene_ID", "File_name", "Sample_name", "Specific_Contig", 
                          "Contig_ORF_Name", "Gene_Name", "Minor_mobileOG_Category", "Pident")]
colnames(int_sub) <- colnames(arg_sub)
int_sub$ORF <- as.numeric(sapply(strsplit(int_sub$ORF, "_"), '[', 4))


# plasmids (takes a min)
plasORFs.ls <- list()
for (i in unique(plas_info$qseqid)) {
  temp <- plas_info[plas_info$qseqid == i,]
  plasORFs.ls[[i]] <- paste0(i, "__", 1:length(temp$qseqid))
}

plasORFs <- data.frame(Contig_ORF = unlist(plasORFs.ls))
plasORFs$qseqid <- sapply(strsplit(plasORFs$Contig_ORF, "__"), '[', 1)

plas_sub <- merge(plasORFs, plas_info, by = "qseqid")
plas_sub <- plas_sub[c("Gene_ID", "File_name", "Sample_name", "qseqid", 
                        "Contig_ORF", "sseqid", "Category", "pident")]
colnames(plas_sub) <- colnames(arg_sub)
plas_sub$ORF <- as.numeric(sapply(strsplit(plas_sub$ORF, "__"), '[', 2))
rm(temp, plasORFs, plasORFs.ls, i)


# taxa 
tax_info <- merge(tax_info, info[1:2], by = "File_name")
tax_sub <- tax_info[c("Gene_ID", "File_name", "Sample_name", "Contig", "Score",
                       "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subgroup")]
colnames(tax_sub)[5] <- ("Bitscore")


# combine contigs and ORFs where necessary - not needed for taxa
arg_sub$Contig_ORF <- paste0(arg_sub$Contig, "__", arg_sub$ORF)
int_sub$Contig_ORF <- paste0(int_sub$Contig, "__", int_sub$ORF)
plas_sub$Contig_ORF <- paste0(plas_sub$Contig, "__", plas_sub$ORF)




############
### save ###
############

int_relabun <- data.frame(File_name = rownames(int_relabun), int_relabun)
arg_relabun <- data.frame(File_name = rownames(arg_relabun), arg_relabun)
plas_relabun <- data.frame(File_name = rownames(plas_relabun), plas_relabun)
tax_relabun <- data.frame(File_name = rownames(tax_relabun), tax_relabun)

int_good_relabun <- data.frame(File_name = rownames(int_good_relabun), int_good_relabun)
arg_good_relabun <- data.frame(File_name = rownames(arg_good_relabun), arg_good_relabun)
plas_good_relabun <- data.frame(File_name = rownames(plas_good_relabun), plas_good_relabun)
tax_good_relabun <- data.frame(File_name = rownames(tax_good_relabun), tax_good_relabun)

write.csv(int_relabun, "./RData/05_output/05_int_relabun_all.csv", row.names = F, na = "")
write.csv(arg_relabun, "./RData/05_output/05_ARGs_relabun_all.csv", row.names = F, na = "")
write.csv(plas_relabun, "./RData/05_output/05_plasmid_relabun_all.csv", row.names = F, na = "")
write.csv(tax_relabun, "./RData/05_output/05_taxa_relabun_all.csv", row.names = F, na = "")

write.csv(int_good_relabun, "./RData/05_output/05_int_relabun_good.csv", row.names = F, na = "")
write.csv(arg_good_relabun, "./RData/05_output/05_ARGs_relabun_good.csv", row.names = F, na = "")
write.csv(plas_good_relabun, "./RData/05_output/05_plasmid_relabun_good.csv", row.names = F, na = "")
write.csv(tax_good_relabun, "./RData/05_output/05_taxa_relabun_good.csv", row.names = F, na = "")

write.csv(int_sub, "./RData/05_output/05_int_geneInfo.csv", row.names = F, na = "")
write.csv(arg_sub, "./RData/05_output/05_ARGs_geneInfo.csv", row.names = F, na = "")
write.csv(plas_sub, "./RData/05_output/05_plasmid_geneInfo.csv", row.names = F, na = "")
write.csv(tax_sub, "./RData/05_output/05_taxa_geneInfo.csv", row.names = F, na = "")

save.image("./RData/05_output/05_env.RData")
