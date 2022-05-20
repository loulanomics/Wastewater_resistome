#######################################
### Explain ARG abundance with ASVs and
### environmental metadata
### Lou LaMartina, May 17, 2022
#######################################


setwd("~/Desktop/Lab/Projects/Resistance")

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(scales)

load("./RData/03_env.RData")


#################
### load data ###
#################

# # v4v5 counts
# month_asvs <- read.csv("./Data/Communities/MonthSeries_V4V5_counts.csv")
# rownames(month_asvs) <- month_asvs$Sample_name
# month_asvs <- month_asvs[-1]
# 
# cities_asvs <- read.csv("./Data/Communities/Cities_V4V5_counts.csv")
# rownames(cities_asvs) <- cities_asvs$Sample_name
# cities_asvs <- cities_asvs[-1]
# 
# neigh_asvs <- read.csv("./Data/Communities/Neighborhoods_V4V5_counts.csv")
# rownames(neigh_asvs) <- neigh_asvs$Sample_name
# neigh_asvs <- neigh_asvs[-1]
# 
# 
# # combine
# counts_v4v5 <- rbind(month_asvs, cities_asvs, neigh_asvs)
# rm(month_asvs, cities_asvs, neigh_asvs)
# 
# 
# # v4 counts
# counts_v4 <- read.csv("./Data/Communities/WeekSeries_V4_counts.csv")
# rownames(counts_v4) <- counts_v4$Sample_name
# counts_v4 <- counts_v4[-1]
# 
# 
# # convert to relative abundances
# relabun_v4v5 <- counts_v4v5 / rowSums(counts_v4v5)
# relabun_v4 <- counts_v4 / rowSums(counts_v4)
# 
# 
# # sample info
# time_info <- read.csv("./Data/Metadata/TimeSeries_sample_info.csv")
# cities_info <- read.csv("./Data/Metadata/Cities_sample_info.csv")
# neigh_info <- read.csv("./Data/Metadata/Neighborhoods_sample_info.csv")
# 
# 
# # taxonomy
# taxa_v4v5 <- read.csv("./Data/Communities/All_V4V5_taxonomy.csv")
# taxa_v4 <- read.csv("./Data/Communities/WeekSeries_V4_taxonomy.csv")
# 
# 
# # human asvs from lamartina 2021
# human_asvs <- read.csv("./Data/Communities/V4V5_bodysite_associations.csv")
# 
# 
# # ddPCR
# args <- read.csv("./RData/02_ddPCR_conversions_uniq.csv")
# 
# 
# # create arg abundance matrix
# args.d <- dcast(Sample_name ~ Target, value.var = "copiespermL", data = args)
# rownames(args.d) <- args.d$Sample_name
# args.d <- args.d[-1]
# 
#
# save env
# save.image("./RData/03_env.RData")





#################
### ARG ~ ASV ###
#################

# subset to samples that we have ddpcr and seq data for
arg_relabun <- relabun_v4v5[rownames(relabun_v4v5) %in% rownames(args.d),]
arg_relabun <- arg_relabun[colSums(arg_relabun) > 0]
arg_asv <- args.d[rownames(args.d) %in% rownames(arg_relabun),]


# get taxa that occur in at least 2 samples (10404)
mintax <- c()
for (i in colnames(arg_relabun)) {
  len <- length(which(arg_relabun[[i]] > 0))
  if (len > 1) {
  mintax[i] <- i
  }
}
arg_relabun <- arg_relabun[colnames(arg_relabun) %in% mintax]
arg_relabun <- arg_relabun[order(rownames(arg_relabun)),]
identical(rownames(arg_asv), rownames(arg_relabun))
arg_asv <- cbind(arg_asv, arg_relabun)
arg_asv[is.na(arg_asv)] <- 0


# spearman correlation
arg_asv.cor <- melt(cor(arg_asv, method = "spearman"))


# subset arg ~ asv
arg_asv.cor$Var1 <- as.character(arg_asv.cor$Var1)
arg_asv.cor$Var2 <- as.character(arg_asv.cor$Var2)
arg_asv.cor <- arg_asv.cor[arg_asv.cor$Var1 %in% colnames(args.d),]
arg_asv.cor <- arg_asv.cor[arg_asv.cor$Var2 %in% colnames(arg_relabun),]
colnames(arg_asv.cor) <- c("ARG", "ASV", "rho")


# get minimum correlations
length(unique(arg_asv.cor[arg_asv.cor$rho < -0.5 | arg_asv.cor$rho > 0.5, "ASV"]))
arg_asv.min <- arg_asv.cor[arg_asv.cor$rho < -0.1 | arg_asv.cor$rho > 0.1,]
arg_asv.min <- merge(arg_asv.min, taxa_v4v5, by = "ASV")


# human?
arg_asv.min <- merge(human_asvs[c(1,6)], arg_asv.min, by = "ASV", all.y = T)
arg_asv.min$Source[is.na(arg_asv.min$Source)] <- "Nonhuman"
#write.csv(arg_asv.df, "./RData/03_ARG_v4v5_corr.csv", row.names = F, na = "")
arg_asv.min$Source[arg_asv.min$Source != "Nonhuman"] <- "Human"


# add phylum
arg_asv.min$xaxis <- paste(arg_asv.min$Phylum, arg_asv.min$ASV)
arg_asv.min <- arg_asv.min[arg_asv.min$ARG == "intI1",]


# ok to combine same names? yes
cbind(aggregate(rho ~ Name, sd, data = arg_asv.min),
      aggregate(rho ~ Name, mean, data = arg_asv.min))


heat.plot <-
  ggplot(arg_asv.min, aes(x = ARG, y = Name, fill = rho)) +
  geom_tile() +
  facet_grid(Phylum ~ ., scales = "free", space = "free", switch = "y") +
  scale_fill_gradientn(colors = c("#A50026", "#F46D43", "white", "#74ADD1", "#313695"),
                       values = rescale(c(min(arg_asv.min$rho), -0.5, 0, 0.5, max(arg_asv.min$rho))),
                       breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits = c(-0.6,0.6)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, color = "black", face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        legend.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8, color = "black", face = "bold"),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black", size = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_colorbar(barheight = 0.5, barwidth = 8, title.position = "top", title.hjust = 0.5)) +
  labs(y = "Phylum of significant ASVs", fill = "Spearman rho score\nASV relative abundance\nto intI1 concentration")
heat.plot

#ggsave("./Plots/heatmap.pdf", plot = heat.plot, device = "pdf", width = 5, height = 10, units = "in")


# top rho ASVs
top <- sort(c(arg_asv.min$ASV[order(arg_asv.min$rho, decreasing = T)][1:3], # "V4V500014" "V4V500066" "V4V500071" hi rho
         arg_asv.min$ASV[order(arg_asv.min$rho, decreasing = F)][1:3])) # "V4V500244" "V4V500242" "V4V500325" lo rho

arg_asv_relabun <- arg_relabun[colnames(arg_relabun) %in% top]
arg_asv_relabun <- data.frame(t(arg_asv_relabun))
arg_asv_relabun <- arg_asv_relabun / rowSums(arg_asv_relabun)
arg_asv_relabun <- data.frame(t(arg_asv_relabun))

arg_asv_relabun <- data.frame(Sample_name = rownames(arg_asv_relabun), arg_asv_relabun)
arg_asv_relabun <- merge(data.frame(Sample_name = rownames(args.d), intI1 = args.d$intI1),
                         arg_asv_relabun, by = "Sample_name")

arg_asv_relabun.m <- melt(arg_asv_relabun, id.vars = c("Sample_name", "intI1"), 
                          variable.name = "ASV", value.name = "Relabun")
#arg_asv_relabun.m <- merge(arg_asv_relabun.m, taxa_v4v5[c("Phylum", "Genus", "ASV")], by = "ASV")
arg_asv_relabun.m$Label <- paste0(arg_asv_relabun.m$Phylum, " (", arg_asv_relabun.m$Genus, ")")

names(top) <- paste0(taxa_v4v5$Phylum[taxa_v4v5$ASV %in% top], " (", taxa_v4v5$Genus[taxa_v4v5$ASV %in% top], ")")


reg.plot <-
  ggplot(arg_asv_relabun.m, aes(x = log10(intI1), y = Relabun, color = ASV)) +
  geom_point(size = 1) +
  geom_smooth(se = F) +
  scale_color_manual(values = c(brewer.pal(11, "RdYlBu")[11:9], brewer.pal(11, "RdYlBu")[1:3]),
                     labels = names(top)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8, face = "italic"),
        strip.text = element_blank(),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(size = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(x = "Log10 intI1 concentration (copies per mL)", color = "ASV phylum (genus)",
       y = "ASV relative abundances\n(sample sum = 1, ASV sum = 1)")
reg.plot

ggsave("./Plots/regress lines.pdf", plot = reg.plot, device = "pdf", width = 8, height = 4, units = "in")




###########
### fit ###
###########

# https://www.frontiersin.org/articles/10.3389/fmicb.2021.679805/full
# Pearson
# between environmental factors (TN, TP, Cl–, Ca2+, and Mg2+) 
# with the relative abundance of bacterial community and ARGs. 
#
# Non-metric multidimensional scaling (NMDS) 
# difference of bacteria between different sampling sites
#
# Redundancy analysis (RDA) was employed to assess the effects of 
# environmental factors and ARGs on the bacterial community. 
#
# The co-occurrence between abundance of ARGs and bacterial taxa 
# was analyzed using network analysis based on the Pearson correlation
#
# RDA, Mantel test, Pearson correlation, and heatmap were performed in RStudio







