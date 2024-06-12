library (dplyr)
library(RColorBrewer)
library(reshape)
library(ggpubr)
library(pheatmap)
library(tidyverse)

#Read in the meta data and modify it for tximport
stats <- read.csv("~/Google Drive/My Drive/vsg/vsgseq2/data/meta.csv")
orthogroups <- read.csv("~/pkgs/analyse/orthofinder/protein/all_samples/OrthoFinder/Results_May31/Orthogroups/Orthogroups.csv")
orthogroups_long <- orthogroups %>% separate_rows(Name, sep = ",")
my.files <-  list.files(list.dirs(path = "/Users/goldriev/pkgs/analyse/all_data_results", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)

tpm <- lapply(Sys.glob(my.files), read.table, header = TRUE)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

names <- sapply(strsplit(my.files, split="\\/"), function(x)x[7])
names <- gsub("_quant", "", names)
names(tpm) <- names

tpm_per <- function(DF) {
  DF$Percent <- (DF$TPM /(sum(DF$TPM)))*100
  DF <- DF %>% select(Name, Percent)
  return(DF)
}

tpm2 <- lapply(tpm, tpm_per)
wide_tpm <- dplyr:::reduce(tpm2, full_join, by = "Name")
names <- append(names, 'Name', after = 0)
colnames(wide_tpm) <- names

tpm_select <- function(DF) {
  DF <- DF %>% select(Name, TPM)
  return(DF)
}

tpm_selected <- lapply(tpm, tpm_select)
wide_tpm_selected <- dplyr:::reduce(tpm_selected, full_join, by = "Name")
colnames(wide_tpm_selected) <- names

#write.csv(wide_tpm_selected, '~/Desktop/wide_tpm.csv', row.names=FALSE)

keep_tpm <- wide_tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 10, Name))

keep_tpm <- keep_tpm %>% mutate(Name = if_else(Names == "TRUE", Name, "zVSGome"))

keep_tpm <- subset(keep_tpm, select = -c(Names) )

long <- melt(keep_tpm, id.vars = c("Name"), variable.name = "isolate", value.name = "P")
long <- merge(x=stats, y=long, by.x="isolate", by.y="variable")[]
long <- merge(x = long, y = orthogroups_long, by.x = "Name", by.y = "Name", all.x = TRUE)
long <- long %>% replace_na(list(col1 = "Other", col2 = "Other"))

png("~/Desktop/orthogroup.png", units="in", width=15, height=15, res=300)
ggboxplot(mouse_comb, x = "type", y = "value",
          add.params = list(size = 2), color = "stage",
          add = "jitter") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap("Orthogroup", scales = 'free_x') 
dev.off()

cow_comb <- long[long$host == 'cow', ]
cow_comb <- cow_comb[cow_comb$bleed == 'small', ]
mouse_comb <- long[long$host == 'mouse', ]
mouse_comb <- mouse_comb[mouse_comb$type == 'WT', ]
cow_comb$stage <- as.numeric(cow_comb$stage)

png("~/Desktop/wt_mouse.png", units="in", width=1, height=6, res=300)
ggbarplot(long, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(cow_comb$Name)))
coul <- append(coul, "black")
cow_comb$stage <- as.numeric(cow_comb$stage)

png("~/Desktop/cow_small.png", units="in", width=20, height=10, res=300)
ggbarplot(cow_comb, x = "stage", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(~type, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

png("~/Desktop/pacbio_combined_mouse.png", units="in", width=12, height=12, res=300)
ggbarplot(mouse_comb, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage + tech, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

#Read in the meta data and modify it for tximport

my.files <-  list.files(list.dirs(path = "/Users/goldriev/pkgs/analyse/full_cds/mouse", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)
tpm <- lapply(Sys.glob(my.files), read.table, header = TRUE)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

names <- sapply(strsplit(my.files, split="\\/"), function(x)x[8])
names <- gsub("_quant", "", names)
names(tpm) <- names

tpm_per <- function(DF) {
  DF$Percent <- (DF$TPM /(sum(DF$TPM)))*100
  DF <- DF %>% select(Name, Percent)
  return(DF)
}

tpm2 <- lapply(tpm, tpm_per)
wide_tpm <- dplyr:::reduce(tpm2, full_join, by = "Name")
names <- append(names, 'Name', after = 0)
colnames(wide_tpm) <- names

keep_tpm <- wide_tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 10, Name))

keep_tpm <- keep_tpm %>% mutate(Name = if_else(Names == "TRUE", Name, "zVSGome"))

keep_tpm <- subset(keep_tpm, select = -c(Names) )

long <- melt(keep_tpm, id.vars = c("Name"), variable.name = "isolate", value.name = "P")
long <- merge(x=stats, y=long, by.x="isolate", by.y="variable")[]

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$Name)))
coul <- append(coul, "black")

png("~/Desktop/mouse.png", units="in", width=12, height=12, res=300)
ggbarplot(long, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

#
my.files <-  list.files(list.dirs(path = "/Users/goldriev/pkgs/analyse/full_cds/cow", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)
tpm <- lapply(Sys.glob(my.files), read.table, header = TRUE)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

names <- sapply(strsplit(my.files, split="\\/"), function(x)x[8])
names <- gsub("_quant", "", names)
names(tpm) <- names

tpm_per <- function(DF) {
  DF$Percent <- (DF$TPM /(sum(DF$TPM)))*100
  DF <- DF %>% select(Name, Percent)
  return(DF)
}

tpm2 <- lapply(tpm, tpm_per)
wide_tpm <- dplyr:::reduce(tpm2, full_join, by = "Name")
names <- append(names, 'Name', after = 0)
colnames(wide_tpm) <- names

keep_tpm <- wide_tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 10, Name))

keep_tpm <- keep_tpm %>% mutate(Name = if_else(Names == "TRUE", Name, "zVSGome"))

keep_tpm <- subset(keep_tpm, select = -c(Names) )

long <- melt(keep_tpm, id.vars = c("Name"), variable.name = "isolate", value.name = "P")
long <- merge(x=stats, y=long, by.x="isolate", by.y="variable")[]

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$Name)))
coul <- append(coul, "black")
long$test <- as.numeric(long$test)

png("~/Desktop/cow.png", units="in", width=12, height=12, res=300)
ggbarplot(long, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()


filter_tpm <- filter(keep_tpm, Name != "zVSGome")
row.names(filter_tpm) <- filter_tpm$'Name'
filter_tpm <- filter_tpm[,-1]
mat <- filter_tpm - rowMeans(filter_tpm)
mat <- log(mat)
mat[is.na(mat)] = 0
breaksList = seq(-7, 7)

pheatmap(mat, show_rownames=FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), breaks = breaksList)

# QC stats

assembly <- subset(stats, select = c("isolate", "type", "stage","assembled_transcripts", "orf", "vsg_blast", "not_blast", "cd.hit", "parasites"))
long <- melt(assembly, id.vars = c("isolate","type", "stage", "parasites"), variable.name = "transcripts", value.name = "transcripts")
mouse_stats <- stats[stats$host == 'mouse', ]

a <- ggboxplot(mouse_stats, x = "type", y = "mapping", color = "stage", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Mapping rate")

b <- ggboxplot(mouse_stats, x = "type", y = "vsgs", color = "stage", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("VSG count")

c <- ggboxplot(mouse_stats, x = "type", y = "log(parasites/vsgs)", color = "stage", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("log(parasite count/ VSG count")

png("~/Desktop/vsgseq2_summary.png", units="in", width=12, height=12, res=300)
ggarrange(a,b,c, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

png("~/Desktop/vsg_summary.png", units="in", width=12, height=10, res=300)
ggarrange(f, c, d, e, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
