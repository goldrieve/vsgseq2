library (dplyr)
library(RColorBrewer)
library(reshape)
library(ggpubr)
library(pheatmap)
library(tidyverse)

#Read in the meta data and modify it for tximport
stats <- read.csv("~/Google Drive/My Drive/vsg/vsgseq2/meta.csv")
orthogroups <- read.csv("~/Desktop/Orthogroups.csv")
orthogroups_long <- orthogroups %>% separate_rows(Name, sep = ",")

my.files <-  list.files(list.dirs(path = "/Users/goldriev/pkgs/analyse/salmon", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)
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

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$Name)))
coul <- append(coul, "black")

cow_comb <- long[long$host == 'cow', ]
mouse_comb <- long[long$host == 'mouse', ]
cow_comb$test <- as.numeric(cow_comb$test)

ggbarplot(long, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)

png("~/Desktop/combined_cow.png", units="in", width=12, height=12, res=300)
ggbarplot(cow_comb, x = "test", y = "value", color = "Orthogroup", fill = "Orthogroup", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

png("~/Desktop/combined_mouse.png", units="in", width=12, height=12, res=300)
ggbarplot(mouse_comb, x = "test", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(type + host ~ stage, scales = 'free_x', ncol = 2) +
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

a <- ggboxplot(long, x = "variable", y = "value", color = "type", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Transcript count") +
  scale_fill_manual(values = c("lightgrey", "white")) 

b <- ggboxplot(long, x = "variable", y = "log(value)", color = "type", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Transcript count (log)") +
  scale_fill_manual(values = c("lightgrey", "white")) 

vsgs <- filter(long,variable == 'cd.hit')

c <- ggboxplot(vsgs, x = "type", y = "log(value)", color = "type", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Transcript count (log)") +
  scale_fill_manual(values = c("lightgrey", "white")) 

d <- ggboxplot(vsgs, x = "type", y = "log(parasites/value)", color = "type", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Parasite count/ transcript count (log)") +
  scale_fill_manual(values = c("lightgrey", "white")) 

mapping <- subset(stats, select = c("isolate", "type", "stage","vsgseq2_mapping", "VSGSeq_mapping", "parasites"))
long <- melt(mapping, id.vars = c("isolate","type", "stage", "parasites"), variable.name = "transcripts", value.name = "transcripts")

e <- ggboxplot(long, x = "type", y = "value", color = "variable", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Mapping rate (%)") +
  scale_fill_manual(values = c("lightgrey", "white")) 

f <- ggboxplot(vsgs, x = "type", y = "(parasites)", color = "type", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Parasite count") +
  scale_fill_manual(values = c("lightgrey", "white")) +
  scale_y_continuous(label=scientific_10) 

png("~/Desktop/vsg_transcripts.png", units="in", width=12, height=12, res=300)
ggarrange(a,b, ncol = 1, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

png("~/Desktop/vsg_summary.png", units="in", width=12, height=10, res=300)
ggarrange(f, c, d, e, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
