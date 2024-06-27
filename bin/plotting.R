library(dplyr)
library(RColorBrewer)
library(reshape)
library(ggpubr)
library(pheatmap)
library(tidyverse)

#Read in the meta data and modify it for tximport
stats <- read.csv("Google Drive/My Drive/vsg/vsgseq2/data/meta.csv")
orthogroups <- read.csv("~/pkgs/analyse/all_data_results_no_vsgdb/protein/OrthoFinder/Results_Jun18/Orthogroups/Orthogroups.csv")
orthogroups_long <- orthogroups %>% separate_rows(Orthogroup, sep = ",")
my.files <-  list.files(list.dirs(path = "Desktop/minion", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)

long_data <- orthogroups %>% 
  gather(key = "sample", value = "value", -Orthogroup)

long_data <- long_data %>% 
  separate_rows(value, sep = ", ")

long_data$value <- sub("_1$", "", long_data$value)

tpm <- lapply(Sys.glob(my.files), read.table, header = TRUE)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

names <- sapply(strsplit(my.files, split="\\/"), function(x)x[4])
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
long <- merge(x = long, y = long_data, by.x = "Name", by.y = "value", all.x = TRUE)
long <- long %>% replace_na(list(col1 = "Other", col2 = "Other"))

png("~/Desktop/minion.png", units="in", width=10, height=8, res=300)
ggbarplot(long, x = "isolate", y = "value",
          add.params = list(size = 2), color = "Name", fill = "Name") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap("old_stage", scales = 'free_x') +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

cow_comb <- long[long$host == 'cow', ]
cow_comb <- cow_comb[cow_comb$bleed == 'small', ]
mouse_comb <- long[long$host == 'mouse', ]
mouse_comb <- mouse_comb[mouse_comb$type == 'WT', ]
cow_comb$stage <- as.numeric(cow_comb$stage)


a <- ggbarplot(long, x = "sample", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(~ cow, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

b <- ggbarplot(long, x = "sample", y = "value", color = "Name", fill = "Name", legend = "right") +
  facet_wrap(~ pool, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

png("~/Desktop/minion.png", units="in", width=10, height=10, res=300)
ggarrange(a,b, ncol = 1, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$Orthogroup)))
coul <- append(coul, "black")


long <- long %>%
  mutate(numeric_value = as.numeric(x_axis))

# Step 2: Arrange the dataframe by the new column in descending order
long <- long %>%
  arrange((numeric_value))


a <- ggbarplot(long, x = "x_axis", y = "value.x", color = "Orthogroup", fill = "Orthogroup", legend = "right") +
  facet_wrap(type~old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)

b <- ggbarplot(long, x = "x_axis", y = "value", color = "Orthogroup", fill = "Orthogroup", legend = "right") +
  facet_wrap(type~old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)

ggarrange(a,b, ncol = 1, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))

png("~/Desktop/orthogroup_major.png", units="in", width=15, height=15, res=300)
a
dev.off()

png("~/Desktop/orthogroup_all.png", units="in", width=15, height=15, res=300)
b
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

coul <- brewer.pal(6,"Set3") 
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

stats <- filter(stats, isolate != "908_d11")
png("~/Desktop/vsg_summary.png", units="in", width=12, height=10, res=300)
ggboxplot(stats, x = "x_axis", y = "log(parasites/vsgs)", color = "old_stage", fill = "old_stage", legend = "right") +
  facet_wrap(type ~ old_stage, scales = 'free_x', ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("VSG count")
dev.off()

c <- ggboxplot(mouse_stats, x = "type", y = "(parasites/vsgs)", color = "stage", fill = "stage", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("log(parasite count/ VSG count")

png("~/Desktop/vsgseq2_summary.png", units="in", width=12, height=12, res=300)
ggarrange(a,b,c, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()

png("~/Desktop/vsg_summary.png", units="in", width=12, height=10, res=300)
ggarrange(f, c, d, e, ncol = 2, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()
