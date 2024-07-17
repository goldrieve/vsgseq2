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
tpm <- read.csv("pkgs/vsgseq2/results/final/summary/tpm.csv", check.names = FALSE)
numeric_columns <- sapply(tpm[-1], is.numeric)
column_sums <- colSums(tpm[-1][, numeric_columns])
tpm <- tpm
tpm[-1][, numeric_columns] <- tpm[-1][, numeric_columns] / column_sums * 100

long_data <- orthogroups %>% 
  gather(key = "sample", value = "value", -Orthogroup)

long_data <- long_data %>% 
  separate_rows(value, sep = ", ")

long_data$value <- sub("_1$", "", long_data$value)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

tpm$id <- tpm$vsg_id

tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 10, vsg_id))

tpm <- tpm %>% mutate(vsg_id = if_else(Names == "TRUE", vsg_id, "other"))
tpm <- subset(tpm, select = -c(Names) )
tpm$vsg_id <- factor(tpm$vsg_id, levels = c(unique(tpm$vsg_id[tpm$vsg_id != "other"]), "other"))

long <- melt(tpm, id.vars = c("id", "vsg_id"), variable.name = "isolate", value.name = "P")
long <- merge(x=stats, y=long, by.x="isolate", by.y="variable")[]
long <- merge(x = long, y = long_data, by.x = "id", by.y = "value", all.x = TRUE)
long <- long %>% replace_na(list(col1 = "Other", col2 = "Other"))

long$value_log <- log(1+ long$value)
cow_long <- long[long$host == 'cow', ]
mouse_long <- long[long$host == 'mouse', ]
cow_long$stage <- as.numeric(cow_long$old_stage)

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$vsg_id)))
coul <- append(coul, "black")

cow_long <- cow_long %>%
  mutate(numeric_value = as.numeric(x_axis))

# Step 2: Arrange the dataframe by the new column in descending order
cow_long <- cow_long %>%
  arrange((numeric_value))

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/cow.png", units="in", width=15, height=10, res=300)
ggbarplot(cow_long, x = "x_axis", y = "value", color = "vsg_id", fill = "vsg_id", legend = "right") +
  facet_wrap(~ old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  xlab ("DPI") +
  ylab ("TPM")
dev.off()

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/mouse.png", units="in", width=15, height=10, res=300)
ggbarplot(mouse_long, x = "x_axis", y = "value", color = "vsg_id", fill = "vsg_id", legend = "right") +
  facet_wrap(type ~ old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  xlab ("Sample") +
  ylab ("TPM")
dev.off()

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$Orthogroup)))
coul <- append(coul, "black")

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/cow_orthogroup_colour.png", units="in", width=15, height=10, res=300)
ggbarplot(cow_long, x = "x_axis", y = "value", color = "Orthogroup", fill = "Orthogroup", legend = "right") +
  facet_wrap(type~old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/mouse_orthogroup_colour.png", units="in", width=15, height=10, res=300)
ggbarplot(mouse_long, x = "x_axis", y = "value", color = "Orthogroup", fill = "Orthogroup", legend = "right") +
  facet_wrap(type~old_stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

coul <- brewer.pal(5,"Set3") 

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/cow_orthogroup.png", units="in", width=15, height=10, res=300)
ggdotplot(cow_long, x = "x_axis", y = "value_log", color = "type",
       add = c("scatter"), facet.by = c("Orthogroup")) +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

coul <- brewer.pal(3,"Set3")
mouse_long$old_stage <- factor(mouse_long$old_stage, levels = c("Early", "Late"))

png("Google Drive/My Drive/vsg/vsgseq2/report/figures/mouse_orthogroup.png", units="in", width=15, height=10, res=300)
ggline(mouse_long, x = "old_stage", y = "value_log", color = "type",
       add = c("boxplot"), facet.by = "Orthogroup") +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul)
dev.off()

cow_filtered <- cow_long[cow_long$value > 0.1, ]

top <- cow_long %>%
  group_by(Orthogroup, x_axis) %>%
  slice(which.max(value))

cow_long$x_axis <- as.numeric(as.character(cow_long$x_axis))

dodge <- position_dodge(width = 0.9)

ggplot(cow_filtered, aes(x = factor(x_axis), y = (value+1), fill = type)) +
  geom_point(aes(color = type), alpha = 0.1)+
  geom_smooth(aes(group = type, color = type), method = "loess", se = FALSE) +  # Corrected color mapping
  theme_minimal() +
  labs(title = "Violin Plot of Values by Group",
       x = "Group",
       y = "Value (log-transformed)") +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~Orthogroup)

ggplot(cow_filtered, aes(x = x_axis, y = (value))) +
  geom_point(aes(color = type), alpha = 0.1)+
  geom_smooth(method= lm, se = F, color = "black") + 
  facet_wrap(~Orthogroup) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(cow_filtered, aes(x = x_axis, y = log(value), colour = type)) +
  geom_point(alpha = 0.1)+
  geom_smooth(method= lm, se = F) + 
  facet_wrap(~Orthogroup) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

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
