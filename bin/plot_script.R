# Load packages
pacman::p_load(dplyr, RColorBrewer, reshape, ggpubr)

# Read in data
tpm <- read.csv("~/pkgs/vsgseq2/data/reads/results/tutorial/summary/tpm/cluster_tpm.csv")
meta <- read.csv("~/pkgs/vsgseq2/data/reads/meta.csv")
count <- read.csv("~/pkgs/vsgseq2/data/reads/results/tutorial/summary/vsgs/vsg_count.csv")

# Filter rows in tpm where any numeric column has a value greater than 100000
keep_tpm <- tpm %>%
  mutate(Names = if_any(where(is.numeric), ~ . > 100000, cluster))
keep_tpm$cluster <- as.character(keep_tpm$cluster)

# Replace 'vsg_id' with 'other' where 'Names' is not TRUE
keep_tpm <- keep_tpm %>% mutate(cluster = if_else(Names == "TRUE", cluster, "other"))

# Remove 'Names' column
keep_tpm <- subset(keep_tpm, select = -c(Names) )

# Reshape data from wide to long format
long <- melt(keep_tpm, id.vars = c("cluster"), variable.name = "isolate", value.name = "P")

# Merge 'long' and 'meta' data frames
long <- merge(x=meta, y=long, by.x="isolate", by.y="variable")[]

# Define color palette
coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(-1 + length(unique(long$cluster)))
coul <- append(coul, "black")

# Create bar plot
a <- ggbarplot(long, x = "isolate", y = "value",
               add.params = list(size = 2), color = "cluster", fill = "cluster") +
  facet_wrap(~ stage, scales = 'free_x', ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) +
  ylab("TPM")

# Merge 'count' and 'meta' data frames

count <- merge(x=meta, y=count, by.x="isolate", by.y="sample")[]

# Create box plot
b <- ggboxplot(count, x = "stage", y = "VSGs", legend = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("VSG count")

# Save plots to a PNG file
png("figures/vsg_summary.png", units="in", width=10, height=12, res=300)
ggarrange(a, b, ncol = 1, nrow = 2,  common.legend = F, legend="right", align = c("hv"), labels = "auto", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
dev.off()