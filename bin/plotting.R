library (dplyr)
#Read in the meta data and modify it for tximport

my.files <-  list.files(list.dirs(path = "/Users/goldriev/pkgs/vsgseq2/results", full.names = TRUE, recursive = TRUE), pattern = "quant.sf", full.names = TRUE)
tpm <- lapply(Sys.glob(my.files), read.table, header = TRUE)


names(tpm) <- sapply(strsplit(my.files, split="\\/"), function(x)x[7])
tpm_per <- function(DF) {
  DF$Percent <- (DF$TPM /(sum(DF$TPM)))*100
  DF <- DF %>% select(Name, Percent)
  return(DF)
}

tpm2 <- lapply(tpm, tpm_per)
wide_tpm <- dplyr:::reduce(tpm2, full_join, by = "Name")
wide_tpm$sum <- rowSums(wide_tpm[,-1])
wide_tpm$V5 <- wide_tpm$sum / (sum(wide_tpm$sum))*100
wide_tpm$names <- ifelse(wide_tpm$V5 < 1, 'other',
                          ifelse(wide_tpm$V5 >= 1, wide_tpm$Name, wide_tpm$Name))
wide_tpm <- subset(wide_tpm, select = -c(Name,sum,V5) )
long <- melt(wide_tpm, id.vars = c("names"), variable.name = "Percent", value.name = "P")

coul <- brewer.pal(11,"Set3") 

ggbarplot(long, x = "Percent", y = "P", color = "names", fill = "names") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = coul) +
  scale_color_manual(values = coul) 
