library(tximport)

txi.salmon <- tximport(salmon_dirs, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
