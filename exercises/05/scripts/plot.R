library(ggplot2)

gene <- read.csv2(snakemake@input[[1]], dec = ",")
gene$id <- as.numeric(gene$id)

pdf(file = snakemake@output[[1]])
ggplot(gene, aes(x = id)) +
    geom_density(alpha=.2, fill="#FF6666") +
    ggtitle(paste0("Plot "), snakemake@wildcards[["file"]])     
dev.off()
