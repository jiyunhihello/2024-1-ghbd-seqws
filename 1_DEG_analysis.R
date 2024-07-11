library(DESeq2)
library(dplyr)

# 1 - Import libraries and data
raw_counts <- read.table("gene_fC.txt", sep='\t', header=T) %>% tibble::column_to_rownames("Geneid")

meta <- read.table("meta_12.txt", sep='\t', hedaer=T)
meta$treatment <- as.factor(meta$treatment)

# 2 - Create the DESeqDataSet object
raw_counts_filtered <- raw_counts[, colnames(raw_counts) %in% meta$Run]
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, colData = meta, design = ~ treatment)

# 3 - Run the DE analysis
dds <- DESeq(dds)
res <- results(dds)

# 4 - Extract the table of differentially expressed genes
diff_genes = res %>% as.data.frame() %>% tibble::rownames_to_column("genes") %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange), desc(padj))
