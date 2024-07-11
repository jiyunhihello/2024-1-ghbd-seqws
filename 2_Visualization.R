
# 1 - Visualize the results with Volcano plot
data$Significance <- ifelse(data$padj < 0.05 & abs(data$log2FoldChange) > 2, "Pass both p-value & Log2 fold change", ifelse(data$padj < 0.05, "Pass p-value cutoff", "Not significant"))

library(ggplot2)
library(plotly)

volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue), text = rownames(data))) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "Treatment_vs_Control", x = "Log2 fold change", y = "-Log10 p-value")

interactive_plot <- ggplotly(volcano_plot, tooltip = "text")

interactive_plot

# 2 - Visualize the results with Heatmap
library(pheatmap)
library(dplyr)

count_data <- total[, 8:19]
rownames(count_data) <- total$Row.names

meta$treatment <- ifelse(meta$treatment == "Clarithromycin", "Treatment", "Control")

anno_col_info <- meta %>% select(Run, treatment) %>% tibble::column_to_rownames("Run")

anno_info_colors <- list(treatment = c(Treatment = "lightgrey", Control = "black"))

pheatmap(count_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = anno_col_info,
         annotation_colors = anno_info_colors,
         main = "Clustering of differentially expressed genes")
