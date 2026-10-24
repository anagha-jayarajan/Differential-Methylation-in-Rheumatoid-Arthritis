# ---------------------------
# 1. Load libraries
# ---------------------------
library(minfi)
library(ChAMP)
library(sesame)
library(sesameData)
library(ggplot2)
library(pheatmap)

# 2. Load Data
# ---------------------------
rgSet <- read.metharray.exp(base = getwd())

targets <- read.csv("Sample_sheet.csv", stringsAsFactors = FALSE)

colnames(targets)[colnames(targets) == "sample_name"] <- "Sample_Name"
pData(rgSet)$Sample_Group <- targets$Sample_Group
pData(rgSet)$Sample_Name <- targets$Sample_Name

# 3. Quality Control
# ---------------------------
qc <- getQC(preprocessRaw(rgSet))
pdf("QC_plot_minfi.pdf")
plotQC(qc)

# 4. Filtering
# ---------------------------
detP <- detectionP(rgSet)
keep <- rowMeans(detP) < 0.01
rgSet_filtered <- rgSet[keep, ]

# 5. Normalization
# ---------------------------
mSet <- preprocessFunnorm(rgSet_filtered)
beta_values <- getBeta(mSet)

# 6. ChAMP (Using minfi output)
# ---------------------------
cat("Running ChAMP pipeline...\n")
champ_beta <- beta_values
champ_pd <- data.frame(
  Sample_Group = targets$Sample_Group
)
rownames(champ_pd) <- colnames(champ_beta)

# Normalize (ChAMP)
champ_norm <- champ.norm(
  beta = champ_beta,
  method = "BMIQ"
)
champ_dmp <- champ.DMP(
  beta = champ_norm,
  pheno = champ_pd$Sample_Group
)

# Save results
champ_results <- champ_dmp$Disease_to_Control
write.csv(champ_results, "Differential_Methylation_ChAMP.csv")

group <- factor(targets$Sample_Group,
                levels = c("Control", "Disease"))

pca <- prcomp(t(beta_values), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = group
)
pdf("PCA_plot.pdf")
ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse() +
  theme_classic()
table(group)
ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  theme_classic(base_size = 14) +
  labs(title = "PCA of DNA Methylation")
# heatmap
top_probes <- rownames(champ_results)[1:50]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(beta_values)
pdf("Heatmap.pdf")
pheatmap(beta_values[top_probes, ],
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE)
volcano_df <- champ_results
# Add significance column
volcano_df$Significant <- "Not Significant"
volcano_df$Significant[
  volcano_df$adj.P.Val < 0.05 & abs(volcano_df$logFC) > 0.2
] <- "Significant"
pdf("Volcano_plot.pdf", width = 7, height = 6)
ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano Plot (ChAMP DMPs)",
    x = "Log Fold Change",
    y = "-log10(adj.P.Val)"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

# Sesame Pipeline
# ---------------------------
cat("Running sesame pipeline...\n")
idat_dir <- getwd()
sdfs <- openSesame(idat_dir)
sesame_beta <- do.call(cbind, lapply(sdfs, getBetas))                       
sesame_beta <- do.call(cbind, sdfs)
sesame_beta <- sdfs
targets <- read.csv("Sample_sheet.csv", stringsAsFactors = FALSE)
targets$Sample_Group <- tolower(targets$Sample_Group)
targets <- targets[match(colnames(sesame_beta), targets$Sample_Name), ]
group1 <- sesame_beta[, targets$Sample_Group == "disease"]
group2 <- sesame_beta[, targets$Sample_Group == "control"]
mean1 <- rowMeans(group1, na.rm = TRUE)
mean2 <- rowMeans(group2, na.rm = TRUE)
delta_beta <- mean1 - mean2
pvals <- apply(sesame_beta, 1, function(x) {
  
  group1 <- x[targets$Sample_Group == "disease"]
  group2 <- x[targets$Sample_Group == "control"]
  
  # Remove NA values
  group1 <- group1[!is.na(group1)]
  group2 <- group2[!is.na(group2)]
  
  # Check if enough samples exist
  if(length(group1) < 2 || length(group2) < 2) {
    return(NA)
  }
  
  t.test(group1, group2)$p.value
})
adj_pvals <- p.adjust(pvals, method = "BH")

targets$Sample_Group <- tolower(targets$Sample_Group)
targets <- targets[match(colnames(sesame_beta), targets$Sample_Name), ]
pvals <- apply(sesame_beta, 1, function(x) {
  
  group1 <- x[targets$Sample_Group == "disease"]
  group2 <- x[targets$Sample_Group == "control"]
  
  # Remove NA values
  group1 <- group1[!is.na(group1)]
  group2 <- group2[!is.na(group2)]
  
  # Skip if not enough samples
  if(length(group1) < 2 || length(group2) < 2) {
    return(NA)
  }
  
  t.test(group1, group2)$p.value
})
adj_pvals <- p.adjust(pvals, method = "BH")
table(targets$Sample_Group)
results <- data.frame(
  CpG = rownames(sesame_beta),
  p_value = pvals,
  adj_p_value = adj_pvals
)
write.csv(results, "sesame_DMP_manual.csv", row.names = FALSE)
# heatmap
top_probes <- rownames(champ_results)[1:50]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(beta_values)
pdf("Heatmap.pdf", width = 10, height = 8)
pheatmap(beta_values[top_probes, ],
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 Differentially Methylated Probes")
dev.off()
volcano_df <- champ_results
volcano_df$Significant <- "Not Significant"
volcano_df$Significant[
  volcano_df$adj.P.Val < 0.05 & abs(volcano_df$logFC) > 0.2
] <- "Significant"
pdf("Volcano_plot.pdf", width = 7, height = 6)
p <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano Plot (ChAMP DMPs)",
    x = "Log Fold Change",
    y = "-log10(adj.P.Val)"
  ) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")
print(p)
dev.off()

# heatmap
top_probes <- rownames(champ_results)[1:50]
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(beta_values)
pdf("Heatmap.pdf", width = 10, height = 8)
pheatmap(beta_values[top_probes, ],
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 Differentially Methylated Probes")
dev.off()
