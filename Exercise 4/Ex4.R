library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(gplots)
library(purrr)
library(xtable)
library(DESeq2)

s_counts <- read.table("./Data/16s_counts.txt", comment.char="", quote="", sep="\t", header= TRUE)
s_annotation <- read.table("./Data/16s_annotation.txt", comment.char="", quote="", sep="\t", header= TRUE)


colSums(s_counts)

# Filtering data
s_counts <- s_counts[rowSums(s_counts)>=5,]


#### PCA #####
pca=prcomp(t(s_counts))
summary(pca)
data <- data.frame(pca$x)
# PC1 and PC2 acount for most of the variation
ggplot(data, aes(y = PC2, x = PC1)) +
  geom_point(aes(color = c("High", "High", "High", "Low", "Low", "Low")), size = 4) +
  labs(color = "Oiliness")
ggsave("./Figures/pca_Oiliness.png", height = 6, width = 8, units = "cm")

# Log-transformed
pca=prcomp(t(log(1+s_counts)))
summary(pca)
data <- data.frame(pca$x)
# PC1 and PC2 acount for most of the variation
ggplot(data, aes(y = PC2, x = PC1)) +
  geom_point(aes(color = c("High", "High", "High", "Low", "Low", "Low")), size = 4) +
  labs(color = "Oiliness")
ggsave("./Figures/pca_log_Oiliness.png", height = 6, width = 8, units = "cm")




##### Rarefaction ####
RareFaction = function(OTU_name,count,depth){
  reads <- rep(OTU_name, times = count)
  reads.sample<-sample(reads,size=depth, replace=FALSE)
  counts.sample<-as.data.frame(table(reads.sample))
return(counts.sample)
  }

depth <- min(colSums(s_counts))
HC1_rf <- RareFaction(rownames(s_counts), s_counts$HC1, depth)
HC2_rf <- RareFaction(rownames(s_counts), s_counts$HC2, depth)
HC3_rf <- RareFaction(rownames(s_counts), s_counts$HC3, depth)
LC1_rf <- RareFaction(rownames(s_counts), s_counts$LC1, depth)
LC2_rf <- RareFaction(rownames(s_counts), s_counts$LC2, depth)
LC3_rf <- RareFaction(rownames(s_counts), s_counts$LC3, depth)

colnames(HC1_rf) <- c("index", "HC1")
colnames(HC2_rf) <- c("index", "HC2")
colnames(HC3_rf) <- c("index", "HC3")
colnames(LC1_rf) <- c("index", "LC1")
colnames(LC2_rf) <- c("index", "LC2")
colnames(LC3_rf) <- c("index", "LC3")

rf_counts <- purrr::reduce(
  list(HC1_rf, HC2_rf, HC3_rf, LC1_rf, LC2_rf, LC3_rf),
  function(x, y) merge(x, y, by = "index", all = TRUE)
)
rf_counts[is.na(rf_counts)] <- 0
rownames(rf_counts) <- rf_counts$index
rf_counts <- rf_counts[,-1]





#### Diversity indices ####

# Diversity
diveristy <- as.data.frame(colSums(rf_counts>0))
diveristy$Sample <- rownames(evenness)
colnames(diveristy) <- c("Diversity", "Sample")
ggplot(diveristy, aes(x = Sample, y = Diversity)) +
  geom_bar(stat = "identity")
ggsave("./Figures/diversity.png")

# Evenness
shannon_sum <- -rf_counts/depth * log(rf_counts/depth)
shannon_sum[is.na(shannon_sum)] <- 0
evenness <- as.data.frame(colSums(shannon_sum))
evenness$Sample <- rownames(evenness)
colnames(evenness)[1] <- "Evenness"

ggplot(evenness, aes(x = Sample, y = Evenness)) +
  geom_bar(stat = "identity")
ggsave("./Figures/evenness.png")


#### Linear Count Model ####
design.matrix <- data.frame(exposure = c(1,1,1,0,0,0))
counts.ds <- DESeqDataSetFromMatrix(countData = rf_counts, design.matrix, design = ~exposure)

res.ds <- DESeq(counts.ds)
results.ds <- results(res.ds, independentFiltering = FALSE, cooksCutoff = FALSE)
results.ds$OTU.ID <- rownames(results.ds)


# Merge DESeq results with annotation
results <- purrr::reduce(
  list(as.data.frame(results.ds), s_annotation),
  function(x, y) merge(x, y, by = "OTU.ID", all = TRUE)
)

# Print annotation for latex
results <- results[order(results$padj),]
subset_results <- results[1:10, c("padj","log2FoldChange", "Class", "Order", "Family", "Genus")]
latex_table <- xtable(subset_results, 
                      caption = "Subset of Results", 
                      label = "tab:subset_results")
print(latex_table, type = "latex", include.rownames = FALSE, scientific = TRUE)


