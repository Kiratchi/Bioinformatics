library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(gplots)



x = read.table("./Data/counts_matrix.txt")

metadata = read.table("./Data/metadata.txt", sep="\t", header=TRUE)
metadata$diagnosis <- relevel( as.factor(metadata$diagnosis), ref = "Not IBD" )
metadata$Sex <- as.factor(metadata$Sex)

gene_annotation = read.table("./Data/geneAnnotation.txt", sep="\t", header=TRUE)
rownames(gene_annotation) <- gene_annotation$ensembl_gene_id

head(metadata)
summary(metadata)

table(metadata$Sex)
table(metadata$diagnosis)

gene2586=x["ENSG00000002586",]
gene4809=x["ENSG00000004809",]


nr_patients = length(x)


# Filtering x (with genes expressed at non 0 level for +25% of samples)
index_keep = (rowSums(x>0) / nr_patients ) >= 0.25
x.filtered <- x[index_keep,]

# Calculating CPM and log-transforming it
CPM <- sweep(x.filtered+1, 2, 10^6/colSums(x.filtered), "*") 
logCPM = log(CPM)
logx.filtered <- log(x.filtered+1)


# Box-plot normalized data
logCPM_long <- pivot_longer(logCPM, cols = everything(), names_to = "person_id", values_to = "expression")
head(logCPM_long)
ggplot(data=logCPM_long, aes(x=person_id, y = expression)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Expression by Person - Normalized", 
       x = "Person ID", 
       y = "log(Expression)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5))
ggsave("./Figures/boxplot_normalized.png",  width = 20, height = 6, units = "cm")

# Box-plot non-normalized data = x.filtered
logx.filtered_long <- pivot_longer(logx.filtered, cols = everything(), names_to = "person_id", values_to = "expression")
head(logx.filtered_long)
ggplot(data=logx.filtered_long, aes(x=person_id, y = expression)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Expression by Person - Non Normalized", 
       x = "Person ID", 
       y = "log(Expression)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5))
ggsave("./Figures/boxplot_non-normalized.png",  width = 20, height = 6, units = "cm")

# Scatter plot sample 1 and 41
ggplot(data=logCPM , aes(x=logCPM[,1] , y=logCPM[,41] )) + geom_point(size = 0.5) +
  labs( x = "log(expression) for Sample 1", 
       y = "log(expression) for Sample 41" ) +
  geom_abline(aes(intercept = 0, slope = 1), color = "orange", linewidth = 0.5)
ggsave("./Figures/scatterplot_41_1.png", width = 8, height = 8, units = "cm")


gene1 <- data.frame(
          logCPM_gene1 = t(logCPM[1,]),
          diagnosis = relevel( as.factor(metadata$diagnosis), ref = "Not IBD" ),
          age = metadata$age.at.diagnosis,
          gender = as.factor(metadata$Sex)
          )

colnames(gene1)[1] <- "logCPM_gene1"


# Creating a boxplot for gene 1 to see if there is a 
# difference between having the disease or not
ggplot(data = gene1, aes(x = diagnosis, y = logCPM_gene1, fill = diagnosis)) + 
  geom_boxplot() + theme(legend.position = "none")
ggsave("./Figures/boxplot_gene1.png")

fit1 <- lm( logCPM_gene1 ~ diagnosis, data=gene1)
summary(fit1)

fit2 <- lm( logCPM_gene1 ~ diagnosis + age + gender, data=gene1)
summary(fit2)


nr_genes <- dim(logCPM)[1]

fit1_df <- data.frame(
        p_value = rep(0, nr_genes),
        coef = rep(0, nr_genes)
) 

fit2_df <- data.frame(
  diagnosis_p_value = rep(0, nr_genes),
  age_p_value = rep(0, nr_genes),
  gender_p_value = rep(0, nr_genes),
  diagnosis_coef = rep(0, nr_genes),
  age_coef = rep(0, nr_genes),
  gender_coef = rep(0, nr_genes)
) 

for(i in 1:nr_genes){
  fit1 <- lm( t(logCPM[i,]) ~ metadata$diagnosis)
  fit1_df[i,] <- c(summary(fit1)$coef[2,4], summary(fit1)$coef[2,1])
  
  
  fit2 <- lm( t(logCPM[i,]) ~ metadata$diagnosis + metadata$age.at.diagnosis 
              + metadata$Sex )
  fit2_df[i,] <- c(summary(fit2)$coef[2,4], summary(fit2)$coef[3,4], summary(fit2)$coef[4,4],
                   summary(fit2)$coef[2,1], summary(fit2)$coef[3,1], summary(fit2)$coef[4,1])
  
}
rownames(fit1_df) <- rownames(logCPM)
rownames(fit2_df) <- rownames(logCPM)


# Adjusting p-values
fit1_df$p_value_adj <- p.adjust(fit1_df$p_value, method="fdr")
fit2_df$diagnosis_p_value_adj <-p.adjust(fit2_df$diagnosis_p_value, method="fdr")
fit2_df$age_p_value_adj <- p.adjust(fit2_df$age_p_value, method="fdr")
fit2_df$gender_p_value_adj <- p.adjust(fit2_df$gender_p_value, method="fdr")

fit1_df <- fit1_df[order(fit1_df$p_value_adj), ]
ggplot(data = fit1_df, aes(x = 1:nr_genes, y = -log10(p_value_adj))) + 
  geom_point(size=0.7) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "Significant Gene Order", 
       y = "-log10(p_adj)")
ggsave("./Figures/fit1_p_values.png")

fit2_df <- fit2_df[order(fit1_df$p_value_adj), ]
ggplot(data = fit2_df, aes(x = 1:nr_genes, y = -log10(diagnosis_p_value_adj))) + 
  geom_point(size=0.7) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs( x = "Significant Gene Order", 
       y = "-log10(p_adj)")
ggsave("./Figures/fit2_diag_p_values.png")

fit2_df <- fit2_df[order(fit2_df$age_p_value_adj), ]
ggplot(data = fit2_df, aes(x = 1:nr_genes, y = -log10(age_p_value_adj))) + 
  geom_point(size=0.1) + 
  labs(title = "", 
       x = "Gene", 
       y = "-log10(p_adj)")
ggsave("./Figures/fit2_age_p_values.png")

ggplot(data = fit2_df, aes(x = 1:nr_genes, y = -log10(gender_p_value_adj))) + 
  geom_point(size=0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "", 
       x = "Gene", 
       y = "-log10(p_adj)")
ggsave("./Figures/fit2_gender_p_values.png")



print(sum(fit1_df$p_value_adj<0.05))
print(sum(fit1_df$p_value_adj<0.05 & fit1_df$coef>0))
print(sum(fit1_df$p_value_adj<0.05 & fit1_df$coef<0))


print(sum(fit2_df$diagnosis_p_value_adj<0.05))
print(sum(fit2_df$diagnosis_p_value_adj<0.05 & fit2_df$diagnosis_coef>0))
print(sum(fit2_df$diagnosis_p_value_adj<0.05 & fit2_df$diagnosis_coef<0))

print(sum(fit2_df$age_p_value_adj<0.05))
print(sum(fit2_df$gender_p_value_adj<0.05))

fit2_df <- fit2_df[order(fit2_df$age_p_value_adj), ]
head(fit2_df)

fit2_df <- fit2_df[order(fit2_df$gender_p_value_adj), ]
head(fit2_df)


fit1_df <- fit1_df[order(fit1_df$p_value_adj), ]
most_significant_genes <- rownames(fit2_df)[1:5]
print(gene_annotation[most_significant_genes, c("ensembl_gene_id", "description")])


###### Heirarchical clustering ####
fit2_df <- fit2_df[order(fit2_df$diagnosis_p_value_adj), ]
genes_sig <- rownames(fit2_df[1:100,])
logCPM_sig <- as.matrix(logCPM[genes_sig,])

mypalette = brewer.pal(11,"RdYlBu")
morecols = colorRampPalette(mypalette)
mycols=rev(morecols(255))
column_cols=c("#F8766D","#00BFC4")[metadata$diagnosis]



pdf("./Figures/top100sigGenesHeatmap.pdf",height=9,width=10)
heatmap.2(logCPM_sig,trace="none",col=mycols,main="The 100 most significant genes",ColSideColors=column_cols)
dev.off()



##### PCA ####
pca=prcomp(t(logCPM))
summary(pca)


data <- data.frame(pca$x, Sex = metadata$Sex, Diagnosis = metadata$diagnosis)

# Create the plot
ggplot(data, aes(x = PC1, y = PC2, shape = Sex, color = Diagnosis)) +
  geom_point(size = 4) +
  labs(x = "PCA 1", y = "PCA 2") + 
  scale_shape_manual(values = c(15, 8))
ggsave("./Figures/pca_12_diagnosis.png", height = 15, width = 12, units = "cm")

ggplot(data, aes(x = PC1, y = PC3, shape = Sex, color = Diagnosis)) +
  geom_point(size = 4) +
  labs(x = "PCA 1", y = "PCA 3") + 
  scale_shape_manual(values = c(15, 8))
ggsave("./Figures/pca_13_diagnosis.png", height = 15, width = 12, units = "cm")

ggplot(data, aes(x = PC2, y = PC3, shape = Sex, color = Diagnosis)) +
  geom_point(size = 4) +
  labs(x = "PCA 2", y = "PCA 3") + 
  scale_shape_manual(values = c(15, 8))
ggsave("./Figures/pca_23_diagnosis.png", height = 15, width = 12, units = "cm")


