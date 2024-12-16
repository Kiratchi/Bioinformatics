library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(gplots)

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

#Rarefaction

RareFaction = function(OTU_name,count,depth){
  reads <- rep(OTU_name, times = count)
  reads.sample<-sample(reads,size=depth, replace=FALSE)
counts.sample<-as.data.frame(table(reads.sample))
return(counts.sample)
  }

RareFaction(rownames(s_counts), s_counts$HC2, 10)
