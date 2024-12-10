library(ggplot2)
library(tidyr)

x = read.table("./Data/counts_matrix.txt")
metadata = read.table("./Data/metadata.txt", sep="\t", header=TRUE)

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
CPM <- sweep(x.filtered, 2, 10^6/colSums(x.filtered), "*") 
logCPM = log(CPM+1)
logx.filtered <- log(x.filtered+1)


# Box-plot normalized data
logCPM_long <- pivot_longer(logCPM, cols = everything(), names_to = "person_id", values_to = "expression")
head(logCPM_long)
ggplot(data=logCPM_long, aes(x=person_id, y = expression)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Expression by Person - Normalized", 
       x = "Person ID", 
       y = "log(Expression)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("./Figures/boxplot_normalized.png")

# Box-plot non-normalized data = x.filtered
logx.filtered_long <- pivot_longer(logx.filtered, cols = everything(), names_to = "person_id", values_to = "expression")
head(logx.filtered_long)
ggplot(data=logx.filtered_long, aes(x=person_id, y = expression)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Expression by Person - Non Normalized", 
       x = "Person ID", 
       y = "log(Expression)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("./Figures/boxplot_non-normalized.png")

# Scatter plot sample 1 and 41
ggplot(data=logCPM , aes(x=logCPM[,1] , y=logCPM[,41] )) + geom_point(size = 0.5) +
  labs(title = "Scatter plot of sample 41 against sample 1", 
       x = "log(expression) for Sample 1", 
       y = "log(expression) for Sample 41" ) +
  geom_abline(aes(intercept = 0, slope = 1), color = "orange", linewidth = 1.5)
ggsave("./Figures/scatterplot_41_1.png")
