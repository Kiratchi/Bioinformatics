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

# Plot
logCPM_long <- pivot_longer(logCPM, cols = everything(), names_to = "person_id", values_to = "expression")
head(logCPM_long)
ggplot(data=logCPM_long,aes(x=person_id, y = expression)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Expression by Person", 
       x = "Person ID", 
       y = "Expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
