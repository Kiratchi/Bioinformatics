library(MASS)
data(mammals)

ls(mammals)
head(mammals)
rownames(mammals)

body_mean = mean(mammals[,1])
body_median = median(mammals[,1])
body_variance = var(mammals[,1])
body_std = sqrt(body_variance)

brain_mean = mean(mammals[,2])
brain_median = median(mammals[,2])
brain_variance = var(mammals[,2])
brain_std = sqrt(body_variance)

brain_to_body_ratio=(mammals[,2]/1000/mammals[,1])
mammals2=cbind(mammals,brain_to_body_ratio)
colnames(mammals2) =c(colnames((mammals)), "Brain/Body")
mammals.order=order(mammals2[,3], decreasing=TRUE)
mammals2.sorted=mammals2[mammals.order,]

# African elephant lowest Ground squirrel highest
