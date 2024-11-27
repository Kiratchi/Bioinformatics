#rm(list=ls())

statistics_of_Bwt = function(matrix){
  mean = mean(matrix)
  median = median(matrix)
  var = var(matrix)
  sd = sqrt(var)
  
  return(c(mean, median, sd, var))
}



library(MASS)
data(cats)


head(cats)
summary(cats)
rownames(summary)



statistics = statistics_of_Bwt(cats[,2])
Bwt_mean = statistics[1]
Bwt_median = statistics[2]
Bwt_sd = statistics[3]
Bwt_var = statistics[4]

Bwt_Hwt_cor=cor(cats[,2],cats[,3])

female = cats[,1]=="F"
cats.female = cats[female,]
cats.male = cats[!female,]

statistics = statistics_of_Bwt(cats.female[,2])
Bwt_mean.female = statistics[1]
Bwt_median.female = statistics[2]
Bwt_sd.female = statistics[3]
Bwt_var.female = statistics[4]

statistics = statistics_of_Bwt(cats.male[,2])
Bwt_mean.male = statistics[1]
Bwt_median.male = statistics[2]
Bwt_sd.male = statistics[3]
Bwt_var.male = statistics[4]

var.test(cats.female[,2],cats.male[,2])
# There is a significant difference in variance of Bwt between female and male cats

t.test(cats.female[,2],cats.male[,2], var.equal = FALSE)
# There is a significant difference in mean Bwt between female and male cats

wilcox.test(cats.female[,2],cats.male[,2])


#### Plotting ####

plot(cats[,2], cats[,3], col="blue",pch=10, main="Cat body-hear weight",
     xlab="Body weight", ylab="Heart weight", cex=1.5, cex.main=1.5, 
     xlim=c(2,2.5),ylim=c(6,15))

hist(cats[,2], breaks=20, main="Cat Heart Weight", xlab = "Heart weight")

pdf(file="cats.pdf",width = 9, height = 8)
barplot(cats[1:10,3], main="10 Cats Heart Weight", ylab = "Heart weight")
dev.off()
