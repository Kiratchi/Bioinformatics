library(ggplot2)

load("./Data/genome1.rdata")
load("./Data/genome2.rdata")
load("./Data/genome3.rdata")

#class(genome1)
#head(genome1)
#dim(genome1)


#class(reference)
#head(reference)
#length(reference)

genome1 = genome1[1:1000,]
reference = reference[1:1000]
# genome1[1,2:5] = c(0,0,0,0)

coverage <- rowSums(genome1[, 2:5])
print(min(coverage))

#head(coverage)
#summary(coverage)


coverage_x_axis = seq(1,1000)
coverage_data = data.frame(Index = coverage_x_axis, Coverage = coverage[coverage_x_axis])
ggplot(data=coverage_data, aes(x = Index, y = Coverage)) + geom_point() + labs(x="Position")
# ggplot(data=coverage_data, aes(x = Index, y = Coverage)) + geom_point()


genome1.length = nrow(genome1)
reference_indices <- match(reference, colnames(genome1)[-1])
matches = genome1[cbind(1:genome1.length, reference_indices+1)]
mismatch_values = coverage - matches
mismatch_indices = which(mismatch_values != 0)
mismatches <- data.frame(
  Index = mismatch_indices,
  Mismatch = mismatch_values[mismatch_indices],
  Proportion = mismatch_values[mismatch_indices] / coverage[mismatch_indices]
)

dim(mismatches)
#ggplot(data = mismatches, aes(x = Index, y = Proportion)) + geom_point()
#ggsave("./Figures/proportion_missing_g1.png")


#Binomial test
mismatches$Pvalue <- numeric(nrow(mismatches))


Binom_calc = function(x,n){
  bin_result <- binom.test(x,n,p=0.1,alternative = "greater")
  return(bin_result$p.value)
}


for (i in 1:nrow(mismatches)){
    #bin_result <- binom.test(mismatches$Mismatch[i],coverage[mismatches$Index[i]],p=0.1,alternative = "greater")
    mismatches$Pvalue[i] <- Binom_calc(mismatches$Mismatch[i],coverage[mismatches$Index[i]])
    #mismatches$Pvalue[i] <- bin_result$p.value
}

plot(mismatches$Index,mismatches$Pvalue)

mismatches <- mismatches[order(mismatches$Pvalue),]

######################


genome1.length = nrow(genome1)
matches = vector(length = genome1.length)

mismatches = data.frame(Index=integer(),Mismatch = integer(),Proportion=numeric())
#vector(length = genome1.length)

for (pos in 1:genome1.length){
  matches[pos]=genome1[pos,reference[pos]]
  current_mismatch=coverage[pos]-matches[pos]
  
  if ( current_mismatch != 0){
    mismatches <- rbind(mismatches,data.frame(Index=pos, Mismatch = current_mismatch, Proportion = current_mismatch/coverage[pos]))
  }
  
}





