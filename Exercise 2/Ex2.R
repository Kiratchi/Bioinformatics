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

genome1 <- genome1[1:1000,]
genome2 <- genome2[1:1000,]
genome3 <- genome3[1:1000,]
reference <- reference[1:1000]


Function_Doing_Stuffs(genome1,reference,1)
Function_Doing_Stuffs(genome2,reference,2)
Function_Doing_Stuffs(genome3,reference,3)

Plot_Coverage = function(coverage,genome_number){
  coverage_x_axis <- seq(1,1000)
  coverage_data <- data.frame(Index = coverage_x_axis, Coverage = coverage[coverage_x_axis])
  ggplot(data=coverage_data, aes(x = Index, y = Coverage)) + geom_area() + labs(x="Position")
  filename <- sprintf("./Figures/coverage%d.png", genome_number)
  ggsave(filename)
}

Binom_calc = function(x,n){
  bin_result <- binom.test(x,n,p=0.1,alternative = "greater")
  return(bin_result$p.value)
}


Function_Doing_Stuffs = function(genomex,reference,genome_number){
  coverage <- rowSums(genomex[, 2:5]) 
  Plot_Coverage(coverage,genome_number)
  
  genomex.length <- nrow(genomex)
  
  reference_indices <- match(reference, colnames(genomex)[-1])
  matches = genomex[cbind(1:genomex.length, reference_indices+1)] 
  
  mismatch_values <- coverage - matches
  mismatch_indices <- which(mismatch_values != 0)
  
  mismatches <- data.frame(
    Index = mismatch_indices,
    Mismatch = mismatch_values[mismatch_indices],
    Proportion = mismatch_values[mismatch_indices] / coverage[mismatch_indices],
    Coverage = coverage[mismatch_indices]
  )
  
  # Plot proportion of mismatches
  ggplot(data = mismatches, aes(x = Index, y = Proportion)) + geom_point()
  filename <- sprintf("./Figures/proportion_missing%d.png", genome_number)
  ggsave(filename)
  
  #Binomial test
  mismatches$Pvalue <- numeric(nrow(mismatches))
  for (i in 1:nrow(mismatches)){
    mismatches$Pvalue[i] <- Binom_calc(mismatches$Mismatch[i],mismatches$Coverage[i])
    
  }
  
  #Plot P-values of mismatches
  ggplot(data=mismatches, aes(x = Index, y = Pvalue)) + geom_point() 
  filename <- sprintf("./Figures/pvalue_plot%d.png", genome_number)
  ggsave(filename)
  
  mismatches <- mismatches[order(mismatches$Pvalue),]
  
  
}








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





