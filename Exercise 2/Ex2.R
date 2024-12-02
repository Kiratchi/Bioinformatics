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

#genome1 <- genome1[1:1000,]
#genome2 <- genome2[1:1000,]
#genome3 <- genome3[1:1000,]
#reference <- reference[1:1000]



Plot_Coverage = function(coverage,genome_number){
  coverage_x_axis <- seq(1,1000)
  coverage_data <- data.frame(Index = coverage_x_axis, Coverage = coverage[coverage_x_axis])
  ggplot(data=coverage_data, aes(x = Index, y = Coverage)) + 
    geom_area() + 
    labs(x="Index")
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
  ggplot(data = mismatches, aes(x = Index, y = Proportion)) + 
    geom_point() +
    xlim(0,4641652)
  filename <- sprintf("./Figures/proportion_missing%d.png", genome_number)
  ggsave(filename)
  
  #Binomial test
  mismatches$Pvalue <- numeric(nrow(mismatches))
  for (i in 1:nrow(mismatches)){
    mismatches$Pvalue[i] <- Binom_calc(mismatches$Mismatch[i],mismatches$Coverage[i])
    
  }
  
  #Plot P-values of mismatches
  ggplot(data=mismatches, aes(x = Index, y = log(Pvalue))) + 
    geom_point(aes(color = log(Pvalue) < -10)) +  
    scale_color_manual(values = c("black", "darkred")) +
    guides(color = "none") +
    geom_hline(yintercept = -10, linetype = "dashed", color = "darkred") +
    xlim(0,4641652) +
    ylim(-43,0)
  filename <- sprintf("./Figures/pvalue_plot%d.png", genome_number)
  ggsave(filename)
  
  mismatches <- mismatches[order(mismatches$Pvalue),]
  print(sprintf("Genome %d", genome_number))
  print(mismatches[1:10,])
  
}


Function_Doing_Stuffs(genome1,reference,1)
Function_Doing_Stuffs(genome2,reference,2)
Function_Doing_Stuffs(genome3,reference,3)




###########










