library(ggplot2)

load("./Data/genome1.rdata")
load("./Data/genome2.rdata")
load("./Data/genome3.rdata")

class(genome1)
head(genome1)
dim(genome1)


class(reference)
head(reference)
length(reference)

#genome1 = genome1[1:1000,]
#reference = reference[1:1000]

# faster then apply since vectorized
coverage <- rowSums(genome1[, 2:5])

head(coverage)
summary(coverage)


coverage_x_axis = seq(1000,2000)
coverage_data = data.frame(Index = coverage_x_axis, Coverage = coverage[coverage_x_axis])
ggplot(data=coverage_data, aes(x = Index, y = Coverage)) + geom_area() + labs(x="Position")
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


ggplot(data = mismatches, aes(x = Index, y = Proportion)) + geom_point()
dim(mismatches)
ggsave("proportion_missing_g1.png", path = "./Figures")

ggsave("./Figures/proportion_missing_g1.png")


######################


genome1.length = nrow(genome1)
matches = vector(length = genome1.length)

mismatches = data.frame(Index=integer(),Mismatch = integer(),Proportion=numeric())
#vector(length = genome1.length)

for (pos in 1:genome1.length){
  matches[pos]=genome1[pos,reference[pos]]
  current_mismatch=coverage[pos]-matches[pos]
  
  if ( current_mismatch != 0){
    mismatches <- rbind(mismatches,data.frame(Index=pos, Mismatch = current_mismatch, Proportion=current_mismatch/coverage[pos]))
  }
  
}





