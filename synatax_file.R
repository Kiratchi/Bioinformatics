O # Generic object
M # Generic matrix
v # Generic vector
D # Generic data


head(O)  #Shows the beginning of most objects
tail(O)  #Shows the end of most objects
summary(O) #Shows a summary 


M[,1]   #Take the first column of matrix M

#### MATRICES ####
rownames() # shows row names
colnames() #Shoes the column names

apply(M,1,sum) # Does the function sum over all rows of matrix z 
apply(M,1,sum) # Does the function sum over all columns of matrix z 

#### DATA MANIPULTAION ####
cbind(M,v) #Add column vector v to end of M
order(v,decreasing=TRUE) #Sort vector v in decreasing and returns the order they should have

#### STATISTICS ####
mean(D) 
median(D) 
var(D)  
sd(D)

t.test(D,D) #Does a t-test
var.test(D,D) #Does a F-test
wilcox.test(D,D) #Does a wilcox Rank sum test


#### PLOTTING ####

plot(d1, d2) # Creates a scatter plot of d2 over d1
# pch : which char character to use
# main : title text
# xlab : x lable text
# ylab : y lable text
# xlim : x min and x max in vector
# ylim : y min and y max in vector
# cex : size of points  "cex.X" to edit "X"
# col : color of points "col.X" to edit "X"

hist(D) #Creates a histogram
# break : number of bars

barplot(D) #Creates a barplot


pdf(file="X.pdf") # Opens a new pdf named "X"
dev.off() #Closes the most recent open graphic 
graphics.off #Cloese all graphics

windows()
layout(matrix(1:2,nrow=1, ncol=2)) # Creates a subfigure where the 2 subsequent
                                   # figures will be saved
