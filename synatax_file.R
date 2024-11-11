O # Generic object
M # Generic matrix
v # Generic vector
D # Generic data
df #Generic data frame




#### IO ####
print()
write.table(data_file,file="data.txt",sep="\t") #Save the Data file with tab separation 
read.table("https://tinyurl.com/y8av4xfn",header=TRUE, sep="\t")

#### Objects ####
# All variables are objects
ls() #List all objects
rm(x) #deletes x
rm(list=ls()) #removes all objects

head(O)  #Shows the beginning of most objects
tail(O)  #Shows the end of most objects
summary(O) #Shows a summary 
dim(o) #Shows the number of rows, columns etc in a vector

#### LOOPS ####
for(i in 1:trials){
  # write your code here
}



#### VECTORS ####
# Vectors contain only one data type

c(2,3,4) #Combine elements into vector
numeric(3) #Creates 3 element zero vector

seq(1, 100, 9) # Create vector 1 to 100 with step length 9
1:100 # =seq(1, 1, 100)

rep(1:3, t=3, e=4) # take the input vector 1:3 and take each element 3 times
# And the entire vector 4 times


x[1:2] #take value 1:2 from x
x[-1:2] #take all values from x except 1:2



#### MATRICES ####
# 2d matrices with only one data type

matrix(,nrow=2,ncol=3) # matrix with 2 rows and 3 columns

A %*% B #Matrix multiplication
t() #Transpose matrix


apply(M,1,sum) # Does the function sum over all rows of matrix z 
apply(M,1,sum) # Does the function sum over all columns of matrix z 



#### DATA FRAMES ####
#

data.frame(name_a = v1, name_b = v2) #Creates a data frame

rownames() # shows row names
colnames() #Shows the column names
rbind() # Add new row

data_frame$name_a #Accesses column name_a
data_frame[,1] #Access column 1

#### DATA MANIPULTAION ####
cbind(M,v) #Add column vector v to end of M
order(v,decreasing=TRUE) #Sort vector v in decreasing and returns the order they should have



#### PSEUDO-RANDOM NUMBERS ####
set.seed(123) #Sets the seed to 123



#### STATISTICS ####
mean(D) 
median(D) 
var(D)  
sd(D)

t.test(D,D) #Does a t-test
var.test(D,D) #Does a F-test
wilcox.test(D,D) #Does a wilcox Rank sum test




#### LINEAR REGRESSION ####
lm(y~x,data=df) #Creates a model E(y) = B0 + B1x from x & y in df
lm(y~x-1,data=df) #Same but B0=0
summary(mymodel) #Can summaries mymodel
resid(mymodel) #List of residuals
predict(mymodel, newdata) #Show model prediction for some data


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

abline(mymodel) #Plot fitted regression line

hist(D) #Creates a histogram
# break : number of bars

barplot(D) #Creates a barplot


pdf(file="X.pdf") # Opens a new pdf named "X"
dev.off() #Closes the most recent open graphic 
graphics.off #Cloese all graphics

windows()
layout(matrix(1:2,nrow=1, ncol=2)) # Creates a subfigure where the 2 subsequent
                                   # figures will be saved



