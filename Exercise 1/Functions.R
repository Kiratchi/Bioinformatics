
sum_of_integers = function(n){
  return((n+1)*n/2) 
  
}



fibonacci = function(n){
  if (n == 0) {
    return(0)}
  if (n == 1) {
    return(1)}
  if (n == 2) {
    return(2)}
  
  fib = numeric(n)
  fib[1] = 0
  fib[2] = 1
  for (i in 3:n){
    fib[i] = fib[i-1] + fib[i-2]
  }
  return(fib)
}

  
Euclidean_distance = function(x,y){
  sum2 = 0
  for (i in 1:length(x)){
    sum2 = sum2 + (x[i]-y[i])^2
  }
  return(sqrt(sum2))
}

fibonacci(15)
sum_of_integers(4)
Euclidean_distance(c(2,5,6,7,3),c(5,2,8,4,4))

