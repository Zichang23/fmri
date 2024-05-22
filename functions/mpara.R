#create function mpara to calculate edges between two corresponding hemisphere regions (L and R)
mpara <- function(x){#x:connectivity matrix
  p <- ncol(x[[1]])
  
  #subset edges for certain brain regions
  x1 <- lapply(1:length(x), function(y)t(x[[y]])[lower.tri(x[[y]])])#convert off-diagonal part of each connectivity matrix into a vector, do this for each frequency
  aa <- combn(1:p, 2)
  
  #count number of edges for each pair of hemisphere regions
  func0 <- function(z){
    y2 <- ifelse(aa[1,] %in% z:(z+1) & aa[2,] %in% z:(z+1) & aa[1,]!=aa[2,], 1, 0)
    x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
    edge0 <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
    names(edge0) <- 1:length(x)
    return(edge0)
  }
  
  result0 <- lapply(seq(1,107,by=2), function(z)func0(z))
  result <- Reduce('+', result0)
  #https://datascienceparichay.com/article/r-list-sum/
  #sum up number of edges for each frequency
  return(result)
}
