#create function medge to calculate edges
mnet <- function(x){#x:connectivity matrix
  p <- ncol(x[[1]])
  
  #subset edges for certain brain regions
  x1 <- lapply(1:length(x), function(y)t(x[[y]])[lower.tri(x[[y]])])#convert off-diagonal part of each connectivity matrix into a vector, do this for each frequency
  aa <- combn(1:p, 2)#combn(x, m): Generate all combinations of the elements of x taken m at a time
  
  #generate the number of edges for VN1
  y2 <- ifelse(aa[1,] %in% 55:56 & aa[2,] %in% 55:56, 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.vn1 <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))#sum all edges that both brain region below to this network, here VN1
  names(edge.vn1) <- 1:length(x)
  
  #generate the number of edges for VN2
  y2 <- ifelse(aa[1,] %in% 51:52 & aa[2,] %in% 51:52, 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.vn2 <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.vn2) <- 1:length(x)
  
  #generate the number of edges for AN
  y2 <- ifelse(aa[1,] %in% 81:82 & aa[2,] %in% 81:82, 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.an <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.an) <- 1:length(x)
  
  #generate the number of edges for SMN
  y2 <- ifelse(aa[1,] %in% c(1,7,57,70) & aa[2,] %in% c(1,7,57,70), 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.smn <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.smn) <- 1:length(x)
  
  #generate the number of edges for SRN
  y2 <- ifelse(aa[1,] %in% c(7,8,23,24,31,32,35,36) & aa[2,] %in% c(7,8,23,24,31,32,35,36), 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.srn <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.srn) <- 1:length(x)
  
  #generate the number of edges for DAN
  y2 <- ifelse(aa[1,] %in% c(3,7,8,61,62) & aa[2,] %in% c(3,7,8,61,62), 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.dan <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.dan) <- 1:length(x)
  
  #generate the number of edges for VAN
  y2 <- ifelse(aa[1,] %in% c(4,8,60,62,86,90) & aa[2,] %in% c(4,8,60,62,86,90), 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.van <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.van) <- 1:length(x)
  
  #generate the number of edges for DMN
  y2 <- ifelse(aa[1,] %in% c(25,26,35,36,37,38,61,62,89,90) & aa[2,] %in% c(25,26,35,36,37,38,61,62,89,90) , 1, 0)
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))
  edge.dmn <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==1]==1))
  names(edge.dmn) <- 1:length(x)
  
  result <- list("VN1"=edge.vn1, "VN2"=edge.vn2, "AN"=edge.an, "SMN"=edge.smn, "SRN"=edge.srn, "DAN"=edge.dan, "VAN"=edge.van, "DMN"=edge.dmn)#https://www.geeksforgeeks.org/r-lists/
  return(result)
}
