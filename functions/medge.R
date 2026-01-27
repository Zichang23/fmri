#create function medge to calculate edges (assume connectivity matrices are symmetric)
medge <- function(x){#x:connectivity matrix
  p <- ncol(x[[1]])
  
  #generate the total number of edges
  edge.all <- sapply(1:length(x), function(y)sum(x[[y]]==1))/2
  names(edge.all) <- 1:length(x)#change label for each element
  
  #subset edges for certain brain regions
  x1 <- lapply(1:length(x), function(y)x[[y]][upper.tri(x[[y]])])#convert off-diagnal part of each connectivity matrix into a vector, do this for each frequency, https://stackoverflow.com/questions/70267887/extracting-upper-off-diagonal-elements-of-square-matrix-in-row-order
  aa <- combn(1:p, 2)#https://stat.ethz.ch/R-manual/R-devel/library/utils/html/combn.html
  y1 <- ifelse(aa %% 2 == 0, 1, 0)
  # y2 <- ifelse(aa[1,] %in% 109:116 & aa[2,] %in% 109:116, 5, ifelse(aa %in% 109:116, 4,0))
  y2 <- ifelse(aa[1,] %in% 109:116 & aa[2,] %in% 109:116, 5, ifelse((aa[1,] %in% 109:116 | (aa[2,] %in% 109:116)), 4, ifelse(y1[1,]==0 & y1[2,]==0, 2, ifelse(y1[1,]==1 & y1[2,]==1,3,4))))#2:edges in Left hemisphere, 3:edges in right hemisphere, 4:edges in inter-hemisphere & hemisphere-vermis, 5:edges in vermis
  x2 <- lapply(1:length(x), function(y)rbind(y2,x1[[y]]))#create matrix with first rows as 'label' for each frequency
  
  #generate the number of edges for left hemisphere
  edge.l <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==2]==1))#https://stackoverflow.com/questions/1923273/counting-the-number-of-elements-with-the-values-of-x-in-a-vector
  #https://ademos.people.uic.edu/Chapter4.html
  names(edge.l) <- 1:length(x)
  
  #generate the number of edges for right hemisphere
  edge.r <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==3]==1))
  names(edge.r) <- 1:length(x)
  
  #generate the number of edges for inter-hemisphere and hemisphere-vermis
  edge.inter <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==4]==1))
  names(edge.inter) <- 1:length(x)
  
  #generate the number of edges for vermis
  edge.v <- sapply(1:length(x), function(y)sum(x2[[y]][2,x2[[y]][1,]==5]==1))
  names(edge.v) <- 1:length(x)
  result <- list("Total"=edge.all, "Left"=edge.l, "Right"=edge.r, "Inter"=edge.inter, "Vermis"=edge.v)#https://www.geeksforgeeks.org/r-lists/
  return(result)
}
