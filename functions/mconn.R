library(astsa)

#create function mconn to calculate connectivity matrix
mconn <-  function(x, alpha=0.05, s, tt){
  #x: time series matrix, s: segments, t: partition time point 
  p <- ncol(x)
  n <- nrow(x)
  m <- matrix(NA, ncol = p, nrow = p)
  
  ff <-  function(x){
    m <- floor(sqrt(dim(x)[1]))
    mspec  <-  mvspec(x, spans=c(2*m+1,2*m+1), kernel="daniell", plot=F)
    Fstat <- apply(mspec$coh, 2, function(y)(y/(1-y))*(2*m))
    pval <- 1-pf(Fstat, 2, 2*(2*m+1) - 2)
    pval1 <- sapply(1:ncol(pval), function(x) p.adjust(pval[,x], method = "bonferroni", n=2*length(pval[,x])))
    Conn <- ifelse(pval1 < alpha, 1, 0)
    nn <- nrow(Conn)
    
    fun <- function(fx){
      # Create a vector with appropriate length
      vec <- Conn[fx,]  # Length of upper diagonal (excluding main diagonal)
      # Get indices of upper diagonal (excluding main diagonal)
      upper_idx <- which(col(m11) > row(m11), arr.ind = TRUE)
      # Fill in column-wise
      m11[upper_idx[order(upper_idx[,2], upper_idx[,1]), ]] <- vec
      # Copy upper triangle to lower triangle to make it symmetric
      m11[lower.tri(m11)] <- t(m11)[lower.tri(m11)]
      diag(m11) <- 0
      return(m11)
    }
    
    result <- lapply(1:nn, fun)
    return(result)
  }
  
  
  if(s==0){
    ff(x)
  }
  
  else if(s==1){
    t1 <- tt+1
    m1 <- ff(x[1:tt,])
    m2 <- ff(x[t1:n,])
    l1 <- length(1:tt)
    l2 <- length(t1:n)
    freq1 <- (1:floor(l1/2))/l1
    freq2 <- (1:floor(l2/2))/l2
    comp <- which.min(c(length(freq1),length(freq2)))
    
    if(comp==1){
      combos <- abs(outer(freq2, freq1, "-"))
      freq.min <- apply(combos,2,which.min)
      m3 <- m2[c(freq.min)]
      f1 <- function(x){
        m <- matrix(1, ncol = p, nrow = p)
        m[!(m1[[x]]==1 & m3[[x]]==1)] <- 0
        return(m)
      }
      result <- lapply(1:length(freq1), f1)
      return(result)
    }
    
    else if(comp==2){
      combos <- abs(outer(freq1, freq2, "-"))
      freq.min <- apply(combos,2,which.min)
      m3 <- m1[c(freq.min)]
      f2 <- function(x){
        m <- matrix(1, ncol = p, nrow = p)
        m[!(m2[[x]]==1 & m3[[x]]==1)] <- 0
        return(m)
      }
      result <- lapply(1:length(freq2), f2)
      return(result)
    }
  }
  else if(s==2){
    if(s==2 & length(tt)==1){
      print("Wrong input: if s=2, tt should be a vector with 2 elements")
    }else{
    t1 <- tt[1]+1
    t2 <- tt[2]+1
    m1 <- ff(x[1:tt[1],])
    m2 <- ff(x[t1:tt[2],])
    m3 <- ff(x[t2:n,])
    l1 <- length(1:tt[1])
    l2 <- length(t1:tt[2])
    l3 <- length(t2:n)
    freq1 <- (1:floor(l1/2))/l1
    freq2 <- (1:floor(l2/2))/l2
    freq3 <- (1:floor(l3/2))/l3
    comp <- which.min(c(length(freq1),length(freq2),length(freq3)))
    
    if(comp==1){
      combo1 <- abs(outer(freq2, freq1, "-"))
      combo2 <- abs(outer(freq3, freq1, "-"))
      freq.min1 <- apply(combo1,2,which.min)
      freq.min2 <- apply(combo2,2,which.min)
      m2.new <- m2[c(freq.min1)]
      m3.new <- m3[c(freq.min2)]
      f3 <- function(x){
        m <- matrix(1, ncol = p, nrow = p)
        m[!(m1[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1)] <- 0
        return(m)
      }
      result <- lapply(1:length(freq1), f3)
      return(result)
    }
    
    else if(comp==2){
      combo1 <- abs(outer(freq1, freq2, "-"))
      combo2 <- abs(outer(freq3, freq2, "-"))
      freq.min1 <- apply(combo1,2,which.min)
      freq.min2 <- apply(combo2,2,which.min)
      m1.new <- m1[c(freq.min1)]
      m3.new <- m3[c(freq.min2)]
      f4 <- function(x){
        m <- matrix(1, ncol = p, nrow = p)
        m[!(m1.new[[x]]==1 & m2[[x]]==1 & m3.new[[x]]==1)] <- 0
        return(m)
      }
      result <- lapply(1:length(freq2), f4)
      return(result)
    }
    
    else if(comp==3){
      combo1 <- abs(outer(freq1, freq3, "-"))
      combo2 <- abs(outer(freq2, freq3, "-"))
      freq.min1 <- apply(combo1,2,which.min)
      freq.min2 <- apply(combo2,2,which.min)
      m1.new <- m1[c(freq.min1)]
      m2.new <- m2[c(freq.min2)]
      f5 <- function(x){
        m <- matrix(1, ncol = p, nrow = p)
        m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3[[x]]==1)] <- 0
        return(m)
      }
      result <- lapply(1:length(freq3), f5)
      return(result)
    }
    }
  }
}
