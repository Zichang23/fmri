---
title: "Example"
# author: "Zichang Xiang"
# date: "2024-05-07"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## Introduction

In this section, we will show you how to use our functions to calculate the number of edge for certain combination of frequencies and networks. 

### Data

The data we will use is `Caltech_0051475_rois_aal.1D`, which is a rs-fMRI dataset for a control subject from the Autism Brain Image Exchange (ABIDE) repository. For more information about this repository, please visit [ABIDE Preprocessed](preprocessed-connectomes-project.org/abide/).


### Functions

Here we will use two functions, `mconn2` and `mnet`.

`mconn2` is a function to calculate connectivity matrix across certain frequencies.

`mnet` is a function to calculate number of edges for all brain networks.

Below are all the arguments in `mconn2`. `x` is the dataset in format ".D". alpha is the significance level, the default value is 0.05. `s` is the number of stationary time series within the dataset `x`. `tt` is the time point for each stationary time series segment. `freq0` are the specific frequencies we would like to focus.

```{}
mconn2(x, alpha = 0.05, s, tt, freq0)
```

```{r,echo=FALSE}
mconn2 <-  function(x, alpha=0.05, s, tt, freq0){
  #x: time series matrix (no need to standardize), s: segments, t: partition time point , freq0: specific frequencies
  p <- ncol(x)
  n <- nrow(x)
  m <- matrix(NA, ncol = p, nrow = p)
  
  ff <-  function(x){
    L <- floor(sqrt(dim(x)[1]))
    mspec  <-  mvspec(x, spans=L, fast =F, kernel="daniell", plot=F)
    Fstat <- apply(mspec$coh, 2, function(y)(y/(1-y))*(L-1))
    Fq <-  qf(1-alpha, 2, 2*L-2)
    Conn <- ifelse(Fstat > Fq, 1, 0)
    nn <- nrow(Conn)
    fun <- function(x){
      m[lower.tri(m)] <- Conn[x,] 
      m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
      diag(m) <- 0
      return(m)
    }
    result <- lapply(1:nn, fun)
    return(result)
  }
  if(s==0){
    # if(s==0 && tt==0){stop("Wrong input: if s=0, please type tt=0")}else{ff(x)}
    result <- ff(x)
    freq1 <- (1:floor(n/2))/n
    combo0 <- abs(outer(freq1, freq0, "-"))
    freq.min0 <- apply(combo0,2,which.min)
    m0 <- result[c(freq.min0)]
    return(m0)
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
      combo0 <- abs(outer(freq1, freq0, "-"))
      freq.min0 <- apply(combo0,2,which.min)
      m0 <- result[c(freq.min0)]
      return(m0)
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
      combo0 <- abs(outer(freq2, freq0, "-"))
      freq.min0 <- apply(combo0,2,which.min)
      m0 <- result[c(freq.min0)]
      return(m0)
    }
  }
  else if(s==2){
    if(s==2 & length(tt)!=2){
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
        combo0 <- abs(outer(freq1, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
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
        combo0 <- abs(outer(freq2, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
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
        combo0 <- abs(outer(freq3, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
    }
    
  }
  else if(s==3){
    if(s==3 & length(tt)!=3){
      print("Wrong input: if s=3, tt should be a vector with 3 elements")
    }else{
      t1 <- tt[1]+1
      t2 <- tt[2]+1
      t3 <- tt[3]+1
      m1 <- ff(x[1:tt[1],])
      m2 <- ff(x[t1:tt[2],])
      m3 <- ff(x[t2:tt[3],])
      m4 <- ff(x[t3:n,])
      l1 <- length(1:tt[1])
      l2 <- length(t1:tt[2])
      l3 <- length(t2:tt[3])
      l4 <- length(t3:n)
      freq1 <- (1:floor(l1/2))/l1
      freq2 <- (1:floor(l2/2))/l2
      freq3 <- (1:floor(l3/2))/l3
      freq4 <- (1:floor(l4/2))/l4
      comp <- which.min(c(length(freq1),length(freq2),length(freq3),length(freq4)))
      if(comp==1){
        combo1 <- abs(outer(freq2, freq1, "-"))
        combo2 <- abs(outer(freq3, freq1, "-"))
        combo3 <- abs(outer(freq4, freq1, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        
        m2.new <- m2[c(freq.min1)]
        m3.new <- m3[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        
        f6 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1 & m4.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq1), f6)
        combo0 <- abs(outer(freq1, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==2){
        combo1 <- abs(outer(freq1, freq2, "-"))
        combo2 <- abs(outer(freq3, freq2, "-"))
        combo3 <- abs(outer(freq4, freq2, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m3.new <- m3[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        
        f7 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2[[x]]==1 & m3.new[[x]]==1 & m4.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq2), f7)
        combo0 <- abs(outer(freq2, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==3){
        combo1 <- abs(outer(freq1, freq3, "-"))
        combo2 <- abs(outer(freq2, freq3, "-"))
        combo3 <- abs(outer(freq4, freq3, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m2.new <- m2[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        
        f8 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3[[x]]==1 & m4.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq3), f8)
        combo0 <- abs(outer(freq3, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==4){
        combo1 <- abs(outer(freq1, freq4, "-"))
        combo2 <- abs(outer(freq2, freq4, "-"))
        combo3 <- abs(outer(freq3, freq4, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m2.new <- m2[c(freq.min2)]
        m3.new <- m3[c(freq.min3)]
        
        f9 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1 & m4[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq4), f9)
        combo0 <- abs(outer(freq4, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
    }
  } 
  else if(s==4){
    if(s==4 & length(tt)!=4){
      print("Wrong input: if s=4, tt should be a vector with 4 elements")
    }else{
      t1 <- tt[1]+1
      t2 <- tt[2]+1
      t3 <- tt[3]+1
      t4 <- tt[4]+1
      
      m1 <- ff(x[1:tt[1],])
      m2 <- ff(x[t1:tt[2],])
      m3 <- ff(x[t2:tt[3],])
      m4 <- ff(x[t3:tt[4],])
      m5 <- ff(x[t4:n,])
      
      l1 <- length(1:tt[1])
      l2 <- length(t1:tt[2])
      l3 <- length(t2:tt[3])
      l4 <- length(t3:tt[4])
      l5 <- length(t4:n)
      
      freq1 <- (1:floor(l1/2))/l1
      freq2 <- (1:floor(l2/2))/l2
      freq3 <- (1:floor(l3/2))/l3
      freq4 <- (1:floor(l4/2))/l4
      freq5 <- (1:floor(l5/2))/l5
      
      comp <- which.min(c(length(freq1),length(freq2),length(freq3),length(freq4),length(freq5)))
      if(comp==1){
        combo1 <- abs(outer(freq2, freq1, "-"))
        combo2 <- abs(outer(freq3, freq1, "-"))
        combo3 <- abs(outer(freq4, freq1, "-"))
        combo4 <- abs(outer(freq5, freq1, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        freq.min4 <- apply(combo4,2,which.min)
        
        m2.new <- m2[c(freq.min1)]
        m3.new <- m3[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        m5.new <- m5[c(freq.min4)]
        
        f10 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1 & m4.new[[x]]==1 & m5.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq1), f10)
        combo0 <- abs(outer(freq1, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==2){
        combo1 <- abs(outer(freq1, freq2, "-"))
        combo2 <- abs(outer(freq3, freq2, "-"))
        combo3 <- abs(outer(freq4, freq2, "-"))
        combo4 <- abs(outer(freq5, freq2, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        freq.min4 <- apply(combo4,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m3.new <- m3[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        m5.new <- m5[c(freq.min4)]
        
        f11 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2[[x]]==1 & m3.new[[x]]==1 & m4.new[[x]]==1 & m5.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq2), f11)
        combo0 <- abs(outer(freq2, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==3){
        combo1 <- abs(outer(freq1, freq3, "-"))
        combo2 <- abs(outer(freq2, freq3, "-"))
        combo3 <- abs(outer(freq4, freq3, "-"))
        combo4 <- abs(outer(freq5, freq3, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        freq.min4 <- apply(combo4,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m2.new <- m2[c(freq.min2)]
        m4.new <- m4[c(freq.min3)]
        m5.new <- m5[c(freq.min4)]
        
        f12 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3[[x]]==1 & m4.new[[x]]==1 & m5.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq3), f12)
        combo0 <- abs(outer(freq3, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==4){
        combo1 <- abs(outer(freq1, freq4, "-"))
        combo2 <- abs(outer(freq2, freq4, "-"))
        combo3 <- abs(outer(freq3, freq4, "-"))
        combo4 <- abs(outer(freq5, freq4, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        freq.min4 <- apply(combo4,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m2.new <- m2[c(freq.min2)]
        m3.new <- m3[c(freq.min3)]
        m5.new <- m5[c(freq.min4)]
        
        f13 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1 & m4[[x]]==1 & m5.new[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq4), f13)
        combo0 <- abs(outer(freq4, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
      
      else if(comp==5){
        combo1 <- abs(outer(freq1, freq5, "-"))
        combo2 <- abs(outer(freq2, freq5, "-"))
        combo3 <- abs(outer(freq3, freq5, "-"))
        combo4 <- abs(outer(freq4, freq5, "-"))
        
        freq.min1 <- apply(combo1,2,which.min)
        freq.min2 <- apply(combo2,2,which.min)
        freq.min3 <- apply(combo3,2,which.min)
        freq.min4 <- apply(combo4,2,which.min)
        
        m1.new <- m1[c(freq.min1)]
        m2.new <- m2[c(freq.min2)]
        m3.new <- m3[c(freq.min3)]
        m4.new <- m4[c(freq.min4)]
        
        f14 <- function(x){
          m <- matrix(1, ncol = p, nrow = p)
          m[!(m1.new[[x]]==1 & m2.new[[x]]==1 & m3.new[[x]]==1 & m4.new[[x]]==1 & m5[[x]]==1)] <- 0
          return(m)
        }
        result <- lapply(1:length(freq5), f14)
        combo0 <- abs(outer(freq5, freq0, "-"))
        freq.min0 <- apply(combo0,2,which.min)
        m0 <- result[c(freq.min0)]
        return(m0)
      }
    }
  } 
}
```

Below is the argument in `mnet`. `x` is the connectivity matrix, which can be get from `mconn2`. We will first apply `mconn2` to the dataset `Caltech_0051475_rois_aal.1D` to get a connectivity matrix for each specific  frequncy. Then we will plug the connectivity matrix into `mnet` to get the number of edge for each network.

```{}
mnet(x)
```

```{r,echo=FALSE}
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

```



```{r,message=FALSE}
#read in data
ctrl <- read.table("Caltech_0051475_rois_aal.1D", header = F)
dim(ctrl)
head(ctrl[,1:8])
#standardize the dataset ctrl (optional step)
library(dplyr)
ctrl1 <- ctrl %>% mutate_all(~(scale(.) %>% as.vector))
#select specific frequencies
freq0 <- seq(0.01,0.1,0.01)
freq0
#apply functions to the data
library(astsa)
myconn <- mconn2(x=ctrl1, alpha=0.05, s=1, tt=41, freq0=freq0)
#check type of the output from function mconn2
class(myconn)
#we get a connectivity matrix under each frequency
length(myconn)
#check the dimension of each connecitivty matrix, which should equal 116. This is the 116 altas used in the dataset ctrl.
dim(myconn[[1]])
mynet <- mnet(myconn)
class(mynet)
mynet
```
