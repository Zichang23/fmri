#create function mconn to calculate connectivity matrix
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

library(astsa)
#read in data
x <- read.table("Leuven_2_0050732_rois_aal.1D", header = F)
#selected frequencies
freq0 <- seq(0.01,0.1,0.016)
#apply the function to the data
mconn2(x=x, alpha=0.01, s=1, tt=82, freq0=freq0)
