library(iZID)

#step 1: check distribution for data
check.dist <- function(x1,x2){#x1,x2 are vectors
  d1=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zip')
  d2=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zinb')
  d3=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='Poisson')
  d4=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='nb')
  pval1 <- c(d1$pvalue, d2$pvalue, d3$pvalue, d4$pvalue)
  llik1 <- round(c(d1$mle_ori[3], d2$mle_ori[4], d3$mle_ori[2], d4$mle_ori[3]),5)
  dist1 <- c("zip", "zinb", "Poisson", "nb")
  id <- 1:4
  dat1 <- data.frame(cbind(id, dist1, pval1, llik1))
  dat2 <- dat1[dat1$pval1 >=0.05,]
  dat3 <- dat2[dat2$llik1==max(as.numeric(dat2$llik1)),]
  dat4 <- dat3[dat3$pval1==max(as.numeric(dat3$pval1)),]
  avg1 <- round(mean(x1),4)
  var1 <- round(var(x1),4)
  result1 <- cbind(dat4[,1:2], avg1, var1)
  d5=dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='zip')
  d6=dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='zinb')
  d7=dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='Poisson')
  d8=dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='nb')
  pval2 <- c(d5$pvalue, d6$pvalue, d7$pvalue, d8$pvalue)
  llik2 <- round(c(d5$mle_ori[3], d6$mle_ori[4], d7$mle_ori[2], d8$mle_ori[3]),5)
  dist2 <- c("zip", "zinb", "Poisson", "nb")
  dat5 <- data.frame(cbind(id, dist2, pval2, llik2))
  dat6 <- dat5[dat5$pval2 >=0.05,]
  dat7 <- dat6[dat6$llik2==max(as.numeric(dat6$llik2)),]
  dat8 <- dat7[dat7$pval2 == max(as.numeric(dat7$pval2)),]
  avg2 <- round(mean(x2),4)
  var2 <- round(var(x2),4)
  result2 <- cbind(dat8[,1:2], avg2, var2)
  names(result1) <- c("ID","Distribution", "Mean", "Variance")
  names(result2) <- c("ID","Distribution", "Mean", "Variance")
  result <-rbind(result1, result2)
  rownames(result) <-c("data1", "data2")
  return(result)
}
aa <- check.dist(x1,x2)


#step 2: compare distribution for ASD and control groups
comp.dist <- function(x1, x2, B=500){#x1,x2 are vectors
  aa <- check.dist(x1,x2)
  dist1 <- as.numeric(aa[1,1])
  dist2 <- as.numeric(aa[2,1])
  if(dist1 !=dist2){
    avg1 <- round(mean(x1),4)
    avg2 <- round(mean(x2),4)
    result <- data.frame(mean1=avg1, mean2=avg2)
    return(result)
    
  }else if(dist1 == dist2){
    if (dist1==1&dist2==1){#1: zip
      # d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      # d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zip')
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='zip')
      N <- d1$N+d2$N
      n <- d1$N
      r <- n/N
      lambda1 <- d1$mle_ori[1]
      omega1 <- d1$mle_ori[2]
      lambda2 <- d2$mle_ori[1]
      omega2 <- d2$mle_ori[2]
      p1 <- omega1+(1-omega1)*exp(-lambda1)
      p2 <- omega2+(1-omega2)*exp(-lambda2)
      pbar <- r*p1+(1-r)*p2
      X1_sq <- (p1-p2)^2/(pbar*(1-pbar)*(1/r+1/(1-r))/N)
      n1_nz <- N*r*(1-p1) #nz: nonzero
      n2_nz <- N*(1-r)*(1-p2)
      N_nz <- n1_nz+n2_nz
      x1_nz <- x1[x1!=0]
      x2_nz <- x2[x2!=0]
      y1_gt_y2 <- length(x1[x1>x2])#gt:greater than
      y1_eq_y2 <- length(x1[x1==x2])#eq: equal
      tie <- unname(table(x1[x1==x2])[table(x1[x1==x2]) >1])
      #https://stat.ethz.ch/pipermail/r-help/2011-November/296349.html
      U <- (y1_gt_y2+0.5*y1_eq_y2-0.5*n1_nz*n2_nz)/sqrt((n1_nz*n2_nz*(n1_nz+n2_nz+1)/12)*((N_nz)^3-N_nz-sum(tie^3-tie)))
      X2_sq <- X1_sq+U^2
      pval <- pchisq(q = X2_sq, df = 2, lower.tail = FALSE)
      
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
      
    }else if(dist1==2&dist2==2){#2: zinb
      # d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      # d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zinb')
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri='zinb')
      N <- d1$N+d2$N
      n <- d1$N
      r <- n/N
      r.1 <- d1$mle_ori[1]
      p.1 <- d1$mle_ori[2]
      r.2 <- d2$mle_ori[1]
      p.2 <- d2$mle_ori[2]
      omega1 <- d1$mle_ori[3]
      omega2 <- d2$mle_ori[3]
      p1 <- omega1+(1-omega1)*choose(r.1-1,0)*p.1^r.1
      p2 <- omega2+(1-omega2)*choose(r.2-1,0)*p.2^r.2
      pbar <- r*p1+(1-r)*p2
      X1_sq <- (p1-p2)^2/(pbar*(1-pbar)*(1/r+1/(1-r))/N)
      n1_nz <- N*r*(1-p1) #nz: nonzero
      n2_nz <- N*(1-r)*(1-p2)
      N_nz <- n1_nz+n2_nz
      x1_nz <- x1[x1!=0]
      x2_nz <- x2[x2!=0]
      y1_gt_y2 <- length(x1[x1>x2])#gt:greater than
      y1_eq_y2 <- length(x1[x1==x2])#eq: equal
      tie <- unname(table(x1[x1==x2])[table(x1[x1==x2]) >1])
      U <- (y1_gt_y2+0.5*y1_eq_y2-0.5*n1_nz*n2_nz)/sqrt((n1_nz*n2_nz*(n1_nz+n2_nz+1)/12)*((N_nz)^3-N_nz-sum(tie^3-tie)))
      X2_sq <- X1_sq+U^2
      pval <- pchisq(q = X2_sq, df =  2, lower.tail = FALSE)
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
      
    }else if(dist1==3&dist2==3){#3: Poisson
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      n <- d1$N
      lambda1 <- d1$mle_ori[1]
      set.seed(1)
      dat1 <- matrix(rpois(n*B, lambda1), nrow=n, ncol=B)
      boot.avg1 <- colMeans(dat1)
      lambda2 <- d2$mle_ori[1]
      dat2 <- matrix(rpois(n*B, lambda2), nrow=n, ncol=B)
      boot.avg2 <- colMeans(dat2)
      pval = 0
      if( (mean(x1)>=quantile(boot.avg2,probs = c(0.025)) & mean(x1)<=quantile(boot.avg2,probs = c(0.975)) ) || 
          (mean(x2)>=quantile(boot.avg1,probs = c(0.025)) & mean(x2)<=quantile(boot.avg1,probs = c(0.975))) ){
        pval = 1
      }
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
      
    }else if(dist1==4&dist2==4){#4: Neg Bin
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      n <- d1$N
      r1 <- d1$mle_ori[1]
      p1 <- d1$mle_ori[2]
      set.seed(1)
      dat1 <- matrix(rnbinom(n*B, r1, p1), nrow=n, ncol=B)
      boot.avg1 <- colMeans(dat1)
      r2 <- d2$mle_ori[1]
      p2 <- d2$mle_ori[2]
      dat2 <- matrix(rnbinom(n*B, r2, p2), nrow=n, ncol=B)
      boot.avg2 <- colMeans(dat2)
      pval = 0
      if( (mean(x1)>=quantile(boot.avg2,probs = c(0.025)) & mean(x1)<=quantile(boot.avg2,probs = c(0.975)) ) || 
          (mean(x2)>=quantile(boot.avg1,probs = c(0.025)) & mean(x2)<=quantile(boot.avg1,probs = c(0.975))) ){
        pval = 1
      }
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
    }
  }
}
