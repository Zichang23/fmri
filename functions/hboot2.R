library(iZID)

#function to check dist for x1 or x2
check.dist <- function(x1, x2, alpha=0.05, B=100){#x1,x2 are vectors
  # d1=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zip')
  # d2=dis.kstest(x1,nsim=100,bootstrap=TRUE,distri='zinb')
  d3=dis.kstest(x1,nsim=B,bootstrap=TRUE,distri='Poisson')
  d4=dis.kstest(x1,nsim=B,bootstrap=TRUE,distri='nb')
  pval1 <- c(d3$pvalue, d4$pvalue)
  llik1 <- round(c(d3$mle_ori[2], d4$mle_ori[3]),5)
  dist1 <- c("Poisson", "nb")
  id <- 1:2
  dat1 <- data.frame(cbind(id, dist1, pval1, llik1))
  dat2 <- dat1[dat1$pval1 >=alpha,]
  dat3 <- dat2[dat2$llik1==max(as.numeric(dat2$llik1)),]
  dat4 <- dat3[dat3$pval1==max(as.numeric(dat3$pval1)),]
  avg1 <- round(mean(x1),4)
  var1 <- round(var(x1),4)
  result1 <- cbind(dat4[,1:2], avg1, var1)
  d7=dis.kstest(x2,nsim=B,bootstrap=TRUE,distri='Poisson')
  d8=dis.kstest(x2,nsim=B,bootstrap=TRUE,distri='nb')
  pval2 <- c(d7$pvalue, d8$pvalue)
  llik2 <- round(c(d7$mle_ori[2], d8$mle_ori[3]),5)
  dist2 <- c("Poisson", "nb")
  dat5 <- data.frame(cbind(id, dist2, pval2, llik2))
  dat6 <- dat5[dat5$pval2 >=alpha,]
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

#function to compare ASD gorup (x1) and control group (x2)
comp.dist1 <- function(x1, x2, alpha=0.05, B=500){#x1,x2 are vectors
  aa <- check.dist1(x1,x2)
  dist1 <- as.numeric(aa[1,1])
  dist2 <- as.numeric(aa[2,1])
  if(dist1 !=dist2){
    avg1 <- round(mean(x1),4)
    avg2 <- round(mean(x2),4)
    result <- data.frame(mean1=avg1, mean2=avg2)
    return(result)
    
  }else if(dist1 == dist2){
    if(dist1==1&dist2==1){#1: Poisson
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
      # boot.stat <- abs(boot.avg1-boot.avg2)
      # pval <- mean(boot.stat >= abs(mean(x1)-mean(x2)))
      hyp <- ifelse(pval>alpha, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
      
    }else if(dist1==2&dist2==2){#2: Neg Bin
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
      # boot.stat <- abs(boot.avg1-boot.avg2)
      # pval <- mean(boot.stat >= abs(mean(x1)-mean(x2)))
      hyp <- ifelse(pval>alpha, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
    }
  }
}
