# install.packages("iZID")
library(iZID)

#step 1
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

check.dist(x1,x2)

#step 2

#run functions to generate sample means from ZIP and ZNB
# generates one i.i.d. sample of size n from ZiIP
# input: estimated lambda, estimated mixing probability, sample size
# output: a sample of size n
boot.zip = function(lambda.est,prob.est,n=100,sn){
  
  zip.sample = rep(0,n)
  
  set.seed(sn)
  unif = runif(n,0,1)
  
  for(i in 1:n){
    
    if(unif[i]<prob.est) zip.sample[i] = 0
    if(unif[i]>=prob.est) zip.sample[i] = rpois(1,lambda.est) 
    
  } # end loop over i
  
  return(zip.sample)  
}

# obtain B bootstrap samples of the mean of ZIP
boot.zip.mean = function(lambda.est,prob.est,n=100,B=500,sn){
  
  boot.zip.mean.sample = rep(0,B)
  for(i in 1:B){
    set.seed(i*sn)
    boot.zip.mean.sample[i] = mean(boot.zip(lambda.est,prob.est,n,i*sn))
  }
  
  outl = list(boot.zip.mean.sample , 
              quantile(boot.zip.mean.sample,probs = c(0.025)),
              quantile(boot.zip.mean.sample,probs = c(0.975 )) )
  return(outl)
} 

# generates one i.i.d. sample of size n from ZiNB
# input: estimated r, estimated p, estimated mixing probability, sample size
# output: a sample of size n
boot.znb = function(r.est,p.est,prob.est,n=100,sn){
  
  znb.sample = rep(0,n)
  
  set.seed(sn)
  unif = runif(n,0,1)
  
  for(i in 1:n){
    
    if(unif[i]<prob.est) znb.sample[i] = 0
    if(unif[i]>=prob.est) znb.sample[i] = rnbinom(1,size=r.est,prob=p.est) 
    
  } # end loop over i
  
  return(znb.sample)  
} # end function boot.znb()
# obtain B bootstrap samples of the mean of ZNB
boot.znb.mean = function(r.est,p.est,prob.est,n=100,B=500,sn){
  
  boot.znb.mean.sample = rep(0,B)
  for(i in 1:B){
    set.seed(i*sn)
    boot.znb.mean.sample[i] = mean(boot.znb(r.est,p.est,prob.est,n,i*sn))
  }
  
  outl = list(boot.znb.mean.sample , 
              quantile(boot.znb.mean.sample,probs = c(0.025)),
              quantile(boot.znb.mean.sample,probs = c(0.975 )) )
  return(outl)
} # end function boot.znb.mean

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
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      n <- d1$N
      lambda1 <- d1$mle_ori[1]
      prob1 <- d1$mle_ori[2]
      boot.avg1 <- boot.zip.mean(lambda.est=lambda1,prob.est=prob1,n,B,sn=1)
      lambda2 <- d2$mle_ori[1]
      prob2 <- d2$mle_ori[2]
      boot.avg2 <- boot.zip.mean(lambda.est=lambda2,prob.est=prob2,n,B,sn=1)
      
      pval = 0
      if( (mean(x1)>=boot.avg2[[2]] & mean(x1)<=boot.avg2[[3]] ) || 
          (mean(x2)>=boot.avg1[[2]] & mean(x2)<=boot.avg1[[3]]) ){
        pval = 1
      }
      
      # boot.stat <- abs(boot.avg1[[1]]-boot.avg2[[1]])
      # pval <- sum(boot.stat >= abs(mean(x1)-mean(x2)))/B
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
      
    }else if(dist1==2&dist2==2){#2: zinb
      d1 <- dis.kstest(x1,nsim=100,bootstrap=TRUE,distri=aa[1,2])
      d2 <- dis.kstest(x2,nsim=100,bootstrap=TRUE,distri=aa[2,2])
      n <- d1$N
      r1 <- d1$mle_ori[1]
      p1 <- d1$mle_ori[2]
      prob1 <- d1$mle_ori[3]
      boot.avg1 <- boot.znb.mean(r.est=r1,p.est=p1,prob.est=prob1,n,B,sn=1)
      r2 <- d2$mle_ori[1]
      p2 <- d2$mle_ori[2]
      prob2 <- d2$mle_ori[3]
      boot.avg2 <- boot.znb.mean(r.est=r1,p.est=p2,prob.est=prob2,n,B,sn=1)
      pval = 0
      if( (mean(x1)>=boot.avg2[[2]] & mean(x1)<=boot.avg2[[3]] ) || 
          (mean(x2)>=boot.avg1[[2]] & mean(x2)<=boot.avg1[[3]]) ){
        pval = 1
      }
      # boot.stat <-abs(boot.avg1[[1]]-boot.avg2[[1]])
      # pval <- mean(boot.stat >= abs(mean(x1)-mean(x2)))
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
      # boot.stat <- abs(boot.avg1-boot.avg2)
      # pval <- mean(boot.stat >= abs(mean(x1)-mean(x2)))
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
      # boot.stat <- abs(boot.avg1-boot.avg2)
      # pval <- mean(boot.stat >= abs(mean(x1)-mean(x2)))
      hyp <- ifelse(pval>0.05, "H0", "H1")
      avg1 <- round(mean(x1),4)
      avg2 <- round(mean(x2),4)
      result <- data.frame(hyp, avg1, avg2)
      colnames(result) <- c("Hypothesis", "Mean 1", "Mean 2")
      return(result)
    }
  }
}

#apply the function to real data
set.seed(2024)
#data (network: SRN, age>=18, freq=0.01)
x1 = c(1,5,4,0,8,0,2,1,5,2,5,1,0,4,4,3,0,2,3,1,4,3,2,5,3,3,3,2,0,1)
x2 = c(1,3,3,1,2,1,7,4,0,3,7,3,15,6,1,2,2,1,4,3,0,3,4,2,16,3,3,6,0,2)
#apply the function
comp.dist(x1,x2)
