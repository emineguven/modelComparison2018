
#EG 20160403 time 12:56 pm
#EG is still trying to prove 
#if the log likelihood function difference of Weibull and Gompertz is always positive
# this is just to simulate if it is really a fact we usually get the positive values or not.
# by simulation I found when x>0, beta and lambda>0.3 
#the difference of log likelihood functions is always positive

#generate postive random variables for gompertz and weibull
#x has to be positive

#x = rnorm(n=N, mean=20, sd=3)
#x <- rnorm(N, mean = 0, sd = 1)

#log normal random numbers
#popSize = 1000
#x = rlnorm(n=N, mean=3, sd=1)
#noise<-sd(v)
#when k=0.1 which means parameter range is [0.1,1] log-likelihood difference is always negative
#when k=0.2 which means parameter range is [0.2,2] log-likelihood difference is always negative
#when k=0.3                                [0.3,3]                                     negative
#when k=0.4 the result start to have more on positive values see figure 

i=10;k=0.9;N=100


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz) for parameters beta,lambda>0.3
# I found between 0 and 0.3 log likelihood difference is negative 
# difference loglikelihood is always negative for 
#restrict beta and lambda to be greater than 0.3 


difference.log=array(NA, dim = c(i,i,i,i))
                     

#std.gompertzR=array(NA, dim = c(i,i)) 

std.weibullR = array (NA,dim=c(i,i))
                   

#random.generator = function(beta,alpha){
  
 # U<-runif(N)
  #inverse.gomp.CDF <- (1/k*beta)*log(1 -(beta/alpha)*log(1-U))
  
  #std=sd(inverse.gomp.CDF)
  
  #return(std)
#}




weibull.random.generator <- function(lambda,v)
{
  U <- runif(N)
  inverse.weibull.CDF <- lambda*k*(-log(1-U))^(1/v*k)
  
  std=sd(inverse.weibull.CDF)
  
  return(std)
}




find.bigger.log = function(lambda,v,beta,alpha){
  
  #unknown parameter coeff ,k tests 

  U<-runif(N)
  inverse.weibull.CDF <- lambda*k*(-log(1-U))^(1/v*k)
  
  difference<-N^2*(k*v-1)*log(k*(lambda+v))*(log(sum(inverse.weibull.CDF))-lambda*sum(inverse.weibull.CDF)^k*v)-
    (N*k*beta*log(k*beta)*sum(inverse.weibull.CDF)*(1-exp(k*alpha*sum(inverse.weibull.CDF))))
  
  return(difference)
 
}






for (lambda in c(1:i)){
  
  for (v in c(1:i)){
    
    for (beta in c(1:i)){
    
    for (alpha in c(1:i)){
    
    #the log likelihood difference is stored in output matrix
    difference.log[lambda,v,beta,alpha]<-find.bigger.log(lambda,v,beta,alpha)
    
  }
  
    }}}

for (lambda in c(1:i)){
  
  for (v in c(1:i)){
    
    #standard deviation of random numbers from Gompertz dist.
    std.weibullR[lambda,v]<-weibull.random.generator(lambda,v)
    
  }
  
}


  DF <- data.frame(noise = as.vector(list(std.weibullR)), difference = as.vector(list(difference.log)))
  
  heatmap_matrix<-matrix(difference.log,nrow=i^2,ncol=i^2)
  
  is.na(heatmap_matrix) <- sapply(heatmap_matrix, is.infinite)
  heatmap_matrix<-heatmap_matrix[, colSums(is.na(heatmap_matrix)) < nrow(heatmap_matrix) * 0.001]
  log.heat.matrix<-log(heatmap_matrix)
  log.heat.matrix<-log.heat.matrix[, colSums(is.na(log.heat.matrix)) < nrow(log.heat.matrix) * 0.001]
  df<-as.data.frame(log.heat.matrix)

  write.csv(df, "test_differenceWeibull.csv")
  #data<-read.csv("test_difference_infsRemoved.csv")



 
  library(gplots)
  library(ColorPalette)
  
  data_mat<-scale(as.matrix(df))



  
  
  jpeg(file = "~/github/model.comparison/plotsWeibull/Heatplot_Weibull_[0.9,9].jpeg")
new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb") 


heatmap(data_mat,Rowv=NA, Colv=NA, ylab= "difference of log likelihood W-G parameter range is [0.9,9]",
        xlab     = "std of random numbers from weibull dist.",
        main     = "log-likelihood vs std",symm=FALSE,col=new.palette(20))
dev.off()


 
#heatmap.2(data_mat,symm=FALSE,dendrogram = "none",
          #scale      = "none",
          #trace      = "none",
          #key        = TRUE,
          #labRow     = NA,  
          #ylab     = "difference of log likelihood W-G parameter range is [0.4,5]",
          #xlab     = "std of random numbers from weibull dist.",
          #main     = "log-likelihood vs std",col=new.palette(20))
  
  

 dev.off()
  
  
#DF<-as.matrix(DF)
#write.csv(DF, "test_difference.csv")
#data<-read.csv("test_difference_infsRemoved.csv")


#heatmap(heatmap_matrix)







