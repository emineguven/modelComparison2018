
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

i=10;k=1.5;N=100


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz) for parameters beta,lambda>0.3
# I found between 0 and 0.3 log likelihood difference is negative 
# difference loglikelihood is always negative for 
#restrict beta and lambda to be greater than 0.3 


difference.log=array(NA, dim = c(i,i,i,i))
                     #dimnames=list( d1=c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10"),
                                    #,"L11","L12","L13","L14","L15","L16","L17","L18","L19","L20",
                                    #"L21","L22","L23","L24","L25","L26","L27","L28","L29","L30"), 
                                    #d2=c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),
                                    #"v11","v12","v13","v14","v15","v16","v17","v18","v19","v20",
                                    # "v21","v22","v23","v24","v25","v26","v27","v28","v29","v30"), 
                                    #d3=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10"),
                                    #d4=c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10")))
#"a11","a12","a13","a14","a15","a16","a17","a18","a19","a20",
#  "a21","a22","a23","a24","a25","a26","a27","a28","a29","a30")))

std.gompertzR=array(NA, dim = c(i,i)) 
                    #dimnames=list( s1=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10"),
                                  # s2=c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10")))
#"a11","a12","a13","a14","a15","a16","a17","a18","a19","a20",
#  "a21","a22","a23","a24","a25","a26","a27","a28","a29","a30")))

random.generator = function(beta,alpha){
  
  U<-runif(N)
  inverse.gomp.CDF <- (1/k*beta)*log(1 -(beta/alpha)*log(1-U))
  
  std=sd(inverse.gomp.CDF)
  
  return(std)
}

find.bigger.log = function(lambda,v,beta,alpha){
  
  #unknown parameter coeff ,k tests 

  U<-runif(N)
  inverse.gomp.CDF <- (1/k*beta)*log(1 -(beta/alpha)*log(1-U))
  
  difference<-N^2*(k*v-1)*log(k*(lambda+v))*(log(sum(inverse.gomp.CDF))-lambda*sum(inverse.gomp.CDF)^k*v)-
    (N*k*beta*log(k*beta)*sum(inverse.gomp.CDF)*(1-exp(k*alpha*sum(inverse.gomp.CDF))))
  
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

for (beta in c(1:i)){
  
  for (alpha in c(1:i)){
    
    #standard deviation of random numbers from Gompertz dist.
    std.gompertzR[beta,alpha]<-random.generator(beta,alpha)
    
  }
  
}


  DF <- data.frame(noise = as.vector(list(std.gompertzR)), difference = as.vector(list(difference.log)))
  
  heatmap_matrix<-matrix(difference.log,nrow=i^2,ncol=i^2)
  
  is.na(heatmap_matrix) <- sapply(heatmap_matrix, is.infinite)
  heatmap_matrix<-heatmap_matrix[, colSums(is.na(heatmap_matrix)) < nrow(heatmap_matrix) * 0.001]
  log.heat.matrix<-log(heatmap_matrix)
  log.heat.matrix<-log.heat.matrix[, colSums(is.na(log.heat.matrix)) < nrow(log.heat.matrix) * 0.001]
  df<-as.data.frame(log.heat.matrix)

  write.csv(df, "test_differenceGompertz.csv")


  library(gplots)
  library(ColorPalette)
  
  data_mat<-scale(as.matrix(df))
  # following code limits the lowest and highest color to 5%, and 95% of your range, respectively
  quantile.range <- quantile(data_mat, probs = seq(0, 1, 0.01))
  palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  
  # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
  
  jpeg(file = "~/Documents/github/model.comparison/plots/Heatplot_[1.5,15].jpeg")
  
  
  heatmap.2(data_mat,dendrogram = "none",
            scale      = "none",
            trace      = "none",
            key        = TRUE,
            labRow     = NA,  
            ylab     = "difference of log likelihood W-G parameter range is [1.5,15]",
            xlab     = "std of random numbers from gompertz dist.",
            main     = "log-likelihood vs std",
            col    = color.palette,
            breaks = palette.breaks) 
  
  dev.off()
  
  
#DF<-as.matrix(DF)
#write.csv(DF, "test_difference.csv")
#data<-read.csv("test_difference_infsRemoved.csv")


#heatmap(heatmap_matrix)







