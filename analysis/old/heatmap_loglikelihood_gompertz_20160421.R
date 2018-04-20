
#EG 20160403 time 12:56 pm

#if the log likelihood function difference of Weibull and Gompertz is always positive
# this is just to simulate if it is really a fact we usually get the positive values or not.
#Weibull model->lambda is a scale v is a shape parameter >0
#Gompertz Model-> beta is a scale alpha is a shape parameter>0 
# by simulation I found when x>0, and all parameters are less than 0.3 the difference is negative
#the difference of log likelihood functions has some positive values for parameters>0.3

#when k=0.1 which means parameter range is [0.1,1] log-likelihood difference is always negative
#when k=0.2 which means parameter range is [0.2,2] log-likelihood difference is always negative
#when k=0.3                                [0.3,3]                                     negative
#when k=0.4 the result start to have more on positive values see figure 

i=2;k=0.4;N=100


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz,X) for parameters Weibull model->lambda is a scale v is a shape parameter >0
   #                                                Gompertz Model-> beta is a scale alpha is a shape parameter>0 
# I found between 0 and 0.3 log likelihood difference is negative 
# difference loglikelihood is always negative for 
#restrict beta and lambda to be greater than 0.3 

# initialize log-likelihood difference and 
#For additive Gaussian noise e ~ N (0, sigma^2) with known variance sigma^2
delta.likelihood=array(NA, dim = c(i,i,i,i))
noise=array(NA, dim = c(i,i)) 
                   
#sd of gaussian noise function
 calculate.noise = function(beta,alpha){
  
  U<-runif(N)
  #prediction
  inverse.gomp.CDF <- (1/k*beta)*log(1 -(beta/alpha)*log(1-U)) 
  gaussian<-rnorm(N, mean = 0, sd = 1)
  #observation
  X<- inverse.gomp.CDF +gaussian
  #noise
  noise=sd(gaussian)
  
  return(noise)
}

 # find delta log-likelihood for each parameter sets
 find.bigger.log = function(lambda,v,beta,alpha){
  
  U<-runif(N)
  #generate gompertz random numbers by using inverse CDF
  inverse.gomp.CDF <- (1/k*beta)*log(1 -(beta/alpha)*log(1-U)) 
  X<- inverse.gomp.CDF +rnorm(N)
  delta.likelihood<-N^2*(k*v-1)*log(k*(lambda+v))*(log(sum(X))-lambda*sum(X)^k*v)- #Weibull -
    (N*k*beta*log(k*beta)*sum(X)*(1-exp(k*alpha*sum(X))))                          #Gompertz
  
  return(delta.likelihood)
 
}






for (lambda in c(1:i)){
  
  for (v in c(1:i)){
    
    for (beta in c(1:i)){
    
    for (alpha in c(1:i)){
    
    #the log likelihood difference is stored in output matrix
    delta.likelihood[lambda,v,beta,alpha]<-find.bigger.log(lambda,v,beta,alpha)
    
  }
  
    }}}

for (beta in c(1:i)){
  
  for (alpha in c(1:i)){
    
    #standard deviation of random numbers.
    noise[beta,alpha]<-calculate.noise(beta,alpha)
    
  }
  
}

  #noise and delta.likelihood row data frame
  DF <- data.frame(noise = as.vector(list(noise)), difference = as.vector(list(delta.likelihood)))
  
  #convert data frame to a matrix for heatmap use
  heatmap_matrix<-matrix(delta.likelihood,nrow=i^2,ncol=i^2)
  
  # omit -inf values to do the heat map
  #first replace -inf and/or inf values with NA
  is.na(heatmap_matrix) <- sapply(heatmap_matrix, is.infinite)
  
  #omit columns with NA values which were -inf and inf initially 
  #because heatmap.2() does hiearchical clustering where no NA,NAN and -inf or inf values are accepted
  heatmap_matrix<-heatmap_matrix[, colSums(is.na(heatmap_matrix)) < nrow(heatmap_matrix) * 0.001]
  
  #take the log of heat map because in row data frame DF there were too big powers of e
  log.heat.matrix<-log(heatmap_matrix)
  
  #just to double check there is no NA/NAN or inf/-inf values
  log.heat.matrix<-log.heat.matrix[, colSums(is.na(log.heat.matrix)) < nrow(log.heat.matrix) * 0.001]
  df<-as.data.frame(log.heat.matrix)

  write.csv(df, "test_differenceGompertz.csv")


  library(gplots)
  library(ColorPalette)
  
  #again convert log likelihood to matrix and scale it
  data_mat<-scale(as.matrix(df))
  # following code limits the lowest and highest color to 5%, and 95% of your range, respectively
  quantile.range <- quantile(data_mat, probs = seq(0, 1, 0.01))
  palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  
  # use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
  
  jpeg(file = "~/github/model.comparison/plotsGompertz/Heatplot_[0.4,0.8].jpeg")
  
  
  heatmap.2(data_mat,dendrogram = "none",
            scale      = "none",
            trace      = "none",
            key        = TRUE,
            labRow     = NA,  
            ylab     = "delta_log-likelihood of W-G, parameter range is [0.4,0.8]",
            xlab     = "noise",
            main     = "noise vs delta_log-likelihood",
            col    = color.palette,
            breaks = palette.breaks) 
  
 dev.off()
  







