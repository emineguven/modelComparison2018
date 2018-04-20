
#EG 20160403 time 12:56 pm

#if the log likelihood function difference of Weibull and Gompertz is always positive
#how does the noise effect the loglikelihood 
#Weibull model->lambda is a scale v is a shape parameter >0
#Gompertz Model-> beta is a scale alpha is a shape parameter>0 
# by simulation I found when x>0, and all parameters are less than 0.3 the difference is negative
#the difference of log likelihood functions has some positive values for parameters>0.3


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz,X) for parameters Weibull model->lambda is a scale v is a shape parameter >0
   #                                                Gompertz Model-> beta is a scale alpha is a shape parameter>0 
# I found between 0 and 0.3 log likelihood difference is negative 
# difference loglikelihood is always negative for 
#restrict beta and lambda to be greater than 0.3 

# initialize log-likelihood difference and 
#For additive Gaussian noise e ~ N (0, sigma^2) with known variance sigma^2
#delta.likelihood=list()
#noise=list()
                   
#sd of gaussian noise function
#max sd would be = 4*mean(inverse.gomp.CDF)
#min sd would be mean(inverse.gomp.CDF)




lambda=0.25
v=0.1
#beta= 0.1#G G= (0.1,0.25)
#alpha=0.001 #R  R= (0.001,0.1)

N=100 # population size
#U<-runif(N)
#generate gompertz random numbers by using inverse CDF
#prediction
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
  
  #(1/(k*beta))*log(1 -(beta/alpha)*log(1-U)) 
  
  x.gompertz = inverse.gomp.CDF(beta,alpha, x.uniform)
  return(x.gompertz)
}
#gompertz.random<-rgompertz(beta,alpha,N)


average.lifespan=mean(gompertz.random)
round.lifespan=ceiling(average.lifespan)

noise<- list()
delta.LL<-list()





 calculate.noise = function(sd){
   
  gaussian<-rnorm(N, mean = 0, sd=sd )
  #observation
  X<- gompertz.random +gaussian
  #noise
  noise=sd(gaussian)
  
  return(noise)
}

 # find delta log-likelihood for each parameter sets
 
# s  <- (1:10) / 10 
 for (beta in seq(from=0.05 , to=0.3, length.out = 2)){
   for (alpha in seq(from=0.0001, to=0.1, length.out=1)){
  #for (i in seq(round.lifespan,3*round.lifespan,by=1)){
     for (sd in seq(round.lifespan,2*round.lifespan,by=2)){
       
       
       #generate gompertz random numbers (lifespan) by using inverse CDF
       #prediction
       
       gompertz.random<-rgompertz(beta,alpha,N)
    
    #standard deviation of random numbers.
    sd.gaussian=calculate.noise(sd)
    gaussian<-rnorm(N, mean = 0, sd=sd.gaussian)
    
    #simulate lifespan plus normal noises
    X<- gompertz.random + gaussian
    
    power.X<-X^v
    power.X<-power.X[complete.cases(power.X)]
    
    delta.likelihood<- N*log(lambda*v)+(v-1)*log(sum(X))+(-lambda)*sum(power.X)- #Weibull -
      (N*log(beta)+alpha*sum(X)+beta/alpha*sum(1-exp(alpha*X)))                 #Gompertz
     
    LL_change<-data.frame(cbind(delta.likelihood)) 
    
    #calculate LL and noise change
    noise[[length(noise)+1]] = sd.gaussian        
    delta.LL[[length(delta.LL)+1]] = delta.likelihood
    
     }
   }
 }
 
 #put the noise into data frame
 df.noise=data.frame(cbind(noise))
 
 #put the LL change into data frame
 df.delta=data.frame(cbind(delta.LL))  
 
 mat<-data.frame(noise=cbind(noise),beta=cbind(seq(from=0.05 , to=0.3, length.out = 2)),delta=cbind(delta.LL))
 
 mat <- mat[order(mat$beta),]
 
 row.names(mat) <- mat$noise
 
 heatmap_matrix <- data.matrix(mat)
 
 new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb") 
 
 
 heatmap(heatmap_matrix,Rowv=NA, Colv=NA, ylab= "difference of log likelihood W-G parameter range is [0.9,9]",
         xlab     = "std of random numbers from weibull dist.",
         main     = "log-likelihood vs std",symm=FALSE,col=new.palette(20))
 
 
# after the loops,  generate matrices for two parameter conditioned on the third parameter for heatmaps


 
 
 





