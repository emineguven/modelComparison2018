
#EG 20160503 


#how does the noise effect the loglikelihood 
#Weibull model->lambda is a scale, v is a shape parameter >0
#Gompertz Model-> beta is a scale, alpha is a shape parameter>0 


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz,X) for parameters Weibull model->lambda is a scale v is a shape parameter >0
#                                                   Gompertz Model-> beta is a scale alpha is a shape parameter>0 

#For additive Gaussian noise e ~ N (0, sigma^2) with known variance sigma^2

                   
#sd of gaussian noise function
#max sd would be = 4*mean(inverse.gomp.CDF)
#min sd would be mean(inverse.gomp.CDF)



#fix the paramaters of weibull model
lambda=0.25
v=0.1
#initialize
#test G and R in nested for loops
beta= 0.1 #G G= (0.1,0.25)
alpha=0.0001 #R  R= (0.001,0.1)

N=100 # population size

#generate gompertz random numbers by using inverse CDF
#prediction
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
    x.gompertz = inverse.gomp.CDF(alpha,beta, x.uniform)
  return(x.gompertz)
}

#generate gompertz random numbers (lifespan) by using inverse CDF
#prediction
gompertz.random<-rgompertz(alpha,beta,N)
average.lifespan=mean(gompertz.random)
round.lifespan=ceiling(average.lifespan)

#put the noise values from loops into noise and delta.LL in a list
noise<- list()
delta.LL<-list() 
#beta.seq<-list()
alpha.seq<-list()

# create a function that calculates noise of lifespan
 calculate.noise = function(i){
   
   #lifespan
  gaussian<-rnorm(N, mean = 0, sd=i )
  
  #observation
  #X<- gompertz.random +gaussian to be used in simulation later.
  #noise
  noise=sd(gaussian)
  
  return(noise)
}

 # find delta log-likelihood for each parameter sets
 for (beta in seq(from=0.05 , to=0.3, length.out = 1)){
   for (alpha in seq(from=0.0001, to=0.1, length.out=1500)){ #fix alpha or in other words R shape parameter
  #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
     for (i in seq(from=0, to=5, by=1)){ 
       
       
       #generate gompertz random numbers (lifespan) by using inverse CDF
       #prediction
       gompertz.random<-rgompertz(alpha,beta,N)
       average.lifespan=mean(gompertz.random)
       round.lifespan=ceiling(average.lifespan)
    
    #standard deviation of lifespan.
    sd.gaussian=calculate.noise(i)
    gaussian<-rnorm(N, mean = 0, sd=sd.gaussian)
    
    #simulate lifespan plus normal noises
    #observation
    X<- gompertz.random + gaussian
    
    #this is to fix a bug problem here power.X is actually the same thing as X=observation
    power.X<-X^v
    power.X<-power.X[complete.cases(power.X)]
    
    delta.likelihood<- N*log(lambda*v)+(v-1)*log(sum(X))+(-lambda)*sum(power.X)- #Weibull -
      (N*log(beta)+alpha*sum(X)+beta/alpha*sum(1-exp(alpha*X)))                 #Gompertz
     
    #calculate LL and noise change
    noise[[length(noise)+1]] = sd.gaussian        
    delta.LL[[length(delta.LL)+1]] = delta.likelihood
    #beta.seq[[length(beta.seq)+1]]=beta
    #switch to alpha.seq when for fixed beta
    alpha.seq[[length(alpha.seq)+1]]=alpha
    
     }
   }
 }
 
 #put the noise into data frame
 df.noise=data.frame(cbind(noise))
 
 #put the LL change into data frame
 df.delta=data.frame(cbind(delta.LL))  
 
 #put everything in one data frame
 mat.data3<-data.frame(beta=cbind(alpha.seq), noise=cbind(noise),delta=cbind(delta.LL))

 
 
 #data<-read.csv("only_noise_data.csv")
 
 #row.names(data) <- data$X
 #heatmap_matrix <- data.matrix(data)
 #heatmap(heatmap_matrix)
 
# new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb") 
 
 
#heatmap(heatmap_matrix,Rowv=NA, Colv=NA, ylab= "difference of log likelihood W-G parameter range is [0.9,9]",
        # xlab     = "std of random numbers from weibull dist.",
        # main     = "log-likelihood vs std",symm=FALSE)
 
 
# after the loops,  generate matrices for two parameter conditioned on the third parameter for heatmaps


 
 
 





