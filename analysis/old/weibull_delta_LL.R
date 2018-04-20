require(flexsurv)

#lambda=0.025
#v=0.001
#test G and R in nested for loops
beta= 0.034  #G G= (0.1,0.25)
alpha=0.01 #R  R= (0.001,0.1)

N=100 # population size

#generate gompertz random numbers by using inverse CDF
#prediction
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
  x.gompertz = inverse.gomp.CDF(alpha,beta, x.uniform)
  return(x.gompertz)
}


rweibull= function(lambda,v,n)
{
  x.uniform= runif(n)
  inverse.wei.CDF=function(lambda,v,y) { lambda*(-log(1-y))^(1/v)}
  x.weibull=inverse.wei.CDF(lambda,v,x.uniform)
  return(x.weibull)
}

#create a function that calculates noise of lifespan
calculate.noise = function(i){
  
  #lifespan
  gaussian<-rnorm(N, mean = 0, sd=i)
  #observation
  #X<- gompertz.random +gaussian to be used in simulation later.
  #noise
  noise=sd(gaussian)
  
  return(noise)
}




#generate gompertz random numbers (lifespan) by using inverse CDF
#prediction
gompertz.random<-rgompertz(alpha,beta,N)
average.lifespan=mean(gompertz.random)

#put the noise values from loops into noise and delta.LL in a list
sderr<- list()
Delta_LL<-list() 
G<-list()
R<-list()
LWei<-list()
LGomp<-list()
MeanLF<-list()
sdLS<-list()
Delta_LL.flex<-list()
LWei.flex<-list()
LGomp.flex<-list()
G.flex.estimated<-list()
R.flex.estimated<-list()
LLG.par<-list()
LLR.par<-list()
v.flex.estimated<-list()
lambda.flex.estimated<-list()

for (beta in c(0.05,0.08, 0.1,0.15,0.17, 0.2, 0.25)){
  for (alpha in c(1E-3, 0.002, 0.005,0.008, 0.01,0.03, 0.05)){ #fix alpha or in other words R shape parameter
    #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
    for (i in c(0, 0.5, 1,2,3,4, 5)){ 
      #for (i in c(0)){ 
      
      #generate gompertz random numbers (lifespan) by using inverse CDF
      #prediction
      gompertz.random<-rgompertz(alpha,beta,N)
      average.lifespan=mean(gompertz.random)
      
      
      #results$sderr[count] = i;
      #results$R[count] = alpha
      #results$G[count] = beta;
      MeanLF[[length(MeanLF)+1]]=average.lifespan
      
      
      
      
      sd.gaussian=calculate.noise(i)
      gaussian<-rnorm(N, mean = 2*average.lifespan, sd=sd.gaussian)
      #standard deviation of lifespan.
      sd.lifespan=sd(gompertz.random)
      sdLS[[length(sdLS)+1]] =sd.lifespan
      
      
      #simulate lifespan plus normal noises
      #observation
      #gompertz.random<-rgompertz(alpha,beta,100)
      #x<- gompertz.random + gaussian
      
      
      lifespan<- gompertz.random +gaussian
      
      #Log likelihood function for the Weibull model
      
      
      weib.likl<-function(param,y){
        theta<-exp(param[1])
        gamma<-exp(param[2])
        delta=1;
        
        y=lifespan[!is.na(lifespan)]
        logl<-sum(delta*(log(gamma) + gamma*log(theta) + (gamma-1)*log(y) -
                           (theta*y)^gamma )) -sum((1-delta)*(theta*y)^gamma)
       
        return(-logl)
          }
      weib=optim(log(c(0.03,0.01)),weib.likl,y=lifespan)
      
      weib$value
      LWei[[length(LWei)+1]] = weib$value
    }
  }
}