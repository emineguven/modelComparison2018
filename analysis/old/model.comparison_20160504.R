
#title: "GaussianNoises and Delta_LL of model comparison project"
#author: "e guven"
#date: "May 3, 2016"
#output: html_document


#EG started this on 20160403 
#how does the noise effect the loglikelihood 
#Weibull model->lambda is a scale, v is a shape parameter >0
#Gompertz Model-> beta(G) is a scale, alpha(R) is a shape parameter>0 

#difference function of loglikelihood function of gompertz and weibull p.d.fs
#test if L(Weibull,X)>L(Gompertz,X) for parameters Weibull model->lambda is a scale, v is a shape parameter >0
#Gompertz Model-> beta is a scale, alpha is a shape parameter>0 

#For additive Gaussian noise e ~ N (0, sigma^2) with known variance sigma^2
#sd of gaussian noise function
#max sd would be = 3*mean(inverse.gomp.CDF)
#min sd would be mean(inverse.gomp.CDF)



#fix the paramaters of weibull model

lambda=0.25
v=0.1
#test G and R in nested for loops
beta= 0.05  #G G= (0.1,0.25)
alpha=0.001 #R  R= (0.001,0.1)

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
  gaussian<-rnorm(N, mean = 0, sd=i )
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
      round.lifespan=ceiling(average.lifespan)

#put the noise values from loops into noise and delta.LL in a list
sderr<- list()
Delta_LL<-list() 
G<-list()
R<-list()
LWei<-list()
LGomp<-list()
MeanLF<-list()
sdLS<-list()

for (beta in c(0.05, 0.1, 0.15, 0.2, 0.25)){
  for (alpha in c(1E-4, 0.002, 0.008, 0.01, 0.05)){ #fix alpha or in other words R shape parameter
    #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
    for (i in c(0, 0.5, 1, 2, 5)){ 
      
      
      #generate gompertz random numbers (lifespan) by using inverse CDF
      #prediction
      gompertz.random<-rgompertz(alpha,beta,N)
      average.lifespan=mean(gompertz.random)
      
      
      #results$sderr[count] = i;
      #results$R[count] = alpha
      #results$G[count] = beta;
      MeanLF[[length(MeanLF)+1]]=average.lifespan
     
      
      
      
      sd.gaussian=calculate.noise(i)
      gaussian<-rnorm(N, mean = 0, sd=sd.gaussian)
      #standard deviation of lifespan.
      sd.lifespan=sd(gompertz.random)
      sdLS[[length(sdLS)+1]] =sd.lifespan
      
      
      #simulate lifespan plus normal noises
      #observation
      X<- gompertz.random + gaussian
      #Y<-weibull.random+gaussian
      
      #this is to fix a bug problem here power.X is actually the same thing as X=observation
      power.X<-X^v
      power.X<-power.X[complete.cases(power.X)]
      
      LL_wei <- function(x)
      { wei<-N*log(lambda*v)+(v-1)*log(sum(X))+(-lambda)*sum(power.X) }   #Weibull 
      
      max_LL_wei<-optim(initial_param_Wei<- c(lambda,v), LL_wei, hessian=TRUE) 
      
      LWei[[length(LWei)+1]] =max_LL_wei$value
      
      
      LL_gomp <- function (x)
        {gomp<- (N*log(beta)+alpha*sum(X)+beta/alpha*sum(1-exp(alpha*X))) }  #Gompertz
      
      max_LL_gomp<-optim(initial_param_gomp<- c(beta,alpha), LL_gomp, hessian=TRUE) 
      
      LGomp[[length(LGomp)+1]] =max_LL_gomp$value
      
      delta.likelihood.wei<-max_LL_wei$value-max_LL_gomp$value
      
      #power.Y<-Y^v
      #power.Y<-power.Y[complete.cases(power.Y)]
      
       #delta.likelihood.gomp<- N*log(lambda*v)+(v-1)*log(sum(Y))+(-lambda)*sum(power.Y)- #Weibull -
       # (N*log(beta)+alpha*sum(Y)+beta/alpha*sum(1-exp(alpha*Y)))                 #Gompertz
    
      
      #calculate LL and noise change
      sderr[[length(sderr)+1]] = i       
      Delta_LL[[length(Delta_LL)+1]] = delta.likelihood.wei
      G[[length(G)+1]]=beta
      #switch to alpha.seq when for fixed beta
      R[[length(R)+1]]=alpha
      
      #results$LGomp[count]=deltalikelihood.gomp
      
      #todo use flexsurv to calculate the LL
    
    }
  }
}

results = data.frame(sderr=cbind(sderr), R=cbind(R), G=cbind(G),Delta_LL=cbind(Delta_LL) ,LWei=cbind(LWei), 
                 LGomp=cbind(LGomp), MeanLF=cbind(MeanLF), sdLS=cbind(sdLS))





