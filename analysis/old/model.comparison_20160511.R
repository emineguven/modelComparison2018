
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

for (beta in c(0.05,0.08, 0.1,0.12,0.15,0.17, 0.2, 0.25)){
  for (alpha in c(1E-4,1E-3, 0.002, 0.005,0.008, 0.01,0.03, 0.05)){ #fix alpha or in other words R shape parameter
    #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
    #for (i in c(0, 0.5, 1, 2,3,4, 5)){ 
    for (i in c(0)){ 
      
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
      #gompertz.random<-rgompertz(alpha,beta,100)
      #x<- gompertz.random + gaussian
      
      
      lifespan<- gompertz.random +gaussian
      LL_wei <- function (param,lifespan)
        
      {   
        my.data<-lifespan[!is.na(lifespan)]
        
        lambda<-param[1];v=param[2]
        if( lambda< 0 ) { lambda= 1E-10 }
        
        l.x<-log(my.data)
        l.x<-l.x[complete.cases(l.x)]
        power.x<-(my.data/lambda)^v
        power.x<-power.x[complete.cases(power.x)]
        
        
        wei<- N*(log((v/lambda)^v))+(v-1)*sum(l.x)-sum(power.x) 
        #print(param)
        return(-wei)}  #Weibull
      
      max_LL_wei<- optim ( c(3, 0.03), fn=LL_wei, lifespan=lifespan ,lower=c(1E-10, 1E-5), method="L-BFGS-B");
      max_LL_wei$value
      LWei[[length(LWei)+1]] = max_LL_wei$value
      
      
      
      #beta=0.05; alpha=0.02
      LL_gomp <- function (param,lifespan){
        
        my.data = lifespan[!is.na(lifespan)]
        beta<-param[1];alpha=param[2];
        if( alpha< 0 ) { alpha = 1E-10 }
        gomp<-N*log(beta)+alpha*sum(my.data)+(beta/alpha)*(N-sum(exp(alpha*my.data)))
        #print (param ); #trace the convergence
        return(-gomp)    # because optim seems to maximize 
        
      }  #Gompertz
      
      
      max_LL_gomp<-optim(param<-c(0.034,0.01), fn=LL_gomp, lifespan=lifespan,lower=c(1E-2, 1E-5), method="L-BFGS-B")
      max_LL_gomp$value
      LGomp[[length(LGomp)+1]] = max_LL_gomp$value
      
      # store R and G estimation from optim of ll functions in Gompertz
      LLG.par[[length(LLG.par)+1]] =max_LL_gomp$par[1]
      LLR.par[[length(LLR.par)+1]]=max_LL_gomp$par[2]
      
      delta.likelihood.wei<-(max_LL_wei$value-max_LL_gomp$value)
      
      #calculate LL and noise change
      sderr[[length(sderr)+1]] = i       
      Delta_LL[[length(Delta_LL)+1]] = delta.likelihood.wei
      G[[length(G)+1]]=beta
      #switch to alpha.seq when for fixed beta
      R[[length(R)+1]]=alpha
      
      #todo use flexsurv to calculate the LL
      
      #flexsurv only works with positive variables.
      #fix gaussian std to 0
      
      gaussian.flex= rnorm(N, mean = 0, sd=0)
      X.flex= gompertz.random +gaussian.flex
      
      
      fitGomp = flexsurvreg(formula = Surv(X.flex) ~ 1, dist="gompertz")
      fitWei = flexsurvreg(formula = Surv(X.flex) ~ 1, dist="weibull")
      
      
      LWei.flex[[length(LWei.flex)+1]]=fitWei$loglik
      
      LGomp.flex[[length(LGomp.flex)+1]]=fitGomp$loglik
      
      param.Gomp<-fitGomp$res; R.flex<-param.Gomp[1]; G.flex<-param.Gomp[2];
      
      R.flex.estimated[[length(R.flex.estimated)+1]]<-R.flex
      G.flex.estimated[[length(G.flex.estimated)+1]]<-G.flex
      
      param.Wei<-fitWei$res; v.flex<-param.Wei[1]; lambda.flex<-param.Wei[2];
      
      v.flex.estimated[[length(v.flex.estimated)+1]]<-v.flex; 
      lambda.flex.estimated[[length(lambda.flex.estimated)+1]]<-lambda.flex
      
      
      
      
      delta_flexsurv=fitWei$loglik-fitGomp$loglik 
      
      #fitWei$loglik
      
      Delta_LL.flex[[length(Delta_LL.flex)+1]]=delta_flexsurv
      
    }
  }
}

results = data.frame(cbind(sderr), cbind(R),cbind(LLR.par),cbind(R.flex.estimated),cbind(G),cbind(LLG.par),cbind(G.flex.estimated),cbind(Delta_LL) , cbind(Delta_LL.flex),
                     cbind(LWei),cbind(LWei.flex), cbind(LGomp),cbind(LGomp.flex), cbind(MeanLF), cbind(sdLS))



results_mat<-as.matrix(results)



dLL<-unlist(results$Delta_LL )
dLL.flex<-unlist(results$Delta_LL.flex)
LLGomp<-unlist(results$LGomp)
LLGomp.flex<- unlist(results$LGomp.flex)
LLWei<- unlist(results$LWei)
LLWei.flex<- unlist(results$LWei.flex)
simulated.G<-unlist(results$G)
estimated.G.flex<-unlist(results$G.flex.estimated)
simulated.R<-unlist(results$R)
estimated.R.flex<-unlist(results$R.flex.estimated)

estimatedLL.G<-unlist(results$LLG.par)
estimatedLL.R<-unlist(results$LLR.par)


#LLGomp.flex
#delta.output<-data.frame(dLL,dLL.flex,LLGomp,LLGomp.flex,LLWei,LLWei.flex)

summary( lm( dLL~ dLL.flex))
summary( lm( LLGomp ~ LLGomp.flex))
summary( lm( LLWei ~ LLWei.flex))
summary(lm(LLWei~LLGomp))
summary(lm(LLWei.flex~LLGomp.flex))
summary(lm(simulated.G~estimated.G.flex))
summary(lm(simulated.R~estimated.R.flex))
summary(lm(simulated.R~estimatedLL.R))
summary(lm(simulated.G~estimatedLL.G))

plot( simulated.R ~ estimated.R.flex,pch=19)#ylim=c(-400, -500), xlim=c(-2, 500)) 
model<-lm(simulated.R~estimated.R.flex)
par(new=TRUE)
abline(model,col='red')

model<-lm(simulated.R~estimatedLL.R)
plot(simulated.R~estimatedLL.R,pch=19) #ylim=c(-400, -500), xlim=c(-2, 500))
par(new=TRUE)
abline(model,col='red')


plot( simulated.G ~ estimated.G.flex,pch=19)#ylim=c(-400, -500), xlim=c(-2, 500)) 
model<-lm(simulated.G~estimated.G.flex)
par(new=TRUE)
abline(model,col='red')

model<-lm(simulated.G~estimatedLL.G)
plot(simulated.G~estimatedLL.G,pch=19)#ylim=c(-400, -500), xlim=c(-2, 500))
par(new=TRUE)
abline(model,col='red')

