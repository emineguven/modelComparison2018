
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



#find delta log-likelihood for each parameter sets



#generate gompertz random numbers (lifespan) by using inverse CDF
      #prediction
      gompertz.random<-rgompertz(alpha,beta,N)
      average.lifespan=mean(gompertz.random)
      round.lifespan=ceiling(average.lifespan)

results = data.frame(t(c(NA, NA, NA, NA, NA, NA, NA)))
names(results) = c("sderr", "R", "G", "LWei", "LGomp", "MeanLF", "sdLS")

count = 1; 
for (beta in c(0.05, 0.1, 0.15, 0.2, 0.25)){
  for (alpha in c(1E-4, 1E-3, 0.002, 0.005, 0.008, 0.01, 0.05)){ #fix alpha or in other words R shape parameter
    #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
    for (i in c(0, 0.5, 1, 2, 5)){ 
      
      #generate gompertz random numbers (lifespan) by using inverse CDF
      #prediction
      gompertz.random<-rgompertz(alpha,beta,N)
      average.lifespan=mean(gompertz.random)
      sd.lifespan=sd(gompertz.random)
      
      results$sderr[count] = i;
      results$R[count] = alpha
      results$G[count] = beta;
      results$MeanLF[count] = average.lifespan
     
      
      
      
      sd.gaussian=calculate.noise(i)
      gaussian<-rnorm(N, mean = 0, sd=sd.gaussian)
      #standard deviation of lifespan.
      results$sdLS[count]=sd.lifespan
      
      #simulate lifespan plus normal noises
      #observation
      X<- gompertz.random + gaussian
      #Y<-weibull.random+gaussian
      
      #this is to fix a bug problem here power.X is actually the same thing as X=observation
      power.X<-X^v
      power.X<-power.X[complete.cases(power.X)]
      
      
      
      delta.likelihood.wei<- N*log(lambda*v)+(v-1)*log(sum(X))+(-lambda)*sum(power.X)- #Weibull -
        (N*log(beta)+alpha*sum(X)+beta/alpha*sum(1-exp(alpha*X)))                 #Gompertz
      
      #power.Y<-Y^v
      #power.Y<-power.Y[complete.cases(power.Y)]
      
       #delta.likelihood.gomp<- N*log(lambda*v)+(v-1)*log(sum(Y))+(-lambda)*sum(power.Y)- #Weibull -
       # (N*log(beta)+alpha*sum(Y)+beta/alpha*sum(1-exp(alpha*Y)))                 #Gompertz
      
      
      results$LWei[count] =delta.likelihood.wei
      #results$LGomp[count]=deltalikelihood.gomp
      
      #todo use flexsurv to calculate the LL
      
      
    }
  }
}

#generate matrix for noise and alpha
results = mat.data
results.sub = results[results$beta==0.05, ]
names(results.sub) = c("beta",'alpha','noise','dLL')

R.els = unlist( unique(results.sub$alpha))
colnum = length(R.els)

tmp = unlist( unique(results.sub$noise))
noise.els = tmp[order(tmp)]
rownum = length(noise.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = noise.els
colnames(mat) = R.els

for( i in 1:rownum){ #rows for noise
  for( j in 1:colum) #cols for alpha
    mat[i,j] = results$dLL[results$alpha == R.els[i] & results$noise== noise.els[j]]
  
}


#todo countour map or heatmap



#put the noise into data frame
df.noise=data.frame(cbind(noise))

#put the LL change into data frame
df.delta=data.frame(cbind(delta.LL))  

#put everything in one data frame
mat.data<-data.frame(G=cbind(beta.seq),R=cbind(alpha.seq),noise=cbind(noise),delta=cbind(delta.LL))
mat.data


