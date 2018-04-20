
#EG 20160403 time 12:56 pm
#EG is still trying to prove 
#if the log likelihood function difference of Weibull and Gompertz is always positive
# this is just to simulate if it is really a fact we usually get the positive values or not.
# by simulation I found when x>0, beta and lambda>0.3 
#the difference of log likelihood functions is always positive

#generate postive random variables for gompertz and weibull
#x has to be positive
N = 100
#x = rnorm(n=N, mean=20, sd=3)
#x <- rnorm(N, mean = 0, sd = 1)

#log normal random numbers
#popSize = 1000
x = rlnorm(n=N, mean=3, sd=1)


#difference function of loglikelihood function of gompertz and weibull p.d.fs
# test if L(Weibull,X)>L(Gompertz) for parameters beta,lambda>0.3
# I found between 0 and 0.3 log likelihood difference is negative 
# difference loglikelihood is always negative for 
#restrict beta and lambda to be greater than 0.3 


find.bigger.log = function(lambda,v,beta,alpha){
  
  #unknown parameter coeff ,k tests for lambda and v,z tests for alpha
  #adjust z as negative or positive for alpha
  k=0.1
  
  
  N^2*(k*v-1)*log(k*(lambda+v))*(log(sum(x))-lambda*sum(x)^k*v)-
    (N*k*beta*log(k*beta)*sum(x)*(1-exp(k*alpha*sum(x))));
  
 
  }

i = 10
difference.log=array(NA, dim = c(i,i,i,i),dimnames=list( d1=c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10"),
                                    d2=c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10"),
                                    d3=c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10"),
                                    d4=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10" )))
                                    
                                    #"L11","L12","L13","L14","L15","L16","L17","L18","L19","L20",
                                    #"L21","L22","L23","L24","L25","L26","L27","L28","L29","L30"), 
                                    #"v11","v12","v13","v14","v15","v16","v17","v18","v19","v20",
                                    # "v21","v22","v23","v24","v25","v26","v27","v28","v29","v30"), 
                                    #"a11","a12","a13","a14","a15","a16","a17","a18","a19","a20",
                                    #"a21","a22","a23","a24","a25","a26","a27","a28","a29","a30"),
                                         #"b11","b12","b13","b14","b15","b16","b17","b18","b19","b20",
                                         #"b21","b22","b23","b24","b25","b26","b27","b28","b29","b30")))

for (lambda in c(1:i)){
  
  for (v in c(1:i)){
    
      for (beta in c(1:i)){ 
        
        for (alpha in c(1:i)){
    
    #the log likelihood difference is stored in output matrix
    difference.log[lambda,v,beta,alpha]<-find.bigger.log(lambda,v,beta,alpha)
    
  }
  
}}}



#after entries [10,10] the log difference gives lambda,v,alpha =0.3
#k*10:100 = (0.3,0.31,0.32....,1)
#count how many negative values
#count.negative<-sum(difference.log[10:100,10:100]<0)

#my.data<-as.data.frame(as.table(difference.log)) 

write.csv(difference.log, "test_difference.zero&one.csv")

mat_data <- read.csv("test_difference.zero&one.csv", comment.char="#")

rnames <- mat_data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(mat_data[,2:ncol(mat_data)]) # transform column 2-: into a matrix
rownames(mat_data) <- rnames      # assign row names

mat_data[is.infinite(mat_data)] <- NA

#mydata <- apply( mat_data, 2, as.numeric )
#mydata[is.na(mydata)] <- 0

library(RColorBrewer)

 rc <- rainbow(nrow(mat_data), start = 0, end = .3)
 cc <- rainbow(ncol(mat_data), start = 0, end = .3)
 heatmap(mat_data, Rowv=NA, Colv=NA, col = cm.colors(256), scale = "column", margins=c(5,10), revC=T)





#dev.off()    