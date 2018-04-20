

require(flexsurv)
require(gplots)




gamma=5 #shape
theta=25   #scale
#test G and R in nested for loops
R=0.01
G=0.01
#R should be in [0, 0.05], and G should be [0.05, 0.3].  
N=2000 # population size





# random Gompertz random generator function

#where alpha and beta are 2 parameters, N is number of population
#generate gompertz random numbers by using inverse CDF
#generate random number with a given distribution of Gompertz
#prediction
#rgompertz = function(N,G, R){
  ### I have a function to generate random-numbers from 
  ### (using the 'inversion method')
  
#  return( (log(1-(G/R*log(1-runif(N)))))/G)
#}

rweibull <- function(N,gamma,theta){
  ### I have a function to generate random-numbers from 
  ### (using the 'inversion method')
  
  X <- theta*(-log(1-runif(N)))^(1/gamma)
  return(X)
}



set.seed(123)
#generate gompertz random numbers (lifespan) 
#prediction
weibull.random<-rweibull(N,gamma,theta)
average.lifespan=mean(weibull.random)

scale=0.01

#generate gaussion random numbers 
set.seed(123)
gaussian<-rnorm(N, mean = 0, sd=sd(weibull.random)*scale)


lifespan = weibull.random+ gaussian
summary(lifespan)

lifespan[lifespan < 0] <- 0
lifespan<-floor(lifespan+0.5)

lifespan <- lifespan[ lifespan != 0 ]




calculate.s = function( lifespan ){
  myData = sort( lifespan[!is.na(lifespan)] );
  tmpGC = table( myData )
  for( i in 2:length(tmpGC)) {
    tmpGC[i] = tmpGC[i-1] + tmpGC[i]        }    
  tot = length(myData)
  tmpGC = tmpGC / tot; 
  s = 1 - tmpGC
  #list( s=s, t=unique(my.data));
  ret = data.frame( cbind(s, unique(myData)));
  names(ret) = c("s", "t");
  ret;
}


GC = calculate.s(lifespan)
plot(GC$s ~ GC$t)



#3)calculate the mortality rates he mortality rate change over time, and then plot log(mortality rate) ~ time.
#Gompertz model should give a linear form.  
#For Weibull model, log(moretality rate) ~ log(time) give the linear form.

#HQin's calculate mortality rate function to calculate the rate of mortality over time
calculate.mortality.rate = function( lifespan ){
  GC = calculate.s(lifespan)
  GC$ds=NA; GC$dt=NA
  #first point
  GC$dt[1] = GC$t[2]
  GC$ds[1] = 1 - GC$s[1]
  GC$mortality.rate[1] = GC$ds[1] / GC$dt[1]
  
  for( j in 2:length(GC[,1])) {
    GC$ds[j] =  GC$s[j-1] - GC$s[j] 
    GC$dt[j] = -GC$t[j-1] + GC$t[j]
    GC$mortality.rate[j] = GC$ds[j] / ( GC$s[j] * GC$dt[j])
  }
  return(GC)
} #end of calculate.mortality.rate()

GC = calculate.mortality.rate(lifespan)



GC = calculate.s(round( lifespan, digits=1))
head(GC)
GC$ds=NA; GC$dt=NA
GC$dt[1] = GC$t[2] #20130321 correct a bug GC$s -> GC$t
GC$ds[1] = 1 - GC$s[1]
GC$mortality.rate[1] = GC$ds[1] / GC$dt[1]

for( j in 2:length(GC[,1])) {
  GC$ds[j] =  GC$s[j-1] - GC$s[j] 
  GC$dt[j] = -GC$t[j-1] + GC$t[j]
  GC$mortality.rate[j] = GC$ds[j] / ( GC$s[j] * GC$dt[j])
}
plot( GC$s ~ GC$t)
plot( GC$mortality.rate ~ GC$t, typ='l', log='y' )

summary(lm(log10(GC$mortality.rate[1:(length(GC$mortality.rate)-1)]) ~ GC$t[1:(length(GC$t)-1)]))
#then plot log(mortality rate) ~ time.
#For Weibull model, log(moretality rate) ~ log(time) give the linear form.

pdf(paste("plots/","Weibull.semi.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ GC$t, type='l') #linear for Gompertz, semi-log plot
text(10,-1,"R2= 0.91")
dev.off()


summary(lm(log10(GC$mortality.rate[1:(length(GC$mortality.rate)-1)]) ~ log10(GC$t[1:(length(GC$t)-1)])))
pdf(paste("plots/","Weibull.log.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ log10(GC$t), type='l' ,main="my weibull random generator" ) #linear for Weibull, log-log plot
text(0.75,-1,"R2=0.97")
dev.off()








# initialize variables 


#generate gompertz random numbers (lifespan) 
#prediction
weibull.random<-rweibull(N,gamma,theta)
average.lifespan=mean(weibull.random)

#create a data frame for fit parameters
fit_names = c( 'sd.rweibull','current.seed','scale','sderr','Delta_LL','gamma','theta','LWei','LGomp','MeanLF','Delta_LL.flex','LWei.flex','LGomp.flex','G.flex.estimated','R.flex.estimated','LLtheta.par','LLgamma.par','gamma.flex.estimated','theta.flex.estimated')
fit= data.frame(matrix(, nrow=length(fit_names), ncol=N)) #set up a skeleton table

names(fit) = c( "sd.rweibull","current.seed","scale","sderr","Delta_LL","theta","gamma","LWei","LGomp","MeanLF","Delta_LL.flex","LWei.flex","LGomp.flex","G.flex.estimated","R.flex.estimated","LLtheta.par","LLgamma.par","gamma.flex.estimated","theta.flex.estimated")




## introduce Weibull model log-likelihood function


#log likelihood function for the Weibull model
#s = exp(-( theta*x)^gamma)
#f = gamma * theta^gamma *x^(gamma-1); 

LL_wei<-function(param,y){
  theta<- exp(param[1]) #take exponential to avoid NaNs when taking log(theta)
  gamma<- exp(param[2]) # avoid NaNs when taking log(gamma)
  
  data=lifespan[!is.na(lifespan)]
  #log_s = -( theta*x)^gamma
  #log_f = log(gamma)+ gamma*log(theta)+ (gamma-1)*log(x) ; 
  #w.lh = sum(log_s)  + sum(log_f);
  
  w.lh<- sum(log(gamma) + gamma*log(theta) + (gamma-1)*log(data)) - sum(
    (theta*data)^gamma)
  
  return(-w.lh)
}
# take log(param) since you take exponential above to avoid NaN values above


##introduce Gompertz model log-likelihood function



#log likelihood function of gompertz model
#s = exp( (R/G) *(1 - exp(G* data)) )  
#f = R * exp( G * data ); 

LL_gomp<- function( param, y ) {
  
  G = exp(param[1]); R = exp(param[2]); 
  data = lifespan[!is.na(lifespan)];
  #log_s = (R/G) *(1 - exp(G* data))
  
  #log_f = log(R) +  G * data ; 
  
  g.lh = sum((R/G) *(1 - exp(G* data)))  + sum(log(R) +  G * data );
  
  return(- g.lh) 
}



## simulate for parameters R, G and scale=i to search effect of noise on delta likl





for( seed in c(12345,20160711, -1881, 9999.1234,300045,50758,-10000,74562,-92345,25434)) {
  set.seed(seed)
  current.seed<-seed
  for (gamma_in in c(1,2,3,4,5,6,10,15,30)){  #gamma_in is shape parameter
    for (theta_in in c(25,26,27,28,29,30,35,40,50)){ # theta_in is scale param.
      
      for (i in c(seq(0,1,by=0.1), seq(1.5,3,by=0.5))){ 
        
        
        #generate gompertz random numbers (lifespan) 
        #prediction
       # gamma_in=2
       # theta_in=1000
       # i=0
        weibull.random<-rweibull(N, gamma_in, theta_in)
        
        lifespan= weibull.random + rnorm(N, mean=0, sd=sd(weibull.random)*i)
        
        lifespan[lifespan < 0] <- 0
        lifespan<-floor(lifespan+0.5)
        lifespan <- lifespan[ lifespan != 0 ]
        
        summary(lifespan)
        
        scale=i
        sdrweibull<-sd(weibull.random)
        
        average.lifespan = mean(lifespan)
        #store average.lifespan into MeanLF list
        MeanLF = average.lifespan
        
        #Log likelihood function for the Weibull model
        
        
        #calculate noise change
        
        sderr = sd(weibull.random)*i   
        
        #G=G_in
        #R=R_in
        
        theta = theta_in
        gamma = gamma_in
        
        weib=optim(log(c(5,0.02)),LL_wei,y=lifespan)
        LWei = - weib$value
        
        # store gamma and theta estimation from optim of likl functions in Gompertz
        LLtheta.par =weib$par[1]
        LLgamma.par=weib$par[2]
        
        gomp<-optim(param<-log(c(G,R)), fn=LL_gomp, y=lifespan)
        LGomp = -gomp$value
        
        # store R and G estimation from optim of likl functions in Gompertz
        LLG.par =gomp$par[1]
        LLR.par=gomp$par[2]
        
        delta.likelihood<- -weib$value-(-gomp$value)
        
        
        #Delta_LL[[length(Delta_LL)+1]] = delta.likelihood.wei
        
        
        #flexsurv to calculate the log-likelihood value for both models
        #flexsurv only works with positive variables.
        
        fitGomp = flexsurvreg(formula = Surv(lifespan) ~ 1, dist="gompertz")
        fitWei = flexsurvreg(formula = Surv(lifespan) ~ 1, dist="weibull")
        
        
        LWei.flex=fitWei$loglik
        
        LGomp.flex=fitGomp$loglik
        
        param.Gomp<-fitGomp$res; R.flex<-param.Gomp[2]; G.flex<-param.Gomp[1];
        
        R.flex.estimated<-R.flex
        G.flex.estimated<-G.flex
        
        param.Wei<-fitWei$res; theta.flex<-param.Wei[1]; gamma.flex<-param.Wei[2];
        
        gamma.flex.estimated<-gamma.flex; 
        theta.flex.estimated<-theta.flex
        
        delta_flexsurv=fitWei$loglik-fitGomp$loglik 
        
        #fitWei$loglik
        
        Delta_LL.flex=delta_flexsurv
        
        
        #sim_names = c( "sd.rweibull","scale","sderr","Delta_LL","G","R","LWei","LGomp","MeanLF","Delta_LL.flex","LWei.flex","LGomp.flex","G.flex.estimated","R.flex.estimated","LLtheta.par","LLgamma.par","gamma.flex.estimated","theta.flex.estimated")
        fit = rbind(fit,c(sdrweibull,current.seed,scale,sderr,delta.likelihood,theta,gamma,LWei,LGomp,MeanLF,Delta_LL.flex,LWei.flex,LGomp.flex,G.flex.estimated,R.flex.estimated,LLtheta.par,LLgamma.par,gamma.flex.estimated,theta.flex.estimated))
        #write.csv(fit, file="Results.csv", row.names=F)
        fit= fit[!is.na(fit[,1]), ]
      }
    }
  }
}

fit<- fit[!is.na(names(fit))]
write.csv(fit, file="ResultsWei.csv", row.names=F)

fit<-read.csv(file="ResultsWei.csv")
#fit[c(1:length(fit[fit$current.seed==12345,5])),]

summary( lm( fit$Delta_LL~ fit$Delta_LL.flex))
summary( lm( fit$LGomp ~ fit$LGomp.flex))
summary( lm( fit$LWei ~ fit$LWei.flex))

summary(lm(fit$theta~fit$theta.flex.estimated))
summary(lm(fit$gamma~fit$gamma.flex.estimated))


pdf(paste("plots/","simulated.theta.vs.estimated.theta.flex_batch.pdf", sep=''), width=5, height=5)
plot(fit$theta~fit$theta.flex.estimated)
dev.off()

pdf(paste("plots/","simulated.gamma.vs.estimated.gamma.flex_batch.pdf", sep=''), width=5, height=5)
plot(fit$gamma~fit$gamma.flex.estimated)
dev.off()




Delta_LL=list()
for (i in 1:length(fit[fit$current.seed==12345,5])){
  Delta_LL[length(Delta_LL)+1]<-(fit[fit$current.seed==12345,5][i]+fit[fit$current.seed==20160711,5][i]+
                                   fit[fit$current.seed==-1881,5][i]+fit[fit$current.seed==9999.1234,5][i]+fit[fit$current.seed== 300045,5][i]
                                 +fit[fit$current.seed==50758,5][i]+fit[fit$current.seed==-10000,5][i]+fit[fit$current.seed== 74562,5][i]+
                                   fit[fit$current.seed== -92345,5][i]+fit[fit$current.seed==25434,5][i])/10
  
}

fit_sub<-fit[1:length(fit[fit$current.seed==12345,5]),]

Delta_LL<-unlist(Delta_LL)
fit_sub$Delta_LL<-Delta_LL

#heatmap for fixed G








#heatmap for fixed G




n.col=256 # number of colors
cm = redblue(n.col) # red-white-blue colormap
mmx = min(abs(min(fit_sub$Delta_LL)),abs(max(fit_sub$Delta_LL))) # find min value, or you can set to a number
colbr <- c(seq(-mmx/2,mmx/2, len=length(cm)+1)) # vector of symmetric breaks

gamma.els = unlist( unique(fit_sub$gamma))
colnum = length(gamma.els)

tmp = unlist( unique(fit_sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = scale.els
colnames(mat) = gamma.els
for (k  in c(25,26,27,28,29,30,35,40,50)){ 
  
 
  
  data = fit_sub[fit_sub[,6]==k, 5]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  rownames(heat_mat) = scale.els
  colnames(heat_mat) = gamma.els
  #rownames(heat_mat) <- scale.els
  #colnames(heat_mat) <- gamma.els
  
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plotsWei/heatplot_zero_noise_G",k,".jpeg",sep="") 
  #par(mfrow=c(2,1)) 
  pdf(paste("plots/",k, "_fixed_theta.pdf", sep=''))
  
  if (k==25){
  
  heatmap.2(heat_mat, scale="none",col = cm, breaks=colbr, notecol="black",
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=TRUE,symkey=F, density.info="none",
            xlab= expression( ~gamma ~"parameters"),
            ylab= expression( ~eta ~"parameters"),margins=c(4,4), main = bquote(paste(~delta[LL], " for fixed", ~ theta==.(k))),par(cex.main=.5),srtCol=315,
            srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,10),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,15),
            sepcolor="black")
  
  }else{
  dev.off()
  pdf(paste("plots/",k, "_fixed_theta.pdf", sep=''))
  heatmap.2(heat_mat, scale="none",col = cm, breaks=colbr, notecol="black",
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=FALSE,symkey=F, density.info="histogram",
            xlab= expression( ~gamma ~"parameters"),
            ylab= expression( ~eta ~"parameters"),margins=c(4,4), main = bquote(paste(~delta[LL], " for fixed", ~ theta==.(k))),
            par(cex.main=.5),srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,10),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,15),
            sepcolor="black")
  
  
  #jpeg(paste("plotsWei/",k, ".fixed_G.dLLhist.jpg", sep=''))
  #hist(heat_mat,  main = bquote(paste("R vs. scale" ~ G==.(k))))
  #dev.off()
  
}

dev.off()
}

# heatmap for fixed scale


theta.els = unlist( unique(fit_sub$theta))
colnum = length(theta.els)

gamma.els=unlist(unique(fit_sub$gamma))
rownum = length(gamma.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = gamma.els
colnames(mat) = theta.els


for (eta in c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
              1.5, 2.0, 2.5,3.0)){
  
  data = fit_sub[ fit_sub[, 3]== eta, 5]
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  
  
  rownames(heat_mat) = gamma.els
  colnames(heat_mat) = theta.els
  
  heat_mat<-t(heat_mat)
  
  #rownames(heat_mat) <- gamma.els
  #colnames(heat_mat) <- theta.els
  
  library(gplots)
  
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plotsWei/heatplot_zero_noise_G",k,".jpeg",sep="") 
  #par(mfrow=c(2,1)) 
  pdf(paste("plots/",eta, "_fixed_eta.pdf", sep=''))
  
  if(eta==0){
  
  heatmap.2(heat_mat,col = cm, breaks=colbr,scale="none", density.info="none",notecol="black",
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=TRUE,
            xlab= expression( ~gamma ~"parameters"), 
          ylab =expression(~theta ~"parameters"),margins=c(4,4),main = bquote(paste(~delta[LL], " for fixed", ~ eta==.(eta))),
          par(cex.main=.5),srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
          sepwidth=c(0.001,0.001),
          colsep=c(0,9),
          #sepwidth=c(0.05,0.05),
          rowsep=c(0,9),
          sepcolor="black")
  
  }else{
    dev.off()
    pdf(paste("plots/",eta, "_fixed_eta.pdf", sep=''))
    heatmap.2(heat_mat,col = cm, breaks=colbr,scale="none", density.info="histogram",notecol="black",
              dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=FALSE,
              xlab= expression( ~gamma ~"parameters"), 
              ylab =expression(~theta ~"parameters"),margins=c(4,4),main = bquote(paste(~delta[LL], " for fixed", ~ eta==.(eta))),
              par(cex.main=.5),srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
              sepwidth=c(0.001,0.001),
              colsep=c(0,9),
              #sepwidth=c(0.05,0.05),
              rowsep=c(0,9),
              sepcolor="black")
    
  }
    dev.off()
} 
  #jpeg(paste("plotsWei/",n, ".fixed_scale.dLLhist.jpg", sep=''))
  #hist(heat_mat, main = bquote(paste("G~R for fixed" ~ scale==.(n))))
  
  #dev.off()
  
}


# heatmap for fixed R


theta.els = unlist( unique(fit_sub$theta))
colnum = length(theta.els)

tmp = unlist( unique(fit_sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
colnames(mat) = theta.els
rownames(mat) = scale.els



for (j in  c(1,2,3,4,5,6,10,15,30)){
  
  data = fit_sub[fit_sub[,7]==j, 5]
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  #rownames(heat_mat) <- scale.els
  
  
  #colnames(heat_mat) <- theta.els
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  colnames(heat_mat) = theta.els
  rownames(heat_mat) = scale.els
  
  #paste(file = "~/github/model.comparison/plotsWei/heatplot_zero_noise_G",k,".jpeg",sep="") 
  
  pdf(paste("plots/",j, "_fixed_gamma.pdf", sep=''))
  
  if (j==1){
  
  heatmap.2(heat_mat,col = cm, breaks=colbr, scale="none", notecol="black", 
            dendrogram='none', density.info="none",Rowv=F, Colv=F,trace='none', key = T,
            xlab     = expression(~theta ~"parameters"),
            ylab     = expression(~eta ~"parameters"),margins=c(4,4),main = bquote(paste(~delta[LL], " for fixed", ~ gamma==.(j))),
            par(cex.main=.5),srtCol=315, srtRow=315,adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,10),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,15),
            sepcolor="black")
  } else{
    dev.off()
    pdf(paste("plots/",j, "_fixed_gamma.pdf", sep=''))
    heatmap.2(heat_mat,col = cm, breaks=colbr, scale="none", notecol="black", 
              dendrogram='none', density.info="histogram",Rowv=F, Colv=F,trace='none', key =FALSE,
              xlab     = expression(~theta ~"parameters"),
              ylab     = expression(~eta ~"parameters"),margins=c(4,4),main = bquote(paste(~delta[LL], " for fixed", ~ gamma==.(j))),
              par(cex.main=.5),srtCol=315, srtRow=315,adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
              sepwidth=c(0.001,0.001),
              colsep=c(0,10),
              #sepwidth=c(0.05,0.05),
              rowsep=c(0,15),
              sepcolor="black")
  }
  
  
  dev.off()
  
  
  
  #jpeg(paste("plotsWei/",j, ".fixed_R.dLLhist.jpg", sep=''))
  #hist(heat_mat,  main = bquote(paste("G vs. scale" ~ R==.(j))))
  #dev.off()
}









