

require(flexsurv)
require(gplots)




gamma=0.01
theta=0.25
#test G and R in nested for loops
R=0.001
G=0.2
#R should be in [0, 0.05], and G should be [0.05, 0.3].  
N=2000 # population size





# random Gompertz random generator function

#where alpha and beta are 2 parameters, N is number of population
#generate gompertz random numbers by using inverse CDF
#generate random number with a given distribution of Gompertz
#prediction
rgompertz = function(N,G, R){
  ### I have a function to generate random-numbers from 
  ### (using the 'inversion method')
  
  return( (log(1-(G/R*log(1-runif(N)))))/G)
}


set.seed(123)
#generate gompertz random numbers (lifespan) 
#prediction
gompertz.random<-rgompertz(N,0.05,0.001)
average.lifespan=mean(gompertz.random)

scale=1

#generate gaussion random numbers 
set.seed(123)
gaussian<-rnorm(N, mean = 0, sd=sd(gompertz.random)*scale)


lifespan = gompertz.random+ gaussian
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

pdf(paste("plots/","Gompertz.semi.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ GC$t, type='l') #linear for Gompertz, semi-log plot
text(10,0,"R2= 0.89")
dev.off()


summary(lm(log10(GC$mortality.rate[1:(length(GC$mortality.rate)-1)]) ~ log10(GC$t[1:(length(GC$t)-1)])))
pdf(paste("plots/","Weibull.log.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ log10(GC$t), type='l'  ) #linear for Weibull, log-log plot
text(0.75,0,"R2=0.79")
dev.off()








# initialize variables 


#generate gompertz random numbers (lifespan) 
#prediction
gompertz.random<-rgompertz(N,G,R)
average.lifespan=mean(gompertz.random)

#create a data frame for fit parameters
fit_names = c( 'sd.rgompertz','current.seed','scale','sderr','Delta_LL','G','R','LWei','LGomp','MeanLF','Delta_LL.flex','LWei.flex','LGomp.flex','G.flex.estimated','R.flex.estimated','LLG.par','LLR.par','gamma.flex.estimated','theta.flex.estimated')
fit= data.frame(matrix(, nrow=length(fit_names), ncol=N)) #set up a skeleton table

names(fit) = c( "sd.rgompertz","current.seed","scale","sderr","Delta_LL","G","R","LWei","LGomp","MeanLF","Delta_LL.flex","LWei.flex","LGomp.flex","G.flex.estimated","R.flex.estimated","LLG.par","LLR.par","gamma.flex.estimated","theta.flex.estimated")




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





for( seed in c(12345, 20160711, -1881, 9999.1234,300045,50758,-10000,74562,-92345,25434)) {
  set.seed(seed)
  current.seed<-seed
  for (R_in in c(0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.04)){ #fix alpha or in other words R shape parameter
    for (G_in in c(0.05,0.08, 0.1,0.12,0.15,0.17, 0.2, 0.25)){
      
      for (i in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,2,3)){ 
        
        
        #generate gompertz random numbers (lifespan) 
        #prediction
        
        gompertz.random<-rgompertz(N, G_in, R_in)
        
        lifespan= gompertz.random + rnorm(N, mean=0, sd=sd(gompertz.random)*i)
        
        lifespan[lifespan < 0] <- 0
        lifespan<-floor(lifespan+0.5)
        lifespan <- lifespan[ lifespan != 0 ]
        
        summary(lifespan)
        
        scale=i
        sdrgompertz<-sd(gompertz.random)
        
        average.lifespan=mean(lifespan)
        #store average.lifespan into MeanLF list
        MeanLF=average.lifespan
        
        #Log likelihood function for the Weibull model
        
        
        #calculate noise change
        
        sderr = sd(gompertz.random)*i   
        
        G=G_in
        R=R_in
        
        weib=optim(log(c(3,0.03)),LL_wei,y=lifespan)
        LWei = - weib$value
        
        gomp<-optim(param<-log(c(G_in,R_in)), fn=LL_gomp, y=lifespan)
        LGomp= -gomp$value
        
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
        
        param.Wei<-fitWei$res; gamma.flex<-param.Wei[1]; theta.flex<-param.Wei[2];
        
        gamma.flex.estimated<-gamma.flex; 
        theta.flex.estimated<-theta.flex
        
        delta_flexsurv=fitWei$loglik-fitGomp$loglik 
        
        #fitWei$loglik
        
        Delta_LL.flex=delta_flexsurv
        
        
        #sim_names = c( "sd.rgompertz","scale","sderr","Delta_LL","G","R","LWei","LGomp","MeanLF","Delta_LL.flex","LWei.flex","LGomp.flex","G.flex.estimated","R.flex.estimated","LLG.par","LLR.par","gamma.flex.estimated","theta.flex.estimated")
        fit = rbind(fit,c(sdrgompertz,current.seed,scale,sderr,delta.likelihood,G,R,LWei,LGomp,MeanLF,Delta_LL.flex,LWei.flex,LGomp.flex,G.flex.estimated,R.flex.estimated,LLG.par,LLR.par,gamma.flex.estimated,theta.flex.estimated))
        #write.csv(fit, file="Results.csv", row.names=F)
        fit= fit[!is.na(fit[,1]), ]
      }
    }
  }
}

fit<- fit[!is.na(names(fit))]
write.csv(fit, file="Results.csv", row.names=F)

fit<-read.csv(file="Results.csv")


summary( lm( fit$Delta_LL~ fit$Delta_LL.flex))
summary( lm( fit$LGomp ~ fit$LGomp.flex))
summary( lm( fit$LWei ~ fit$LWei.flex))

summary(lm(fit$G~fit$G.flex.estimated))
summary(lm(fit$R~fit$R.flex.estimated))


pdf(paste("plots/","simulated.G.vs.estimated.G.flex_batch.pdf", sep=''), width=5, height=5)
plot(fit$G~fit$G.flex.estimated)
dev.off()

pdf(paste("plots/","simulated.R.vs.estimated.R.flex_batch.pdf", sep=''), width=5, height=5)
plot(fit$R~fit$R.flex.estimated)
dev.off()




#fit<-fit[c(1:length(fit[fit$current.seed==12345,5])),]

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

pdf(paste("plots/", "Histogram_simulated_data.pdf", sep=''))
par(mfrow=c(3,1)) 
hist(as.numeric(fit_sub$MeanLF),xlab="Mean Lifespan",main="Simulated Data",
     breaks=length(as.numeric(fit_sub$MeanLF))/5,col="gray",xlim=c(0,100),ylim=c(0,60))

mtext('Total experiments=704',line=-2,cex=0.5)
box(lty='1373',col='blue')

hist(as.numeric(fit_sub$G.flex.estimated),xlab="G shape parameters",main="",
     breaks=length(as.numeric(fit_sub$G.flex.estimated))/5,col="blue",xlim=c(0,0.35),ylim=c(0,40))

box(lty='1373',col='blue')

hist(log(as.numeric(fit_sub$R.flex.estimated)),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(fit_sub$R.flex.estimated))/5,col="green",xlim=c(-7.15,-2.9),ylim=c(0,40))
box(lty='1373',col='blue')
dev.off()


pdf(paste("plots/", "Histogram_Delta_LL_simulated_data.pdf", sep=''))

fifty.percent1<-quantile(fit_sub$Delta_LL, c(.50)) 


hist(fit_sub$Delta_LL,breaks=length(fit_sub$Delta_LL)/5,col="blue",xlim=c(-150,50),ylim=c(0,20),xlab=expression(paste( ~delta[LL])),main="Simulated Data Delta_LL Histogram")


mtext('Total experiments=704',line=-3,at=-100,side=3,cex=0.8)
mtext('median=-23.29',line=-4,at=-100,side=3,cex=0.8)
mtext('mean=-34.55',line=-5,at=-100,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent1, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-12,at=15,cex=0.8)

dev.off()


pdf(paste("plots/", "Histogram_Delta_LL_flexsurv_simulated_data.pdf", sep=''))

fifty.percent2<-quantile(fit_sub$Delta_LL.flex, c(.50)) 


hist(fit_sub$Delta_LL.flex,breaks=length(fit_sub$Delta_LL.flex)/5,col="blue",xlim=c(-160,50),ylim=c(0,40),xlab=expression(paste( ~delta[LL])),main="Simulated Data Delta_LL.flexsurv Histogram")


mtext('Total experiments=704',line=-3,at=-100,side=3,cex=0.8)
mtext('median=-23.81',line=-4,at=-100,side=3,cex=0.8)
mtext('mean=-34.09',line=-5,at=-100,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent2, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-11.5,at=30,cex=0.8)

dev.off()


gompShape<-as.numeric(fit_sub$G.flex.estimated)
gompRate<-as.numeric(fit_sub$R.flex.estimated)

CV_gompShape<-sd(gompShape)/mean(gompShape)
CV_gompRate<-sd(gompRate)/mean(gompRate)

sd.Rls<-as.numeric(fit_sub$MeanLF)
avg.Rls<-as.numeric(fit_sub$sderr)

sd.Rls<-sd.Rls[!is.na(sd.Rls)]
avg.Rls<-avg.Rls[!is.na(avg.Rls)]

CV_sd.Rls<- sd(sd.Rls)/mean(sd.Rls)
CV_avg.Rls<-sd(avg.Rls)/mean(avg.Rls)

CV_names = c('CV of mean(LS)','CV of sd(LS)', 'CV of gompShape','CV of gompRate')
CV_results = c( CV_avg.Rls,CV_sd.Rls, CV_gompShape,CV_gompRate)
CV.df<-data.frame(CV_names,CV_results)
CV.df



n.col=256 # number of colors
cm = redblue(n.col) # red-white-blue colormap
mmx = min(abs(min(fit_sub$Delta_LL)),abs(max(fit_sub$Delta_LL))) # find min value, or you can set to a number
colbr <- c(seq(-mmx/2,mmx/2, len=length(cm)+1)) # vector of symmetric breaks

R.els = unlist( unique(fit_sub$R))
colnum = length(R.els)

tmp = unlist( unique(fit_sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = scale.els
colnames(mat) = R.els
for (k in c(0.05,0.08, 0.1,0.12,0.15,0.17, 0.2, 0.25)){
  
  data = fit_sub[fit_sub[,6]==k, 5]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  
  rownames(heat_mat) <- scale.els
  colnames(heat_mat) <- R.els
  
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  #par(mfrow=c(2,1)) 
  
  pdf(paste("plots/",k, "_hist_.pdf", sep=''))
  
  hist(heat_mat,xlab=expression(paste( ~delta[LL])),main=expression(paste( ~delta[LL]," of Gompertz lifespans")),
       ylim=c(0,12),breaks=length(heat_mat),col="grey")
  dev.off()
  
  pdf(paste("plots/",k, "_fixed_G.pdf", sep=''))
  
  if (k == 0.05){
  
  
  #cellnote=hM,
  heatmap.2(heat_mat, scale="none",col = cm, breaks=colbr, notecol="black",
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=TRUE,symkey=F, density.info="none",
            xlab     = "R parameters",
            ylab     = expression( ~eta ~"parameters"), margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed", ~ G==.(k))),par(cex.main=.5),
            srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,8),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,11),
            sepcolor="black")
 
  } else {
    dev.off()
    pdf(paste("plots/",k, "_fixed_G.pdf", sep=''))
    heatmap.2(heat_mat, scale="none",col = cm, breaks=colbr, notecol="black",
              dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=F,symkey=F, density.info="histogram",
              xlab     = "R parameters",
              ylab     = expression( ~eta ~"parameters"), margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed", ~ G==.(k))),par(cex.main=.5),
              srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
              sepwidth=c(0.001,0.001),
              colsep=c(0,8),
              #sepwidth=c(0.05,0.05),
              rowsep=c(0,11),
              sepcolor="black")
    
  }
  dev.off()
} 

# heatmap for fixed scale


G.els = unlist( unique(fit_sub$G))
colnum = length(G.els)

R.els=unlist(unique(fit_sub$R))
rownum = length(R.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = R.els
colnames(mat) = G.els


for (eta in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,2,3) ){
  
  data = fit_sub[fit_sub[,3]==0, 5]
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  heat_mat<-t(heat_mat)
  
  
  
  rownames(heat_mat) <- R.els
  colnames(heat_mat) <- G.els
  
  library(gplots)
  
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  #par(mfrow=c(2,1)) 
  
  pdf(paste("plots/",eta, "_hist_.pdf", sep=''))
  
  hist(heat_mat,xlab=expression(paste( ~delta[LL])),main=expression(paste( ~delta[LL]," of Gompertz lifespans")),
       ylim=c(0,12),breaks=length(heat_mat),col="grey")
  dev.off()
  pdf(paste("plots/",eta, "_fixed_eta.pdf", sep=''))
  
  if (eta==0){
  heatmap.2(heat_mat,col = cm, breaks=colbr,scale="none", density.info="none",notecol="black",
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=TRUE,
            xlab     = "G parameters",
            ylab     = "R parameters", margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed",~eta==.(eta))),par(cex.main=.5),
            srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
  sepwidth=c(0.001,0.001),
  colsep=c(0,8),rowsep=c(0,8),sepcolor="black")
  #dev.off()
  } else{
   dev.off()
    pdf(paste("plots/",eta, "_fixed_eta.pdf", sep=''))
    heatmap.2(heat_mat,col = cm, breaks=colbr,scale="none", density.info="none",notecol="black",
              dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=F,
              xlab     = "G parameters",
              ylab     = "R parameters", margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed", ~eta==.(eta))),par(cex.main=.5),
              srtCol=315, srtRow=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
              sepwidth=c(0.001,0.001),
              colsep=c(0,8),rowsep=c(0,8),sepcolor="black")
    
    
   #dev.off()
  }
  #jpeg(paste("plots/",n, ".fixed_scale.dLLhist.jpg", sep=''))
  #hist(heat_mat, main = bquote(paste("G~R for fixed" ~ scale==.(n))))
  
  dev.off()
  
}


# heatmap for fixed R


G.els = unlist( unique(fit_sub$G))
colnum = length(G.els)

tmp = unlist( unique(fit_sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
colnames(mat) = G.els
rownames(mat) = scale.els



for (j in c(0.001,0.002,0.003,0.005,0.01,0.02,0.03,0.04) ){
  
  data = fit_sub[fit_sub[,7]==j, 5]
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  rownames(heat_mat) <- scale.els
  
  
  colnames(heat_mat) <- G.els
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  
  
  pdf(paste("plots/",j, "_hist_.pdf", sep=''))
  
  hist(heat_mat,xlab=expression(paste( ~delta[LL])),main=expression(paste( ~delta[LL]," of Gompertz lifespans")),
       ylim=c(0,12),breaks=length(heat_mat),col="grey")
  dev.off()
  pdf(paste("plots/",j, "_fixed_R.pdf", sep=''))
  
  if(j==0.0001){
  
  heatmap.2(heat_mat,col = cm, breaks=colbr, scale="none", notecol="black",
            dendrogram='none', density.info="none",Rowv=FALSE, Colv=FALSE,trace='none', key = T,
            xlab     = "G parameters",
            ylab     = expression( ~eta ~"parameters"), margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed",~ R==.(j))),par(cex.main=.5),
            srtCol=315, srtRow=315,adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,8),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,11),
            sepcolor="black")
  
  } else{
  dev.off()
  pdf(paste("plots/",j, "_fixed_R.pdf", sep=''))
  heatmap.2(heat_mat,col = cm, breaks=colbr, scale="none", notecol="black",
            dendrogram='none', density.info="none",Rowv=FALSE, Colv=FALSE,trace='none', key = F,
            xlab     = "G parameters",
            ylab     = expression( ~eta ~"parameters"), margins=c(5,5),main = bquote(paste(~delta[LL], " for fixed", ~ R==.(j))),par(cex.main=.5),
            srtCol=315, srtRow=315,adjCol = c(0,1),cexRow=0.8,cexCol=0.8,
            sepwidth=c(0.001,0.001),
            colsep=c(0,8),
            #sepwidth=c(0.05,0.05),
            rowsep=c(0,11),
            sepcolor="black")
  
  }
  dev.off()
  
  #jpeg(paste("plots/",j, ".fixed_R.dLLhist.jpg", sep=''))
  #hist(heat_mat,  main = bquote(paste("G vs. scale" ~ R==.(j))))
  #dev.off()
}


gompShape<-as.numeric(new_fit$gompShape)
gompRate<-as.numeric(new_fit$gompRate)

CV_gompShape<-sd(gompShape)/mean(gompShape)
CV_gompRate<-sd(gompRate)/mean(gompRate)

sd.Rls<-as.numeric(new_fit$sd.ls)
avg.Rls<-as.numeric(new_fit$AvgLS)

sd.Rls<-sd.Rls[!is.na(sd.Rls)]
avg.Rls<-avg.Rls[!is.na(avg.Rls)]

CV_sd.Rls<- sd(sd.Rls)/mean(sd.Rls)
CV_avg.Rls<-sd(avg.Rls)/mean(avg.Rls)

CV_names = c('CV of mean(LS)','CV of sd(LS)', 'CV of gompShape','CV of gompRate')
CV_results = c( CV_avg.Rls,CV_sd.Rls, CV_gompShape,CV_gompRate)
CV.df<-data.frame(CV_names,CV_results)
CV.df






