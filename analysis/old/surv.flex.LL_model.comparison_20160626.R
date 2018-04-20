#title: "GaussianNoises and Delta_LL of model comparison project"
#author: "e guven"
#date: "May 3, 2016"
#output: html_document


#EG started this on 20160403 
#how does the noise effect the loglikelihood 
#Weibull model->theta is a scale, gamma is a shape parameter >0
#Gompertz Model-> beta(G) is a scale, alpha(R) is a shape parameter>0 

#difference function of loglikelihood function of gompertz and weibull p.d.fs
#test if L(Weibull,X)>L(Gompertz,X) for parameters Weibull model->theta is a scale, gamma is a shape parameter >0
#Gompertz Model-> beta is a scale, alpha is a shape parameter>0 

#For additive Gaussian noise e ~ N (0, sigma^2) with known variance sigma^2
#sd of gaussian noise function
#max sd would be = 3*mean(inverse.gomp.CDF)
#min sd would be mean(inverse.gomp.CDF)



#fix the paramaters of weibull model

require(flexsurv)

gamma=0.01
theta=0.25
#test G and R in nested for loops
R=0.01
G=0.2
#R should be in [0, 0.05], and G should be [0.05, 0.3].  
N=2000 # population size


#rgompertz(alpha,beta,N) gives random Gompertz numbers from inverse CDF of Gompertz
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
gompertz.random<-rgompertz(N,0.05,0.05)
average.lifespan=mean(gompertz.random)

scale=0

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

summary(lm(log10(GC$mortality.rate[1:80]) ~ GC$t[1:80]))
#then plot log(mortality rate) ~ time.
#For Weibull model, log(moretality rate) ~ log(time) give the linear form.

pdf(paste("plots/","Gompertz.semi.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ GC$t, type='l') #linear for Gompertz, semi-log plot
text(10,0,"R2= 0.9677")
dev.off()

summary(lm(log10(GC$mortality.rate[1:80]) ~ log10(GC$t[1:80])))
pdf(paste("plots/","Weibull.log.log.plot.batch.pdf", sep=''), width=5, height=5)
plot( log10(GC$mortality.rate) ~ log10(GC$t), type='l'  ) #linear for Weibull, log-log plot
text(0.75,0,"R2=0.89")
dev.off()










#put the variable values from loops into lists
sderr<- list()

G<-list()
R<-list()

sd.rgompertz<-list()
scale<-list()
MeanLF<-list()
sdLS<-list()
Delta_LL.flex<-list()
LWei.flex<-list()
LGomp.flex<-list()
G.flex.estimated<-list()
R.flex.estimated<-list()

gamma.flex.estimated<-list()
theta.flex.estimated<-list()

N=2000

  for (G_in in c(0.05,0.08, 0.1,0.12,0.15,0.17, 0.2, 0.25)){
    for (R_in in c(1E-3, 0.002,0.003, 0.005,0.008, 0.01,0.03, 0.05)){ #fix alpha or in other words R shape parameter
    #for (sd in seq(round.lifespan,3*round.lifespan,by=1)){
    #for (i in c(0,0.09,0.05,0.1,0.2,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7)){ 
      for (i in c(0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3,0.5,1,1.5,2,3,4,5)){ 
      
      #generate gompertz random numbers (lifespan) 
      #prediction
        require(flexsurv)
      set.seed(123)
     # gompertz.random<-rgompertz(N,R_in,G_in)
     
     
      
      gompertz.random<-rgompertz(N, G_in, R_in)
     
     
     #lifespan = gompertz.random + rnorm(N, mean=0, sd=sd(gompertz.random)*i
             
     summary(lifespan)
             
      
      scale[[length(scale)+1]]=i
      
      
      sdrgompertz<-sd(gompertz.random)
      
      sd.rgompertz[[length(sd.rgompertz)+1]]=sdrgompertz
      
      #generate gaussion random numbers 
      #gaussian<-rnorm(N, mean = 0, 1)*i
      
      #lifespan =  gompertz.random+gaussian
      
      average.lifespan=mean(gompertz.random)
      #store average.lifespan into MeanLF list
      MeanLF[[length(MeanLF)+1]]=average.lifespan
      
      #Log likelihood function for the Weibull model
    
      
      #calculate LL and noise change
      sderr[[length(sderr)+1]] = sd(gompertz.random)*i   
      
      G[[length(G)+1]]=G_in
      #switch to alpha.seq when for fixed beta
      R[[length(R)+1]]=R_in
      
      #todo use flexsurv to calculate the LL
      
      #flexsurv only works with positive variables.
      #fix gaussian std to 0
      
      
      X.flex= gompertz.random + rnorm(N, mean=0, sd=sd(gompertz.random)*i)
      
      X.flex[X.flex < 0] <- 0
      X.flex<-floor(X.flex+0.5)
     
      X.flex <- X.flex[ X.flex != 0 ]
      
      
      fitGomp = flexsurvreg(formula = Surv(X.flex) ~ 1, dist="gompertz")
      fitWei = flexsurvreg(formula = Surv(X.flex) ~ 1, dist="weibull")
      
      
      LWei.flex[[length(LWei.flex)+1]]=fitWei$loglik
      
      LGomp.flex[[length(LGomp.flex)+1]]=fitGomp$loglik
      
      param.Gomp<-fitGomp$res; R.flex<-param.Gomp[2]; G.flex<-param.Gomp[1];
      
      R.flex.estimated[[length(R.flex.estimated)+1]]<-R.flex
      G.flex.estimated[[length(G.flex.estimated)+1]]<-G.flex
      
      param.Wei<-fitWei$res; gamma.flex<-param.Wei[1]; theta.flex<-param.Wei[2];
      
      gamma.flex.estimated[[length(gamma.flex.estimated)+1]]<-gamma.flex; 
      theta.flex.estimated[[length(theta.flex.estimated)+1]]<-theta.flex
      
      delta_flexsurv=fitWei$loglik-fitGomp$loglik 
      
      #fitWei$loglik
      
      Delta_LL.flex[[length(Delta_LL.flex)+1]]=delta_flexsurv
      
    }
  }
}




results = data.frame(cbind(sderr), cbind(R),cbind(R.flex.estimated),cbind(G),cbind(G.flex.estimated) , cbind(Delta_LL.flex),
                    cbind(LWei.flex),cbind(LGomp.flex),cbind(MeanLF),cbind(sd.rgompertz),cbind(scale))

results_mat<-as.matrix(results)

write.csv(results_mat,file="Results.csv")

#write.csv(results_mat,file="noise_zero.csv")


dLL.flex<-unlist(results$Delta_LL.flex)

LLGomp.flex<- unlist(results$LGomp.flex)

LLWei.flex<- unlist(results$LWei.flex)
simulated.G<-unlist(results$G)
estimated.G.flex<-unlist(results$G.flex.estimated)
simulated.R<-unlist(results$R)
estimated.R.flex<-unlist(results$R.flex.estimated)
scale<-unlist(scale)





summary(lm(simulated.G~estimated.G.flex))
summary(lm(simulated.R~estimated.R.flex))


pdf(paste("plots/","simulated.G.vs.estimated.G.flex.batch.pdf", sep=''), width=5, height=5)
plot(simulated.G~estimated.G.flex)
text(0.06,0.20,"R^2=0.71")
dev.off()

pdf(paste("plots/","simulated.R.vs.estimated.R.flex.batch.pdf", sep=''), width=5, height=5)
plot(simulated.R~estimated.R.flex)
text(0.01,0.04,"R^2=0.94")
dev.off()


results.sub<-data.frame(cbind(sderr),cbind(R),cbind(G),cbind(Delta_LL.flex),cbind(scale))

results_sub_mat<-as.matrix(results.sub)

write.csv(results_sub_mat,file="Results.sub.csv")



R.els = unlist( unique(results.sub$R))
colnum = length(R.els)

tmp = unlist( unique(results.sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = scale.els
colnames(mat) = R.els
for (k in  c(0.05,0.08, 0.1,0.12,0.15,0.17, 0.2, 0.25)){
  
  #R in c(1E-3, 0.002, 0.005,0.008, 0.01,0.03, 0.05)
  data = results.sub[results.sub[,3]==k, 4]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  rownames(heat_mat) <- scale.els
  colnames(heat_mat) <- R.els
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  jpeg(paste("plots.flex/",k, ".fixed.G.jpg", sep=''))
  
  #paste(“myplot_”, i, “.jpeg”, sep=””)
  
  heatmap.2(heat_mat, cellnote=hM,col = cm.colors(256), scale="none", notecol="black",  margins=c(5,10),
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
            xlab     = "R parameters",
            ylab     = "scale", main = bquote(paste("R~scale for fixed" ~ G==.(k))),par(cex.main=.5), adjCol = c(1,1),cexRow=0.8,cexCol=0.8)
  
  
  dev.off()
  
}

G.els = unlist( unique(results.sub$G))
colnum = length(G.els)

R.els=unlist(unique(results.sub$R))
rownum = length(R.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = R.els
colnames(mat) = G.els
for (n in c(0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3,0.5,1,1.5,2,3,4,5) ){
  #for (n in c(0)){
  data = results.sub[results.sub[,5]==n, 4]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  rownames(heat_mat) <- R.els
  
  colnames(heat_mat) <- G.els
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  jpeg(paste("plots.flex/",n, ".fixed_scale.jpg", sep=''))
  
  #paste(“myplot_”, i, “.jpeg”, sep=””)
  
  heatmap.2(heat_mat, cellnote=hM,col = cm.colors(256), scale="none", notecol="black",  margins=c(5,10),
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
            xlab     = "G parameters",
            ylab     = "R parameters", main = bquote(paste("G~R for fixed"~ scale==.(n))),par(cex.main=.5),adjCol=c(1,1),cexRow=0.8,cexCol=0.8)
  
  
  dev.off()
  
}

G.els = unlist( unique(results.sub$G))
colnum = length(G.els)

tmp = unlist( unique(results.sub$scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
colnames(mat) = G.els
rownames(mat) = scale.els

#R should be in [0, 0.05], and G should be [0.05, 0.3].  
for (j in c(1E-3, 0.002, 0.003, 0.005,0.008, 0.01,0.03, 0.05) ){
  data = results.sub[results.sub[,2]==j, 4]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  colnames(heat_mat) <- G.els
  
  rownames(heat_mat) <- scale.els
  library(gplots)
  hM <- format(round(heat_mat, 1))
  #data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  jpeg(paste("plots.flex/",j, ".fixed_R.jpg", sep=''))
  
  #paste(“myplot_”, i, “.jpeg”, sep=””)
  
  heatmap.2(heat_mat, cellnote=hM,col = cm.colors(256), scale="none", notecol="black",  margins=c(5,10),
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
            xlab     = "G parameters",
            ylab     = "scale", main = bquote(paste("G~scale for fixed" ~ R==.(j))),par(cex.main=.5),adjCol=c(1,1),cexRow=0.8,cexCol=0.8)
  
  
  dev.off()
  
}



