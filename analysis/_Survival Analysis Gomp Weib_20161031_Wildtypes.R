
library('varhandle')
library('flexsurv')
library('stringr')



##load original data for analyzing

new_rls<-read.csv("rls.csv")


my.data=new_rls

#my.data<-my.data[complete.cases(my.data$ref_lifespan_mean),]

# uniquely determine the ref_lifespans of each experiments individually
big_data<-my.data[!duplicated(my.data$ref_lifespans),] 
#row 94 ref_lifespan is blank


#big_data at column 26 has ref_lifespan_mean with NA values to be ignored
#big_data= big_data[!is.na(big_data[,26]), ]

#big_data = big_data[complete.cases(big_data$ref_lifespan_mean),]
ref_lifespan_mean<-list()



for (k in 1:(length(sort(unique(big_data$ref_lifespans))))){
#k=94
f1<-big_data$ref_lifespans
f<-unfactor(f1[k])
ref_lifespan_single<-as.numeric(unlist(str_split(f, ",")))
ref_lifespan_single_mean<-mean(ref_lifespan_single)
ref_lifespan_mean[[length(ref_lifespan_mean)+1]]<-ref_lifespan_single_mean

#avg.ref.ls<-round(mean(c(ref_lifespan_single)))
#x<-big_data$id
#ref_lifespan_mean<-round(big_data$ref_lifespan_mean[k])


}


big_data$ref_lifespan_mean<-unlist(ref_lifespan_mean)



fit_names = c( 'popsize','genotype','media','mat','sd.ls','medianLS','Delta_LL','AvgLS','gompLogLik','gompRate','gompShape','weibLogLik','weibScale','weibShape')
fit = data.frame(matrix(, nrow=length(fit_names), ncol=length(fit_names))) #set up a skeleton table
names(fit) = c( "popsize","genotype","media","mat","sd.ls","medianLS","Delta_LL","AvgLS","gompLoglik","gompRate","gompShape","weibLogLik","weibScale","weibShape")


lifespansChar_set<-big_data$set_lifespans
lifespansChar_ref<-big_data$ref_lifespans

rgompertz = function(N,G, R){
  ### I have a function to generate random-numbers from 
  ### (using the 'inversion method')
  
  return( (log(1-(G/R*log(1-runif(N)))))/G)
}


new_fit$CV_exp<-new_fit$sd.ls/new_fit$AvgLS

CV_exp_sim<-list()
sd_sim<-list()
mean_sim<-list()

R_squared<-list()

lifespansTemp_list<-list()

for (genotype_in in c(1:length(big_data$ref_name))){
  
  #genotype=117
  
  
  lifespansChar_set1<-unfactor(lifespansChar_set[genotype_in])
  #if lifespansChar_set1=="NA"
  lifespansTemp_set =as.numeric(unlist(str_split(lifespansChar_set1, ",")))

 
  lifespansChar_ref1<-unfactor(lifespansChar_ref[genotype_in])
  lifespansTemp_ref =as.numeric(unlist(str_split(lifespansChar_ref1, ",")))
  
  lifespansTemp<-c(lifespansTemp_set,lifespansTemp_ref)  
  
  
                       
  
  lifespansTemp[lifespansTemp < 0] <- 0
  lifespansTemp <-floor(lifespansTemp+0.5)
  lifespansTemp <- lifespansTemp[ lifespansTemp != 0 ]
  
  lifespansTemp_list[[length(lifespansTemp_list)+1]]<-lifespansTemp
  
  pop_size<-length(lifespansTemp)

  
  lifespanGomp = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
  lifespanWeib = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'weibull')  
  
## Fill in added columns in the controlConditions table with data from gompertz and weibull ##fitting. Columns are named according to respective variables
  media<-big_data$ref_media
  media=unfactor(media[genotype_in])
  mat<-big_data$ref_mating_type
  mat=unfactor(mat[genotype_in])
  avgLS = mean(lifespansTemp)
  StddevLS = sd(lifespansTemp)
  medianLS = median(lifespansTemp)
  gompShape = lifespanGomp$res[1,1]
  gompRate = lifespanGomp$res[2,1]
  gompLogLik = lifespanGomp$loglik
  gompAIC = lifespanGomp$AIC
  
  weibShape = lifespanWeib$res[1,1]
  weibScale = lifespanWeib$res[2,1]
  weibLogLik = lifespanWeib$loglik
  weibAIC = lifespanWeib$AIC   
  
  delta_LL = lifespanWeib$loglik- lifespanGomp$loglik
  
  genotype<-big_data$ref_name
  
  genotype=unfactor(genotype[genotype_in])
  
#fit_names = c( 'sd.ls','current.seed','i','medianLS','Delta_LL','AvgLS','gompLogLik','gompRat#e','gompShape','weibLogLik','weibScale','weibShape')
  
  fit = rbind(fit,c(pop_size,genotype,media,mat,StddevLS,medianLS,delta_LL,avgLS,gompLogLik,gompRate,gompShape,weibLogLik,weibScale,weibShape))
  #write.csv(fit, file="Results.csv", row.names=F)
  #names(fit) = c( "genotype","media","mat","sd.ls","medianLS","Delta_LL","AvgLS","gompLoglik","gompRate","gompShape","weibLogLik","weibScale","weibShape")
  
  
#fit= fit[!is.na(fit[,1]), ]
  
}


new_fit<-fit[(length(fit_names)+1):length(fit[,1]),]

write.csv(new_fit,"conditionsWeibRedo_emine_10312016.csv")


###read the analyzed data

new_fit<-read.csv("conditionsWeibRedo_emine_10312016.csv")

new_fit<-new_fit[-c(1),-c(21)]

r.square = function( lifespan.observed, lifespan.fitted){

  # compute coefficient of determination
  if (class(lifespan.fitted) == "numeric") {
    return(cor(lifespan.observed, lifespan.fitted)^2)
  } else {
    R2 = double(ncol(lifespan.observed))
    for (ic in 1:ncol(lifespan.observed))
      R2[ic] = cor(lifespan.observed[,ic], lifespan.fitted[,ic])^2
    return(R2)
  }
}



new_fit$CV_exp<-new_fit$sd.ls/new_fit$AvgLS

CV_exp_fit<-list()
sd_fit<-list()
mean_fit<-list()

R_squared<-list()

#lifespansTemp_list


for (i in 1:length(new_fit[,1])){
  
  lifespan_fit<- exp(new_fit$gompRate[i]/new_fit$gompShape[i] *(1-exp(new_fit$gompShape[i] * lifespansTemp_list[[i]])))
  
  #lifespan_fit<-rgompertz(new_fit$popsize[i],new_fit$gompShape[i],new_fit$gompRate[i])
  #lifespan_fit<-floor(lifespan_fit+0.5)
  
  #avg.lifespansTemp<-round(mean(lifespansTemp_list[[i]]))
  lifespan.observed<-lifespansTemp_list[[i]]
  
  
  #cook's distance test
  #ols<-lm(lifespan.observed~lifespan_fit)
  #library(influence.ME)
  #infl <- influence(ols, obs = TRUE)
  #opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
  #plot(ols, las = 1)
  #d1 <- cooks.distance(ols)
  #d1<-d1[d1<3*mean(d1)]
  #plot(d1,ylab="Cooks distance")
  

  
  R_squ<- r.square(lifespan.observed,lifespan_fit)
  R_squared[[length(R_squared)+1]]<-R_squ
  
  
  
  sd_fit[[length(sd_fit)+1]]<-sd(lifespan_fit)
  mean_fit[[length(mean_fit)+1]]<-mean(lifespan_fit)
  CV_fit<-sd(lifespan_fit)/mean(lifespan_fit)
  CV_exp_fit[[length(CV_exp_fit)+1]]<-CV_fit
  
}


R_squared<-unlist(R_squared)
R_squared<-R_squared[!is.na(R_squared)]

summary(R_squared)

pdf(paste("plots/", "Histogram_empirical_vs_fitted_R_squared.pdf", sep=''))
hist(unlist(R_squared),breaks=length(R_squared)/5,xlim=c(0.75,1),ylim=c(0,150))


sd_fit<-unlist(sd_fit)
mean_fit<-unlist(mean_fit)

pdf(paste("plots/", "Histogram_fitted_data_summarize.pdf", sep=''))
par(mfrow=c(1,2)) 
hist(sd_fit,breaks=length(sd_fit)/5)
hist(mean_fit,breaks=length(mean_fit)/5)

#calculate CV's
CV_exp_fit<-unlist(CV_exp_fit)
pdf(paste("plots/", "Histogram_empirical_vs_fitted_data_CV.pdf", sep=''))
par(mfrow=c(1,2))
hist(CV_exp_fit,main="CV of fitted data",breaks=length(CV_exp_fit)/5,xlim=c(0,1),ylim=c(0,120),col="green")
hist(new_fit$CV_exp,main="CV of fitted data",breaks=length(CV_exp_fit)/5,xlim=c(0.2,0.6),ylim=c(0,120),col="blue")



CV_ratio<-new_fit$CV_exp/CV_exp_fit
CV_ratio<-CV_ratio[!is.na(CV_ratio)]
summary(CV_ratio)
hist(CV_ratio,breaks=length(CV_ratio)/5,xlim=c(0,2))






pdf(paste("plots/", "Histogram_empirical_data_model_params.pdf", sep=''))

par(mfrow=c(3,1)) 
hist(new_fit$AvgLS,xlab="Mean Lifespan",main="Empirical Data",
     breaks=length(new_fit$AvgLS)/10,col="gray",xlim=c(10,40),ylim=c(0,70))
mtext('Total experiments=5303',line=-2,cex=0.5)
box(lty='1373',col='blue')



hist(new_fit$gompShape,xlab="G shape parameters",main="",
     breaks=length(new_fit$gompShape)/10,col="blue",xlim=c(0,0.25),ylim=c(0,200))

box(lty='1373',col='blue')

hist(log(new_fit$gompRate),xlab="log(R) rate parameters",main="",
     breaks=length(new_fit$gompRate)/10,col="green",xlim=c(-8,-3),ylim=c(0,150))

box(lty='1373',col='blue')
dev.off()

Delta_LL<-unlist(new_fit$Delta_LL)

pdf(paste("plots/", "Histogram_Delta_LL_empirical_data_all.pdf", sep=''))

fifty.percent<-quantile(Delta_LL, c(.50)) 

hist(Delta_LL,breaks=length(Delta_LL)/10,col="blue",xlim=c(-50,120),ylim=c(0,1200),xlab=expression(paste( ~delta[LL])),main="Histogram of Delta_LL for Experimental lifespan data")
mtext('Total experiments=5303',line=-3,at=70,side=3,cex=0.8)
mtext('median=5.22',line=-4,at=70,side=3,cex=0.8)
mtext('mean=10.76',line=-5,at=70,side=3,cex=0.8)
box(lty='1373',col='blue')



abline(v = fifty.percent, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-20,at=1000,cex=0.8)



#hist(new_fit$avgLS)
dev.off()


WT.BY4742<- new_fit[new_fit$genotype=="BY4742",]

#WT.BY4742.temp<-WT.BY4742[WT.BY4742$temp==30,]
WT.BY4742.media<-WT.BY4742[WT.BY4742$media=="YPD",]

WT.BY4742.media<-WT.BY4742.media[WT.BY4742.media$mat=="MATalpha",]

WT.BY4742.media= WT.BY4742.media[!is.na(WT.BY4742.media[,1]), ]




pdf(paste("plots/", "Histogram_empirical_data_WT_BY4742_model_params.pdf", sep=''))

par(mfrow=c(3,1)) 
hist(as.numeric(WT.BY4742.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4742",
     breaks=length(as.numeric(WT.BY4742.media$AvgLS))/5,col="gray",xlim=c(10,40),ylim=c(0,40))
mtext('Total_experiments=2108',line=-3,cex=0.5)
box(lty='1373',col='blue')


hist(as.numeric(WT.BY4742.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$gompShape))/5,col="blue",xlim=c(0,0.16),ylim=c(0,70))

box(lty='1373',col='blue')

hist(log10(as.numeric(WT.BY4742.media$gompRate)),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(WT.BY4742.media$gompRate))/10,col="green",xlim=c(-3,-1),ylim=c(0,150))

box(lty='1373',col='blue')
dev.off()

pdf(paste("plots/", "Histogram_Delta_LL_empirical_data_WT_BY4742.pdf", sep=''))

fifty.percent2<-quantile(WT.BY4742.media$Delta_LL, c(.50)) 


hist(WT.BY4742.media$Delta_LL,breaks=length(WT.BY4742.media$Delta_LL)/5,col="blue",xlim=c(-50,120),ylim=c(0,600),xlab=expression(paste( ~delta[LL])),main="Histogram of Delta_LL for WT BY4742")


mtext('Total experiments=2108',line=-3,at=70,side=3,cex=0.8)
mtext('median=4.342',line=-4,at=70,side=3,cex=0.8)
mtext('mean=10.1',line=-5,at=70,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent2, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-20.5,at=500,cex=0.8)


dev.off()






WT.BY4741<- new_fit[new_fit$genotype=="BY4741",]



#WT.BY4742.temp<-WT.BY4742[WT.BY4742$temp==30,]
WT.BY4741.media<-WT.BY4741[WT.BY4741$media=="YPD",]

WT.BY4741.media<-WT.BY4741.media[WT.BY4741.media$mat=="MATa",]

WT.BY4741.media= WT.BY4741.media[!is.na(WT.BY4741.media[,1]), ]






pdf(paste("plots/", "Histogram_empirical_data_WT_BY4741_model_params.pdf", sep=''))
par(mfrow=c(3,1)) 
hist(as.numeric(WT.BY4741.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4741",
     breaks=length(as.numeric(WT.BY4741.media$AvgLS))/5,col="gray",xlim=c(10,40),ylim=c(0,60))

mtext('Total experiments=381',line=-2,cex=0.5)
box(lty='1373',col='blue')

hist(as.numeric(WT.BY4741.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4741.media$gompShape))/5,col="blue",xlim=c(0.02,0.16),ylim=c(0,40))

box(lty='1373',col='blue')

hist(log10(as.numeric(WT.BY4741.media$gompRate)),xlab="log(R) rate parameters",main="",
     breaks=length(as.numeric(WT.BY4741.media$gompRate))/5,col="green",xlim=c(-3,-1),ylim=c(0,40))
box(lty='1373',col='blue')
dev.off()


pdf(paste("plots/", "Histogram_Delta_LL_empirical_data_WT_BY4741.pdf", sep=''))

fifty.percent1<-quantile(WT.BY4741.media$Delta_LL, c(.50)) 


hist(WT.BY4741.media$Delta_LL,breaks=length(WT.BY4741.media$Delta_LL)/5,col="blue",xlim=c(-10,60),ylim=c(0,80),xlab=expression(paste( ~delta[LL])),main="Histogram of Delta_LL for WT BY4741 ")


mtext('Total experiments=381',line=-3,at=30,side=3,cex=0.8)
mtext('median=3.6',line=-4,at=30,side=3,cex=0.8)
mtext('mean=6.025',line=-5,at=30,side=3,cex=0.8)
box(lty='1373',col='blue')


abline(v = fifty.percent1, lty = 1,col="red")
mtext('50% quantile',col="red",side=4,line=-24,at=60,cex=0.8)

dev.off()


WT.BY4743<- new_fit[new_fit$genotype=="BY4743",]

WT.BY4743.media<-WT.BY4743[WT.BY4743$media=="YPD",]

WT.BY4743.media<-WT.BY4743.media[WT.BY4743.media$mat=="diploid",]

WT.BY4743.media= WT.BY4743.media[!is.na(WT.BY4743.media[,1]), ]



#pdf(paste("plots/", "Histogram_empirical_data_WT_BY4743_model_params.pdf", sep=''))
par(mfrow=c(3,1)) 
hist(as.numeric(WT.BY4743.media$AvgLS),xlab="Mean Lifespan",main="Wild type: BY4743",
     breaks=length(as.numeric(WT.BY4743.media$AvgLS)),col="gray",xlim=c(25,45),ylim=c(0,30))


hist(as.numeric(WT.BY4743.media$gompShape),xlab="G shape parameters",main="",
     breaks=length(as.numeric(WT.BY4743.media$gompShape)),col="blue",xlim=c(0.05,0.25),ylim=c(0,30))

hist(as.numeric(WT.BY4743.media$gompRate),xlab="R rate parameters",main="",
     breaks=length(as.numeric(WT.BY4743.media$gompRate)),col="green",xlim=c(0,0.01),ylim=c(0,40))
#dev.off()



###check coefficent of the variation: CV
###Because small CV indicate less noisy and more robustness

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


gompShape_BY4742<-as.numeric(WT.BY4742.media$gompShape)
gompRate_BY4742<-as.numeric(WT.BY4742.media$gompRate)

CV_gompShape_BY4742<-sd(gompShape_BY4742)/mean(gompShape_BY4742)
CV_gompRate_BY4742<-sd(gompRate_BY4742)/mean(gompRate_BY4742)

sd.Rls.BY4742<-as.numeric(WT.BY4742.media$sd.ls)
avg.Rls.BY4742<-as.numeric(WT.BY4742.media$AvgLS)

sd.Rls.BY4742<-sd.Rls.BY4742[!is.na(sd.Rls.BY4742)]
avg.Rls.BY4742<-avg.Rls.BY4742[!is.na(avg.Rls.BY4742)]

CV_sd.Rls.BY4742<- sd(sd.Rls.BY4742)/mean(sd.Rls.BY4742)
CV_avg.Rls.BY4742<-sd(avg.Rls.BY4742)/mean(avg.Rls.BY4742)


CV_names = c('CV of mean(LS)','CV of sd(LS)', 'CV of gompShape','CV of gompRate')
CV_results = c( CV_avg.Rls.BY4742,CV_sd.Rls.BY4742, CV_gompShape_BY4742,CV_gompRate_BY4742)
CV.df.BY4742<-data.frame(CV_names,CV_results)
CV.df.BY4742



gompShape_BY4741<-as.numeric(WT.BY4741.media$gompShape)
gompRate_BY4741<-as.numeric(WT.BY4741.media$gompRate)

CV_gompShape_BY4741<-sd(gompShape_BY4741)/mean(gompShape_BY4741)
CV_gompRate_BY4741<-sd(gompRate_BY4741)/mean(gompRate_BY4741)

sd.Rls.BY4741<-as.numeric(WT.BY4741.media$sd.ls)
avg.Rls.BY4741<-as.numeric(WT.BY4741.media$AvgLS)

sd.Rls.BY4741<-sd.Rls.BY4741[!is.na(sd.Rls.BY4741)]
avg.Rls.BY4741<-avg.Rls.BY4741[!is.na(avg.Rls.BY4741)]

CV_sd.Rls.BY4741<- sd(sd.Rls.BY4741)/mean(sd.Rls.BY4741)
CV_avg.Rls.BY4741<-sd(avg.Rls.BY4741)/mean(avg.Rls.BY4741)


CV_names = c('CV of mean(LS)','CV of sd(LS)', 'CV of gompShape','CV of gompRate')
CV_results = c( CV_avg.Rls.BY4741,CV_sd.Rls.BY4741, CV_gompShape_BY4741,CV_gompRate_BY4741)
CV.df.BY4741<-data.frame(CV_names,CV_results)
CV.df.BY4741





#10/13/2016
#next thing to do compare all CV values both from simulation and seperately in WILD TYPE data


#How do noises influence likelihood surface of Gompertz and Weibull model? 

WT.BY4742.media$CV<-WT.BY4742.media$sd.ls/WT.BY4742.media$AvgLS

lm( WT.BY4742.media$gompLoglik ~ WT.BY4742.media$CV )

summary(lm( WT.BY4742.media$gompLoglik ~ WT.BY4742.media$CV ))


# We expect that larger CV (more noisy) 
#data will bring down LLH of both Gomeprtz and Weibull models. 
#However, Weibull’s decreasing LLH maybe slower than Gompertz. 


plot(WT.BY4742.media$gompLoglik ~ WT.BY4742.media$CV,col='green',ylab='' )
par(new=TRUE)
plot(WT.BY4742.media$weibLogLik ~ WT.BY4742.media$CV,col='red',ylab='Loglik values')

plot(WT.BY4742.media$gompLoglik ~ WT.BY4742.media$CV )
m<-lm( WT.BY4742.media$gompLoglik ~ WT.BY4742.media$CV )
abline(m,col='red',lty=2,lwd=2)


summary(lm( WT.BY4742.media$weibLogLik ~ WT.BY4742.media$CV))

plot(WT.BY4742.media$weibLogLik ~ WT.BY4742.media$CV )
m<-lm( WT.BY4742.media$weibLogLik ~ WT.BY4742.media$CV )
abline(m,col='red',lty=2,lwd=2)




# We expect that larger CV (more noisy) 
#data will bring down LLH of both Gomeprtz and Weibull models. 
#However, Weibull’s decreasing LLH maybe slower than Gompertz. 










