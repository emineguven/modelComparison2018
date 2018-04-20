

new_rls<-read.csv("rls 2016-08-02.csv")
my.data=new_rls

#my.data<-my.data[complete.cases(my.data$ref_lifespan_mean),]

# uniquely determine the ref_lifespans of each experiments individually
big_data<-my.data[!duplicated(my.data$ref_lifespans),] 
#row 94 ref_lifespan is blank

big_data<-big_data[-c(94,56,85,94,99,109),]

#big_data at column 26 has ref_lifespan_mean with NA values to be ignored
#big_data= big_data[!is.na(big_data[,26]), ]

#big_data = big_data[complete.cases(big_data$ref_lifespan_mean),]
ref_lifespan_mean<-list()
for (k in 1:(length(sort(unique(big_data$ref_lifespans))))){
  #k=94
  f1<-big_data$ref_lifespans
  f<-unfactor(f1[k])
  ref_lifespan_single<-as.numeric(unlist(strsplit(f, ",")))
  ref_lifespan_single_mean<-mean(ref_lifespan_single)
  ref_lifespan_mean[[length(ref_lifespan_mean)+1]]<-ref_lifespan_single_mean
  
  #avg.ref.ls<-round(mean(c(ref_lifespan_single)))
  #x<-big_data$id
  #ref_lifespan_mean<-round(big_data$ref_lifespan_mean[k])
  
  
}
big_data$ref_lifespan_mean<-unlist(ref_lifespan_mean)



fit_names = c( 'genotype','media','mat','sd.ls','medianLS','Delta_LL','AvgLS','gompLogLik','gompRate','gompShape','weibLogLik','weibScale','weibShape')
fit = data.frame(matrix(, nrow=length(fit_names), ncol=length(fit_names))) #set up a skeleton table
names(fit) = c( "genotype","media","mat","sd.ls","medianLS","Delta_LL","AvgLS","gompLoglik","gompRate","gompShape","weibLogLik","weibScale","weibShape")



for (genotype_in in c(1:length(big_data$set_name))){
  
  #genotype=117
  
  lifespansChar_set<-big_data$set_lifespans
  lifespansChar_set1<-unfactor(lifespansChar_set[genotype_in])
  lifespansTemp_set =as.numeric(unlist(strsplit(lifespansChar_set1, ",")))
  
  
  
  lifespansChar_ref<-big_data$ref_lifespans
  lifespansChar_ref1<-unfactor(lifespansChar_ref[genotype_in])
  lifespansTemp_ref =as.numeric(unlist(strsplit(lifespansChar_ref1, ",")))
  
  lifespansTemp<-c(lifespansTemp_set,lifespansTemp_ref)             
  
  
  lifespansTemp[lifespansTemp < 0] <- 0
  lifespansTemp <-floor(lifespansTemp+0.5)
  lifespansTemp <- lifespansTemp[ lifespansTemp != 0 ]
  
  lifespanGomp = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
  lifespanWeib = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'weibull')  
  
  ### Fill in added columns in the controlConditions table with data from gompertz and weibull fitting. Columns are named according to respective variables
  media<-big_data$set_media
  media=unfactor(media[genotype_in])
  mat<-big_data$set_mating_type
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
  
  genotype<-big_data$set_genotype
  
  genotype=unfactor(genotype[genotype_in])
  
  #fit_names = c( 'sd.ls','current.seed','i','medianLS','Delta_LL','AvgLS','gompLogLik','gompRate','gompShape','weibLogLik','weibScale','weibShape')
  
  fit = rbind(fit,c(genotype,media,mat,StddevLS,medianLS,delta_LL,avgLS,gompLogLik,gompRate,gompShape,weibLogLik,weibScale,weibShape))
  #write.csv(fit, file="Results.csv", row.names=F)
  #names(fit) = c( "genotype","media","mat","sd.ls","medianLS","Delta_LL","AvgLS","gompLoglik","gompRate","gompShape","weibLogLik","weibScale","weibShape")
  
  
  #fit= fit[!is.na(fit[,1]), ]
  
}

new_fit<-fit[length(fit_names)+1:length(fit),]

#hist(new_fit$avgLS)








