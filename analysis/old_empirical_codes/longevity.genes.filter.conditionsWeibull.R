
library('survival')
library('KMsurv')
library('data.table')
library('RSQLite')
library('flexsurv')
library('muhaz')

#########################################################################################################
#########################################################################################################

#####################Declare names of files and folders. Move working directory to desired folder. Establish connection to database
#folderLinearHazard = '//udrive.uw.edu/udrive/Kaeberlein Lab/RLS Survival Analysis 2/Linear Hazard'
Folder = "~/Desktop/emprical data/9.ken"

setwd(Folder)

setwd(Folder)
tb<-read.csv("conditionsWeibRedo_qin.csv") ##### read conditions table from a CSV file

tb_filter<-tb[tb$media == "YPD",] 
tb_filter<-tb_filter[tb_filter$n > 50,]

tb_filter_matalpha <- tb_filter[tb_filter$mat == "MATalpha",]

MAT_alpha<-length(tb_filter_matalpha[,1])

#tb_filter_matalpha[grep("[a-z][a-z][a-z][1-9]", tb_filter_matalpha$genotype, perl = TRUE, ignore.case=T),]

tb_filter_3L1D<-tb_filter_matalpha[grep("^[a-z]{3}[1-9]$", tb_filter_matalpha$genotype, perl = TRUE),]

MAT_alpha_genotype<-length(tb_filter_3L1D[,1])

tb_filter_WT<-tb_filter_matalpha[grep("[A-Z]{2}[1-9][1-9][1-9][1-9]$",tb_filter_matalpha$genotype,perl=T,ignore.case=F),]

MAT_alpha_WT<-length(tb_filter_WT[,1])

filtered.matalpha<-rbind(tb_filter_WT,tb_filter_3L1D)


pdf(paste("plots/", "G_MATalpha_WT_quantiles.pdf", sep=''))
hist(filtered.matalpha$gompShape,breaks=length(filtered.matalpha),xlim=c(0,0.3),ylim=c(0,250),xlab="G parameters",main="MATalpha gompertz shape parameter hist")
arrows(filtered.matalpha[1:4,]$gompShape,135,filtered.matalpha[1:4,]$gompShape,2,col="red" )
text(filtered.matalpha[1,]$gompShape,145,"BY4742",col="red",pos = 3, cex = 0.8,srt=90)
#text(filtered.matalpha[3,]$gompShape,152,"BY4742",col="red",pos = 4, cex = 0.6,srt=90)
text(filtered.matalpha[2,]$gompShape,145,"BY4742",col="red",pos = 3, cex = 0.8, srt = 90)
text(filtered.matalpha[4,]$gompShape,145,"BY4742",col="red",pos = 3, cex = 0.8, srt = 90)
#text(filtered.matalpha[4,]$gompShape,55,"BY4742",pos = 4, cex = 0.6, srt = 90)
#lines(density(filtered.mata.data$gompShape),na.rm=TRUE)
quantiles.matalpha<-quantile(filtered.matalpha$gompShape, c(.05,.10,.15,.85,.90,.95))
#lines(A)
quantiles.matalpha<-unname(quantiles.matalpha)
percent.5<-quantiles.matalpha[[1]]
percent.10<-quantiles.matalpha[[2]]
percent.15<-quantiles.matalpha[[3]]
percent.85<-quantiles.matalpha[[4]]
percent.90<-quantiles.matalpha[[5]]
percent.95<-quantiles.matalpha[[6]]

segments(percent.5,150,percent.5,0, col= "magenta",lwd=2,lty=2)
text(percent.5,150,"5%",pos = 3, cex = 0.8, srt = 90)

segments(percent.10,160,percent.10,0, col= "magenta",lwd=2,lty=2)
text(percent.10,160,"10%",pos = 3, cex = 0.8, srt = 90)

segments(percent.15,190,percent.15,0, col= "magenta",lwd=2,lty=2)
text(percent.15,190,"15%",pos = 4, cex = 0.8, srt = 90)

segments(percent.85,190,percent.85,0, col= "royalblue",lwd=2,lty=2)
text(percent.85,190,"85%",pos = 3, cex = 0.8, srt = 90)

segments(percent.90,160,percent.90,0, col= "royalblue",lwd=2,lty=2)
text(percent.90,160,"90%",pos = 3, cex = 0.8, srt = 90)

segments(percent.95,150,percent.95,0, col= "royalblue",lwd=2,lty=2)
text(percent.95,150,"95%",pos = 3, cex = 0.8, srt = 90)

dev.off()

#conditions_names = c('MAT_a','MAT_a_WT','MAT_a_percent5','MAT_a_percent10',
#                     'MAT_a_percent15','MAT_a_percent85','MAT_a_percent90','MAT_a_percent95',
 #                    'MAT_alpha','MAT_alpha_WT','MAT_alpha_percent5','MAT_alpha_percent10',
 #                    'MAT_alpha_percent15','MAT_alpha_percent85','MAT_alpha_percent90','MAT_alpha_percent95')
 

  
  MAT_alpha_percent5<- filtered.matalpha[filtered.matalpha$gompShape<percent.5,]
  MAT_alpha_percent5<-length(MAT_alpha_percent5[,1])
  MAT_alpha_percent10<- filtered.matalpha[filtered.matalpha$gompShape<percent.10,]
  MAT_alpha_percent10<-length(MAT_alpha_percent10[,1])
  MAT_alpha_percent15<- filtered.matalpha[filtered.matalpha$gompShape<percent.15,]
  MAT_alpha_percent15<-length(MAT_alpha_percent15[,1])
  MAT_alpha_percent85<- filtered.matalpha[filtered.matalpha$gompShape<percent.85,]
  MAT_alpha_percent85<-length(MAT_alpha_percent85[,1])
  MAT_alpha_percent90<- filtered.matalpha[filtered.matalpha$gompShape<percent.90,]
  MAT_alpha_percent90<-length(MAT_alpha_percent90[,1])
  MAT_alpha_percent95<- filtered.matalpha[filtered.matalpha$gompShape<percent.95,]
  MAT_alpha_percent95<-length(MAT_alpha_percent95[,1])




tb_filter_mata <- tb_filter[tb_filter$mat == "MATa",]

MAT_a<-length(tb_filter_mata[,1])

#tb_filter_mata[grep("[a-z][a-z][a-z][1-9]", tb_filter_mata$genotype, perl = TRUE, ignore.case=T),]

tb_filter_3L1D<-tb_filter_mata[grep("^[a-z]{3}[1-9]$", tb_filter_mata$genotype, perl = TRUE),]

MAT_a_genotype<-length(tb_filter_3L1D[,1])

tb_filter_WT<-tb_filter_mata[grep("^BY474",tb_filter_mata$genotype,perl=T,ignore.case=F),]

MAT_a_WT<-length(tb_filter_WT[,1])

filtered.mata<-rbind(tb_filter_WT,tb_filter_3L1D)

pdf(paste("plots/", "G_MATa_WT_quantiles.pdf", sep=''))

hist(filtered.mata$gompShape,breaks=length(filtered.mata),xlim=c(0,0.3),ylim=c(0,70),xlab="G parameters",main="MATa gompertz shape parameter hist")
arrows(filtered.mata[1:2,]$gompShape,40,filtered.mata[1:2,]$gompShape,2,col="red" )
text(filtered.mata[1,]$gompShape,48,"BY4741", col="red", pos = 2, cex = 0.8, srt = 90)
text(filtered.mata[2,]$gompShape,41,"BY4741", col="red",pos = 4, cex = 0.8, srt = 90)
#lines(density(filtered.mata.data$gompShape),na.rm=TRUE)
quantiles.mata<-quantile(filtered.mata$gompShape, c(.05,.10,.15,.85,.90,.95))
#lines(A)
quantiles.mata<-unname(quantiles.mata)
percent.5<-quantiles.mata[[1]]
percent.10<-quantiles.mata[[2]]
percent.15<-quantiles.mata[[3]]
percent.85<-quantiles.mata[[4]]
percent.90<-quantiles.mata[[5]]
percent.95<-quantiles.mata[[6]]

segments(percent.5,45,percent.5,0, col= "magenta",lwd=2,lty=2)
text(percent.5,45,"5%",pos = 3, cex = 0.8, srt = 90)

segments(percent.10,50,percent.10,0, col= "magenta",lwd=2,lty=2)
text(percent.10,47,"10%",pos = 3, cex = 0.8, srt = 90)

segments(percent.15,60,percent.15,0, col= "magenta",lwd=2,lty=2)
text(percent.15,60,"15%",pos = 3, cex = 0.8, srt = 90)

segments(percent.85,60,percent.85,0, col= "royalblue",lwd=2,lty=2)
text(percent.85,60,"85%",pos = 3, cex = 0.8, srt = 90)

segments(percent.90,50,percent.90,0, col= "royalblue",lwd=2,lty=2)
text(percent.90,47,"90%",pos = 3, cex = 0.8, srt = 90)

segments(percent.95,45,percent.95,0, col= "royalblue",lwd=2,lty=2)
text(percent.95,45,"95%",pos = 3, cex = 0.8, srt = 90)

dev.off()


#for (current.percent in c(percent.5,percent.10,percent.15,percent.85,percent.90,percent.95)){
#tmp =  filtered.mata[filtered.mata$gompShape<current.percent,]

#genotype.name<-tmp$genotype

#genotype.number<-length(genotype.name)

#}

MAT_a_percent5<- filtered.mata[filtered.mata$gompShape<percent.5,]
MAT_a_percent5<-length(MAT_a_percent5[,1])
MAT_a_percent10<- filtered.mata[filtered.mata$gompShape<percent.10,]
MAT_a_percent10<-length(MAT_a_percent10[,1])
MAT_a_percent15<- filtered.mata[filtered.mata$gompShape<percent.15,]
MAT_a_percent15<-length(MAT_a_percent15[,1])
MAT_a_percent85<- filtered.mata[filtered.mata$gompShape<percent.85,]
MAT_a_percent85<-length(MAT_a_percent85[,1])
MAT_a_percent90<- filtered.mata[filtered.mata$gompShape<percent.90,]
MAT_a_percent90<-length(MAT_a_percent90[,1])
MAT_a_percent95<- filtered.mata[filtered.mata$gompShape<percent.95,]
MAT_a_percent95<-length(MAT_a_percent95[,1])


#set up a skeleton table

condition_names = c("MAT_a","MAT_a_genotype","MAT_a_WT","MAT_a_percent5","MAT_a_percent10",
                    "MAT_a_percent15","MAT_a_percent85","MAT_a_percent90","MAT_a_percent95",
                    "MAT_alpha","MAT_alpha_genotype","MAT_alpha_WT","MAT_alpha_percent5","MAT_alpha_percent10",
                    "MAT_alpha_percent15","MAT_alpha_percent85","MAT_alpha_percent90","MAT_alpha_percent95")



conditions = c(MAT_a,MAT_a_genotype,MAT_a_WT,MAT_a_percent5,MAT_a_percent10,
               MAT_a_percent15,MAT_a_percent85,MAT_a_percent90,MAT_a_percent95,
               MAT_alpha,MAT_alpha_genotype,MAT_alpha_WT,MAT_alpha_percent5,MAT_alpha_percent10,
               MAT_alpha_percent15,MAT_alpha_percent85,MAT_alpha_percent90,MAT_alpha_percent95)


df=data.frame(condition_names,conditions)


write.csv(df,"filtered_condLongevity.csv")








break()

llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H2  = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb$wt, lifespan2=tb$wtCR );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );


