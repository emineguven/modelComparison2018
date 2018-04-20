
#####################################################################################################
#################install and declare needed packaqes for survival analysis

#install.packages('survival') 
#install.packages('KMsurv')
#install.packages('data.table')
#install.packages('RSQLite')

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
File = 'rls.db'
setwd(Folder)
drv <- SQLite()
con <- dbConnect(drv, dbname = File)

##########################################################################################################
##########################################################################################################

######################Create a complete table of conditions that have data in the RLS database. 
######################(Each row in this table will be a unique combination of genotype, mating type, media, and temperature

###get unique conditions from "set" columns (experimental treatments/mutations)
conditions1 = dbGetQuery(con, "
  SELECT DISTINCT set_genotype as genotype, set_mating_type as mat, set_media as media, set_temperature as temp
  FROM result
  WHERE pooled_by = 'file' 	
  ")
### rows in the result table of the database are not mutually exclusive. In some rows, data has sometimes been pooled by genotype, background strain, etc.   

###get unique conditions from "reference" columns (these columns are the control conditions/lifespan results for each row of experimental lifespan results. Rows are not mutually exclusive (1 control to many experimental conditions)
conditions2 = dbGetQuery(con, "
  SELECT DISTINCT ref_genotype as genotype, ref_mating_type as mat, ref_media as media, ref_temperature as temp
  FROM result 
  WHERE pooled_by = 'file'
  ")

####combine and take unique conditions from these two
conditions = rbind(conditions1, conditions2)
conditions = unique(conditions)
conditions = conditions[complete.cases(conditions),]
row.names(conditions) = NULL													###renumber the rows. important because future processes will refer to a unique condition by its row number in this table
controlConditions = conditions[conditions$genotype %in% c('BY4741', 'BY4742', 'BY4743','sir2','tor1'),]			###create a table of conditions that have WT genotypes	

###Add columns to the conditions data frame. These columns will be filled in by their respective variable: ie. gompertz shape/rate of the lifespan data associated with a given conditions (genotype, mating type, media, temp)
controlConditions$n = apply(controlConditions, 1, function(row) 0)
controlConditions$avgLS = apply(controlConditions, 1, function(row) 0)
controlConditions$StddevLS = apply(controlConditions, 1, function(row) 0)
controlConditions$medianLS = apply(controlConditions, 1, function(row) 0)
controlConditions$gompShape = apply(controlConditions, 1, function(row) 0)
controlConditions$gompRate = apply(controlConditions, 1, function(row) 0)
controlConditions$gompLogLik = apply(controlConditions, 1, function(row) 0)
controlConditions$gompAIC = apply(controlConditions, 1, function(row) 0)
controlConditions$weibShape = apply(controlConditions, 1, function(row) 0)
controlConditions$weibScale = apply(controlConditions, 1, function(row) 0)
controlConditions$weibLogLik = apply(controlConditions, 1, function(row) 0)
controlConditions$weibAIC = apply(controlConditions, 1, function(row) 0)
controlConditions$delta_LL = apply(controlConditions, 1, function(row) 0)
#controlConditions$noise   =  apply(controlConditions, 1, function(row) 0)



fit_names = c( 'genotype','sd.ls','current.seed','medianLS','Delta_LL','AvgLS','gompLogLik','gompRate','gompShape','weibLogLik','weibScale','weibShape')
fit = data.frame(matrix(, nrow=length(fit_names), ncol=2000)) #set up a skeleton table
names(fit) = c( "genotype","sd.ls","current.seed","medianLS","Delta_LL","AvgLS","gompLoglik","gompRate","gompShape","weibLogLik","weibScale","weibShape")





#############################################################################################################
#############################################################################################################

#########################################Loop through the conditions to get lifespan data 
for (r in 1:length(controlConditions$genotype)) {
  #i=0.2
  
  
  genotypeTemp = controlConditions$genotype[r]
  mediaTemp = controlConditions$media[r]
  temperatureTemp = controlConditions$temp[r]
  matTemp = controlConditions$mat[r]
  
  conditionName = apply(controlConditions[r,1:4], 1, paste, collapse=" ") 				#### create a string to name a possible output file		
  conditionName = gsub("[[:punct:]]", "", conditionName) 						#### remove special characters from the name (ie. quotations marks, slashes, etc.)
  
  genotypeTemp = gsub("'", "''", genotypeTemp)
  genotypeTemp = gsub('"', '""', genotypeTemp)
  
  ##### Query the database to take data (including lifespan data) from every mutually exclusive row (pooled by file, not genotype, not background, etc). 
  ##### There will often be multiple rows (representing different experiments) for each unique condition. This analysis pools lifespan data from multiple experiments if the controlConditions are all the same
  queryStatementSet = paste(
    "SELECT * ",
    "FROM result ",
    "WHERE pooled_by = 'file' AND set_genotype = '", genotypeTemp, "' AND set_mating_type = '", matTemp, "' AND set_media = '", mediaTemp, "' AND set_temperature = '", temperatureTemp, "'", 
    sep = ""
  )
  
  queryStatementRef = paste(		### both the reference and set columns will be searched for matching conditions
    "SELECT * ",
    "FROM result ",
    "WHERE pooled_by = 'file' AND ref_genotype = '", genotypeTemp, "' AND ref_mating_type = '", matTemp, "' AND ref_media = '", mediaTemp, "' AND ref_temperature = '", temperatureTemp, "'", 
    sep = ""
  )
  
  dataListSet = dbGetQuery(con, queryStatementSet)
  dataListRef = dbGetQuery(con, queryStatementRef)
  lifespansChar = unique(c(dataListSet$set_lifespans, dataListRef$ref_lifespans)) ### combine lifespan values for a given condition into a single data structure. (problems of having non-mutually exclusive rows in the ref_lifespans column are overcome by only taking unique groups of lifespans. The assumption is that no two experiments produced identical lifespan curvs)
  
  ##### Database codes the lifespan data for each experiment as a string. So, lifespanChar is a vector of strings
  ##### Convert lifespanChar into a single vector of integers lifespansTemp
  lifespansTemp = c()
  if (length(lifespansChar) > 0) {
    for (s in 1:length(lifespansChar)) {
      lifespansNum = as.numeric(unlist(strsplit(lifespansChar[s], ",")))
      lifespansNum = lifespansNum[!is.na(lifespansNum)]
      if (length(lifespansNum) > 0) {
        lifespansTemp = c(lifespansTemp, lifespansNum)
      } 
    }
  }
  
  controlConditions$n[r] = length(lifespansTemp)		#### record number of individuals 
  
  
  ##### Now that there is a single vector of lifespans for a condition. Do the analysis
  if (length(lifespansTemp) > 50) {
    
    
    
    lifespanGomp = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
    lifespanWeib = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'weibull')  
    
    ### Fill in added columns in the controlConditions table with data from gompertz and weibull fitting. Columns are named according to respective variables
    controlConditions$avgLS[r] = mean(lifespansTemp)
    controlConditions$StddevLS[r] = sd(lifespansTemp)
    controlConditions$medianLS[r] = median(lifespansTemp)
    controlConditions$gompShape[r] = lifespanGomp$res[1,1]
    controlConditions$gompRate[r] = lifespanGomp$res[2,1]
    controlConditions$gompLogLik[r] = lifespanGomp$loglik
    controlConditions$gompAIC[r] = lifespanGomp$AIC
    
    controlConditions$weibShape[r] = lifespanWeib$res[1,1]
    controlConditions$weibScale[r] = lifespanWeib$res[2,1]
    controlConditions$weibLogLik[r] = lifespanWeib$loglik
    controlConditions$weibAIC[r] = lifespanWeib$AIC   
    
    controlConditions$delta_LL[r] = lifespanWeib$loglik- lifespanGomp$loglik
    
    #controlConditions<-controlConditions[controlConditions$n > 50,]
    
    #controlConditions$noise[r]=i
    
    for( seed in c(12345)){#}, 20160711, -1881, 9999.1234,300045,50758,-10000,74562,-92345,25434)) {
      set.seed(seed)
      current.seed<-seed
     for (genotype_in in sort(unique(controlConditions$genotype))){
      
          lifespansTemp= lifespansTemp
          
          lifespansTemp[lifespansTemp < 0] <- 0
          lifespansTemp <-floor(lifespansTemp+0.5)
          lifespansTemp <- lifespansTemp[ lifespansTemp != 0 ]
          
          lifespanGomp = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'gompertz') ### Use the flexsurvreg package to fit lifespan data to gompertz or weibull distribution
          lifespanWeib = flexsurvreg(formula = Surv(lifespansTemp) ~ 1, dist = 'weibull')  
          
          ### Fill in added columns in the controlConditions table with data from gompertz and weibull fitting. Columns are named according to respective variables
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
          
          
          genotype=genotype_in
          
          #fit_names = c( 'sd.ls','current.seed','i','medianLS','Delta_LL','AvgLS','gompLogLik','gompRate','gompShape','weibLogLik','weibScale','weibShape')
          
          fit = rbind(fit,c(genotype,StddevLS,current.seed,medianLS,delta_LL,avgLS,gompLogLik,gompRate,gompShape,weibLogLik,weibScale,weibShape))
          #write.csv(fit, file="Results.csv", row.names=F)
          
          fit= fit[!is.na(fit[,1]), ]
          
        }
      }
    } s
  }

new_fit<-fit[,1:12]
write.csv(new_fit, file = "fit.csv", row.names=F)

genotype.els = unlist( unique(new_fit$genotype))
colnum = length(genotype.els)

tmp = unlist( unique(new_fit$noise_scale))
scale.els = tmp[order(tmp)]
rownum = length(scale.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = scale.els
colnames(mat) = genotype.els

sort(new_fit[new_fit$genotype[,1]])
data = new_fit[new_fit[,1]=="BY4741", 6 ]

data<-unlist(data)

heat_mat<-matrix(data,ncol=colnum,nrow=rownum)


#data<-read.csv("conditionsWeibGomp_LL_guven_.csv")


