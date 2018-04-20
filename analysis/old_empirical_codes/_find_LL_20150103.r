setwd("~/projects/0.network.aging.prj/9.ken")
rm(list=ls())

list.files(pattern='csv')
#get sce ORF and names
tb.sce = read.csv("SceORF_name.csv",header=F)



tb = read.csv("rls.csv", colClasses=c("integer", rep("character",8),
                                      rep("numeric",5), rep("character",8),
                                      rep("numeric",5),'character',
                                      rep('numeric',3), 'character')
              )
str(tb)

tb$setRLS_vs_ref_RLS = tb$set_lifespan_mean / tb$ref_lifespan_mean
hist(tb$setRLS_vs_ref_RLS, br = 30)
summary(tb$ranksum_p)
tb$LL = ifelse( tb$ranksum_p <0.05 & tb$setRLS_vs_ref_RLS>1, 'LL', 'NLL' )
table(tb$LL)

sort( unique(toupper(tb$set_genotype[tb$LL=='LL'])) 

# need to pick up single-gene mutants from tb$set_genotype 
tb$single_gene_mutant = NA;
for( i in 1:length(tb[,1])) {
  
  m
  regmatches(x, m)
   
}
#tb$single_gene_mutant[grep( "[\\s|\\/]+", tb$set_genotype)] = 2
table(tb$single_gene_mutant)
summary(tb$single_gene_mutant)

tb.s = tb[tb$single_gene_mutant == 1, ] 
#tb.s contains only single element
sort( unique(toupper(tb.s$set_genotype[tb.s$LL=='LL'])) )
head(tb.s)




