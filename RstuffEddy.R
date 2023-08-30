#example spearman rank analyses eddy 2022 sep
#set working directory to where you want to run things
setwd("~/OneDrive - University of Otago/CRISPR masters project - Anya/R_stuff/")
#read in table (check on excell sheet as to how I set things up for ease of use notice the removal of weird characters etc)
test<-read.table('example_input.txt',header=T,row.names=NULL,sep='\t')
head(test)
meta_80<-read.table('meta_OTU_80%.txt',header=T,row.names=NULL,sep='\t')
meta_98<-read.table('meta_OTU_98%.txt',header=T,row.names=NULL,sep='\t')
crispr_80<-read.table('crispr_OTU_80%.txt',header=T,row.names=NULL,sep='\t')
crispr_98<-read.table('crispr_OTU_98%.txt',header=T,row.names=NULL,sep='\t')

#what are you reads here?
#Can you do normalised read counts
S#normalised read counts = number of reads per sample/divided by the total number of reads per sample (exluding negative controls)

#spearmank rank test#
#real easy test to run you have two variables an x and a y 
#make sure you are changing both the x and y variables in your tests
#here I am testing MC1 of crispr 80 against the biomass estimates for MCI
ccorr
?cor.test
sp_crispr_80_sm <- cor.test(x=crispr_80$MC1, y=crispr_80$MC1_biomass_dryweight, method = 'spearman')
sp_crispr_80_sm
#which gives:
#Spearman's rank correlation rho
#
#data:  crispr_80$MC1 and crispr_80$MC1_biomass_dryweight
#S = 66, p-value = 0.7131
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.1785714 
#here there is nosignificance (p-value >>0.05), so there is no relationship between biomass and reads

#we dont have the smaple sizes for pearson (spearman is less impacted to outliers)
sp_crispr_80_pr <- cor.test(x=crispr_80$MC1, y=crispr_80$MC1_biomass_dryweight, method = 'pearson')
sp_crispr_80_pr
#another option if kendall which can handle smaller sample sizes, it similar to spearman 
sp_crispr_80 <- cor.test(x=crispr_80$MC1, y=crispr_80$MC1_biomass_dryweight, method = 'kendall')
sp_crispr_80_pr
#also none significant

#second test is a linear model - which should give basically the same answer as pearson/spearman
#which should be fine 
lm_test<-lm(MC1~MC1_biomass_dryweight,data=crispr_80)
summary(lm_test)
#the benefit being that you can graph it out:
library(tidyverse)
#draw it out:
ggplot(crispr_80, aes(x=MC1_biomass_dryweight,y=MC1)) + geom_point(color='BLACK', size = 4) + 
  geom_smooth(method=lm,color='RED')+
  theme_bw()  

#to get the equation/r2 on the graph use this library:
library(ggpmisc)
ggplot(crispr_80, aes(x=MC1_biomass_dryweight,y=MC1)) + geom_point(color='BLACK', size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()  

#create a table that is:
#test spearmanrho spearmanp-value kendaltau kendallp-value lm_eq lm_r2


#now this is the bit Im not 100% about
crispr_80_MC1<-cbind(crispr_80$Species,crispr_80$MC1,crispr_80$MC1_biomass_dryweight,'MC1')
crispr_80_MC1
colnames(crispr_80_MC1) <- c('Species','Read_count','biomass','test')

crispr_80_MC2<-cbind(crispr_80$Species,crispr_80$MC2,crispr_80$MC2_biomass_dryweight,'MC2')
crispr_80_MC2
colnames(crispr_80_MC2) <- c('Species','Read_count','biomass','test')

crispr_80_MC3<-cbind(crispr_80$Species,crispr_80$MC3,crispr_80$MC3_biomass_dryweight,'MC3')
crispr_80_MC3
colnames(crispr_80_MC3) <- c('Species','Read_count','biomass','test')

crispr_80_MC4<-cbind(crispr_80$Species,crispr_80$MC4,crispr_80$MC4_biomass_dryweight,'MC4')
crispr_80_MC4
colnames(crispr_80_MC4) <- c('Species','Read_count','biomass','test')

crispr_80_MC5<-cbind(crispr_80$Species,crispr_80$MC5,crispr_80$MC5_biomass_dryweight,'MC5')
crispr_80_MC5
colnames(crispr_80_MC5) <- c('Species','Read_count','biomass','test')

test<-rbind(crispr_80_MC1,crispr_80_MC2,crispr_80_MC3,crispr_80_MC4,crispr_80_MC5)
test<-as.data.frame(test)



#this part not working yet dont run this yet
lm_test<-lm(Read_count~biomass,data=test)
summary(lm_test)
library(tidyverse)
#draw it out:
ggplot(test, aes(x=biomass,y=Read_count)) + geom_point(color='BLACK', size = 4) + 
  geom_smooth(method=lm,color='RED')+
  theme_bw()  

ggplot(test, aes(x=biomass,y=Read_count)) + geom_point(color='BLACK', size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()  

#probs need to be a GLM with biomass~ read count + species + test

