#example spearman rank analyses eddy 2022 sep
#set working directory to where you want to run things
setwd('~/OneDrive - University of Otago/CRISPR masters project - Anya/R_stuff/')
#read in table (check on excell sheet as to how I set things up for ease of use notice the removal of weird characters etc)
test<-read.table('example_input.txt',header=T,row.names=NULL,sep='\t')
head(test)

#spearmank rank test#
#real easy test to run you have two variables an x and a y 
corr <- cor.test(x=test$MC1_biomass_dryweight, y=test$MC1_readcounts_metabar, method = 'spearman')
corr <- cor.test(x=test$MC1_biomass_dryweight, y=test$MC1_readcounts_metabar, method = 'pearson')
corr
#so spearman rank is a non-parametric test while a pearson is a parametric test....Im currently having a debatein my mind as to which is better. But this allows you to do both while I mull that one. 

#linear model#
#install tidyverse if required if installed skip
#install.packages('tidyverse')
library(tidyverse)
#draw it out:
ggplot(test, aes(x=MC1_biomass_dryweight, y=MC1_readcounts_metabar)) + geom_point(color='BLACK', size = 4) + 
geom_smooth(method=lm,color='RED')+
theme_bw()  
#fit the linear model:
lm_test<-lm(MC1_biomass_dryweight~MC1_readcounts_metabar,data=test)
summary(lm_test)
#here is a good explanation of variables and what to look for https://www.statology.org/lm-function-in-r/
#now for whatever reason its a faff to get the r2 values in a ggplot (which look nicer than a standard R plot) but it is posible https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph 

#and that is it, easy right!