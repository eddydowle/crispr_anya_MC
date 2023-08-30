setwd("~/OneDrive - University of Otago/CRISPR masters project - Anya/R_stuff/")
#now switching to 98% OTU tables
#OTU tables are:
#Species  MC1 MC2 MC3 MC4 MC5 MC1_biomass_dryweight .... MC1_biomass_dryweight

#species are species names of animals present in a sample
#MC1 MC2 .. MC5 are read counts per sample for each species
#e.g. mock community 1 (MC1) has standardised read count (read count/total read count) of 0.22 for species Aoteapsyche
#MC1_biomass_dryweights are the biomass estimates for each of the species in the various pools
#e.g. MC1_biomass_dryweights 7.44 grams of Aoteasphye was put into the MC1 pool
#the same pools were used in both methods so dry weights are identical for each method
#two methods were used to generate read counts: metabarcoding and crispr enrichment
#neither worked particularily well
#going back to generate raw genetic estimates from shot gun sequencing to test dry-weight to raw genetic material estimates is next step in lab (probably)

normmeta_80<-read.table('meta_normalisedreads_98%.txt',header=T,row.names=NULL,sep='\t')
normcrispr_80_COIonly <-read.table('crispr_normalisedreads_98%.txt',header=T,row.names=NULL,sep='\t')

# convert the crspr data first
crispr <- normcrispr_80_COIonly %>% 
  # clean up species names
  mutate(Species = str_trim(Species)) %>% 
  # make long
  pivot_longer(-Species) %>% 
  mutate(
    # extract sample number
    sample = as.integer(str_extract(name, "[0-9]")),
    # change name to reflect data kind
    name = case_when(
      str_detect(name, "biomass") ~ "biomass",
      TRUE ~ "read_count"
    )
  ) %>% 
  # pivot wider
  pivot_wider() %>% 
  # add method tag
  bind_cols(method = "Cas_enriched")
# do the same for meta, and combine
dat <- normmeta_80 %>% 
  mutate(Species = str_trim(Species)) %>% 
  pivot_longer(-Species) %>% 
  mutate(
    sample = as.integer(str_extract(name, "[0-9]")),
    name = case_when(
      str_detect(name, "biomass") ~ "biomass",
      TRUE ~ "read_count"
    )
  ) %>% 
  pivot_wider() %>% 
  bind_cols(method = "Metabarcoding") %>% 
  bind_rows(crispr)

#something weird going on with austrosimulium merge 
dat$Species <- gsub("Austrosimulium sp.", "Austrosimulium sp", dat$Species)

# data consistency: species look OK after trimming above
unique(dat$Species)

# eyeball the data with plots
#by species
dat %>% 
  ggplot(aes(read_count, biomass, colour = method)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", se = F) +
  facet_wrap(~ Species, scales = "free")
# make a more convenient wide format for regression analysis
dat_wide <- dat %>% 
  pivot_wider(names_from = method, values_from = read_count) %>%
  # data is non-negative, derived from counts, and on plots shows heteroscedasticity
  # so log-transforming biomass and read counts for lm, thus avoiding the need to mess
  # with GLM functions
  mutate(across(c(biomass, Metabarcoding, Cas_enriched), log))


x<-dat %>% 
  pivot_wider(names_from = method, values_from = read_count)
x$log_biomass<-log(x$biomass)
# now we can make linear models, allowing for interaction between species and read
# count predictors of biomass
lm_meta_sp <- lm(biomass ~ Metabarcoding * Species, dat_wide)
# some significant interaction terms
summary(lm_meta_sp)
plot_model(lm_meta_sp,type = "pred", terms = c("Metabarcoding", "Species"))
library(interactions)
interact_plot(lm_meta_sp, pred = Metabarcoding, modx = Species, plot.points = TRUE)
qqnorm(residuals(lm_meta_sp))

# the same for crispr
lm_crispr_sp <- lm(biomass ~ Cas_enriched * Species, dat_wide)
plot_model(lm_crispr_sp,type = "pred", terms = c("Cas_enriched", "Species"))
interact_plot(lm_crispr_sp, pred = Cas_enriched, modx = Species, plot.points = TRUE)
summary(lm_crispr_sp)
qqnorm(residuals(lm_crispr_sp))
plot(lm_meta_sp)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(lm_meta_sp)
tab_model(lm_crispr_sp)

#comparison of models?
#anova(lm_meta_sp, lm_crispr_sp) no cause different underlying datasets

#within mock community test (no difference between 98% and 80% ID thresholds and none are significant in pearson or kendall)
#are switching to 98% OTU tables but cant be bother renaming files in R
#are sample numbers are not great for this 7 species per community
# comparison of meta or cas to biomass for each mock community. This stuff I understand, nothing was significant
#metabarcoding MC1
normmeta_80_sm <- cor.test(x=normmeta_80$MC1, y=normmeta_80$MC1_biomass_dryweight, method = 'spearman')
normmeta_80_sm
normmeta_80_k <- cor.test(x=normmeta_80$MC1, y=normmeta_80$MC1_biomass_dryweight, method = 'kendall')
normmeta_80_k
lm_normmeta_80<-lm(MC1~MC1_biomass_dryweight,data=normmeta_80)
summary(lm_normmeta_80)
library(ggpmisc)
ggplot(normmeta_80, aes(x=MC1_biomass_dryweight,y=MC1)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 
#crispr MC1
normcrispr_80_sm <- cor.test(x=normcrispr_80_COIonly$MC1, y=normcrispr_80_COIonly$MC1_biomass_dryweight, method = 'spearman')
normcrispr_80_sm
normcrispr_80_k <- cor.test(x=normcrispr_80_COIonly$MC1, y=normcrispr_80_COIonly$MC1_biomass_dryweight, method = 'kendall')
normcrispr_80_k
lm_normcrispr_80<-lm(MC1~MC1_biomass_dryweight,data=normcrispr_80_COIonly)
summary(lm_normcrispr_80)
ggplot(normcrispr_80_COIonly, aes(x=MC1_biomass_dryweight,y=MC1)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()

#metabarcoding MC2
normmeta_80_MC2_sm <- cor.test(x=normmeta_80$MC2, y=normmeta_80$MC2_biomass_dryweight, method = 'spearman')
normmeta_80_MC2_sm
normmeta_80_MC2_k <- cor.test(x=normmeta_80$MC2, y=normmeta_80$MC2_biomass_dryweight, method = 'kendall')
normmeta_80_MC2_k
lm_normmeta_80_MC2<-lm(MC2~MC2_biomass_dryweight,data=normmeta_80)
summary(lm_normmeta_80_MC2)
ggplot(normmeta_80, aes(x=MC2_biomass_dryweight,y=MC2)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 
#crispr MC2
normcrispr_80_MC2_sm <- cor.test(x=normcrispr_80_COIonly$MC2, y=normcrispr_80_COIonly$MC2_biomass_dryweight, method = 'spearman')
normcrispr_80_MC2_sm
normcrispr_80_MC2_k <- cor.test(x=normcrispr_80_COIonly$MC2, y=normcrispr_80_COIonly$MC2_biomass_dryweight, method = 'kendall')
normcrispr_80_MC2_k
lm_normcrispr_80_MC2<-lm(MC2~MC2_biomass_dryweight,data=normcrispr_80_COIonly)
summary(lm_normcrispr_80_MC2)
ggplot(normcrispr_80_COIonly, aes(x=MC2_biomass_dryweight,y=MC2)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 

#metabarcoding MC3
normmeta_80_MC3_sm <- cor.test(x=normmeta_80$MC3, y=normmeta_80$MC3_biomass_dryweight, method = 'spearman')
normmeta_80_MC3_sm
normmeta_80_MC3_k <- cor.test(x=normmeta_80$MC3, y=normmeta_80$MC3_biomass_dryweight, method = 'kendall')
normmeta_80_MC3_k
lm_normmeta_80_MC3<-lm(MC3~MC3_biomass_dryweight,data=normmeta_80)
summary(lm_normmeta_80_MC3)
ggplot(normmeta_80, aes(x=MC3_biomass_dryweight,y=MC3)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 
#crispr MC3
normcrispr_80_MC3_sm <- cor.test(x=normcrispr_80_COIonly$MC3, y=normcrispr_80_COIonly$MC3_biomass_dryweight, method = 'spearman')
normcrispr_80_MC3_sm
normcrispr_80_MC3_k <- cor.test(x=normcrispr_80_COIonly$MC3, y=normcrispr_80_COIonly$MC3_biomass_dryweight, method = 'kendall')
normcrispr_80_MC3_k
lm_normcrispr_80_MC3<-lm(MC3~MC3_biomass_dryweight,data=normcrispr_80_COIonly)
summary(lm_normcrispr_80_MC3)
ggplot(normcrispr_80_COIonly, aes(x=MC3_biomass_dryweight,y=MC3)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()

#metabarcoding MC4
normmeta_80_MC4_sm <- cor.test(x=normmeta_80$MC4, y=normmeta_80$MC4_biomass_dryweight, method = 'spearman')
normmeta_80_MC4_sm
normmeta_80_MC4_k <- cor.test(x=normmeta_80$MC4, y=normmeta_80$MC4_biomass_dryweight, method = 'kendall')
normmeta_80_MC4_k
lm_normmeta_80_MC4<-lm(MC4~MC4_biomass_dryweight,data=normmeta_80)
summary(lm_normmeta_80_MC4)
ggplot(normmeta_80, aes(x=MC4_biomass_dryweight,y=MC4)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 
#crispr MC4
normcrispr_80_MC4_sm <- cor.test(x=normcrispr_80_COIonly$MC4, y=normcrispr_80_COIonly$MC4_biomass_dryweight, method = 'spearman')
normcrispr_80_MC4_sm
normcrispr_80_MC4_k <- cor.test(x=normcrispr_80_COIonly$MC4, y=normcrispr_80_COIonly$MC4_biomass_dryweight, method = 'kendall')
normcrispr_80_MC4_k
lm_normcrispr_80_MC4<-lm(MC4~MC4_biomass_dryweight,data=normcrispr_80_COIonly)
summary(lm_normcrispr_80_MC4)
ggplot(normcrispr_80_COIonly, aes(x=MC4_biomass_dryweight,y=MC4)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()

#metabarcoding MC5
normmeta_80_MC5_sm <- cor.test(x=normmeta_80$MC5, y=normmeta_80$MC5_biomass_dryweight, method = 'spearman')
normmeta_80_MC5_sm
normmeta_80_MC5_k <- cor.test(x=normmeta_80$MC5, y=normmeta_80$MC5_biomass_dryweight, method = 'kendall')
normmeta_80_MC5_k
lm_normmeta_80_MC5<-lm(MC5~MC5_biomass_dryweight,data=normmeta_80)
summary(lm_normmeta_80_MC5)
ggplot(normmeta_80, aes(x=MC5_biomass_dryweight,y=MC5)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw() 
#crispr MC5
normcrispr_80_MC5_sm <- cor.test(x=normcrispr_80_COIonly$MC5, y=normcrispr_80_COIonly$MC5_biomass_dryweight, method = 'spearman')
normcrispr_80_MC5_sm
normcrispr_80_MC5_k <- cor.test(x=normcrispr_80_COIonly$MC5, y=normcrispr_80_COIonly$MC5_biomass_dryweight, method = 'kendall')
normcrispr_80_MC5_k
lm_normcrispr_80_MC5<-lm(MC5~MC5_biomass_dryweight,data=normcrispr_80_COIonly)
summary(lm_normcrispr_80_MC5)
ggplot(normcrispr_80_COIonly, aes(x=MC5_biomass_dryweight,y=MC5)) + geom_point(aes(color = Species), size = 4) + 
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  theme_bw()

#could do a dis-simularity matrix analysis kulczynski Dissimilarity Index or something
#not sure if there is any point




