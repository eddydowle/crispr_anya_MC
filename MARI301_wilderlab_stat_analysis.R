##################################################
## DATA ANALYSIS: MARI301 - Wilderlab eDNA Data ##
##################################################

###########################
## Importing data into R ##
###########################

## LOAD PACKAGES
library('vegan')
library('tidyverse')
library('BiodiversityR')
library('ggsci')
library('readxl')
library('dplyr')
library('reshape2')
library('viridis')
library('rstatix')
library('labdsv')
library('gplots')

## SET WORKING DIRECTORY
setwd("~/Dropbox (Otago University)/BiomeGene/")

## FORMAT DATA
raw_data <- read_excel('WLJ603943.xlsx', sheet = 3)
frequency_table <- read_excel('WLJ603943.xlsx', sheet = 2)
metadata <- read_excel('WLJ603943.xlsx', range = "metadata!A19:J34") %>% mutate(Location = case_when(startsWith(ClientSampleID, 'Saltmarsh') ~ 'Saltmarsh',
                                                                                                     startsWith(ClientSampleID, 'Careys') ~ 'Careys Bay',
                                                                                                     startsWith(ClientSampleID, 'Cockle') ~ 'Cockle Bed',
                                                                                                     startsWith(ClientSampleID, 'Open') ~ 'Open Water',
                                                                                                     startsWith(ClientSampleID, 'Har') ~ 'Harrington Point'))

# remove unused columns and set rownames
#frequency_table <- frequency_table %>% filter(!is.na(ScientificName)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName') %>% select (-c(Rank, TaxID, CommonName, Group))
#raw_data <- raw_data %>% remove_rownames %>% column_to_rownames(var = 'Sequence') %>% select (-c(Target, ScientificName, Rank, TaxID, CommonName, Group)) 
#metadata <- metadata %>% select (-c(DeploymentDuration, ClientNotes))


##################################
## Preliminary data exploration ##
##################################

## RAREFACTION CURVES: determine if sufficient sequencing is achieved.
# prepare dataframe for analysis
rarefaction_df <- raw_data %>% select (-c(Target, ScientificName, Rank, TaxID, CommonName, Group, Sequence)) %>% t()
# identify lowest number of reads for samples
raremax_rarefaction <- min(rowSums(rarefaction_df))
# generate graph
rarecurve(rarefaction_df, step = 100, sample = raremax_rarefaction, col = 'blue', cex = 0.6,label=FALSE)
?rarecurve
getAnywhere(rarecurve)
body(rarecurve)
# what is the meaning of the "Warning message"?

#my input from python code
x<-read.table('BCI_mycode_rarefraction.txt',sep='\t',header=T)
m.av<-rollmean(x$x, 3,fill = list(NA, NULL, NA))
x$movave<-m.av

x %>%ggplot(aes(x=y, y=x, group=Sample)) + geom_line() + 
  geom_line(aes(movave, group=Sample),color="red")
x<-read.table('wilderlab_rarefraction_martyversion.txt',sep='\t',header=T)

x %>%
  ggplot( aes(x=x, y=y, group=Sample, color=Sample)) +
  geom_line()

x<-read.table('wilderlab_rarefraction.txt',sep='\t',header=T)

x %>%
  ggplot( aes(x=x, y=y, group=Sample, color=Sample)) +
  geom_line()


#BCI
BCI<-read.table('BCI.csv',sep='\t',header=T)
raremax_rarefaction <- min(rowSums(BCI))
# generate graph
rarecurve(BCI, step = 100, sample = raremax_rarefaction, col = 'blue', cex = 0.6,label=FALSE)



## SPECIES ACCUMULATION CURVES: determine if sufficient sampling was conducted.
# prepare dataframe for analysis
speccurve <- frequency_table %>% filter(!is.na(ScientificName)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName') %>% select (-c(Rank, TaxID, CommonName, Group)) %>% t() 
# set theme of figure
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())
# calculate values
Accum.1 <- accumcomp(speccurve, y=metadata, factor='Location', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)
# plot graph
plotgg1 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE) + 
  BioR.theme +
  scale_colour_npg() +
  labs(x = "Sample", y = "Number of species", colour = "source", shape = "source")
plotgg1


## TAXON ABUNDANCE: determine the most abundant taxa in our data.
# prepare dataframe for analysis
taxonab <- frequency_table %>% select(-c(TaxID, Rank, CommonName, Group)) %>% mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% arrange(desc(Total)) %>% slice_max(n = 20, order_by = Total) %>% select(-c(Total)) %>% melt()
# plot data
ggplot(taxonab, aes(fill = ScientificName, y = value, x = variable)) + 
  geom_bar(position = 'stack', stat = 'identity')
# plot relative abundance
ggplot(taxonab, aes(fill = ScientificName, y = value, x = variable)) + 
  geom_bar(position = 'fill', stat = 'identity')
# filter data for only species and plot again
speciesab <- frequency_table %>% filter(Rank == 'species') %>% select(-c(TaxID, Rank, CommonName, Group)) %>% mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% arrange(desc(Total)) %>% slice_max(n = 20, order_by = Total) %>% select(-c(Total)) %>% melt()
ggplot(speciesab, aes(fill = ScientificName, y = value, x = variable)) + 
  geom_bar(position = 'fill', stat = 'identity')
# combine data per site and plot again
var_names <- metadata %>% select(ClientSampleID, UID) %>% deframe()
siteab <- frequency_table %>% rename(!!!var_names) %>% filter(Rank == 'species') %>% select(-c(TaxID, Rank, CommonName, Group)) %>% mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% arrange(desc(Total)) %>% slice_max(n = 20, order_by = Total) %>% select(-c(Total)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName')
names(siteab) <- substring(names(siteab), 1, 10)
sitesum <- melt(t(rowsum(t(siteab), group = colnames(siteab), na.rm = TRUE)))
ggplot(sitesum, aes(fill = Var1, y = value, x = Var2)) + 
  geom_bar(position = 'fill', stat = 'identity')
# plot per species for easy visualization
ggplot(sitesum, aes(fill = Var2, y = value, x = Var2)) +
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_viridis(discrete = T, option = 'E') +
  ggtitle('Graphing data per species') +
  facet_wrap(~Var1) +
  theme(legend.position = 'none') +
  xlab('')

## VENN DIAGRAM: display overlap of taxa between locations
# prepare dataframe for analysis
venntable <- frequency_table %>% filter(!is.na(ScientificName)) %>% select (-c(Rank, TaxID, CommonName, Group)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName') %>% t() %>% as.data.frame() %>% rownames_to_column() #%>% group_by(metadata$Location) %>% summarise(across(everything(), sum), .groups = 'drop') %>% as.data.frame()
# sum values per location
venntable['Location'] <- metadata$Location
venntable <- venntable %>% select(-c(rowname)) %>% group_by(Location) %>% summarise(across(everything(), sum)) %>% as.data.frame() %>% remove_rownames %>% column_to_rownames(var = 'Location')
# list which species at which location
venntable <- replace(venntable, venntable >0, NA)
careysbay <- colnames(venntable[1,sapply(venntable[1,], is.na)])
cocklebed <- colnames(venntable[2,sapply(venntable[2,], is.na)])
harringtonpoint <- colnames(venntable[3,sapply(venntable[3,], is.na)])
openwater <- colnames(venntable[4,sapply(venntable[4,], is.na)])
saltmarsh <- colnames(venntable[5, sapply(venntable[5,], is.na)])
# create venn diagram
venn.table <- venn(list(careysbay = careysbay, cocklebed = cocklebed, harringtonpoint = harringtonpoint, openwater = openwater, saltmarsh = saltmarsh))

##############################
## Alpha Diversity Analysis ##
##############################

## SPECIES RICHNESS COMPARISON
# calculate the number of taxa detected per sample and group per sampling location
var_names <- metadata %>% select(ClientSampleID, UID) %>% deframe()
specrich <- frequency_table %>% rename(!!!var_names) %>% filter(!is.na(ScientificName)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName') %>% select (-c(Rank, TaxID, CommonName, Group))
specrich <- (specrich > 0) *1L
specrichtable <- setNames(nm = c('colname', 'SpeciesRichness'), stack(colSums(specrich))[2:1])
specrichtable <- transform(specrichtable, group = substr(colname, 1, 10))
# test assumptions of statistical test, first normal distribution, next homoscedasticity
histogram(~ SpeciesRichness | group, data = specrichtable, layout = c(5,1))
bartlett.test(SpeciesRichness ~ group, data = specrichtable)
# due to low sample number per location and non-normal distribution, run Welch's ANOVA and visualize using boxplot
oneway.test(SpeciesRichness ~ group, data = specrichtable, var.equal = FALSE)
specrichtable %>% games_howell_test(SpeciesRichness ~ group)
boxplot(SpeciesRichness ~ group, data = specrichtable, names = c('Careys Bay', 'Cockle bed', 'Harrington', 'Open water', 'Saltmarsh'), ylab = 'Species Richness')


##############################
## Beta Diversity Analysis ##
##############################

## ORDINATION ANALYSIS
# prepare the dataframe
var_names <- metadata %>% select(ClientSampleID, UID) %>% deframe()
orditable <- frequency_table %>% rename(!!!var_names) %>% filter(!is.na(ScientificName)) %>% remove_rownames %>% column_to_rownames(var = 'ScientificName') %>% select (-c(Rank, TaxID, CommonName, Group)) %>% t()
orditable <- (orditable > 0) *1L
# run analysis
dis <- vegdist(orditable, method = 'jaccard')
groups <- factor(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)), labels = c('Saltmarsh', 'Careys Bay', 'Cockle bed', 'Open Water', 'Harrington Point'))
mod <- betadisper(dis, groups)
# plot ordination and eigenvalues
plot(mod, hull = FALSE, ellipse = TRUE)
ordination_mds <- wcmdscale(dis, eig = TRUE)
ordination_eigen <- ordination_mds$eig
ordination_eigenvalue <- ordination_eigen/sum(ordination_eigen)
ordination_eigen_frame <- data.frame(Inertia = ordination_eigenvalue*100, Axes = c(1:14))
plot(ordination_mds)
eigenplot <- ggplot(data = ordination_eigen_frame, aes(x = factor(Axes), y = Inertia)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
  theme_classic() +
  xlab("Axes") +
  ylab("Inertia %") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
eigenplot

## PERMANOVA
source <- adonis(orditable ~ Location, data = metadata, by = 'terms')
print(source$aov.tab)

## PERMDISP
boxplot(mod)
set.seed(25)
permutest(mod)

## INDICATOR SPECIES ANALYSIS
specrichtable$groupnumeric <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)
ISA <- indval(orditable, clustering = as.numeric(specrichtable$groupnumeric))
summary(ISA)
