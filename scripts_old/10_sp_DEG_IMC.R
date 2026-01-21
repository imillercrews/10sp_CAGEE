#### 10 sp bird analysis 
### DEG analysis
# R 4.2.1

## setwd
setwd("/geode2/home/u040/imillerc/Quartz/10_sp")

#### load libraries ####
library(DESeq2)
library(edgeR)
library(ape)
# library(RedRibbon)
library(ComplexUpset)
library(tidyverse)

#### load data ####
### load normalized counts file
data = readxl::read_excel('CAGEE/data/normalizedCounts.xlsx')

### format sample names 
## needs to be in format: species_tissue_sample
## here use the format: species_sex_sample
sample_ids = data_frame(sample_name = data %>% 
                          dplyr::select(-c(GeneName,
                                           GeneDescription,
                                           ZebraFinchProteinID)) %>% 
                          colnames())

# replace tree swallow name from "KYTS" to "TS"
sample_ids = sample_ids %>% 
  mutate(sample_id = str_replace(sample_name,
                                 "KYTS",
                                 "TS"))

# create sex column 
sample_ids = sample_ids %>% 
  mutate(sex = case_when(
    grepl("F", sample_id) ~ "F",
    grepl("M", sample_id) ~ "M"
  ))
# 
# #  seperate into multiple columns
# sample_ids = sample_ids %>% 
#   separate(sample_id,
#            into = c(NA,""))

## extract species data
sample_ids = sample_ids %>% 
  mutate(species = str_sub(sample_id,
                           1,
                           2)
  )

## create long format dataframe
data.long = data %>% 
  dplyr::select(-c(ZebraFinchProteinID,
                   GeneDescription)) %>% 
  t()

# rename columns
colnames(data.long) = data.long[1,]

# remove first row
data.long = data.long %>% 
  as.data.frame() %>% 
  slice(-1)

# convert to numeric
data.long = data.long %>% 
  mutate_all(as.numeric)


# add sample id column
data.long = data.long %>% 
  rownames_to_column("sample_name") 

# combine with sample name
data.long = full_join(sample_ids,
                      data.long)

## Calculate phylogenetic distance
# load in species tree
tree = read.tree('CAGEE/data/10_species_birds_edit.nwk')

# calculate species pairwise distance  
data.dist = cophenetic.phylo(tree) 

# pivot to long pairwise format
data.dist = data.frame(col=colnames(data.dist)[col(data.dist)], 
                       row=rownames(data.dist)[row(data.dist)], 
                       dist=c(data.dist))

# divide distance by 2 to get divergence time
data.dist = data.dist %>% 
  mutate(dist = dist/2)

# get combination of every clade
data.species.pair = combn(unique(data.dist$col),
                   2) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('comp') %>% 
  mutate(comp = as.numeric(comp)) 

# get distance for every pair
data.species.pair.dist = data.species.pair %>% 
  left_join(data.dist %>% 
  na.omit() %>%
  filter(dist != 0)%>% 
    rename(V1 = col) %>% 
    rename(V2 = row) )

# ## graph tree
# library(ggtree)
p = ggtree(tree) +
  geom_nodelab()+
  geom_tiplab() +
  theme_tree2()+
    scale_x_continuous(labels = abs)
revts(p)

  
  
  

# #### BlueBird vs Tree swallow QC ####
# #### create pvclust matrix for all samples
# ### just use blue birds and tree swallows
# ## check data for outliers
# ## create tmp dataframe
# tmp.data = data.long %>% 
#                         filter(species %in% c('BB',
#                                               'TS')) %>% 
#                         select(-c(sample_name,
#                                   sex,
#                                   species)) %>% 
#                         column_to_rownames('sample_id') %>% 
#                         t()
# 
# # calculate the variance for each tissue
# tmp.rv <- rowVars(tmp.data)
# 
# 
# #graph variance for each gene
# data.frame(variance = tmp.rv) %>% 
#   ggplot(aes(variance)) +
#   geom_histogram() +
#   theme_classic()
# ggsave('DEG/figures/BBvsTS/Gene variance histogram.png')
# 
# # select the 1000 genes by variance
# tmp.select <- order(tmp.rv, decreasing=TRUE)[seq_len(min(1000, length(tmp.rv)))]
# 
# # only use selected genes
# tmp.data = tmp.data[tmp.select,]
# 
# ## run pvclust
# tmp.p = pvclust::pvclust(tmp.data,
#                          method.dist="cor",
#                          method.hclust="average",
#                          nboot=1000,
#                          quiet = F)
# # graph
# png("DEG/figures/BBvsTS/pvclust BB vs TS sex.png",
#       height = 10,
#       width = 10,
#       unit = 'in',
#       res = 480)
#   plot(tmp.p)
#   dev.off()
#   
# #### perform a PCA on the data for the selected genes
# pca <- prcomp(t(tmp.data))
#   
# # the contribution to the total variance for each component
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
#   
# # compare PC to behavior data
# # PC 1, 2, 3
# pca.data = sample_ids %>%
#   select(-c(sample_name)) %>% 
#     right_join(pca$x %>% 
#                 as.data.frame() %>% 
#                 rownames_to_column('sample_id')) 
#   
# #graph PCA var
# data.frame(percentVar = percentVar,
#              PC = seq(1:25)) %>% 
#     mutate(PC = paste('PC',
#                       PC,
#                       sep ='')) %>% 
#     ggplot(aes(x = reorder(PC,
#                            -percentVar),
#                y = percentVar)) +
#     geom_point() +
#     geom_segment(aes(x=reorder(PC,
#                                -percentVar), 
#                      xend=reorder(PC,
#                                   -percentVar), 
#                      y=0, 
#                      yend=percentVar)) +
#     theme_classic()
# ggsave('DEG/figures/BBvsTS/PC variance.png')
#   
#   
# #plot pca
# # PC 1 vs PC 2
# pca.data %>%
#     ggplot(aes(x = PC1,
#                y = PC2,
#                color = species,
#                shape = sex)) +
#     geom_point(size=5) +
#     xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
#     ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
#     coord_fixed() +
#     theme_classic()+ 
#     scale_color_manual(values = c('black',
#                                   'red'))
# ggsave('DEG/figures/BBvsTS/PCA 1 vs 2.png',
#        height = 10,
#        width = 10)
#   
# # PC 1 vs PC 2
# pca.data %>%
#   ggplot(aes(x = PC3,
#              y = PC2,
#              color = species,
#              shape = sex)) +
#   geom_point(size=5) +
#   xlab(paste0("PC3: ",round(100*percentVar[3]),"% variance")) +
#   ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
#   coord_fixed() +
#   theme_classic()+ 
#   scale_color_manual(values = c('black',
#                                 'red'))
# ggsave('DEG/figures/BBvsTS/PCA 3 vs 2.png',
#        height = 10,
#        width = 10)
#   
#   
#   
# #### BlueBird vs Tree swallow DEG analysis ####
# # https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# 
# #### pre-processing
# ## create counts data frame
# # just use bluebird and treeswallow
# counts = data.long %>% 
#   filter(species %in% c('BB',
#                         'TS')) %>% 
#   select(-c(sample_name,
#             species,
#             sex)) %>% 
#   column_to_rownames('sample_id') %>% 
#   t()
# 
# ## create DGElist object
# d0 = DGEList(counts)
# 
# ## calculate nomalization factors
# d0 = calcNormFactors(d0)
# 
# ## filter low-expressed genes
# # set cut off
# cutoff <- 5
# # which genes to drop
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# # remove low expressed genes
# d <- d0[-drop,] 
# # number of genes left
# dim(d) 
# 
# ## create sample information
# sample.df = data.frame(sample_id = colnames(counts))
# # add sample data on species and sex
# sample.df = sample.df %>% 
#   left_join(sample_ids )
# # create list of sex
# sample.sex = sample.df %>% 
#   pull(sex)
# # create list of species 
# sample.sex = sample.df %>% 
#   pull(sex)
# # create list of species
# sample.species = sample.df %>% 
#   pull(species)
# # create list of species
# sample.group = sample.df %>% 
#   mutate(group = paste(species,
#                        sex,
#                        sep = '.')) %>% 
#   pull(group)
# 
# 
# ### graph MDS
# png("DEG/figures/BBvsTS/MDS BB vs TS.png",
#     height = 10,
#     width = 10,
#     unit = 'in',
#     res = 480)
# plotMDS(d, 
#         col = as.numeric(as.factor(sample.group)))
# dev.off()
# 
# 
# 
# #### voom transformation and calculation of variance weight
# #### by species 
# ## specify model to be fitted
# mm <- model.matrix(~0+sample.group)
# 
# ### voom transformation
# png("DEG/figures/BBvsTS/voom BB vs TS.png",
#     height = 10,
#     width = 10,
#     unit = 'in',
#     res = 480)
# y <- voom(d, 
#           mm, 
#           plot = T)
# dev.off()
# 
# ## fitting ilnear models in limma
# fit = lmFit(y,
#             mm)
# 
# ### make contrasts
# ## BB
# contr.BB = makeContrasts(sample.groupBB.M - sample.groupBB.F,
#                       levels = colnames(coef(fit)))
# 
# # Estimate contrast for each gene
# tmp.BB <- contrasts.fit(fit, 
#                              contr.BB)
# 
# # Empirical Bayes smoothing 
# tmp.BB <- eBayes(tmp.BB)
# 
# # What genes are most differentially expressed?
# top.table.BB <- topTable(tmp.BB, sort.by = "P", n = Inf)
# 
# # number of DE genes
# length(which(top.table.BB$adj.P.Val < 0.05))
# 
# # graph DEGs
# # volcano plots
# # BB
# top.table.BB %>% 
#   mutate(color = ifelse(adj.P.Val < 0.05,
#                         'sig',
#                         'not_sig')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log(adj.P.Val),
#              color = color)) +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlim(-8,8) +
#   ggtitle('BB male vs female') 
# ggsave('DEG/figures/BBvsTS/BB M vs F volcano.png')
# 
# ## TS
# contr.TS = makeContrasts(sample.groupTS.M - sample.groupTS.F,
#                          levels = colnames(coef(fit)))
# 
# # Estimate contrast for each gene
# tmp.TS <- contrasts.fit(fit, 
#                         contr.TS)
# 
# # Empirical Bayes smoothing 
# tmp.TS <- eBayes(tmp.TS)
# 
# # What genes are most differentially expressed?
# top.table.TS <- topTable(tmp.TS, sort.by = "P", n = Inf)
# 
# # number of DE genes
# length(which(top.table.TS$adj.P.Val < 0.05))
# 
# # graph DEGs
# # volcano plots
# # TS
# top.table.TS %>% 
#   mutate(color = ifelse(adj.P.Val < 0.05,
#                         'sig',
#                         'not_sig')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log(adj.P.Val),
#              color = color)) +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlim(-9,9) +
#   ggtitle('TS male vs female') 
# ggsave('DEG/figures/BBvsTS/TS M vs F volcano.png')
# 
# #### interaction
# ## specify model to be fitted
# # interaction of sex and species
# mm.int <- model.matrix(~sample.species*sample.sex)
# 
# ### voom transformation
# png("DEG/figures/BBvsTS/voom BB vs TS interaction.png",
#     height = 10,
#     width = 10,
#     unit = 'in',
#     res = 480)
# yint <- voom(d, 
#           mm.int, 
#           plot = T)
# dev.off()
# 
# ## fitting ilnear models in limma
# fit.int = lmFit(yint,
#             mm.int)
# 
# ### make contrasts
# ## species
# tmp.species = contrasts.fit(fit.int,
#                               coef = 2)
# 
# # Empirical Bayes smoothing 
# tmp.species <- eBayes(tmp.species)
# 
# # What genes are most differentially expressed?
# top.table.species <- topTable(tmp.species, sort.by = "P", n = Inf)
# 
# # number of DE genes
# length(which(top.table.species$adj.P.Val < 0.05))
# 
# # graph DEGs
# # volcano plots
# # set x limits
# x.species = top.table.species %>% 
#   pull(logFC) %>% 
#   abs() %>% 
#   max()
# 
# # species
# top.table.species %>% 
#   mutate(color = ifelse(adj.P.Val < 0.05,
#                         'sig',
#                         'not_sig')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log(adj.P.Val),
#              color = color)) +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlim(-x.species,x.species) +
#   ggtitle('Species BB vs TS') 
# ggsave('DEG/figures/BBvsTS/species BB vs TS volcano.png')
# 
# ## sex
# tmp.sex = contrasts.fit(fit.int,
#                             coef = 3)
# 
# # Empirical Bayes smoothing 
# tmp.sex <- eBayes(tmp.sex)
# 
# # What genes are most differentially expressed?
# top.table.sex <- topTable(tmp.sex, sort.by = "P", n = Inf)
# 
# # number of DE genes
# length(which(top.table.sex$adj.P.Val < 0.05))
# 
# # graph DEGs
# # volcano plots
# # set x limits
# x.sex = top.table.sex %>% 
#   pull(logFC) %>% 
#   abs() %>% 
#   max()
# 
# # sex
# top.table.sex %>% 
#   mutate(color = ifelse(adj.P.Val < 0.05,
#                         'sig',
#                         'not_sig')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log(adj.P.Val),
#              color = color)) +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlim(-x.sex,x.sex) +
#   ggtitle('Sex BB vs TS') 
# ggsave('DEG/figures/BBvsTS/Sex M vs F volcano.png')
# 
# ## interaction
# tmp.interaction = contrasts.fit(fit.int,
#                         coef = 4)
# 
# # Empirical Bayes smoothing 
# tmp.interaction <- eBayes(tmp.interaction)
# 
# # What genes are most differentially expressed?
# top.table.interaction <- topTable(tmp.interaction, sort.by = "P", n = Inf)
# 
# # number of DE genes
# length(which(top.table.interaction$adj.P.Val < 0.05))
# 
# # graph DEGs
# # volcano plots
# # set x limits
# x.interaction = top.table.interaction %>% 
#   pull(logFC) %>% 
#   abs() %>% 
#   max()
# 
# # interaction
# top.table.interaction %>% 
#   mutate(color = ifelse(adj.P.Val < 0.05,
#                         'sig',
#                         'not_sig')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log(adj.P.Val),
#              color = color)) +
#   geom_hline(yintercept = -log(0.05)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlim(-x.interaction,x.interaction) +
#   ggtitle('interaction species and sex') 
# ggsave('DEG/figures/BBvsTS/interaction species and sex.png')
# 
# #### RRHO2 BB vs TS ####
# ## compare concordance between bluebirds and tree swallow male and female analysis
# ## use updated RedRibbon
# # devtools::install_github("antpiron/RedRibbon")
# library(RedRibbon)
# 
# ## Create RedRibbon object
# # use log fold change for BB and TS
# df.rr = top.table.BB %>% 
#   select(logFC) %>% 
#   rownames_to_column('id') %>% 
#   dplyr::rename(a = logFC) %>% 
#   full_join( top.table.TS %>% 
#                select(logFC) %>% 
#                rownames_to_column('id') %>% 
#                dplyr::rename(b = logFC))
# 
# # create RR object
# rr <- RedRibbon(df.rr, 
#                 enrichment_mode="hyper-two-tailed")
# 
# ## Run the overlap using evolutionary algorithm,
# ## computing permutation adjusted p-value for the four quadrants
# quad <- quadrants(rr, 
#                   algorithm="ea", 
#                   permutation=TRUE, 
#                   whole=FALSE)
# 
# ## Plots the results
# ggRedRibbon(rr, 
#             quadrants=quad) + 
#   coord_fixed(ratio = 1, 
#               clip = "off")
# ggsave('DEG/figures/BBvsTS/RRHO BB vs TS across sex.png')
# 
# ## Get the list of genes for each quadrant
# df.rr.results = df.rr[quad$upup$positions,] %>% 
#   mutate(position = 'uu') %>% 
#   rbind(df.rr[quad$downdown$positions,] %>% 
#           mutate(position = 'dd')) %>% 
#   rbind(df.rr[quad$updown$positions,] %>% 
#           mutate(position = 'ud')) %>% 
#   rbind(df.rr[quad$downup$positions,] %>% 
#           mutate(position = 'du')) 
# 
# # check if these are Z genes?
# # load gene chromsome position
# data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')
# 
# # combine with up-up genes
# df.rr.results = df.rr.results %>% 
#   left_join(data.gene.chromosome %>% 
#               select(c(Symbol,
#                        chromosome_name)) %>% 
#             dplyr::rename(id = Symbol))
# 
# # graph
# df.rr.results %>% 
#   group_by(chromosome_name,
#            position) %>% 
#   summarise(count = n()) %>% 
#   mutate(chromosome_name = ifelse(is.na(chromosome_name),
#                                   'A',
#                                   chromosome_name)) %>% 
#   ggplot(aes(x = reorder(chromosome_name,
#                      -count),
#          y = count,
#          color = position)) +
#   geom_point() +
#   theme_classic()
# ggsave('DEG/figures/BBvsTS/RRHO chromosomes BB vs TS across sex.png')
# 
# # graph
# # filter
# df.rr.results %>% 
#   group_by(chromosome_name,
#            position) %>% 
#   summarise(count = n()) %>% 
#   mutate(chromosome_name = ifelse(is.na(chromosome_name),
#                                   'A',
#                                   chromosome_name)) %>% 
#   filter(chromosome_name != 'A') %>% 
#   filter(count > 20) %>% 
#   ggplot(aes(x = reorder(chromosome_name,
#                          -count),
#              y = count,
#              color = position)) +
#   geom_point() +
#   theme_classic()
# ggsave('DEG/figures/BBvsTS/RRHO chromosomes BB vs TS across sex reduce.png')
# 
# 
# 
# 
# 
#   
#### all species sex DEG analysis ####
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

#### pre-processing
## create counts data frame
counts = data.long %>% 
  select(-c(sample_name,
            species,
            sex)) %>% 
  column_to_rownames('sample_id') %>% 
  t()

## create DGElist object
d0 = DGEList(counts)

## calculate nomalization factors
d0 = calcNormFactors(d0)

## filter low-expressed genes
# set cut off
cutoff <- 5
# which genes to drop
drop <- which(apply(cpm(d0), 1, max) < cutoff)
# remove low expressed genes
d <- d0[-drop,] 
# number of genes left
dim(d) 

## create sample information
sample.df = data.frame(sample_id = colnames(counts))
# add sample data on species and sex
sample.df = sample.df %>% 
  left_join(sample_ids )
# create list of sex
sample.sex = sample.df %>% 
  pull(sex)
# create list of species 
sample.sex = sample.df %>% 
  pull(sex)
# create list of species
sample.species = sample.df %>% 
  pull(species)
# create list of species
sample.group = sample.df %>% 
  mutate(group = paste(species,
                       sex,
                       sep = '.')) %>% 
  pull(group)


### graph MDS
png("DEG/figures/All_sex/MDS all_sex.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
plotMDS(d, 
        col = as.numeric(as.factor(sample.group)))
dev.off()



#### voom transformation and calculation of variance weight
#### compare each species by sex
## specify model to be fitted
mm <- model.matrix(~0+sample.group)

### voom transformation
png("DEG/figures/All_sex//voom all_sex.png",
    height = 10,
    width = 10,
    unit = 'in',
    res = 480)
y <- voom(d, 
          mm, 
          plot = T)
dev.off()

## fitting linear models in limma
fit = lmFit(y,
            mm)

### make contrasts
## all
contr.all = makeContrasts(sample.groupBB.M - sample.groupBB.F,
                          sample.groupBS.M - sample.groupBS.F,
                          sample.groupCW.M - sample.groupCW.F,
                          sample.groupET.M - sample.groupET.F,
                          sample.groupHS.M - sample.groupHS.F,
                          sample.groupHW.M - sample.groupHW.F,
                          sample.groupPW.M - sample.groupPW.F,
                          sample.groupRO.M - sample.groupRO.F,
                          sample.groupTS.M - sample.groupTS.F,
                          sample.groupYW.M - sample.groupYW.F,
                         levels = colnames(coef(fit)))

# Estimate contrast for each gene
tmp.cont <- contrasts.fit(fit, 
                          contr.all)

# Empirical Bayes smoothing 
tmp.cont <- eBayes(tmp.cont)

### What genes are most differential expressed?
## create top table with every contrast
# create empty dataframe
top.table = data.frame(GeneName = as.character(),
                       logFC = as.numeric(),
                       AveExpr = as.numeric(),
                       t = as.numeric(),
                       P.Value = as.numeric(),
                       adj.P.Val = as.numeric(),
                       B = as.numeric(),
                       contrast = as.character())

# for each of 10 species and contrasts
for (i in colnames(tmp.cont$coefficients)) {
  # get information for each contrast
top.table.tmp <- topTable(tmp.cont,
                         sort.by = "P",
                      coef = i,
                         n = Inf) %>% 
  rownames_to_column('GeneName') %>% 
  mutate(contrast = i)
# combine with top.table
top.table = top.table %>% 
  rbind(top.table.tmp)
}

# convert contrast names to species
top.table = top.table %>%
  separate_wider_position(contrast,
                          c(drop = 12,
                            species = 2),
                          too_many = 'drop',
                          cols_remove = F) %>% 
  select(-c(drop))

### save results 
write.csv(top.table,
          file = 'DEG/top.table.all.species.csv')

# load results
top.table = read.csv('DEG/top.table.all.species.csv')

#### graph DEGs all species ####
## stats
top.table %>% 
  filter(adj.P.Val < 0.05) %>% 
  mutate(direction = ifelse(logFC > 0,
                            'male_bias',
                            'female_bias')) %>% 
  select(species,
         direction) %>% 
  table() %>% 
  as.data.frame() %>% 
  group_by(direction) %>% 
  summarise(mean = mean(Freq),
            sd = sd(Freq))
  
  
  
# number of DE genes
# top.table %>% 
#   filter(adj.P.Val < 0.05) %>% 
#   mutate(direction = ifelse(logFC > 0,
#                             'male_bias',
#                             'female_bias')) %>% 
#   select(species,
#          direction) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Freq = ifelse(direction == 'female_bias',
#                        -Freq,
#                        Freq)) %>% 
#   ggplot(aes(y = reorder(species,
#                          -Freq),
#              x = Freq,
#              fill = direction)) +
#   geom_vline(xintercept = 0) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   scale_fill_manual(values = c('black',
#                                'red')) +
#   ylab('Species') +
#   xlab('Number of DEG')
# ggsave('DEG/figures/All_sex/Number of DEG per species.png')

## graph DEGs for all species
# volcano plots
# loop across all species 
for (i in unique(top.table$species)) {
top.table %>% 
    filter(species == i) %>% 
  mutate(color = ifelse(adj.P.Val < 0.05,
                        'sig',
                        'not_sig')) %>% 
  ggplot(aes(x = logFC,
             y = -log(adj.P.Val),
             color = color)) +
  geom_hline(yintercept = -log(0.05)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle(paste0(i,
                 ' male vs female'))
ggsave(paste0('DEG/figures/All_sex/volcano_plots/',
              i,
              ' M vs F volcano.png'))

}


### Compare DEG with Z chromosome
## load gene chromosome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

## add Z to all species DEG
# paper
top.table %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Symbol,
                     chromosome_name) %>% 
              dplyr::rename(GeneName = Symbol) %>% 
              dplyr::rename(chr = chromosome_name) %>% 
              mutate(chr = ifelse(chr == 'Z',
                                  'Z',
                                  'A'),
                     chr = ifelse(is.na(chr),
                                  'A',
                                  chr))) %>% 
  filter(adj.P.Val < 0.05) %>% 
  mutate(direction = ifelse(logFC > 0,
                            'male_bias',
                            'female_bias')) %>% 
  select(species,
         direction,
         chr) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(Freq = ifelse(direction == 'female_bias',
                       -Freq,
                       Freq)) %>% 
  left_join(data.frame(species = c("ET",
                                   "HS",
                                   "PW",
                                   "YW",
                                   "BB",
                                   "RO",
                                   "HW",
                                   "CW",
                                   "TS",
                                   "BS"),
                       order = c(1:10))) %>% 
  ggplot(aes(y = reorder(species,
                         -order),
             x = Freq)) +
  geom_bar(aes(group = chr,
               fill = chr),
           stat = 'identity',
           # fill = 'lightgrey',
           # color = 'black'
           ) +
  geom_vline(xintercept = 0) +
  theme_classic(base_size = 8) +
  ylab('Species') +
  xlab('Number of DEGs') +
  scale_fill_manual(values = c('lightgrey',
                               'red'))+ 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.9, 
                            0.5))
ggsave('DEG/figures/All_sex/Number of DEG per species.pdf',
       height = 3.9,
       width = 3.1,
       units = 'in',
       dpi = 320)

## stats
top.table %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Symbol,
                     chromosome_name) %>% 
              dplyr::rename(GeneName = Symbol) %>% 
              dplyr::rename(chr = chromosome_name) %>% 
              mutate(chr = ifelse(chr == 'Z',
                                  'Z',
                                  'A'),
                     chr = ifelse(is.na(chr),
                                  'A',
                                  chr))) %>% 
  filter(adj.P.Val < 0.05) %>% 
  mutate(direction = ifelse(logFC > 0,
                            'male_bias',
                            'female_bias')) %>% 
  select(species,
         direction,
         chr) %>% 
  table() %>% 
  as.data.frame() %>% 
  group_by(direction,
           chr) %>% 
  summarise(mean = mean(Freq),
            sd = sd(Freq))

### graph upset R
## create dataframe with variables:
# Gene name, sig_TS (1 or 0), sig_BB (1 or 0), etc.., Z chromosome (Z or A)
top.table.upset =
  top.table %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Symbol,
                     chromosome_name) %>% 
              dplyr::rename(GeneName = Symbol) %>% 
              dplyr::rename(chr = chromosome_name) %>% 
              mutate(chr = ifelse(chr == 'Z',
                                  'Z',
                                  'A'),
                     chr = ifelse(is.na(chr),
                                  'A',
                                  chr))) %>% 
  mutate(sig = ifelse(adj.P.Val < 0.05,
                      1,
                      0)) %>% 
    pivot_wider(id_cols = c('GeneName',
                            'chr'),
                names_from = species,
                values_from = sig) %>% 
  na.omit() %>% 
  filter_at(vars(-GeneName,
                 -chr), 
            any_vars(. != 0)) # remove rows with all zeros

## make one for male biased DEGs and female biased DEGs
# male
top.table.upset.m =
  top.table %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Symbol,
                     chromosome_name) %>% 
              dplyr::rename(GeneName = Symbol) %>% 
              dplyr::rename(chr = chromosome_name) %>% 
              mutate(chr = ifelse(chr == 'Z',
                                  'Z',
                                  'A'),
                     chr = ifelse(is.na(chr),
                                  'A',
                                  chr))) %>% 
  mutate(sig = ifelse(adj.P.Val < 0.05,
                      1,
                      0),
         sig.male = ifelse(logFC>0,
                                sig,
                                0)) %>% 
  pivot_wider(id_cols = c('GeneName',
                          'chr'),
              names_from = species,
              values_from = sig.male) %>% 
  na.omit() %>% 
  filter_at(vars(-GeneName,
                 -chr), 
            any_vars(. != 0)) # remove rows with all zeros

# male
top.table.upset.f =
  top.table %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Symbol,
                     chromosome_name) %>% 
              dplyr::rename(GeneName = Symbol) %>% 
              dplyr::rename(chr = chromosome_name) %>% 
              mutate(chr = ifelse(chr == 'Z',
                                  'Z',
                                  'A'),
                     chr = ifelse(is.na(chr),
                                  'A',
                                  chr))) %>% 
  mutate(sig = ifelse(adj.P.Val < 0.05,
                      1,
                      0),
         sig.female = ifelse(logFC<0,
                           sig,
                           0)) %>% 
  pivot_wider(id_cols = c('GeneName',
                          'chr'),
              names_from = species,
              values_from = sig.female) %>% 
  na.omit() %>% 
  filter_at(vars(-GeneName,
                 -chr), 
            any_vars(. != 0)) # remove rows with all zeros

# species
species.list = top.table %>% 
  pull(species) %>% 
  unique()

## graph upset plot
# png('DEG/figures/All_sex/upset/All DEGs upset vs Z.png')
# upset(top.table.upset,
#       species.list,
#       min_size = 10,
#       base_annotations=list(
#         'All DEGs'=intersection_size(
#           # counts=FALSE,
#           mapping=aes(fill = chr)
#         ) + scale_fill_manual(values = c(
#           'A'='black',
#           'Z'='red'
#         ))
#       ),
#       width_ratio=0.1,
#       wrap = T,
#       set_sizes = (
#         upset_set_size(
#           geom = geom_bar(
#             aes(fill = chr)) 
#         )+ scale_fill_manual(values = c(
#           'A'='black',
#           'Z'='red'
#         ),
#         guide = 'none')
#       ))
# dev.off()
# 
# 
# # male only DEGs
# png('DEG/figures/All_sex/upset/Male DEGs upset vs Z.png')
# upset(top.table.upset.m,
#       species.list,
#       min_size = 5,
#       base_annotations=list(
#         'Male DEGs'=intersection_size(
#           # counts=FALSE,
#           mapping=aes(fill = chr)
#         ) + scale_fill_manual(values = c(
#           'A'='black',
#           'Z'='red'
#         ))
#       ),
#       width_ratio=0.1,
#       wrap = T,
#       set_sizes = (
#         upset_set_size(
#           geom = geom_bar(
#             aes(fill = chr)) 
#         )+ scale_fill_manual(values = c(
#           'A'='black',
#           'Z'='red'
#         ),
#         guide = 'none')
#       ))
# dev.off()

# males
# paper
png('DEG/figures/All_sex/upset/Male DEGs upset vs Z.png',
    height = 7,
    width = 10.85,
    units = 'in',
    res = 320,
    pointsize = 24)
upset(top.table.upset.m,
      species.list,
      min_size = 5,
      base_annotations=list(
        'Male DEGs'=intersection_size(
          # counts=FALSE,
          mapping=aes(fill = chr)
        ) + scale_fill_manual(values = c(
          'A'='black',
          'Z'='red'
        ), guide = 'none')
      ),
      width_ratio=0.1,
      wrap = T,
      set_sizes = (
        upset_set_size(
          geom = geom_bar(
            aes(fill = chr)) 
        )+ scale_fill_manual(values = c(
          'A'='black',
          'Z'='red'
        ),
        guide = 'none') 
      ))
dev.off()

# female only DEGs
png('DEG/figures/All_sex/upset/Female DEGs upset vs Z.png',
    height = 7,
    width = 10.85,
    units = 'in',
    res = 320)
upset(top.table.upset.f,
      species.list,
      min_size = 4,
      base_annotations=list(
        'Female DEGs'=intersection_size(
          # counts=FALSE,
          mapping=aes(fill = chr)
        ) + scale_fill_manual(values = c(
          'A'='black',
          'Z'='red'
        ), guide = 'none')
      ),
      width_ratio=0.1,
      wrap = T,
      set_sizes = (
        upset_set_size(
          geom = geom_bar(
            aes(fill = chr)) 
        )+ scale_fill_manual(values = c(
          'A'='black',
          'Z'='red'
        ),
        guide = 'none')
      ))
dev.off()





# #### RRHO2 ####
# ## compare concordance between all pairs of clades
# # use DEG of cavity nesting vs not-cavity nesting species 
# 
# # create empty dataframe
# df.rr.results = data.frame(species.a = as.character(),
#                            species.b = as.character(),
#                            comp = as.numeric(),
#                            dist = as.numeric(),
#                            uu = as.numeric(),
#                            dd = as.numeric(),
#                            up = as.numeric(),
#                            du = as.numeric())
# 
# 
# ## loop through every species pair
# for (i in data.species.pair.dist$comp) {
#   # get first species for comparison
#   species.a = data.species.pair.dist %>% 
#     filter(comp == i) %>% 
#     pull(V1)
#   # get second species for comparison
#   species.b = data.species.pair.dist %>% 
#     filter(comp == i) %>% 
#     pull(V2)
#   
#   ## Create RedRibbon object
#   # use log fold change for two speciess
#   # species a
#   df.rr.a = top.table %>% 
#     filter(species == species.a) %>% 
#     select(GeneName,
#            logFC) %>% 
#     rename(id = GeneName) %>% 
#     rename(a = logFC)
#   
#   # species b
#   df.rr.b = top.table %>% 
#     filter(species == species.b) %>% 
#     select(GeneName,
#            logFC) %>% 
#     rename(id = GeneName) %>% 
#     rename(b = logFC)
#   
#   # combine 
#   df.rr = df.rr.a %>% 
#     full_join(df.rr.b)
#   
#   ## create RR object
#   rr <- RedRibbon(df.rr, 
#                   enrichment_mode="hyper-two-tailed")
#   
#   
#   ## Run the overlap using evolutionary algorithm,
#   ## computing permutation adjusted p-value for the four quadrants
#   quad <- quadrants(rr, 
#                     algorithm="ea", 
#                     permutation=TRUE, 
#                     whole=FALSE)
#   
#   
#   ## Plots the results
#   ggRedRibbon(rr, 
#               quadrants=quad) + 
#     coord_fixed(ratio = 1, 
#                 clip = "off")
#   ggsave(paste0('DEG/figures/All_sex/RRHO/heatmaps/RRHO ',
#                 species.a,
#                 ' vs ',
#                 species.b,
#                 '.png'))
#   
#   
#   ## Get the list of genes for each quadrant
#   tmp.results = data.frame(species.a = species.a,
#                            species.b = species.b,
#                            comp = i,
#                            dist = data.species.pair.dist %>% 
#                              filter(comp == i) %>% 
#                              pull(dist),
#                            uu = quad$upup$positions %>%
#                              length(),
#                            dd = quad$downdown$positions %>%
#                              length(),
#                            ud = quad$updown$positions %>%
#                              length(),
#                            du = quad$downup$positions %>%
#                              length())
#   
#   # combine with results 
#   df.rr.results = df.rr.results %>% 
#     rbind(tmp.results)
#   
# }
# 
# ## save RRHO comparison
# write.csv(df.rr.results,
#           file = 'DEG/df.rr.results.all.species.csv')
# 
# ### compare distance with RRHO results
# # concordance
# df.rr.results %>% 
#   mutate(concordance = uu + dd) %>% 
#   ggplot(aes(x = dist,
#              y = concordance)) +
#   geom_point() +
#   # geom_smooth(method = 'lm') +
#   theme_classic()
# ggsave('DEG/figures/All_sex/RRHO/Distance vs concordance.png')
# 
# # disconcordance
# df.rr.results %>% 
#   mutate(discordance = ud + du,
#          discordance = -discordance) %>% 
#   ggplot(aes(x = dist,
#              y = discordance)) +
#   geom_point() +
#   # geom_smooth(method = 'lm') +
#   theme_classic()
# ggsave('DEG/figures/All_sex/RRHO/Distance vs discordance.png')
# 
# # ### if wanted to pull gene names from RedRibbon
# # # convert quad output to dataframe
# # # lists information for quadrant but also gene position in list
# # df = quad %>% 
# #   plyr::ldply(data.frame)
# # 
# # # add gene name and logfc in species A and B
# # df = df %>% 
# #   left_join(df.rr %>% 
# #               rownames_to_column('positions') %>% 
# #               mutate(positions = as.numeric(positions)))
# 
# 

#### format data for paper supplemental ####
### take sex difference DEG output
## add chromosomal data
# load DEG results
top.table = read.csv('DEG/top.table.all.species.csv') %>% 
  dplyr::select(-c(X))

# load chromosome data
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

# combine DEG results with chromosome data
top.table.chr = top.table %>% 
  left_join(data.gene.chromosome %>% 
  dplyr::select(Symbol,
                chromosome_name) %>% 
  dplyr::rename(GeneName = Symbol) %>% 
    mutate(chromosome_name = ifelse(is.na(chromosome_name),
                                    'A',
                                    chromosome_name)))

## save file
write.csv(top.table.chr,
          'DEG/Dataset_S2.csv',
          row.names = F)

