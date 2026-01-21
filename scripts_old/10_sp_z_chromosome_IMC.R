#### 10 sp bird analysis 
### Z chromosome 
# R 4.2.1

## help making trees: yulab-smu.top/treedata-book/

## setwd
setwd("/geode2/home/u040/imillerc/Quartz/10_sp")

#### load libraries ####
## install changepoint.np
# install.packages('changepoint.np')

# load library
library(changepoint.np)
library(gtools)
library(rstatix)
library(boot)
library(ape)
library(tidytree)
library(ggtree)
library(tidyverse)

#### load data ####
### load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

### load ratio data
data.wide.format.z = read.delim('CAGEE/data/normalizedCounts_ratio_z.tsv')

## combine ratio data with chromosome data
data.wide.format.chr = data.wide.format.z %>% 
  left_join(data.gene.chromosome,
            by = c('GeneName' = 'Symbol'))

## load gametolog data
data.gametolog = read_tsv('Gene_list/10sp_ensembl_gametologs_edit.tsv') %>% 
  select(-c(Ensembl.GeneIDs)) %>% 
  rename(GeneName = Symbol) %>% 
  distinct()

# # check chromosome position of gametologs
data.wide.format.chr %>%
  inner_join(data.gametolog ) %>% 
  pull(SAMPLETYPE) 
# most on Z 
# 3 have no chromosome assigned

# count
data.wide.format.chr %>%
  inner_join(data.gametolog ) %>% 
  pull(chromosome_name) %>% 
  length()
# 35

## save gametolog for CAGEE
# change SAMPLETYPE to compare gametolog
# make sure to add in genes that were unassigned to Z
data.wide.format.chr.gametolog = data.wide.format.chr %>% 
  left_join(data.gametolog) %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'Z',
                             SAMPLETYPE)) %>% 
  filter(SAMPLETYPE == 'Z') %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'ZG',
                             SAMPLETYPE))

# remove extra columns
data.wide.format.chr.gametolog = data.wide.format.chr.gametolog[,-c(14:23)]

# save file
write_delim(data.wide.format.chr.gametolog,
            'CAGEE/data/normalizedCounts_ratio_gametolog.tsv',
            delim = '\t')

## median gene expression for gametologs
### load ratio data
data.wide.format.median.z = read.delim('CAGEE/data/normalizedCounts_median_z.tsv')

# combine median data with chromosome data
data.wide.format.median.z.chr = data.wide.format.median.z %>% 
  left_join(data.gene.chromosome,
            by = c('GeneName' = 'Symbol'))

# combine with gametolog data
data.wide.format.median.chr.gametolog = data.wide.format.median.z.chr %>% 
  left_join(data.gametolog) %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'Z',
                             SAMPLETYPE)) %>% 
  filter(SAMPLETYPE == 'Z') %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'ZG',
                             SAMPLETYPE))

# remove extra columns
data.wide.format.median.chr.gametolog = data.wide.format.median.chr.gametolog[,-c(14:23)]

# save file
write_delim(data.wide.format.median.chr.gametolog,
            'CAGEE/data/normalizedCounts_median_gametolog.tsv',
            delim = '\t')

### format z chromsome dataframe
## reduce columns, select z chromsome
data.z = data.wide.format.chr %>% 
  dplyr::select(-c(GeneDescription,
                   SAMPLETYPE,
                   NCBI.GeneID,
                   Description,
                   Transcripts,
                   Ensembl.GeneIDs,
                   external_gene_name)) %>% 
  filter(chromosome_name == 'Z')

## pivot longer  
data.z.long = data.z %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression')

#### create gametolog upsetR across birds ####
### load in bird gametologs
data.all.bird.game = read.csv('Gene_list/Gametologs_Xu_etal_NatEcolEvol_2019.csv')

## create presence/absence matrix for each gene
data.all.bird.game = data.all.bird.game %>% 
  select(-c(order,
            Strata)) %>% 
  mutate(across(everything(), 
                as.character)) %>% 
  pivot_longer(cols = house.sparrow:chicken,
               names_to = 'species',
               values_to = 'present') %>% 
  mutate(present = ifelse(is.na(present),
                          0,
                          1)) %>% 
  pivot_wider(names_from = species,
              values_from = present) %>% 
  column_to_rownames('gene')
  
### create Upset plot
library(UpSetR)

## graph all birds and all intersects
png('z_chromosome/figures/gametolog/upset all birds.png',
    height = 5,
    width = 10,
    units = 'in',
    res = 480)
upset(data.all.bird.game,
      nsets = ncol(data.all.bird.game),
      order.by = 'freq',
      decreasing = T,
      nintersects = NA,
      sets = colnames(data.all.bird.game),
      sets.bar.color = c(rep('black', 
                             8),
                         'red',
                         rep('black',
                             3)))
dev.off()

## graph subset of birds and all intersects
png('z_chromosome/figures/gametolog/upset song birds.png',
    height = 5,
    width = 10,
    units = 'in',
    res = 480)
upset(data.all.bird.game,
      nsets = ncol(data.all.bird.game),
      order.by = 'freq',
      decreasing = T,
      nintersects = NA,
      sets = c("house.sparrow",
               "common.canary",
               "medium.ground.finch",
               "collared.flycatcher",
               "ground.tit"),
      sets.bar.color = c('red',
                         rep('black',
                             4)))
dev.off()

#### calculate stats Z chromosome ####
#### find change points in Z
### create wide format z chr data
# each row is a species
# 2 genes have same start position: CDC14B & PRXL2C
# remove CDC14B
data.z.long.wide = data.z.long %>% 
  filter(GeneName != 'CDC14B') %>% 
  arrange(start_position) %>% 
  pivot_wider(id_cols = species,
              names_from = start_position,
              values_from = gene_expression) %>% 
  column_to_rownames('species') %>% 
  as.matrix()


### run across each species
data.change.point = cpt.np(data.z.long.wide)
 
## create dataframe of results
# create empty dataframe
data.change.point.results = NULL

# loop through species
for (i in 1:nrow(data.z.long.wide)) {
  # create species temporary dataframe
tmp = data.frame(species = rownames(data.z.long.wide)[i],
                                       change.points = data.change.point[[i]]@cpts)
# combine with other species
data.change.point.results = data.change.point.results %>% 
  rbind(tmp)
}


## combine with position and gene name
data.change.point.results = data.change.point.results %>% 
  left_join(data.z.long.wide %>% 
              colnames() %>% 
              as.data.frame() %>% 
              rename(position = '.') %>% 
              rownames_to_column('change.points') %>% 
              mutate(position = as.integer(position),
                     change.points = as.integer(change.points)))

### graph change points
data.change.point.results %>% 
  ggplot(aes(x = position/1000000,
             y = species)) +
  geom_rect(aes(xmin = 25,
                xmax = 35,
                ymin = 0,
                ymax = 11,
                fill = 'red')) +
  geom_rect(aes(xmin = 72.5,
                xmax = 73.5,
                ymin = 0,
                ymax = 11,
                fill = 'red')) +
  geom_point() +
  theme_classic() +
  xlim(0,80)+
  xlab('Z chr position') +
  ggtitle('Change points Z chromosome') +
  theme(legend.position = "none")
ggsave('z_chromosome/figures/Change points Z chr.png')

### calculate rolling average
# moving window average across 30 genes with shift of 1 gene
data.z.long.run = data.z.long %>% 
  arrange(start_position) %>% 
  group_by(species) %>% 
  mutate(gene_expression_running = running(gene_expression,
                                           width = 30,
                                           fun = mean,
                                           align = 'left',
                                           allow.fewer = TRUE,
                                           by = 1)) 


## graph running average
# loop through all species

for (i in unique(data.z.long.run$species)) {

# get box heights set to max and min value
# max
max.v = data.z.long.run %>% 
  filter(species == i) %>% 
  pull(gene_expression) %>% 
  max() %>% 
  log()
# min
min.v = data.z.long.run %>% 
  filter(species == i) %>% 
  pull(gene_expression) %>% 
  min() %>% 
  log()

# add change points and MHM windows
data.z.long.run %>% 
  filter(species == i) %>% 
  ggplot() +
  annotate("rect",
           xmin = 25,
                xmax = 35,
                ymin = min.v,
                ymax = max.v,
                fill = 'red') +
  annotate("rect",
          xmin = 72.5,
                xmax = 73.5,
                ymin = min.v,
                ymax = max.v,
                fill = 'red') +
  geom_point(aes(x = start_position/1000000,
                 y = log(gene_expression))) +
  geom_line(aes(x = start_position/1000000,
                y = log(gene_expression_running)),
            linewidth = 2) +
  geom_vline(xintercept = data.change.point.results %>%
               filter(species == i) %>% 
               pull(position)/1000000,
             linetype = 'dashed') +
  theme_classic() +
  xlim(0,80)+
  xlab('Z chr position') +
  ylab('Male:Female ratio') +
  ggtitle(i) + 
  theme(legend.position = "none")
ggsave(paste0('z_chromosome/figures/z_chr_windows/',
              i, 
              ' window Z chr.png'))
}


#### calculate median for each species
data.z.median.MHM1 = data.z.long %>% 
  mutate(MHM1 = ifelse(25 < start_position/1000000 & 35 > start_position/1000000,
                       "MHM1",
                       "other"),
         log_gene_expression = log(gene_expression)) 

## loop through each species
# create empty dataframe
df.MHM1.median = NULL
# run loop
for (i in unique(data.z.long.run$species)) {
  # get list of MHM1 genes
tmp.1 = data.z.median.MHM1 %>% 
  filter(species == i) %>% 
  filter(MHM1 == "MHM1") %>% 
  pull(log_gene_expression)

# get list of other genes
tmp.2 = data.z.median.MHM1 %>% 
  filter(species == i) %>% 
  filter(MHM1 != "MHM1") %>% 
  pull(log_gene_expression)

## MHM1 genes
# set up confidence interval
b.1 <- boot(tmp.1,
          statistic = function(tmp.1,i) median(tmp.1[i]),
          R = 1000)
# get confidence interval
b.1.ci = boot.ci(b.1)
# calculate median
m.1 = median(tmp.1)

## Other Z genes
# set up confidence interval
b.2 <- boot(tmp.2,
                  statistic = function(tmp.2,i) median(tmp.2[i]),
                  R = 1000)
# get confidence interval
b.2.ci = boot.ci(b.2)
# calculate median
m.2 = median(tmp.2)

# run wilcox test
tmp = data.z.median.MHM1 %>% 
  filter(species == i) %>% 
  wilcox_test(log_gene_expression ~ MHM1,
              conf.level = 0.95,
              detailed = T) %>% 
  mutate(species = i,
         MHM1.median = m.1,
         MHM1.ci.low = b.1.ci[["normal"]][2],
         MHM1.ci.high = b.1.ci[["normal"]][3],
         Other.median = m.2,
         Other.ci.low = b.2.ci[["normal"]][2],
         Other.ci.high = b.2.ci[["normal"]][3])

# add to dataframe
df.MHM1.median = df.MHM1.median %>% 
  rbind(tmp)

}

## save statistics 
write.csv(df.MHM1.median,
          'z_chromosome/figures/z_chr_windows/df.MHM1.median.csv')
# # load
# df.MHM1.median = read.csv('z_chromosome/figures/z_chr_windows/df.MHM1.median.csv')

# graph median 
df.MHM1.median %>% 
  select(species,
         MHM1.median,
         Other.median) %>% 
  pivot_longer(cols = -c(species),
               names_to = 'type',
               values_to = 'median') %>% 
  ggplot(aes(x = type,
             y = median)) +
  geom_boxplot() +
  geom_line(aes(group = species,
                color = species)) +
  geom_point(aes(color = species)) +
  theme_classic()
ggsave('z_chromosome/figures/z_chr_windows/MHM1 vs Other median.png')


#### graph Z chromosome gametolog sex ratio distribution ####
### format dataframe for Z
## reduce columns, select z chromsome
data.z = data.wide.format.chr %>% 
  dplyr::select(-c(GeneDescription,
                   SAMPLETYPE,
                   NCBI.GeneID,
                   Description,
                   Transcripts,
                   external_gene_name)) %>% 
  left_join(data.gametolog)  %>% 
  mutate(chromosome_name = ifelse(!is.na(gametolog),
                             'Z',
                             chromosome_name)) %>% 
  filter(chromosome_name == 'Z') 

## pivot longer  
data.z.long = data.z %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression')

### format dataframe for all
## reduce columns, select z chromsome
data.long = data.wide.format.chr %>% 
  dplyr::select(-c(GeneDescription,
                   SAMPLETYPE,
                   NCBI.GeneID,
                   Description,
                   Transcripts,
                   external_gene_name))%>% 
  mutate(chromosome_type = case_when(chromosome_name == 'Z' ~ 'Z',
                                  TRUE ~ 'Autosome')) %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression') 

# ### graph ratio of Z for each species
# # boxplot
# data.z.long %>% 
#   ggplot(aes(y = log(gene_expression),
#              x = species,
#              fill = as.factor(gametolog))) +
#   geom_boxplot(color = 'black') +
#   theme_classic() +
#   ylim(-1,1)+
#   ylab('log Male:Female gene expression ratio') +
#   scale_fill_manual(values = c('#9A1B1F')) +
#   labs(fill = 'Gametolog')
# ggsave('z_chromosome/figures/gametolog/gametolog species boxplot.png',
#        height = 10,
#        width = 10)
# 
# # ridge plot
# data.z.long%>% 
#   ggplot(aes(x = log(gene_expression),
#              y = species,
#              fill = as.factor(gametolog))) +
#   geom_vline(xintercept = 0.69,
#              linetype = 'dashed')+
#   geom_vline(xintercept = -0.69,
#              linetype = 'dashed')+
#   ggridges::geom_density_ridges() +
#   theme_classic() + 
#   theme(text = element_text(size = 24),
#         axis.text.y =  element_text(size = 16),
#         axis.text.x =  element_text(size = 16),
#         legend.position = 'null') +
#   xlab('log Male:Female gene expression ratio') +
#   xlim(-1,1) +
#   scale_fill_manual(values = c('#9A1B1F')) +
#   labs(fill = 'Gametolog')
# ggsave('z_chromosome/figures/gametolog/gametolog species ridge.png',
#        height = 10,
#        width = 10)
# 
# # just Z
# data.z.long%>% 
#   ggplot(aes(x = log(gene_expression),
#              y = species)) +
#   ggridges::geom_density_ridges(fill = '#9A1B1F') +
#   geom_vline(xintercept = 0.69,
#              linetype = 'dashed')+
#   geom_vline(xintercept = -0.69,
#              linetype = 'dashed')+
# 
#   geom_vline(xintercept = 0)+
#   theme_classic() + 
#   theme(text = element_text(size = 24),
#         axis.text.y =  element_text(size = 16),
#         axis.text.x =  element_text(size = 16),
#         legend.position = 'null') +
#   xlab('log Male:Female gene expression ratio') +
#   xlim(-1,1) 
# ggsave('z_chromosome/figures/Z species ridge all.png',
#        height = 10,
#        width = 10)
# 
# # Z vs Autosomes
# data.long %>% 
#   ggplot(aes(x = log(gene_expression),
#              y = species,
#              fill = chromosome_type)) +
#   ggridges::geom_density_ridges() +
#   geom_vline(xintercept = 0.69,
#              linetype = 'dashed')+
#   geom_vline(xintercept = -0.69,
#              linetype = 'dashed')+
#   
#   geom_vline(xintercept = 0)+
#   theme_classic() + 
#   theme(text = element_text(size = 24),
#         axis.text.y =  element_text(size = 16),
#         axis.text.x =  element_text(size = 16),
#         legend.position = 'null'
#         ) +
#   xlab('log Male:Female gene expression ratio') +
#   xlim(-1,1) +
#   scale_fill_manual(values = c('Z' = '#9A1B1F',
#                                'Autosome' = 'grey'))
# ggsave('z_chromosome/figures/Z species ridge all.png',
#        height = 10,
#        width = 10)

# Z vs Autosomes
data.long %>% 
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
  ggplot(aes(x = log(gene_expression),
             y = reorder(species,
                         -order),
             fill = chromosome_type)) +
  ggridges::geom_density_ridges() +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  
  geom_vline(xintercept = 0)+
  theme_classic(base_size = 8) + 
  theme(
    # text = element_text(size = 24),
    #     axis.text.y =  element_text(size = 16),
    #     axis.text.x =  element_text(size = 16),
        legend.position = 'null'
  ) +
  xlab('log Male:Female gene expression ratio') +
  ylab('species') +
  xlim(-1,1) +
  scale_fill_manual(values = c('Z' = 'red',
                               # 'Z' = '#9A1B1F',
                               'Autosome' = 'grey'))
ggsave('z_chromosome/figures/Z species ridge all paper.pdf',
       height = 4,
       width = 3.1,
       units = 'in',
       dpi = 320)


### graph species average
# boxplot
data.z.long %>% 
  group_by(Ensembl.GeneIDs,
           gametolog) %>% 
  summarise(avg_gene_expression = mean(gene_expression)) %>% 
  ggplot(aes(y = log(avg_gene_expression),
             fill = as.factor(gametolog))) +
  geom_boxplot(color = 'black') +
  theme_classic() +
  ylab('log Male:Female gene expression ratio') +
  scale_fill_manual(values = c('#9A1B1F')) +
  labs(fill = 'Gametolog')
ggsave('z_chromosome/figures/gametolog/average gametolog species boxplot.png',
       height = 10,
       width = 10)

# ridge plot
data.z.long %>% 
  group_by(Ensembl.GeneIDs,
           gametolog) %>% 
  summarise(avg_gene_expression = mean(gene_expression)) %>% 
  ggplot(aes(x = log(avg_gene_expression),
             fill = as.factor(gametolog))) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  geom_density() +
  theme_classic() + 
  theme(text = element_text(size = 24),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null') +
  xlab('log Male:Female gene expression ratio') +
  xlim(-1,1) +
  scale_fill_manual(values = c('#9A1B1F')) +
  labs(fill = 'Gametolog')
ggsave('z_chromosome/figures/gametolog/average gametolog species ridge.png',
       height = 10,
       width = 10)

# ridge plot
data.z.long %>% 
  group_by(Ensembl.GeneIDs) %>% 
  summarise(avg_gene_expression = mean(gene_expression)) %>% 
  ggplot(aes(x = log(avg_gene_expression))) +
  geom_density(fill = '#9A1B1F') +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = 0)+
  theme_classic() + 
  theme(text = element_text(size = 24),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null') +
  xlab('Average M:F (log)') +
  xlim(-1,1)  
ggsave('z_chromosome/figures/Z average species ridge.png',
       height = 10,
       width = 10)

## compare to quantile data
data.quantile = read_tsv('CAGEE/data/normalizedCounts_ratio_quantile_ts.tsv')

# get quantile list for genes
data.quantile.list = data.quantile %>% 
  select(GeneDescription,
         GeneName,
         SAMPLETYPE) %>% 
  separate_wider_delim(SAMPLETYPE,
                       delim = '.',
                       names = c('chr',
                                 'quantile'),
                       too_many = 'merge') %>% 
  filter(chr == 'Z')
  
# add gametologs
data.quantile.list.gametolog = data.quantile.list %>% 
  left_join(data.wide.format.chr %>% 
  select(GeneDescription,
         GeneName,
         Ensembl.GeneIDs)) %>% 
  left_join(data.gametolog ) %>% 
  mutate(gametolog = ifelse(is.na(gametolog),
                            0,
                            gametolog),
         count = 1) 
  
 # graph counts per quantile 
data.quantile.list.gametolog %>% 
  group_by(gametolog,
           quantile) %>% 
  summarise(count = sum(count)) %>% 
  filter(gametolog == 1) %>% 
  ggplot(aes(x = quantile,
             y = count)) +
  geom_hline(yintercept = 35/5,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  ylim(0,12) +
  ggtitle('Tree swallow gametolog quantile') +
  ylab('Number of gametolog genes')
ggsave('z_chromosome/figures/gametolog/TS gametolog quantile.png') 
  
### compare gametolog expression level  
data.long.format.median.chr.gametolog = data.wide.format.median.chr.gametolog %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'median_gene_expression')
  
## graph boxplot comparing expression
data.long.format.median.chr.gametolog  %>%
  ggplot(aes(y = log(median_gene_expression),
             x = species,
             fill = as.factor(SAMPLETYPE))) +
  geom_boxplot(color = 'black') +
  theme_classic() +
  ylab('log median gene expression') +
  scale_fill_manual(values = c('#9A1B1F',
                               'grey')) +
  labs(fill = 'Gametolog')
ggsave('z_chromosome/figures/gametolog/gametolog species median boxplot.png',
       height = 10,
       width = 10)





#### graph ratio of gene expression across z chromosome ####
### graph along chromosome
## include gene length
data.z.long %>% 
  mutate(log.gene.expression = log(gene_expression),
         log.gene.expression.min = ifelse(log.gene.expression > 0,
                                          log.gene.expression,
                                          log.gene.expression - 0.1),
         log.gene.expression.max = ifelse(log.gene.expression > 0,
                                          log.gene.expression + 0.1,
                                          log.gene.expression),
         log.direction = ifelse(log.gene.expression > 0,
                                'Male_bias',
                                'Female_bias')) %>% 
  ggplot(aes(xmin = start_position,
             xmax = end_position,
             ymin = log.gene.expression.min,
             ymax = log.gene.expression.max,
             fill = log.direction)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.5,
             linetype = 'dashed') +
  geom_rect() +
  theme_classic() +
  facet_grid(species~.) +
  xlab('Z chromosome position')+
  ylab('log Male:Female gene expression') +
  theme(legend.position = 'null')
ggsave('z_chromosome/figures/gene expression ratio z chromosome across species length.png',
       height = 10,
       width = 5)

## genes as points
data.z.long %>% 
  mutate(log.gene.expression = log(gene_expression),
         log.direction = ifelse(log.gene.expression > 0,
                                'Male_bias',
                                'Female_bias')) %>% 
  ggplot(aes(x = start_position,
             xmax = end_position,
             y = log.gene.expression,
             color = log.direction)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.69,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  facet_grid(species~.) +
  xlab('Z chromosome position')+
  ylab('log Male:Female gene expression') +
  theme(legend.position = 'null')
ggsave('z_chromosome/figures/gene expression ratio z chromosome across species.png',
       height = 10,
       width = 5)

## genes as points
# add 2x line for double dose of z chromosome
data.z.long %>% 
  mutate(log.gene.expression = log(gene_expression),
         log.direction = case_when(log.gene.expression > 0.69 ~ "Male_bias",
                                   log.gene.expression < 0 ~ "Female_bias",
                                   TRUE ~ "none")) %>% 
  ggplot(aes(x = start_position,
             xmax = end_position,
             y = log.gene.expression,
             color = log.direction)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.69,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  facet_grid(species~.) +
  xlab('Z chromosome position')+
  ylab('log Male:Female gene expression') +
  theme(legend.position = 'null') +
  scale_color_manual(values = c('red',
                                'blue',
                                'grey'))
ggsave('z_chromosome/figures/gene expression ratio z chromosome across species color.png',
       height = 10,
       width = 5)
  
  
#### graph every chromosome ####
### create list of chromosomes
# remove NA
# remove RRCB
chr.list = data.wide.format.chr %>% 
  filter(!is.na(chromosome_name)) %>%
  filter(!grepl("RRCB",
                chromosome_name)) %>% 
  pull(chromosome_name) %>% 
    unique()
  
## create color list
cols.value = c('Male_bias' = "red",
               'Female_bias' = "blue",
               'none' = "grey")

### create graph across every chromosome of gene expression ratio
for (i in chr.list) {
## reduce columns, select chromosome, pivot longer
tmp = data.wide.format.chr %>% 
  dplyr::select(-c(GeneDescription,
                   SAMPLETYPE,
                   NCBI.GeneID,
                   Description,
                   Transcripts,
                   Ensembl.GeneIDs,
                   external_gene_name)) %>% 
  filter(chromosome_name == i) %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression')

# add 2x line for genes above 2x
tmp %>% 
  mutate(log.gene.expression = log(gene_expression),
         log.direction = case_when(log.gene.expression > 0.69 ~ "Male_bias",
                                   log.gene.expression < -0.69 ~ "Female_bias",
                                   TRUE ~ "none")) %>% 
  ggplot(aes(x = start_position,
             xmax = end_position,
             y = log.gene.expression,
             color = log.direction)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.69,
             linetype = 'dashed') +
  geom_hline(yintercept = -0.69,
             linetype = 'dashed') +
  geom_point() +
  theme_classic() +
  facet_grid(species~.) +
  xlab(paste(i,
             ' chromosome position',
             sep = '')) +
  ylab('log Male:Female gene expression') +
  theme(legend.position = 'null') +
  scale_color_manual(values = cols.value)
ggsave(paste('z_chromosome/figures/chr_position_sex_ratio/',
       i,
       ' chromosome gene expression ratio across species color.png',
       sep = ''),
       height = 10,
       width = 10)
}

### create ridge plot for each species
## create list of species
species.list = data.wide.format.chr %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression') %>% 
  pull(species) %>% 
  unique()

## create species long format data
# remove extra chromosome info
# convert ratio to log
# remove chromsome 16
# add chromosome size
data.long.format.chr = data.wide.format.chr %>% 
  dplyr::select(-c(GeneDescription,
                   SAMPLETYPE,
                   NCBI.GeneID,
                   Description,
                   Transcripts,
                   Ensembl.GeneIDs,
                   external_gene_name)) %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'gene_expression') %>% 
  filter(chromosome_name %in% chr.list) %>% 
  filter(chromosome_name != '16') %>% 
  mutate(log.gene.expression = log(gene_expression)) %>% 
  group_by(chromosome_name) %>% 
  mutate(chromosome_count = n())

## graph ridge plot of sex ratio for each chromosome for each species
for (i in species.list) {
  ## reduce columns, select chromosome, pivot longer
  tmp = data.long.format.chr %>% 
    filter(species == i) 
  
  # add 2x line for genes above 2x
  tmp %>% 
    droplevels() %>% 
    mutate(color = ifelse(chromosome_name == 'Z',
                          'Z',
                          'Auto')) %>% 
    ggplot(aes(x = log.gene.expression,
               y = reorder(chromosome_name,
                           -chromosome_count),
               fill = color)) +
    geom_vline(xintercept = 0.69,
               linetype = 'dashed')+
    geom_vline(xintercept = -0.69,
               linetype = 'dashed')+
    ggridges::geom_density_ridges() +
    theme_classic() +
    theme(legend.position = 'null') +
    xlab('log Male:Female gene expression') +
    ylab(paste(i,
               ' Chromosomes',
               sep = '')) +
    xlim(-1,1)+
    scale_fill_manual(values = c('grey',
                                 '#9A1B1F'))
  ggsave(paste('z_chromosome/figures/species_chr_ridge/',
               i,
               ' chromosome gene expression ratio ridgeplot.png',
               sep = ''),
         height = 10,
         width = 10)
}
  
### graph ridge for poster
## reduce columns, select chromosome, pivot longer
tmp.poster = data.long.format.chr %>% 
  filter(species == 'TS') 

# add 2x line for genes above 2x
tmp.poster %>% 
  droplevels() %>% 
  mutate(color = ifelse(chromosome_name == 'Z',
                        'Z',
                        'Auto')) %>% 
  ggplot(aes(x = log.gene.expression,
             y = reorder(chromosome_name,
                         -chromosome_count),
             fill = color)) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  theme(text = element_text(size = 24),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null') +
  xlab('log Male:Female gene expression ratio') +
  ylab(paste('Tree Swallow',
             ' Chromosomes',
             sep = '')) +
  xlim(-1,1) +
  scale_fill_manual(values = c('grey',
                               '#9A1B1F'))
ggsave('z_chromosome/figures/TS chromosome gene expression ratio ridgeplot poster.pdf',
       height = 5,
       width = 5)

## graph ancestral sex ratio
# get ancestral state at root 
ratio.ancestral.state.df = read_tsv('CAGEE/ratio/ratio_cagee_outs/ancestral_states.tab') %>% 
  dplyr::select(TranscriptID,
                '<11>') %>% 
  rename(root = '<11>') %>% 
  rename(GeneName = TranscriptID)

# graph 
data.long.format.chr %>% 
  droplevels() %>% 
  select(GeneName,
         chromosome_name,
         chromosome_count) %>% 
  distinct() %>% 
  left_join(ratio.ancestral.state.df) %>% 
  mutate(color = ifelse(chromosome_name == 'Z',
                        'Z',
                        'Auto'),
         log.root = log(root)) %>% 
  ggplot(aes(x = log.root,
             y = reorder(chromosome_name,
                         -chromosome_count),
             fill = color)) +
  geom_vline(xintercept = 0)+
  # geom_vline(xintercept = 0.69,
  #            linetype = 'dashed')+
  # geom_vline(xintercept = -0.69,
  #            linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  theme(
    # text = element_text(size = 24),
    #     axis.text.y =  element_text(size = 16),
    #     axis.text.x =  element_text(size = 16),
        legend.position = 'null') +
  xlab('log Male:Female gene expression ratio') +
  ylab('Chromosomes') +
  xlim(-0.38,0.78) +
  scale_fill_manual(values = c('grey',
                               '#9A1B1F'))
ggsave('z_chromosome/figures/Ancestral gene expression ratio chromosome ridgeplot.png',
       height = 5,
       width = 5)

# subset chromosomes 
data.long.format.chr %>% 
  droplevels() %>% 
  select(GeneName,
         chromosome_name,
         chromosome_count) %>% 
  distinct() %>% 
  left_join(ratio.ancestral.state.df) %>% 
  mutate(color = ifelse(chromosome_name == 'Z',
                        'Z',
                        'Auto'),
         log.root = log(root)) %>% 
  filter(chromosome_name %in% c('4','Z')) %>% 
  ggplot(aes(x = log.root,
             y = reorder(chromosome_name,
                         -chromosome_count),
             fill = color)) +
  geom_vline(xintercept = 0)+
  # geom_vline(xintercept = 0.69,
  #            linetype = 'dashed')+
  # geom_vline(xintercept = -0.69,
  #            linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  theme(
    text = element_text(size = 24),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
    legend.position = 'null') +
  xlab('log Male:Female gene expression ratio') +
  ylab('Chromosomes') +
  xlim(-0.38,0.78) +
  scale_fill_manual(values = c('grey',
                               '#9A1B1F'))
ggsave('z_chromosome/figures/Ancestral gene expression ratio chromosome ridgeplot Z and 4.png')
 
#### graph bar graph of number of sex biased genes
### set bias to be above or below 50%
data.long.format.chr = data.long.format.chr %>% 
  mutate(bias = case_when(log.gene.expression > 0.4055 ~ 'Male_bias',
                          log.gene.expression < -0.4055 ~ 'Female_bias',
                          TRUE ~ 'none'))


## get count data
data.long.format.chr.count = data.long.format.chr %>% 
  select(bias,
         chromosome_name,
         species) %>% 
  table() %>% 
  as.data.frame()

## graph boxplot
# all chromosomes
data.long.format.chr.count %>% 
  ggplot(aes(x = reorder(chromosome_name,
                         -Freq),
             y = Freq,
             fill = bias)) +
  geom_boxplot(position = 'dodge') +
  theme_classic() +
  scale_fill_manual(values = cols.value) +
  ylab('Number of genes') +
  xlab('Chromosome')
ggsave('z_chromosome/figures/bias per chromosome.png')

# just sex bias
data.long.format.chr.count %>% 
  filter(bias != 'none') %>% 
  ggplot(aes(x = reorder(chromosome_name,
                         -Freq),
             y = Freq,
             fill = bias)) +
  geom_boxplot(position = 'dodge') +
  theme_classic() +
  scale_fill_manual(values = cols.value)+
  ylab('Number of genes') +
  xlab('Chromosome')
ggsave('z_chromosome/figures/bias per chromosome sex bias.png')

# remove z chromosome
data.long.format.chr.count %>% 
  filter(bias != 'none') %>% 
  filter(chromosome_name != 'Z') %>% 
  ggplot(aes(x = reorder(chromosome_name,
                         -Freq),
             y = Freq,
             fill = bias)) +
  geom_boxplot(position = 'dodge') +
  theme_classic() +
  scale_fill_manual(values = cols.value)+
  ylab('Number of genes') +
  xlab('Chromosome')
ggsave('z_chromosome/figures/bias per chromosome sex bias no z.png')

























#### test in other species ####
### load bird data from: https://doi.org/10.1016/j.jgg.2020.10.005
## chicken 
data.chicken = read.csv('z_chromosome/Other_birds/Copy of Supplemental Table S3 - Gallus Gallus.csv')
# add species label
data.chicken = data.chicken %>% 
  mutate(species = 'chicken')

## turkey
data.turkey = read.csv('z_chromosome/Other_birds/Copy of Supplemental Table S3 - Meleagris gallopavo.csv')
# add species label
data.turkey = data.turkey %>% 
  mutate(species = 'turkey')

## mouse
data.mouse = read.csv('z_chromosome/Other_birds/Copy of Supplemental Table S3 - Mus musculus.csv')
# add species label
data.mouse = data.mouse %>% 
  mutate(species = 'mouse')


### pivot longer to separate out sex and species 
## chicken
data.chicken.long = data.chicken %>% 
  pivot_longer(cols = -c(species,
                         Gene_id,
                         Chr),
               names_to = 'type',
               values_to = 'expression') %>% 
  separate_wider_delim(type,
                       delim = '_',
                       names = c('sex',
                                 'tissue'),
                       too_many = "merge")

## turkey
data.turkey.long = data.turkey %>% 
  pivot_longer(cols = -c(species,
                         Gene_id,
                         Chr),
               names_to = 'type',
               values_to = 'expression') %>% 
  separate_wider_delim(type,
                       delim = '_',
                       names = c('sex',
                                 'tissue'),
                       too_many = "merge")

## mouse
data.mouse.long = data.mouse %>% 
  pivot_longer(cols = -c(species,
                         Gene_id,
                         Chr),
               names_to = 'type',
               values_to = 'expression') %>% 
  separate_wider_delim(type,
                       delim = '_',
                       names = c('sex',
                                 'tissue'),
                       too_many = "merge")

### calculate male to female ratio
## chicken 
# pivot wider
# remove na to keep tissue with both sexes
# remove zero counts 
# add chromosome count per tissue
data.chicken.long.ratio = data.chicken.long %>% 
  pivot_wider(names_from = sex,
              values_from = expression) %>% 
  na.omit() %>% 
  mutate(ratio = male/female,
         log.ratio = log(ratio)) %>% 
  filter(log.ratio != Inf) %>% 
  filter(log.ratio != -Inf) %>% 
  group_by(Chr,
           tissue) %>% 
  mutate(Chromosome.count = n())
  
## turkey 
# pivot wider
# remove na to keep tissue with both sexes
# remove zero counts 
# add chromosome count per tissue
data.turkey.long.ratio = data.turkey.long %>% 
  pivot_wider(names_from = sex,
              values_from = expression) %>% 
  na.omit() %>% 
  mutate(ratio = male/female,
         log.ratio = log(ratio)) %>% 
  filter(log.ratio != Inf) %>% 
  filter(log.ratio != -Inf) %>% 
  group_by(Chr,
           tissue) %>% 
  mutate(Chromosome.count = n()) 
  

## mouse 
# pivot wider
# remove na to keep tissue with both sexes
# remove zero counts 
# add chromosome count per tissue
data.mouse.long.ratio = data.mouse.long %>% 
  pivot_wider(names_from = sex,
              values_from = expression) %>% 
  na.omit() %>% 
  mutate(ratio = male/female,
         log.ratio = log(ratio)) %>% 
  filter(log.ratio != Inf) %>% 
  filter(log.ratio != -Inf) %>% 
  group_by(Chr,
           tissue) %>% 
  mutate(Chromosome.count = n()) 


### graph ridgeplot
## chicken
# all chromosome
data.chicken.long.ratio %>%
  ggplot(aes(x = log.ratio,
             y = reorder(Chr,
                         -Chromosome.count),
             fill = tissue)) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  xlab('log Male:Female gene expression') +
  ylab('Chicken Chromosomes') 
ggsave('z_chromosome/Other_birds/figures/chicken chromosome gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)

# z chromosome
data.chicken.long.ratio %>%
  filter(Chr == 'Z') %>% 
  ggplot(aes(x = log.ratio,
             y = reorder(tissue,
                         -Chromosome.count))) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  geom_text(aes(x = -5,
                y = tissue,
                label = Chromosome.count,
                hjust = 1)) +
  theme_classic() +
  xlim(-5,5) +
  xlab('Z chromosome log Male:Female gene expression') +
  ylab('Chicken tissue') 
ggsave('z_chromosome/Other_birds/figures/chicken z gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)

## turkey
# all chromosome
data.turkey.long.ratio %>%
  ggplot(aes(x = log.ratio,
             y = reorder(Chr,
                         -Chromosome.count),
             fill = tissue)) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  xlab('log Male:Female gene expression') +
  ylab('Turkey Chromosomes') 
ggsave('z_chromosome/Other_birds/figures/turkey chromosome gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)

# z chromosome
data.turkey.long.ratio %>%
  filter(Chr == 'Z') %>% 
  ggplot(aes(x = log.ratio,
             y = reorder(tissue,
                         -Chromosome.count))) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  geom_text(aes(x = -5,
                y = tissue,
                label = Chromosome.count,
                hjust = 1)) +
  theme_classic() +
  xlim(-5,5) +
  xlab('Z chromosome log Male:Female gene expression') +
  ylab('Turkey tissue') 
ggsave('z_chromosome/Other_birds/figures/turkey z gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)

## mouse
# all chromosome
data.mouse.long.ratio %>%
  ggplot(aes(x = log.ratio,
             y = reorder(Chr,
                         -Chromosome.count),
             fill = tissue)) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  theme_classic() +
  xlab('log Male:Female gene expression') +
  ylab('mouse Chromosomes') 
ggsave('z_chromosome/Other_birds/figures/mouse chromosome gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)

# x chromosome
data.mouse.long.ratio %>%
  filter(Chr == 'X') %>% 
  ggplot(aes(x = log.ratio,
             y = reorder(tissue,
                         -Chromosome.count))) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed')+
  geom_vline(xintercept = -0.69,
             linetype = 'dashed')+
  ggridges::geom_density_ridges() +
  geom_text(aes(x = -5,
                y = tissue,
                label = Chromosome.count,
                hjust = 1)) +
  theme_classic() +
  xlim(-5,5) +
  xlab('X chromosome log Male:Female gene expression') +
  ylab('mouse tissue') 
ggsave('z_chromosome/Other_birds/figures/mouse x gene expression ratio ridgeplot.png',
       height = 10,
       width = 10)










