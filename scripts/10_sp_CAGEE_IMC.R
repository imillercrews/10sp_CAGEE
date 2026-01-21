#### 10 sp bird analysis 
### CAGEE on sex difference
# R 4.2.1

## help making trees: yulab-smu.top/treedata-book/

## setwd
setwd("/geode2/home/u040/imillerc/Quartz/10_sp")

#### load libraries ####
library(ape)
library(tidytree)
library(ggtree)
library(tidyverse)
library(gghalves)
#### convert counts file to CAGEE input file ####
### load gene chromesome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

### load normalized counts file
## downloaded from CGB sharepoint 
#https://indiana.sharepoint.com/:x:/r/sites/msteams_863e21/Shared%20Documents/General/BL-OVPR-CGB/Labs/Rosvall/GSF2959_CavityNesting/foldchange_species/normalizedCounts.xlsx?d=webbe6cccfd6048258adf97282db1df10&csf=1&web=1&e=NWlnku
# saved as txt file
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

### calculate median for every gene in each sex in each species 
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

## calculate median for each sex and species
data.long.median = data.long %>% 
  dplyr::select(-c(sample_name,
                   sample_id)) %>% 
  group_by(species,
           sex) %>% 
  summarize_if(is.numeric,
            median)

# convert to true long
data.long.median = data.long.median %>% 
  pivot_longer(cols = -c(species,
                         sex),
               names_to = 'GeneName',
               values_to = 'median')

# pivot wider for CAGEE
data.wide.median = data.long.median %>% 
  pivot_wider(names_from = 'species',
              values_from = 'median')

## calculate median for each species
data.long.median.sp = data.long %>% 
  dplyr::select(-c(sample_name,
                   sample_id,
                   sex)) %>% 
  group_by(species) %>% 
  summarize_if(is.numeric,
               median)

# convert to true long
data.long.median.sp = data.long.median.sp %>% 
  pivot_longer(cols = -c(species),
               names_to = 'GeneName',
               values_to = 'median')

# pivot wider for CAGEE
data.long.median.sp = data.long.median.sp %>% 
  pivot_wider(names_from = 'species',
              values_from = 'median')

## calculate ratio for each species
# males/females
# add 1
data.long.ratio = data.long.median %>% 
  group_by(species,
           GeneName) %>% 
  summarize(ratio = (median[sex == 'M']+1)/(median[sex == 'F']+1))

# pivot wider for CAGEE
data.wide = data.long.ratio %>% 
  pivot_wider(names_from = 'species',
              values_from = 'ratio')

### formatting columns for use with CAGEE
## median species and sex
# add gene description
data.wide.format.median = full_join(data %>% 
                               dplyr::select(c(GeneDescription,
                                               GeneName)),
                             data.wide.median)
# rename SAMPLETYPE
data.wide.format.median = data.wide.format.median %>% 
  dplyr::rename(SAMPLETYPE = sex)

## median species
# add gene description
data.wide.format.median.sp = full_join(data %>% 
                                      dplyr::select(c(GeneDescription,
                                                      GeneName)),
                                    data.long.median.sp)

## ratio
# add gene description
data.wide.format = full_join(data %>% 
                        dplyr::select(c(GeneDescription,
                                        GeneName)),
                      data.wide)

### add z chromosome information
## make z chromosome gene list
# anything not on the z chromosome is A
data.gene.chromosome.z = data.gene.chromosome %>% 
  dplyr::select(Symbol,
                chromosome_name) %>% 
  mutate(SAMPLETYPE = case_when(chromosome_name == 'Z' ~ 'Z',
                                chromosome_name != 'Z' ~ 'A',
                                is.na(chromosome_name) ~ 'A',
                                TRUE ~ 'A')) %>% 
  dplyr::select(-c(chromosome_name)) %>% 
  dplyr::rename(GeneName = Symbol)

## add z chromosome to other data frames
# ratio
data.wide.format.z = data.wide.format %>% 
  left_join(data.gene.chromosome.z) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

# median sex and species 
# combine sex and chromosome
data.wide.format.median.z = data.wide.format.median %>% 
  left_join(data.gene.chromosome.z %>% 
              rename(SAMPLETYPE.Z = SAMPLETYPE)) %>% 
  mutate(SAMPLETYPE.Z = ifelse(is.na(SAMPLETYPE.Z),
                             'A',
                             SAMPLETYPE.Z)) %>% 
  mutate(SAMPLETYPE = paste0(SAMPLETYPE,
                             SAMPLETYPE.Z)) %>% 
  select(-c(SAMPLETYPE.Z)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

# median species
data.wide.format.median.sp.z = data.wide.format.median.sp %>% 
  left_join(data.gene.chromosome.z) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

## create random chromosomes for testing ratio data
# remove Z genes 
data.wide.format.rand = data.wide.format.z %>% 
  filter(SAMPLETYPE != 'Z') %>% 
  mutate(SAMPLETYPE = sample(25,
                             size = 10268,
                             replace = T),
         SAMPLETYPE = as.character(SAMPLETYPE))

### save file
## expression
# ratio
write_delim(data.wide.format,
            'CAGEE/data/normalizedCounts_ratio.tsv',
            delim = '\t')

# median sex and species
write_delim(data.wide.format.median,
            'CAGEE/data/normalizedCounts_median_sex.tsv',
            delim = '\t')

# median species
write_delim(data.wide.format.median.sp,
            'CAGEE/data/normalizedCounts_median.tsv',
            delim = '\t')

## z chromosome
# ratio
write_delim(data.wide.format.z,
            'CAGEE/data/normalizedCounts_ratio_z.tsv',
            delim = '\t')

# median sex and species
write_delim(data.wide.format.median.z,
            'CAGEE/data/normalizedCounts_median_sex_z.tsv',
            delim = '\t')

# median species
write_delim(data.wide.format.median.sp.z,
            'CAGEE/data/normalizedCounts_median_z.tsv',
            delim = '\t')

# ratio
write_delim(data.wide.format.rand,
            'CAGEE/data/normalizedCounts_ratio_z_rand.tsv',
            delim = '\t')

#### calculate quantile for ratio data ####
### load data 
## load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

## load ratio data
data.wide.format = read_delim(
            'CAGEE/data/normalizedCounts_ratio.tsv',
            delim = '\t')

## make all chromosome gene list
data.gene.chromosome.all = data.gene.chromosome %>% 
  dplyr::select(Symbol,
                chromosome_name) %>% 
  dplyr::rename( SAMPLETYPE = chromosome_name) %>% 
  dplyr::rename(GeneName = Symbol)


## add all chromosome to other data frames
# ratio
# assign genes with no chromosome data to 'A'
data.wide.format.all = data.wide.format %>% 
  left_join(data.gene.chromosome.all) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.all.quantile = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                    quantile(ratio, 
                                             probs = seq(0, 1, 0.2)), 
                                    rightmost.closed = TRUE)) 

# ## check that quantiles worked
# # count table
# data.wide.format.all.quantile %>%
#   filter(SAMPLETYPE == '4') %>%
#   select(quantile) %>%
#   table()
# # graph histogram
# data.wide.format.all.quantile %>%
#   filter(SAMPLETYPE == '4') %>%
#   filter(species =='TS') %>%
#   ggplot(aes(log(ratio),
#              fill = as.factor(quantile))) +
#   geom_histogram()+
#   theme_classic()
# # graph scatterplot
# data.wide.format.all.quantile %>%
#   filter(SAMPLETYPE == '4') %>%
#   filter(species =='TS') %>%
#   ggplot(aes(x = log(ratio),
#              y = quantile,
#              fill = as.factor(quantile))) +
#   geom_boxplot() +
#   theme_classic()

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.all.quantile.ts = data.wide.format.all.quantile %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.all.quantile.ts = data.wide.format %>% 
  left_join(data.wide.format.all.quantile.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 


# ## save data
# write_delim(data.wide.format.all.quantile.ts,
#   'CAGEE/data/normalizedCounts_ratio_quantile_ts.tsv',
#   delim = '\t')

### create shifted towards male biased dataset
## convert each ratio 
data.wide.format.all.quantile.shifted = data.wide.format.all.quantile %>% 
  mutate(ratio = ifelse(SAMPLETYPE == 'Z',
                        ratio,
                        exp(log(ratio)+0.5)
                        ))


# ## check that quantiles worked
# # count table
# data.wide.format.all.quantile.shifted %>%
#   filter(SAMPLETYPE == '4') %>%
#   select(quantile) %>%
#   table()
# # graph histogram
# data.wide.format.all.quantile.shifted %>%
#   filter(SAMPLETYPE == 'Z') %>%
#   filter(species =='TS') %>%
#   ggplot(aes(log(ratio),
#              fill = as.factor(quantile))) +
#   geom_histogram()+
#   theme_classic()
# # graph scatterplot
# data.wide.format.all.quantile.shifted %>%
#   filter(SAMPLETYPE == '4') %>%
#   filter(species =='TS') %>%
#   ggplot(aes(x = log(ratio),
#              y = quantile,
#              fill = as.factor(quantile))) +
#   geom_boxplot() +
#   theme_classic()

## create tree swallow list
data.wide.format.all.quantile.shifted.ts = data.wide.format.all.quantile.shifted %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## convert shifted data into wide format
data.wide.format.shifted = data.wide.format.all.quantile.shifted %>% 
  ungroup() %>% 
  select(-c(SAMPLETYPE,
            quantile)) %>% 
  pivot_wider(names_from = 'species',
              values_from = 'ratio')


## combine quantile with ratio data
data.wide.format.all.quantile.shifted.ts = data.wide.format.shifted %>% 
  left_join(data.wide.format.all.quantile.shifted.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

# ## save data
# write_delim(data.wide.format.all.quantile.shifted.ts,
#             'CAGEE/data/normalizedCounts_ratio_quantile_ts_shifted.tsv',
#             delim = '\t')

# list of chromosomes above threshold
data.wide.format.all %>% 
  select(SAMPLETYPE) %>% 
  group_by(SAMPLETYPE) %>% 
  summarise(count = n()) %>% 
  arrange(count) %>% 
  filter(count > 25) %>% 
  pull(SAMPLETYPE)

# pull z genes 
data.wide.format.z.quantile = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  filter(SAMPLETYPE == 'Z')

#### calculate chr Z quantile ratio for all bird species ####
# need to run previous chunk: "calculate quantile for ratio data"

### center Z chromosome data
## convert to log space and center on 0
# convert back to ratio 
data.wide.format.z.quantile.center = data.wide.format.z.quantile %>% 
  mutate(ratio = ratio %>% 
           log() %>% 
           scale(scale = FALSE,
                 center = TRUE) %>% 
           exp(),
         ratio = as.numeric(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.center = data.wide.format.z.quantile.center %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

### loop through each bird species and save
for (i in unique(data.wide.format.z.quantile.center$species)) {
# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
tmp.bird = data.wide.format.z.quantile.center %>% 
  filter(species == i) %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
tmp = data.wide.format.z.quantile.center %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(tmp.bird) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## save data
write_delim(tmp,
            paste0('CAGEE/data/species_z/',
            i,
            '_normalizedCounts_ratio_quantile_sim_center.tsv'),
            delim = '\t')

# remove files
rm(tmp)
rm(tmp.bird)

}

#### calculate chr 4 quantile ratio for all bird species ####
# need to run previous chunk: "calculate quantile for ratio data"

# pull chr 4 genes 
data.wide.format.4.quantile = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  filter(SAMPLETYPE == '4')

## assign quantiles per chromosome of sex ratio
# want to categori4e genes into groups per chromosome
data.wide.format.4.quantile = data.wide.format.4.quantile %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

### loop through each bird species and save
for (i in unique(data.wide.format.4.quantile$species)) {
  # pivot wider
  # create dataframe of just Tree Swallow (TS) quantile counts
  tmp.bird = data.wide.format.4.quantile %>% 
    filter(species == i) %>% 
    ungroup() %>% 
    mutate(SAMPLETYPE = paste(SAMPLETYPE,
                              quantile,
                              sep = '.')) %>% 
    dplyr::select(-c(species,
                     ratio,
                     quantile))
  
  ## combine quantile with ratio data
  tmp = data.wide.format.4.quantile %>% 
    pivot_wider(id_cols = c(GeneDescription,
                            GeneName),
                names_from = 'species',
                values_from = 'ratio') %>% 
    left_join(tmp.bird) %>% 
    relocate(SAMPLETYPE) %>% 
    relocate(GeneName) %>% 
    relocate(GeneDescription)
  
  ## save data
  write_delim(tmp,
              paste0('CAGEE/data/species_4/',
                     i,
                     '_normalizedCounts_ratio_quantile_sim_center.tsv'),
              delim = '\t')
  
  # remove files
  rm(tmp)
  rm(tmp.bird)
  
}

#### calculate rank for ratio data ####
### load data 
## load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

## load ratio data
data.wide.format = read_delim(
  'CAGEE/data/normalizedCounts_ratio.tsv',
  delim = '\t')

## make all chromosome gene list
data.gene.chromosome.all = data.gene.chromosome %>% 
  dplyr::select(Symbol,
                chromosome_name) %>% 
  dplyr::rename( SAMPLETYPE = chromosome_name) %>% 
  dplyr::rename(GeneName = Symbol)


## add all chromosome to other data frames
# ratio
# assign genes with no chromosome data to 'A'
data.wide.format.all = data.wide.format %>% 
  left_join(data.gene.chromosome.all) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

## assign rank per chromosome of sex ratio
# want to categorize genes into groups per chromosome and species
# calculate average rank of gene across all species
data.wide.format.all.rank = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(Rank = dense_rank(desc(ratio))) %>% 
  ungroup() %>% 
  group_by(SAMPLETYPE,
           GeneDescription) %>% 
  mutate(Rank.avg = mean(Rank))

### graph average rank per chromosome
## graph average rank per gene for all chromosomes
data.wide.format.all.rank %>% 
  group_by(SAMPLETYPE) %>% 
  mutate(Chr.count = n()) %>% 
  filter(Chr.count >= 200) %>% 
  select(GeneDescription,
         SAMPLETYPE,
         GeneName,
         Rank.avg,
         Chr.count) %>% 
  distinct() %>% 
  mutate(Rank.avg.center = scale(Rank.avg)) %>% 
  ggplot(aes(x = Rank.avg.center,
             y = reorder(SAMPLETYPE,
                         -Chr.count))) +
  ggridges::geom_density_ridges() +
  theme_classic() +
  ylab('Chromosome') +
  xlab('Scaled average rank per gene across all species')
ggsave('CAGEE/global_figures/chromsome_evolution/rank average per gene across chromosome.png')

# poster
## graph average rank per gene for Z, 4, 7
data.wide.format.all.rank %>%
  filter(SAMPLETYPE %in% c('Z',
                           '4',
                           '7',
                           '6',
                           '1A')) %>% 
  group_by(SAMPLETYPE) %>% 
  mutate(Chr.count = n()) %>% 
  filter(Chr.count >= 200) %>% 
  select(GeneDescription,
         SAMPLETYPE,
         GeneName,
         Rank.avg,
         Chr.count) %>% 
  distinct() %>% 
  mutate(Rank.avg.center = scale(Rank.avg)) %>% 
  mutate(color = ifelse(SAMPLETYPE == 'Z',
                        'Z',
                        'Auto')) %>% 
  ggplot(aes(x = -Rank.avg.center,
             y = reorder(SAMPLETYPE,
                         -Chr.count),
             fill = color)) +
  ggridges::geom_density_ridges() +
  theme_classic() +
  ylab('Chromosome') +
  xlab('Scaled average rank per gene across all species') +
  theme(title = element_text(size = 28),
        axis.title.x =  element_text(size = 28),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null') +
  scale_fill_manual(values = c('grey',
                               '#9A1B1F'))
ggsave('CAGEE/global_figures/chromsome_evolution/rank average per gene across z vs 4_7_6_1A.pdf',
       width = 7.5,
       height = 4.5)


### compare average ratio rank to gene expression

# load median sex and species
data.wide.format.median = read_delim('CAGEE/data/normalizedCounts_median_sex.tsv',
            delim = '\t')

## calculate average per gene for males and females across all species
data.wide.format.median.long = data.wide.format.median %>% 
  pivot_longer(cols = BB:YW,
               names_to = 'species',
               values_to = 'expression') %>% 
  rename(Sex = SAMPLETYPE)

# add chromosome data
data.wide.format.median.long = data.wide.format.median.long %>% 
  left_join(data.gene.chromosome.all) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 


## assign rank gene expression per chromosome for each sex 
# want to categorize genes into groups per chromosome and species
# calculate average rank of gene across all species
data.wide.format.median.long.rank = data.wide.format.median.long %>% 
  group_by(SAMPLETYPE,
           species,
           Sex) %>% 
  mutate(Rank = dense_rank(desc(expression))) %>% 
  ungroup() %>% 
  group_by(SAMPLETYPE,
           GeneDescription,
           Sex) %>% 
  mutate(Rank.avg.exp = mean(Rank))

# reduce 
data.wide.format.median.long.rank.reduce = data.wide.format.median.long.rank %>% 
  select(-c(expression,
            Rank,
            species)) %>% 
  distinct() %>% 
  pivot_wider(names_from = Sex,
              values_from = Rank.avg.exp)
  


## combine ratio rank and expression rank data 
data.wide.format.all.rank.avg = data.wide.format.all.rank %>% 
  select(-c(species,
            ratio,
            Rank)) %>% 
  distinct() %>% 
  full_join(data.wide.format.median.long.rank.reduce)



### graph average rank vs expression for each chromosome 
for (i in unique(data.wide.format.all.rank.avg$SAMPLETYPE)) {
data.wide.format.all.rank.avg %>% 
  filter(SAMPLETYPE == i) %>% 
  ungroup() %>% 
  mutate(Rank.avg.center = scale(Rank.avg)) %>% 
  mutate(Rank.avg.expression = ifelse(Rank.avg.center > 0,
                                      M,
                                      F),
         Rank.avg.expression.center = scale(Rank.avg.expression)) %>% 
  ggplot(aes(y = Rank.avg.expression.center,
             x = Rank.avg.center)) + 
  geom_point() +
  theme_classic() +
  labs(x = 'Average rank M:F ratio scaled',
       y = 'Average rank expression scaled',
       title = i) +
    xlim(-5,5) +
    ylim(-4,4)
ggsave(paste0('CAGEE/global_figures/chromsome_evolution/average_rank/',
              i,
              ' avg rank expression vs average rank ratio.png'))

}



#### simulate z chromosome data ####
### create simulated z chromosome data
## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.tmp = data.wide.format.z.quantile %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.tmp.ts = data.wide.format.z.quantile.tmp %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.tmp = data.wide.format.z.quantile.tmp %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.tmp.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.tmp %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.tmp %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

#### simulate z chromosome ratios
### uniform distribution
## pull from uniform distribution centered around 0
# use a range of 50% male to female bias (from 0.66667 to 1.5)
data.wide.format.z.quantile.sim = data.wide.format.z.quantile %>% 
  mutate(ratio = runif(nrow(data.wide.format.z.quantile),
                       min = 0.66667,
                       max = 1.5))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.sim = data.wide.format.z.quantile.sim %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.sim.ts = data.wide.format.z.quantile.sim %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.sim = data.wide.format.z.quantile.sim %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.sim.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.sim %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.sim %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.sim,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim.tsv',
            delim = '\t')

### simulate z chromosome ratios
## pull from uniform distribution centered around exp(0.5)
# use a range similar to that of 50% male to female bias
# shift max and min of previous by 0.5
data.wide.format.z.quantile.sim_0.5 = data.wide.format.z.quantile %>% 
  mutate(ratio = runif(nrow(data.wide.format.z.quantile),
                       min = exp(log(0.66667)+0.5),
                       max = exp(log(1.5)+0.5)))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.sim_0.5 = data.wide.format.z.quantile.sim_0.5 %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.sim.ts_0.5 = data.wide.format.z.quantile.sim_0.5 %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.sim_0.5 = data.wide.format.z.quantile.sim_0.5 %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.sim.ts_0.5) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.sim_0.5 %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.sim_0.5 %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.sim_0.5,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_0.5.tsv',
            delim = '\t')

### normal distribution
## calculate standard deviation for autosomes
# use average standard deviation from autosomes 4 and 7 from Tree swallow
# need the log scale sd 
data.wide.format.all  %>% 
  filter(SAMPLETYPE %in% c('4',
                           '7')) %>% 
  select(SAMPLETYPE,
         TS) %>% 
  group_by(SAMPLETYPE) %>% 
  summarise(sd = sd(log(TS))) %>% 
  pull(sd) %>% 
  mean()
# 0.1645212

## pull from normal distribution centered around 1:1
# use standard deviation from tree swallow chromosome 4 and 7
# convert back to ratio 
data.wide.format.z.quantile.sim.norm = data.wide.format.z.quantile %>% 
  mutate(ratio = rnorm(nrow(data.wide.format.z.quantile),
                       mean = 0,
                       sd = 0.1645212),
         ratio = exp(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.sim.norm = data.wide.format.z.quantile.sim.norm %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.sim.norm.ts = data.wide.format.z.quantile.sim.norm %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.sim.norm = data.wide.format.z.quantile.sim.norm %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.sim.norm.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.sim.norm %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.sim.norm %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.sim.norm,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim.norm.tsv',
            delim = '\t')

## pull from normal distribution centered around 0.5 log scale
# use standard deviation from tree swallow chromosome 4 and 7
# convert back to ratio 
data.wide.format.z.quantile.sim.norm_0.5 = data.wide.format.z.quantile %>% 
  mutate(ratio = rnorm(nrow(data.wide.format.z.quantile),
                       mean = 0.5,
                       sd = 0.1645212),
         ratio = exp(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.sim.norm_0.5 = data.wide.format.z.quantile.sim.norm_0.5 %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.sim.norm_0.5.ts = data.wide.format.z.quantile.sim.norm_0.5 %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.sim.norm_0.5 = data.wide.format.z.quantile.sim.norm_0.5 %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.sim.norm_0.5.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.sim.norm_0.5 %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.sim.norm_0.5 %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.sim.norm_0.5,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_0.5.norm.tsv',
            delim = '\t')

### center Z chromosome data
## convert to log space and center on 0
# convert back to ratio 
data.wide.format.z.quantile.center = data.wide.format.z.quantile %>% 
  mutate(ratio = ratio %>% 
           log() %>% 
           scale(scale = FALSE,
                 center = TRUE) %>% 
           exp(),
         ratio = as.numeric(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.center = data.wide.format.z.quantile.center %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.center.ts = data.wide.format.z.quantile.center %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.center = data.wide.format.z.quantile.center %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.center.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.center %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.center %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.center,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_center.tsv',
            delim = '\t')

### center and scale variance Z chromosome data
## calculate standard deviation of all chromosomes
data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  # filter(SAMPLETYPE == '4') %>% 
  group_by(SAMPLETYPE) %>% 
  summarise(tmp = sd(log(ratio))) %>% view()

# standard deviation
# chr 4 = 0.1760972
# chr 7 = 0.1663503
# average of chr 4 & 7 = 0.1712238
# chr z = 0.3059818

## convert to log space and center on 0
# convert back to ratio 
data.wide.format.z.quantile.scale = data.wide.format.z.quantile %>% 
  mutate(ratio = ratio %>% 
           log() %>% 
           scale(scale = FALSE,
                 center = TRUE),
         ratio = ratio*0.1712238/0.3059818,
         ratio = exp(ratio),
         ratio = as.numeric(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.scale = data.wide.format.z.quantile.scale %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.scale.ts = data.wide.format.z.quantile.scale %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.scale = data.wide.format.z.quantile.scale %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.scale.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

# ## check that quantiles worked
# # count table
# data.wide.format.z.quantile.scale %>%
#   select(SAMPLETYPE) %>%
#   table()
# # graph histogram
# data.wide.format.z.quantile.scale %>%
#   ggplot(aes(log(TS),
#              fill = SAMPLETYPE)) +
#   geom_histogram()+
#   theme_classic()

## save data
write_delim(data.wide.format.z.quantile.scale,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_scale.tsv',
            delim = '\t')

### Remove gametologs Z chromosome data
## load gametolog data
data.gametolog = read_tsv('Gene_list/10sp_ensembl_gametologs_edit.tsv') %>% 
  select(-c(Ensembl.GeneIDs)) %>% 
  rename(GeneName = Symbol) %>% 
  distinct()

## convert to log space and center on 0
# convert back to ratio 
data.wide.format.z.quantile.game = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  left_join(data.gametolog) %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'Z',
                             SAMPLETYPE)) %>% 
  filter(SAMPLETYPE == 'Z') %>% 
  filter(is.na(gametolog)) %>% 
  select(-c(gametolog)) 

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.game = data.wide.format.z.quantile.game %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.game.ts = data.wide.format.z.quantile.game %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.game = data.wide.format.z.quantile.game %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.game.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.wide.format.z.quantile.game %>%
  select(SAMPLETYPE) %>%
  table()
# graph histogram
data.wide.format.z.quantile.game %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic()

## save data
write_delim(data.wide.format.z.quantile.game,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_game.tsv',
            delim = '\t')

### Remove gametologs, center Z chromosome data
## load gametolog data
data.gametolog = read_tsv('Gene_list/10sp_ensembl_gametologs_edit.tsv') %>% 
  select(-c(Ensembl.GeneIDs)) %>% 
  rename(GeneName = Symbol) %>% 
  distinct()

## convert to log space and center on 0
# convert back to ratio 
data.wide.format.z.quantile.center.game = data.wide.format.all %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName,
                         SAMPLETYPE),
               names_to = 'species',
               values_to = 'ratio') %>% 
  left_join(data.gametolog) %>% 
  mutate(SAMPLETYPE = ifelse(!is.na(gametolog),
                             'Z',
                             SAMPLETYPE)) %>% 
  filter(SAMPLETYPE == 'Z') %>% 
  filter(is.na(gametolog)) %>% 
  select(-c(gametolog)) %>% 
   mutate(ratio = ratio %>% 
           log() %>% 
           scale(scale = FALSE,
                 center = TRUE) %>% 
           exp(),
         ratio = as.numeric(ratio))

## assign quantiles per chromosome of sex ratio
# want to categorize genes into groups per chromosome
data.wide.format.z.quantile.center.game = data.wide.format.z.quantile.center.game %>% 
  group_by(SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE)) 

# pivot wider
# create dataframe of just Tree Swallow (TS) quantile counts
data.wide.format.z.quantile.center.game.ts = data.wide.format.z.quantile.center.game %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.wide.format.z.quantile.center.game = data.wide.format.z.quantile.center.game %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.wide.format.z.quantile.center.game.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.wide.format.z.quantile.center.game %>%
  select(SAMPLETYPE) %>%
  table()
# graph histogram
data.wide.format.z.quantile.center.game %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic()

## save data
write_delim(data.wide.format.z.quantile.center.game,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_center_game.tsv',
            delim = '\t')


#### graph simulated data
data.wide.format.z.quantile.sim %>%
  mutate(bias = 'No bias',
         simulation = 'uniform') %>% 
  rbind(data.wide.format.z.quantile.sim_0.5 %>%
          mutate(bias = 'Male bias',
                 simulation = 'uniform')) %>% 
  rbind(data.wide.format.z.quantile.sim.norm %>%
          mutate(bias = 'No bias',
                 simulation = 'normal')) %>% 
  rbind(data.wide.format.z.quantile.sim.norm_0.5 %>%
          mutate(bias = 'Male bias',
                 simulation = 'normal')) %>% 
  rbind(data.wide.format.z.quantile.center %>%
          mutate(bias = 'Real',
                 simulation = 'center')) %>% 
  rbind(data.wide.format.z.quantile.tmp %>%
          mutate(bias = 'Real',
                 simulation = 'none')) %>% 
  rbind(data.wide.format.z.quantile.scale %>%
          mutate(bias = 'Real',
                 simulation = 'scale')) %>% 
  rbind(data.wide.format.z.quantile.game %>%
          mutate(bias = 'Gametolog removed',
                 simulation = 'none')) %>% 
  rbind(data.wide.format.z.quantile.center.game %>%
          mutate(bias = 'Gametolog removed',
                 simulation = 'center')) %>% 
 ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_vline(xintercept = 0) +
  geom_histogram()+
  theme_classic() +
  facet_grid(bias + simulation~.) 
ggsave('CAGEE/global_figures/chromsome_evolution/histogram simulated z.png')

data.wide.format.z.quantile.center %>%
          mutate(bias = 'Real',
                 simulation = 'center') %>% 
  rbind(data.wide.format.z.quantile.tmp %>%
          mutate(bias = 'Real',
                 simulation = 'none')) %>% 
  rbind(data.wide.format.z.quantile.scale %>%
          mutate(bias = 'Real',
                 simulation = 'scale')) %>% 
  rbind(data.wide.format.z.quantile.game %>%
          mutate(bias = 'Gametolog removed',
                 simulation = 'none')) %>% 
  rbind(data.wide.format.z.quantile.center.game %>%
          mutate(bias = 'Gametolog removed',
                 simulation = 'center')) %>% 
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_vline(xintercept = 0) +
  geom_histogram()+
  theme_classic() +
  facet_grid(simulation + bias~.) 
ggsave('CAGEE/global_figures/chromsome_evolution/histogram simulated z real.png')

#### simulate evolution z chromosome ####
#### create root values for CAGEE simulation
### pull random tree swallow data to treat as root of tree for simulation
# needs to be in format: 
# "tab-separated file specifying value(col 1) and root distribution weight (col 2)"
## uniform distribution
# no bias
data.wide.format.z.quantile.sim %>% 
  select(TS) %>% 
  mutate(weight = 1) %>%
  write_delim('CAGEE/data/sim/normalizedCounts_ratio_sim_ts_root.tsv',
              delim = '\t',
              col_names = F)

# male bias
data.wide.format.z.quantile.sim %>% 
  select(TS) %>% 
  mutate(weight = 1) %>%
  write_delim('CAGEE/data/sim/normalizedCounts_ratio_sim_ts_root_0.5.tsv',
              delim = '\t',
              col_names = F)

## normal distribution
# no bias
data.wide.format.z.quantile.sim.norm %>% 
  select(TS) %>% 
  mutate(weight = 1) %>%
  write_delim('CAGEE/data/sim/normalizedCounts_ratio_sim_ts_root_norm.tsv',
              delim = '\t',
              col_names = F)

# male bias
data.wide.format.z.quantile.sim.norm_0.5 %>% 
  select(TS) %>% 
  mutate(weight = 1) %>%
  write_delim('CAGEE/data/sim/normalizedCounts_ratio_sim_ts_root_norm_0.5.tsv',
              delim = '\t',
              col_names = F)


#### run these root distribution data sets in CAGEE simulation


#### calculate simulated evolution z chromosome ratio quantiles ####
### load simulated data
### note: use comment to remove first 4 rows from txt file (get rid of metadata in hashtag rows)
## uniform distribution
# no bias
data.evo.z.uni = read_delim('CAGEE/sim/evo/ratio_z_sim_ts_cagee_outs/simulation.txt',
            delim = '\t',
            comment = '#')
# male bias
data.evo.z.uni.0.5 = read_delim('CAGEE/sim/evo/ratio_z_sim_0.5_ts_cagee_outs/simulation.txt',
                            delim = '\t',
                            comment = '#')
## normal distribution
# no bias
data.evo.z.norm = read_delim('CAGEE/sim/evo/ratio_z_sim_norm_ts_cagee_outs/simulation.txt',
                            delim = '\t',
                            comment = '#')
# male bias
data.evo.z.norm.0.5 = read_delim('CAGEE/sim/evo/ratio_z_sim_norm_0.5_ts_cagee_outs/simulation.txt',
                            delim = '\t',
                            comment = '#')
  
#### assign quantiles per chromosome of sex ratio
### uniform distribution
## no bias
data.evo.z.uni =  data.evo.z.uni %>% 
  pivot_longer(cols = -c(DESC,
                         GENE_ID),
               names_to = 'species',
               values_to = 'ratio') %>% 
  mutate(SAMPLETYPE = 'Z') %>% 
  group_by(
    # SAMPLETYPE,
           species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE
                                 )) %>% 
  rename(GeneDescription = DESC) %>% 
  rename(GeneName = GENE_ID)  %>% 
  ungroup()

# create dataframe of just Tree Swallow (TS) quantile counts
data.evo.z.uni.ts = data.evo.z.uni %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.evo.z.uni = data.evo.z.uni %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.evo.z.uni.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.evo.z.uni %>%
  select(SAMPLETYPE) %>%
  table()

# graph histogram
data.evo.z.uni %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic() 

## save data
write_delim(data.evo.z.uni,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_evo.tsv',
            delim = '\t')

### uniform distribution
## male bias
data.evo.z.uni.0.5 =  data.evo.z.uni.0.5 %>% 
  pivot_longer(cols = -c(DESC,
                         GENE_ID),
               names_to = 'species',
               values_to = 'ratio') %>% 
  mutate(SAMPLETYPE = 'Z') %>% 
  group_by(
    # SAMPLETYPE,
    species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE
  )) %>% 
  rename(GeneDescription = DESC) %>% 
  rename(GeneName = GENE_ID)  %>% 
  ungroup()

# create dataframe of just Tree Swallow (TS) quantile counts
data.evo.z.uni.0.5.ts = data.evo.z.uni.0.5 %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.evo.z.uni.0.5 = data.evo.z.uni.0.5 %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.evo.z.uni.0.5.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.evo.z.uni.0.5 %>%
  select(SAMPLETYPE) %>%
  table()

# graph histogram
data.evo.z.uni.0.5 %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic() 

## save data
write_delim(data.evo.z.uni.0.5,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_0.5_evo.tsv',
            delim = '\t')

### normal distribution
## no bias
data.evo.z.norm =  data.evo.z.norm %>% 
  pivot_longer(cols = -c(DESC,
                         GENE_ID),
               names_to = 'species',
               values_to = 'ratio') %>% 
  mutate(SAMPLETYPE = 'Z') %>% 
  group_by(
    # SAMPLETYPE,
    species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE
  )) %>% 
  rename(GeneDescription = DESC) %>% 
  rename(GeneName = GENE_ID)  %>% 
  ungroup()

# create dataframe of just Tree Swallow (TS) quantile counts
data.evo.z.norm.ts = data.evo.z.norm %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.evo.z.norm = data.evo.z.norm %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.evo.z.norm.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.evo.z.norm %>%
  select(SAMPLETYPE) %>%
  table()

# graph histogram
data.evo.z.norm %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic() 

## save data
write_delim(data.evo.z.norm,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_norm_evo.tsv',
            delim = '\t')

### normal distribution
## male bias
data.evo.z.norm.0.5 =  data.evo.z.norm.0.5 %>% 
  pivot_longer(cols = -c(DESC,
                         GENE_ID),
               names_to = 'species',
               values_to = 'ratio') %>% 
  mutate(SAMPLETYPE = 'Z') %>% 
  group_by(
    # SAMPLETYPE,
    species) %>% 
  mutate(quantile = findInterval(ratio, 
                                 quantile(ratio, 
                                          probs = seq(0, 1, 0.2)), 
                                 rightmost.closed = TRUE
  )) %>% 
  rename(GeneDescription = DESC) %>% 
  rename(GeneName = GENE_ID)  %>% 
  ungroup()

# create dataframe of just Tree Swallow (TS) quantile counts
data.evo.z.norm.0.5.ts = data.evo.z.norm.0.5 %>% 
  filter(species == 'TS') %>% 
  ungroup() %>% 
  mutate(SAMPLETYPE = paste(SAMPLETYPE,
                            quantile,
                            sep = '.')) %>% 
  dplyr::select(-c(species,
                   ratio,
                   quantile))

## combine quantile with ratio data
data.evo.z.norm.0.5 = data.evo.z.norm.0.5 %>% 
  pivot_wider(id_cols = c(GeneDescription,
                          GeneName),
              names_from = 'species',
              values_from = 'ratio') %>% 
  left_join(data.evo.z.norm.0.5.ts) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription)

## check that quantiles worked
# count table
data.evo.z.norm.0.5 %>%
  select(SAMPLETYPE) %>%
  table()

# graph histogram
data.evo.z.norm.0.5 %>%
  ggplot(aes(log(TS),
             fill = SAMPLETYPE)) +
  geom_histogram()+
  theme_classic() 

## save data
write_delim(data.evo.z.norm.0.5,
            'CAGEE/data/normalizedCounts_ratio_quantile_ts_sim_norm_0.5_evo.tsv',
            delim = '\t')


#### load simulated ratio data ####
### load results
## uniform
# results.txt z chromosome simulated
data.cagee.results.z.sim = read.delim('CAGEE/ratio/ratio_z_sim_quantile_ts_cagee_outs/results.txt')

# results.txt z chromosome simulated
# bias
data.cagee.results.z.sim.0.5 = read.delim('CAGEE/ratio/ratio_z_sim_0.5_quantile_ts_cagee_outs/results.txt')

## norm
# results.txt z chromosome simulated
data.cagee.results.z.sim.norm = read.delim('CAGEE/ratio/ratio_z_sim_norm_quantile_ts_cagee_outs/results.txt')

# results.txt z chromosome simulated
# bias
data.cagee.results.z.sim.0.5.norm = read.delim('CAGEE/ratio/ratio_z_sim_norm_0.5_quantile_ts_cagee_outs/results.txt')

## centered z chromosome
# results.txt z chromosome simulated
data.cagee.results.z.sim.center = read.delim('CAGEE/ratio/ratio_z_sim_center_quantile_ts_cagee_outs/results.txt')
## centered and scaled z chromosome
# results.txt z chromosome simulated
data.cagee.results.z.sim.scale = read.delim('CAGEE/ratio/ratio_z_sim_scale_quantile_ts_cagee_outs/results.txt')

### evolution
## uniform
# results.txt z chromosome simulated
data.cagee.results.z.sim.evo = read.delim('CAGEE/sim/ratio/ratio_z_sim_evo_quantile_ts_cagee_outs/results.txt')

# results.txt z chromosome simulated
# bias
data.cagee.results.z.sim.evo.0.5 = read.delim('CAGEE/sim/ratio/ratio_z_sim_evo_0.5_quantile_ts_cagee_outs/results.txt')
## norm
# results.txt z chromosome simulated
data.cagee.results.z.sim.evo.norm = read.delim('CAGEE/sim/ratio/ratio_z_sim_norm_evo_quantile_ts_cagee_outs/results.txt')
# results.txt z chromosome simulated
# bias
data.cagee.results.z.sim.evo.norm.0.5 = read.delim('CAGEE/sim/ratio/ratio_z_sim_evo_norm_0.5_quantile_ts_cagee_outs/results.txt')

## results.txt 4 chromosome
# shifted
data.cagee.results.4.shifted = read.delim('CAGEE/ratio/ratio_4_quantile_shifted_ts_cagee_outs/results.txt')

## results.txt 7 chromosome
# shifted
data.cagee.results.7.shifted = read.delim('CAGEE/ratio/ratio_7_quantile_shifted_ts_cagee_outs/results.txt')

## results.txt Z chromosome
# gametolog
data.cagee.results.z.sim.game = read.delim('CAGEE/ratio/ratio_z_sim_game_quantile_ts_cagee_outs/results.txt')

## results.txt z chromosome
# centered, gametolog
data.cagee.results.z.sim.center.game = read.delim('CAGEE/ratio/ratio_z_sim_center_game_quantile_ts_cagee_outs/results.txt')



## uniform
## convert to table
# z chromosome simulated
data.cagee.results.z.sim = data.cagee.results.z.sim %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome simulated
# bias
data.cagee.results.z.sim.0.5 = data.cagee.results.z.sim.0.5 %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))


## normal
## convert to table
# z chromosome simulated
data.cagee.results.z.sim.norm = data.cagee.results.z.sim.norm %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome simulated
# bias
data.cagee.results.z.sim.0.5.norm = data.cagee.results.z.sim.0.5.norm %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome center
# bias
data.cagee.results.z.sim.center = data.cagee.results.z.sim.center %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))
# z chromosome center
# bias
data.cagee.results.z.sim.scale = data.cagee.results.z.sim.scale %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome 
# gametolog
data.cagee.results.z.sim.game = data.cagee.results.z.sim.game %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))


# z chromosome center
# gametolog
data.cagee.results.z.sim.center.game = data.cagee.results.z.sim.center.game %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## evolution
## uniform
## convert to table
# z chromosome simulated
data.cagee.results.z.sim.evo = data.cagee.results.z.sim.evo %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome simulated
# bias
data.cagee.results.z.sim.evo.0.5 = data.cagee.results.z.sim.evo.0.5 %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))


## normal
## convert to table
# z chromosome simulated
data.cagee.results.z.sim.evo.norm = data.cagee.results.z.sim.evo.norm %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# z chromosome simulated
# bias
data.cagee.results.z.sim.evo.norm.0.5 = data.cagee.results.z.sim.evo.norm.0.5 %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# 4 chromosome shifted
data.cagee.results.4.shifted = data.cagee.results.4.shifted %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# 7 chromosome shifted
data.cagee.results.7.shifted = data.cagee.results.7.shifted %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## combine results
# create quantile column
data.cagee.results.sim.table = data.cagee.results.z.sim %>% 
  filter(grepl("Sigma",
               CAGEE.1.1)) %>% 
  mutate(bias = 'No bias',
         simulation = 'uniform') %>% 
  rbind(data.cagee.results.z.sim.0.5 %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Male bias',
                 simulation = 'uniform'))  %>% 
  rbind(data.cagee.results.z.sim.norm %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'No bias',
                 simulation = 'normal'))  %>% 
  rbind(data.cagee.results.z.sim.0.5.norm %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Male bias',
                 simulation = 'normal'))  %>% 
  rbind(data.cagee.results.z.sim.center %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Z',
                 simulation = 'shifted'))  %>% 
  rbind(data.cagee.results.z.sim.center.game %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'ZG',
                 simulation = 'No bias'))  %>% 
  rbind(data.cagee.results.z.sim.game %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'ZG',
                 simulation = 'center'))  %>% 
  rbind(data.cagee.results.z.sim.scale %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Z.scale',
                 simulation = 'shifted'))  %>% 
  rbind(data.cagee.results.z.sim.evo %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'No bias',
                 simulation = 'evo.uniform')) %>% 
  rbind(data.cagee.results.z.sim.evo.0.5 %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Male bias',
                 simulation = 'evo.uniform')) %>% 
  rbind(data.cagee.results.z.sim.evo.norm %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'No bias',
                 simulation = 'evo.normal')) %>% 
  rbind(data.cagee.results.z.sim.evo.norm.0.5 %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = 'Male bias',
                 simulation = 'evo.normal')) %>% 
  rbind(data.cagee.results.4.shifted %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = '4',
                 simulation = 'shifted'))  %>% 
  rbind(data.cagee.results.7.shifted %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1)) %>% 
          mutate(bias = '7',
                 simulation = 'shifted'))  %>% 
    separate_wider_delim(CAGEE.1.1,
                       delim = 'Sigma for ',
                       names = c(NA,
                                 'chr')) %>% 
  separate_wider_delim(chr,
                       delim = '.',
                       names = c('chr',
                                 'quantile')) %>% 
  mutate(sigma = round(as.numeric(result),
                       digits = 8))

#### graph simulated ratio data ####
### graph results
## sigma rate across quantiles
data.cagee.results.sim.table %>% 
  filter(bias != 'ZG') %>% 
  mutate(id = paste0(simulation,
                     bias)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = simulation,
             shape = bias)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles simulated' ,
          subtitle = '(Tree Swallow quantiles)')+
  scale_color_manual(values = c('red',
                                'black',
                                'blue',
                                'grey',
                                'lightblue'))
ggsave('CAGEE/global_figures/chromsome_evolution/quantile z simulated chromsome.png',
       height = 10,
       width = 10)

## sigma rate across quantiles
data.cagee.results.sim.table %>% 
  filter(bias != 'ZG') %>% 
  mutate(id = paste0(simulation,
                     bias)) %>% 
  filter(simulation != 'uniform') %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = simulation,
             shape = bias)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles simulated' ,
          subtitle = '(Tree Swallow quantiles)')+
  scale_color_manual(values = c('red',
                                'black',
                                'blue',
                                'grey',
                                'lightblue'))
ggsave('CAGEE/global_figures/chromsome_evolution/quantile z simulated chromsome outlier.png',
       height = 10,
       width = 10)

#### graph quantile ratio data ####
### load quantile results chromosome data 
### load results
## results.txt z chromosome
data.cagee.results.z = read.delim('CAGEE/ratio/ratio_z_quantile_ts_cagee_outs/results.txt')

## results.txt 4 chromosome
data.cagee.results.4 = read.delim('CAGEE/ratio/ratio_4_quantile_ts_cagee_outs/results.txt')

## results.txt 7 chromosome
data.cagee.results.7 = read.delim('CAGEE/ratio/ratio_7_quantile_ts_cagee_outs/results.txt')


## convert to table
# z chromosome
data.cagee.results.z = data.cagee.results.z %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# 4 chromosome
data.cagee.results.4 = data.cagee.results.4 %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

# 7 chromosome
data.cagee.results.7 = data.cagee.results.7 %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))



## combine results
# create quantile column
data.cagee.results.table = data.cagee.results.z %>% 
  filter(grepl("Sigma",
               CAGEE.1.1)) %>% 
  rbind(data.cagee.results.7 %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1))) %>% 
    rbind(data.cagee.results.4 %>% 
          filter(grepl("Sigma",
                       CAGEE.1.1))) %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = 'Sigma for ',
                       names = c(NA,
                                 'chr')) %>% 
  separate_wider_delim(chr,
                       delim = '.',
                       names = c('chr',
                                 'quantile')) %>% 
  mutate(sigma = round(as.numeric(result),
                       digits = 8))

### graph results
## sex ratio across chromosomes
# graph histogram
data.wide.format.all.quantile %>%
  filter(SAMPLETYPE == '4' | SAMPLETYPE == 'Z'| SAMPLETYPE == '7') %>%
  filter(species =='TS') %>%
  ggplot(aes(log(ratio),
             fill = as.factor(quantile))) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed') +
  geom_vline(xintercept = -0.69,
             linetype = 'dashed') +
  geom_histogram(binwidth = 0.01)+
  xlab('Gene expression sex ratio log') +
  labs(fill = 'Quantile')+
  theme_classic() +
  ggtitle('Sex ratio gene expression evolution across quantiles' ,
          subtitle = '(Tree Swallow)') +
  facet_grid(SAMPLETYPE ~ .) 
ggsave('CAGEE/global_figures/chromsome_evolution/histogram sex ratio z vs 4_7 chromsome.png')

## sex ratio across chromosomes
# scaled z data
# graph histogram
data.wide.format.all.quantile %>%
  filter(SAMPLETYPE == '4' | SAMPLETYPE == 'Z'| SAMPLETYPE == '7') %>%
  filter(species =='TS') %>%
  rbind(data.wide.format.z.quantile.scale %>% 
  separate_wider_delim(SAMPLETYPE,
                       delim = '.',
                       names = c('SAMPLETYPE',
                                 'quantile')) %>% 
    mutate(SAMPLETYPE = 'Z.scale',
           quantile = as.integer(quantile)) %>% 
  select(GeneDescription,
         GeneName,
         SAMPLETYPE,
         quantile,
         TS) %>% 
    rename(ratio = TS)) %>% 
  ggplot(aes(log(ratio),
             fill = as.factor(quantile))) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed') +
  geom_vline(xintercept = -0.69,
             linetype = 'dashed') +
  geom_histogram(binwidth = 0.01)+
  xlab('Gene expression sex ratio log') +
  labs(fill = 'Quantile')+
  theme_classic() +
  ggtitle('Sex ratio gene expression evolution across quantiles' ,
          subtitle = '(Tree Swallow)') +
  facet_grid(SAMPLETYPE ~ .) 
ggsave('CAGEE/global_figures/chromsome_evolution/histogram sex ratio z scaled vs 4_7 chromsome.png')

## poster
p = data.wide.format.all.quantile %>%
  filter(SAMPLETYPE == '4' | SAMPLETYPE == 'Z') %>%
  filter(species =='TS') %>%
  ggplot(aes(log(ratio),
             alpha = as.factor(quantile),
             fill = SAMPLETYPE)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.69,
             linetype = 'dashed') +
  geom_vline(xintercept = -0.69,
             linetype = 'dashed') +
  geom_histogram(binwidth = 0.01)+
  xlab('log Male:Female gene expression ratio') +
  ylab('Gene Count') +
  ggtitle('Tree Swallow quantiles') +
  # labs(fill = 'Quantile')+
  theme_classic() +
  facet_grid(SAMPLETYPE ~ .)+
  xlim(-1.25,1.25)+
  theme(title = element_text(size = 28),
        axis.title.x =  element_text(size = 28),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        strip.text.y = element_text(size = 28),
        legend.position = 'null') +
  scale_fill_manual(values = c('black',
                               '#9A1B1F')) +
  scale_alpha_manual(values = c(1,
                                0.5,
                                1,
                                0.5,
                                1))
ggsave('CAGEE/global_figures/chromsome_evolution/histogram sex ratio z vs 4 chromsome poster.pdf',
       p,
       width = 7.5,
       height = 4.5)


## sigma rate across quantiles
data.cagee.results.table %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = chr,
             color = chr)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
          subtitle = '(Tree Swallow quantiles)')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile z vs 4_7 chromsome.png')

### combine with simulation data
data.cagee.results.sim.table %>% 
  filter(simulation != 'evo.normal') %>% 
  filter(simulation != 'evo.uniform') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'shifted',
                       'simulation')) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>% 
  mutate(label = ifelse(quantile == 5,
                        id,
                        NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
          subtitle = '(Tree Swallow quantiles); random simulation') +
  scale_color_manual(values = c('red',
                                'black',
                                'grey')) 
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome.png',
       height = 10,
       width = 10)

# remove uniform
data.cagee.results.sim.table %>% 
  filter(simulation != 'evo.normal') %>% 
  filter(simulation != 'evo.uniform') %>% 
  filter(simulation != 'uniform') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'shifted',
                       'simulation')) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>% 
  mutate(label = ifelse(quantile == 5,
                        id,
                        NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
          subtitle = '(Tree Swallow quantiles); random simulation') +
  scale_color_manual(values = c('red',
                                'black',
                                'grey')) 
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome outlier.png',
       height = 10,
       width = 10)

## just graph evolution simulations
data.cagee.results.sim.table %>% 
  filter(simulation != 'normal') %>% 
  filter(simulation != 'uniform') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'shifted',
                       'simulation'),
         # id = ifelse(simulation == 'shifted',
         #               paste(),
         #               id)
         ) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>% 
  mutate(label = ifelse(quantile == 5,
                        id,
                        NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type,
             shape = chr)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
          subtitle = '(Tree Swallow quantiles); evolution simulation') +
  scale_color_manual(values = c('red',
                                'black',
                                'grey')) 
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome evo.png',
       height = 10,
       width = 10)

## just graph z
data.cagee.results.sim.table %>% 
  filter(chr != 'normal') %>% 
  filter(simulation != 'uniform') %>% 
  filter(bias != 'No bias') %>% 
  filter(bias != 'Male bias') %>% 
  filter(chr == 'Z') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'shifted',
                       'simulation'),
         # id = ifelse(simulation == 'shifted',
         #               paste(),
         #               id)
  ) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>% 
  mutate(label = ifelse(quantile == 5,
                        id,
                        NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type,
             shape = chr)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2) +
  theme_classic() +
  ylab('sigma (gene expression ratio evolution)') +
  xlab('Gene expression sex ratio quantile') +
  ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
          subtitle = '(Tree Swallow quantiles); evolution simulation') +
  scale_color_manual(values = c('black',
                                'red',
                                'grey')) 
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome Z.png',
       height = 10,
       width = 10)

# poster
data.cagee.results.sim.table %>% 
  filter(bias != 'ZG') %>% 
  filter(simulation != 'normal') %>% 
  filter(simulation != 'uniform') %>% 
  filter(bias != '4') %>% 
    filter(bias != '7') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'centered',
                       'simulation'),
         # id = ifelse(simulation == 'shifted',
         #               paste(),
         #               id)
  ) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>% 
  mutate(label = id,
         label = ifelse(quantile == 5,
                        label,
                        NA),
         label = ifelse(type == 'simulation',
                        NA,
                        label),
         label = ifelse(id == 'Z shifted' & quantile == 5,
                        'Z*',
                        label)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type,
             shape = chr)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  # ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
  #         subtitle = '(Tree Swallow quantiles); evolution simulation') +
  scale_color_manual(values = c('#9A1B1F',
                                'black',
                                'grey')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome evo poster.pdf',
       height = 7.5,
       width = 7.5)

# poster
data.cagee.results.sim.table %>% 
  filter(simulation != 'normal') %>% 
  filter(simulation != 'uniform') %>% 
  filter(simulation != 'evo.normal') %>% 
  filter(simulation != 'evo.uniform') %>% 
filter(bias != 'Z.scale') %>% 
  filter(bias != '4') %>% 
  filter(bias != '7') %>% 
  mutate(id = paste(bias,
                    simulation,
                    sep = ' '),
         type = ifelse(simulation == 'shifted',
                       'centered',
                       'simulation'),
         # id = ifelse(simulation == 'shifted',
         #               paste(),
         #               id)
  ) %>% 
  select(-c(bias,
            simulation)) %>% 
  rbind(data.cagee.results.table %>% 
          mutate(id = chr,
                 type = 'real')) %>%
  mutate(label = id,
         label = ifelse(quantile == 5,
                        label,
                        NA),
         label = ifelse(type == 'simulation',
                        NA,
                        label),
         label = ifelse(id == 'Z shifted' & quantile == 5,
                        'Z*',
                        label),
         label = ifelse(id == 'ZG center' & quantile == 5,
                        'ZG*',
                        label),
         label = ifelse(id == 'ZG No bias' & quantile == 5,
                        'ZG',
                        label)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = id,
             color = type,
             shape = chr)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  # ggtitle('Rate of sex ratio gene expression evolution across sex ratio quantiles' ,
  #         subtitle = '(Tree Swallow quantiles); evolution simulation') +
  scale_color_manual(values = c('#9A1B1F',
                                'black',
                                'grey')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real vs simulated chromosome evo gametolog.png',)

### create average quantile ranking per gene across species 
# also create average ratio per gene
data.wide.format.all.quantile.avg = data.wide.format.all.quantile %>% 
  group_by(GeneDescription,
           GeneName,
           SAMPLETYPE) %>% 
  summarise(ratio.avg = mean(ratio),
            quantile.avg = mean(quantile))

## quantile
# quantile across each chromsome 
data.wide.format.all.quantile.avg %>% 
  group_by(SAMPLETYPE) %>% 
  mutate(Chr.count = n()) %>%
  filter(Chr.count >= 50) %>% 
  ggplot(aes(x = quantile.avg,
         y = reorder(SAMPLETYPE,
                   -Chr.count))) +
  ggridges::geom_density_ridges() +
  theme_classic() +
  ylab('Chromosome') +
  xlab('Average quantile score across species per gene')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile average per gene across chromosome.png')

# quantile across each chromsome 
data.wide.format.all.quantile.avg %>% 
  group_by(SAMPLETYPE) %>% 
  mutate(Chr.count = n()) %>%
  filter(Chr.count >= 50) %>% 
  filter(SAMPLETYPE %in% c('4',
                           '7',
                           'Z')) %>% 
  ggplot(aes(x = quantile.avg,
             y = reorder(SAMPLETYPE,
                         -Chr.count))) +
  ggridges::geom_density_ridges() +
  theme_classic() +
  ylab('Chromosome') +
  xlab('Average quantile score across species per gene')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile average per gene across z vs 4_7.png')

# compare Z to autosomes
data.wide.format.all.quantile.avg %>% 
  mutate(Chr.z = ifelse(SAMPLETYPE=='Z',
                        'Z chromosome',
                        'Autosome')) %>% 
  group_by(Chr.z) %>% 
  mutate(Chr.count = n()) %>%
  ggplot(aes(x = quantile.avg,
             y = reorder(Chr.z,
                         -Chr.count))) +
  ggridges::geom_density_ridges() +
  theme_classic() +
  ylab('') +
  xlab('Average quantile score across species per gene')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile average per gene Z vs Autosome.png')

### combine quantile and ratio data
## add TS quantile
data.wide.format.all.quantile.rank.avg = data.wide.format.all.quantile.avg %>% 
  right_join(data.wide.format.all.rank %>% 
  group_by(SAMPLETYPE) %>% 
  mutate(Chr.count = n())) %>% 
  left_join(data.wide.format.all.quantile) %>% 
  left_join(data.wide.format.all.quantile.ts %>% 
              select(GeneDescription,
                     GeneName,
                     SAMPLETYPE) %>% 
              separate_wider_delim(SAMPLETYPE,
                                   delim = '.',
                                   names = c('SAMPLETYPE',
                                             'TS.quantile'),
                                   too_many = 'merge'))

## graph Z chromosome synteny style graph
# quantile 1
# data.wide.format.all.quantile.rank.avg %>% 
#   filter(SAMPLETYPE == 'Z') %>% 
#   group_by(GeneDescription) %>% 
#   ggplot(aes(x = -Rank,
#              y = species,
#              color = TS.quantile == 1,
#              group = GeneDescription)) +
#   geom_point() +
#   geom_line() +
#   theme_classic() +
#   scale_color_manual(values = c('grey',
#                                 'red'))
for (i in 1:5) {
  data.wide.format.all.quantile.rank.avg %>% 
    filter(SAMPLETYPE == 'Z') %>% 
    group_by(GeneDescription) %>% 
    filter(TS.quantile != i) %>% 
    ggplot(aes(x = -Rank,
               y = species,
               group = GeneDescription)) +
    geom_point(color = 'grey') +
    geom_line(color = 'grey') +
    theme_classic() +
    geom_point(data = data.wide.format.all.quantile.rank.avg %>% 
                 filter(SAMPLETYPE == 'Z') %>% 
                 group_by(GeneDescription) %>% 
                 filter(TS.quantile == i), 
               aes(x = -Rank,
                   y = species,
                   group = GeneDescription),
               color = 'red',
               inherit.aes = F) +
    geom_line(data = data.wide.format.all.quantile.rank.avg %>% 
                filter(SAMPLETYPE == 'Z') %>% 
                group_by(GeneDescription) %>% 
                filter(TS.quantile == i), 
              aes(x = -Rank,
                  y = species,
                  group = GeneDescription),
              color = 'red',
              inherit.aes = F) +
    ggtitle(paste('TS quantile',
                  i,
                  sep = ' '))
  ggsave(paste('CAGEE/global_figures/chromsome_evolution/TS_quantile/Z/TS quantile synteny Z ',
               i,
               '.png',
               sep = ''),
         height = 10,
         width = 10)
}


## graph 4 chromosome synteny style graph
for (i in 1:5) {
  data.wide.format.all.quantile.rank.avg %>% 
    filter(SAMPLETYPE == '4') %>% 
    group_by(GeneDescription) %>% 
    filter(TS.quantile != i) %>% 
    ggplot(aes(x = -Rank,
               y = species,
               group = GeneDescription)) +
    geom_point(color = 'grey') +
    geom_line(color = 'grey') +
    theme_classic() +
    geom_point(data = data.wide.format.all.quantile.rank.avg %>% 
                 filter(SAMPLETYPE == '4') %>% 
                 group_by(GeneDescription) %>% 
                 filter(TS.quantile == i), 
               aes(x = -Rank,
                   y = species,
                   group = GeneDescription),
               color = 'red',
               inherit.aes = F) +
    geom_line(data = data.wide.format.all.quantile.rank.avg %>% 
                filter(SAMPLETYPE == '4') %>% 
                group_by(GeneDescription) %>% 
                filter(TS.quantile == i), 
              aes(x = -Rank,
                  y = species,
                  group = GeneDescription),
              color = 'red',
              inherit.aes = F) +
    ggtitle(paste('TS quantile',
                  i,
                  sep = ' '))
  ggsave(paste('CAGEE/global_figures/chromsome_evolution/TS_quantile/4/TS quantile synteny 4 ',
               i,
               '.png',
               sep = ''),
         height = 10,
         width = 10)
}





#### graph all quantile ratio data ####
### load data 
## load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

## load ratio data
data.wide.format = read_delim(
  'CAGEE/data/normalizedCounts_ratio.tsv',
  delim = '\t')

## make all chromosome gene list
data.gene.chromosome.all = data.gene.chromosome %>% 
  dplyr::select(Symbol,
                chromosome_name) %>% 
  dplyr::rename( SAMPLETYPE = chromosome_name) %>% 
  dplyr::rename(GeneName = Symbol)


## add all chromosome to other data frames
# ratio
# assign genes with no chromosome data to 'A'
data.wide.format.all = data.wide.format %>% 
  left_join(data.gene.chromosome.all) %>% 
  mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                             'A',
                             SAMPLETYPE)) %>% 
  relocate(SAMPLETYPE) %>% 
  relocate(GeneName) %>% 
  relocate(GeneDescription) 

### loop through results in 'CAGEE/ratio/all'
## get list of files
quantile.files = data.frame(files = list.files('CAGEE/ratio/all/'))
# get list of chromosome names
quantile.files = quantile.files %>% 
  separate_wider_delim(cols = 'files',
                       delim = '_',
                       names = c(NA,
                                 'Chromosome',
                                 NA,
                                 NA,
                                 NA,
                                 NA),
                       cols_remove = F)

#### create loop to load all results into one dataframe
## create empty dataframe
data.cagee.results.table.all = data.frame(chr = character(),
                                          quantile = character(),
                                          results = character(),
                                          sigma = numeric())
  
## run loop

for (i in quantile.files %>% filter(Chromosome != c('W')) %>% pull(Chromosome)) {
  ### load results
  tmp = read.delim(paste0('CAGEE/ratio/all/ratio_',
                          i,
                          '_quantile_ts_cagee_outs/results.txt'))
  
  ## convert to table
  tmp = tmp %>%
    mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                       "Sigma2: ",
                                       "")) %>%
    separate_wider_delim(CAGEE.1.1,
                         delim = ': ',
                         names = c("CAGEE.1.1",
                                   "result"),
                         too_few = 'align_end',
                         too_many = "merge") %>%
    mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                              "Attempts",
                              CAGEE.1.1))
  
  ## combine results
  # create quantile column
  data.cagee.results.table.all = tmp %>% 
    filter(grepl("Sigma",
                 CAGEE.1.1)) %>% 
    separate_wider_delim(CAGEE.1.1,
                         delim = 'Sigma for ',
                         names = c(NA,
                                   'chr')) %>% 
    separate_wider_delim(chr,
                         delim = '.',
                         names = c('chr',
                                   'quantile')) %>% 
    mutate(sigma = round(as.numeric(result),
                         digits = 8)) %>% 
    rbind(data.cagee.results.table.all)
    
    # remove temp file
    rm(tmp)
  
}

#### graph all quantile ratio
# combine with count data
p = data.cagee.results.table.all %>% 
  # filter(chr != 'A') %>% 
  left_join(data.wide.format.all %>% 
              select(SAMPLETYPE) %>% 
              group_by(SAMPLETYPE) %>% 
              summarise(count = n()) %>% 
              arrange(count) %>% 
              filter(count > 25) %>% 
              rename(chr = SAMPLETYPE)) %>% 
  mutate(color = ifelse(chr == 'Z',
                        'Z',
                        'A'),
         alpha = count/max(count)) %>% 
  mutate(label = case_when(chr == 'Z' & quantile == 5 ~ 'Z',
                           chr == 'A' & quantile == 5 ~ 'A',
                           TRUE ~ NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = chr,
             color = color,
             shape = color,
             label = label)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5,
             aes(alpha = count)) +
  geom_line(size = 2,
            aes(alpha = count)) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  scale_color_manual(values = c('black',
                                '#9A1B1F')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real chromosome.pdf',
       p,
       height = 7.5,
       width = 7.5)

# outlier
p = data.cagee.results.table.all %>% 
  filter(chr != 'A') %>%
  left_join(data.wide.format.all %>% 
              select(SAMPLETYPE) %>% 
              group_by(SAMPLETYPE) %>% 
              summarise(count = n()) %>% 
              arrange(count) %>% 
              filter(count > 25) %>% 
              rename(chr = SAMPLETYPE)) %>% 
  mutate(color = ifelse(chr == 'Z',
                        'Z',
                        'A'),
         alpha = count/max(count)) %>% 
  mutate(label = case_when(chr == 'Z' & quantile == 5 ~ 'Z',
                           # chr == 'A' & quantile == 5 ~ 'A',
                           TRUE ~ NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = chr,
             color = color,
             shape = color,
             label = label)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5,
             aes(alpha = count)) +
  geom_line(size = 2,
            aes(alpha = count)) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  scale_color_manual(values = c('black',
                                '#9A1B1F')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real chromosome outlier.pdf',
       p,
       height = 7.5,
       width = 7.5)

# outlier
p = data.cagee.results.table.all %>% 
  filter(chr != 'A') %>%
  left_join(data.wide.format.all %>% 
              select(SAMPLETYPE) %>% 
              group_by(SAMPLETYPE) %>% 
              summarise(count = n()) %>% 
              arrange(count) %>% 
              filter(count > 25) %>% 
              rename(chr = SAMPLETYPE)) %>% 
  filter(count > 200) %>% 
  mutate(color = ifelse(chr == 'Z',
                        'Z',
                        'A'),
         alpha = count/max(count)) %>% 
  mutate(label = case_when(chr == 'Z' & quantile == 5 ~ 'Z',
                           # chr == 'A' & quantile == 5 ~ 'A',
                           TRUE ~ NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = chr,
             color = color,
             shape = color,
             label = label)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5,
             aes(alpha = count)) +
  geom_line(size = 2,
            aes(alpha = count)) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  scale_color_manual(values = c('black',
                                '#9A1B1F')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real chromosome outlier threshold.pdf',
       p,
       height = 7.5,
       width = 7.5)

## add simulated data 
# need from above code: data.cagee.results.sim.table
# create reduced version
data.cagee.results.sim.table.reduce = data.cagee.results.sim.table %>% 
  filter(chr =='Z') %>% 
  filter(simulation != 'normal') %>% 
  filter(simulation != 'uniform') %>% 
  mutate(simulation = ifelse(simulation =='shifted',
                             '*',
                             simulation)) %>% 
  mutate(chr = paste0(bias,
                      simulation)) %>% 
  select(-c(bias,
            simulation)) %>% 
  mutate(count = 404) %>% 
  mutate(color = ifelse(chr == 'Z*',
                        'Z',
                        'sim'),
         alpha = 1)

# remove outlier and threshold chromosomes at 250
p = data.cagee.results.table.all %>% 
  filter(chr != 'A') %>%
  left_join(data.wide.format.all %>% 
              select(SAMPLETYPE) %>% 
              group_by(SAMPLETYPE) %>% 
              summarise(count = n()) %>% 
              arrange(count) %>% 
              filter(count > 25) %>% 
              rename(chr = SAMPLETYPE)) %>% 
  filter(count > 250) %>% 
  mutate(color = ifelse(chr == 'Z' | chr == 'Z*',
                        'Z',
                        'A')) %>% 
  mutate(alpha = count/max(count)) %>%
  full_join(data.cagee.results.sim.table.reduce) %>% 
  mutate(label = case_when(chr == 'Z' & quantile == 5 ~ 'Z',
                           chr == 'Z*' & quantile == 5 ~ 'Z*',
                           TRUE ~ NA)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = chr,
             color = color,
             shape = color,
             label = label)) +
  geom_hline(yintercept = 0.000696) +
  geom_hline(yintercept = 0.00038) +
  geom_point(size = 5,
             aes(alpha = alpha)) +
  geom_line(size = 2,
            aes(alpha = alpha)) +
  ggrepel::geom_label_repel(aes(label = label),
                            nudge_x = 0.2,
                            size = 16) +
  theme_classic() +
  ylab('sigma (ratio evolution rate)') +
  xlab('Gene expression sex ratio quantile') +
  scale_color_manual(values = c('black',
                                '#B27C2B',
                                '#9A1B1F')) +
  theme(axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text.y =  element_text(size = 16),
        axis.text.x =  element_text(size = 16),
        legend.position = 'null')
ggsave('CAGEE/global_figures/chromsome_evolution/quantile real chromosome outlier threshold evo poster.pdf',
       p,
       height = 7.5,
       width = 7.5)

#### graph results Z quantile all birds ####
### create data frame with quantile CAGEE results from each bird
## make dummy data frame
data.cagee.results.z.all.bird.table = data.frame(CAGEE.1.1 = as.character(),
                                                 results = as.character(),
                                                 species = as.character())


for (i in unique(data.wide.format.z.quantile$species)) {
  ## results.txt for each bird
  tmp = read.delim(paste0('CAGEE/ratio/ratio_species_z/',
  i,
  '/results.txt'))
  
  ## convert to table
  # z chromosome
  tmp.tbl = tmp %>%
    mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                       "Sigma2: ",
                                       "")) %>%
    separate_wider_delim(CAGEE.1.1,
                         delim = ': ',
                         names = c("CAGEE.1.1",
                                   "result"),
                         too_few = 'align_end',
                         too_many = "merge") %>%
    mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                              "Attempts",
                              CAGEE.1.1)) %>% 
    filter(grepl("Sigma",
                 CAGEE.1.1)) %>% 
    mutate(species = i)
  
  ## combine results
  # create quantile column
  data.cagee.results.z.all.bird.table = data.cagee.results.z.all.bird.table  %>% 
    rbind(tmp.tbl) 
  
  # remove
  rm(tmp)
  rm(tmp.tbl)
}



## create quantile column
data.cagee.results.z.all.bird.table = data.cagee.results.z.all.bird.table %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = 'Sigma for ',
                       names = c(NA,
                                 'chr')) %>% 
  separate_wider_delim(chr,
                       delim = '.',
                       names = c('chr',
                                 'quantile')) %>% 
  mutate(sigma = round(as.numeric(result),
                       digits = 8),
         quantile = as.integer(quantile))

### graph results
# label every species
data.cagee.results.z.all.bird.table %>% 
  mutate(species.color = ifelse(species == 'TS',
                                'TS',
                                'Other')) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = species,
             color = species)) +
  geom_line(size = 2)+
  geom_point(size = 3)+
  theme_classic() +
  scale_colour_brewer(palette = "Spectral")
ggsave('CAGEE/global_figures/chromsome_evolution/quantile/all birds quantile sigma.png')

# just color tree swallow
data.cagee.results.z.all.bird.table %>% 
  mutate(species.color = ifelse(species == 'TS',
                                'TS',
                                'Other')) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = species,
             color = species.color)) +
  geom_line(size = 2)+
  geom_point(size = 3)+
  theme_classic()+
  scale_color_manual(values = c('black',
                                '#9A1B1F'))
ggsave('CAGEE/global_figures/chromsome_evolution/quantile/all birds quantile sigma TS.png')


#### graph results chr 4 quantile all birds ####
### create data frame with quantile CAGEE results from each bird
## make dummy data frame
data.cagee.results.4.all.bird.table = data.frame(CAGEE.1.1 = as.character(),
                                                 results = as.character(),
                                                 species = as.character())


for (i in unique(data.wide.format.4.quantile$species)) {
  ## results.txt for each bird
  tmp = read.delim(paste0('CAGEE/ratio/ratio_species_4/',
                          i,
                          '/results.txt'))
  
  ## convert to table
  # 4 chromosome
  tmp.tbl = tmp %>%
    mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                       "Sigma2: ",
                                       "")) %>%
    separate_wider_delim(CAGEE.1.1,
                         delim = ': ',
                         names = c("CAGEE.1.1",
                                   "result"),
                         too_few = 'align_end',
                         too_many = "merge") %>%
    mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                              "Attempts",
                              CAGEE.1.1)) %>% 
    filter(grepl("Sigma",
                 CAGEE.1.1)) %>% 
    mutate(species = i)
  
  ## combine results
  # create quantile column
  data.cagee.results.4.all.bird.table = data.cagee.results.4.all.bird.table  %>% 
    rbind(tmp.tbl) 
  
  # remove
  rm(tmp)
  rm(tmp.tbl)
}



## create quantile column
data.cagee.results.4.all.bird.table = data.cagee.results.4.all.bird.table %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = 'Sigma for ',
                       names = c(NA,
                                 'chr')) %>% 
  separate_wider_delim(chr,
                       delim = '.',
                       names = c('chr',
                                 'quantile')) %>% 
  mutate(sigma = round(as.numeric(result),
                       digits = 8),
         quantile = as.integer(quantile))

### graph results
# label every species
data.cagee.results.4.all.bird.table %>% 
  mutate(species.color = ifelse(species == 'TS',
                                'TS',
                                'Other')) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = species,
             color = species)) +
  geom_line(size = 2)+
  geom_point(size = 3)+
  theme_classic() +
  scale_colour_brewer(palette = "Spectral")
ggsave('CAGEE/global_figures/chromsome_evolution/quantile/all birds quantile sigma chr 4.png')

# just color tree swallow
data.cagee.results.4.all.bird.table %>% 
  mutate(species.color = ifelse(species == 'TS',
                                'TS',
                                'Other')) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = species,
             color = species.color)) +
  geom_line(size = 2)+
  geom_point(size = 3)+
  theme_classic()+
  scale_color_manual(values = c('black',
                                '#9A1B1F'))
ggsave('CAGEE/global_figures/chromsome_evolution/quantile/all birds quantile sigma chr 4 TS.png')

### compare chromsome 4 to Z
# just color tree swallow
data.cagee.results.z.all.bird.table %>% 
  mutate(chr = 'Z') %>% 
  rbind(data.cagee.results.4.all.bird.table %>% 
          mutate(chr = '4')) %>% 
  mutate(species.color = ifelse(species == 'TS',
                                'TS',
                                'Other'),
         species.color = paste(species.color,
                               chr),
         ID = paste(species,
                    chr)) %>% 
  ggplot(aes(x = quantile,
             y = sigma,
             group = ID,
             color = species.color)) +
  geom_line(size = 2)+
  geom_point(size = 3)+
  theme_classic()+
  scale_color_manual(values = c('grey',
                                'black',
                                '#9A1B1F',
                                '#9A1B1F'))
ggsave('CAGEE/global_figures/chromsome_evolution/quantile/all birds quantile sigma chr 4 vs Z TS.png')


#### graph ratios histogram ####
### load data
data.wide.format = read.delim('CAGEE/normalizedCounts_ratio.tsv')

## convert to long format
data.long = data.wide.format %>% 
  pivot_longer(cols = -c(GeneDescription,
                         GeneName),
               names_to = 'species',
               values_to = 'sex_ratio')
# add log 
data.long = data.long %>% 
  mutate(sex_ratio_log = log(sex_ratio))

## graph histogram
# all
data.long %>% 
  ggplot(aes(sex_ratio_log)) +
  geom_histogram(bins = 200) +
  theme_classic()
ggsave('CAGEE/histogram log sex ratio.png')

# filter middle
data.long %>% 
  mutate(keep = case_when(sex_ratio_log > 1 ~ 1,
                          sex_ratio_log < -1 ~ 1)) %>% 
  filter(keep == 1) %>% 
  ggplot(aes(sex_ratio_log)) +
  geom_histogram(bins = 200) +
  theme_classic()
ggsave('CAGEE/histogram log sex ratio filter.png')

#### simulation results for n_gamma_cats ####
#### load and combine data for gene expression and ratio simulations
### gene expression
# make unique transcript ID
ngamma_sim_gene_exp = read_tsv('CAGEE/n_gamma_simulation/median_100_0.005/simulation.txt',
                               skip = 4) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_0.01/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_0.015/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_0.02/simulation.txt',
                 skip = 4)) %>% 
  mutate(GENE_ID = paste0(GENE_ID,
                          '_',
                          DESC))
## save file
write_tsv(ngamma_sim_gene_exp,
          'CAGEE/n_gamma_simulation/median_100_combined.tsv')

## large sigma
# make unique transcript ID
ngamma_sim_gene_exp_large = read_tsv('CAGEE/n_gamma_simulation/median_100_0.001/simulation.txt',
                               skip = 4) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_0.01/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_0.1/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/median_100_1/simulation.txt',
                 skip = 4)) %>% 
  mutate(GENE_ID = paste0(GENE_ID,
                          '_',
                          DESC))
## save file
write_tsv(ngamma_sim_gene_exp_large,
          'CAGEE/n_gamma_simulation/median_100_combined_large.tsv')

### ratio
# make unique transcript ID
ngamma_sim_gene_ratio = read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.0003/simulation.txt',
                               skip = 4) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.0006/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.0009/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.0012/simulation.txt',
                 skip = 4)) %>% 
  mutate(GENE_ID = paste0(GENE_ID,
                          '_',
                          DESC))
## save file
write_tsv(ngamma_sim_gene_exp,
          'CAGEE/n_gamma_simulation/ratio_100_combined.tsv')

## large sigma
# make unique transcript ID
ngamma_sim_gene_ratio_large = read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.0001/simulation.txt',
                                 skip = 4) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.001/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.01/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/ratio_100_0.1/simulation.txt',
                 skip = 4)) %>% 
  mutate(GENE_ID = paste0(GENE_ID,
                          '_',
                          DESC))
## save file
write_tsv(ngamma_sim_gene_ratio_large,
          'CAGEE/n_gamma_simulation/ratio_100_combined_large.tsv')



#### load n_gamma_cats data ####
# remove NA
# pivot long
# only keep category with the top liklihood score per gene

#### 4 categories
### simulation
## median
# 4 categories 
ngamma_median_4_sim = read_tsv('CAGEE/n_gamma_simulation/median_sim_100_cagee_outs_n_gamma_cats//category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'sim',
         gamma.count = 4) 

# large sigma
ngamma_median_4_sim_large = read_tsv('CAGEE/n_gamma_simulation/median_sim_100_large_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'sim.large',
         gamma.count = 4) 

## ratio
# 4 categories 
ngamma_ratio_4_sim = read_tsv('CAGEE/n_gamma_simulation/ratio_sim_100_cagee_outs_n_gamma_cats//category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'sim',
         gamma.count = 4) 

# large sigma
ngamma_ratio_4_sim_large = read_tsv('CAGEE/n_gamma_simulation/ratio_sim_100_large_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'sim.large',
         gamma.count = 4) 

### transcriptome
## median
# 4 categories 
ngamma_median_4 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 4) 

## ratio
# 4 categories 
ngamma_ratio_4 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 4)

### Z chromosome
## median
# 4 categories 
ngamma_median_z_4 = read_tsv('CAGEE/median_z/median_z_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'Z',
         gamma.count = 4)


## ratio
# 4 categories 
ngamma_ratio_z_4 = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats/category_likelihoods.txt')%>% 
  select(-c('...6')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'Z',
         gamma.count = 4)

#### 10 categories
### transcriptome
## median
# 10 categories 
ngamma_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 

## ratio
# 10 categories 
ngamma_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10)

### Z chromosome
## median
# 10 categories 
ngamma_median_z_10 = read_tsv('CAGEE/median_z/median_z_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'Z',
         gamma.count = 10)


## ratio
# 10 categories 
ngamma_ratio_z_10 = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats_10/category_likelihoods.txt')%>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = max(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'Z',
         gamma.count = 10)

### combine into one data frame
ngamma.df = rbind(ngamma_median_4,
                  ngamma_ratio_4) %>% 
  rbind(ngamma_median_z_4) %>% 
  rbind(ngamma_ratio_z_4) %>% 
  rbind(ngamma_median_10) %>% 
  rbind(ngamma_ratio_10) %>% 
  rbind(ngamma_median_z_10) %>% 
  rbind(ngamma_ratio_z_10) %>% 
  rbind(ngamma_median_4_sim) %>% 
  rbind(ngamma_ratio_4_sim) %>% 
  rbind(ngamma_ratio_4_sim_large) %>% 
  rbind(ngamma_median_4_sim_large)

## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma.df.table = ngamma.df %>% 
  select(-c(transcript,
            liklihood,
            keep,
            top.cat)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat = as.numeric(as.character(gamma.cat))) %>% 
  filter(Freq > 0) %>% 
  group_by(type,
           gene.list,
           gamma.count) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(percent = 100*Freq/total)

#### graph n_gamma_cats data ####
### graph genes per category
ngamma.df.table %>% 
  ggplot(aes(x = gamma.cat,
             y = percent,
             group = gene.list,
             color = gene.list)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  facet_grid(gamma.count ~ type)
ggsave('CAGEE/global_figures/n_gamma/Percent of genes per sigma category.png')

## graph gene categories and gamma.cat
ngamma_ratio_4_sim %>% 
  rbind(ngamma_median_4_sim) %>% 
  rbind(ngamma_median_4_sim_large) %>% 
  rbind(ngamma_ratio_4_sim_large) %>% 
  separate_wider_delim(cols = transcript,
                       delim = 'SIG',
                       names = c(NA,
                                 'sigma'),
                       cols_remove = F) %>% 
  mutate(sigma = as.numeric(sigma)) %>% 
  select(sigma, 
         gamma.cat,
         type,
         gene.list) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>%
  mutate(gamma.cat = as.numeric(as.character(gamma.cat)),
         gamma.cat.type = ifelse(gamma.cat >= 1,
                                 'high',
                                 'low')) %>% 
  ggplot(aes(x = sigma,
             y = Freq,
             fill = gamma.cat.type)) + 
  geom_bar(position="stack", 
           stat="identity") +
  theme_classic() +
  facet_grid(gene.list~type,
             scales = 'free') +
  xlab('simulated sigma') 
ggsave('CAGEE/global_figures/n_gamma/Simulation percent of genes per sigma category.png')

#### load n_gamma_cats data low ####
### assign to lowest liklihood

# remove NA
# pivot long
# only keep category with the top liklihood score per gene

#### 4 categories
### simulation
## median
# 4 categories 
ngamma_low_median_4_sim = read_tsv('CAGEE/n_gamma_simulation/median_sim_100_cagee_outs_n_gamma_cats//category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'sim',
         gamma.count = 4) 

# large sigma
ngamma_low_median_4_sim_large = read_tsv('CAGEE/n_gamma_simulation/median_sim_100_large_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'sim.large',
         gamma.count = 4) 

## ratio
# 4 categories 
ngamma_low_ratio_4_sim = read_tsv('CAGEE/n_gamma_simulation/ratio_sim_100_cagee_outs_n_gamma_cats//category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'sim',
         gamma.count = 4) 

# large sigma
ngamma_low_ratio_4_sim_large = read_tsv('CAGEE/n_gamma_simulation/ratio_sim_100_large_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'sim.large',
         gamma.count = 4) 

### transcriptome
## median
# 4 categories 
ngamma_low_median_4 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 4) 

## ratio
# 4 categories 
ngamma_low_ratio_4 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 4)

### Z chromosome
## median
# 4 categories 
ngamma_low_median_z_4 = read_tsv('CAGEE/median_z/median_z_cagee_outs_n_gamma_cats/category_likelihoods.txt') %>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'Z',
         gamma.count = 4)


## ratio
# 4 categories 
ngamma_low_ratio_z_4 = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats/category_likelihoods.txt')%>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'Z',
         gamma.count = 4)

#### 10 categories
### transcriptome
## median
# 10 categories 
ngamma_low_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  dplyr::select(-c('...12')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 

## ratio
# 10 categories 
ngamma_low_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  dplyr::select(-c('...12')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10)

### Z chromosome
## median
# 10 categories 
ngamma_low_median_z_10 = read_tsv('CAGEE/median_z/median_z_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  dplyr::select(-c('...12')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'Z',
         gamma.count = 10)


## ratio
# 10 categories 
ngamma_low_ratio_z_10 = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats_10/category_likelihoods.txt')%>% 
  dplyr::select(-c('...12')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'Z',
         gamma.count = 10)

### combine into one data frame
ngamma_low.df = rbind(ngamma_low_median_4,
                  ngamma_low_ratio_4) %>% 
  rbind(ngamma_low_median_z_4) %>% 
  rbind(ngamma_low_ratio_z_4) %>% 
  rbind(ngamma_low_median_10) %>% 
  rbind(ngamma_low_ratio_10) %>% 
  rbind(ngamma_low_median_z_10) %>% 
  rbind(ngamma_low_ratio_z_10) %>% 
  rbind(ngamma_low_median_4_sim) %>% 
  rbind(ngamma_low_ratio_4_sim) %>% 
  rbind(ngamma_low_ratio_4_sim_large) %>% 
  rbind(ngamma_low_median_4_sim_large)

## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma_low.df.table = ngamma_low.df %>% 
  dplyr::select(-c(transcript,
            liklihood,
            keep,
            top.cat)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat = as.numeric(as.character(gamma.cat))) %>% 
  filter(Freq > 0) %>% 
  group_by(type,
           gene.list,
           gamma.count) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(percent = 100*Freq/total)

#### graph n_gamma_cats data low ####
### graph genes per category
# ngamma_low.df.table %>% 
#   ggplot(aes(x = gamma.cat,
#              y = percent,
#              group = gene.list,
#              color = gene.list)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   facet_grid(gamma.count ~ type)
# ggsave('CAGEE/global_figures/n_gamma/Percent of genes per sigma category low.png')

## graph gene categories and gamma.cat
# ngamma_low_ratio_4_sim %>% 
#   rbind(ngamma_low_median_4_sim) %>% 
#   rbind(ngamma_low_median_4_sim_large) %>% 
#   rbind(ngamma_low_ratio_4_sim_large) %>% 
#   separate_wider_delim(cols = transcript,
#                        delim = 'SIG',
#                        names = c(NA,
#                                  'sigma'),
#                        cols_remove = F) %>% 
#   mutate(sigma = as.numeric(sigma)) %>% 
#   select(sigma, 
#          gamma.cat,
#          type,
#          gene.list) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq != 0) %>% 
#   group_by(gene.list,
#             type) %>% 
#   mutate(gamma.cat.type = dense_rank(gamma.cat),
#          sigma.type = dense_rank(sigma)) %>%
#   ggplot(aes(x = sigma.type,
#              y = Freq,
#              fill = gamma.cat.type)) + 
#   geom_bar(position="stack", 
#            stat="identity") +
#   theme_classic() +
#   facet_grid(gene.list~type,
#              scales = 'free') +
#   xlab('simulated sigma') 
# ggsave('CAGEE/global_figures/n_gamma/Simulation percent of genes per sigma category_low.png')

# paper
ngamma_low_ratio_4_sim %>% 
  rbind(ngamma_low_median_4_sim) %>% 
  rbind(ngamma_low_median_4_sim_large) %>% 
  rbind(ngamma_low_ratio_4_sim_large) %>% 
  separate_wider_delim(cols = transcript,
                       delim = 'SIG',
                       names = c(NA,
                                 'sigma'),
                       cols_remove = F) %>% 
  mutate(sigma = as.numeric(sigma)) %>% 
  select(sigma, 
         gamma.cat,
         type,
         gene.list) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>% 
  group_by(gene.list,
           type) %>% 
  mutate(gamma.cat.type = dense_rank(gamma.cat),
         sigma.type = dense_rank(sigma)) %>%
  ggplot(aes(x = sigma.type,
             y = Freq,
             fill = gamma.cat.type)) + 
  geom_bar(position="stack", 
           stat="identity") +
  theme_classic(base_size = 8) +
  facet_grid(gene.list~type,
             scales = 'free') +
  xlab('simulated sigma') +
  # theme(legend.key.size = unit(.2, 'in')) +
  theme(legend.position = 'none')
ggsave('CAGEE/global_figures/n_gamma/Simulation percent of genes per sigma category_low.pdf',
       height = 3.25,
       width = 3.25,
       units = 'in',
       dpi = 320)

## compare known sigma to calculated rates
# pull sigma from results files
# ngamma_low.df.table %>% 
#   filter(gene.list %in% c('sim',
#                           'sim.large')) %>% 
#   select(gamma.cat, 
#          type,
#          gene.list) %>% 
#   distinct() %>%  
#   mutate(rate = case_when(
#     gene.list == 'sim' & type == 'median' ~ 0.012539,
#     gene.list == 'sim' & type == 'ratio' ~ 0.000926,
#     gene.list == 'sim.large' & type == 'median' ~ 0.250032,
#     gene.list == 'sim.large' & type == 'ratio' ~ 0.044866
#   ),
#   rate.gamma = gamma.cat*rate) %>% 
#   cbind(data.frame(rate.sim = c(0.005, 0.010, 0.015, 0.020,
#                                 0.0003, 0.0006, 0.0009, 0.0012,
#                                 0.001, 0.010, 0.1, 1,
#                                 0.0001, 0.001, 0.01, 0.1))) %>% 
#   ggplot(aes(y = rate.gamma,
#              x = rate.sim,
#              group = gene.list,
#              color = gene.list)) +
#   geom_abline(slope = 1,
#               intercept = 0)+
#     geom_point() +
#     geom_line() + 
#     facet_wrap(~gene.list + type,
#           scales = 'free') +
#   theme_classic()
# ggsave('CAGEE/global_figures/n_gamma/Simulation rate vs gamma rate low.png')

# paper
ngamma_low.df.table %>% 
  filter(gene.list %in% c('sim',
                          'sim.large')) %>% 
  select(gamma.cat, 
         type,
         gene.list) %>% 
  distinct() %>%  
  mutate(rate = case_when(
    gene.list == 'sim' & type == 'median' ~ 0.012539,
    gene.list == 'sim' & type == 'ratio' ~ 0.000926,
    gene.list == 'sim.large' & type == 'median' ~ 0.250032,
    gene.list == 'sim.large' & type == 'ratio' ~ 0.044866
  ),
  rate.gamma = gamma.cat*rate) %>% 
  cbind(data.frame(rate.sim = c(0.005, 0.010, 0.015, 0.020,
                                0.0003, 0.0006, 0.0009, 0.0012,
                                0.001, 0.010, 0.1, 1,
                                0.0001, 0.001, 0.01, 0.1))) %>% 
  ggplot(aes(y = rate.gamma,
             x = rate.sim,
             group = gene.list)) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = 'dashed',
              color = 'darkgrey')+
  geom_point() +
  geom_line() + 
  facet_wrap(~gene.list + type,
             scales = 'free') +
  theme_classic(base_size = 8) +
  theme(legend.position = 'none')
ggsave('CAGEE/global_figures/n_gamma/Simulation rate vs gamma rate low.pdf',
       height = 3.25,
       width = 3.25,
       units = 'in',
       dpi = 320)

#### compare gamma cats with ancestral state ####
#### z
### load data
## ratio
# # 10 categories 
# ngamma_low_ratio_z_10 = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats_10/category_likelihoods.txt')%>% 
#   dplyr::select(-c('...12')) %>% 
#   rename(transcript = "Transcript ID") %>% 
#   pivot_longer(cols = -c(transcript),
#                values_to = "liklihood",
#                names_to = 'gamma.cat') %>% 
#   mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
#   group_by(transcript) %>%
#   mutate(top.cat = min(liklihood)) %>% 
#   ungroup() %>% 
#   mutate(keep = ifelse(top.cat == liklihood,
#                        1,
#                        0)) %>% 
#   filter(keep == 1) %>% 
#   mutate(type = 'ratio',
#          gene.list = 'Z',
#          gamma.count = 10)
# 
# # get sigma 
# ngamma_low_ratio_z_10.sigma = read.delim('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats_10/results.txt') %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>% 
#   filter(CAGEE.1.1 == 'Sigma2') %>% 
#   separate_wider_delim(result,
#                        names = c('result',
#                                  'value'),
#                        delim = ': ') %>% 
#   pull(value) %>% 
#   as.numeric()
# 
# ## check tree for root node
# ape::read.tree('CAGEE/data/10_species_birds_edit.nwk') %>% 
# plot()
# edgelabels(cex = 0.8) 
# nodelabels(cex = 0.8)
# dev.off()
# 
# ## ancestral state
# # remove autosomes 
# #select root node: 11
# ratio.ancestral.state.df = read_tsv('CAGEE/ratio/ratio_z_cagee_outs/credible_intervals.tab') %>% 
#   filter(SAMPLETYPE == 'Z') %>% 
#   select(TranscriptID,
#          '<11>') %>% 
#   rename(root = '<11>') %>% 
#   mutate(root = substr(root, 2, nchar(root))) %>% 
#   mutate(root = substr(root, 0, nchar(root)-1)) %>% 
#   separate_wider_delim(cols = root,
#                        delim = '-',
#                        names = c('min.root',
#                                  'max.root')) %>% 
#   mutate(min.root.log = log(as.numeric(min.root)),
#          max.root.log = log(as.numeric(max.root)))
# 
# ## graph distribution of min and max root ratio values
# ratio.ancestral.state.df %>% 
#   
# 
# # min
# hist(ratio.ancestral.state.df$min.root.log)
# # max
# hist(ratio.ancestral.state.df$max.root.log)

### load data
## ratio
ngamma_low_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  dplyr::select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10)

# get sigma 
ngamma_low_ratio_10.sigma = read.delim('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

## get ancestral state at root 
ratio.z.ancestral.state.df = read_tsv('CAGEE/ratio/ratio_z_cagee_outs_n_gamma_cats_10/ancestral_states.tab') %>% 
  filter(SAMPLETYPE == 'Z') %>% 
  dplyr::select(TranscriptID,
           '<11>') %>% 
  rename(root = '<11>') %>% 
  rename(transcript = TranscriptID)

## median
ngamma_low_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 

# get sigma 
ngamma_low_median_10.sigma = read.delim('CAGEE/median/median_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

# #### graph
# ### compare ancestral state values with sigma categories
# # use log of root value
# ngamma_low_ratio_10 %>% 
#   left_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = log.root,
#              group = gamma.cat)) +
#   # geom_hline(yintercept = log(2),
#   #            linetype = 'dashed') +
#   geom_hline(yintercept = 0.4252187,
#                         linetype = 'dashed') +
#   geom_hline(yintercept = 0) +
#   geom_violin() +
#   geom_boxplot() +
#   stat_summary(fun.y=mean, 
#                geom="point", 
#                shape=20, 
#                size=2, 
#                color="red", 
#                fill="red") +
#   theme_classic() +
#   geom_text(data = ngamma_low_ratio_10 %>% 
#               group_by(gamma.cat) %>% 
#               summarise(count = n()),
#             aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#                 label = paste0(count),
#                 y = 0.85)) +
#   xlab('Sigma category') +
#   ylab('Sex ratio ancestral state (log)') +
#   ggtitle('Z chromosome sigma categories vs ancestral state') +
#   scale_y_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome sigma categories vs ancestral state.png')


### statistics
library(nptest)
# Test whether each group differs from average
# could add * to boxplots

# # t-test
# ngamma_low_ratio_z_10.pvalue = ngamma_low_ratio_10 %>%
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>%
#   group_by(gamma.cat) %>%
#   summarise(P = t.test(log.root, mu = 0.4252187)$p.value,
#             Sig = ifelse(P < 0.05, "*", NA)) %>% 
#   mutate(p.adj = p.adjust(P,
#                           n = 10),
#          Sig.adj = ifelse(p.adj < 0.05, "*", NA))
# # doesn't work as data isn't normal...
# 
# # wilcox test
# ngamma_low_ratio_10 %>%
#   right_join(ratio.z.ancestral.state.df)  %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>%
# group_by(gamma.cat) %>%
#   summarise(P = wilcox.test(log.root, mu = 0.4252187,
#                             correct = F)$p.value,
#             Sig = ifelse(P < 0.05, "*", "ns")) %>%
#   mutate(p.adj = p.adjust(P,
#                           n = 10),
#          Sig.adj = ifelse(p.adj < 0.05, "*", "ns"))
# # doesn't work as there are ties in data...

## wilcoxon signed rank test with permutation
# ratio
ngamma_low_ratio_z_10.pvalue = ngamma_low_ratio_10 %>% 
  right_join(ratio.z.ancestral.state.df)  %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>%
  group_by(gamma.cat) %>%
  summarise(P = np.loc.test(log.root, 
                            mu = mean.log.root,
                            R = 10000,
                            parallel = T,
                            median.test = T,
                            symmetric = T)$p.value,
            Sig = ifelse(P < 0.05, "*", NA)) %>% 
  mutate(p.adj = p.adjust(P,
                          method = 'fdr',
                          n = 10),
         Sig.adj = ifelse(p.adj < 0.05, "*", NA))

# median
ngamma_low_median_z_10.pvalue = ngamma_low_median_10 %>% 
  right_join(ratio.z.ancestral.state.df)  %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>%
  group_by(gamma.cat) %>%
  summarise(P = np.loc.test(log.root, 
                            mu = mean.log.root,
                            R = 10000,
                            parallel = T,
                            median.test = T,
                            symmetric = T)$p.value,
            Sig = ifelse(P < 0.05, "*", NA)) %>% 
  mutate(p.adj = p.adjust(P,
                          method = 'fdr',
                          n = 10),
         Sig.adj = ifelse(p.adj < 0.05, "*", NA))


## ordinal regression
# ratio
ngamma_low_ratio_z_10.pvalue.ord = ordinal::clm(gamma.cat.value ~ log.root,
                    data = ngamma_low_ratio_10 %>% 
                      right_join(ratio.z.ancestral.state.df)  %>%
                      mutate(log.root = log(root),
                             gamma.cat.value = gamma.cat*ngamma_low_ratio_10.sigma,
                             gamma.cat.value = factor(gamma.cat.value, order= TRUE)
                             )) 

summary(ngamma_low_ratio_z_10.pvalue.ord)

# median
ngamma_low_median_z_10.pvalue.ord = ordinal::clm(gamma.cat.value ~ log.root,
                                                data = ngamma_low_median_10 %>% 
                                                  right_join(ratio.z.ancestral.state.df)  %>%
                                                  mutate(log.root = log(root),
                                                         gamma.cat.value = gamma.cat*ngamma_low_median_10.sigma,
                                                         gamma.cat.value = factor(gamma.cat.value, order= TRUE)
                                                  )) 

summary(ngamma_low_median_z_10.pvalue.ord)





# # graph with significance
# ngamma_low_ratio_z_10 %>% 
#   left_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#              y = log.root,
#              group = gamma.cat)) +
#   # geom_hline(yintercept = log(2),
#   #            linetype = 'dashed') +
#   geom_hline(yintercept = 0.4252187,
#              linetype = 'dashed') +
#   geom_hline(yintercept = 0) +
#   geom_violin() +
#   geom_boxplot() +
#   geom_text(data = ngamma_low_ratio_z_10.pvalue,
#             aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#                 label = Sig.adj,
#                 y = 0.85)) +
#   stat_summary(fun.y=mean, 
#                geom="point", 
#                shape=20, 
#                size=2, 
#                color="red", 
#                fill="red") +
#   theme_classic() +
#   geom_text(data = ngamma_low_ratio_z_10 %>% 
#               group_by(gamma.cat) %>% 
#               summarise(count = n()),
#             aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#                 label = paste0(count),
#                 y = -0.4)) +
#   xlab('Sigma category') +
#   ylab('Sex ratio ancestral state (log)') +
#   ggtitle('Z chromosome sigma categories vs ancestral state') +
#   scale_y_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome sigma categories vs ancestral state sig.png')
# 
# ## for presentation
# # graph with significance
# ngamma_low_ratio_z_10 %>% 
#   left_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#              y = log.root,
#              group = gamma.cat)) +
#   # geom_hline(yintercept = log(2),
#   #            linetype = 'dashed') +
#   geom_hline(yintercept = 0) +
#   geom_violin() +
#   # geom_boxplot() +
#   geom_text(data = ngamma_low_ratio_z_10.pvalue,
#             aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#                 label = Sig.adj,
#                 y = 0.85)) +
#   stat_summary(fun.y=mean, 
#                geom="point", 
#                shape=20, 
#                size=2, 
#                color="red", 
#                fill="red") +
#   theme_classic() +
#   geom_text(data = ngamma_low_ratio_z_10 %>% 
#               group_by(gamma.cat) %>% 
#               summarise(count = n()),
#             aes(x = gamma.cat*ngamma_low_ratio_z_10.sigma,
#                 label = paste0(count),
#                 y = -0.4)) +
#   geom_hline(yintercept = 0.4252187,
#              linetype = 'dashed') +
#   xlab('Sigma category') +
#   ylab('Sex ratio ancestral state (log)') +
#   ggtitle('Z chromosome sigma categories vs ancestral state') +
#   scale_y_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome sigma categories vs ancestral state sig presentation.png')

## for paper
## graph with significance
## ratio
# ngamma_low_ratio_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = log.root,
#              group = gamma.cat)) +
#   geom_hline(yintercept = 0) +
#   geom_violin() + 
#   stat_summary(fun.y=mean, 
#                geom="point", 
#                shape=18,
#                size=3, 
#                color="red") +
#   geom_text(data = ngamma_low_ratio_z_10.pvalue,
#             aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#                 label = Sig.adj,
#                 y = 0.85)) +
#   theme_classic() +
#   geom_hline(yintercept = 0.4252187,
#              linetype = 'dashed') +
#   xlab('Ratio sigma category') +
#   ylab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome ratio sigma categories vs ancestral state') +
#   scale_y_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ratio sigma categories vs ancestral state sig paper.png')

# flip axis
# boxplot
ngamma_low_ratio_10 %>% 
  right_join(ratio.z.ancestral.state.df) %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = 0.4252187,
             linetype = 'dashed') +
  geom_hline(yintercept = 0) +
  geom_boxplot(position = position_dodge2(preserve = "single"),
               width = 0.00003,
               fill = 'lightgrey') +
  theme_classic() +
  xlab('Ratio sigma category') +
  ylab('Sex ratio ancestral state (ln)') +
  ggtitle('Z chromosome ratio sigma categories vs ancestral state') +
  scale_y_continuous(breaks = c(-0.4,
                                0,
                                0.4,
                                0.8))+
  coord_flip()
ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ratio sigma categories vs ancestral state sig flip paper.png')


# paper
ngamma_low_ratio_10 %>% 
  right_join(ratio.z.ancestral.state.df) %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(position = position_dodge2(preserve = "single"),
               width = 0.00003,
               fill = 'red') +
  geom_hline(yintercept = 0.4252187,
             linetype = 'dashed') +
  theme_classic(base_size = 8) +
  xlab('Ratio sigma category') +
  ylab('Sex ratio ancestral state (ln)') +
  scale_y_continuous(breaks = c(-0.4,
                                0,
                                0.4,
                                0.8))+
  coord_flip()
ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ratio sigma categories vs ancestral state sig flip paper.pdf',
       height = 4,
       width = 6.5,
       units = 'in',
       dpi = 320)

# paper
# ordered
ngamma_low_ratio_10 %>% 
  right_join(ratio.z.ancestral.state.df) %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root),
         sigma_value =formatC(gamma.cat*ngamma_low_ratio_10.sigma, 
                              format = "e", 
                              digits = 2)) %>% 
  ggplot(aes(x = reorder(sigma_value,
                         gamma.cat),
             y = log.root,
             group = gamma.cat*ngamma_low_ratio_10.sigma
             )) +
  geom_hline(yintercept = 0) +
  geom_boxplot(position = position_dodge2(preserve = "single"),
               fill = 'red') +
  geom_hline(yintercept = 0.4252187,
             linetype = 'dashed') +
  theme_classic(base_size = 8) +
  xlab('Ratio sigma category') +
  ylab('Sex ratio ancestral state (ln)') +
  scale_y_continuous(breaks = c(-0.4,
                                0,
                                0.4,
                                0.8)) +
  # scale_x_discrete(breaks = c(1e-05,
  #                               0.00041,
  #                               0.00183)) +
  coord_flip()
ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ratio sigma categories vs ancestral state sig flip order paper.pdf',
       height = 4,
       width = 6.5,
       units = 'in',
       dpi = 320)

# 
# # flip axis
# # open circle jitter
# ngamma_low_ratio_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_ratio_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0,
#              color = 'grey50') +
#   geom_jitter(width = 0.006,
#               height = 0.0004,
#               shape = 21,
#               color = 'red') +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed',
#              color = 'grey50') +
#   ylab('Ratio sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome ratio ancestral state vs ratio sigma category') +
#   coord_fixed(0.8/0.0020)
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs ratio sigma category open.png')
# 
# # 
# # # add area
# # ngamma_low_ratio_10 %>% 
# #   right_join(ratio.z.ancestral.state.df) %>%
# #   mutate(log.root = log(root),
# #          mean.log.root = mean(log.root)) %>% 
# #   ggplot(aes(y = gamma.cat*ngamma_low_ratio_10.sigma,
# #              x = log.root,
# #              group = gamma.cat)) +
# #   geom_vline(xintercept = 0) +
# #   geom_count(aes(color = after_stat(n), 
# #                  size = after_stat(n))) +
# #   guides(color = 'legend') +
# #   # scale_x_binned(n.breaks = 10) +
# #   theme_classic() +
# #   geom_vline(xintercept = 0.4252187,
# #              linetype = 'dashed') +
# #   ylab('Ratio sigma category') +
# #   xlab('Sex ratio ancestral state (ln)') +
# #   ggtitle('Z chromosome ratio sigma categories vs ancestral state') +
# #   # scale_x_continuous(breaks = c(-0.4,
# #   #                               0,
# #   #                               0.4,
# #   #                               0.8)) +
# #   scale_color_gradient(low = 'red1',
# #                        high = 'red4')
# # ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs ratio sigma category area.png')
# # 
# 
# ## median
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_median_10.sigma,
#              y = log.root,
#              group = gamma.cat)) +
#   geom_hline(yintercept = 0) +
#   geom_violin() + 
#   stat_summary(fun.y=mean, 
#                geom="point", 
#                shape=18,
#                size=3, 
#                color="red") +
#   geom_text(data = ngamma_low_median_z_10.pvalue,
#             aes(x = gamma.cat*ngamma_low_median_10.sigma,
#                 label = Sig.adj,
#                 y = 0.85)) +
#   theme_classic() +
#   geom_hline(yintercept = 0.4252187,
#              linetype = 'dashed') +
#   xlab('Median sigma category') +
#   ylab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome median sigma categories vs ancestral state') +
#   # coord_fixed(0.04/0.8) +
#   scale_y_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome sigma categories vs ancestral state sig paper.png')
# 
# # flip axis
# # open circle jitter
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_median_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0,
#              color = 'grey50') +
#   geom_jitter(width = 0.006,
#               height = 0.0004,
#               shape = 21,
#               color = 'red') +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed',
#              color = 'grey50') +
#   ylab('Median sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome ratio ancestral state vs median sigma category') +
#   scale_x_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8)) +
#   coord_fixed(0.8/0.04)
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs median sigma category open.png')
# # 
# # open circle jitter alpha
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_median_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0,
#              color = 'grey50') +
#   geom_jitter(width = 0.006,
#               height = 0.0004,
#               shape = 21,
#               color = 'red',
#               fill = "red",
#               alpha = 0.5) +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed',
#              color = 'grey50') +
#   ylab('Median sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome median sigma categories vs ancestral state') +
#   scale_x_continuous(breaks = c(-0.4,
#                                 0,
#                                 0.4,
#                                 0.8))
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs median sigma category open alpha.png')
# 
# 
#  # add area
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_median_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0) +
#   geom_count(aes(color = after_stat(n), 
#                  size = after_stat(n))) +
#   guides(color = 'legend') +
#   # scale_x_binned(n.breaks = 10) +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed') +
#   ylab('Median sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome median sigma categories vs ancestral state') +
#   # scale_x_continuous(breaks = c(-0.4,
#   #                               0,
#   #                               0.4,
#   #                               0.8)) +
#   scale_color_gradient(low = 'red1',
#                        high = 'red4')
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs median sigma category area.png')
# 
# # add alpha
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_median_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0) +
#   geom_count(aes(alpha = after_stat(n),
#              color = 'red')) +
#   scale_alpha(range = c(0.5,1)) +
#   guides(color = 'legend') +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed') +
#   ylab('Median sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome median sigma categories vs ancestral state') 
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs median sigma category area alpha.png')
# 
# # alpha
# ngamma_low_median_10 %>% 
#   right_join(ratio.z.ancestral.state.df) %>%
#   mutate(log.root = log(root),
#          mean.log.root = mean(log.root)) %>% 
#   ggplot(aes(y = gamma.cat*ngamma_low_median_10.sigma,
#              x = log.root,
#              group = gamma.cat)) +
#   geom_vline(xintercept = 0) +
#   geom_point(alpha = 0.25,
#              color = 'red')+
#   # geom_count(aes(alpha = after_stat(n),
#   #            color = 'red')) +
#   # scale_alpha(range = c(0.5,1)) +
#   guides(color = 'legend') +
#   theme_classic() +
#   geom_vline(xintercept = 0.4252187,
#              linetype = 'dashed') +
#   ylab('Median sigma category') +
#   xlab('Sex ratio ancestral state (ln)') +
#   ggtitle('Z chromosome median sigma categories vs ancestral state') 
# ggsave('CAGEE/global_figures/n_gamma/Z/Z chromosome ancestral state vs median sigma category alpha.png')




#### all genes
### load data
## ratio
# 10 categories 
ngamma_low_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  dplyr::select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10)

# save data
write.csv(ngamma_low_ratio_10,
          file = 'CAGEE/csv_files/ngamma_low_ratio_10.csv')

# get sigma 
ngamma_low_ratio_10.sigma = read.delim('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

# ## check tree for root node
# ape::read.tree('CAGEE/data/10_species_birds_edit.nwk') %>% 
#   plot()
# edgelabels(cex = 0.8) 
# nodelabels(cex = 0.8)
# dev.off()

## get ancestral state at root 
ratio.ancestral.state.df = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/ancestral_states.tab') %>% 
  dplyr::select(TranscriptID,
         '<11>') %>% 
  rename(root = '<11>') %>% 
  rename(transcript = TranscriptID)

#### graph
### compare ancestral state values with sigma categories
# use log of root value
ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df) %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = -0.01590597,
             linetype = 'dashed') +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() +
  geom_text(data = ngamma_low_ratio_10 %>% 
              group_by(gamma.cat) %>% 
              summarise(count = n()),
            aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
                label = paste0(count),
                y = 0.85)) +
  xlab('Sigma category') +
  ylab('Sex ratio ancestral state (log)') +
  ggtitle('All genes sigma categories vs ancestral state') 
ggsave('CAGEE/global_figures/n_gamma/Z/All genes sigma categories vs ancestral state.png')

# filter 
ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df) %>%
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = -0.01590597,
             linetype = 'dashed') +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() +
  geom_text(data = ngamma_low_ratio_10 %>% 
              group_by(gamma.cat) %>% 
              summarise(count = n()),
            aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
                label = paste0(count),
                y = -0.25)) +
  xlab('Sigma category') +
  ylab('Sex ratio ancestral state (log)') +
  ggtitle('All genes sigma categories vs ancestral state') +
  ylim(-0.25,
       0.25)
ggsave('CAGEE/global_figures/n_gamma/Z/All genes sigma categories vs ancestral state filter.png')


# remove Z genes 
ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df) %>% 
  filter(!transcript %in% ratio.z.ancestral.state.df$transcript) %>% 
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = -0.01590597,
             linetype = 'dashed') +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() +
  geom_text(data = ngamma_low_ratio_10 %>% 
              group_by(gamma.cat) %>% 
              summarise(count = n()),
            aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
                label = paste0(count),
                y = -0.25)) +
  xlab('Sigma category') +
  ylab('Sex ratio ancestral state (log)') +
  ggtitle('All genes sigma categories vs ancestral state, no Z') +
  ylim(-0.25,
       0.25)
ggsave('CAGEE/global_figures/n_gamma/Z/All genes sigma categories vs ancestral state filter no Z.png')


### statistics
library(nptest)
# Test whether each group differs from average
# could add * to boxplots
# wilcoxon signed rank test with permutation
# no Z
ngamma_low_ratio_10.pvalue = ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df)  %>%
  filter(!transcript %in% ratio.z.ancestral.state.df$transcript) %>% 
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  group_by(gamma.cat) %>%
  summarise(P = np.loc.test(log.root, 
                            mu = -0.033238,
                            R = 1000,
                            parallel = T,
                            median.test = T,
                            symmetric = T)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns")) %>% 
  mutate(p.adj = p.adjust(P,
                          n = 10),
         Sig.adj = ifelse(p.adj < 0.05, "*", "ns"))

# get all stats
tmp = ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df)  %>%
  filter(!transcript %in% ratio.z.ancestral.state.df$transcript) %>% 
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  group_by(gamma.cat) %>%
  summarise(estimate = np.loc.test(log.root, 
                            mu = -0.033238,
                            R = 1000,
                            parallel = T,
                            median.test = T,
                            symmetric = T)$estimate,
            estimate.exp = exp(estimate),
            mu.exp = exp(-0.033238),
            t.statistic = np.loc.test(log.root, 
                                   mu = -0.033238,
                                   R = 1000,
                                   parallel = T,
                                   median.test = T,
                                   symmetric = T)$statistic,
            P = np.loc.test(log.root, 
                            mu = -0.033238,
                            R = 1000,
                            parallel = T,
                            median.test = T,
                            symmetric = T)$p.value,
            P.adj = p.adjust(P,
                             n = 10)) 


# graph with significance 
# filter 
# no Z
ngamma_low_ratio_10 %>% 
  left_join(ratio.ancestral.state.df) %>%
  filter(!transcript %in% ratio.z.ancestral.state.df$transcript) %>% 
  mutate(log.root = log(root),
         mean.log.root = mean(log.root)) %>% 
  ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
             y = log.root,
             group = gamma.cat)) +
  geom_hline(yintercept = -0.01590597,
             linetype = 'dashed') +
  geom_hline(yintercept = 0) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = ngamma_low_ratio_10.pvalue, 
             aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
                 label = Sig.adj,
                 y = 0.25)) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() +
  geom_text(data = ngamma_low_ratio_10 %>% 
              group_by(gamma.cat) %>% 
              summarise(count = n()),
            aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
                label = paste0(count),
                y = -0.25)) +
  xlab('Sigma category') +
  ylab('Sex ratio ancestral state (log)') +
  ggtitle('All genes sigma categories vs ancestral state, no Z') +
  ylim(-0.25,
       0.25)
ggsave('CAGEE/global_figures/n_gamma/Z/All genes sigma categories vs ancestral state filter no Z sig.png')


#### compare sigma category for ratio and median ####
### create reduced dataframe with ranked gamma categories
## ratio
ngamma_low_ratio_10.reduce = ngamma_low_ratio_10 %>% 
  dplyr::select(transcript,
                gamma.cat) %>% 
  rename(gamma.cat.ratio = gamma.cat) 

# rank gamma.cat
ngamma_low_ratio_10.ranks = ngamma_low_ratio_10.reduce %>% 
  dplyr::select(gamma.cat.ratio) %>% 
  distinct() %>% 
  dplyr::arrange(gamma.cat.ratio) %>% 
  rownames_to_column('rank.ratio')  

# combine data
# add sigma
ngamma_low_ratio_10.reduce = ngamma_low_ratio_10.reduce %>% 
  left_join(ngamma_low_ratio_10.ranks) %>% 
  mutate(sigma.ratio = gamma.cat.ratio*ngamma_low_ratio_10.sigma)

## 10 categories
## transcriptome
# median
# 10 categories 
ngamma_low_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 

# save data
write.csv(ngamma_low_median_10,
          file = 'CAGEE/csv_files/ngamma_low_median_10.csv')


## median
# get sigma 
ngamma_low_median_10.sigma = read.delim('CAGEE/median/median_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

ngamma_low_median_10.reduce = ngamma_low_median_10 %>% 
  dplyr::select(transcript,
                gamma.cat) %>% 
  rename(gamma.cat.median = gamma.cat) 

# rank gamma.cat
ngamma_low_median_10.ranks = ngamma_low_median_10.reduce %>% 
  dplyr::select(gamma.cat.median) %>% 
  distinct() %>% 
  dplyr::arrange(gamma.cat.median) %>% 
  rownames_to_column('rank.median') 

# combine data
# add sigma
ngamma_low_median_10.reduce = ngamma_low_median_10.reduce %>% 
  left_join(ngamma_low_median_10.ranks) %>% 
  mutate(sigma.median = gamma.cat.median*ngamma_low_median_10.sigma)

### compare distribution of rate category
## ratio
ngamma_low_median_10.reduce %>% 
  full_join(ngamma_low_ratio_10.reduce) %>%
  ggplot(aes(x = sigma.ratio,
             y = sigma.median,
             group = gamma.cat.ratio)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/compare/ratio vs median sigma cat.png')

# rank
ngamma_low_median_10.reduce %>% 
  full_join(ngamma_low_ratio_10.reduce) %>% 
  ggplot(aes(x = as.numeric(rank.ratio),
             y = as.numeric(rank.median),
             group = gamma.cat.ratio)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/compare/rank ratio vs median sigma cat.png')

## median
ngamma_low_median_10.reduce %>% 
  full_join(ngamma_low_ratio_10.reduce) %>%
  ggplot(aes(y = sigma.ratio,
             x = sigma.median,
             group = gamma.cat.median)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/compare/median vs ratio sigma cat.png')

# rank
ngamma_low_median_10.reduce %>% 
  full_join(ngamma_low_ratio_10.reduce) %>% 
  ggplot(aes(y = as.numeric(rank.ratio),
             x = as.numeric(rank.median),
             group = gamma.cat.median)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, 
               geom="point", 
               shape=20, 
               size=2, 
               color="red", 
               fill="red") +
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/compare/rank median vs ratio sigma cat.png')


#### gamma cats Z enrichment ####
#### median 
## median
# 10 categories 
ngamma_low_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 

# get sigma 
ngamma_low_median_10.sigma = read.delim('CAGEE/median/median_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

## load Z chromosome position
## load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

# add chromsome position to ngamma_low_median_10
ngamma_low_median_10_pos = ngamma_low_median_10 %>% 
  left_join(data.gene.chromosome %>% 
              select(Symbol,
                     chromosome_name) %>% 
              rename(transcript = Symbol))
  
  
## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma_low_median_10_pos.table = ngamma_low_median_10_pos %>% 
  select(-c(transcript,
            liklihood,
            keep,
            top.cat)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat = as.numeric(as.character(gamma.cat))) %>% 
  filter(Freq > 0) %>% 
  group_by(type,
           gene.list,
           gamma.count,
           chromosome_name) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(percent = 100*Freq/total)

## add mean and SE
# filter microchromosomes 
ngamma_low_median_10_pos.table.stat = ngamma_low_median_10_pos.table %>% 
  mutate(Chr = ifelse(chromosome_name == 'Z',
                      'Z',
                      'Auto')) %>% 
  filter(total > 26) %>%  
  group_by(Chr,
           gamma.cat) %>% 
  mutate(Percent.cat = mean(percent),
         SE.cat = sd(percent)/sqrt(length((percent)))) %>% 
  mutate(gamma.cat.value = gamma.cat*ngamma_low_median_10.sigma)
  

#### graph n_gamma_cats data low position
### graph genes per category per chromsoome
# ngamma_low_median_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_median_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlab('Sigma category')
# ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name.png')
# 
# ## filter for chrosmomes with genes > 150
# ngamma_low_median_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   filter(total > 26) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_median_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category')
# ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name filter.png')
# 
# ## add label
# ngamma_low_median_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   filter(total > 26) %>% 
#   mutate(label = ifelse(gamma.cat == max(ngamma_low_median_10_pos.table$gamma.cat),
#                         as.character(chromosome_name),
#                         NA)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_median_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   ggrepel::geom_label_repel(aes(label = label),
#                             max.overlaps = 50,
#                             color = 'orange') +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category')
# ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name filter label.png')
# 
# ### graph autosomes vs Z
# ngamma_low_median_10_pos.table.stat %>%
#   dplyr::select(gamma.cat.value,
#                 Percent.cat,
#                 Chr,
#                 SE.cat) %>% 
#   ggplot(aes(x = gamma.cat.value,
#              y = Percent.cat,
#              color = Chr)) +
#   geom_errorbar(aes( ymin = Percent.cat - SE.cat,
#                      ymax = Percent.cat + SE.cat)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category') +
#   ggtitle('Percent of genes per sigma category, standard error')
# ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name filter Autosomes.png')

### statistics
## one sample two tailed t test
## create empty dataframe
ngamma_low_median_10_stat = data.frame()

## run ttest for each rate category
for (i in unique(ngamma_low_median_10_pos.table.stat$gamma.cat.value)) {
  # subset data
  tmp = ngamma_low_median_10_pos.table.stat %>%
    filter(gamma.cat.value == i) 
  
  # get Z value
  tmpz = tmp %>% 
    filter(chromosome_name == 'Z')
  
  # get autosome z value
  tmpa = tmp %>% 
    filter(chromosome_name != 'Z')
  
  # run one sample ttest
  tmpy = t.test(tmpa$percent,
                mu = tmpz$percent,
                alternative = 'two.sided')
  
  # create tmp dataframe of results
  tmp.df = data.frame(
    gamma.cat.value = i,
    type = tmpz$type,
    gamma.count = tmpz$gamma.count,
    auto.percent.mean = unique(tmpa$Percent.cat),
    auto.percent.SE = unique(tmpa$SE.cat),
    z.percent = tmpz$percent,
    t.statistic = tmpy$statistic,
    df = tmpy$parameter,
    p.value = tmpy$p.value,
    alternative = tmpy$alternative,
    method = tmpy$method
  ) 
  
  # combine results
  ngamma_low_median_10_stat =  ngamma_low_median_10_stat %>% 
    rbind(tmp.df)
  
}
 
## adjust pvalue for multiple testing
ngamma_low_median_10_stat = ngamma_low_median_10_stat %>% 
  mutate(p.value.adj = p.adjust(p.value,
                                method = 'fdr',
                                n = length(unique(ngamma_low_median_10_pos.table.stat$gamma.cat.value))))
 
## add significance to graph
# ngamma_low_median_10_pos.table.stat %>%
#   dplyr::select(gamma.cat.value,
#                 Percent.cat,
#                 Chr,
#                 SE.cat) %>% 
#   ggplot(aes(x = gamma.cat.value,
#              y = Percent.cat,
#              color = Chr)) +
#   geom_errorbar(aes( ymin = Percent.cat - SE.cat,
#                      ymax = Percent.cat + SE.cat)) +
#   geom_text(data = ngamma_low_median_10_stat %>% 
#                mutate(sig = case_when(p.value.adj <= 0.001 ~ '***',
#                                       p.value.adj > 0.001 & p.value.adj <= 0.01 ~ '**',
#                                       p.value.adj > 0.01 & p.value.adj <= 0.05 ~ '*',
#                                    TRUE ~ NA),
#                       Percent = case_when(auto.percent.mean >= z.percent ~ auto.percent.mean,
#                                           auto.percent.mean < z.percent ~ z.percent)),
#              aes(x = gamma.cat.value,
#                  y = Percent + 0.9,
#                  label = sig),
#              color = 'black',
#             size = 5) + 
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   labs(x = 'Median gamma category',
#        y = 'Percent of genes',
#        caption = 'FDR < *0.05, **0.01, ***0.001') + 
#   ggtitle('Percent of genes per sigma category, standard error')
# ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name filter Autosomes sig.png')

# paper
ngamma_low_median_10_pos.table.stat %>%
  dplyr::select(gamma.cat.value,
                Percent.cat,
                Chr,
                SE.cat) %>% 
  ggplot(aes(x = gamma.cat.value,
             y = Percent.cat,
             color = Chr)) +
  geom_errorbar(aes( ymin = Percent.cat - SE.cat,
                     ymax = Percent.cat + SE.cat),
                size = 0.15) +
  geom_text(data = ngamma_low_median_10_stat %>% 
              mutate(sig = case_when(p.value.adj <= 0.001 ~ '***',
                                     p.value.adj > 0.001 & p.value.adj <= 0.01 ~ '**',
                                     p.value.adj > 0.01 & p.value.adj <= 0.05 ~ '*',
                                     TRUE ~ NA),
                     Percent = case_when(auto.percent.mean >= z.percent ~ auto.percent.mean,
                                         auto.percent.mean < z.percent ~ z.percent)),
            aes(x = gamma.cat.value,
                y = Percent + 0.9,
                label = sig),
            color = 'black',
            size = 2) + 
  geom_line() +
  geom_point(size = 0.75) +
  theme_classic(base_size = 8) +
  scale_color_manual(values = c('black',
                                'red'))+
  labs(x = 'Median gamma category',
       y = 'Percent of genes',
       # caption = 'FDR < *0.05, **0.01, ***0.001'
  ) +
  theme(legend.position = 'null')
ggsave('CAGEE/global_figures/n_gamma/Z/Percent of genes per sigma category low chromosome_name filter Autosomes sig.pdf',
       height = 1.34,
       width = 3.1,
       units = 'in',
       dpi = 320)
#### ratio 
## ratio
# 10 categories 
ngamma_low_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10) 

# get sigma 
ngamma_low_ratio_10.sigma = read.delim('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/results.txt') %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  filter(CAGEE.1.1 == 'Sigma2') %>% 
  pull(result) %>% 
  as.numeric()

## load Z chromosome position
## load gene chromsome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

# add chromsome position to ngamma_low_ratio_10
ngamma_low_ratio_10_pos = ngamma_low_ratio_10 %>% 
  left_join(data.gene.chromosome %>% 
              select(Symbol,
                     chromosome_name) %>% 
              rename(transcript = Symbol))


## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma_low_ratio_10_pos.table = ngamma_low_ratio_10_pos %>% 
  select(-c(transcript,
            liklihood,
            keep,
            top.cat)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat = as.numeric(as.character(gamma.cat))) %>% 
  filter(Freq > 0) %>% 
  group_by(type,
           gene.list,
           gamma.count,
           chromosome_name) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(percent = 100*Freq/total)

## add mean and SE
# filter microchromosomes 
ngamma_low_ratio_10_pos.table.stat = ngamma_low_ratio_10_pos.table %>% 
  mutate(Chr = ifelse(chromosome_name == 'Z',
                      'Z',
                      'Auto')) %>% 
  filter(total > 26) %>%  
  group_by(Chr,
           gamma.cat) %>% 
  mutate(Percent.cat = mean(percent),
         SE.cat = sd(percent)/sqrt(length((percent)))) %>% 
  mutate(gamma.cat.value = gamma.cat*ngamma_low_ratio_10.sigma)

#### graph n_gamma_cats data low position
# ### graph genes per category per chromsoome
# ngamma_low_ratio_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red')) +
#   xlab('Sigma category') +
#   ggtitle('Sex ratio')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name.png')
# 
# ## filter for chrosmomes with genes > 150
# ngamma_low_ratio_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   filter(total > 150) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category')+
#   ggtitle('Sex ratio')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter.png')
# 
# ## add label
# ngamma_low_ratio_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   filter(total > 150) %>% 
#   mutate(label = ifelse(gamma.cat == max(ngamma_low_ratio_10_pos.table$gamma.cat),
#                         as.character(chromosome_name),
#                         NA)) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = percent,
#              group = chromosome_name,
#              color = Chr)) +
#   geom_line() +
#   geom_point() +
#   ggrepel::geom_label_repel(aes(label = label),
#                             max.overlaps = 50,
#                             color = 'orange') +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category') +
#   ggtitle('Sex ratio')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter label.png')
# 
# ### graph autosomes vs Z
# # no Z
# ngamma_low_ratio_10_pos.table %>% 
#   mutate(Chr = ifelse(chromosome_name == 'Z',
#                       'Z',
#                       'Auto')) %>% 
#   filter(total > 150) %>% 
#   group_by(gamma.cat,
#            Chr) %>% 
#   summarise(Percent = mean(percent),
#             SE = sd(percent)/sqrt(length((percent)))) %>% 
#   ggplot(aes(x = gamma.cat*ngamma_low_ratio_10.sigma,
#              y = Percent,
#              color = Chr)) +
#   geom_errorbar(aes( ymin = Percent - SE,
#                      ymax = Percent + SE)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'white'))+
#   xlab('Sigma category') +
#   ggtitle('Sex ratio percent of genes per sigma category, standard error')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter Autosomes no z.png')
# 
# # with Z
# ngamma_low_ratio_10_pos.table.stat %>%
#   dplyr::select(gamma.cat.value,
#                 Percent.cat,
#                 Chr,
#                 SE.cat) %>% 
#   ggplot(aes(x = gamma.cat.value,
#              y = Percent.cat,
#              color = Chr)) +
#   geom_errorbar(aes( ymin = Percent.cat - SE.cat,
#                      ymax = Percent.cat + SE.cat)) +
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   xlab('Sigma category') +
#   ggtitle('Sex ratio percent of genes per sigma category, standard error')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter Autosomes.png')

### statistics
## one sample two tailed t test
## create empty dataframe
ngamma_low_ratio_10_stat = data.frame()

## run ttest for each rate category
for (i in unique(ngamma_low_ratio_10_pos.table.stat$gamma.cat.value)) {
  # subset data
  tmp = ngamma_low_ratio_10_pos.table.stat %>%
    filter(gamma.cat.value == i) 
  
  # get Z value
  tmpz = tmp %>% 
    filter(chromosome_name == 'Z')
  
  # get autosome z value
  tmpa = tmp %>% 
    filter(chromosome_name != 'Z')
  
  # run one sample ttest
  tmpy = t.test(tmpa$percent,
                mu = tmpz$percent,
                alternative = 'two.sided')
  
  # create tmp dataframe of results
  tmp.df = data.frame(
    gamma.cat.value = i,
    type = tmpz$type,
    gamma.count = tmpz$gamma.count,
    auto.percent.mean = unique(tmpa$Percent.cat),
    auto.percent.SE = unique(tmpa$SE.cat),
    z.percent = tmpz$percent,
    t.statistic = tmpy$statistic,
    df = tmpy$parameter,
    p.value = tmpy$p.value,
    alternative = tmpy$alternative,
    method = tmpy$method
  ) 
  
  # combine results
  ngamma_low_ratio_10_stat =  ngamma_low_ratio_10_stat %>% 
    rbind(tmp.df)
  
}

## adjust pvalue for multiple testing
ngamma_low_ratio_10_stat = ngamma_low_ratio_10_stat %>% 
  mutate(p.value.adj = p.adjust(p.value,
                                method = 'fdr',
                                n = length(unique(ngamma_low_ratio_10_pos.table.stat$gamma.cat.value))))

## add significance to graph
# ngamma_low_ratio_10_pos.table.stat %>%
#   dplyr::select(gamma.cat.value,
#                 Percent.cat,
#                 Chr,
#                 SE.cat) %>% 
#   ggplot(aes(x = gamma.cat.value,
#              y = Percent.cat,
#              color = Chr)) +
#   geom_errorbar(aes( ymin = Percent.cat - SE.cat,
#                      ymax = Percent.cat + SE.cat)) +
#   geom_text(data = ngamma_low_ratio_10_stat %>% 
#               mutate(sig = case_when(p.value.adj <= 0.001 ~ '***',
#                                      p.value.adj > 0.001 & p.value.adj <= 0.01 ~ '**',
#                                      p.value.adj > 0.01 & p.value.adj <= 0.05 ~ '*',
#                                      TRUE ~ NA),
#                      Percent = case_when(auto.percent.mean >= z.percent ~ auto.percent.mean,
#                                          auto.percent.mean < z.percent ~ z.percent)),
#             aes(x = gamma.cat.value,
#                 y = Percent + 1.1,
#                 label = sig),
#             color = 'black',
#             size = 5) + 
#   geom_line() +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('black',
#                                 'red'))+
#   labs(x = 'Ratio gamma category',
#        y = 'Percent of genes',
#        caption = 'FDR < *0.05, **0.01, ***0.001') + 
#   ggtitle('Sex ratio percent of genes per sigma category, standard error')
# ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter Autosomes sig.png')

# paper
ngamma_low_ratio_10_pos.table.stat %>%
  dplyr::select(gamma.cat.value,
                Percent.cat,
                Chr,
                SE.cat) %>% 
  ggplot(aes(x = gamma.cat.value,
             y = Percent.cat,
             color = Chr)) +
  geom_errorbar(aes( ymin = Percent.cat - SE.cat,
                     ymax = Percent.cat + SE.cat),
                size = 0.15) +
  geom_text(data = ngamma_low_ratio_10_stat %>% 
              mutate(sig = case_when(p.value.adj <= 0.001 ~ '***',
                                     p.value.adj > 0.001 & p.value.adj <= 0.01 ~ '**',
                                     p.value.adj > 0.01 & p.value.adj <= 0.05 ~ '*',
                                     TRUE ~ NA),
                     Percent = case_when(auto.percent.mean >= z.percent ~ auto.percent.mean,
                                         auto.percent.mean < z.percent ~ z.percent)),
            aes(x = gamma.cat.value,
                y = Percent + 1.1,
                label = sig),
            color = 'black',
            size = 2) + 
  geom_line() +
  geom_point(size = 0.75) +
  theme_classic(base_size = 8) +
  scale_color_manual(values = c('black',
                                'red'))+
  labs(x = 'Ratio gamma category',
       y = 'Percent of genes',
       # caption = 'FDR < *0.05, **0.01, ***0.001'
       ) +
  theme(legend.position = 'null')
ggsave('CAGEE/global_figures/n_gamma/Z/ratio Percent of genes per sigma category low chromosome_name filter Autosomes sig.pdf',
       height = 1.34,
       width = 3.1,
       units = 'in',
       dpi = 320)


#### save stat results together
write.csv(rbind(ngamma_low_median_10_stat, 
                ngamma_low_ratio_10_stat),
          file = 'CAGEE/global_figures/n_gamma/Z/ngamma_z_stats.csv',
          row.names = FALSE)



#### compare ratio to gene expression evolution
## create dataframe of ratio
ngamma_low_df = ngamma_low_ratio_10 %>% 
  dplyr::select(transcript,
                gamma.cat) %>% 
  mutate(gamma.cat.ratio = gamma.cat*ngamma_low_ratio_10.sigma) %>% 
  select(-gamma.cat) %>% 
  full_join(ngamma_low_median_10 %>% 
              dplyr::select(transcript,
                            gamma.cat) %>% 
              mutate(gamma.cat.median = gamma.cat*ngamma_low_median_10.sigma) %>% 
              dplyr::select(-gamma.cat))

## save for supplemental
# rename variables
write.csv(ngamma_low_df %>% 
            dplyr::rename('median_sigma' = 'gamma.cat.median')%>% 
            dplyr::rename('M:F_sigma' = 'gamma.cat.ratio'),
          'CAGEE/global_figures/n_gamma/compare/Dataset_S3.csv',
          row.names = F)

## stats
ngamma_low_df_sum = ngamma_low_df %>% 
  dplyr::select(-transcript) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat.ratio = as.numeric(as.character(gamma.cat.ratio))) %>% 
  mutate(gamma.cat.median = as.numeric(as.character(gamma.cat.median)))

# add values
ngamma_low_df_sum = ngamma_low_df %>% 
  dplyr::select(gamma.cat.ratio) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename(ratio.count = Freq) %>% 
  mutate(gamma.cat.ratio = as.numeric(as.character(gamma.cat.ratio))) %>% 
  full_join(ngamma_low_df_sum) %>% 
  full_join(ngamma_low_df %>% 
              dplyr::select(gamma.cat.median) %>% 
              table() %>% 
              as.data.frame() %>% 
              dplyr::rename(median.count = Freq) %>% 
              mutate(gamma.cat.median = as.numeric(as.character(gamma.cat.median)))) 

# test for over-representation
ngamma_low_df_sum = ngamma_low_df_sum %>% 
  mutate(enrichment.p.val = phyper(Freq-1,
                                   median.count,
                                   10677-median.count,
                                   ratio.count,
                                   lower.tail = F)) %>% 
  mutate(enrichment.p.val.fdr = p.adjust(enrichment.p.val,
                                         method = 'fdr',
                                         n = 100)) %>% 
  mutate(sig = case_when(enrichment.p.val.fdr < 0.05 & enrichment.p.val.fdr > 0.01 ~ "*",
                         enrichment.p.val.fdr < 0.01 & enrichment.p.val.fdr > 0.001~ "**",
                         enrichment.p.val.fdr < 0.001 ~ "***",
                         TRUE ~ NA))
  
  
  
### graph
## heatmap
# ngamma_low_df_sum %>% 
#   ggplot(aes(y = as.numeric(as.character(gamma.cat.ratio)),
#              x = as.numeric(as.character(gamma.cat.median)),
#              label = sig,
#              fill = Freq)) +
#   geom_tile(color = 'darkgrey') +
#   geom_text(size = 2) +
#   labs(x = 'Median gamma category',
#        y = 'Ratio gamma category',
#        caption = 'FDR < *0.05, **0.01, ***0.001') +
#   scale_fill_gradient2(low = 'white',
#                        high = 'red',
#                        limits = c(0,400)) + 
#   theme_classic() +
#   coord_fixed(ratio = 0.05/0.002)
# ggsave('CAGEE/global_figures/n_gamma/compare/heatmap median vs ratio sigma cat.png')

# paper
ngamma_low_df_sum %>% 
  ggplot(aes(y = as.numeric(as.character(gamma.cat.ratio)),
             x = as.numeric(as.character(gamma.cat.median)),
             label = sig,
             fill = Freq)) +
  geom_tile(color = 'black') +
  geom_text(size = 1,
            vjust = 0.75) +
  labs(x = 'Median gamma category',
       y = 'Ratio gamma category',
       # caption = 'FDR < *0.05, **0.01, ***0.001'
       ) +
  scale_fill_gradient2(low = 'white',
                       high = 'red',
                       limits = c(0,400)) + 
  theme_classic(base_size = 8) +
  coord_fixed(ratio = 0.05/0.002) + 
  theme(legend.key.size = unit(0.075,
                               "in"))
  # theme(legend.key.height = unit(0.3, "in"),
  #   legend.key.width = unit(0.5, "in"))
ggsave('CAGEE/global_figures/n_gamma/compare/heatmap median vs ratio sigma cat.pdf',
       height = 2.6895,
       width = 3.25,
       units = 'in',
       dpi = 320)


#### set up GO terms  genes ####
#### all genes
#### load data
### ratio
## 10 categories 
ngamma_low_ratio_10 = read_tsv('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'ratio',
         gene.list = 'all',
         gamma.count = 10)


# get top category value
ngamma_low_ratio_10.max = ngamma_low_ratio_10 %>% 
  pull(gamma.cat) %>% 
  max()

# get bottom category value
ngamma_low_ratio_10.min = ngamma_low_ratio_10 %>% 
  pull(gamma.cat) %>% 
  min()

## filter down to list of genes in fastest rate category
# get list of transcripts in fastest rate category
ngamma_low_ratio_10_fast_genes = ngamma_low_ratio_10 %>% 
  filter(gamma.cat == ngamma_low_ratio_10.max) %>% 
  pull(transcript)

## filter down to list of genes in fastest rate category
# get list of transcripts in fastest rate category
ngamma_low_ratio_10_slow_genes = ngamma_low_ratio_10 %>% 
  filter(gamma.cat == ngamma_low_ratio_10.min) %>% 
  pull(transcript)

### median
## 10 categories 
ngamma_low_median_10 = read_tsv('CAGEE/median/median_cagee_outs_n_gamma_cats_10/category_likelihoods.txt') %>% 
  select(-c('...12')) %>% 
  rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(type = 'median',
         gene.list = 'all',
         gamma.count = 10) 


# get top category value
ngamma_low_median_10.max = ngamma_low_median_10 %>% 
  pull(gamma.cat) %>% 
  max()

# get bottom category value
ngamma_low_median_10.min = ngamma_low_median_10 %>% 
  pull(gamma.cat) %>% 
  min()

## filter down to list of genes in fastest rate category
# get list of transcripts in fastest rate category
ngamma_low_median_10_fast_genes = ngamma_low_median_10 %>% 
  filter(gamma.cat == ngamma_low_median_10.max) %>% 
  pull(transcript)


## filter down to list of genes in slowest rate category
ngamma_low_median_10_slow_genes = ngamma_low_median_10 %>% 
  filter(gamma.cat == ngamma_low_median_10.min) %>% 
  pull(transcript)

### combine into one dataframe
## fast
ngamma_10_fast_genes = data.frame(transcript = ngamma_low_median_10_fast_genes,
                                  type = 'median') %>% 
  rbind(
  data.frame(transcript = ngamma_low_ratio_10_fast_genes,
        type = 'ratio')
  )


## slow
ngamma_10_slow_genes = data.frame(transcript = ngamma_low_median_10_slow_genes,
                                  type = 'median') %>% 
  rbind(
    data.frame(transcript = ngamma_low_ratio_10_slow_genes,
               type = 'ratio')
  )


### get list of Z genes
## load gene chromosome position
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

# add chromosome to fast gene data
ngamma_10_fast_genes = ngamma_10_fast_genes %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Ensembl.GeneIDs,
                     external_gene_name,
                     chromosome_name) %>% 
              rename(transcript = external_gene_name) %>% 
              distinct())

# add chromosome to slow gene data
ngamma_10_slow_genes = ngamma_10_slow_genes %>% 
  left_join(data.gene.chromosome %>% 
              dplyr::select(Ensembl.GeneIDs,
                     external_gene_name,
                     chromosome_name) %>% 
              rename(transcript = external_gene_name) %>% 
              distinct())

#### create Zebrafinch Go database ####
# ### use makeOrgPackage in AnnotationForge
# # https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
# ## load library
# library(AnnotationForge)
# 
# # use NCBI to get data for zebrafinch
# # takes a long time to download
# makeOrgPackageFromNCBI(version = "0.1",
#                        author = "Isaac Miller-Crews <imillerc@iu.edu>",
#                        maintainer = "Isaac Miller-Crews <imillerc@iu.edu",
#                        outputDir = "./CAGEE/data/",
#                        tax_id = "59729",
#                        genus = "Taeniopygia",
#                        species = "guttata")
# 
# # can rebuild database if NCBI.sqlite is already downloaded in working directory
# makeOrgPackageFromNCBI(version = "0.1",
#                        author = "Isaac Miller-Crews <imillerc@iu.edu>",
#                        maintainer = "Isaac Miller-Crews <imillerc@iu.edu",
#                        outputDir = "./CAGEE/data/",
#                        tax_id = "59729",
#                        genus = "Taeniopygia",
#                        species = "guttata",
#                        rebuildCache = F) # stops files from being redownloaded, need files in output directory
# 
# ## Makes an organism package for Zebra Finch data.frames:
# finchFile <- system.file("extdata","finch_info.txt",
#                          package="AnnotationForge")
# finch <- read.table(finchFile,
#                     sep="\t")
# 
# ## Now prepare some data.frames
# fSym <- finch[,c(2,3,9)]
# fSym <- fSym[fSym[,2]!="-",]
# fSym <- fSym[fSym[,3]!="-",]
# colnames(fSym) <- c("GID","SYMBOL","GENENAME")
# 
# fChr <- finch[,c(2,7)]
# fChr <- fChr[fChr[,2]!="-",]
# colnames(fChr) <- c("GID","CHROMOSOME")
# 
# finchGOFile <- system.file("extdata","GO_finch.txt",
#                            package="AnnotationForge")
# fGO <- read.table(finchGOFile,sep="\t")
# fGO <- fGO[fGO[,2]!="",]
# fGO <- fGO[fGO[,3]!="",]
# colnames(fGO) <- c("GID","GO","EVIDENCE")
# 
# ## Then call the function
# makeOrgPackage(gene_info=fSym,
#                chromosome=fChr,
#                go=fGO,
#                version="0.1",
#                maintainer="Isaac Miller-Crews <imillerc@iu.edu>",
#                author="Isaac Miller-Crews <imillerc@iu.edu>",
#                outputDir = "./CAGEE/NCBI/",
#                tax_id="59729",
#                genus="Taeniopygia",
#                species="guttata",
#                goTable="go")
# 
# ## then you can call install.packages based on the return value
# install.packages("./CAGEE/NCBI/org.Tguttata.eg.db", repos=NULL)

#### make own Zebra Finch reference
## load in all gene data
data = readxl::read_excel('CAGEE/data/normalizedCounts.xlsx')

## subset to gene columns
data.genes = data %>% 
  select("ZebraFinchProteinID",
         "GeneName",
         "GeneDescription")


### use ensemble to get gene GO terms, chromosome, and evidence
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl")

# List available datasets (species)
listDatasets(ensembl) %>% 
  View()

# Select a dataset
# Zebra finch
ensembl.zf <- useDataset("tguttata_gene_ensembl",
                      mart = ensembl)

# List available attributes 
listAttributes(ensembl.zf) %>% 
  View()

# Retrieve gene information based on a filter, such as a list of gene symbols
attributes <- c("ensembl_gene_id", 
                "external_gene_name", 
                "chromosome_name",
                "go_id",
                "go_linkage_type",
                "refseq_peptide_predicted")

# Get the data from GeneName
data.genes.name <- getBM(attributes = attributes,
                   filters = "external_gene_name", 
                   values = data.genes %>% 
                     pull(GeneName), 
                   mart = ensembl.zf)

# Get the data from protein
data.genes.protein <- getBM(attributes = attributes,
                         filters = "refseq_peptide_predicted", 
                         values = data.genes %>% 
                           dplyr::select(ZebraFinchProteinID) %>% 
                           separate_wider_delim(ZebraFinchProteinID,
                                                names = c('ZebraFinchProteinID',
                                                          NA),
                                                delim = '.') %>% 
                           pull(ZebraFinchProteinID), 
                         mart = ensembl.zf)

# combine protein and genename results
data.genes.name.protein = full_join(data.genes.name,
                                    data.genes.protein)



### Makes an organism package for Zebra Finch data.frames:
## Need dataframe of "GID","SYMBOL","GENENAME"
data.genes.sym = data.genes
colnames(data.genes.sym) = c("GID",
                         "SYMBOL",
                         "GENENAME")

## Need dataframe of "GID","CHROMOSOME"
# add GID to chromosome data
# replace NA with A
data.genes.chr = data.genes %>% 
  left_join(data.genes.name %>% 
              dplyr::select(external_gene_name,
                     chromosome_name) %>% 
              rename(GeneName = external_gene_name) %>% 
              distinct()) %>% 
  dplyr::select(ZebraFinchProteinID,
                chromosome_name) %>% 
  mutate(chromosome_name = ifelse(is.na(chromosome_name),
                                  'A',
                                  chromosome_name))

colnames(data.genes.chr) = c("GID",
                             "CHROMOSOME")

## Need dataframe of "GID","GO","EVIDENCE"
# remove NA and blank
data.genes.go = data.genes %>% 
  left_join(data.genes.name %>% 
              dplyr::select(external_gene_name,
                            go_id,
                            go_linkage_type) %>% 
              rename(GeneName = external_gene_name) %>% 
              distinct()) %>% 
  dplyr::select(ZebraFinchProteinID,
                go_id,
                go_linkage_type) %>% 
  na.omit() %>% 
  filter(go_id != "")

colnames(data.genes.go) = c("GID",
                             "GO",
                             "EVIDENCE")

## Then call the function
makeOrgPackage(gene_info=data.genes.sym,
               chromosome=data.genes.chr,
               go=data.genes.go,
               version="0.1",
               maintainer="Isaac Miller-Crews <imillerc@iu.edu>",
               author="Isaac Miller-Crews <imillerc@iu.edu>",
               outputDir = "./CAGEE/GO_Ref/",
               tax_id="59729",
               genus="Taeniopygia",
               species="guttata2",
               goTable="go")


## then you can call install.packages based on the return value
install.packages("./CAGEE/GO_Ref/org.Tguttata2.eg.db", repos=NULL)



### could also try downloading from annotation hub
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("taeniopygia guttata","orgdb"))
orgdb <- hub[["AH108527"]]


#### calculate GO terms  genes ####

## load GO libraries 
# load libraries 
library(enrichplot)
library(clusterProfiler)
library(org.Tguttata2.eg.db) # need to make database in code above
library(AnnotationDbi)

#### Ratio
### fast genes
## all
# can use gene symbol, don't need ensembl ids
ngamma_low_ratio_10_fast_genes_all <- enrichGO(gene = ngamma_10_fast_genes %>% 
                                             filter(type == 'ratio') %>% 
                                             pull(transcript) %>% 
                                             unique(), 
                                           OrgDb = "org.Tguttata2.eg.db",
                                           keyType = "SYMBOL", 
                                           ont = "BP",
                                           pvalueCutoff = 0.1,
                                           pAdjustMethod = 'fdr',
                                           minGSSize = 10,
                                           maxGSSize = 500
                                           )


# graph GO results
# create plot
ngamma_low_ratio_10_fast_genes_all.fit <- plot(barplot(ngamma_low_ratio_10_fast_genes_all,
                                                   showCategory = 15))

ngamma_low_ratio_10_fast_genes_all.fit

# save plot
png("./CAGEE/global_figures/GO/GO ratio fast genes all barplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(ngamma_low_ratio_10_fast_genes_all.fit)
dev.off()

# create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

# graph tree
png("./CAGEE/global_figures/GO/GO ratio fast genes all treeplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
ngamma_low_ratio_10_fast_genes_all %>% 
  pairwise_termsim() %>% 
  treeplot() +
  ggtitle('GO ratio fast genes')
dev.off()

### slow genes
## all
# can use gene symbol, don't need ensembl ids
ngamma_low_ratio_10_slow_genes_all <- enrichGO(gene = ngamma_10_slow_genes %>% 
                                                 filter(type == 'ratio') %>% 
                                                 pull(transcript) %>% 
                                                 unique(), 
                                               OrgDb = "org.Tguttata2.eg.db",
                                               keyType = "SYMBOL", 
                                               ont = "BP",
                                               pvalueCutoff = 0.1,
                                               pAdjustMethod = 'fdr',
                                               minGSSize = 10,
                                               maxGSSize = 500
)


# graph GO results
# create plot
ngamma_low_ratio_10_slow_genes_all.fit <- plot(barplot(ngamma_low_ratio_10_slow_genes_all,
                                                       showCategory = 15))

ngamma_low_ratio_10_slow_genes_all.fit

# save plot
png("./CAGEE/global_figures/GO/GO ratio slow genes all barplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(ngamma_low_ratio_10_slow_genes_all.fit)
dev.off()

# create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

# graph tree
png("./CAGEE/global_figures/GO/GO ratio slow genes all treeplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
ngamma_low_ratio_10_slow_genes_all %>% 
  pairwise_termsim() %>% 
  treeplot() +
  ggtitle('GO ratio slow genes')
dev.off()

#### median
### fast genes
## all
# can use gene symbol, don't need ensembl ids
ngamma_low_median_10_fast_genes_all <- enrichGO(gene = ngamma_10_fast_genes %>% 
                                                 filter(type == 'median') %>% 
                                                 pull(transcript) %>% 
                                                 unique(), 
                                               OrgDb = "org.Tguttata2.eg.db",
                                               keyType = "SYMBOL", 
                                               ont = "BP",
                                               pvalueCutoff = 0.1,
                                               readable = T,
                                               pAdjustMethod = 'fdr',
                                               minGSSize = 10,
                                               maxGSSize = 500
                                               )


# graph GO results
# create plot
ngamma_low_median_10_fast_genes_all.fit <- plot(barplot(ngamma_low_median_10_fast_genes_all,
                                                       showCategory = 15))

ngamma_low_median_10_fast_genes_all.fit

# save plot
png("./CAGEE/global_figures/GO/GO median fast genes all barplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(ngamma_low_median_10_fast_genes_all.fit)
dev.off()

# create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

# graph tree
png("./CAGEE/global_figures/GO/GO median fast genes all treeplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
ngamma_low_median_10_fast_genes_all %>% 
  pairwise_termsim() %>% 
  treeplot()+
  ggtitle('GO median fast genes')
dev.off()

### slow genes
## all
# can use gene symbol, don't need ensembl ids
ngamma_low_median_10_slow_genes_all <- enrichGO(gene = ngamma_10_slow_genes %>% 
                                                  filter(type == 'median') %>% 
                                                  pull(transcript) %>% 
                                                  unique(), 
                                                OrgDb = "org.Tguttata2.eg.db",
                                                keyType = "SYMBOL", 
                                                ont = "BP",
                                                pvalueCutoff = 0.1,
                                                readable = T,
                                                pAdjustMethod = 'fdr',
                                                minGSSize = 10,
                                                maxGSSize = 500
)


# graph GO results
# create plot
ngamma_low_median_10_slow_genes_all.fit <- plot(barplot(ngamma_low_median_10_slow_genes_all,
                                                        showCategory = 15))

ngamma_low_median_10_slow_genes_all.fit

# save plot
png("./CAGEE/global_figures/GO/GO median slow genes all barplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
print(ngamma_low_median_10_slow_genes_all.fit)
dev.off()

# create GO term tree
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
# download enrichplot
# remotes::install_github("GuangchuangYu/enrichplot")

# graph tree
png("./CAGEE/global_figures/GO/GO median slow genes all treeplot.png", 
    res = 300, 
    width = 10, 
    height = 10,
    units = 'in')
ngamma_low_median_10_slow_genes_all %>% 
  pairwise_termsim() %>% 
treeplot() +
  ggtitle('GO median slow genes')
dev.off()


#### simulation results for multi sigma n_gamma_cats ####
### load and combine data for gene expression and ratio simulations
## gene expression simulations
# make unique transcript ID
ngamma_sim_gene_exp_multi = read_tsv('CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_100_0.01/simulation.txt',
                               skip = 4) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_100_0.005/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_100_0.001/simulation.txt',
                 skip = 4)) %>% 
  rbind(read_tsv('CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_100_0.0001/simulation.txt',
                 skip = 4)) %>% 
  mutate(GENE_ID = paste0(GENE_ID,
                          '_',
                          DESC))
## save file
write_tsv(ngamma_sim_gene_exp_multi,
          'CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_100_combined.tsv')


### load CAGEE results
# 
ngamma_multi_sigma = read_tsv('CAGEE/n_gamma_simulation/multi_sigma/multi_sigma_cagee_outs_n_gamma_cats/category_likelihoods.txt')%>% 
  dplyr::select(-c('...6')) %>% 
  dplyr::rename(transcript = "Transcript ID") %>% 
  pivot_longer(cols = -c(transcript),
               values_to = "liklihood",
               names_to = 'gamma.cat') %>% 
  mutate(gamma.cat = as.numeric(gamma.cat)) %>% 
  group_by(transcript) %>%
  mutate(top.cat = min(liklihood)) %>% 
  ungroup() %>% 
  mutate(keep = ifelse(top.cat == liklihood,
                       1,
                       0)) %>% 
  filter(keep == 1)

## get count of genes per gamma category per rate
# calculate percentage
ngamma_multi_sigma.table = ngamma_multi_sigma %>% 
  separate_wider_delim(cols = transcript,
                       delim = '_',
                       names = c(NA,
                                 'rate')) %>% 
  separate_wider_delim(cols = rate,
                       delim = 'SIG',
                       names = c(NA,
                                 'rate')) %>% 
  dplyr::select(c(rate,
            gamma.cat)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(gamma.cat = as.numeric(as.character(gamma.cat))) 

## graph results
# graph distribution of transcripts across categories
ngamma_multi_sigma.table %>% 
  ggplot(aes(x = rate,
             y = Freq,
             group = rate,
             fill = as.factor(gamma.cat))) +
  geom_bar(stat='identity') +
  theme_classic() +
  scale_fill_manual(values = c('grey75',
                                  'grey50',
                                  'grey25',
                                  'grey1')) +
  ggtitle('Multiple sigma sim') +
  labs(fill = 'Gamma Category') 
ggsave('CAGEE/global_figures/n_gamma/multi_sigma_sim/Percent of genes per sigma category.png')

#### create tree from newick file ####
### load tree 
tree_10sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

png('CAGEE/Tree_10_species_birds_edit.png')
plot(tree_10sp)
edgelabels(cex = 0.8)
nodelabels(cex = 0.8)
dev.off()

# regular tree
p = tree_10sp %>%
  ggtree(size = 4) +
  # geom_nodelab(aes(label = paste0("-",
  #                                 Decrease)),
  #              geom = 'text',
  #              node = 'all',
  #              nudge_x = -3,
  #              # hjust = 1,
  #              nudge_y = -0.25,
  #              color = 'darkred',
  #              size = 12) +
  # geom_tiplab(aes(label = species),
  #             size = 16) +
theme_tree2() +
  scale_x_continuous(labels = abs) +
  theme(axis.text.x = element_text(size = 18))
revts(p)
ggsave('CAGEE/global_figures/tree/tree poster.pdf',
       height = 15,
       width = 4)

# regular tree
# labels
p = tree_10sp %>%
  ggtree(size = 4) +
theme_tree2() +
  geom_tiplab() +
  scale_x_continuous(labels = abs) +
  theme(axis.text.x = element_text(size = 18))
revts(p)
ggsave('CAGEE/global_figures/tree/tree poster labels.pdf',
       height = 15,
       width = 4)

#### create poster graphic of sigma results ####
### load in data
data.poster = read.csv('CAGEE/CAGEE_models_summary.csv')

#### compare Z vs autosomes
## graph sigma for sex across chromsomes
# median
p1 = data.poster %>%
  filter(model == 1) %>% 
  filter(data.type == 'median') %>% 
  filter(type %in% c('A',
                     'Z',
                     "MA",
                     "MZ",
                     "FA",
                     "FZ")) %>% 
  mutate(sigma.type = case_when(
    str_detect(type, "M") ~ 'males',
    str_detect(type, "F") ~ 'females',
    TRUE ~ "Overall"
  )) %>%  
  mutate(chr = case_when(
    str_detect(type, "A") ~ 'A',
    str_detect(type, "Z") ~ 'Z',
  )) %>% 
  mutate(sigma.type = fct_relevel(sigma.type,
                                  "males",
                                  "Overall",
                                  "females")) %>% 
  ggplot(aes(x = chr,
             y= sigma,
             shape =sigma.type,
             color = chr)) +
  geom_point(size = 5,
             position=position_dodge(width = .5)) +
  theme_classic() +
  theme(axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24),
        axis.title.x =  element_text(size = 28),
        axis.text.x =  element_text(size = 24)) +
  ggtitle('Species median') +
  scale_color_manual(values = c("black",
                                "red")) +
  scale_shape_manual(values = c(15,
                                16,
                                17)) +
  expand_limits(y = 0)

# sex ratio
p2 = data.poster %>%
  filter(model == 1) %>% 
  filter(data.type != 'median') %>%
  filter(type != 'All') %>%
  na.omit() %>% 
  mutate(chr = case_when(
    str_detect(type, "A") ~ 'A',
    str_detect(type, "Z") ~ 'Z',
  )) %>% 
  ggplot(aes(x = chr,
             y= sigma,
             color = chr)) +
  geom_point(size = 5) +
  theme_classic() +
  theme(axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24),
        axis.title.x =  element_text(size = 28),
        axis.text.x =  element_text(size = 24))+
  ggtitle('M:F') +
  scale_color_manual(values = c("black",
                                "red")) +
  theme(legend.position = "none") +
    expand_limits(y = 0)
# graph together
p2+p1
ggsave('CAGEE/global_figures/Z vs A across sex sigma poster.pdf',
       width = 10,
       height = 5)

#### all genes
### median gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             # group = model
             )) +
  # geom_line(linewidth = 2) +
  geom_point(size = 5,
             aes(color = color)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0)
ggsave('CAGEE/global_figures/cavity models median gene expression poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             # group = model
  )) +
  # geom_line(linewidth = 2) +
  geom_point(size = 5,
             aes(color = color)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0)+
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') 
ggsave('CAGEE/global_figures/cavity models median gene expression facultative.pdf',
       width = 6.5,
       height = 5)

### ratio gene expression
## compare sigma rates across models 
## create table with sigma values for each model
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             # group = model
  )) +
  # geom_line(linewidth = 2) +
  geom_point(size = 5,
             aes(color = color)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0)
ggsave('CAGEE/global_figures/cavity models ratio gene expression poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             # group = model
  )) +
  # geom_line(linewidth = 2) +
  geom_point(size = 5,
             aes(color = color)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0)+
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') 
ggsave('CAGEE/global_figures/cavity models ratio gene expression facultative.pdf',
       width = 6.5,
       height = 5)

#### Male vs female all genes
### median gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("M",
                     "F",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                           "M",
                           "All",
                           "F")) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(15,
                                16,
                                17)) 
ggsave('CAGEE/global_figures/cavity models median gene expression sex poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("M",
                     "F",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "M",
                            "All",
                            "F")) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(15,
                                16,
                                17)) 
ggsave('CAGEE/global_figures/cavity models median gene expression sex facultative.pdf',
       width = 6.5,
       height = 5)


#### Z vs autosomes
### median gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("A",
                     "Z",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All",
                            "Z")) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16,
                                90)) 
ggsave('CAGEE/global_figures/cavity models median gene expression Z poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("A",
                     "Z",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All",
                            "Z")) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16,
                                90)) 
ggsave('CAGEE/global_figures/cavity models median gene expression Z facultative.pdf',
       width = 6.5,
       height = 5)

### ratio gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("A",
                     "Z",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All",
                            "Z")) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16,
                                90)) 
ggsave('CAGEE/global_figures/cavity models ratio gene expression Z poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("A",
                     "Z",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All",
                            "Z")) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16,
                                90)) 
ggsave('CAGEE/global_figures/cavity models ratio gene expression Z facultative.pdf',
       width = 6.5,
       height = 5)



#### No Z just autosomes
### median gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("A",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All")) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16)) 
ggsave('CAGEE/global_figures/cavity models median gene expression no Z poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type == 'median') %>% 
  filter(type %in% c("A",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All")) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16)) 
ggsave('CAGEE/global_figures/cavity models median gene expression no Z facultative.pdf',
       width = 6.5,
       height = 5)

### ratio gene expression
## compare sigma rates across models 
## graph sigma for every parameter in each model
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("A",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All")) %>% 
  filter(model != 4) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16)) 
ggsave('CAGEE/global_figures/cavity models ratio gene expression no Z poster.pdf',
       width = 6,
       height = 4.5)

# facultative
data.poster %>%
  filter(data.type != 'median') %>% 
  filter(type %in% c("A",
                     "All"))  %>% 
  mutate(color = case_when(model == 3 & model.pos == 2 ~ 'cavity',
                           model == 3 & model.pos == 3 ~ 'flexible',
                           model == 4 & model.pos == 2 ~ 'cavity',
                           model == 4 & model.pos == 3 ~ 'flexible',
                           model == 2 & model.pos == 2 ~ 'split',
                           TRUE ~ 'internal')) %>% 
  mutate(type = fct_relevel(type,
                            "A",
                            "All")) %>% 
  ggplot(aes(x = parameter,
             y= sigma,
             shape = type,
             group = type)) +
  geom_point(size = 5,
             aes(color = color),
             position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =  element_text(size = 28),
        axis.text.y =  element_text(size = 24)) +
  scale_x_continuous(breaks = seq(0,
                                  9,
                                  by = 1)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 3.5,
             linetype = 'dashed') +
  geom_vline(xintercept = 6.5,
             linetype = 'dashed') +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black',
                                'darkgrey')) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(65,
                                16)) 
ggsave('CAGEE/global_figures/cavity models ratio gene expression Z facultative.pdf',
       width = 6.5,
       height = 5)

#### graph sample gamma categories presentation ####
### graph example normal distributions of sigma 
### show multiple rate categories by dividing data
### set up distribution
delta     <- 0.001 
z.df     <- data.frame(x = seq(from=-3, to=3, by=delta))
z.df$density <- dnorm(z.df$x)
z.df$qt1  <- 1
z.df$qt3  <- cut(pnorm(z.df$x),breaks=3,labels=F)
z.df$qt10  <- cut(pnorm(z.df$x),breaks=10,labels=F)

## k = 1
# just distribution
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt1,fill=qt1),color="black", fill = 'white')+
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k1 distribution.png')

# color quantiles
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt1,fill=qt1),color="black", fill = 'darkgrey')+
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k1 distribution color.png')

# add mean value
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt1,fill=qt1),color="black", fill = 'darkgrey')+
  theme_classic() +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) 
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k1 distribution color mean.png')


## k = 3
# color quantiles
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt3,fill=qt3),color="black")+
  scale_fill_gradient2(midpoint=median(unique(z.df$qt3)), guide="none",
                       low = 'grey25',
                       mid = 'grey',
                       high = 'grey25') +
  theme_classic() 
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k3 distribution color.png')

# add mean value
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt3,fill=qt3),color="black")+
  theme_classic() +
  scale_fill_gradient2(midpoint=median(unique(z.df$qt3)), guide="none",
                       low = 'grey25',
                       mid = 'grey',
                       high = 'grey25') +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt3 == 1) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt3 == 2) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt3 == 3) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) 
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k3 distribution color mean.png')

## k = 10
# color quantiles
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt10,fill=qt10),color="black")+
  scale_fill_gradient2(midpoint=median(unique(z.df$qt10)), guide="none",
                       low = 'grey25',
                       mid = 'grey',
                       high = 'grey25') +
  theme_classic() +
  geom_segment(aes(x = 0,
                   y = 0,
                   yend = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 xend = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(x) %>% 
                   median()))
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k10 distribution color.png')

# add mean value
ggplot(z.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,y=density,group=qt10,fill=qt10),color="black")+
  theme_classic() +
  scale_fill_gradient2(midpoint=median(unique(z.df$qt10)), guide="none",
                       low = 'grey25',
                       mid = 'grey',
                       high = 'grey25') +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 1) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 2) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 3) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 4) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 5) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 6) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 7) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 8) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 9) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_point(aes(y = z.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(density) %>% 
                   max(),
                 x = z.df %>% 
                   filter(qt10 == 10) %>% 
                   pull(x) %>% 
                   median()),
             size = 3) +
  geom_segment(aes(x = 0,
                   y = 0,
                   yend = z.df %>% 
                     filter(qt1 == 1) %>% 
                     pull(density) %>% 
                     max(),
                   xend = z.df %>% 
                     filter(qt1 == 1) %>% 
                     pull(x) %>% 
                     median()))
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k10 distribution color mean.png')


### graph example gamma distributions of sigma 
### show multiple rate categories by dividing data
### set up distribution
delta     <- 0.01 
zg.df     <- data.frame(x = seq(from=0, to=10, by=delta))
zg.df$density <- dgamma(zg.df$x,
                        shape = 2,
                        scale = 2)
zg.df$qt1  <- 1
# zg.df$qt3  <- cut(pnorm(zg.df$x),breaks=3,labels=F)
zg.df$qt10  <- cut(pgamma(zg.df$x,
                          shape = 2,
                          scale = 2),
                   breaks=10,
                   labels=F)
zg.df$x = zg.df$x/200

## k = 1
# color quantiles
# add mean value
ggplot(zg.df,
       aes(x=x,
           y=density)) + 
  geom_area(aes(x=x,
                y=density,
                # group=qt1,
                # fill=qt1
                ),
            color="black", 
            fill = 'darkgrey')+
  theme_classic() +
  geom_point(aes(y = 0.011,
                 x = (zg.df %>% 
                   filter(qt1 == 1) %>% 
                   pull(x) %>% 
                   median())),
             size = 1,
             shape = 25,
             fill = 'black')  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab('sigma') +
  xlim(0,0.05)+
  ylim(0,0.19) +
  theme_classic(base_size = 8)+
  theme(axis.text = element_text(colour = "black"))
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k1 distribution color mean gamma.pdf',
       width = 2.2,
       height = 1.2,
       units = 'in',
       dpi = 300)


## k = 10
# color quantiles
# add mean value
ggplot(zg.df,aes(x=x,y=density)) + 
  geom_area(aes(x=x,
                y=density,
                ),
            fill = 'darkgrey',
            color = 'black'
            )+
  theme_classic() +
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 1) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 1) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)),
               linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 2) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 2) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 3) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 3) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 4) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 4) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 5) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 5) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 6) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 6) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 7) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 7) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 8) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 8) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 9) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 9) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
  geom_segment(aes(y = 0,
                   yend = zg.df %>% 
                     filter(qt10 == 10) %>% 
                     filter(x == max(x)) %>% 
                     pull(density),
                   x = zg.df %>% 
                     filter(qt10 == 10) %>% 
                     filter(x == max(x)) %>% 
                     pull(x)), linewidth = 0.25)+
    geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 1) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 2) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 3) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 4) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 5) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 6) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 7) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 8) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 9) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  geom_point(aes(y = 0.011,
                 x = zg.df %>% 
                   filter(qt10 == 10) %>% 
                   pull(x) %>% 
                   median()),
             size = 1,
             shape = 25,
             fill = 'black') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab('sigma') +
  xlim(0,0.05) +
  ylim(0,0.19) +
  theme_classic(base_size = 8)+
  theme(axis.text = element_text(colour = "black"))
ggsave('CAGEE/global_figures/n_gamma/example_presentation/k10 distribution color mean gamma.pdf',
       width = 2.2,
       height = 1.2,
       units = 'in',
       dpi = 720)


#### Credible changes ratio ####
### want to compare genes that changed in sex ratio
## load data
# get credible changes for each gene
data.cagee.changes.ratio = read.delim('CAGEE/ratio/ratio_cagee_outs/change.tab')

## convert data so that each gene is indicated as change (1) or not change (0)
# credible change genes have "*" 
data.cagee.changes.ratio.sig = data.cagee.changes.ratio %>% 
  column_to_rownames('TranscriptID') %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "\\*",
                                     negate = TRUE),
                          0))) %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "\\*"),
                          1))) %>% 
  mutate_all((function(x) as.numeric(as.character(x)))) %>% 
  as.matrix()

# get direction of change by assigning -1 and 1 to each value
data.cagee.changes.ratio.direction = data.cagee.changes.ratio %>% 
  column_to_rownames('TranscriptID') %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "-",
                                     negate = TRUE),
                          1))) %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "-"),
                          -1))) %>% 
  mutate_all((function(x) as.numeric(as.character(x)))) %>% 
  as.matrix()

# combine significance and direction
data.cagee.changes.ratio.direction.sig = data.cagee.changes.ratio.direction*data.cagee.changes.ratio.sig

## get summary for each gene
# number of credible changes 
# number of changes that are positive or negative
data.cagee.changes.ratio.sum = data.cagee.changes.ratio.sig %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  dplyr::rename(Num_credible_changes = '.') %>% 
  rownames_to_column('TranscriptID') %>% 
  full_join(data.cagee.changes.ratio.direction.sig %>% 
              as.data.frame() %>% 
              mutate(
                positive = rowSums(across(everything(), ~ .x == 1)),
                negative = rowSums(across(everything(), ~ .x == -1)),
              ) %>% 
              select(positive,
                     negative)%>% 
              rownames_to_column('TranscriptID')) %>% 
  mutate(concordance = case_when(Num_credible_changes == positive ~ 'concordant',
                                 Num_credible_changes == negative ~ 'concordant',
                                 TRUE ~ 'discordant'))

# number of genes with changes
# 5308
data.cagee.changes.ratio.sum %>% 
  filter(Num_credible_changes >= 1) %>% 
  nrow()

# number of genes with multiple changes
# 3244
data.cagee.changes.ratio.sum %>% 
  filter(Num_credible_changes >= 2) %>% 
  nrow()

# number of genes with multiple changes
# 266
data.cagee.changes.ratio.sum %>% 
  filter(Num_credible_changes >= 2) %>% 
  filter(concordance == 'concordant') %>% 
  nrow()

### graph number of changes
# data.cagee.changes.ratio.sum %>% 
#   group_by(Num_credible_changes) %>% 
#   summarise(Total_credible_changes = n()) %>% 
#   ggplot(aes(y = Total_credible_changes,
#              x = as.factor(Num_credible_changes),
#              label = Total_credible_changes)) +
#   geom_vline(xintercept = 1.5,
#              linetype = 'dashed') +
#   geom_bar(stat = 'identity') +
#   geom_text(vjust = -0.5) +
#   theme_classic() +
#   ggtitle("Number of genes with credible change in M:F across branches") +
#   xlab('Number of branches') +
#   ylab('Number of genes')
# ggsave("CAGEE/global_figures/credible_change/ratio credible changes hist.png")

# graph concordance
data.cagee.changes.ratio.sum %>% 
  group_by(Num_credible_changes,
           concordance) %>% 
  summarise(Total_credible_changes = n()) %>% 
  mutate(concordance = ifelse(Num_credible_changes <= 1,
                              'NA',
                              concordance)) %>% 
  ggplot(aes(y = Total_credible_changes,
             x = as.factor(Num_credible_changes),
             label = Total_credible_changes,
             fill = concordance)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_bar(stat = 'identity',
           position = 'dodge',
           color = 'black') +
  geom_text(vjust = -0.5,
            position = position_dodge(width = 1)) +
  theme_classic() +
  ggtitle("Number of genes with credible change in M:F across branches") +
  xlab('Number of branches') +
  ylab('Number of genes') +
  scale_fill_manual(values = c('grey',
                               'black',
                               'white'))
ggsave("CAGEE/global_figures/credible_change/ratio credible changes concordance hist.png")

# paper
data.cagee.changes.ratio.sum %>% 
  group_by(Num_credible_changes,
           concordance) %>% 
  summarise(Total_credible_changes = n()) %>% 
  mutate(concordance = ifelse(Num_credible_changes <= 1,
                              'NA',
                              concordance)) %>% 
  ggplot(aes(y = Total_credible_changes,
             x = as.factor(Num_credible_changes),
             label = Total_credible_changes,
             fill = concordance)) +
  geom_vline(xintercept = 1.5,
             linetype = 'dashed') +
  geom_bar(stat = 'identity',
           position = 'dodge',
           color = 'black') +
  geom_text(vjust = -0.5,
            position = position_dodge(width = 1)) +
  # ggtitle("Number of genes with credible change in M:F across branches") +
  xlab('Number of branches') +
  ylab('Number of genes') +
  scale_fill_manual(values = c('grey',
                               'black',
                               'white')) +
  theme_classic(base_size = 8) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.90,
                                   0.5))
ggsave("CAGEE/global_figures/credible_change/ratio credible changes concordance hist.pdf",
       height = 4,
       width = 6.5,
       dpi = 320,
       units = 'in')



## summarize across cavity and flexible nesting 
# filter to external branches
data.cagee.changes.ratio.direction.sig.cavity =
  data.cagee.changes.ratio.direction.sig %>% 
    as.data.frame() %>% 
  rownames_to_column('TranscriptID') %>% 
  pivot_longer(cols = -TranscriptID,
               names_to = 'branch',
               values_to = 'direction.sig') %>% 
  filter(!str_detect(branch,
                    'X.')) %>% 
    separate_wider_delim(branch,
                         delim = '.',
                         names = c('species',
                                   NA,
                                   NA)) %>% 
    mutate(cavity = case_when(species %in% c("ET","PW","BB","HW","TS") ~ 'obligate',
                              species %in% c("HS","YW","RO","CW","BS") ~ 'flexible',
                              TRUE ~ NA)) %>% 
    select(TranscriptID,
             cavity,
           direction.sig) %>% 
    table() %>% 
    as.data.frame() %>% 
    mutate(direction.sig = case_when(direction.sig == -1 ~ 'negative',
                                     direction.sig == 1 ~ 'positive',
                                     direction.sig == 0 ~ 'none')) %>% 
    pivot_wider(id_cols = c('TranscriptID',
                            'cavity'),
                names_from = 'direction.sig',
                values_from = 'Freq') %>% 
    mutate(concordance = case_when(negative > 0 & positive > 0 ~ 'discordant',
                                   negative > 1 & positive == 0 ~ 'neg.concordant',
                                   negative == 0 & positive > 1 ~ 'pos.concordant',
                                   # negative == 2 & positive == 0 ~ 'neg.concordant.2',
                                   # negative == 0 & positive == 2 ~ 'pos.concordant.2',
                                   # negative > 2 & positive == 0 ~ 'neg.concordant.3',
                                   # negative == 0 & positive > 2 ~ 'pos.concordant.3',
                                   negative == 1 & positive == 0 ~ 'neg.single',
                                   negative == 0 & positive == 1 ~ 'pos.single',
                                   TRUE ~ 'none'),
           discordant.ratio = paste0(as.character(positive),'/',as.character(negative))) 

## graph
# concordance
data.cagee.changes.ratio.direction.sig.cavity %>% 
  group_by(cavity,
           concordance) %>% 
  summarise(count = n()) %>% 
ggplot(aes(x = reorder(concordance,
                       -count),
           y = count,
           group = cavity)) +
  geom_bar(stat = 'identity',
           position = position_dodge2('single', 
                                      width = 0.9),
           aes(fill = cavity)) +
  theme_classic() +
  geom_label(aes(label = count),
             vjust = -0.05,
             position = position_dodge2('single', 
                                        width = 0.9)) +
  xlab('credible changes ratio')
ggsave('CAGEE/global_figures/credible_change/ratio credible changes concordance cavity.png')



# discordance ratio
data.cagee.changes.ratio.direction.sig.cavity %>% 
  group_by(cavity,
           discordant.ratio) %>% 
  summarise(count = n()) %>% 
  filter(discordant.ratio != '0/0') %>% 
  ggplot(aes(x = as.character(discordant.ratio),
             y = count,
             group = cavity)) +
  geom_bar(stat = 'identity',
           position = position_dodge2('single', 
                                      width = 0.9),
           aes(fill = cavity)) +
  theme_classic() +
  geom_label(aes(label = count),
             vjust = -0.05,
             position = position_dodge2('single', 
                                        width = 0.9)) +
  xlab('Discordance ratio (pos/neg)')+
  theme_classic()
ggsave('CAGEE/global_figures/credible_change/ratio credible changes discordance cavity.png')


















#### Credible changes median ####
### want to compare genes that changed in sex ratio
## load data
# get credible changes for each gene
data.cagee.changes.median = read.delim('CAGEE/median_cavity/median_cavity_cagee_outs/change.tab')

## convert data so that each gene is indicated as change (1) or not change (0)
# credible change genes have "*" 
data.cagee.changes.median.sig = data.cagee.changes.median %>% 
  column_to_rownames('TranscriptID') %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "\\*",
                                     negate = TRUE),
                          0))) %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "\\*"),
                          1))) %>% 
  mutate_all((function(x) as.numeric(as.character(x)))) %>% 
  as.matrix()

# get direction of change by assigning -1 and 1 to each value
data.cagee.changes.median.direction = data.cagee.changes.median %>% 
  column_to_rownames('TranscriptID') %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "-",
                                     negate = TRUE),
                          1))) %>% 
  mutate(across(everything(),
                ~ replace(.,  
                          str_detect(.,
                                     "-"),
                          -1))) %>% 
  mutate_all((function(x) as.numeric(as.character(x)))) %>% 
  as.matrix()

# combine significance and direction
data.cagee.changes.median.direction.sig = data.cagee.changes.median.direction*data.cagee.changes.median.sig

## get summary for each gene
# number of credible changes 
# number of changes that are positive or negative
data.cagee.changes.median.sum = data.cagee.changes.median.sig %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  dplyr::rename(Num_credible_changes = '.') %>% 
  rownames_to_column('TranscriptID') %>% 
  full_join(data.cagee.changes.median.direction.sig %>% 
              as.data.frame() %>% 
              mutate(
                positive = rowSums(across(everything(), ~ .x == 1)),
                negative = rowSums(across(everything(), ~ .x == -1)),
              ) %>% 
              select(positive,
                     negative)%>% 
              rownames_to_column('TranscriptID')) %>% 
  mutate(concordance = case_when(Num_credible_changes == positive ~ 'concordant',
                                 Num_credible_changes == negative ~ 'concordant',
                                 TRUE ~ 'discordant'))

# number of genes with changes
# 2779
data.cagee.changes.median.sum %>% 
  filter(Num_credible_changes >= 1) %>% 
  nrow()

# number of genes with multiple changes
# 2211
data.cagee.changes.median.sum %>% 
  filter(Num_credible_changes >= 2) %>% 
  nrow()

# number of genes with multiple changes
# 49
data.cagee.changes.median.sum %>% 
  filter(Num_credible_changes >= 2) %>% 
  filter(concordance == 'concordant') %>% 
  nrow()

### graph number of changes
# data.cagee.changes.median.sum %>% 
#   group_by(Num_credible_changes) %>% 
#   summarise(Total_credible_changes = n()) %>% 
#   ggplot(aes(y = Total_credible_changes,
#              x = as.factor(Num_credible_changes),
#              label = Total_credible_changes)) +
#   geom_vline(xintercept = 1.5,
#              linetype = 'dashed') +
#   geom_bar(stat = 'identity') +
#   geom_text(vjust = -0.5) +
#   theme_classic() +
#   ggtitle("Number of genes with credible change in M:F across branches") +
#   xlab('Number of branches') +
#   ylab('Number of genes')
# ggsave("CAGEE/global_figures/credible_change/median credible changes hist.png")

# graph concordance
# data.cagee.changes.median.sum %>% 
#   group_by(Num_credible_changes,
#            concordance) %>% 
#   summarise(Total_credible_changes = n()) %>% 
#   mutate(concordance = ifelse(Num_credible_changes <= 1,
#                               'NA',
#                               concordance)) %>% 
#   ggplot(aes(y = Total_credible_changes,
#              x = as.factor(Num_credible_changes),
#              label = Total_credible_changes,
#              fill = concordance)) +
#   geom_vline(xintercept = 1.5,
#              linetype = 'dashed') +
#   geom_bar(stat = 'identity',
#            position = 'dodge',
#            color = 'black') +
#   geom_text(vjust = -0.5,
#             position = position_dodge(width = 1)) +
#   theme_classic() +
#   ggtitle("Number of genes with credible change in M:F across branches") +
#   xlab('Number of branches') +
#   ylab('Number of genes') +
#   scale_fill_manual(values = c('grey',
#                                'black',
#                                'white'))
# ggsave("CAGEE/global_figures/credible_change/median credible changes concordance hist.png")

## summarize across cavity and flexible nesting 
# filter to external branches
data.cagee.changes.median.direction.sig.cavity =
  data.cagee.changes.median.direction.sig %>% 
  as.data.frame() %>% 
  rownames_to_column('TranscriptID') %>% 
  pivot_longer(cols = -TranscriptID,
               names_to = 'branch',
               values_to = 'direction.sig') %>% 
  filter(!str_detect(branch,
                     'X.')) %>% 
  separate_wider_delim(branch,
                       delim = '.',
                       names = c('species',
                                 NA,
                                 NA)) %>% 
  mutate(cavity = case_when(species %in% c("ET","PW","BB","HW","TS") ~ 'obligate',
                            species %in% c("HS","YW","RO","CW","BS") ~ 'flexible',
                            TRUE ~ NA)) %>% 
  select(TranscriptID,
         cavity,
         direction.sig) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(direction.sig = case_when(direction.sig == -1 ~ 'negative',
                                   direction.sig == 1 ~ 'positive',
                                   direction.sig == 0 ~ 'none')) %>% 
  pivot_wider(id_cols = c('TranscriptID',
                          'cavity'),
              names_from = 'direction.sig',
              values_from = 'Freq') %>% 
  mutate(concordance = case_when(negative > 0 & positive > 0 ~ 'discordant',
                                 negative > 1 & positive == 0 ~ 'neg.concordant',
                                 negative == 0 & positive > 1 ~ 'pos.concordant',
                                 # negative == 2 & positive == 0 ~ 'neg.concordant.2',
                                 # negative == 0 & positive == 2 ~ 'pos.concordant.2',
                                 # negative > 2 & positive == 0 ~ 'neg.concordant.3',
                                 # negative == 0 & positive > 2 ~ 'pos.concordant.3',
                                 negative == 1 & positive == 0 ~ 'neg.single',
                                 negative == 0 & positive == 1 ~ 'pos.single',
                                 TRUE ~ 'none'),
         discordant.ratio = paste0(as.character(positive),'/',as.character(negative))) 

## graph
# concordance
data.cagee.changes.median.direction.sig.cavity %>% 
  group_by(cavity,
           concordance) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = reorder(concordance,
                         -count),
             y = count,
             group = cavity)) +
  geom_bar(stat = 'identity',
           position = position_dodge2('single', 
                                      width = 0.9),
           aes(fill = cavity)) +
  theme_classic() +
  geom_label(aes(label = count),
             vjust = -0.05,
             position = position_dodge2('single', 
                                        width = 0.9)) +
  xlab('credible changes median')
ggsave('CAGEE/global_figures/credible_change/median credible changes concordance cavity.png')



# discordance median
data.cagee.changes.median.direction.sig.cavity %>% 
  group_by(cavity,
           discordant.ratio) %>% 
  summarise(count = n()) %>% 
  filter(discordant.ratio != '0/0') %>% 
  ggplot(aes(x = as.character(discordant.ratio),
             y = count,
             group = cavity)) +
  geom_bar(stat = 'identity',
           position = position_dodge2('single', 
                                      width = 0.9),
           aes(fill = cavity)) +
  theme_classic() +
  geom_label(aes(label = count),
             vjust = -0.05,
             position = position_dodge2('single', 
                                        width = 0.9)) +
  xlab('Discordance median (pos/neg)')+
  theme_classic()
ggsave('CAGEE/global_figures/credible_change/median credible changes discordance cavity.png')


## compare median and ratio
# need to run ratio section before this
data.cagee.changes.medianvsratio = 
  data.cagee.changes.median.direction.sig.cavity %>% 
  dplyr::select(TranscriptID,
                cavity,
                concordance) %>% 
  dplyr::rename(concordance.median = concordance) %>% 
  full_join(data.cagee.changes.ratio.direction.sig.cavity %>% 
              dplyr::select(TranscriptID,
                            cavity,
                            concordance) %>% 
              dplyr::rename(concordance.ratio = concordance))
  
## graph upset plot
library(ggupset)
# graph median
 p = data.cagee.changes.medianvsratio %>%
  filter(concordance.median %in% c('neg.concordant',
                                   'pos.concordant')) %>% 
  mutate(concordance.median = paste0(cavity,
                                     '_',
                                     concordance.median)) %>% 
  group_by(TranscriptID) %>%
  summarize(concordance.median = list(concordance.median)) %>%
  ggplot(aes(x = concordance.median)) +
  geom_bar() +
  geom_text(stat='count', 
            aes(label=after_stat(count)),
            vjust=1,
            color = 'white') +
  scale_x_upset() +
  theme_classic()
png('CAGEE/global_figures/credible_change/median credible changes cavity upset.png')
p
dev.off()

# graph ratio
p = data.cagee.changes.medianvsratio %>%
  filter(concordance.ratio %in% c('neg.concordant',
                                   'pos.concordant')) %>% 
  mutate(concordance.ratio = paste0(cavity,
                                     '_',
                                     concordance.ratio)) %>% 
  group_by(TranscriptID) %>%
  summarize(concordance.ratio = list(concordance.ratio)) %>%
  ggplot(aes(x = concordance.ratio)) +
  geom_bar() +
  geom_text(stat='count', 
            aes(label=after_stat(count)),
            vjust=1,
            color = 'white') +
  scale_x_upset() +
  theme_classic()
png('CAGEE/global_figures/credible_change/ratio credible changes cavity upset.png')
p
dev.off()
# graph ratio vs median
p = data.cagee.changes.medianvsratio %>%
  mutate(concordance.ratio = paste0(cavity,
                                    '_ratio_',
                                    concordance.ratio)) %>% 
  mutate(concordance.median = paste0(cavity,
                                    '_median_',
                                    concordance.median)) %>% 
  dplyr::select(-c(cavity)) %>% 
  pivot_longer(cols = -TranscriptID,
               names_to = 'type',
               values_to = 'concordance.all') %>% 
  filter(str_detect(concordance.all,
                    '.concordant')) %>% 
  group_by(TranscriptID) %>%
  summarize(concordance.all = list(concordance.all)) %>%
  ggplot(aes(x = concordance.all)) +
  geom_bar() +
  geom_text(stat='count', 
            aes(label=after_stat(count)),
            vjust=1,
            color = 'white') +
  scale_x_upset() +
  theme_classic()
png('CAGEE/global_figures/credible_change/all credible changes cavity upset.png')
p
dev.off()

# graph ratio vs median
# reduce
p = data.cagee.changes.medianvsratio %>%
  mutate(concordance.ratio = paste0(cavity,
                                    '_ratio_',
                                    concordance.ratio)) %>% 
  mutate(concordance.median = paste0(cavity,
                                     '_median_',
                                     concordance.median)) %>% 
  dplyr::select(-c(cavity)) %>% 
  pivot_longer(cols = -TranscriptID,
               names_to = 'type',
               values_to = 'concordance.all') %>% 
  filter(str_detect(concordance.all,
                    '.concordant')) %>% 
  group_by(TranscriptID) %>%
  summarize(concordance.all = list(concordance.all)) %>%
  ggplot(aes(x = concordance.all)) +
  geom_bar() +
  geom_text(stat='count', 
            aes(label=after_stat(count)),
            vjust=1,
            color = 'white') +
  scale_x_upset(n_intersections = 12) +
  theme_classic()





#### CAGEE results median ####
### load result data
## load clad results
data.cagee.clade.results = read.delim('CAGEE/median/median_cagee_outs_2//clade_results.tab')
# separate data
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

## load results
data.cagee.results = read.delim('CAGEE/median/median_cagee_outs_2/results.txt')
# separate data
# remove Sigma2 line
data.cagee.results = data.cagee.results %>% 
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%  
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
# change to get tree
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[5,2],
                                          ";",
                                          sep = ''))

## add edge length
# load species tree
tree_10sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

# add edge length
tree.cagee$edge.length = tree_10sp$edge.length

# convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

## graph
p = tree.cagee %>% 
  ggtree() +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               # hjust = 1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               # hjust = 1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE male to female median",
          subtitle = paste(data.cagee.results[3,1],
                           data.cagee.results[3,2] %>% 
                             as.numeric() %>% 
                             round(5),
                           "; ",
                           data.cagee.results[4,1],
                           data.cagee.results[4,2] %>% 
                             as.numeric() %>% 
                             round(5),
                           "; ",
                           "-lnl:",
                           data.cagee.results[2,2]))+ 
  scale_x_continuous(labels = abs) 
revts(p)
ggsave('CAGEE/median/median_cagee_outs_2/figures/tree clade results median.png',
       height = 10,
       width = 10)



#### CAGEE results median all ####
### load result data
## load clad results
data.cagee.clade.results = read.delim('CAGEE/median/median_cagee_outs_all/clade_results.tab')
# separate data
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

## load results
data.cagee.results = read.delim('CAGEE/median/median_cagee_outs_all/results.txt')
# separate data
# remove Sigma2 line
data.cagee.results = data.cagee.results %>% 
  # mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
  #                                    "Sigma2: ",
  #                                    "")) %>%  
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
# change to get tree
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[4,2],
                                          ";",
                                          sep = ''))

## add edge length
# load species tree
tree_10sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

# add edge length
tree.cagee$edge.length = tree_10sp$edge.length

# convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

## graph
p = tree.cagee %>% 
  ggtree() +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               # hjust = 1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               # hjust = 1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE all median",
          subtitle = paste(data.cagee.results[3,1],
                           data.cagee.results[3 ,2] %>% 
                             as.numeric() %>% 
                             round(5),
                           "; ",
                           "-lnl:",
                           data.cagee.results[2,2])) +  
  scale_x_continuous(labels = abs) 
revts(p)
ggsave('CAGEE/median/median_cagee_outs_all/figures/tree clade results median all.png',
       height = 10,
       width = 10)

# poster
# # to run ggimage need to start in geode2 terminal
# # run the following to open rstudio
# # >module load imagemagick/imagemagick-7.1.1
# # >module load rstudio
# # >rstudio
# library(ggimage)
# # create bar graphs for each branch of changed genes
# bars2 <- nodebar(as.data.frame(tree.cagee@data), 
#                  cols=4:5, 
#                  position='dodge')
# # add color, scale, and flip bar graphs
# bars3 <- lapply(bars2,
#                function(g) g+scale_fill_manual(values = c('darkred',
#                                                           'darkgreen'))+
#                  # scale_y_continuous(limits = c(0,
#                  #                               1800)) +
#                  scale_y_reverse() +
#                  coord_flip(ylim = c(1800,
#                                      0)) 
#                  # scale_x_reverse()
#                )
# # create tree
# p = tree.cagee %>%
#   ggtree(size = 2) +
#   # geom_nodelab(aes(label = paste0("+",
#   #                                 Increase)),
#   #              geom = 'text',
#   #              node = 'all',
#   #              nudge_x = -3,
#   #              # hjust = 1,
#   #              nudge_y = 0.25,
#   #              color = 'darkgreen',
#   #              size = 12) +
#   # geom_nodelab(aes(label = paste0("-",
#   #                                 Decrease)),
# #              geom = 'text',
# #              node = 'all',
# #              nudge_x = -3,
# #              # hjust = 1,
# #              nudge_y = -0.25,
# #              color = 'darkred',
# #              size = 12) +
# # geom_tiplab(aes(label = species),
# #             size = 16) +
# theme_tree2() +
#   ggtitle("CAGEE all median",
#           subtitle = paste(data.cagee.results[3,1],
#                            data.cagee.results[3 ,2] %>%
#                              as.numeric() %>%
#                              round(5),
#                            "; ",
#                            "-lnl:",
#                            data.cagee.results[2,2]))+
#   scale_x_continuous(labels = abs)
# revts(p)
# 
# # for each bar graph need to adjust starting position
# p + geom_inset(bars3[9:10], 
#           x='branch',
#           height = 0.05, 
#           width=0.25,
#           vjust = 0) + 
#   geom_inset(bars3[7:8], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -0.45,
#              vjust = 0)+ 
#   geom_inset(bars3[5:6], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -6.25,
#              vjust = .05)+ 
#   geom_inset(bars3[3:4], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -5.45,
#              vjust = .05)+ 
#   geom_inset(bars3[1:2], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -5.5,
#              vjust = .05)+ 
#   geom_inset(bars3[11:12], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -1.65,
#              vjust = .05)+ 
#   geom_inset(bars3[19], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -5.95,
#              vjust = .05)+ 
#   geom_inset(bars3[18], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -5.50,
#              vjust = .05)+ 
#   geom_inset(bars3[17], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = 1.90,
#              vjust = .05)+ 
#   geom_inset(bars3[16], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -0.095,
#              vjust = .05)+ 
#   geom_inset(bars3[15], 
#              x='branch',
#              height = 0.05, 
#              width=0.25,
#              hjust = -0.90,
#              vjust = 0.05)





#### CAGEE results cavity  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/median_cavity_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/median_cavity_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species_Data/data/Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/median_cavity_cagee_outs/figures/tree clade results median cavity.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)

# create bar graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  ),
  size = 5) +
  geom_rect(aes(xmin = x-Increase/150,
                xmax = x,
                ymin = y,
                ymax = y+0.25),
            fill = 'grey75',
            color = 'black') +
  geom_rect(aes(xmin = x-Decrease/150,
                xmax = x,
                ymin = y,
                ymax = y-0.25),
            fill = 'grey25',
            color = 'black') +
  # geom_nodelab(aes(label = paste0("-",
  #                                 Decrease)),
  #              geom = 'text',
  #              node = 'all',
  #              nudge_x = -3,
  #              # hjust = 1,
  #              nudge_y = -0.25,
  #              color = 'darkred',
  #              size = 12) +
  # geom_tiplab(aes(label = species),
  #             size = 16) +
theme_tree2() +
  ggtitle("CAGEE results cavity facultative",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result) ,
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) +
  scale_x_continuous(labels = abs) +
  theme(axis.text.x = element_text(size = 18),
        legend.position = "none")+
  scale_color_manual(values = c('#A96B51',
                                '#34B881',
                                'black'))
revts(p)  +
  scale_x_reverse() 
ggsave('CAGEE/global_figures/tree/tree clade results median cavity poster.pdf',
       height = 10,
       width = 7.5)



#### CAGEE results cavity facultative  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/median_cavity_facultative_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/median_cavity_facultative_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_facultative.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity facultative",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/median_cavity_facultative_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)





#### CAGEE results cavity 2p  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/median_cavity_2p_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/median_cavity_2p_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_2p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity 2p",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/median_cavity_2p_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)





#### CAGEE results cavity 1p  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/median_cavity_1p_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/median_cavity_1p_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_1p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity 1p",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/median_cavity_1p_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)


# create bar graph
p = tree.cagee %>%
  ggtree(size = 4) +
  geom_rect(aes(xmin = x-Increase/150,
                xmax = x,
                ymin = y,
                ymax = y+0.25),
            fill = 'grey75',
            color = 'black') +
  geom_rect(aes(xmin = x-Decrease/150,
                xmax = x,
                ymin = y,
                ymax = y-0.25),
            fill = 'grey25',
            color = 'black') +
  # geom_nodelab(aes(label = paste0("-",
  #                                 Decrease)),
  #              geom = 'text',
  #              node = 'all',
  #              nudge_x = -3,
  #              # hjust = 1,
  #              nudge_y = -0.25,
  #              color = 'darkred',
  #              size = 12) +
  # geom_tiplab(aes(label = species),
  #             size = 16) +
theme_tree2() +
  ggtitle("CAGEE all median",
          subtitle = paste(data.cagee.results[3,1],
                           data.cagee.results[3 ,2] %>%
                             as.numeric() %>%
                             round(5),
                           "; ",
                           "-lnl:",
                           data.cagee.results[2,2]))+
  scale_x_continuous(labels = abs) +
  theme(axis.text.x = element_text(size = 18))
revts(p)
ggsave('CAGEE/global_figures/tree/tree clade results median all poster.pdf',
       height = 10,
       width = 7.5)


#### CAGEE results ratio ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

# ### Optional: add sigma tree to CAGEE tree
# ## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# ## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    # colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/ratio/ratio_cagee_outs/figures/tree clade results ratio.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)


#### CAGEE results cavity ratio  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_cavity_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_cavity_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity sex ratio",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/ratio/ratio_cavity_cagee_outs/figures/tree clade results ratio cavity.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)





#### CAGEE results cavity facultative ratio  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_cavity_facultative_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_cavity_facultative_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_facultative.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity facultative sex ratio",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/ratio/ratio_cavity_facultative_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)





#### CAGEE results cavity 2p ratio  ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_cavity_2p_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_cavity_2p_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_2p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

## load in species 
data.species.cheat = readxl::read_xlsx('Species cheat sheet.xlsx')

# rename KYTS to TS
# add cavity nesting label
data.species.cheat = data.species.cheat %>% 
  mutate(Abbreviation = str_replace(Abbreviation,
                                    "KYTS",
                                    "TS")) 

## add branch colors for the 3 sigmas
tree.cagee = tree.cagee %>% 
  full_join(data.species.cheat %>% 
              dplyr::select(Abbreviation,
                            `Nesting strategy`),
            by = join_by(species == Abbreviation))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species,
                  color = `Nesting strategy`)) +
  theme_tree2() +
  ggtitle("CAGEE results cavity 2p sex ratio",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/ratio/ratio_cavity_2p_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

### note for mult-itissue analysis need to add sigma values as table separately 
# ## create sigma table
# sigma.table = data.cagee.results %>% 
#   filter(grepl("Sigma",
#                CAGEE.1.1)) 
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>% 
#   ggtree(aes(colour = sigma)) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>% 
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                              pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# multiplot(revts(p), 
#           p1, 
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)





#### CAGEE results median z chromosome ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_z/median_z_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_z/median_z_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

# note: for multi-tissue analysis run this instead
# see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

# ### Optional: add sigma tree to CAGEE tree
# ## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# ## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
# ### create tree graph
# p = tree.cagee %>% 
#   ggtree(aes(
#     colour = sigma #Optional
#   )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle =  paste(data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Sigma2") %>% 
#                               pull(CAGEE.1.1),
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Sigma2") %>% 
#                               pull(result),
#                             ";",
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                               pull(CAGEE.1.1),
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                               pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# revts(p)
# 
# ## save graph
# ggsave('CAGEE/median_z/median_z_cagee_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)

## note for mult-itissue analysis need to add sigma values as table separately
## create sigma table
sigma.table = data.cagee.results %>%
  filter(grepl("Sigma",
               CAGEE.1.1))

# graph sigma table
library(ggpubr)
p1 = ggpubr::ggtexttable(sigma.table,
                 rows = NULL)

## graph CAGEE results
p = tree.cagee %>%
  ggtree(aes(
    # colour = sigma
             )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle = paste(data.cagee.results %>%
                             filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
                             pull(CAGEE.1.1),
                           data.cagee.results %>%
                             filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
                             pull(result))) +
  scale_x_continuous(labels = abs)

## reorder timeline
multiplot(revts(p),
          p1,
          widths = c(8,2),
          ncol = 2)

## save graph
ggsave('CAGEE/median_z/median_z_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)







#### CAGEE results median 2p 4 gamma cats ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median/median_cagee_outs_2p_n_gamma_cats/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median/median_cagee_outs_2p_n_gamma_cats/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_2p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median/median_cagee_outs_2p_n_gamma_cats/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)








#### CAGEE results ratio ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_z_cagee_outs//clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_z_cagee_outs//results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

# note: for multi-tissue analysis run this instead
# see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

# ### Optional: add sigma tree to CAGEE tree
# ## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# ## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

# #### graph tree with results 
# ### create tree graph
# p = tree.cagee %>% 
#   ggtree(aes(
#     # colour = sigma #Optional
#   )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle =  paste(data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Sigma2") %>% 
#                               pull(CAGEE.1.1),
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Sigma2") %>% 
#                               pull(result),
#                             ";",
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                               pull(CAGEE.1.1),
#                             data.cagee.results %>% 
#                               filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
#                               pull(result))) + 
#   scale_x_continuous(labels = abs) 
# 
# ## reorder timeline
# revts(p)
# 
# ## save graph
# ggsave('CAGEE/ratio/ratio_cagee_outs/figures/tree clade results ratio.png',
#        height = 10,
#        width = 10)

## note for mult-itissue analysis need to add sigma values as table separately
## create sigma table
sigma.table = data.cagee.results %>%
  filter(grepl("Sigma",
               CAGEE.1.1))

# graph sigma table
library(ggpubr)
p1 = ggpubr::ggtexttable(sigma.table,
                 rows = NULL)

## graph CAGEE results
p = tree.cagee %>%
  ggtree(aes(
    # colour = sigma
    )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle = paste(data.cagee.results %>%
                             filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
                             pull(CAGEE.1.1),
                           data.cagee.results %>%
                             filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
                             pull(result))) +
  scale_x_continuous(labels = abs)

## reorder timeline
p2 = ggarrange(revts(p),
          p1,
          widths = c(8,2),
          ncol = 2,
          nrow = 1)
p2
## save graph
ggsave('CAGEE/ratio/ratio_z_cagee_outs//figures/tree clade results.png',
       p2,
       height = 10,
       width = 10)

#### CAGEE results cavity 4p ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/median_cavity_4p_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/median_cavity_4p_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_4p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/median_cavity_4p_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)








#### CAGEE results cavity ratio 4p ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/ratio/ratio_cavity_4p_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/ratio/ratio_cavity_4p_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity_4p.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/ratio/ratio_cavity_4p_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)








#### CAGEE results freerate ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median/median_freerate_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median/median_freerate_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# works with free rate model
data.cagee.results = data.cagee.results %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## Optional: add free rate 
tree.cagee = left_join(tree.cagee,
                       data.cagee.results %>% 
                         mutate(CAGEE.1.1 = as.integer(CAGEE.1.1)) %>% 
  rename(node = CAGEE.1.1))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    # colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = result),
               geom = 'text',
               node = 'all') + # optional for freerate
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median/median_freerate_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)


## graph free rate by branch length
# add free rate results to child branches
# just species nodes
as_tibble(tree.cagee) %>% 
  filter(is.na(result)) %>%
  left_join(as_tibble(tree.cagee) %>% 
              select(result, 
                     node) %>% 
              rename(parent = node) %>% 
              rename(freerate = result) %>% 
              na.omit()) %>% 
  select(branch.length,
         freerate) %>% 
  distinct() %>% 
  ggplot(aes(x = branch.length,
             y = as.numeric(freerate))) +
  geom_smooth(method = 'lm') +
  geom_point() +
  theme_classic() +
  ylab('Sigma freerate')
ggsave('CAGEE/median/median_freerate_cagee_outs/figures/Free rate sigma by branch length species.png')

# all nodes
as_tibble(tree.cagee) %>% 
  left_join(as_tibble(tree.cagee) %>% 
              select(result, 
                     node) %>% 
              rename(parent = node) %>% 
              rename(freerate = result) %>% 
              na.omit()) %>% 
  select(branch.length,
         freerate) %>% 
  distinct() %>% 
  ggplot(aes(x = branch.length,
             y = as.numeric(freerate))) +
  geom_smooth(method = 'lm') +
  geom_point() +
  theme_classic() +
  ylab('Sigma freerate')
ggsave('CAGEE/median/median_freerate_cagee_outs/figures/Free rate sigma by branch length.png')




#### CAGEE results median females ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/sex/median_cavity_1p_sex_F_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/sex/median_cavity_1p_sex_F_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

# ### Optional: add sigma tree to CAGEE tree
# ## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# ## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    # colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results female",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "F") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/sex/median_cavity_1p_sex_F_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)








#### CAGEE results median males ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/sex/median_cavity_1p_sex_M_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/sex/median_cavity_1p_sex_M_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
# tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

# ### Optional: add sigma tree to CAGEE tree
# ## create dummy tree
# tree.cagee.sigma = tree.cagee
# # add sigma
# tree.cagee.sigma$edge.length = tree_sigma$edge.length
# # convert to tibble
# # rename to sigma
# tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
#   dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

# ## Optional: add sigma
# tree.cagee = full_join(tree.cagee,
#                        tree.cagee.sigma) %>% 
#   mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    # colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results male",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "M") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/sex/median_cavity_1p_sex_M_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)








#### CAGEE results median females cavity ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/sex/median_cavity_sex_F_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/sex/median_cavity_sex_F_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>%
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>%
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results female",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "F") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/sex/median_cavity_sex_F_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)









#### CAGEE results median males cavity ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE/median_cavity/sex/median_cavity_sex_M_cagee_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE/median_cavity/sex/median_cavity_sex_M_cagee_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('CAGEE/data/10_species_birds_edit.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('CAGEE/data/10_species_birds_edit_cavity.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
# data.cagee.results = data.cagee.results %>% 
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end') %>% 
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
data.cagee.results = data.cagee.results %>%
  mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
                                     "Sigma2: ",
                                     "")) %>%
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end',
                       too_many = "merge") %>%
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>%
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>%
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results male",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "M") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE/median_cavity/sex/median_cavity_sex_M_cagee_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)









#### CAGEE results Automate ####
### load clade results
# 'clade_results.tab'
data.cagee.clade.results = read.delim('CAGEE_outs/clade_results.tab')

### load results
# results.txt
data.cagee.results = read.delim('CAGEE_outs/results.txt')

### load species tree
tree_sp = ape::read.tree('species_tree_used_for_CAGEE.nwk')

### optional: load sigma tree
tree_sigma = ape::read.tree('sigma_tree_used_for_CAGEE.nwk')

### clade results formatting
## separate data into columns 
# make species variable
# make tree label variable 
data.cagee.clade.results = data.cagee.clade.results %>% 
  separate_wider_delim(X.Taxon_ID,
                       delim = "<",
                       cols_remove = F,
                       names = c('species',"label")) %>% 
  mutate(label = paste0("<",
                        label))

### results formatting
## separate data into columns
# remove attempt row
data.cagee.results = data.cagee.results %>% 
  separate_wider_delim(CAGEE.1.1,
                       delim = ': ',
                       names = c("CAGEE.1.1",
                                 "result"),
                       too_few = 'align_end') %>% 
  mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                            "Attempts",
                            CAGEE.1.1))

## note: for multi-tissue analysis run this instead
## see special graphing instructions below too
# data.cagee.results = data.cagee.results %>%
#   mutate(CAGEE.1.1 = str_replace_all(CAGEE.1.1,
#                                      "Sigma2: ",
#                                      "")) %>%
#   separate_wider_delim(CAGEE.1.1,
#                        delim = ': ',
#                        names = c("CAGEE.1.1",
#                                  "result"),
#                        too_few = 'align_end',
#                        too_many = "merge") %>%
#   mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
#                             "Attempts",
#                             CAGEE.1.1))

### create CAGEE tree
# graph labels from CAGEE results
# add semicolon at end
tree.cagee =  ape::read.tree(text = paste(data.cagee.results[which(data.cagee.results$CAGEE.1.1=="IDs of Nodes"),2],
                                          ";",
                                          sep = ''))


### add edge length to CAGEE tree from species tree
tree.cagee$edge.length = tree_sp$edge.length

### Optional: add sigma tree to CAGEE tree
## create dummy tree
tree.cagee.sigma = tree.cagee
# add sigma
tree.cagee.sigma$edge.length = tree_sigma$edge.length
# convert to tibble
# rename to sigma
tree.cagee.sigma = as_tibble(tree.cagee.sigma) %>% 
  dplyr::rename(sigma = branch.length)

### add data results data to CAGEE tree
## convert tree to tibble
tree.cagee = as_tibble(tree.cagee)

## add data results data to CAGEE tree
tree.cagee = full_join(tree.cagee,
                       data.cagee.clade.results,
                       by = "label")

## Optional: add sigma
tree.cagee = full_join(tree.cagee,
                       tree.cagee.sigma) %>% 
  mutate(sigma = as.character(sigma))


## if you want to color specific branches you can make a variable here
## use label for any specific branch
## can use species.x to add metadata on those branches and species label tips
## then add as an aesthetic to ggtree: ggtree(aes(colour = new_variable))

# convert back to tree
tree.cagee = as.treedata(tree.cagee)

#### graph tree with results 
### create tree graph
p = tree.cagee %>% 
  ggtree(aes(
    colour = sigma #Optional
  )) +
  geom_nodelab(aes(label = paste0("+",Increase)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = 0.2,
               color = 'darkgreen') +
  geom_nodelab(aes(label = paste0("-",Decrease)),
               geom = 'text',
               node = 'all',
               nudge_x = -1,
               nudge_y = -0.2,
               color = 'darkred') +
  geom_tiplab(aes(label = species)) +
  theme_tree2() +
  ggtitle("CAGEE results",
          subtitle =  paste(data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Sigma2") %>% 
                              pull(result),
                            ";",
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(CAGEE.1.1),
                            data.cagee.results %>% 
                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>% 
                              pull(result))) + 
  scale_x_continuous(labels = abs) 

## reorder timeline
revts(p)

## save graph
ggsave('CAGEE_outs/figures/tree clade results.png',
       height = 10,
       width = 10)

# ## note for mult-itissue analysis need to add sigma values as table separately
# ## create sigma table
# sigma.table = data.cagee.results %>%
#   filter(grepl("Sigma",
#                CAGEE.1.1))
# 
# # graph sigma table
# library(ggpubr)
# p1 = ggpubr::ggtexttable(sigma.table,
#                  rows = NULL)
# 
# ## graph CAGEE results
# p = tree.cagee %>%
#   ggtree(aes(
#              colour = sigma #Optional
#              )) +
#   geom_nodelab(aes(label = paste0("+",Increase)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = 0.2,
#                color = 'darkgreen') +
#   geom_nodelab(aes(label = paste0("-",Decrease)),
#                geom = 'text',
#                node = 'all',
#                nudge_x = -1,
#                nudge_y = -0.2,
#                color = 'darkred') +
#   geom_tiplab(aes(label = species)) +
#   theme_tree2() +
#   ggtitle("CAGEE results",
#           subtitle = paste(data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(CAGEE.1.1),
#                            data.cagee.results %>%
#                              filter(CAGEE.1.1 == "Final Likelihood (-lnL)") %>%
#                              pull(result))) +
#   scale_x_continuous(labels = abs)
# 
# ## reorder timeline
# multiplot(revts(p),
#           p1,
#           widths = c(8,2),
#           ncol = 2)
# 
# ## save graph
# ggsave('CAGEE_outs/figures/tree clade results.png',
#        height = 10,
#        width = 10)







