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
### load gene chromosome position
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


#### create random datasets for bootstrapping ####
### create 100 dataframes with randomly chosen male and female from each species 
### load gene chromosome position
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


### randomly sample a male and female from each species
for (i in 1:100) {
  # get list of random names
  sample_ids.random = sample_ids %>% 
    group_by(sex,
             species) %>% 
    sample_n(1) %>% 
    pull(sample_id)
  
  # subset data
  tmp = data.long %>% 
    filter(sample_id %in% sample_ids.random)
  
  # create long format
  tmp = tmp %>% 
    dplyr::select(-c(sample_name,
                     sample_id)) %>% 
    pivot_longer(cols = -c(sex,
                           species),
                 names_to = 'GeneName',
                 values_to = 'counts')
  
  # create wide format data
  # pivot wider for CAGEE
  tmp = tmp %>% 
    pivot_wider(names_from = 'species',
                values_from = 'counts')
  
  # add gene description
  tmp = full_join(data %>% 
                    dplyr::select(c(GeneDescription,
                                    GeneName)),
                  tmp)
  
  ## add z chromosome to other data frames
  tmp = tmp %>% 
    left_join(data.gene.chromosome.z) %>% 
    mutate(SAMPLETYPE = ifelse(is.na(SAMPLETYPE),
                               'A',
                               SAMPLETYPE),
           SAMPLETYPE = paste0(sex,
                               SAMPLETYPE)) %>% 
    dplyr::select(-c(sex)) %>% 
    relocate(SAMPLETYPE) %>% 
    relocate(GeneName) %>% 
    relocate(GeneDescription) 
  
  # save
  write_delim(tmp,
              paste0('CAGEE/data/random/median/normalizedCounts_median_sex_z_rand',
                     i,
                     '.tsv'),
              delim = '\t')
  
}

## create bash script to run CAGEE across all random files
# get list of files
files_in_folder <- list.files("CAGEE/data/random/median/")


## bootstrapping
## compare MZ vs FZ vs MA vs FA
## two tissue approach with MZ vs FZ vs MA vs FA
for (i in files_in_folder) {
  
  tmp = paste0("module load boost
  module load eigen
  module load intel
               ./CAGEE/CAGEE/build/cagee --tree ./CAGEE/data/10_species_birds_edit.nwk --infile ./CAGEE/data/random/median/", 
               i, 
               " --output_prefix CAGEE/median/median_cagee_outs_random/", 
               i,
               "_outs --sample_group MA --sample_group FA --sample_group MZ --sample_group FZ --cores 32")
  
  
  # tmp = paste0("./CAGEE/build/cagee --tree ./data/10_species_birds_edit.nwk --infile ./data/random/median/", 
  #        i, 
  #        " --output_prefix median/", 
  #        i,
  #        "_outs --sample_group MA --sample_group FA --sample_group MZ --sample_group FZ --cores 32")
  # 
  # cat(tmp)
  system(tmp)
}

#### Graph bootstrapping results ####   
## load model results
data.poster = read.csv('CAGEE/CAGEE_models_summary.csv')

### pull sigma and likelihood results from 100 bootstraps
# get list of files
files_in_folder <- list.files("CAGEE/data/random/median/")


## bootstrapping
## compare MZ vs FZ vs MA vs FA
## two tissue approach with MZ vs FZ vs MA vs FA
# create empty data frame
median_random_df = data.frame()

for (i in files_in_folder) {
  # get output results
  # get likelihood
  tmp = read.delim(paste0("./CAGEE/median/median_cagee_outs_random/", 
                         i,
                         "_outs/results.txt")) %>% 
    separate_wider_delim(CAGEE.1.1,
                         delim = 'Sigma for ',
                         names = c(NA,
                                   "CAGEE.1.1"),
                         too_few = 'align_end',
                         too_many = "merge") %>% 
    separate_wider_delim(CAGEE.1.1,
                         delim = ': ',
                         names = c("parameter",
                                   "sigma"),
                         too_few = 'align_end',
                         too_many = "merge") %>% 
    filter(parameter %in% c('MA',
                            'FA',
                            'MZ',
                            'FZ')) %>% 
    mutate(file = i,
           sigma = as.numeric(sigma)) %>% 
    full_join(read.delim(paste0("./CAGEE/median/median_cagee_outs_random/", # get likelihood
                          i,
                          "_outs/results.txt")) %>% 
    separate_wider_delim(CAGEE.1.1,
                         delim = 'Sigma for ',
                         names = c(NA,
                                   "CAGEE.1.1"),
                         too_few = 'align_end',
                         too_many = "merge") %>% 
    separate_wider_delim(CAGEE.1.1,
                         delim = ': ',
                         names = c("parameter",
                                   "likelihood"),
                         too_few = 'align_end',
                         too_many = "merge") %>% 
    filter(parameter %in% c('Final Likelihood (-lnL)')) %>% 
    mutate(file = i,
           likelihood = as.numeric(likelihood)) %>% 
    dplyr::select(-c(parameter)))
  
  # add to dataframe
  median_random_df = median_random_df %>% 
    rbind(tmp)
  
}

## get stats from sigma distribution
median_random_df_stat = median_random_df %>% 
  group_by(parameter) %>% 
  summarize(mean_val = mean(sigma),
            sd_val = sd(sigma),
            n_val = n(),
            sem = sd_val / sqrt(n_val),
            critical_t = qt(0.975, df = n_val - 1), # For 95% CI
            moe = critical_t * sem,
            lower_ci = mean_val - moe,
            upper_ci = mean_val + moe)

# add overall results
median_random_df_stat = median_random_df_stat %>% 
  mutate(type = 'bootstrap') %>% 
  full_join(
    data.poster %>%
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
      dplyr::select(type,
                    sigma) %>% 
      dplyr::rename(mean_val = sigma) %>% 
      dplyr::rename(parameter = type) %>% 
      mutate(type = 'overall')
  )

### graph results
## sigma distribution
median_random_df %>% 
  ggplot(aes(x = sigma)) +
  geom_histogram() +
  facet_wrap(.~parameter,
             nrow = 2,
             scales = 'free') +
  theme_classic()
ggsave('CAGEE/global_figures/bootstrap/sex_chrom_4parameter/Sigma distribution.png')


## likelihood distribution
median_random_df %>% 
  dplyr::select(file,
                likelihood) %>% 
  distinct() %>% 
  ggplot(aes(x = likelihood)) +
  geom_histogram()  +
  theme_classic()
ggsave('CAGEE/global_figures/bootstrap/sex_chrom_4parameter/Likelihood distribution.png')
  
## sigma 95% interval
median_random_df_stat %>% 
  mutate(parameter = factor(parameter,
                            levels = c("FA",
                                       "MA",
                                       "A",
                                       "FZ",
                                       "MZ",
                                       "Z"))) %>% 
  ggplot(aes(y = parameter)) +
  geom_errorbar(aes(xmin = lower_ci,
                    xmax = upper_ci),
                width = 0.1) +
  geom_point(aes(x = mean_val,
                 color = type))  +
  theme_classic() +
  xlab('sigma')
ggsave('CAGEE/global_figures/bootstrap/sex_chrom_4parameter/Compare sigma all.png')

# just bootstrap interval
median_random_df_stat %>% 
  filter(type == 'bootstrap') %>% 
  mutate(parameter = factor(parameter,
                levels = c("FA",
                           "MA", 
                           "FZ",
                           "MZ"))) %>% 
  ggplot(aes(y = parameter)) +
  geom_errorbar(aes(xmin = lower_ci,
                   xmax = upper_ci),
                width = 0.1) +
  geom_point(aes(x = mean_val))  +
  theme_classic() +
  xlab('sigma')
ggsave('CAGEE/global_figures/bootstrap/sex_chrom_4parameter/Compare sigma bootstrap.png')


  













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



#### load sim n_gamma_cats data low ####
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


### combine into one data frame
ngamma_low.sim.df = 
  rbind(ngamma_low_median_4_sim) %>% 
  rbind(ngamma_low_ratio_4_sim) %>% 
  rbind(ngamma_low_ratio_4_sim_large) %>% 
  rbind(ngamma_low_median_4_sim_large)

## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma_low.sim.df.table = ngamma_low.sim.df %>% 
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

#### graph sim n_gamma_cats data low ####
### graph genes per category
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
# paper
ngamma_low.sim.df.table %>% 
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

#### load n_gamma_cats data low ####
### assign to lowest liklihood

# remove NA
# pivot long
# only keep category with the top liklihood score per gene

#### load results from gamma category testing
### both ratio and median
# make empty dataframe
ngamma_low.df = data.frame()

### loop through multiple gamma categories
for (i in c(2,4,6,8,10)) {
  
  # combine data
  ngamma_low.df = ngamma_low.df %>% 
    rbind(read_tsv(paste0('CAGEE/median/median_cagee_outs_n_gamma_cats_',i,'/category_likelihoods.txt')) %>% #median
            dplyr::select(-last_col()) %>% 
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
                   gamma.count = i,
                   sigma.all = read.delim(paste0('CAGEE/median/median_cagee_outs_n_gamma_cats_',i,'/results.txt')) %>% # add sigma value
                     separate_wider_delim(CAGEE.1.1,
                                          delim = ': ',
                                          names = c("CAGEE.1.1",
                                                    "result"),
                                          too_few = 'align_end',
                                          too_many = "merge") %>%
                     mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                                               "Attempts",
                                               CAGEE.1.1)) %>% 
                     filter(CAGEE.1.1 == "Sigma2") %>% 
                     pull(result) %>% 
                     as.numeric()) %>% 
            group_by(transcript) %>% 
            slice(1) %>% 
            ungroup()) %>% 
    rbind(read_tsv(paste0('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_',i,'/category_likelihoods.txt')) %>% #ratio 
            dplyr::select(-last_col()) %>% 
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
                   gamma.count = i,
                   sigma.all = read.delim(paste0('CAGEE/ratio/ratio_cagee_outs_n_gamma_cats_',i,'/results.txt')) %>% # add sigma value
                     separate_wider_delim(CAGEE.1.1,
                                          delim = ': ',
                                          names = c("CAGEE.1.1",
                                                    "result"),
                                          too_few = 'align_end',
                                          too_many = "merge") %>%
                     mutate(CAGEE.1.1 = ifelse(is.na(CAGEE.1.1),
                                               "Attempts",
                                               CAGEE.1.1)) %>% 
                     filter(CAGEE.1.1 == "Sigma2") %>% 
                     pull(result) %>% 
                     as.numeric())%>% 
            group_by(transcript) %>% 
            slice(1) %>% 
            ungroup()) 
  
}

## calculate sigma
ngamma_low.df = ngamma_low.df %>% 
  mutate(sigma = sigma.all * gamma.cat) %>% 
  dplyr::select(-c(liklihood,
                   top.cat,
                   keep))

# add k = 1 for median and sigma
ngamma_low.df = ngamma_low.df %>% 
  rbind(data.frame(transcript = unique(ngamma_low.df$transcript),
             gamma.cat = 1,
             type = 'median',
             gene.list = 'all',
             gamma.count = 1,
             sigma.all = 0.011384,
             sigma = 0.011384)) %>% 
  rbind(data.frame(transcript = unique(ngamma_low.df$transcript),
                   gamma.cat = 1,
                   type = 'ratio',
                   gene.list = 'all',
                   gamma.count = 1,
                   sigma.all = 0.000381,
                   sigma = 0.000381))

## get count of genes per gamma category per type, gene.list, and gamma.count
# calculate percentage
ngamma_low.df.table = ngamma_low.df %>% 
  dplyr::select(c(sigma,
                  type,
                  gene.list,
                  gamma.count)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  group_by(type,
           gene.list,
           gamma.count) %>% 
  mutate(total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(percent = 100*Freq/total)

#### graph n_gamma_cats data low ####
### graph sigma per category
ngamma_low.df %>%
  dplyr::select(gamma.count,
                sigma.all,
                type) %>% 
  distinct() %>% 
  ggplot(aes(x = gamma.count,
             y = sigma.all)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  xlab('k categories') +  
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5,
                                                    nint = TRUE)) +
  facet_wrap(. ~ type,
             scales = "free_y") 
ggsave('CAGEE/global_figures/n_gamma/compare/Sigma across k categories.png',
       width = 6,
       height = 3)

### graph genes per category
ngamma_low.df.table %>%
  ggplot(aes(x = as.numeric(as.character(sigma)),
             y = percent,
             group = gamma.count,
             color = gamma.count)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  facet_grid(. ~ type,
             scales = 'free_x'
             ) +
  xlab('sigma') +
  labs(color = 'k = ')
ggsave('CAGEE/global_figures/n_gamma/compare/Percent of genes per sigma category low.png')


# ### compare number of gamma categories 
# ## heat map of k = 4 vs k = 10
# ## real
# # median
# ngamma_low.df %>% 
#   filter(type == 'median') %>% 
#   filter(gamma.count == '4' | gamma.count == '10') %>%
#   mutate(sigma = round(sigma, 5)) %>% 
#   pivot_wider(id_cols = 'transcript',
#               names_from = 'gamma.count',
#               names_prefix = 'k.',
#               values_from = sigma) %>% 
#   column_to_rownames('transcript') %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0 ) %>% 
#   ggplot(aes(x = as.numeric(as.character(k.10)),
#              y = as.numeric(as.character(k.4)),
#              fill = Freq,
#              label = Freq)) +
#   geom_tile() +
#   geom_text() +
#   theme_classic() +
#   scale_fill_gradient(low = 'white',
#                       high = 'darkred') +
#   xlab('k = 10') +
#   ylab('k = 4') +
#   ggtitle('median') +
#   coord_fixed()
# ggsave('CAGEE/global_figures/n_gamma/compare/median k=4 vs k=10.png')
# 
# # ratio
# ngamma_low.df %>% 
#   filter(type == 'ratio') %>%
#   filter(gamma.count == '4' | gamma.count == '10') %>%
#   mutate(sigma = round(sigma, 5)) %>% 
#   pivot_wider(id_cols = 'transcript',
#               names_from = 'gamma.count',
#               names_prefix = 'k.',
#               values_from = sigma) %>% 
#   column_to_rownames('transcript') %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0 ) %>% 
#   ggplot(aes(x = as.numeric(as.character(k.10)),
#              y = as.numeric(as.character(k.4)),
#              fill = Freq,
#              label = Freq)) +
#   geom_tile() +
#   geom_text() +
#   theme_classic() +
#   scale_fill_gradient(low = 'white',
#                       high = 'darkred') +
#   xlab('k = 10') +
#   ylab('k = 4')+
#   ggtitle('ratio')+
#   coord_fixed()
# ggsave('CAGEE/global_figures/n_gamma/compare/ratio k=4 vs k=10.png')
# 
# ## relative
# # median
# ngamma_low.df %>% 
#   filter(type == 'median') %>% 
#   filter(gamma.count == '4' | gamma.count == '10') %>%
#   mutate(sigma = round(sigma, 5)) %>% 
#   pivot_wider(id_cols = 'transcript',
#               names_from = 'gamma.count',
#               names_prefix = 'k.',
#               values_from = sigma) %>% 
#   column_to_rownames('transcript') %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0 ) %>% 
#   ggplot(aes(x = k.10,
#              y = k.4,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile() +
#   geom_text() +
#   theme_classic() +
#   scale_fill_gradient(low = 'white',
#                       high = 'darkred') +
#   xlab('k = 10') +
#   ylab('k = 4') +
#   ggtitle('median')
# ggsave('CAGEE/global_figures/n_gamma/compare/median k=4 vs k=10 relative.png')
# 
# # ratio
# ngamma_low.df %>% 
#   filter(type == 'ratio') %>%
#   filter(gamma.count == '4' | gamma.count == '10') %>%
#   mutate(sigma = round(sigma, 5)) %>% 
#   pivot_wider(id_cols = 'transcript',
#               names_from = 'gamma.count',
#               names_prefix = 'k.',
#               values_from = sigma) %>% 
#   column_to_rownames('transcript') %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0 ) %>% 
#   ggplot(aes(x = k.10,
#              y = k.4,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile() +
#   geom_text() +
#   theme_classic() +
#   scale_fill_gradient(low = 'white',
#                       high = 'darkred') +
#   xlab('k = 10') +
#   ylab('k = 4')+
#   ggtitle('ratio')
# ggsave('CAGEE/global_figures/n_gamma/compare/ratio k=4 vs k=10 relative.png')

### create  plot of genes across K gamma categories
library(igraph)
library(ggraph)

## median
## create edge list
# create empty dataframe
ngamma_low.df.median.edge = data.frame()
# loop through K values
for (i in c(2,4,6,8)) {

tmp = ngamma_low.df %>% 
  filter(type == 'median') %>%
  dplyr::select(gamma.count,
                sigma,
                transcript) %>% 
  filter(gamma.count %in% c(i,
                            i+2)) %>% 
  mutate(sigma = round(sigma, 5),
         gamma.count = paste0('k_',gamma.count)) %>% 
  pivot_wider(names_from = 'gamma.count',
              values_from = 'sigma') %>% 
  column_to_rownames('transcript') %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0)

# rename columns
colnames(tmp)[1:2] = c('X',
                       'Y')

# add gamma count to sigma
tmp = tmp %>% 
  mutate(X = paste0(i, "_", X),
         Y = paste0(i+2, "_", Y))


# add percentage 
tmp = tmp %>% 
  group_by(X) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total)


# add to dataframe
ngamma_low.df.median.edge = ngamma_low.df.median.edge %>% 
  rbind(tmp)

}

# create nodes list
# add counts per category
ngamma_low.df.median.node = data.frame(name = unique(c(ngamma_low.df.median.edge$X,
                                                       ngamma_low.df.median.edge$Y))) %>% 
  separate_wider_delim(cols = name,
                       delim = '_',
                       names = c('k',
                                 'sigma'),
                      cols_remove = F) %>% 
  dplyr::relocate(name) %>% 
  mutate(k = as.numeric(k),
         sigma = as.numeric(sigma)) %>% 
  left_join(ngamma_low.df.table %>% 
  filter(type == 'median') %>% 
  mutate(name = paste0(gamma.count, 
                       '_',
                       round(as.numeric(as.character(sigma)),
                                               5)),
         percent = round(percent, 0)) %>% 
  dplyr::select(name,
                percent))

# make graph object
ngamma_low.df.median.net = igraph::graph_from_data_frame(d = ngamma_low.df.median.edge,
                               vertices = ngamma_low.df.median.node,
                               directed = F)
  
  
# graph object
# # freq
# ggraph(ngamma_low.df.median.net,
#        layout = 'manual',
#        y = k,
#        x = sigma) +
#   geom_edge_link(aes(edge_width = Freq,
#                      edge_colour = Freq))+
#   geom_node_label(aes(label = paste0(percent,'%'))) +
#   theme_classic() +
#   scale_y_reverse() +
#   scale_edge_color_gradient(low = 'white',
#                        high = 'darkred') +
#   xlab('sigma') +
#   ylab('k categories') +
#   ggtitle('Median')
# ggsave('CAGEE/global_figures/n_gamma/compare/median k category flow.png')

# percent
ggraph(ngamma_low.df.median.net,
       layout = 'manual',
       y = k,
       x = sigma) +
  geom_edge_link(aes(edge_width = Percent,
                     edge_colour = Percent))+
  geom_node_label(aes(label = paste0(percent,'%')),
                  size = 2) +
  theme_classic() +
  scale_y_reverse() +
  scale_edge_color_gradient(low = 'white',
                            high = 'darkred') +
  xlab('sigma') +
  ylab('k categories') +
  ggtitle('Median')
ggsave('CAGEE/global_figures/n_gamma/compare/median k category flow percent.png',
       width = 6,
       height = 3)



## ratio
## create edge list
# create empty dataframe
ngamma_low.df.ratio.edge = data.frame()
# loop through K values
for (i in c(2,4,6,8)) {
  
  tmp = ngamma_low.df %>% 
    filter(type == 'ratio') %>%
    dplyr::select(gamma.count,
                  sigma,
                  transcript) %>% 
    filter(gamma.count %in% c(i,
                              i+2)) %>% 
    mutate(sigma = round(sigma, 5),
           gamma.count = paste0('k_',gamma.count)) %>% 
    pivot_wider(names_from = 'gamma.count',
                values_from = 'sigma') %>% 
    column_to_rownames('transcript') %>% 
    table() %>% 
    as.data.frame() %>% 
    filter(Freq > 0)
  
  # rename columns
  colnames(tmp)[1:2] = c('X',
                         'Y')
  
  # add gamma count to sigma
  tmp = tmp %>% 
    mutate(X = paste0(i, "_", X),
           Y = paste0(i+2, "_", Y))
  
  # add percentage 
  tmp = tmp %>% 
    group_by(X) %>% 
    mutate(Total = sum(Freq)) %>% 
    ungroup() %>% 
    mutate(Percent = 100*Freq/Total)
  
  # add to dataframe
  ngamma_low.df.ratio.edge = ngamma_low.df.ratio.edge %>% 
    rbind(tmp)
  
}

# create nodes list
# add counts per category
ngamma_low.df.ratio.node = data.frame(name = unique(c(ngamma_low.df.ratio.edge$X,
                                                       ngamma_low.df.ratio.edge$Y))) %>% 
  separate_wider_delim(cols = name,
                       delim = '_',
                       names = c('k',
                                 'sigma'),
                       cols_remove = F) %>% 
  dplyr::relocate(name) %>% 
  mutate(k = as.numeric(k),
         sigma = as.numeric(sigma)) %>% 
  left_join(ngamma_low.df.table %>% 
              filter(type == 'ratio') %>% 
              mutate(name = paste0(gamma.count, 
                                   '_',
                                   round(as.numeric(as.character(sigma)),
                                         5)),
                     percent = round(percent, 0)) %>% 
              dplyr::select(name,
                            percent))

# make graph object
ngamma_low.df.ratio.net = igraph::graph_from_data_frame(d = ngamma_low.df.ratio.edge,
                                                         vertices = ngamma_low.df.ratio.node,
                                                         directed = F)


# graph object
# ggraph(ngamma_low.df.ratio.net,
#        layout = 'manual',
#        y = k,
#        x = sigma) +
#   geom_edge_link(aes(edge_width = Freq,
#                      edge_colour = Freq))+
#   geom_node_label(aes(label = paste0(percent,'%'))) +
#   theme_classic() +
#   scale_y_reverse() +
#   scale_edge_color_gradient(low = 'white',
#                             high = 'darkred') +
#   xlab('sigma') +
#   ylab('k categories') +
#   ggtitle('ratio')
# ggsave('CAGEE/global_figures/n_gamma/compare/ratio k category flow.png')

# percent
ggraph(ngamma_low.df.ratio.net,
       layout = 'manual',
       y = k,
       x = sigma) +
  geom_edge_link(aes(edge_width = Percent,
                     edge_colour = Percent))+
  geom_node_label(aes(label = paste0(percent,'%')),
                  size = 2) +
  theme_classic() +
  scale_y_reverse() +
  scale_edge_color_gradient(low = 'white',
                            high = 'darkred') +
  xlab('sigma') +
  ylab('k categories') +
  ggtitle('ratio')
ggsave('CAGEE/global_figures/n_gamma/compare/ratio k category flow percent.png',
       width = 6,
       height = 3)







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







