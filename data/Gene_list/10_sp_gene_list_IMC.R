#### Getting gene names from 10 sp bird data set
### convert NCBI ZebraFinchProteinID to ensembl gene name

## setwd
setwd("/geode2/home/u040/imillerc/Quartz/10_sp")

#### Get ZebraFinchProteinID list from data and save to 10sp_genes.txt
### upload 10sp_genes.txt into NCBI gene table: ncbi.nlm.nih.gov/datasets/gene
## Use 'Upload a text file' tab and in results use the 'select columns' to add:
# GeneID, Symbol, Gene name, Transcripts, Ensembl Gene ID, Proteins
## save all genes in table to: 10sp_NCBI_ensembl.tsv


#### load libraries ####
### install
## biomart
# BiocManager::install("biomaRt")

### load libraries 
library(biomaRt)
library(tidyverse)

#### load data ####
### Protein gene table
data = read.delim('./Gene_list/10sp_NCBI_ensembl.tsv')

#### biomart to get ensembl gene names ####
### selecting biomart database
ensembl = useMart('ensembl')
# list datasets
dataset = listDatasets(ensembl)
# select dataset
# zebrafinch
#Zebra finch genes (bTaeGut1_v1.p)
ensembl.zfinch = useDataset('tguttata_gene_ensembl',
                            mart = ensembl)
# get list of attributes
listAttributes(ensembl.zfinch) %>% 
  View()

# create zfinch attributes 
zfinch.attributes = c('external_gene_name',
                      'ensembl_gene_id')


# identify ensembl gene names
data.gene = getBM(attributes = zfinch.attributes,
                  mart = ensembl.zfinch,
                  values = data$Ensembl.GeneIDs,
                  filter = 'ensembl_gene_id',
                  useCache = FALSE) #useCache has to do with version of R?

### combine gene lists
data.gene.ID = data %>% 
  dplyr::select(Symbol,
                Ensembl.GeneIDs) %>% 
  left_join(data.gene,
            by = c('Ensembl.GeneIDs' = 'ensembl_gene_id'))

### filter down to genes with LOC name
data.gene.ID.filter = data.gene.ID %>% 
  filter(substr(Symbol, 1,3) == "LOC") %>% 
  filter(!is.na(external_gene_name)) %>% 
  filter(external_gene_name != "") 


### save results
write_tsv(data.gene.ID.filter,
          'Gene_list/10sp_NCBI_ensembl_LOC_names.tsv')


#### biomart to get gene chromosome ####
### selecting biomart database
ensembl = useMart('ensembl')
# list datasets
dataset = listDatasets(ensembl)
# select dataset
# zebrafinch
#Zebra finch genes (bTaeGut1_v1.p)
ensembl.zfinch = useDataset('tguttata_gene_ensembl',
                            mart = ensembl)
# get list of attributes
listAttributes(ensembl.zfinch) %>% 
  View()

# create zfinch attributes 
zfinch.attributes = c('external_gene_name',
                      'ensembl_gene_id',
                      'chromosome_name',
                      'start_position',
                      'end_position')


# identify ensembl gene names
data.chromosome = getBM(attributes = zfinch.attributes,
                  mart = ensembl.zfinch,
                  values = data$Ensembl.GeneIDs,
                  filter = 'ensembl_gene_id',
                  useCache = FALSE) #useCache has to do with version of R?

### combine gene lists
data.gene.chromosome = data %>% 
  left_join(data.chromosome %>% 
              distinct(),
            by = c('Ensembl.GeneIDs' = 'ensembl_gene_id'))


### save results
write_tsv(data.gene.chromosome,
          'Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')



#### biomart to get anolis chromosome for Z genes ####
#### load bird Z genes
### load gene chromosome position in birds
data.gene.chromosome = read_tsv('Gene_list/10sp_NCBI_ensembl_chromosome_position.tsv')

## get Z genes
data.gene.chromosome.z = data.gene.chromosome %>% 
  filter(chromosome_name == 'Z')

### selecting biomart database
ensembl.zfinch = useEnsembl(biomart='ensembl', 
                     dataset='tguttata_gene_ensembl', 
                     mirror = "useast")

# get list of attributes
listAttributes(ensembl.zfinch) %>% 
  View()

# create zfinch attributes 
zfinch.attributes = c('acarolinensis_homolog_ensembl_gene',
                      'acarolinensis_homolog_chromosome',
                      'acarolinensis_homolog_associated_gene_name',
                      'ensembl_gene_id')


# identify ensembl gene names
data.chromosome.anolis = getBM(attributes = zfinch.attributes,
                        mart = ensembl.zfinch,
                        values = data.gene.chromosome.z$Ensembl.GeneIDs,
                        filter = 'ensembl_gene_id',
                        useCache = FALSE) #useCache has to do with version of R?

### graph results
data.chromosome.anolis %>% 
  select(acarolinensis_homolog_chromosome) %>% 
  table() %>% 
  data.frame() %>% 
  mutate(Percent = 100*Freq/432) %>% 
  filter(Percent > 3) %>% 
  ggplot(aes(x = reorder(acarolinensis_homolog_chromosome,
                         -Percent),
             y= Percent,
             label = paste0(round(Percent),
                            "%"))) +
  geom_label() +
  theme_classic() +
  xlab('Anolis chromosome') +
  ylab('Percent Z genes') +
  ggtitle('Percent of Bird Z genes on A. carolinesis chromosomes')
ggsave('z_chromosome/Other_birds/figures/Anolis Z chromsome genes.png')

#### biomart to get gametologs from collared flycatcher to zebrafinch ####
### load in gametolog data
data.gametolog.fly = read.csv('Gene_list/Ficedula albicollis_gametologs.csv')


### selecting biomart database
ensembl = useMart('ensembl')
# list datasets
dataset = listDatasets(ensembl)
# select dataset
# zebrafinch
#Zebra finch genes (bTaeGut1_v1.p)
ensembl.fly = useDataset('falbicollis_gene_ensembl',
                            mart = ensembl)
# get list of attributes
listAttributes(ensembl.fly) %>% 
  View()

# create fly attributes 
fly.attributes = c('external_gene_name',
                      'ensembl_gene_id',
                   'tguttata_homolog_ensembl_gene',
                   'tguttata_homolog_associated_gene_name')


# identify ensembl gene names
data.gene.fly.zfinch = getBM(attributes = fly.attributes,
                  mart = ensembl.fly,
                  values = data.gametolog.fly$Gene.ID.of.Z.linked.gametolog_ensembl,
                  filter = 'ensembl_gene_id',
                  useCache = FALSE) #useCache has to do with version of R?

### get flycatcher gametologs without an ensembl match in zebra finch
# will go blast these genes individually to identify orthologs on NCBI
data.gametolog.fly.ens = data.gametolog.fly %>% 
  left_join(data.gene.fly.zfinch %>% 
              rename(Gene.ID.of.Z.linked.gametolog_ensembl = ensembl_gene_id))

## save list
write.csv(data.gametolog.fly.ens,
          'Gene_list/Ficedula albicollis_gametologs_zebrafinch.csv')

### combine gene lists
data.gametolog.zfinch = data %>% 
  select(c(Symbol,
           Description,
           Ensembl.GeneIDs)) %>% 
  distinct() %>% 
  left_join(data.gene.fly.zfinch %>% 
              select(tguttata_homolog_ensembl_gene) %>% 
              filter(tguttata_homolog_ensembl_gene != '') %>% 
              mutate(gametolog = 1) %>% 
              rename(Ensembl.GeneIDs = tguttata_homolog_ensembl_gene) %>% 
              distinct()) %>% 
  filter(gametolog == 1)


### save results
write_tsv(data.gametolog.zfinch,
          'Gene_list/10sp_ensembl_gametologs.tsv')
