#### CAGEE results tree
### graph tree of CAGEE results with
## positive and negative changes listed on each branch
## Sigma2 and final likelihood
# R 4.2.1

### required inputs:
# 'clade_results.tab' from CAGEE
# 'results.txt' from CAGEE
# newick tree file used for CAGEE
# optional: newick sigma tree from CAGEE

## help making trees: yulab-smu.top/treedata-book/

## setwd
setwd("/directory_with_CAGEE_outs/")

#### load libraries ####
library(ape)
library(tidytree)
library(ggtree)
library(tidyverse)

#### load data ####
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

#### graph tree with results ####
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




