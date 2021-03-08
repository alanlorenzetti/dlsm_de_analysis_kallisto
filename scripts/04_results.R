# alorenzetti 2021201

# description ####
# this script will download protein
# product annotation and will report the results
# it is going to download cog from halo_nr_tx
# repo as well and then perform enrichment analysis

# downloading functional annotation from remote resource ####
halo_nr_tx = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv")
nrtx = halo_nr_tx

# parsing to get products for all possible locus_tags
halo_nr_tx = halo_nr_tx %>% 
  dplyr::select(locus_tag,
                product) %>% 
  separate_rows(locus_tag,
                sep = ",")

# writing result files ####
finalRes = list()

# processing raw counts and TPM counts
finalRes[["kallistoEstCounts"]] = totrna %>% 
  mutate(target_id = str_replace(target_id, "\\|.*$", "")) %>% 
  left_join(., halo_nr_tx, by = c("target_id" = "locus_tag")) %>% 
  dplyr::rename(locus_tag = target_id)

finalRes[["kallistoTPM"]] = totrnaTPM %>% 
  mutate(target_id = str_replace(target_id, "\\|.*$", "")) %>% 
  left_join(., halo_nr_tx, by = c("target_id" = "locus_tag")) %>% 
  dplyr::rename(locus_tag = target_id)

# processing DESeqResults
for(i in names(results)){
  finalRes[["DESeq2"]][[i]] = results[[i]] %>% 
    mutate(target_id = str_replace(target_id, "\\|.*$", "")) %>% 
    left_join(., halo_nr_tx, by = c("target_id" = "locus_tag")) %>% 
    dplyr::rename(locus_tag = target_id)
}

# writing files
for(i in names(finalRes)){
  finalRes[[i]] %>% 
    write.xlsx(., file = paste0("results/", i, ".xlsx"))
}

# performing functional enrichment analysis ####
# genes will be considered differentially expressed
# if satisfying the pvalue threshold below
pthr = 0.01

# downloading COG data for halo_nr_tx
cog = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/cog.tsv")

# parsing deseq results to merge cog info
# djusting COG categories
# for unknown instances
inputobj = finalRes$DESeq2$dlsm_vs_dura3 %>%
  left_join(x = .,
            y = cog,
            by = c("locus_tag" = "representative")) %>% 
  mutate(cog_category = case_when(is.na(cog_category) ~ "Undefined",
                                  TRUE ~ as.character(cog_category)))

# adjusting COG functional categories
# there are a few genes having more than
# one functional class; for those I am
# going to keep the first one in order
# to reduce complexity
# I will also give "Mobilome..." class
# to all transposases and ISH elements
inputobj = inputobj %>% 
  mutate(cog_category = str_replace(cog_category, "\\|.*$", "")) %>% 
  mutate(cog_category = case_when(str_detect(string = product,
                                             pattern = "transposase|ISH") ~ "Mobilome: prophages, transposons",
                                  TRUE ~ as.character(cog_category)))

# adding diff express status to the obj
inputobj = inputobj %>% 
  filter(pvalue < pthr) %>% 
  mutate(status = case_when(log2FoldChange < 0 ~ "Downregulated",
                            log2FoldChange > 0 ~ "Upregulated",
                            TRUE ~ "NonDE") )

# preparing background object
# performing the same arbitrary
# modifications above to the 
# cog categories
funcatobj = nrtx %>% 
  left_join(x = .,
            y = cog,
            by = "representative") %>%
  mutate(cog_category = case_when(is.na(cog_category) ~ "Undefined",
                                  TRUE ~ as.character(cog_category))) %>% 
  mutate(cog_category = str_replace(cog_category, "\\|.*$", "")) %>% 
  mutate(cog_category = case_when(str_detect(string = product,
                                             pattern = "transposase|ISH") ~ "Mobilome: prophages, transposons",
                                  TRUE ~ as.character(cog_category)))

# testing for enrichment ####
# lsm interaction as primary clusters
# hypergeometric enrichment test of
# COG var inside lsm interaction status
vars = "cog_category"

# setting qval threshold
qthr = 0.01

# creating list to store results
enrich = list()
  
for(j in inputobj$status %>% unique()){
  curRegGroup = inputobj %>%
    filter(status == j)
  for(k in vars){
    curRegGroupVec = curRegGroup %>% 
      dplyr::select(all_of(k)) %>% 
      unlist(use.names = F)
    
    curRegGroupLvs = curRegGroupVec %>% 
      unique()
    
    for(l in curRegGroupLvs){
      wb = sum(curRegGroupVec == l)
      vecu = funcatobj %>% 
        dplyr::select(all_of(k)) %>% 
        unlist(use.names = F)
      wu = sum(vecu == l)
      bu = sum(vecu != l)
      drawn = curRegGroupVec %>% length()
      
      pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
      
      tib = tibble(regRule = j,
                   criteria = k,
                   level = l,
                   pval = pval)
      enrich = bind_rows(enrich, tib)
    }
  }
}

# correcting pvalues using BH method
# filtering by pval
enrich$qval = p.adjust(enrich$pval)
enrich = enrich[enrich$qval < qthr,]

# enrichement was found for only one class, but they were just
# three genes pertaining to the "Cell cycle ..." category
# I don't think it is worthy of mentioning the analysis