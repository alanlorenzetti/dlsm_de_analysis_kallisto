# alorenzetti 2021201

# description ####
# this script will download protein product annotation
# and will report the results

# downloading functional annotation from remote resource ####
halo_nr_tx = read_tsv(file = "https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv")

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