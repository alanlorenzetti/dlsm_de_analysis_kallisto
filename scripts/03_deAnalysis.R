# alorenzetti 2021201

# description ####
# this script will take data from kallisto
# perform differential expression
# analysis

# loading raw counts ####
# to manually compute TPMs, follow the instructions on the following
# page: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# reading totalrna counts
totrna=read_delim("data/tableEstCounts15.tsv",delim="\t")
totrnaTPM=read_delim("data/tableTpm15.tsv",delim="\t")

# wrangling colnames
# and removing genes with no counts
totrna = totrna %>%
  rename_with(.fn = ~ str_replace(.x, "-", "_TP"),
              .cols = starts_with("d")) %>% 
  filter(totrna %>%
           dplyr::select(starts_with("d")) %>%
           apply(X = ., MARGIN = 1, FUN = function(x) sum(x)) != 0)

# getting tpm results
totrnaTPM = totrnaTPM %>%
  rename_with(.fn = ~ str_replace(.x, "-", "_TP"),
              .cols = starts_with("d")) %>% 
  filter(totrnaTPM %>%
           dplyr::select(starts_with("d")) %>%
           apply(X = ., MARGIN = 1, FUN = function(x) sum(x)) != 0)
  

# correlation matrix shows
# that TP1 and TP2 are highly correlated
# I am going to use them as replicates
totrna %>%
  drop_na() %>% 
  dplyr::select(starts_with("d")) %>% 
  as.matrix() %>% 
  cor(method = "spearman") %>% 
  ggcorrplot(method = "circle")

# dropping TP3 and TP4
# totrna = totrna %>% 
#   dplyr::select(-matches("TP3|TP4"))

# processing raw counts with DESeq2  ####
## building se object
samples = totrna %>% 
  dplyr::select(starts_with("d")) %>% 
  colnames()
tp = totrna %>% 
  dplyr::select(starts_with("d")) %>% 
  colnames() %>% 
  str_replace(., "^.*_(.*)$", "\\1")
strain = totrna %>% 
  dplyr::select(starts_with("d")) %>% 
  colnames() %>% 
  str_replace(., "^(.*)_.*$", "\\1")

colData = data.frame(row.names = samples,
                     timepoint = tp,
                     strain = strain)

assay = dplyr::select(totrna, starts_with("d")) %>%
  as.matrix() %>% 
  round(digits = 0)

totrnaSE = SummarizedExperiment(assay = list(counts=assay),
                                rowData = totrna[,c(1:3)],
                                colData = colData)
rownames(totrnaSE) = rowData(totrnaSE)$target_id

# creating deseq2 objects
# design controls for the effect of timepoint
# https://stats.stackexchange.com/questions/17336/how-exactly-does-one-control-for-other-variables
totrnadds = totrnaSE
totrnadds = DESeqDataSet(totrnadds, design = ~ strain + timepoint)

# doing rlog transformation for distance and PCA
# before eliminating genes with zero counts
# save apart this object for exploratory analysis
rld = rlog(totrnadds)

# removing genes with zero counts and performing DESeq2 analysis
totrnadds = totrnadds[rowSums(counts(totrnadds)) > 1, ]
totrnadds = DESeq(totrnadds)

# result tables for contrast dlsm vs dura3
results = list()
results[["dlsm_vs_dura3"]] =  results(totrnadds,
                                      contrast = c("strain", "dlsm", "dura3"),
                                      alpha = padjthreshold) %>% 
  as_tibble(rownames = "target_id")

# final result tables
resultsSig = list()
resultsSig[["dlsm_vs_dura3"]] = results[["dlsm_vs_dura3"]] %>%
  dplyr::filter(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold) %>% 
  dplyr::filter(padj < padjthreshold)

# exploratory analysis ####
# setting distance matrix for dds
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(rld)
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)

# plotting heatmap 
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# # plotting principal componet analysis 
plotPCA(rld, intgroup = c("strain", "timepoint"))

# heatmap of gene clustering (based on rlog distance) 
# the 40 genes with highest variance between samples
topVarGenes = head(order(rowVars(assay(rld)),decreasing=TRUE), 40)
mat = assay(rld)[ topVarGenes, ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("strain", "timepoint")])
pheatmap(mat, annotation_col=df, width=10, height = 15)
