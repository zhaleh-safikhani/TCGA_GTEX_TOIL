options(stringsAsFactors=TRUE)
tissue <- "Breast"
data_type_index <- 1
data_set_index <- 1

data_type <- c("tpm", "counts")[data_type_index]
dataset <- c("GTEx", "TCGA")[data_set_index]
library(Biobase)
load("/mnt/work1/users/bhklab/Users/zhaleh/kallisto.hg38/Gencode.v23.annotation.RData", verbose=T)

transcript_level_file <- sprintf("/mnt/work1/users/bhklab/Data/%s/temp/%s_%s_%s.RData", dataset, tolower(dataset), tolower(tissue), data_type)
gene_level_file <- sprintf("/mnt/work1/users/bhklab/Data/%s/esets/%s_%s_%s.RData", dataset, tolower(dataset), tolower(tissue), data_type)

ff <- fData(transcript_level)
ee <- exprs(transcript_level)

epsilon <- c(0.001, 1)[data_type_index]
ee_gene_level <- t(sapply(rownames(toil.genes), function(gene){
  isoforms <- ee[which(ff$gene_id== gene),, drop=F]
log2(colSums(2 ^ isoforms - epsilon) + epsilon)}))

rownames(ee_gene_level) <- rownames(toil.genes)
colnames(ee_gene_level) <- colnames(ee)

gene_level <- ExpressionSet(ee_gene_level, annotation = "rnaseq")
fData(gene_level) <- toil.genes
pData(gene_level) <- pData(transcript_level)


##test #1
sum(2 ^ exprs(transcript_level)[which(fData(transcript_level)[,"gene_name"] == "THRA")] - epsilon)
2 ^ exprs(gene_level)[which(fData(gene_level)[,"gene_name"] == "THRA")] - epsilon

##test #2
sapply(1:3, function(x){sum(2 ^ exprs(transcript_level)[,x] - epsilon)})
sapply(1:3, function(x){sum(2 ^ exprs(gene_level)[,x] - epsilon)})

### A gene level eset contains both transcript and gene level
save(transcript_level, gene_level, file=gene_level_file)


