#module load gcc/6.2.0 
#module load R/3.4.0
#R

options(stringsAsFactors=TRUE)
tissue <- "Pancreas"
data_type_index <- 1
data_set_index <- 1
feature_level_index <- 1

data_type <- c("tpm", "counts")[data_type_index]
dataset <- c("GTEx", "TCGA")[data_set_index]
feature_level <- c("gene", "transcript")[feature_level_index]

load("/mnt/work1/users/bhklab/Users/zhaleh/kallisto.hg38/Gencode.v23.annotation.RData", verbose=T)

if(dataset == "GTEx"){
  path <- "/mnt/work1/users/bhklab/Data/GTEx" #"~/Google Drive/Toil_tcga_gtex/gtext_annotation"
  phenotypes <- read.csv(file.path(path, "GTEX_phenotype"), header=T, sep="\t", row.names=1)
}else{
  path <- "/mnt/work1/users/bhklab/Data/TCGA" #"~/Google Drive/Toil_tcga_gtex/tcga_annotation"
  phenotypes <- read.csv(file.path(path, "PANCAN_clinicalMatrix"), header=T, sep="\t", row.names=1)
}
dim(phenotypes)
phenotypes[1:3, 1:3]

ff <- file.path(path, "sample_names")
con <- file(ff, open = "r")
sample.names <- readLines(con, n = 1, warn = FALSE)
close(con)
sample.names <- unlist(strsplit(sample.names, split="\t"))
sample.names <- sample.names[-1]
length(sample.names)
sample.names[1:3]
sample.names[which(duplicated(sample.names))]
length(which(sample.names %in% rownames(phenotypes)))

tissue.specific <- rownames(phenotypes)[which(phenotypes$X_primary_site == tissue)]
xx <- which(sample.names %in% tissue.specific)

ff <- file.path(path, sprintf("%s_Kallisto_%s%s", tolower(dataset), data_type, ifelse(feature_level=="gene", "_genes", "")))
fo <- file.path(path, sprintf("%s_Kallisto_%s%s_%s", tolower(dataset), data_type, ifelse(feature_level=="gene", "_genes", ""), tolower(tissue)))
con <- file(ff, open = "r")
flag <- T
i <- 1
rr <- readLines(con, n = 1, warn = FALSE)
if(feature_level=="gene"){rr <- unlist(strsplit(rr, split=","))}else{rr <- unlist(strsplit(rr, split="\t"))}
if(feature_level=="gene"){yy <- xx}else{yy <- which(rr %in% tissue.specific)}
if(feature_level=="gene"){
  cat(paste(c("sample", rr[yy]), collapse=","), sep="\n", file=fo)
}else{
    cat(paste(c("sample", sample.names[yy]), collapse=","), sep="\n", file=fo)
  }
while(flag){
  rr <- readLines(con, n = 1, warn = FALSE)
  if(length(rr) != 0){
  #if(i != 3){
    if(feature_level=="gene"){
      rr <- unlist(strsplit(rr, split=","))
      cat(paste(c(toil.genes[i, "gene_id"], rr[xx]), collapse=","), sep="\n", file=fo, append=TRUE)
    }else{
      rr <- unlist(strsplit(rr, split="\t"))
      cat(paste(c(rr[1], rr[yy]), collapse=","), sep="\n", file=fo, append=TRUE)
    }
  }else{
    flag <- FALSE
  }
  print(i)
  i <- i + 1
}
close(con)
