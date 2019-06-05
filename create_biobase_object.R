library(Biobase)

expression <- read.csv(file.path(path, sprintf("%s_Kallisto_%s%s_%s", tolower(dataset), data_type, ifelse(feature_level=="gene", "_genes", ""), tolower(tissue))), stringsAsFactor=FALSE, row.names=1)
dim(expression)

eSet <- Biobase::ExpressionSet(as.matrix(expression), annotation=sprintf("%s_Kallisto_Gencode_v23_%s%s_%s", tolower(dataset), data_type, ifelse(feature_level=="gene", "_genes", "_isoforms"), tolower(tissue)))

sampleNames(eSet) <- gsub("[.]", "-", sampleNames(eSet))
pData(eSet) <- phenotypes[sampleNames(eSet),]

if(feature_level=="gene"){
  fData(eSet) <- toil.genes[featureNames(eSet),]
}else{
  fData(eSet) <- toil.transcripts[featureNames(eSet),]
}

save(eSet, file=file.path(path, "temp", sprintf("%s_Kallisto_Gencode_v23_%s%s_%s.RData", tolower(dataset), data_type, ifelse(feature_level=="gene", "_genes", "_isoforms"), tolower(tissue))))
