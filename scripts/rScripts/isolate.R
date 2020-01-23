# -------------------------------------------------------------------------------------- #
#miRNA - Isolate novel and known miRNAs from mirDeep2 output
# -------------------------------------------------------------------------------------- #


# load packages and set empty vectors
library(seqinr)

IDNovel <- c()
sequencesNovel <- c()
importantNovel <- c()
matureCountsNovel <- c()
miRBaseNameNovel <- c()

sequencesKnown <- c()
IDKnown <- c()
importantKnown <- c()
matureCountsKnown <- c()
miRBaseNameKnown <- c()

IDTotal <- c()
sequencesTotal <- c()
matureCountsTotal <- c()
miRBaseNameTotal <- c()

# set names and fix data type prior to analysis
novel <- read.csv(snakemake@input[[1]], header = T, sep = "\t")
novel$provisional.id <- as.character(novel$provisional.id)
novel$consensus.mature.sequence <- as.character(novel$consensus.mature.sequence)

known <- read.csv(snakemake@input[[2]], header = T, sep = "\t")
known$mature.miRBase.miRNA <- as.character(known$mature.miRBase.miRNA)
known$consensus.mature.sequence <- as.character(known$consensus.mature.sequence)
known$tag.id <- as.character(known$tag.id)

# isolate novel miRNAs with score >3 and significant randfold value
for (i in 1:nrow(novel)){
  if ((novel$miRDeep2.score[i] >= 3) & (novel$significant.randfold.p.value[i] == "yes")){
    (importantNovel <- c(importantNovel, i))
  } 
}
 
# isoate novel usable miRNA IDs, sequences, and read counts 
for (n in 1:length(importantNovel)){
  IDNovel <- c(IDNovel, novel$provisional.id[importantNovel[n]])
}

for (n in 1:length(importantNovel)){
  sequencesNovel <- c(sequencesNovel, novel$consensus.mature.sequence[importantNovel[n]])
} 

for (n in 1:length(importantNovel)){
  matureCountsNovel <- c(matureCountsNovel, novel$mature.read.count[importantNovel[n]])
} 

for (n in 1:length(importantNovel)){
	miRBaseNameNovel <- c(miRBaseNameNovel, novel$miRBase.miRNA[importantNovel[n]])
}

# isolate known miRNAs with score >3 and significant randfold value
for (i in 1:nrow(known)){
  if ((known$miRDeep2.score[i] >= 3) & (known$significant.randfold.p.value[i] == "yes")){
    (importantKnown <- c(importantKnown, i))
  } 
}

# isoate usable miRNA IDs, sequences, and read counts 
for (n in 1:length(importantKnown)){
  IDKnown <- c(IDKnown, known$tag.id[importantKnown[n]])
}

for (n in 1:length(importantKnown)){
  sequencesKnown <- c(sequencesKnown, known$consensus.mature.sequence[importantKnown[n]])
} 

for (n in 1:length(importantKnown)){
  matureCountsKnown <- c(matureCountsKnown, known$mature.read.count[importantKnown[n]])
} 

for (n in 1:length(importantKnown)){
	miRBaseNameKnown <- c(miRBaseNameKnown, known$mature.miRBase.miRNA[importantKnown[n]])
}

# combine known and novel read counts into single table
IDTotal <- c(IDNovel, IDKnown)
sequencesTotal <- c(sequencesNovel, sequencesKnown)
matureCountsTotal <- c(matureCountsNovel, matureCountsKnown)
miRBaseNameTotal <- c(miRBaseNameNovel, miRBaseNameKnown)
tableTotal <- data.frame(IDTotal, sequencesTotal, matureCountsTotal, miRBaseNameTotal)
write.csv(tableTotal, file=snakemake@output[[1]], row.names = F)

# write .fa for novel sequences
for (i in 1:length(importantNovel)){ 
  write.fasta(sequencesNovel[i], IDNovel[i], file.out = snakemake@output[[2]], open = "a")
}

# write .fa for known sequences
for (i in 1:length(importantKnown)){ 
  write.fasta(sequencesKnown[i], IDKnown[i], file.out = snakemake@output[[3]], open = "a")
}
