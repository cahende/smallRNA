library(seqinr)

novelMiRNA <- read.fasta(file = snakemake@input[[3]], seqtype = "DNA", as.string = TRUE)
novel <- read.table(file = snakemake@input[[1]], sep = "\t", header = T)
known <- read.table(file = snakemake@input[[2]], sep = "\t", header = T)

novelMiRNA <- data.frame(unlist(novelMiRNA))
novelMiRNA$seq <- novelMiRNA$unlist.novelMiRNA.
novelMiRNA$ID <- row.names(novelMiRNA)
novelMiRNA <- data.frame(novelMiRNA$ID, novelMiRNA$seq)

novelMiRNA$novelMiRNA.seq <- as.character(novelMiRNA$novelMiRNA.seq)
colnames(novelMiRNA) <- c("new.id", "consensus.mature.sequence")

novel$consensus.mature.sequence <- as.character(novel$consensus.mature.sequence)

novelKeep <- subset(novel, novel$miRDeep2.score >= 3 & novel$significant.randfold.p.value == "yes")
novelMorph <- merge(novelKeep, novelMiRNA, by = "consensus.mature.sequence")

novelNew <- data.frame(novelMorph$new.id, novelMorph$consensus.mature.sequence, novelMorph$consensus.precursor.sequence, novelMorph$mature.read.count, novelMorph$total.read.count, novelMorph$loop.read.count, novelMorph$star.read.count)
colnames(novelNew) <- c("id", "consensus.mauture.sequence", "consensus.precursor.sequence", "mature.read.count", "total.read.count", "loop.read.count", "star.read.count")

known$consensus.mature.sequence <- as.character(known$consensus.mature.sequence)

knownKeep <- subset(known, known$miRDeep2.score >= 3 & known$significant.randfold.p.value == "yes")

knownNew <- data.frame(knownKeep$mature.miRBase.miRNA, knownKeep$consensus.mature.sequence, knownKeep$consensus.precursor.sequence, knownKeep$mature.read.count, knownKeep$total.read.count, knownKeep$loop.read.count, knownKeep$star.read.count)
colnames(knownNew) <- c("id", "consensus.mauture.sequence", "consensus.precursor.sequence", "mature.read.count", "total.read.count", "loop.read.count", "star.read.count")

knownUnique <- unique(knownNew)
novelUnique <- unique(novelNew)

write.csv(novelUnique, file = snakemake@output[[1]], row.names = F)
write.csv(knownUnique, file = snakemake@output[[2]], row.names = F)
