library("tidyr")
library("dplyr")

table <- read.table("data/processedData/miRNATotalTargets-parsed.txt", header = T, sep = "")
gtf <- read.table("genomes/anophelesStephensiBaseFeaturesAstel2.3.gtf", header = F, sep = "\t")

gtf <- separate(data = gtf, col = V9, sep = ";", into = c("V9", "V10", "V11"))

table$mirna <- as.character(table$mirna)
gtf$V9 <- as.character(gtf$V9)
gtf$V3 <- as.character(gtf$V3)

miRNA <- c()
ID <- c()
gene <- c()
CDS <- c()
exon <- c()
five_prime_utr <- c()
three_prime_utr <- c()
start_codon <- c()
stop_codon <- c()
transcript <- c()

for (n in 1:nrow(table)){
  for (i in 1:nrow(gtf)){
    if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
        (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "gene")){
      miRNA <- c(miRNA, table$mirna[n])
      ID <- c(ID, gtf$V9[i])
      gene <- c(gene, 1)
      CDS <- c(CDS, 0)
      exon <- c(exon, 0)
      five_prime_utr <- c(five_prime_utr, 0)
      three_prime_utr <- c(three_prime_utr, 0)
      start_codon <- c(start_codon, 0)
      stop_codon <- c(stop_codon, 0)
      transcript <- c(transcript, 0)
        }
      else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
               (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "CDS")){
                 miRNA <- c(miRNA, table$mirna[n])
                 ID <- c(ID, gtf$V9[i])
                 gene <- c(gene, 0)
                 CDS <- c(CDS, 1)
                 exon <- c(exon, 0)
                 five_prime_utr <- c(five_prime_utr, 0)
                 three_prime_utr <- c(three_prime_utr, 0)
                 start_codon <- c(start_codon, 0)
                 stop_codon <- c(stop_codon, 0)
                 transcript <- c(transcript, 0)    
               }
             else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                      (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "exon")){
                        miRNA <- c(miRNA, table$mirna[n])
                        ID <- c(ID, gtf$V9[i])
                        gene <- c(gene, 0)
                        CDS <- c(CDS, 0)
                        exon <- c(exon, 1)
                        five_prime_utr <- c(five_prime_utr, 0)
                        three_prime_utr <- c(three_prime_utr, 0)
                        start_codon <- c(start_codon, 0)
                        stop_codon <- c(stop_codon, 0)
                        transcript <- c(transcript, 0)    
                      }
                      else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                               (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "five_prime_utr")){
                                 miRNA <- c(miRNA, table$mirna[n])
                                 ID <- c(ID, gtf$V9[i])
                                 gene <- c(gene, 0)
                                 CDS <- c(CDS, 0)
                                 exon <- c(exon, 0)
                                 five_prime_utr <- c(five_prime_utr, 1)
                                 three_prime_utr <- c(three_prime_utr, 0)
                                 start_codon <- c(start_codon, 0)
                                 stop_codon <- c(stop_codon, 0)
                                 transcript <- c(transcript, 0)    
                               }
                               else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                                        (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "three_prime_utr")){
                                          miRNA <- c(miRNA, table$mirna[n])
                                          ID <- c(ID, gtf$V9[i])
                                          gene <- c(gene, 0)
                                          CDS <- c(CDS, 0)
                                          exon <- c(exon, 0)
                                          five_prime_utr <- c(five_prime_utr, 0)
                                          three_prime_utr <- c(three_prime_utr, 1)
                                          start_codon <- c(start_codon, 0)
                                          stop_codon <- c(stop_codon, 0)
                                          transcript <- c(transcript, 0)    
                                        }
                                        else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                                                 (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "start_codon")){
                                                   miRNA <- c(miRNA, table$mirna[n])
                                                   ID <- c(ID, gtf$V9[i])
                                                   gene <- c(gene, 0)
                                                   CDS <- c(CDS, 0)
                                                   exon <- c(exon, 0)
                                                   five_prime_utr <- c(five_prime_utr, 0)
                                                   three_prime_utr <- c(three_prime_utr, 0)
                                                   start_codon <- c(start_codon, 1)
                                                   stop_codon <- c(stop_codon, 0)
                                                   transcript <- c(transcript, 0)    
                                                 }
                                                 else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                                                          (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "stop_codon")){
                                                            miRNA <- c(miRNA, table$mirna[n])
                                                            ID <- c(ID, gtf$V9[i])
                                                            gene <- c(gene, 0)
                                                            CDS <- c(CDS, 0)
                                                            exon <- c(exon, 0)
                                                            five_prime_utr <- c(five_prime_utr, 0)
                                                            three_prime_utr <- c(three_prime_utr, 0)
                                                            start_codon <- c(start_codon, 0)
                                                            stop_codon <- c(stop_codon, 1)
                                                            transcript <- c(transcript, 0)    
                                                          }
                                                          else if ((table$Target[n] == gtf$V1[i]) && (table$Subject.Aln.start.[n] >= gtf$V4[i]) && 
                                                                   (table$Subject.Aln.end.[n] <= gtf$V5[i]) && (gtf$V3[i] == "transcript")){
                                                                     miRNA <- c(miRNA, table$mirna[n])
                                                                     ID <- c(ID, gtf$V9[i])
                                                                     gene <- c(gene, 0)
                                                                     CDS <- c(CDS, 0)
                                                                     exon <- c(exon, 0)
                                                                     five_prime_utr <- c(five_prime_utr, 0)
                                                                     three_prime_utr <- c(three_prime_utr, 0)
                                                                     start_codon <- c(start_codon, 0)
                                                                     stop_codon <- c(stop_codon, 0)
                                                                     transcript <- c(transcript, 1)    
                                                                   }
  }
}

miRNATargetAnnotations <- data.frame(miRNA, ID, gene, CDS, exon, five_prime_utr, three_prime_utr, start_codon, stop_codon, transcript)

write.table(miRNATargetAnnotations, file = "targetAnnotationsTest.txt")
