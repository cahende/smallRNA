#!/usr/bin/env/python
shell.prefix("set -o pipefail; ")
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
configfile: "config.yaml"

rule all:
    input:
        expand("data/processedData/virus/{sample}/star_output/{sample}-piRNA-virus.bam", sample=config["SAMPLES"]),
        expand("data/processedData/piRNA/{sample}/star_output/{sample}-piRNA.bam", sample=config["SAMPLES"]),
        expand("data/processedData/totalMap/{sample}/star_output/{sample}.bam", sample=config["SAMPLES"]),
        expand("data/processedData/virusTotal/{sample}/star_output/{sample}-virus.bam", sample=config["SAMPLES"]),
        expand("data/processedData/condensedOutputVirus/{sampleVirus}-novelVirus.csv", sampleVirus=config["SAMPLES_VIRUS"])
rule trim_reads:
    input:
        "data/rawData/{sample}_R1_001.fastq.gz"
    output:
        "data/processedData/trimmed_reads/{sample}_trimmed.fastq.gz"
    log: "logs/trim/{sample}.trim.log"
    shell:
        "conda activate bioinfo;"
        "fastqc {input};"
        "java -jar {config[TRIMMOMATIC]} SE -phred33 -trimlog {log} {input} {output} ILLUMINACLIP:{config[ADAPTERS]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15"

rule decompress_trimmed_reads:
    input:
        "data/processedData/trimmed_reads/{sample}_trimmed.fastq.gz"
    output:
        "data/processedData/trimmed_reads/{sample}_trimmed.fastq"
    shell:
        "gunzip -dc {input} > {output}"

rule convert_to_fasta:
    input:
        "data/processedData/trimmed_reads/{sample}_trimmed.fastq"
    output:
        "data/processedData/mirDeep/{sample}/{sample}.fasta"
    shell:
        "conda activate bioinfo;"
        "sed '/^@/!d;s//>/;N' {input} > {output};"
        "sed 's, ,_,g' -i {output}"
        
rule map:
    input:
        "data/processedData/mirDeep/{sample}/{sample}.fasta"
    output:
        "data/processedData/mirDeep/{sample}/{sample}.arf",
        "data/processedData/mirDeep/{sample}/{sample}-collapsed.fa"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeep/{wildcards.sample}/;"
        "mapper.pl {config[WORKDIR]}/{input} -c -v -i -m -p {config[WORKDIR]}/{config[GENOME]} -s {config[WORKDIR]}/{output[1]} -t {config[WORKDIR]}/{output[0]}"

rule mapVirus:
    input:
        "data/processedData/mirDeep/{sampleVirus}/{sampleVirus}.fasta"
    output:
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}.arf",
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-collapsed.fa"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeep/{wildcards.sampleVirus}/;"
        "mapper.pl {config[WORKDIR]}/{input} -c -v -i -m -p {config[WORKDIR]}/{config[VIRUS_GENOME]} -s {config[WORKDIR]}/{output[1]} -t {config[WORKDIR]}/{output[0]}"

rule mirDeep2:
    input:
        "data/processedData/mirDeep/{sample}/{sample}.arf",
        "data/processedData/mirDeep/{sample}/{sample}-collapsed.fa"
    output:
        "data/processedData/mirDeep/{sample}/{sample}-report.log",
        "data/processedData/mirDeep/{sample}/{sample}-results.csv"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeep/{wildcards.sample}/;"
        "miRDeep2.pl {config[WORKDIR]}/{input[1]} {config[WORKDIR]}/{config[GENOME]} {config[WORKDIR]}/{input[0]} {config[WORKDIR]}/{config[SAME_MATURE]} {config[WORKDIR]}/{config[CLOSE_MATURE]} {config[WORKDIR]}/{config[SAME_PRE]} 2> {config[WORKDIR]}/{output[0]};"
        "mv result*.csv {config[WORKDIR]}/{output[1]}"

rule mirDeep2Virus:
    input:
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}.arf",
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-collapsed.fa",
        "data/processedData/mirDeep/miRNACombined-total.fa"
    output:
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-report.log",
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-results.csv"
    shell:
        "conda activate bioinfo;"
#UNCOMMENT BELOT LINE TO GENERATE VIRUS COMPARATIVE FASTA SEQUENCE
#        "cat {input[2]} {config[SAME_MATURE]} {config[CLOSE_MATURE]} > genomes/virusRefMiRNA.fa;"
        "cd {config[WORKDIR]}/data/processedData/mirDeepVirus/{wildcards.sampleVirus}/;"
        "miRDeep2.pl {config[WORKDIR]}/{input[1]} {config[WORKDIR]}/{config[VIRUS_GENOME]} {config[WORKDIR]}/{input[0]} none {config[WORKDIR]}/genomes/virusRefMiRNA.fa none 2> {config[WORKDIR]}/{output[0]};"
        "cp {config[WORKDIR]}/{output[0]} {wildcards.sampleVirus}-copy.log;"
        "mv result*.csv {config[WORKDIR]}/{output[1]}"

rule splitMirDeep2Results:
    input:
        "data/processedData/mirDeep/{sample}/{sample}-results.csv"
    output:
        "data/processedData/mirDeep/{sample}/results-{sample}-novel.csv",
        "data/processedData/mirDeep/{sample}/results-{sample}-known.csv"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeep/{wildcards.sample}/;"
        "csplit -f {wildcards.sample} {config[WORKDIR]}/{input} /miRDeep2/ '{{*}}';"
        "mv {wildcards.sample}03 {config[WORKDIR]}/{output[0]};"
        "mv {wildcards.sample}05 {config[WORKDIR]}/{output[1]}"
        
rule splitMirDeep2ResultsVirus:
    input:
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-results.csv"
    output:
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-novel.csv",
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-known.csv"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeepVirus/{wildcards.sampleVirus}/;"
        "csplit -f {wildcards.sampleVirus} {config[WORKDIR]}/{input} /miRDeep2/ '{{*}}';"
        "mv {wildcards.sampleVirus}03 {config[WORKDIR]}/{output[0]};"
        "mv {wildcards.sampleVirus}05 {config[WORKDIR]}/{output[1]}"

rule isolateMiRNAs:
    input:
        "data/processedData/mirDeep/{sample}/results-{sample}-novel.csv",
        "data/processedData/mirDeep/{sample}/results-{sample}-known.csv"
    output:
        "data/processedData/mirDeep/{sample}/{sample}-counts.csv",
        "data/processedData/mirDeep/{sample}/{sample}-novel.fa",
        "data/processedData/mirDeep/{sample}/{sample}-known.fa"
    shell:
        "conda activate bioinfo;"
        "Rscript scripts/rScripts/isolate.R"

rule isolateMiRNAsVirus:
    input:
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-novel.csv",
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-known.csv"
    output:
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-counts.csv",
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-novel.fa",
        "data/processedData/mirDeepVirus/{sampleVirus}/{sampleVirus}-known.fa"
    shell:
        "conda activate bioinfo;"
        "Rscript scripts/rScripts/isolate.R"

rule miRNANovelCombine:
    input:
        expand("data/processedData/mirDeep/{samples}/{samples}-novel.fa", samples=config["SAMPLES"])
    output:
        "data/processedData/mirDeep/miRNACombined-novelWithDups-old.fa",
        "data/processedData/mirDeep/miRNACombined-novel-old.fa",
        "data/processedData/mirDeep/miRNACombined-novel.fa"
    shell:
        "conda activate bioinfo;"
        "cat {input} > {output[0]};"
        "seqkit rmdup -s {output[0]} > {output[1]};"
        "awk '/^>/{{print "">as-mir"" ++i; next}}{{print}}' < {output[1]} > {output[2]}"

rule miRNANovelCombineVirus:
    input:
        expand("data/processedData/mirDeepVirus/{samples}/{samples}-novel.fa", samples=config["SAMPLES"])
    output:
        "data/processedData/mirDeepVirus/miRNACombinedVirus-novelWithDups-old.fa",
        "data/processedData/mirDeepVirus/miRNACombinedVirus-novel-old.fa",
        "data/processedData/mirDeepVirus/miRNACombinedVirus-novel.fa"
    shell:
        "conda activate bioinfo;"
        "cat {input} > {output[0]};"
        "seqkit rmdup -s {output[0]} > {output[1]};"
        "awk '/^>/{{print "">mayv-mir"" ++i; next}}{{print}}' < {output[1]} > {output[2]}"

rule miRNAKnownCombine:
    input:
        expand("data/processedData/mirDeep/{samples}/{samples}-known.fa", samples=config["SAMPLES"])
    output:         
        "data/processedData/mirDeep/miRNACombined-knownWithDups.fa",
        "data/processedData/mirDeep/miRNACombined-known.fa"
    shell:
        "conda activate bioinfo;"
        "cat {input} > {output[0]};"
        "seqkit rmdup -n {output[0]}  > {output[1]}"

rule miRNAKnownCombineVirus:
    input:
        expand("data/processedData/mirDeepVirus/{samples}/{samples}-known.fa", samples=config["SAMPLES"])
    output:
        "data/processedData/mirDeepVirus/miRNACombinedVirus-knownWithDups.fa",
        "data/processedData/mirDeepVirus/miRNACombinedVirus-known.fa"
    shell:
        "conda activate bioinfo;"
        "cat {input} > {output[0]};"
        "seqkit rmdup -n {output[0]}  > {output[1]}"

rule miRNATotalCombine:
    input:
        "data/processedData/mirDeep/miRNACombined-novel.fa",
        "data/processedData/mirDeep/miRNACombined-known.fa"
    output:
        "data/processedData/mirDeep/miRNACombined-total.fa"
    shell:
        "cat {input} > {output}"

rule miRNATotalCombineVirus:
    input:
        "data/processedData/mirDeepVirus/miRNACombinedVirus-novel.fa",
        "data/processedData/mirDeepVirus/miRNACombinedVirus-known.fa"
    output:
        "data/processedData/mirDeepVirus/miRNACombinedVirus-total.fa"
    shell:
        "cat {input} > {output}"

rule renameNovelmiRNA:
    input:
        "data/processedData/mirDeep/{sample}/{sample}-novel.fa",
        "data/processedData/mirDeep/miRNACombined-novel.fa"
    output:
        "data/processedData/mirDeep/{sample}/{sample}-novelRenamed.fa"
    shell:
        "grep -f {input[0]} {input[1]} -B1 | grep -v -- ""^--$"" > {output}"


rule combineMiRNAs:
    input:
        "data/processedData/mirDeep/{sample}/{sample}-novelRenamed.fa",
        "data/processedData/mirDeep/{sample}/{sample}-known.fa"
    output:
        "data/processedData/mirDeep/{sample}/{sample}-combined.fa"
    shell:
        "cat {input[0]} {input[1]} > {output}"

rule isolateCondenseReadCounts:
    input:
        "data/processedData/mirDeep/{sample}/results-{sample}-novel.csv",
        "data/processedData/mirDeep/{sample}/results-{sample}-known.csv",
        "data/processedData/mirDeep/miRNACombined-novel.fa"
    output:
        "data/processedData/condensedOutput/{sample}-novel.csv",
        "data/processedData/condensedOutput/{sample}-known.csv"
    script:
        "scripts/rScripts/condenseReads.R"

rule isolateCondenseReadCountsVirus:
    input:
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-novel.csv",
        "data/processedData/mirDeepVirus/{sampleVirus}/results-{sampleVirus}-known.csv",
        "data/processedData/mirDeepVirus/miRNACombinedVirus-novel.fa"
    output:
        "data/processedData/condensedOutputVirus/{sampleVirus}-novelVirus.csv",
        "data/processedData/condensedOutputVirus/{sampleVirus}-knownVirus.csv"
    script:
        "scripts/rScripts/condenseReads.R"

rule IDTargetSitesTotal:
    input:
        "data/processedData/mirDeep/miRNACombined-novel.fa",
        "data/processedData/mirDeep/miRNACombined-known.fa"
    output:
        "data/processedData/mirDeep/miRandaTargetSites-totalNovel.txt",
        "data/processedData/mirDeep/miRandaTargetSites-totalknown.txt"
    shell:
        "conda activate bioinfo;"
        "miranda {input[0]} {config[GENOME]} -out {output[0]} -strict -sc 140 -go -9 -ge -4 -en -20 -quiet;"

rule parseMiRandaTotal:
    input:
        "data/processedData/mirDeep/miRandaTargetSites-totalNovel.txt",
        "data/processedData/mirDeep/miRandaTargetSites-totalknown.txt"
    output:
        "data/processedData/mirDeep/miRandaTargetSites-totalNovel-parsed.txt",
        "data/processedData/mirDeep/miRandaTargetSites-totalKnown-parsed.txt"
    shell:
        "grep -A 1 'Scores for this hit:' {input[0]} | sort | grep '>' > data/processedData/mirDeep/1.tmp;"
        "cat header.txt data/processedData/mirDeep/1.tmp > {output[0]};"
        "grep -A 1 'Scores for this hit:' {input[1]} | sort | grep '>' > data/processedData/mirDeep/2.tmp;"
        "cat header.txt data/processedData/mirDeep/2.tmp > {output[1]};"
        "rm data/processedData/mirDeep/*.tmp"

rule retreiveAnnotationsTotal:
    input:
        "data/processedData/mirDeep/miRandaTargetSites-totalNovel-parsed.txt",
        "data/processedData/mirDeep/miRandaTargetSites-totalKnown-parsed.txt"
    output:
        "data/processedData/mirDeep/miRandaTargetSitesNovel-annotationsRaw.txt",
        "data/processedData/mirDeep/miRandaTargetSitesNovel-annotations.txt",
        "data/processedData/mirDeep/miRandaTargetSitesKnown-annotationsRaw.txt",
        "data/processedData/mirDeep/miRandaTargetSitesKnown-annotations.txt"
    shell:
        "conda activate bioinfo;"
        "awk '{{$2=$3=$5=$6=$7=$8=""; print $0}}' {input[0]} > {output[0]};"
        "bedtools intersect -a {config[BED]} -b {output[0]} -wa -wb > {output[1]};"
        "awk '{{$2=$3=$5=$6=$7=$8=""; print $0}}' {input[1]} > {output[2]};"
        "bedtools intersect -a {config[BED]} -b {output[2]} -wa -wb > {output[3]}"

rule removeMiRNAsFromRawDataAndSizeFilter:
    input:
        "data/processedData/mirDeep/{sample}/{sample}.fasta",
        "data/processedData/mirDeep/miRNACombined-total.fa"
    output:
        "data/processedData/piRNA/{sample}/{sample}-miRNARemoved.fa",
        "data/processedData/piRNA/{sample}/{sample}-FilteredForAlign.fa"
    shell:
        "conda activate bioinfo;"
        "cat {input[0]} | seqkit grep -ivsf <(seqkit seq -s --rna2dna {input[1]}) > {output[0]};"
        "seqkit seq -m 24 -M 30 {output[0]} > {output[1]}"

rule piRNAStarMap:
    input:
        "data/processedData/piRNA/{sample}/{sample}-FilteredForAlign.fa"
    output:
        "data/processedData/piRNA/{sample}/star_output/{sample}-piRNA.bam"
    shell:
        "conda activate bioinfo;"
        "STAR --runThreadN 8 --genomeDir genomes --readFilesIn {input} --sjdbGTFfile {config[GENOME_ANNOTATION]} --outFileNamePrefix  data/processedData/piRNA/{wildcards.sample}/star_output/ --outFilterMismatchNmax 3;"
        "samtools view -S -b data/processedData/piRNA/{wildcards.sample}/star_output/Aligned.out.sam | samtools sort > {output};"
        "samtools index {output}"

rule piRNAStarMapVirus:
    input:
        "data/processedData/piRNA/{sample}/{sample}-FilteredForAlign.fa"
    output:
        "data/processedData/virus/{sample}/star_output/{sample}-piRNA-virus.bam"
    shell:
        "conda activate bioinfo;"
        "STAR --runThreadN 8 --genomeDir genomes/virus --readFilesIn {input} --outFileNamePrefix  data/processedData/virus/{wildcards.sample}/star_output/ --outFilterMismatchNmax 3;"
        "samtools view -S -b data/processedData/virus/{wildcards.sample}/star_output/Aligned.out.sam | samtools sort > {output};"
        "samtools index {output}"

rule smallRNATotalMap:
    input:
        "data/processedData/mirDeep/{sample}/{sample}.fasta"
    output:
        "data/processedData/totalMap/{sample}/star_output/{sample}.bam"
    shell:
        "conda activate bioinfo;"
        "STAR --runThreadN 8 --genomeDir genomes --readFilesIn {input} --sjdbGTFfile {config[GENOME_ANNOTATION]} --outFileNamePrefix data/processedData/totalMap/{wildcards.sample}/star_output/ --outFilterMismatchNmax 3;"     
        "samtools view -S -b data/processedData/totalMap/{wildcards.sample}/star_output/Aligned.out.sam | samtools sort > {output};"
        "samtools index {output}"

rule totalStarMapVirus:
    input:
        "data/processedData/mirDeep/{sample}/{sample}.fasta"
    output:
        "data/processedData/virusTotal/{sample}/star_output/{sample}-virus.bam"
    shell:
        "conda activate bioinfo;"
        "STAR --runThreadN 8 --genomeDir genomes/virus --readFilesIn {input} --outFileNamePrefix  data/processedData/virusTotal/{wildcards.sample}/star_output/ --outFilterMismatchNmax 3;"
        "samtools view -S -b data/processedData/virusTotal/{wildcards.sample}/star_output/Aligned.out.sam | samtools sort > {output};"
        "samtools index {output}"
