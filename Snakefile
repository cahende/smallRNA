#!/usr/bin/env/python
shell.prefix("set -o pipefail; ")
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
configfile: "config.yaml"

rule all:
    input:
        expand("data/processedData/mirDeep/{sample}/{sample}-counts.csv", sample=config["SAMPLES"])

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
        "mv results*.csv {output[1]}"

rule splitMirDeep2Results:
    input:
        "data/processedData/mirDeep/{sample}/{sample}-results.csv"
    output:
        "data/processedData/mirDeep/{sample}/results-{sample}-novel.csv",
        "data/processedData/mirDeep/{sample}/results-{sample}-known.csv"
    shell:
        "conda activate bioinfo;"
        "cd {config[WORKDIR]}/data/processedData/mirDeep/{wildcards.sample}/;"
        "csplit -f {wildcards.sample} {config[WORKDIR]}/{input} '{{*}}';"
        "mv {wildcards.sample}01 {output[0]};"
        "mv {wildcards.sample}03 {output[1]}"

rule isolateMiRNAs:
    input:
        "data/processedData/mirDeep/{sample}/results-{sample}-novel.csv",
        "data/processedData/mirDeep/{sample}/results-{sample}-known.csv"
    output:
        "data/processedData/mirDeep/{sample}/{sample}-counts.csv",
        "data/processedData/mirDeep/{sample}/{sample}-novel.fa",
        "data/processedData/mirDeep/{sample}/{sample}-known.fa"
    script:
        "scripts/rScripts/isolate.R"
        
#Star map with modified paramaters: >=16b matched to the genome, number of mismatches <= 5% of mapped length, i.e. 0MM for 16-19b, 1MM for 20-39b etc, splicing switched off
rule star_map:
    input:
        expand("data/processedData/trimmed_reads/{{sample}}_{read}_paired.fastq", read=["R1", "R2"])
    output:
        directory("data/processedData/aligned_reads/star_output/{sample}/")
    log: "logs/{sample}.map_and_bam.log"
    shell:
        "conda activate bioinfo;"
        "fastqc {input};"
        "STAR --runThreadN 8 --genomeDir genomes --readFilesIn {input} --sjdbGTFfile {config[GENOME_ANNOTATION]} --outFileNamePrefix {output} \
                --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1"

rule sort_bam:
    input:
        "data/processedData/aligned_reads/star_output/{sample}/"
    output:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    log: "logs/{sample}.sort_bam.log"
    shell:
        "module load samtools fastqc;"
        "fastqc {input};"
        "samtools view -Sb {input}/Aligned.out.sam | samtools sort -o {output} --threads 8"

rule quality_filter_reads:
    input:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    output:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.passed.bam"
    log: "logs/{sample}.quality_filter_reads.log"
    shell:
        "module load bamtools fastqc;"
        "fastqc {input};"
        "bamtools filter -mapQuality '>=20' -length '<=40' -in {input} -out {output};"
        "fastqc {output}"

rule index_bam:
    input:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.passed.bam"
    output:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.passed.bam.bai"
    log: "logs/{sample}.index_bam.log"
    shell:
        "module load samtools fastqc;"
        "fastqc {input};"
        "samtools index {input} {output}"

rule fix_mate_pairs:
    input:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.passed.bam"
    output:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.passed.fixed.bam"
    log: "logs/{sample}.fix_mate_pairs.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} FixMateInformation INPUT={input} OUTPUT={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"

rule filter_mapped_and_paired_reads:
    input:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.passed.fixed.bam"
    output:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.passed.fixed.filtered.bam"
    log: "logs/{sample}.filter_mapped_and_paired_reads.log"
    shell:
        "module load bamtools fastqc;"
        "fastqc {input};"
        "bamtools filter -isMapped true -in {input} -out {output}"

rule remove_duplicate_reads:
    input:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.passed.fixed.filtered.bam"
    output:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.passed.fixed.filtered.postdup.bam"
    log: "logs/{sample}.remove_duplicate_reads.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} MarkDuplicates INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE={log}"

rule add_read_groups:
    input:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.passed.fixed.filtered.postdup.bam"
    output:
        "data/processedData/aligned_reads/read_group/{sample}.PE.star.sorted.passed.fixed.filtered.postdup.RG.bam"
    log: "logs/{sample}.add_read_groups.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={wildcards.sample}.PE RGPL=Illumina RGPU=Group1 RGSM={wildcards.sample}.PE"

rule rSubread:
    input:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.passed.fixed.filtered.postdup.RG.bam"
    output:
        "data/processedData/read_counts/{sample}.readCounts.txt"
    log: "logs/{sample}.rSubread.log"
    script:
        "scripts/rScripts/rSubreadFeatureCounts.R"

rule differential_expression:
    input:
        expand("data/processedData/read_counts/{sample}.readCounts.txt", sample=config["SAMPLES"])
    output:
        "data/processedData/differentialExpression/glmControlVsInfectedPValue.txt",
        "data/processedData/differentialExpression/topGOTerms-BP.csv",
        "data/processedData/differentialExpression/topGOTerms-MF.csv",
        "data/processedData/differentialExpression/topGOTerms-CC.csv"
    log: "logs/edgeR.log"
    script:
        "scripts/rScripts/edgeR.R"


