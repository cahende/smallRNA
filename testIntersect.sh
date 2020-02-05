#!/bin/bash

# -------------------------------------------------- #
# Test bedTools intersect command before integration into snakemake pipeline
# -------------------------------------------------- #

conda activate bioinfo

cd ~/$HOME/scratch/stephensiSmallRNA-november2018-mayvInf

gtf2bed genomes/*.gtf > genomes/*.bed

bedtools intersect -a genomes/*.bed -b data.bed
