#!/bin/bash
file=$1
name=$(basename $file .bam)
fa=$2

#samtools stats
samtools stats -@ $3 $file > ${name}_samtools_stats.txt

#samtools flagstat
samtools flagstat -@ $3 $file > ${name}_samtools_flagstats.txt

#calmd
samtools calmd -@ $3 -b $file $fa > "${name}_tag.bam" 
samtools index -@ $3 "${name}_tag.bam"

#gapless
samtools view -h -@ $3 "${name}_tag.bam" | awk '$6 !~ /[IDNH]/ || $1 ~ /^@/' | samtools view -bo "${name}_gapless.bam" -
samtools index -@ $3 "${name}_gapless.bam"
samtools flagstat -@ $3 "${name}_gapless.bam" > "${name}_gapless_samtools_flagstats.txt"

#perfect
samtools view -h -@ $3 "${name}_gapless.bam" | awk 'BEGIN{FS="\t";OFS="\t"} $1 ~ /^@/ {print; next} {for(i=1;i<=NF;i++) if($i ~ /^MD:Z:[0-9]+$/) {print; break}}' | samtools view -bo "${name}_perfect.bam" -
samtools index -@ $3 "${name}_perfect.bam"
samtools flagstat -@ $3 "${name}_perfect.bam" > "${name}_perfect_samtools_flagstats.txt"


