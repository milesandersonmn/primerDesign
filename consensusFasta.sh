#!/bin/bash

conda activate samtools
while read region; do
while read sample; do
fasta=$(echo ${sample}.${region}.fa);
samtools faidx ../S_chilense_reference_rename.fasta $region | bcftools consensus -s $sample ../DFE/Schil_capture_variants_AllChroms.decomposedVariants.filtered.maxMissing.sorted.reheader.vcf.gz -o $fasta --mark-ins lc --mark-del '-' --mark-snv lc;
sed -i "s/^>\+.*/>$sample/g" $fasta;
cat $fasta >> $region.allSamples.fa;
rm -f $fasta;
done < primerSamples.txt;
done < primerLoci.txt




while read region; do
while read sample; do
fasta=$(echo ${sample}.${region}.fa);
samtools faidx ../S_chilense_reference_rename.fasta $region | bcftools consensus -s $sample ../DFE/Schil_capture_variants_AllChroms.decomposedVariants.filtered.maxMissing.sorted.reheader.vcf.gz -o $fasta --mark-ins lc --mark-del '-' --mark-snv lc;
sed -i "s/^>\+.*/>$sample/g" $fasta;
cat $fasta >> $region.allSamples.fa;
done < primerSamples.txt;
done < <(head -n 1 primerLoci.txt)
