#!/bin/bash
set -e

module load dorado
module load modkit
module load samtools
module load modbam2bed

#get user args
POD5_FILE=$1
methyl_context=$2 #5mCG_5hmCG
basecall_model_file=$3 #/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0
reference_genome=$4
modified_calls_bam=$5
modified_calls_sorted_bam=$6
modified_calls_sorted_bam_bai=$7
CpG_5mC_5hmC_merged_bed=$8
CpG_5mC_bed=$9

#sampleName=$(basename -s .pod5 ${POD5_FILE})


#dorado download --model dna_r10.4.1_e8.2_400bps_sup@v4.1.0

echo "running dorado..."
dorado basecaller \
    $basecall_model_file \
    $POD5_FILE/ \
    --reference $reference_genome \
    --modified-bases $methyl_context \
    --verbose | \
    samtools view --threads 32 -O BAM -o $modified_calls_bam  && echo "done!" || echo "failed."
    

#sort bams and index
echo "sorting.."
samtools sort --threads 32 -o $modified_calls_sorted_bam $modified_calls_bam && echo "done." || echo "failed."

echo "indexing.."
samtools index $modified_calls_sorted_bam && echo "done." || echo "failed."



# ## methylation bed generation

modkit pileup $modified_calls_sorted_bam \
$CpG_5mC_bed \
--cpg --ref $reference_genome \
--ignore h

modkit pileup $modified_calls_sorted_bam \
$CpG_5mC_5hmC_merged_bed \
--cpg --ref $reference_genome --combine-mods 