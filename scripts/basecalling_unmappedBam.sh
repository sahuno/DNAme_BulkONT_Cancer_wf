##run 3: run basecaller to unaligned reads. 
dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --verbose --emit-sam > ONTwf_5mC5hmC6mA_D01_5000kHz_withMoves.sam

# sgpu numpyro D01_5000kHzGPU5 /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/basecalling_unmappedBam.sh 5

