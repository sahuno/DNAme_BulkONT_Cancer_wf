#keep things to mm10 if you don't have anyway to build other ref databases dfor mm39
#  --emit-summary doesn't work with `dorado basecaller` yet
# - Note: FASTQ output is not recommended as not all data can be preserved. #why woyl dyou keep fastq file

#run 1: not what needed
# dorado basecaller sup,5mCG_5hmCG@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --emit-moves --output-dir "ONTwf_test_D01_5000kHz" --mm2-opts "-Y"
# sgpu numpyro D01_5000kHz /data1/greenbab/users/ahunos/apps/ont_cancer_pipeline/scripts/test_modbase_calling_commands_dorado.sh 4

#run 2; not what needed
# dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --emit-moves --output-dir "ONTwf_5mC5hmC6mA_D01_5000kHz" --mm2-opts "-Y"
# l /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5

##run 3: this is what i need i think. move tables
# dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --emit-moves --mm2-opts "-Y" --output-dir ONTwf_5mC5hmC6mA_D01_5000kHz_withMoves #--resume-from ONTwf_5mC5hmC6mA_D01_5000kHz/calls_2024-10-15_T16-35-20.sam > ONTwf_5mC5hmC6mA_D01_5000kHz_withMoves.sam

# dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --emit-moves --mm2-opts ""-Y"" --resume-from ONTwf_5mC5hmC6mA_D01_5000kHz/calls_2024-10-15_T16-35-20.sam > ONTwf_5mC5hmC6mA_D01_5000kHz_withMoves.sam
# sgpu numpyro D01_5000kHzGPU5 /data1/greenbab/users/ahunos/apps/ont_cancer_pipeline/scripts/test_modbase_calling_commands_dorado.sh 5


# scontrol show job
# dorado basecaller sup /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --verbose --emit-fastq --emit-moves --output-dir "fastq/D01_5000kHz"



##run 4: this is without move tables
### i want to know how many sam files are written in the --outDir? shouldn't be 2
dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --output-dir "ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves_rerun" --mm2-opts "-Y"
# sgpu numpyro D01_5000kHzG8 /data1/greenbab/users/ahunos/apps/ont_cancer_pipeline/scripts/test_modbase_calling_commands_dorado.sh 8


#sort small file
# samtools view -u /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T03-29-52.sam | samtools sort @24 m 8G - > /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T03-29-52.bam
# samtools index /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T03-29-52.sam
# # sacct -j 6732262 -X --format=JobID,Partition%20,Account,AllocCPUS,State,ExitCode
# # scontrol show job 6732262 

#sort big file
# samtools view -u /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T01-53-48.sam | samtools sort | samtools view -h > /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T01-53-48.bam
# samtools index /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T01-53-48.bam
# # samtools index /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T01-53-48.sam

#index big bam file 
# sambamba index -p /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/ONTwf_test_D01_5000kHz_v2/ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves/calls_2024-10-16_T01-53-48.bam calls_2024-10-16_T01-53-48.bam.bai
# samtools index calls_2024-10-16_T01-53-48.bam
# scontrol show job 6732262  

# sambamba view --sam-input calls_2024-10-16_T01-53-48.sam > calls_2024-10-16_T01-53-48_unsorted.bam
# sambamba sort -p -t 12 calls_2024-10-16_T01-53-48_unsorted.bam calls_2024-10-16_T01-53-48_sorted.bam

# sambamba index -p calls_2024-10-16_T01-53-48_unsorted.bam calls_2024-10-16_T01-53-48_unsorted.bam.bai




# samtools view -bS  calls_2024-10-16_T01-53-48.sam > calls_2024-10-16_T01-53-48_unsorted.bam
# ssend numpyro samtoolsSort /data1/greenbab/users/ahunos/apps/ont_cancer_pipeline/scripts/test_modbase_calling_commands_dorado.sh

#--no-home 
# singularity exec -B /data1/greenbab --nv /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --mm2-opts "-Y" #| samtools sort --threads 12 -o ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves_singularity_sorted.bam -O BAM --write-index


# dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --output-dir "ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves_rerun" --mm2-opts "-Y"
