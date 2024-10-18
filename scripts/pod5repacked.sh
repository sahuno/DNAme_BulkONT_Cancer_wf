pod5 repack /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 \
-o /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5Repacked/D-0-1/sampleRate5000/
# ssend pod5 pod5Repack /data1/greenbab/users/ahunos/apps/ont_cancer_pipeline/scripts/pod5repacked.sh

# mamba activate pod5

# POD5_DEBUG=1 
pod5 subset /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 \
-s /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5view/D-0-1/D-0-1_sampleRateSummary.tsv \
-o /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5Repacked/D-0-1/sampleRate5000subset/ \
--READ_ID_COLUMN read_id
# less /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5view/D-0-1/D-0-1_sampleRateSummary.tsv


# pod5 subset --help
# sampleRate5000subset
# pip install pod5
# pip install --upgrade pod5
