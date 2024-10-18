#i have 1 big pod5 with numerous reads with sampling rates. this messes up base calling. i want to split the pod5 into multiple pod5s based on the sampling rate.

dataset_pod5=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/pod5/D-S-1/D-S-1.pod5
pod5 view $dataset_pod5 --include "read_id, sample_rate" --output D-S-1_sampleRateSummary.tsv


pod5 subset --summary D-S-1_sampleRateSummary.tsv --columns sample_rate --output split_by_sample_rate $dataset_pod5
