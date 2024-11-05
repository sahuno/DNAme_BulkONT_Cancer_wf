from oct 31st
goal:
scalable, platform independent, parallel workflow to preprocess Bulk ONT DNA modification 

nextflow & snakemake workflow 


Question to the experts
how to assess phasing quality
sasha, julio, sobran andrew


todo
add optional arguments to rule shell commands

think about
what should we extract with modkit
modkit stats <input_bedfile> <bam> #get methylation stats per region
modkit stats <input_bedfile> <bam> #get genome methylation calls stats

modkit pileup <bam> <outBedmethyl> #aggregate methylation calls per modification
modkit pileup <bam> <outBedmethyl> #aggregate methylation calls per modification

Gotcha!
sort bedfiles with bedtools for outputs to work with bedtools map
`sortBed -i $bedNotZipped > output.bed`


# tuesdat nov 5
sample witout replacement to get run minimal test
awk 'rand()>0.99' /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5view/D-0-2/D-0-2_sampleRateSummary.tsv > D-0-2_sampleRateSummary_Sampledpercentage.tsv

sed -i '1s/^/read_id\tsample_rate\n/' D-0-2_sampleRateSummary_Sampledpercentage.tsv

less D-0-2_sampleRateSummary_Sampledpercentage.tsv
wc -l D-0-2_sampleRateSummary_Sampledpercentage.tsv


D02=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/pod5/D-0-2/D-0-2.pod5
pod5 subset $D02 --table D-0-2_sampleRateSummary_Sampledpercentage.tsv --columns sample_rate --output sample_rate_subset --template "D-0-2.{sample_rate}Hz.pod5"
POD5_DEBUG=1