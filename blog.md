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
