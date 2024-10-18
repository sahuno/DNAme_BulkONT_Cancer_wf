#name:Samuel ahuno
#date: Sept 19th 2023
#create yaml files

import glob as glob
import yaml as yaml

#get paths
paths_sorted_bams = glob.glob("/lila/data/greenbaum/users/ahunos/TRI_EPI_DIVYA/mod_bases/*/*_sorted.bam",recursive=True)

sampleNames=[paths_sorted_bams[i].split("/")[-2] for i in range(len(paths_sorted_bams))] #get sample names
samples_dict=dict.fromkeys(sampleNames) #make dictionary with sample names as keys

for i in range(len(paths_sorted_bams)):
    samples_dict[paths_sorted_bams[i].split("/")[-2]] = paths_sorted_bams[i]

# dict_samples = {}
# dict_samples[0] = "samples"
# samples_dict.update(dict_samples)
#write output to file.yaml
with open(r'/lila/data/greenbaum/users/ahunos/TRI_EPI_DIVYA/scripts/samples_sorted_bams.yaml', 'w') as file:
    documents = yaml.dump(samples_dict, file)