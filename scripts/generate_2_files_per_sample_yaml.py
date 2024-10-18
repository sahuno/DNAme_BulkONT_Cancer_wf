#name:Samuel ahuno
#date: Sept 19th 2023
#create sample yaml files. 2 files for each each 

import os
import yaml as yaml
import argparse
#inputs
#1. parent dir
#2. extension of files to search for
#4. output file name

# ext = "_5mC.bed"
# directory_to_search = '/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/modkit/D-Q-1_5000_4000/'

#function to find bed files
def find_files(directory, ext):
    files_output = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(ext):
                files_output.append(os.path.join(root, file))
    return files_output

def main(args):
    files_output = find_files(args.bamDir, args.ext_B)
    files_seqeuncing = find_files(args.seqSummmaryDir, args.ext_S)
    
    sampleNames=[files_output[i].split("/")[-2] for i in range(len(files_output))] #get sample names
    samples_dict=dict.fromkeys(sampleNames) #make dictionary with sample names as keys

    sampleNamesSeq=[files_seqeuncing[i].split("/")[-2] for i in range(len(files_seqeuncing))] #get sample names
    samples_dictSeq=dict.fromkeys(sampleNamesSeq) #make dictionary with sample names as keys


# Output: {'a': [1], 'b': [2, 4], 'c': [3, 5], 'd': [6]}



for i in range(len(files_output)):
    samples_dict[files_output[i].split("/")[-2]] = files_output[i]
for i in range(len(files_seqeuncing)):
    samples_dictSeq[files_seqeuncing[i].split("/")[-2]] = files_seqeuncing[i]

from collections import defaultdict

merged_dict = defaultdict(list)
# merged_dict["samples"] = []

for key, value in samples_dict.items():
    merged_dict[key].append(value)
    
for key, value in samples_dictSeq.items():
    merged_dict[key].append(value)

merged_dict2 = dict(merged_dict)
# merged_dict2[1] = "samples"
# merged_dict2[1] = "samples"
# print(dict(merged_dict))
    # dict_samples = {}
    # dict_samples[0] = "samples"
    # samples_dict.update(dict_samples)

if args.outputFile:
    with open(args.outputFile, 'w') as file:
        documents = yaml.dump(merged_dict2, file)
else:
    with open(f"samples_{args.ext}.yaml", 'w') as file:
        documents = yaml.dump(merged_dict2, file)

if __name__ == "__main__":
parser = argparse.ArgumentParser(description="make sample yaml file")
parser.add_argument("-B", "--bamDir", nargs="?", type=str)
parser.add_argument("-S", "--seqSummmaryDir", nargs="?", type=str)
parser.add_argument("--ext_B", nargs="?", default=".bed", type=str)
parser.add_argument("--ext_S", nargs="?", default=".bed", type=str)
parser.add_argument("--outputFile", nargs="?", type=str, help='file name to save to')
args = parser.parse_args()
main(args)

# Run this script with the following command:
# python make_generic_yamlSampleList.py -d /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed --ext "_5mC.bed" --outputFile samples_5mC.yaml

#make both 5mc and 5hmc
# python make_generic_yamlSampleList.py -d /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed --ext "_combined.bed" --outputFile samples_5mC_5hmC.yaml

#make both bamfiles
# python make_generic_yamlSampleList.py -B /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed -S /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed --ext_B "_modBaseCalls_sorted_dup.bam" --ext_S "_seq_summary.txt" --outputFile samples_bams_pycoqc.yaml
# args.bamDir="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed";
# args.seqSummmaryDir="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed" args.ext_B="_modBaseCalls_sorted_dup.bam" args.ext_S="_seq_summary.txt" args.outputFile="samples_bams_pycoqc.yaml"

# find /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed -name *_modBaseCalls_sorted_dup.bam
# find /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed -name *_seq_summary.txt
