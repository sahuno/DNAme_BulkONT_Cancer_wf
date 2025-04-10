#author: samule ahuno 
#make samples_conditions.csv for merging bedmethyl files files
import os
from pathlib import Path
import glob as glob
import pandas as pd
import yaml


DIR = "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/modkit/"
filenames=glob.glob(os.path.join(DIR, "**","*_modpileup_combined.bed.gz"))

samples = {os.path.basename(f).replace("_modpileup_combined.bed.gz", ""): f for f in filenames}
# samples = {os.path.splitext(os.path.basename(f))[0]: f for f in filenames}# samples{0}
samples

# samples = dict(list(samples.items())[:5])
samplesDF = pd.DataFrame(samples.items(), columns=["samples", "path"])

# Load mapping from provided CSV file
mapping_file = "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_DNArecoded.csv"
mapping_df = pd.read_csv(mapping_file)

# Create a mapping dictionary based on the 'samples' column
condition_map = mapping_df.set_index("samples")[["condition", "condition_long"]]

# Extract base sample names (e.g., 'D-0-1' from 'D-0-1_5000_4000')
samplesDF["sample_base"] = samplesDF["samples"].str.extract(r"(D-[^-]+-\d+)")

samplesDF = samplesDF.join(condition_map, on="sample_base")
samplesDF.drop(columns=["sample_base"], inplace=True)

#write to dd
samplesDF.to_csv("samples_per_conditions.csv", index=False, sep="\t")


condition_to_paths = samplesDF.groupby("condition")["path"].apply(list).to_dict()
condition_to_paths


# Save to YAML file
output_dict = {"samples": condition_to_paths}
output_yaml = "samplesByConditions.yaml"

with open(output_yaml, "w") as yaml_file:
    yaml.dump(output_dict, yaml_file, default_flow_style=False, sort_keys=False)

print(f"YAML saved to {output_yaml}")

