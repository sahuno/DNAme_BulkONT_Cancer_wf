#make samplesDotYaml
import os
import glob
import yaml
import argparse
import pandas as pd
DirPath = "/data1/greenbab/projects/tcr_ont_DNAme/data/20231121_jurkat_tcr/E6/20231115_2101_3G_PAS57174_1472b895/pod5_pass/"
# DirPath = "/data1/greenbab/projects/tcr_ont_DNAme/data/20231121_jurkat_tcr/Jurkat/20231115_2101_3E_PAS88082_937c71ac/pod5_pass"

filenames = glob.glob(os.path.join(DirPath, "*.pod5"))
filenamesJk = glob.glob(os.path.join("/data1/greenbab/projects/tcr_ont_DNAme/data/20231121_jurkat_tcr/Jurkat/20231115_2101_3E_PAS88082_937c71ac/pod5_pass", "*.pod5"))

mergedFileNames = filenames + filenamesJk


#get samplee names
samples = {os.path.splitext(os.path.basename(f))[0]: f for f in filenames}# samples{0}
samplesJk = {os.path.splitext(os.path.basename(f))[0]: f for f in filenamesJk}# samples{0}
samplesMerged = {os.path.splitext(os.path.basename(f))[0]: f for f in mergedFileNames}# samples{0}

# samples.keys()[0:3]
samples_subset = dict(list(samples.items())[:5])
samplesMerged_subset = dict(list(samplesMerged.items())[:5])

samplesMergedPD = pd.DataFrame(samplesMerged.items(), columns=["readId", "path"])
samplesMergedPD["sample"] = samplesMergedPD["path"].str.split("/").str[-4]
samplesMergedPD["conditionTumorNormal"] = "cellLines"
samplesMergedPD['patientID'] = samplesMergedPD['sample'].map({
    'E6': 'PID001',
    'Jurkat': 'PID002'
})

#sample a few to test
# sampled_paths = samplesMergedPD.groupby('patientID').sample(n=5, random_state=42)
# sampled_paths = sampled_paths.reset_index(drop=True)
# sampled_paths.to_csv("/data1/greenbab/projects/tcr_ont_DNAme/scripts/configs/sampled_pathsTcells.csv", index=False, sep="\t")

# sampled_paths = samplesMergedPD.reset_index(drop=True)
samplesMergedPD.to_csv("/data1/greenbab/projects/tcr_ont_DNAme/scripts/configs/samples_pathsTcells.csv", index=False, sep="\t")

# samplesMerged_subsetPD[samplesMerged_subsetPD["Path"].str.contains("Jurkat")]
# pd.DataFrame(samples.items(), columns=["readIds", "Path"])

# Wrap samples dictionary under 'samples' key
output_dict = {"samples": samples_subset}



# Save to YAML file
output_yaml = "head_E6_samples.yaml"

with open(output_yaml, "w") as yaml_file:
    yaml.dump(output_dict, yaml_file, default_flow_style=False, sort_keys=False)

print(f"YAML saved to {output_yaml}")




