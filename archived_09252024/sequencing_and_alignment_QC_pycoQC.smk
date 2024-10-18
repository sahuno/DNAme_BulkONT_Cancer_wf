parent_dir = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
# configfile: parent_dir + "/config/samples_mergeBam_modkit.yaml" #human samples
configfile: parent_dir + "/config/samples_bams_pycoqc.yaml" #mouse samples
set_species = "mouse"

#mouse directory
#$ /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit

#launch the pipeline by runnign the command below
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/mergeBams_mdup_modkit.smk --cores 5 -np
##2.  on cluster
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/sequencing_and_alignment_QC_pycoQC.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 1000 --cores all --use-conda --keep-going --forceall -np
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/sequencing_and_alignment_QC_pycoQC.smk --jobs 1000 --cores all --use-conda --keep-going --forceall -np

# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/sequencing_and_alignment_QC_pycoQC.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np

rule all:
    input: 
        expand("results/pycoQC_html/{samples}/{samples}_seq_n_Alignment_pycoQC.html", samples=config["samples"]),
        # expand("results/pycoQC_html/{samples}/{samples}_done.txt", samples=config["samples"]),
        expand("results/pycoQC_json/{samples}/{samples}_seq_n_Alignment_pycoQC.json", samples=config["samples"])
        # expand("results/pycoQC_json/{samples}/{samples}_done.txt", samples=config["samples"])


print(f"config['samples']\n {config["samples"]}")

rule pycoQC_html:
    input:
        bamsAll=lambda wildcards: config["samples"][wildcards.samples]
        # bam1=lambda wildcards: config["samples"][wildcards.samples][0]
    output:
         out_html="results/pycoQC_html/{samples}/{samples}_seq_n_Alignment_pycoQC.html"
        #  done_txt="results/pycoQC_html/{samples}/{samples}_done.txt"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=16
    conda: "pycoQC"
    log:
        "logs/pycoQC_html/{samples}/{samples}.log"
    shell:
        """ 
        echo "pycoQC bams for {wildcards.samples}"
        pycoQC -f {input.bamsAll[1]} -a {input.bamsAll[0]} -o {output.out_html} 2> {log}
        """

rule pycoQC_json:
    input:
        bamsAll=lambda wildcards: config["samples"][wildcards.samples]
        # bam1=lambda wildcards: config["samples"][wildcards.samples][0]
    output:
         out_json="results/pycoQC_json/{samples}/{samples}_seq_n_Alignment_pycoQC.json"
        #  done_txt_json="results/pycoQC_json/{samples}/{samples}_done.txt"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=16
    conda: "pycoQC"
    log:
        "logs/pycoQC_json/{samples}/{samples}.log"
    shell:
        """ 
        echo "pycoQC bams for {wildcards.samples}"
        pycoQC -f {input.bamsAll[1]} -a {input.bamsAll[0]} --json_outfile {output.out_json} 2> {log}
        """
# pycoQC -f ${d01_5000_summary} ${dA1_4000_summary} -a ${d01_5000_chr8_bam} ${DA1_4000_bam} -o D01_5000_comb_chr8_pycoQC_output_withAlignment.html dA1_4000_comb_chr8_pycoQC_output_withAlignment.html


        # if [ ${{input.bamsAll}} -gt 1 ]; then
        #     pycoQC -f {input.bamsAll[1]} -a {input.bamsAll[0]} -o {output.out_html} && touch {output.done_txt} 2> {log}
        # else
        #     touch {output.done_txt} 2> {log}
        # fi        