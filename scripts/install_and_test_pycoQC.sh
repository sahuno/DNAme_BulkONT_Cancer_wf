mamba env remove -n pycoQC
mamba create -n pycoQC python=3.12
conda activate pycoQC

pip install pycoQC
# pip install --index-url https://test.pypi.org/simple/ pycoQC -U

# srun -p componc_cpu --pty --nodes=1 --ntasks-per-node=24 --mem=100G --time=04:00:00 bash


#store reesults here
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/test_dorado_summary
D03_4000_bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mod_bases/D-0-3_4000/D-0-3_4000_modBaseCalls_sorted.bam
D03_4000_summary=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/dorado_summary/D-0-3_4000/D-0-3_4000_seq_summary.txt

OUT_DIR=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/test_dorado_summary

dorado summary ${D03_4000_bam} > ${OUT_DIR}/D-0-3_4000_summary.tsv


pycoQC -f ${D03_4000_summary} -a ${D03_4000_bam} -o D03_4000_pycoQC_output_withAlignment.html

d01_5000_summary=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/dorado_summary/D-0-1_5000/D-0-1_5000_seq_summary.txt
d01_5000_chr8_bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/chr8_markdup_bams/D-0-1_5000_chr8.sorted.bam

dA1_4000_summary=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/dorado_summary/D-A-1_4000/D-A-1_4000_seq_summary.txt
DA1_4000_bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/chr8_markdup_bams/D-A-1_4000/D-A-1_4000_chr8.sorted.bam

#works
# pycoQC -f ${d01_5000_summary} -a ${d01_5000_chr8_bam} -o _chr8_pycoQC_output_withAlignment.html

pycoQC -f ${d01_5000_summary} ${dA1_4000_summary} -a ${d01_5000_chr8_bam} ${DA1_4000_bam} -o D01_5000_dA1_4000_chr8_pycoQC_output_withAlignment.html

pycoQC -f ${d01_5000_summary} ${dA1_4000_summary} -a ${d01_5000_chr8_bam} ${DA1_4000_bam} -o D01_5000_comb_chr8_pycoQC_output_withAlignment.html dA1_4000_comb_chr8_pycoQC_output_withAlignment.html

pycoQC -f /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/dorado_summary/D-0-1_5000/D-0-1_5000_seq_summary.txt -a /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/chr8_markdup_bams/D-0-1_5000_chr8.sorted.bam -o D01_5000_chr8_pycoQC_output_withAlignment.html




pycoQC -f ${D03_4000_summary} -a /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/chr8_markdup_bams/D-0-1_5000_chr8.sorted.bam /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/chr8_markdup_bams/D-A-1_4000/D-A-1_4000_chr8.sorted.bam -o D03_4000_pycoQC_output_withAlignment.html
