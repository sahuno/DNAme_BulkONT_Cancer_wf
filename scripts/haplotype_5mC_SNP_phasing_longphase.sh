## longphase example

lphase=/data1/greenbab/users/ahunos/apps/longphase_linux-x64
bams=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam
save_dir=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/longphase_5mc_haplotypes
ref_file=/data1/greenbab/database/mm10/mm10.fa

$lphase modcall -b $bams -t 8 -o modcall_D-0-1_5000_4000 -r /data1/greenbab/database/mm10/mm10.fa




INPUT_DIR="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/longphase_5mc_haplotypes"      # e.g. /home/user1/output (absolute path needed)
REF_DIR="/data1/greenbab/database/mm10"      # e.g. /home/user1/reference (absolute path needed)
SOFTWARES_DIR="/data1/greenbab/users/ahunos/apps/containers"      # e.g. /home/user1/software (absolute path needed)
PLATFORM='ont_r10_dorado_sup_5khz'
THREADS=12
mkdir -p ${OUTPUT_DIR}


# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR},${REF_DIR},${SOFTWARES_DIR} \
  ${SOFTWARES_DIR}/clairs-to_latest.sif \
  /opt/bin/run_clairs_to \
  --tumor_bam_fn ${INPUT_DIR}/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam \
  --ref_fn ${REF_DIR}/mm10.fa \
  --threads ${THREADS} \
  --platform ${PLATFORM} \
  --output_dir ${OUTPUT_DIR} \
  --ctg_name="chr19" \
  --conda_prefix /opt/micromamba/envs/clairs-to



snpFile=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/longphase_5mc_haplotypes/snv.vcf.gz

$lphase phase \
-s $snpFile \
--mod-file modcall_D-0-1_5000_4000.vcf \
-b $bams \
-r $ref_file \
-t 8 \
-o D-0-1_5000_4000 \
--ont # or --pb for PacBio Hifi

phased_mod=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/longphase_5mc_haplotypes/D-0-1_5000_4000_mod.vcf

# haplotag the phased mod calls; can add svs and s
# $lphase haplotag \
# --mod-file=$phased_mod 
# -r $ref_file \
# -b $bams \
# -t 8 \
# -o D-0-1_5000_4000_haplotagged

#ideally run
$lphase haplotag \
-s phased_snp.vcf \
--sv-file phased_sv.vcf \
--mod-file=$phased_mod 
-r $ref_file \
-b $bams \
-t 8 \
-o D-0-1_5000_4000_haplotagged


#subset the bam file to only include the standard chromosomes
# samtools view -b $bams "chr1" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chrM" "chrX" "chrY"  > D-0-1_5000_4000_modBaseCalls_sorted_dup_standard_chrom.bam
# samtools index D-0-1_5000_4000_modBaseCalls_sorted_dup_standard_chrom.bam



#pipeline for longphase
#1. run modcall
#2. call snps with `clair3_tumor` enable phasing option 
#3. co-pahse snps and 5mc with `longphase phase`
#4. haplotag reads in bams the phased mod calls with `longphase haplotag`


#### try with d02
D02_bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-2_5000_4000/D-0-2_5000_4000_modBaseCalls_sorted_dup.bam
$lphase modcall -b $D02_bam -t 8 -o modcalls/modcall_D-0-2_5000_4000_modcall -r $ref_file

snpFileD02=/data1/greenbab/projects/methylRNA/Methyl2Expression/sandbox/longphase_test/results/call_snps_indels/D-0-2_5000/snv.vcf.gz

$lphase phase \
-s $snpFileD02 \
--mod-file modcalls/modcall_D-0-2_5000_4000_modcall.vcf \
-b $D02_bam \
-r $ref_file \
-t 8 \
-o phased_mod/D-0-2_5000_4000_longphase_phase \
--ont


$lphase haplotag \
-s phased_snp.vcf \
--sv-file phased_sv.vcf \
--mod-file=phased_mod/D-0-2_5000_4000_longphase_phase.vcf \
-r $ref_file \
-b $bams \
-t 8 \
-o D-0-2_5000_4000_haplotagged