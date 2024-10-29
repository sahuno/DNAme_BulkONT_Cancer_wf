longphase modcall -b /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-0-1_4000/D-0-1_4000_modBaseCalls_sorted_dup.bam -t 8 -o results/longphase_modcall/D-0-1_4000/modcall_D-0-1_4000.vcf -r /data1/greenbab/database/mm10/mm10.fa 2> logs/longphase_modcall/D-0-1_4000/D-0-1_4000.log
mv results/longphase_modcall/D-0-1_4000/modcall_D-0-1_4000.vcf.vcf results/longphase_modcall/D-0-1_4000/modcall_D-0-1_4000.vcf

# ssend numpyro lgphase /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/test_longphase.sh
