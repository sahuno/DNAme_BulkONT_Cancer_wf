#Run with
# checkForMethylTags.sh /data1/shahs3/users/schrait/analyses/ont_temp/45665/45665.merged.bam

samtools view @1 | grep -E '\tMN:|\tML:'

# samtools view /data1/shahs3/users/schrait/analyses/ont_temp/45665/45665.merged.bam | grep -E '\tMN:|\tML:'


