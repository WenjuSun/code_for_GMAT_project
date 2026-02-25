# sam to sorted bam and index 
for i in `ls -1 *mESC*HiC*.sam|sort`; do echo $i; samtools view --threads 16 -bS $i | samtools sort --threads 16 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
