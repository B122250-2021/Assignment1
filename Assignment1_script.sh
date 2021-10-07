#B122250

#With one file
#Directory for fastqc results
mkdir fastqc_results

#-f flag used, to specify the input format, in this case fastq
# For future versions, I might want to use the -t flag,to specify the threads
fastqc -f fastq /localdisk/home/data/BPSM/AY21/fastq/100k.C1-1-501_1.fq.gz -o $PWD/fastqc_results

#This script needs a section to produce the fastqc report, would go here

#to allign the reads, I need to make an index of the reference genome in a directory called index_files

mkdir index_files
bowtie2-build /localdisk/home/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz index_files/index_ref_file

#the allignment itself, with output in alligned_reads
mkdir alligned_reads
bowtie2 -x index_files/index_ref_file -1 /localdisk/home/data/BPSM/AY21/fastq/100k.C1-1-501_1.fq.gz -2 /localdisk/home/data/BPSM/AY21/fastq/100k.C1-1-501_2.fq.gz -S alligned_reads/100k.C1-1-501.aligned.sam

#convert output to indexed bam
samtools view -b alligned_reads/100k.C1-1-501.aligned.sam  > alligned_reads/100k.C1-1-501.aligned.bam && samtools sort alligned_reads/100k.C1-1-501.aligned.bam > alligned_reads/100k.C1-1-501.aligned.sorted.bam && samtools index alligned_reads/100k.C1-1-501.aligned.sorted.bam

# sam to bam and then count gene reads using bedtools coverage, which outputs text file into folder called counts

mkdir counts
bedtools bamtobed -i alligned_reads/100k.C1-1-501.aligned.sorted.bam > alligned_reads/100k.C1-1-501.aligned.sorted.bed
bedtools coverage -counts -a /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b alligned_reads/100k.C1-1-501.aligned.sorted.bed > counts/100k.C1-1-501.counts.txt

#Now I want to work on making a file that shows the gene, gene description and then the counts for each file
# make a file called gene_expressions.txt, add headers and the contents	form columns 4,5 in the	bed file
touch gene_expressions.txt
echo "Gene	Gene_description">>gene_expressions.txt
cut -f 4,5 /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed >>gene_expressions.txt

#Now I want to merge the gene count file, but need to sort both
sort -t$'\t' -k1,1 gene_expressions.txt
sort -t$'\t' -k4,4 counts/100k.C1-1-501.counts.txt > counts/100k.C1-1-501.counts.sorted.txt
cut -f 4,6 counts/100k.C1-1-501.counts.sorted.txt > temp_count_file.txt
file="100k.C1-1-501.counts.sorted.txt"
sed -i "1i Gene ""$file" temp_count_file.txt
#unset file
join -t$'\t' -e 'NaN' gene_expressions.txt temp_count_file.txt > gene_expressions_test.txt
