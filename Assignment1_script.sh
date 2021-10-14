#B122250

#With many files
#Directory for fastqc results, index files, alligned_reads and counts
mkdir fastqc_results && mkdir index_files && mkdir alligned_reads && mkdir counts && mkdir bamfiles && mkdir final_counts

#to allign the reads, I need to make an index of the reference genome in a directory called index_files, only need to do this once, so outside$
bowtie2-build --threads 60 /localdisk/home/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz index_files/index_ref_file

#now looping for every file
for file in /localdisk/home/data/BPSM/AY21/fastq/*.fq.gz
do
#-f flag used, to specify the input format, in this case fastq
# For future versions, I might want to use the -t flag,to specify the threads
fastqc -t 60 -f fastq $file -o $PWD/fastqc_results
done

#making a file that has both forwards and reverse reads on the same line, so that I can use a single variable for both
touch forward_and_reverse_reads.txt && touch forward_reads.txt && touch reverse_reads.txt
for file in /localdisk/home/data/BPSM/AY21/fastq/*1.fq.gz
do
echo "${file:37}" >> forward_reads.txt
done
for file in /localdisk/home/data/BPSM/AY21/fastq/*2.fq.gz
do
echo "${file:37}" >> reverse_reads.txt
done
paste -d "" forward_reads.txt reverse_reads.txt > forward_and_reverse_reads.txt

#This script needs a section to produce the fastqc report, would go here

#Now I want to work on making a file that shows the gene, gene description and then the counts for each file. This goes before the loop, will need it later for the final table of counts
# make a file called gene_expressions.txt, add headers and the contents form columns 4,5 in the original bed file
touch genes.txt
echo "Gene	Gene_description">>genes.txt
cut -f 4,5 /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed >>genes.txt
#make a copy for the final gene counts
cp genes.txt final_counts/gene_counts_0.txt

cat forward_and_reverse_reads.txt | while read line
do
input_file_1=$(echo ${line:0:21})
input_file_2=$(echo ${line:21})
sam_name=${line:0:13}
#the allignment itself, with output in alligned_reads
bowtie2 -p 60 -x index_files/index_ref_file -1 /localdisk/home/data/BPSM/AY21/fastq/$input_file_1 -2 /localdisk/home/data/BPSM/AY21/fastq/$input_file_2 -S alligned_reads/$sam_name.aligned.sam
echo "$sam_name" >> sam_filenames.txt
done

i=0
p=1

#convert output to indexed bam
cat sam_filenames.txt | while read samfile
do
samtools view -b alligned_reads/$samfile.aligned.sam  > bamfiles/$samfile.aligned.bam && samtools sort bamfiles/$samfile.aligned.bam > bamfiles/$samfile.aligned.sorted.bam && samtools index bamfiles/$samfile.aligned.sorted.bam

# sam to bam and then count gene reads using bedtools coverage, which outputs text file into folder called counts

bedtools bamtobed -i bamfiles/$samfile.aligned.sorted.bam > bamfiles/$samfile.aligned.sorted.bed
bedtools coverage -counts -a /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b bamfiles/$samfile.aligned.sorted.bed > counts/$samfile.counts.txt

#Now I want to merge the gene count file, but need to sort both, as you need it before the merging
sort -t$'\t' -k1 genes.txt
sort -t$'\t' -k4 counts/$samfile.counts.txt > counts/$samfile.counts.sorted.txt
#cut the gene id and the count columns
cut -f 4,6 counts/$samfile.counts.sorted.txt > temp_count_file.txt
sed -i "1i Gene""	$samfile" temp_count_file.txt
join -t$'\t' -e 'NaN' genes.txt temp_count_file.txt > temp_gene_counts.txt

#I need to make sure that the files are the same length, before joining them

genelines1=$(wc -l  genes.txt | cut -d ' ' -f1)
genelines2=$(wc -l  temp_gene_counts.txt | cut -d ' ' -f1)

if [ $genelines1 != $genelines2 ]
then
   echo "Error in replicate "$samfile"" >> logfile.log
else
   echo "Expression of replicate "$samfile" was reccorded correctly" >> logfile.log
fi 

cut -f 3 temp_gene_counts.txt > temp_gene_counts_single_columns.txt

paste -d'\t' final_counts/gene_counts_$i.txt temp_gene_counts_single_columns.txt > final_counts/gene_counts_$p.txt
((i=i+1))
((p=p+1))
done
#Select the last file, which is the one containing all the gene counts
last_file=$(wc -l  sam_filenames.txt | cut -d ' ' -f1)
mv final_counts/*$last_file.txt replicate_gene_expressions.txt
#Remove messy files and directories
rm -fr forward_and_reverse_reads.txt && rm -fr forward_reads.txt && rm -fr reverse_reads.txt && rm -fr temp_gene_counts_single_columns.txt && rm -fr temp_gene_counts.txt && rm -fr temp_count_file.txt && rm -fr genes.txt && rm -fr sam_filenames.txt
rm -fr counts && rm -fr final_counts
