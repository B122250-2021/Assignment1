#B122250

#With many files
#Directory for fastqc results, index files, alligned_reads and counts
mkdir fastqc_results
mkdir index_files
mkdir alligned_reads
mkdir counts
mkdir bamfiles


#to allign the reads, I need to make an index of the reference genome in a directory called index_files, only need to do this once, so outside$
bowtie2-build /localdisk/home/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz index_files/index_ref_file

#now looping for every file
for file in /localdisk/home/data/BPSM/AY21/fastq/*.fq.gz
do
#-f flag used, to specify the input format, in this case fastq
# For future versions, I might want to use the -t flag,to specify the threads
fastqc -t 6 -f $file -o $PWD/fastqc_results
done

#loop to produce file with filenames
touch filenames.txt
for file in /localdisk/home/data/BPSM/AY21/fastq/*.fq.gz
do
echo "${file:37}" >> filenames.txt
done

#This script needs a section to produce the fastqc report, would go here

cat forward_and_reverse_reads.txt | while read line
do
input_file_1=$(echo ${line:0:21})
input_file_2=$(echo ${line:21})
sam_name=${line:0:13}
#the allignment itself, with output in alligned_reads
bowtie2 -x index_files/index_ref_file -1 /localdisk/home/data/BPSM/AY21/fastq/input_file_1 -2 /localdisk/home/data/BPSM/AY21/fastq/input_file_2 -S alligned_reads/$sam_name.aligned.sam
echo "$sam_name" >> sam_filenames.txt
done

#convert output to indexed bam
cat sam_filenames.txt | while read samfile
do
samtools view -b alligned_reads/$samfile.aligned.sam  > bamfiles/$samfile.aligned.bam && samtools sort bamfiles/$samfile.aligned.bam > bamfiles/$samfile.aligned.sorted.bam && samtools index bamfiles/$samfile.aligned.sorted.bam

# sam to bam and then count gene reads using bedtools coverage, which outputs text file into folder called counts

bedtools bamtobed -i bamfiles/$samfile.aligned.sorted.bam > bamfiles/$samfile.aligned.sorted.bed
bedtools coverage -counts -a /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b bamfiles/$samfile.aligned.sorted.bed > counts/$samfile.counts.txt
done
