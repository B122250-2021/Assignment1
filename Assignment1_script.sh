#B122250
#!/bin/bash
#With many files
#Directory for fastqc results, index files, aligned_reads and counts
mkdir fastqc_results && mkdir index_files && mkdir aligned_reads && mkdir counts && mkdir bamfiles && mkdir final_counts && mkdir fastqc_reports && mkdir groupwise_comparisons

#to align the reads, I need to make an index of the reference genome in a directory called index_files, only need to do this once, so outside$
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

#Making the FASTQC report
counter_0=0
counter_1=1
echo "Sample" >> fastqc_categories.txt && echo "Basic Statistics" >> fastqc_categories.txt && echo "Per base sequence quality" >> fastqc_categories.txt && echo "Per sequence quality scores" >> fastqc_categories.txt && echo "Per base sequence content" >> fastqc_categories.txt && echo "Per sequence GC content" >> fastqc_categories.txt && echo "Per base N content" >> fastqc_categories.txt && echo "Sequence Length Distribution" >> fastqc_categories.txt && echo "Sequence Duplication Levels" >> fastqc_categories.txt && echo "Overrepresented sequences" >> fastqc_categories.txt && echo "Adapter Content" >> fastqc_categories.txt
cp fastqc_categories.txt fastqc_reports/fastqc_report_0.txt
cat forward_reads.txt reverse_reads.txt > all_reads.txt
cd fastqc_results
unzip '*.zip'
cd ..
cat all_reads.txt | while read filename
do
current_file=$(echo ${filename:0:15}) 
echo "$current_file" > temp_fastqc_file.txt
cut -f 1 fastqc_results/"$current_file"_fastqc/summary.txt >> temp_fastqc_file.txt
paste -d'\t' fastqc_reports/fastqc_report_$counter_0.txt temp_fastqc_file.txt > fastqc_reports/fastqc_report_$counter_1.txt
rm -fr temp_fastqc_file.txt
((counter_0=counter_0+1))
((counter_1=counter_1+1))
done
finished_report=$(wc -l  all_reads.txt | cut -d ' ' -f1)
mv fastqc_reports/*$finished_report.txt fastqc_report.txt
rm -fr fastqc_reports && rm -fr all_reads && rm -fr fastqc_categories && rm -fr forward_reads && rm -fr reverse_reads && rm -fr forward_and_reverse_reads


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
#the alignment itself, with output in aligned_reads
bowtie2 -p 60 -x index_files/index_ref_file -1 /localdisk/home/data/BPSM/AY21/fastq/$input_file_1 -2 /localdisk/home/data/BPSM/AY21/fastq/$input_file_2 -S aligned_reads/$sam_name.aligned.sam 2> stats.file
echo "$sam_name" >> sam_filenames.txt
lastline=$(tail -1 stats.file)
success=$(echo ${lastline:0:1})
alignmentrate=$(echo ${lastline:0:6})
if [ "$success" -lt "8" ]
then
   echo "Replicate "$sam_name" has a rate of alignment of only "$alignmentrate". You might want to test for contamination using BLAST" >> logfile.log
fi 
done
i=0
p=1
#convert output to indexed bam
sort -t$'\t' -k1 genes.txt > genes.sorted
cat sam_filenames.txt | while read samfile
do
samtools view -b aligned_reads/$samfile.aligned.sam  > bamfiles/$samfile.aligned.bam && samtools sort bamfiles/$samfile.aligned.bam > bamfiles/$samfile.aligned.sorted.bam && samtools index bamfiles/$samfile.aligned.sorted.bam
# sam to bam and then count gene reads using bedtools coverage, which outputs text file into folder called counts
bedtools bamtobed -i bamfiles/$samfile.aligned.sorted.bam > bamfiles/$samfile.aligned.sorted.bed
bedtools coverage -counts -a /localdisk/home/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -b bamfiles/$samfile.aligned.sorted.bed > counts/$samfile.counts.txt

#Now I want to merge the gene count file, but need to sort both, as you need it before the merging
sort -t$'\t' -k4 counts/$samfile.counts.txt > counts/$samfile.counts.sorted.txt
#cut the gene id and the count columns
cut -f 4,6 counts/$samfile.counts.sorted.txt > temp_count_file.txt
sed -i "1i Gene""	$samfile" temp_count_file.txt
join -t$'\t' -e 'NaN' genes.sorted temp_count_file.txt > temp_gene_counts.txt

#I need to make sure that the files are the same length, before joining them

genelines1=$(wc -l  genes.txt | cut -d ' ' -f1)
genelines2=$(wc -l  temp_gene_counts.txt | cut -d ' ' -f1)

if [ $genelines1 != $genelines2 ]
then
   echo "Error in replicate "$samfile"" >> logfile.log
fi 

cut -f 3 temp_gene_counts.txt > temp_gene_counts_single_columns.txt

paste -d'\t' final_counts/gene_counts_$i.txt temp_gene_counts_single_columns.txt > final_counts/gene_counts_$p.txt
((i=i+1))
((p=p+1))
done
#Select the last file, which is the one containing all the gene counts
last_file=$(wc -l  sam_filenames.txt | cut -d ' ' -f1)
mv final_counts/*$last_file.txt replicate_gene_expressions.txt

mkdir mean_counts
tail -n+2 /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | awk '{print $2}'| sort | uniq > mean_counts/samples.txt
cat mean_counts/samples.txt | while read sample
do
# 0 hours uninduced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | awk -F '\t' '$4 == "0" {print $1}' > mean_counts/$sample.0hr.replicates.txt
   cat mean_counts/$sample.0hr.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.0hr.uninduced.samples
   rm -fr mean_counts/*temp
   rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.0hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.0hr.uninduced.mean
   sed -i "1i""$sample"".0hr.uninduced" mean_counts/$sample.0hr.uninduced.mean
   rm -fr mean_counts/*samples

#24 hours induced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Induced"| awk -F '\t' '$4 == "24" {print $1}' > mean_counts/$sample.24hr.induced.replicates.txt
   cat mean_counts/$sample.24hr.induced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.24hr.induced.samples
   rm -fr mean_counts/*temp
   rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.24hr.induced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.24hr.induced.mean
   sed -i "1i""$sample"".24hr.induced" mean_counts/$sample.24hr.induced.mean
   rm -fr mean_counts/*samples
#24 hours uninduced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Uninduced"| awk -F '\t' '$4 == "24" {print $1}' > mean_counts/$sample.24hr.uninduced.replicates.txt
   cat mean_counts/$sample.24hr.uninduced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.24hr.uninduced.samples
   rm -fr mean_counts/*temp
   rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.24hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.24hr.uninduced.mean
   sed -i "1i""$sample"".24hr.uninduced" mean_counts/$sample.24hr.uninduced.mean
   rm -fr mean_counts/*samples
#48 hours induced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Induced"| awk -F '\t' '$4 == "48" {print $1}' > mean_counts/$sample.48hr.induced.replicates.txt
   cat mean_counts/$sample.48hr.induced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.48hr.induced.samples
   rm -fr mean_counts/*temp
   rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.48hr.induced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.48hr.induced.mean
   sed -i "1i""$sample"".48hr.induced" mean_counts/$sample.48hr.induced.mean
   rm -fr mean_counts/*samples
#48 hours uninduced
   cat /localdisk/home/data/BPSM/AY21/fastq/100k.fqfiles | grep "$sample" | grep "Uninduced"| awk -F '\t' '$4 == "48" {print $1}' > mean_counts/$sample.48hr.uninduced.replicates.txt
   cat mean_counts/$sample.48hr.uninduced.replicates.txt | while read replicate
   do
      column=$(head -1 replicate_gene_expressions.txt | tr '\t' '\n' | cat -n | grep "$replicate"| cut -d$'\t' -f1)
      cut -f $column replicate_gene_expressions.txt > mean_counts/$replicate.temp
   done
   paste mean_counts/*temp > mean_counts/$sample.48hr.uninduced.samples
   rm -fr mean_counts/*temp
   rm -fr mean_counts/*replicates.txt
   tail -n+2 mean_counts/$sample.48hr.uninduced.samples | awk '{sum=0; for (i=1; i<=NF; i++) {sum=sum+$i;} m=sum/NF; print $0, m; }' | awk '{ print $4 }' > mean_counts/$sample.48hr.uninduced.mean
   sed -i "1i""$sample"".48hr.uninduced" mean_counts/$sample.48hr.uninduced.mean
   rm -fr mean_counts/*samples
done
paste mean_counts/*mean > mean_counts/replicate_mean_counts.txt
paste genes.txt mean_counts/replicate_mean_counts.txt > replicate_mean_counts.txt
cp mean_counts/samples.txt samples.txt

#Fold change
#doing all possible comparisons
head -1 replicate_mean_counts.txt | tr "\t" "\n" | grep -v Gene | grep -v Gene_description > all_groups.txt
cat all_groups.txt | while read i
do
  cat all_groups.txt | while read j
  do
    if [ "$i" \< "$j" ]
    then
     echo "$i:$j" >> all_combinations.txt
     echo "$j:$i" >> all_combinations.txt
    fi
  done
done
rm -fr all_groups.txt
#Now to find the fold change for every possible comparison
cat all_combinations.txt | while read STR
do
   comp1=$(echo $STR | cut -f1 -d ':')
   comp2=$(echo $STR | cut -f2 -d ':')

   column1=$(head -1 replicate_mean_counts.txt | tr '\t' '\n' | cat -n | grep $comp1 | cut -d$'\t' -f1)
   column2=$(head -1 replicate_mean_counts.txt | tr '\t' '\n' | cat -n | grep $comp2 | cut -d$'\t' -f1)
   cut -f $column1 replicate_mean_counts.txt > groupwise_comparisons/mean.counts.$comp1 && cut -f $column2 replicate_mean_counts.txt > groupwise_comparisons/mean.counts.$comp2
   paste -d'\t' groupwise_comparisons/mean.counts.$comp1 groupwise_comparisons/mean.counts.$comp2 > groupwise_comparisons/comparison_columns
   sed -i '1s/$/  fold_change/' groupwise_comparisons/comparison_columns
   awk -v OFS='\t' 'NR!=1 {$3 = ($2 != 0) ? sprintf("%.3f", $1 / $2) : "NAN"}1' groupwise_comparisons/comparison_columns > groupwise_comparisons/$comp1.vs.$comp2.single_column
   paste -d'\t' genes.sorted groupwise_comparisons/$comp1.vs.$comp2.single_column | tail -n+2 | sort -t$'\t' -n -k5 -r > groupwise_comparisons/$comp1.vs.$comp2
   sed -i "1i""Gene  Gene_description  ""$comp1""  ""$comp2""  ""fold_change" groupwise_comparisons/$comp1.vs.$comp2
   rm -fr groupwise_comparisons/mean* && rm -fr groupwise_comparisons/comparison* && rm -fr groupwise_comparisons/*single_column
done

#Remove messy files and directories
rm -fr forward_and_reverse_reads.txt && rm -fr forward_reads.txt && rm -fr reverse_reads.txt && rm -fr temp_gene_counts_single_columns.txt && rm -fr temp_gene_counts.txt && rm -fr temp_count_file.txt && rm -fr genes.txt && rm -fr genes.sorted.txt && rm -fr sam_filenames.txt && rm -fr all_combinations.txt
rm -fr counts && rm -fr final_counts && rm -fr mean_counts && rm -fr stats.file && rm -fr all_reads.txt && rm -fr fastqc_categories.txt && rm -fr samples.txt && rm -fr genes.sorted
