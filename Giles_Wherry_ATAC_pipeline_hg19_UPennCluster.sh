#!/bin/bash

### Script can be run straight through, but I have advised certain stopping points for check the files -- some of which is important for parameter adjustments.

### TO DO BEFORE RUNNING THIS SCRIPT!!!!!
# Unzip all using gunzip *.gz

### Move files for processing
# Move bowtie index files (*bt2) into working directory
# Move hg19.chrom.size file into the working directory
# Move wgEncodeDacMapabilityConsensusExcludable.bed into working directory

### Set # of cores and memory
# Fastqc: correct -t # for number of cores
# Bowtie: Correct the number of cores "p"
# All samtools: # Change '@' for number of cores; Change 'm' for amount of memory PER CORE
# Picar tools: Set memory  '-Xmx_g'


if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load python-2.7.5
module load bowtie2-2.1.0
module load samtools-1.1
module add java-sdk-1.7.0
module load picard-tools-1.141
module load FastQC-0.11.2

# PART ONE #



### Alignment
# Mapping atac paired end files to hg19 using bowtie2
# Move bowtie index files into working directory
# Correct the number of cores "p"

for i in *_R1.fastq;
do
 sampleID="${i%_*.*}"
 echo "${sampleID}"
 bowtie2 -k1 -N1 -p16 -x mm9 -1 ${sampleID}_R1.fastq -2 ${sampleID}_R2.fastq -S ${sampleID}.sam 2> ${sampleID}.out
done

# Combine all of the .out files produced from the sam alignment

for i in *.out;
do
	sampleID="${i%.*}"
	echo "${sampleID}" > ${sampleID}_out.txt
	cat $i >> ${sampleID}_out.txt
done

find . -name '*_out.txt' | xargs cat >> combOutSams.txt
rm *_out.txt

# STOPPING POINT -- Check fastqc and alignment reports


######################################################################

# PART TWO #

######################################################################

### Make bam files
# Change '@' for number of cores
# Change 'm' for amount of memory PER CORE


# Make sorted bam
for i in `ls *.sam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	samtools view -bS -@ 16 $i | samtools sort - -m 2G -@ 16 ${sampleID}_sort
done
echo $?
echo 'You have sorted bam files :)!'


# Make flagstat file to check conversion to bam
for i in `ls *_sort.bam`;
do
	echo $i >> combFlagstat_bam.txt
	samtools flagstat $i >> combFlagstat_bam.txt
done

echo $?


# Make index for sorted bam files
for i in `ls *_sort.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'

# Make bam files with only mapped reads
for i in `ls *_sort.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	samtools view -bS -f 2 -@ 16 $i > ${sampleID}_map.bam
done
echo $?
echo 'You have mapped-only bam files!'

# Make index for mapped-only bam files
for i in `ls *_map.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'



# Make bam files without mitochondrial reads

for i in `ls *_map.bam`;
do
    sampleID="${i%*.*}"
    echo "${sampleID}"
    samtools view -b -@ 16 $i chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${sampleID}_deMito.bam
done
echo $?
echo 'You have bams without mitocondrial reads!'


# Make index for deMito bams
for i in `ls *_deMito.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'


######################################################################


### Remove duplicates

# Change 'Xmx_g for amount of available ram

for i in `ls *_deMito.bam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	java -Xmx32g -jar $PICARD_JAR MarkDuplicates INPUT=$i OUTPUT=${sampleID}_dedup.bam METRICS_FILE=${sampleID}_dedupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
done

echo $?



# Combine picard duplicate metrics
for i in `ls *_dedupMetrics.txt`;
do
	sampleID="${i%_*_*_*.*}"
	awk 'FNR==8 {print FILENAME, $0}' ${sampleID}_map_deMito_dedupMetrics.txt >> comb_picardDedupStats1.txt
done

awk '{print $1"\t"$9"\t"$10"\t"$11}' comb_picardDedupStats1.txt > comb_picardDedupStats2.txt


echo $'sampleID\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE' | cat - comb_picardDedupStats2.txt > comb_picardDedupStats.txt
rm comb_picardDedupStats1.txt
rm comb_picardDedupStats2.txt




# Make index for deduped bam files
for i in `ls *_dedup.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'


######################################################################

### Get insert metrics -- fragment size information

# Get insert metrics -- file and histogram
for i in `ls *_dedup.bam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	java -Xmx32g -jar $PICARD_JAR CollectInsertSizeMetrics I=$i OUTPUT=${sampleID}_insertMetrics H=${sampleID}_histo.pdf
done

echo $?

#####################################################################

### Convert bams to beds --> shift reads to reflect center of Tn5 insertion
#(The Tn5 transposon binds as a dimer and inserts two adaptors separated by 9 bp, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the – strand were offset −5 bp.)


for i in `ls *_dedup.bam`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools bamtobed -i  $i | awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 5, $3 + 5, $4, $5, $6; else print $1, $2 - 4, $3 - 4, $4, $5, $6}' > ${sampleID}.bed
done

echo $?
echo 'You have bed files!'


#####################################################################


### Remove blacklisted regions

# Copy mm10.blacklist.bed into working directory

for i in `ls *.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools subtract -a $i -b mm9.blacklist.bed > ${sampleID}_debl.bed
done

echo $?
echo 'You have bed files without blacklist regions!'

rm mm9.blacklist_debl.bed



#####################################################################


### Make bam files and associated indices from deblacklisted bed files

for i in `ls *_debl.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools bedtobam -i $i -g mm9.chrom.sizes | samtools sort - -m 2G -@ 16 ${sampleID}
done

echo $?
echo 'You have deblackedlisted bam files!'

# Make index for deduped bam files
for i in `ls *_debl.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for deduplicated-deblacklist bam files!'


######################################################################


### Make normalized bigwig files

# Must copy hg19.chrom.size file into working directory

for i in `ls *_debl.bam`;
do
 sampleID="${i%.*}"
 lines=$(samtools view -c ${sampleID}.bam);\
 bedtools genomecov -ibam ${sampleID}.bam -bg -scale $(echo "1000000 / ${lines} " | bc -l) -g mm9.chrom.sizes | \
 wigToBigWig -clip stdin mm9.chrom.sizes ${sampleID}.bw 2> ${sampleID}_log
done

echo $?



#######################################################################

### Get map and chromosome numbers with sorted bam and deduped_debl_bams


### With sort_bam

## Create file with mapping numbers

# Create file with row names for mapStats
echo sampleID >> rowNames_mapStats.txt
echo totalReads >> rowNames_mapStats.txt
echo totalMappedReads >> rowNames_mapStats.txt
echo unmappedReads >> rowNames_mapStats.txt
echo pairedMappedReads >> rowNames_mapStats.txt

# Get total reads, total mapped reads, total unmapped reads, properly paired reads -- for each sample
for i in *_sort.bam;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" >> ${sampleID}_mapStats.txt
	samtools view $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -F4 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -f4 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -f2 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt

done

# Combine individual mapStat files into one file
find . -name '*_mapStats.txt' | sort | xargs paste | column -s $'\t' -t > combMapStats.txt




### Create file with chromosome stats

# Create file with row names for chrStats# Get number of reads per chromosome

echo sampleID >> rowNames_chrStats.txt
echo chr1 >> rowNames_chrStats.txt
echo chr2 >> rowNames_chrStats.txt
echo chr3 >> rowNames_chrStats.txt
echo chr4 >> rowNames_chrStats.txt
echo chr5 >> rowNames_chrStats.txt
echo chr6 >> rowNames_chrStats.txt
echo chr7 >> rowNames_chrStats.txt
echo chr8 >> rowNames_chrStats.txt
echo chr9 >> rowNames_chrStats.txt
echo chr10 >> rowNames_chrStats.txt
echo chr11 >> rowNames_chrStats.txt
echo chr12 >> rowNames_chrStats.txt
echo chr13 >> rowNames_chrStats.txt
echo chr14 >> rowNames_chrStats.txt
echo chr15 >> rowNames_chrStats.txt
echo chr16 >> rowNames_chrStats.txt
echo chr17 >> rowNames_chrStats.txt
echo chr18 >> rowNames_chrStats.txt
echo chr19 >> rowNames_chrStats.txt
echo chr20 >> rowNames_chrStats.txt
echo chr21 >> rowNames_chrStats.txt
echo chr22 >> rowNames_chrStats.txt
echo chrX >> rowNames_chrStats.txt
echo chrY >> rowNames_chrStats.txt
echo chrM >> rowNames_chrStats.txt

# Get number of reads per chromosome -- for each samples
for i in *_sort.bam;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}" >> ${sampleID}_chrStats.txt
	samtools view $i chr1 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr2 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr3 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr4 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr5 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr6 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr7 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr8 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr9 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr10 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr11 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr12 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr13 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr14 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr15 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr16 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr17 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr18 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr19 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr20 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr21 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr22 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrX -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrY -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrM -@ 16 | wc -l >> ${sampleID}_chrStats.txt
done

# Combine individual chrStat files into one file
find . -name '*_chrStats.txt' | sort | xargs paste > combChrStats.txt




#####################################################

### With fully processed bam file

echo sampleID >> rowNames_finalMapStats.txt
echo totalReads >> rowNames_finalMapStats.txt


for i in *_map_deMito_dedup_debl.bam;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" >> ${sampleID}_finalMapStats.txt
	samtools view $i | wc -l >> ${sampleID}_finalMapStats.txt
done

# Combine individual mapStat files into one file
find . -name '*_finalMapStats.txt' | sort | xargs paste | column -s $'\t' -t > combFinalMapStats.txt




# STOPPING POINT -- Check all raw metrics and bigwigs generated above -- samples look ok? May want to remove bad samples before proceeding.

######################################################################





# PART THREE #

######################################################################

### Call peaks
## This should be optimized for your samples!!!!!!!!!!!!!!
## Put peaks into igv with bigwig files -- do you like the peak calling??????
## Can adjust the q value
## Can add -broad option and alter --broad cutoff (look at documentation)

# Activate python virtual environment

source /home/jogiles/my_python-2.7.5/bin/activate

for i in `ls *_debl.bed`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	macs2 callpeak -t $i -f BEDPE -g hs -n ${sampleID} -q 0.01
done

echo $?
echo 'You have called peaks!'


# Make a file with the chr, start, end, AND AN EXTRA 4TH COLUMN -- necessary for downstream -- for each sample called peak file
for i in `ls *.narrowPeak`;
do
 sampleID="${i%.*}"
 awk '{print $1"\t"$2"\t"$3"\t"$9}' ${sampleID}.narrowPeak > ${sampleID}.bg
done
echo $?



# Make a file with the chr, start, end (standard bed) -- useful for looking at peak regions in igv
for i in `ls *_map_deMito_dedup_debl_peaks.bg`;
do
 sampleID="${i%_*_*_*_*_*.*}"
 echo $sampleID
 awk '{print $1"\t"$2"\t"$3}' ${sampleID}_map_deMito_dedup_debl_peaks.bg > ${sampleID}_peaks.bed
done
echo $?

# STOPPING POINT -- Check *_peaks.bed files with bigwigs to determine best peak calling paramters.

######################################################################




# PART FOUR #

######################################################################


### Make INDIVIDUAL peak count files

### Make individual _count.bedgraph files with the 4th column being # of tags in each peak
for i in `ls *_debl.bed`;
do
 sampleID="${i%.*}"
 bedtools coverage -counts -a $i -b ${sampleID}_peaks.bg > ${sampleID}_indCounts.bedgraph
done

echo $?

# Make edited file with only 4 columns
for i in `ls *_indCounts.bedgraph`;
do
	sampleID="${i%.*}"
	awk 'BEGIN {OFS="\t"} $1 ~ /chr./ {print $1, $2, $3, $5}' $i > ${sampleID}F.bedgraph
done

echo $?
echo 'The individual coverage files are ready!'


# Remove unformatted file
# Rename formatted file
rm *_indCounts.bedgraph

for i in *_map_deMito_dedup_debl_indCountsF.bedgraph
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	mv $i ${sampleID}_indCounts.bedgraph
done


# Make file with the total number of reads in the individually called peaks --> use to calculate FRiPs

for i in `ls *_indCounts.bedgraph`;
do
 sampleID="${i%_*_*_*_*_*.*}"
 awk -v sampleID="$sampleID" 'BEGIN {OFS="\t"} {s+=$4}END{print sampleID,s}' $i >> combTotalIndCounts.txt
done

# Get number of peaks called for each sample
for i in *.bg;
do
	sampleID="${i%_*_*_*_*_*.*}"
	echo "${sampleID}" >> ${sampleID}_indPeakNumber.txt
	wc -l $i | awk '{print $1}' >> ${sampleID}_indPeakNumber.txt
done

find . -name '*_indPeakNumber.txt' | sort | xargs paste | column -s $'\t' -t > combIndPeakNumber.txt

rm *.bg_indPeakNumber.txt



########################################################################






# PART FIVE #

# This section will need to be repeated for every group comparison -- it is basically creating a consensus list of peak regions (like genes or transcripts for RNA-seq).

########################################################################


### Make UNION peak count files


## Make union of all peaks by combining peak call bed files
bedtools unionbedg -g /home/jogiles/bin/bedtools2/genomes/mouse.mm9.genome -i *.bg > unionbed.text
awk 'BEGIN {OFS="\t"} $1 ~ /chr./ {print $1, $2, $3}' unionbed.text > myunionbed.bed
echo $?

## Merge overlapping peaks
mergeBed -i myunionbed.bed > unionbed_merged.bed
echo $?

## Make  counts bedgraph files with the 4th column being # of tags in each union peak:
for i in `ls *_debl.bed`;
do
 sampleID="${i%.*}"
 bedtools coverage -counts -a $i -b unionbed_merged.bed | bedtools sort -i - > ${sampleID}_counts.bedgraph
done

echo $?
echo 'The union coverage files are ready!'

## Make file with the total number of reads in the union peak list
for i in `ls *_counts.bedgraph`;
do
 echo $i
 sampleID="${i%.*}"
 awk '{s+=$4}END{print FILENAME,s}' ${sampleID}.bedgraph >> combTotalCountUnionPeaks.txt
done

echo $?
echo 'File with total tags in union_merged list is ready!'


# Combine all count bedgraph files into one file with chr, start, end AND headers
for i in *_counts.bedgraph;
do
	sampleID="${i%_*_*_*_*_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" > ${i}_unionTags.txt
	awk '{print $4}' $i >> ${i}_unionTags.txt
done

find . -name '*_unionTags.txt' | sort | xargs paste | column -s $'\t' -t > combUnionReads.txt

echo $'chr\tstart\tend' | cat - unionbed_merged.bed > unionbed_merged2.bed

paste unionbed_merged2.bed combUnionReads.txt > combUnionReadsWithLabels.txt

# Get number of peaks in union&merged list
wc -l unionbed_merged.bed | cut -d' ' -f1 >> combUnionPeakNumber.txt

rm *_unionTags.txt

######################################################################


########################################################################

### Make folders and sort files

mkdir output

mv fastqc output
mkdir output/fastqs
mkdir output/sams
mkdir output/sortBams
mkdir output/mapBams
mkdir output/dedupBams
mkdir output/deblBams
mkdir output/deMitoBams
mkdir output/dedupBeds
mkdir output/deblBeds
mkdir output/peaks
mkdir output/unionBedgraphs
mkdir output/indBedgraphs
mkdir output/bws
mkdir output/histos
mkdir output/combFiles
mkdir output/indStats
mkdir output/insertMetrics
mkdir output/dedupMetrics

mv *fastq output/fastqs
mv *sam output/sams
mv *out output/sams
mv *sort.bam *sort.bam.bai output/sortBams
mv *map.bam *map.bam.bai output/mapBams
mv *dedup.bam *dedup.bam.bai output/dedupBams
mv *map_deMito.bam *map_deMito.bam.bai output/deMitoBams
mv *map_deMito_dedup.bed output/dedupBeds
mv *debl.bed output/deblBeds
mv *debl.bam *debl.bam.bai output/deblBams
mv *_peaks.xls *_peaks.narrowPeak *_summits.bed *_peaks.bg *_peaks.bed output/peaks
mv *_counts.bedgraph output/unionBedgraphs
mv *indCounts.bedgraph output/indBedgraphs
mv *.bw output/bws
mv *_log output/bws

mv *_mapStats.txt output/indStats
mv *_chrStats.txt output/indStats
mv *_indPeakNumber.txt output/indStats

mv *map_deMito_dedupMetrics.txt output/dedupMetrics
mv *map_deMito_dedup_insertMetrics output/insertMetrics

mv *.pdf output/histos
mv output/histos output/combFiles
mv comb* *union* output/combFiles
