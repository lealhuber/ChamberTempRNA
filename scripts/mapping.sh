#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# run in analysis/testsamples/mapping and have a symlink to data there? But then there's a path in f...
# create file with base-filenames without read1 or read2 called filebase and run in same directory
# adjust for paired-end data!!! Test this on 2 samples first or so

## go through each file seperately for mapping
# define path to data directory
dat_dir=/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/illumina_test/trimmed

while read filebase
do
    echo $filebase
    # assemble and define filenames
    fq_read1=$dat_dir/${filebase}_R1.fq.gz
    fq_read2=$dat_dir/${filebase}_R2.fq.gz
    # map with STAR with an output BAM 
    STAR --runThreadN 10 --genomeDir /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/reference/reference_indexed \
        --readFilesIn $fq_read1 $fq_read2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $filebase.
    # makes files sample#.Aligned.sortedByCoord.out.bam
    # Grab all alignments except ones that map to more than 4 locations:
    samtools view -b -q 1 $filebase.Aligned.sortedByCoord.out.bam > $filebase.Aligned.sortedByCoord.q1.out.bam
    # index because htseq is complaining
    samtools index -@ 10 $filebase.Aligned.sortedByCoord.q1.out.bam # generates BAI-format index Aligned.sortedByCoord.q1.out.bam.bai
    #Count the reads
    htseq-count --stranded=no $filebase.Aligned.sortedByCoord.q1.out.bam \
        /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/reference/Struthio_camelus_australis.ASM69896v1.106.gtf \
        > $filebase.count_matrix.txt ;
done < filebase

# make R style file with filenames for DESeqDataSetFromHTSeqCount function
ls *count_matrix.txt | sed 's/.*/"&",/' > matrixFileNames.txt
# same with first part of sample names (ADJUST!!!)
ls *count_matrix.txt | sed 's/\(B...._rep.\)\(_.*\)/"\1",/' > SampleNames.txt
# make file with just file names to copy them to local computer
ls *count_matrix.txt > matrixFileNamesCopy.txt

exit 0