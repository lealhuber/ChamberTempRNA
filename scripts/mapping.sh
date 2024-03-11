#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# run in analysis/testsamples

# go through each file seperately for mapping
for f in *.fq.gz ; do
    # map with STAR with an output BAM 
    STAR --runThreadN 12 --genomeDir /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/reference/reference_indexed \
        --readFilesIn $f --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix STARout/$f.
    # makes files sample#.Aligned.sortedByCoord.out.bam
    # Grab all alignments except ones that map to more than 4 locations:
    samtools view -b -q 1 STARout/$f.Aligned.sortedByCoord.out.bam > STARout/$f.Aligned.sortedByCoord.q1.out.bam
    # index because htseq is complaining
    samtools index -@ 12 STARout/$f.Aligned.sortedByCoord.q1.out.bam # generates BAI-format index Aligned.sortedByCoord.q1.out.bam.bai
    #Count the reads
    htseq-count --stranded=no STARout/$f.Aligned.sortedByCoord.q1.out.bam \
        /faststorage/project/ostrich_thermal/people/leah/broiler_testDGEA/broiler_testdata/reference/ensembl/Gallus_gallus.GRCg6a.106.gtf \
        > $f.count_matrix.txt ;
done

# make R style file with filenames for DESeqDataSetFromHTSeqCount function
ls *count_matrix.txt | sed 's/.*/"&",/' > matrixFileNames.txt
# same with first part of sample names (ADJUST!!!)
ls *count_matrix.txt | sed 's/\(B...._rep.\)\(_.*\)/"\1",/' > SampleNames.txt
# make file with just file names to copy them to local computer
ls *count_matrix.txt > matrixFileNamesCopy.txt

exit 0