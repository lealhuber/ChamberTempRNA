#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=72:00:00

# run in steps/mapping and have a symlink to data there? But then there's a path in f...
# create file with base-filenames without read1 or read2 called filebase and run in same directory
# adjust for paired-end data!!! Test this on 2 samples first or so

# I renamed them all in a previous attempt that otherwise failed with this code:
# for filename in $dat_dir/*.fq.gz; do
#     # rename them, ie remove the _val part
#     baseread=${filename%_val_*.fq.gz}
#     mv $filename $baseread.trimmed.fq.gz

## go through each file seperately for mapping
# define path to data directory
dat_dir=/faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/trimmed

# want to do everything only once of course!! so only take the R1 file
for filename in $dat_dir/*_R1.trimmed.fq.gz; do
    # extract just the sample filebase name without read information
    namenopath=${filename##*/}
    filebase=${namenopath%_*}
    echo $filebase
    # assemble and define filenames for each read pair
    fq_read1=$dat_dir/${filebase}_R1.trimmed.fq.gz
    fq_read2=$dat_dir/${filebase}_R2.trimmed.fq.gz
    # I want to do this only
    # map with STAR with an output BAM 
    STAR --runThreadN 10 --genomeDir /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/HiC_indexed \
        --readFilesIn $fq_read1 $fq_read2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $filebase.
    # makes files sample#.Aligned.sortedByCoord.out.bam
    echo "Star done"
    # Grab all alignments except ones that map to more than 4 locations:
    samtools view -b -q 1 $filebase.Aligned.sortedByCoord.out.bam > $filebase.Aligned.sortedByCoord.q1.out.bam
    # index because htseq is complaining
    samtools index -@ 10 $filebase.Aligned.sortedByCoord.q1.out.bam # generates BAI-format index Aligned.sortedByCoord.q1.out.bam.bai

    ## this didn't work / was added later
    # Check coverage, try to infer about library creation and sequencing
    samtools stats -@ 10 -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta $filebase.Aligned.sortedByCoord.out.bam > $filebase.Aligned.sortedByCoord.out.bam.stats
    plot-bamstats -p plots/ $filebase.Aligned.sortedByCoord.out.bam.stats
done

exit 0