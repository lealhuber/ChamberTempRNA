#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00

#run in steps/mapping

for filename in mapped/L38.Aligned.sortedByCoord.q255.out.bam; do
    ## index (if not already done)
    if ! test -f ${filename}.bai; then
        samtools index -@ 4 $filename
    fi
    ## Count the reads
    # default mode is union and noneunique mode none
    # format bam, stranded no, ordered by position (not by name)
    htseq-count -f bam -s no -r pos $filename \
        /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gff \
        > counts/L38.count_matrix.txt
    echo "htseq-count done"
done

# # make R style file with filenames for DESeqDataSetFromHTSeqCount function
# ls *count_matrix.txt | sed 's/.*/"&",/' > matrixFileNames.txt
# # same with first part of sample names (ADJUST!!!)
# ls *count_matrix.txt | sed 's/\(L.\)\(.Aligned*\)/"\1",/' > SampleNames.txt
# # make file with just file names to copy them to local computer
# ls *count_matrix.txt > matrixFileNamesCopy.txt

exit 0