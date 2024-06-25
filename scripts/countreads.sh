#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00

#run in steps/mapping

for filename in *q1.out.bam; do
    #Count the reads
    htseq-count --stranded=no $filename \
        /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gff \
        > $filename.count_matrix.txt
    echo $filename
    echo "htseq-count done"
done

# make R style file with filenames for DESeqDataSetFromHTSeqCount function
ls *count_matrix.txt | sed 's/.*/"&",/' > matrixFileNames.txt
# same with first part of sample names (ADJUST!!!)
ls *count_matrix.txt | sed 's/\(L.\)\(.Aligned*\)/"\1",/' > SampleNames.txt
# make file with just file names to copy them to local computer
ls *count_matrix.txt > matrixFileNamesCopy.txt

exit 0