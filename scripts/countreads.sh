#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=10:00:00

#run in steps/**/mapping

samtools idxstats L38.Aligned.sortedByCoord.q255.out.bam | awk '{print $1, $2}' > L38.chr_sizes.txt
awk '$2 > 10000000 {print $1}' L38.chr_sizes.txt > L38.large_chromosomes.txt
awk '$2 <= 10000000 {print $1}' L38.chr_sizes.txt > L38.small_chromosomes.txt
while read chr; do
    samtools view -b L38.Aligned.sortedByCoord.q255.out.bam $chr > L38.chunk_$chr.bam
done < L38.large_chromosomes.txt
# Create a single BAM file for all small chromosomes
samtools view -b L38.Aligned.sortedByCoord.q255.out.bam $(cat L38.small_chromosomes.txt | tr '\n' ' ') > L38.chunk_small.bam

# sample L38 keeps running out of memory, so I am removing multiple mapped reads and duplicate reads
while read chr; do
    samtools sort -n -@ 4 -o L38.chunk_${chr}.sortedByName.bam L38.chunk_${chr}.bam &&
    samtools fixmate -r -m -@ 4 L38.chunk_${chr}.sortedByName.bam - |
        samtools sort -@ 4 -o L38.chunk_${chr}.sortedByCoord.bam &&
    samtools markdup -r -s --duplicate-count -@ 4 L38.${chr}.sortedByCoord.bam L38.chunk_${chr}.sortedByCoord.dedup.bam &&
    samtools index -@ 4 L38.chunk_${chr}.sortedByCoord.dedup.bam
    samtools stats -@ 4 -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta \
        L38.chunk_${chr}.sortedByCoord.dedup.bam > L38.chunk_${chr}.sortedByCoord.dedup.bam.stats
done < L38.large_chromosomes.txt

samtools sort -n -@ 4 -o L38.chunk_small.sortedByName.bam L38.chunk_small.bam &&
samtools fixmate -r -m -@ 4 L38.chunk_small.sortedByName.bam - |
    samtools sort -@ 4 -o L38.chunk_small.sortedByCoord.bam &&
samtools markdup -r -s --duplicate-count -@ 4 L38.chunk_small.sortedByCoord.bam L38.chunk_small.sortedByCoord.dedup.bam
samtools index -@ 4 L38.chunk_small.sortedByCoord.dedup.bam
samtools stats -@ 4 -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta \
    L38.chunk_small.sortedByCoord.dedup.bam > L38.chunk_small.sortedByCoord.dedup.bam.stats

samtools merge -@ 4 L38.sortedByCoord.dedup.q255.out.bam L38.chunk_*.sortedByCoord.dedup.bam


for filename in L38.sortedByCoord.dedup.q255.out.bam; do
    ## index (if not already done)
    if ! test -f ${filename}.bai; then
        samtools index -@ 4 $filename
    fi
    ## Count the reads
    # default mode is union and noneunique mode none
    # format bam, stranded no, ordered by position (not by name)
    htseq-count -f bam -s no -r name $filename \
        /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gff \
        > ../counts/L38.count_matrix.txt
    echo "htseq-count done"
done

# make R style file with filenames for DESeqDataSetFromHTSeqCount function
ls *dedup.count_matrix.txt | sed 's/.*/"&",/' > matrixFileNames.txt
# same with first part of sample names (ADJUST!!!)
ls *dedup.count_matrix.txt | sed 's/\(L.\)\(.Aligned*\)/"\1",/' > SampleNames.txt
# make file with just file names to copy them to local computer
ls *dedup.count_matrix.txt > matrixFileNamesCopy.txt

exit 0