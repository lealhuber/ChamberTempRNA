#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=18:00:00

# run in mapping
# QC of mapped samples
for filename in L*.sortedByCoord.out.bam; do
    # samtools stats stuff
    filebase=${filename%.Aligned*}
    plot-bamstats -p plots/${filebase}/ $filebase.Aligned.sortedByCoord.out.bam.stats
    # Calculate the distribution of mismatches across reads. Doesn't work because my BAM files do not have an MD tag :(
    # mismatch_profile.py -i $filename -l 150 -o ../QC/rseqc/${filename}

    # Calculate the distributions of inserted nucleotides across reads
    insertion_profile.py -i $filename -o .../QC/rseqc/${filename} -s "PE"

    # read duplication sequence vs. mapping based
    read_duplication.py -i $filename -o ../QC/rseqc/${filename}

    # need BED format for this !!
    # Calculate inner distance between read pairs.
    inner_distance.py -i $filename -o ../QC/rseqc/${filename} -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.bed

    # this program will compare detected splice junctions to reference gene model
    # splicing annotation is performed in two levels: splice event level and splice junction level.
    junction_annotation.py -i $filename -o ../QC/rseqc/${filename} -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.bed

    # read distribution over genome features like exons, introns etc.
    read_distribution.py -i $filename -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.bed > ../QC/rseqc/${filename}_dist.txt

    FPKM_count.py -i $filename -o ../QC/rseqc/${filename} -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.bed
done