#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=02:00:00

# indexing reference and annotation
genome=/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/HiC_indexed \
    --genomeFastaFiles $genome \
    --sjdbGTFfile /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gtf \
    --sjdbOverhang 149

exit 0