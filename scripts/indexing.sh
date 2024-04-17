#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=02:00:00

# indexing reference and annotation
genome=/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus.20130116.OM.fa

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/reference/reference_indexed \
    --genomeFastaFiles $genome \
    --sjdbGTFfile /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/reference/Struthio_camelus_australis.ASM69896v1.106.gtf \
    --sjdbOverhang 149

exit 0