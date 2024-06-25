#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=36:00:00

# run trim galore! to trim adaptors (will detect them automatically) and reads with <Q10 (good enough for DGE)
trim_galore --paired --fastqc --gzip --cores 8 -q 10 \
    -o /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/trimmed \
    /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/90-964973806/00_fastq/L*.fastq.gz
# -q 10 : Rachel recommended only light quality trimming of <Q10 (instead of default <Q20)
exit 0
