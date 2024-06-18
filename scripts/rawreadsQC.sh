#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=06:00:00

# rename them, ie remove the _001
cd /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/90-964973806/00_fastq
for filename in *.fastq.gz; do
    base=${filename%_001.fastq.gz}
    mv $filename $base.fastq.gz
done

# run in steps/QC
cd /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/steps/QC

# run fastqc on all samples, should make an .html files as report in same folder
fastqc -t 8 /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/90-964973806/00_fastq/*.fastq.gz

# run trim galore! to trim adaptors (will detect them automatically) and reads with <Q10 (good enough for DGE)
trim_galore --paired --fastqc --gzip --cores 8 -q 10 \
    -o /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/trimmed \
    /faststorage/project/ostrich_thermal/BACKUP/ChamberRNAseq/data/90-964973806/00_fastq/*.fastq.gz
# -q 10 : Rachel recommended only light quality trimming of <Q10 (instead of default <Q20)
exit 0