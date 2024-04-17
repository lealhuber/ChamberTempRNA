#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=06:00:00

# run in analysis/testsamples/QCtrimming
# run fastqc on all samples, should make an .html files as report in same folder
fastqc -t 8 ../../../data/00_fastq/*.fastq.gz

# run trim galore! to trim adaptors (will detect them automatically) and reads with <Q10
trim_galore --paired --fastqc --gzip --cores 8 -q 10 \ 
    -o /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/illumina_test/trimmed/ \
    /faststorage/project/ostrich_thermal/BACKUP/RNAseqTest/00_fastq/*.fastq.gz
# -q 10 : Rachel recommended only light quality trimming of <Q10 (instead of default <Q20)
exit 0