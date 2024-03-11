#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=01:00:00

# run in analysis/intermediate

# run fastqc on all samples, should make an .html files as report in same folder
fastqc -t 6 ../../data/illumina_test/*.fq.gz

# run trim galore! to trim adaptors (will detect them automatically) and reads with <Q10
trim_galore --paired --fastqc --gzip --cores 8 -q 10 \
    -o /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/illumina_test/trimmed \
    ../../data/illumina_test/*.fq.gz

exit 0