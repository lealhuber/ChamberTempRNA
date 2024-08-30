#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64g
#SBATCH --time=60:00:00

emapper.py --cpu 4 -i /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.aa \
    --report_orthologs -o Struthio_camelus_HiC

exit 0