## declare the workflow
from gwf import Workflow # type: ignore
import numpy as np # type: ignore

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

#############################
### Process all fastq files 
### with fastqc 
#############################


## Extract all file names from the folder data which end with .fastq.gz
#  and make a list with all the names
import os # module provides functions for interacting with operating system, eg with files
import glob # module for regex patterns for finding file patterns or names

pathname = '/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/Resequenced'
out_dir = 'Reseq'
 #pattern for path, join inserts /
pattern = os.path.join(pathname, '*.fastq.gz')
fastq_files = glob.glob(pattern) #list with names and their path
# print("Fastq files in the data folder")
# print(fastq_files)
fastq_files = [os.path.basename(file) for file in fastq_files] #remove the path from the file
# print("Fastq files after removing path from the name")
# print(fastq_files) #print names


#############################
### Trimming step with Trim Galore
### for each pair of fastq files in parallel
#############################

fastq_basenames = [ FASTQ.split("_0")[0] for FASTQ in fastq_files ] #names without extensions
#print("Just the fastq basenames")
#print(fastq_basenames)

fastq_pairs = np.unique([ FASTQ.split('_')[0] for FASTQ in fastq_basenames ]) #names without pair number and extension
# print(fastq_pairs)

for FASTQ_PAIR in fastq_pairs:
    gwf.target( f"TrimGalore_{FASTQ_PAIR}", #name of the target
           cores=4,
           memory='8gb',
           walltime='06:00:00',	
           inputs= {'FASTQS': [f'{pathname}/{FASTQ_PAIR}_R1_001.fastq.gz',
                               f'{pathname}/{FASTQ_PAIR}_R2_001.fastq.gz']}, 
           outputs=[f'{out_dir}/trimmed/{FASTQ_PAIR}_val_1.fq.gz', 
                    f'{out_dir}/trimmed/{FASTQ_PAIR}_val_2.fq.gz',
                    f'{out_dir}/trimmed/{FASTQ_PAIR}_R1_001.fastq.gz_trimming_report.txt',
                    f'{out_dir}/trimmed/{FASTQ_PAIR}_R2_001.fastq.gz_trimming_report.txt',
                    f'{out_dir}/trimmed/{FASTQ_PAIR}_val_1_fastqc.html',
                    f'{out_dir}/trimmed/{FASTQ_PAIR}_val_2_fastqc.html']) << """
    mkdir -p {out_dir}/trimmed
    trim_galore --paired --fastqc --gzip --length 80 --cores 4 -q 10 --basename {FASTQ_PAIR} \
        -o {out_dir}/trimmed \
        {pathname}/{FASTQ_PAIR}_R1_001.fastq.gz \
        {pathname}/{FASTQ_PAIR}_R2_001.fastq.gz
    """.format(FASTQ_PAIR=FASTQ_PAIR, pathname=pathname, out_dir=out_dir) # type: ignore


#############################
### Mapping in parallel
#############################

for FASTQ_PAIR in fastq_pairs:
    # print(FASTQ_PAIR)
    gwf.target( f"STAR_{FASTQ_PAIR}",
        cores=10,
        memory='32G',
        walltime='72:00:00',
        inputs= {'FASTQ_TRIM': [f'{out_dir}/trimmed/{FASTQ_PAIR}_val_1.fq.gz',
                                f'{out_dir}/trimmed/{FASTQ_PAIR}_val_2.fq.gz']},
        outputs=[f'{out_dir}/mapped/{FASTQ_PAIR}.Aligned.sortedByCoord.out.bam']) << """
    mkdir -p {out_dir}/mapped
    STAR --runThreadN 10 --genomeDir /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/HiC_indexed \
        --readFilesIn {out_dir}/trimmed/{FASTQ_PAIR}_val_1.fq.gz {out_dir}/trimmed/{FASTQ_PAIR}_val_2.fq.gz \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {out_dir}/mapped/{FASTQ_PAIR}.
    """.format(FASTQ_PAIR=FASTQ_PAIR, out_dir=out_dir)  # type: ignore

# STAR mapps the reads and outputs a BAM file that's sorted by coordinates (otherwise all standard settings)


#############################
### Removing duplicates
#############################

for FASTQ_PAIR in fastq_pairs:
    # print(FASTQ_PAIR)
    gwf.target( f"Dedup_{FASTQ_PAIR}",
        cores=4,
        memory='64G',
        walltime='72:00:00',
        inputs= {'BAMs': [f'{out_dir}/mapped/{FASTQ_PAIR}.Aligned.sortedByCoord.out.bam']},
        outputs=[f'{out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam',
                 f'{out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam.bai',
                 f'{out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam.stats'],
        protect=[f'{out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam.stats']) << """
    samtools view -b -q 255 {out_dir}/mapped/{FASTQ_PAIR}.Aligned.sortedByCoord.out.bam |
        samtools sort -n -@ 4 -o {out_dir}/mapped/{FASTQ_PAIR}.sortedByName.bam &&
    samtools fixmate -r -m -@ 4 {out_dir}/mapped/{FASTQ_PAIR}.sortedByName.bam - |
        samtools sort -@ 4 -o {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.bam &&
    samtools markdup -r -s --duplicate-count -@ 4 {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.bam {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam
    samtools index -@ 4 {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam
    samtools stats -@ 4 -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta \
        {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam > {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam.stats
    """.format(FASTQ_PAIR=FASTQ_PAIR, out_dir=out_dir)  # type: ignore

# samtools view grabs all except ones that map (I doublechecked this with STAR manual))
# --> next time grab only reads that mapped to 1 location to increase speed of htseq-count
# using this (curtesy of ChatGPT): samtools view -h {FASTQ_PAIR}.Aligned.sortedByCoord.out.bam | awk '$1 ~ /^@/ || $5 == 255' | samtools view -b > {FASTQ_PAIR}.uniquely_mapped.bam
# this prints bam including header, awk grabs all headers (start with @) and reads with MAPQ=255, than converts back to bam
# or maybe test first whether simply -q 255 works 
# sorting by name necessary for fixmate
# fixmate and sorting by coordinate necessary for markdup
# samtools indexing is needed by htseq
# samtools stats makes some stats

#############################
### Count reads
#############################

for FASTQ_PAIR in fastq_pairs:
    # print(FASTQ_PAIR)
    gwf.target( f"Count_{FASTQ_PAIR}",
        cores=1,
        memory='64gb',
        walltime='10:00:00',
        inputs={'BAMS': [f'{out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam']},
        outputs=[f'{out_dir}/counts/{FASTQ_PAIR}.dedup.count_matrix.txt'],
        protect=[f'{out_dir}/counts/{FASTQ_PAIR}.dedup.count_matrix.txt']) << """
    mkdir -p {out_dir}/counts
    htseq-count -f bam -s no -r pos {out_dir}/mapped/{FASTQ_PAIR}.sortedByCoord.dedup.q255.out.bam \
        /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gff \
        > {out_dir}/counts/{FASTQ_PAIR}.dedup.count_matrix.txt
""".format(FASTQ_PAIR=FASTQ_PAIR, out_dir=out_dir) # type: ignore

# counts genes with default mode union and noneunique mode none
# format bam, not stranded, and ordered by position (not name)

#################################
### Make list of files for DESeq2
#################################

# make R style file with filenames for DESeqDataSetFromHTSeqCount function
# same with first part of sample names (ADJUST!!!)
# make file with just file names to copy them to local computer

# Input file is not perfect yet!!!

# gwf.target("FileList",
#            cores = 1,
#            memory='1gb',
#            walltime='00:05:00',
#            inputs=[f'{out_dir}/counts/L66.count_matrix.txt'],
#            outputs=[f'{out_dir}/matrixFileNames.txt',
#                     f'{out_dir}/SampleNames.txt',
#                     f'{out_dir}/matrixFileNamesCopy.txt']) << """
#     ls {out_dir}/counts/ | sed 's/.*/"&",/' > {out_dir}/matrixFileNames.txt
#     ls {out_dir}/counts/ | sed 's/\\(L.\\)\\(.Aligned*\\)/"\\1",/' > {out_dir}/SampleNames.txt
#     ls {out_dir}/counts/ > {out_dir}/matrixFileNamesCopy.txt
# """.format(out_dir=out_dir) # type: ignore

