## declare the workflow
from gwf import Workflow
import numpy as np

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
out_dir = '/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/steps/QC_Resequenced'
 #pattern for path, join inserts /
pattern = os.path.join(pathname, '*.fastq.gz')
fastq_files = glob.glob(pattern) #list with names and their path
#print("Fastq files in the data folder")
#print(fastq_files)
fastq_files = [os.path.basename(file) for file in fastq_files] #remove the path from the file
#print("Fastq files after removing path from the name")
#print(fastq_files) #print names

for FASTQ in fastq_files:
    ## A target doing FastQC on a file
    FASTQ_BASENAME = FASTQ.split("_")[0,1] #name without extensions and without the _001
    gwf.target( f'FastQC_{FASTQ_BASENAME}', # name of the target
           cores=1,
           memory='8gb',
           walltime='00:05:00',	
           inputs=[f'{pathname}/{FASTQ}'], 
           outputs=[f'{out_dir}/{FASTQ_BASENAME}_fastqc.html', 
                    f'.{out_dir}/{FASTQ_BASENAME}_fastqc.zip']) << """ 
    eval “$(conda shell.bash hook)”
    conda activate QC                
    fastqc {pathname}/{FASTQ}
    """.format(FASTQ=FASTQ, pathname=pathname, outdir=outdir)


#############################
### Trimming step with Trim Galore
### for each pair of fastq files in parallel
#############################

fastq_basenames = [ FASTQ.split("_")[0,1] for FASTQ in fastq_files ] #names without extensions
#print("Just the fastq basenames")
#print(fastq_basenames)

fastq_pairs = np.unique([ FASTQ.split('_')[0] for FASTQ in fastq_basenames ]) #names without pair number and extension


for FASTQ_PAIR in fastq_pairs:
    gwf.target( f"TrimGalore_{FASTQ_PAIR}", #name of the target
           cores=4,
           memory='8gb',
           walltime='00:10:00',	
           inputs= {'FASTQS': [f'{pathname}/{FASTQ_PAIR}_R1_001.fastq.gz', # OBS! most probably needs adjustment!
                               f'{pathname}/{FASTQ_PAIR}_R2_001.fastq.gz'],
                    'MULTIQC': ['results/multiqc_output/multiqc_report.html']}, 
           outputs=[f'{pathname}/trimmed/{FASTQ_PAIR}_1_trimmed.fq.gz', 
                    f'{pathname}/trimmed/{FASTQ_PAIR}_2_trimmed.fq.gz',
                    f'{out_dir}/{FASTQ_PAIR}_1.fastq.gz_trimming_report.txt',
                    f'{out_dir}/{FASTQ_PAIR}_1_val_1_fastqc.html',
                    f'{out_dir}/{FASTQ_PAIR}_2.fastq.gz_trimming_report.txt',
                    f'{out_dir}/{FASTQ_PAIR}_2_val_2_fastqc.html']) << """
    eval “$(conda shell.bash hook)”
    conda activate broiler_test
    mkdir -p {pathname}/trimmed
    trim_galore --paired --fastqc --gzip --length 80 --cores 4 -q 10 \
        -o {pathname}/trimmed \
        {pathname}/{FASTQ_PAIR}_R1_001.fastq.gz \
        {pathname}/{FASTQ_PAIR}_R2_001.fastq.gz
    """.format(FASTQ_PAIR=FASTQ_PAIR, pathname='/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/Resequenced')


#############################
### Mapping in parallel
#############################

for FASTQ in fastq_files:
    FASTQ_BASENAME = FASTQ.split("_")[0] # just sample name now
    gwf.target( f"STAR_{FASTQ_BASENAME}",
        cores=10
        memory='32G',
        walltime='72:00:00',
        inputs= {'FASTQ_TRIM': [f'{pathname}/trimmed/{FASTQ_BASENAME}_R1.trimmed.fq.gz',
                                f'{pathname}/trimmed/{FASTQ_BASENAME}_R2.trimmed.fq.gz']},
        outputs=[f'{pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam',
                 f'{pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam.bai',
                 f'{pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam.stats']) <<< """
    eval “$(conda shell.bash hook)”
    conda activate broiler_test
    mkdir -p {pathname}/mapped
    STAR --runThreadN 10 --genomeDir /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/HiC_indexed \
        --readFilesIn {pathname}/trimmed/{FASTQ_BASENAME}_R1.trimmed.fq.gz {pathname}/trimmed/{FASTQ_BASENAME}_R2.trimmed.fq.gz \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {FASTQ_BASENAME}.
    
    # index because htseq is complaining
    samtools index -@ 10 {pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam # generates BAI-format

    # Check coverage, try to infer about library creation and sequencing
    samtools stats -@ 10 -r /faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC.fasta \
        {pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam > {pathname}/mapped/{FASTQ_BASENAME}.Aligned.sortedByCoord.out.bam.stats
    """.format(FASTQ=FASTQ, pathname=pathname)

