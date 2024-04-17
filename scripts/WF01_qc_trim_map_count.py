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

pathname = '/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/00_fastq'
out_dir = '/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/steps/QC'
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
    FASTQ_BASENAME = FASTQ.split(".")[0] #name without extensions
    gwf.target( f'FastQC_{FASTQ_BASENAME}', # name of the target
           cores=1,
           memory='8gb',
           walltime='00:05:00',	
           inputs=[f'../data/azenta/{FASTQ}'], 
           outputs=[f'{out_dir}/{FASTQ_BASENAME}_fastqc.html', 
                    f'.{out_dir}/{FASTQ_BASENAME}_fastqc.zip']) << """ 
    eval “$(conda shell.bash hook)”
    conda activate QC                
    fastqc ../data/azenta/{FASTQ}
    """.format(FASTQ=FASTQ, out_dir='/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/steps/QC')

#############################
### Wait for all fastqc outputs
### as input for the next step
### where we run multiqc on 
#############################


fastq_basenames = [ FASTQ.split(".")[0] for FASTQ in fastq_files ] #names without extensions
#print("Just the fastq basenames")
#print(fastq_basenames)

gwf.target( 'MultiQC', #name of the target
           cores=1,
           memory='8gb',
           walltime='00:05:00',	
           inputs= {'ZIP': [f'{out_dir}/{FASTQ_BASENAME}_fastqc.zip' for FASTQ_BASENAME in fastq_basenames],
                    'REPORTS': [f'{out_dir}/{FASTQ_BASENAME}_fastqc.html' for FASTQ_BASENAME in fastq_basenames] }, 
           outputs=[f'{out_dir}/multiqc_report.html'],
           # everything except the endpoints (targets which not dependencies of others) will be removed (cleaned),
           # unless you protect it (next line)
           protect=[f'{out_dir}/multiqc_report.html']) << """
    eval “$(conda shell.bash hook)”
    conda activate QC 
    multiqc -f --outdir ../steps/QC \
            ../steps/QC/
    """

#############################
### Trimming step with Trim Galore
### for each pair of fastq files in parallel
#############################


fastq_pairs = np.unique([ FASTQ.split('_')[0] for FASTQ in fastq_basenames ]) #names without pair number and extension


for FASTQ_PAIR in fastq_pairs:
    gwf.target( f"TrimGalore_{FASTQ_PAIR}", #name of the target
           cores=4,
           memory='8gb',
           walltime='00:10:00',	
           inputs= {'FASTQS': [f'{pathname}/{FASTQ_PAIR}_1.fastq.gz', # OBS! most probably needs adjustment!
                               f'{pathname}/{FASTQ_PAIR}_2.fastq.gz'],
                    'MULTIQC': ['results/multiqc_output/multiqc_report.html']}, 
           outputs=[f'../data/azenta/trimmed/{FASTQ_PAIR}_1_trimmed.fq.gz', 
                    f'../data/azenta/trimmed/{FASTQ_PAIR}_2_trimmed.fq.gz',
                    f'{FASTQ_PAIR}_1.fastq.gz_trimming_report.txt',
                    f'{FASTQ_PAIR}_1_val_1_fastqc.html',
                    f'{FASTQ_PAIR}_2.fastq.gz_trimming_report.txt',
                    f'{FASTQ_PAIR}_2_val_2_fastqc.html']) << """
    eval “$(conda shell.bash hook)”
    conda activate broiler_test
    mkdir -p ../data/azenta/trimmed
    trim_galore --paired --fastqc --gzip --length 80 --cores 4 -q 10 \
        -o /faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/azenta/trimmed/ \
        {pathname}/{FASTQ_PAIR}_1.fastq.gz \
        {pathname}/{FASTQ_PAIR}_2.fastq.gz
    """.format(FASTQ_PAIR=FASTQ_PAIR, pathname='/faststorage/project/ostrich_thermal/people/leah/ChamberTempRNA/data/azenta')
