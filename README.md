# ChamberTempRNA: The Effects of Temperature Stress on Gene Regulation of Ostrich Chicks
In this project, I am interested in genetic underpinnings of the chicks’ thermoregulatory responses
to different temperature conditions.
This is the outline of this repository:
## Data
The data for this project is RNA-seq data. Chicks of different age groups were subjected to three temperature treatments - hot, cold and benign - and their blood was sampled. Following RNA extraction and library preparation by PolyA enrichment, the samples were sequenced using Illumina paired-end sequencing, resulting in this many 150 bp reads.
## Analysis
This folder contains scripts for differential expression analysis and WGCNA.
## Results
This folder contains scripts for generation of figures, tables etc. and some summary text files, .htmls and .pdfs
## Steps
This folder contains intermediary files for the different parts of the analysis (QC, mapping)
## Scripts
This folder contains scripts for read QC, trimming, mapping and read counting
### trimming
Trimming was done using TrimGalore, which is a wrapper around cutadapt and fastQC. The minimum length was set to 80 bp and the per-base quality cutoff for read retension was 10. No adaptor sequence was supplied because cutadapt should recognise adapters automatically.
### mapping
STAR was used to map the trimmed and lightly filtered reads to the ostrich reference genome. The resulting BAM files were sorted by coordinate and reads that mapped to more than 4 locations were discarded.
### deduplicating
To asses the influence that putative PCR duplicates have on the sucess of downstream count normalisation and the results of DGE analysis, we did everything downstream both with and without removing duplicate reads. For this purpose, sample L38 was split because it was too big.
### counting
The number of reads that mapped to each gene per sample was counted using htseq-count using the default modes: union (which means that reads do not have to completely overlap with genes in order to be counted for that gene), and non-unique *none* (the read is counted as ambiguous if it aligns to overlapping genes, and as not unique if it aligns to multiple genes)

