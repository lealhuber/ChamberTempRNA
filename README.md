# ChamberTempRNA: The Effects of Temperature Stress on Gene Regulation of Ostrich Chicks
In this project, I am interested in genetic underpinnings of the chicks’ thermoregulatory responses to different temperature conditions. Chicks of two age groups were subjected to three temperature treatments - hot, cold and benign - and their blood was sampled. Following RNA extraction and library preparation by PolyA enrichment, the samples were sequenced using Illumina paired-end sequencing.
The scripts folder contains all the scripts of the analysis. The main script is gwf_trim_map_count.py. Here are some descriptions of some of the steps:
### Trimming
Trimming was done using TrimGalore, which is a wrapper around cutadapt and fastQC. The minimum length was set to 80 bp and the per-base quality cutoff for read retension was 10. No adaptor sequence was supplied because cutadapt should recognise adapters automatically.
### Mapping
STAR was used to map the trimmed and lightly filtered reads to the ostrich reference genome. The resulting BAM files were sorted by coordinate and reads that mapped to more than 4 locations were discarded.
### Deduplicating
To asses the influence that putative PCR duplicates have on the sucess of downstream count normalisation and the results of DGE analysis, we did everything downstream both with and without removing duplicate reads. For this purpose, sample L38 was split because it was too big.
### Counting
The number of reads that mapped to each gene per sample was counted using htseq-count using the default modes: union (which means that reads do not have to completely overlap with genes in order to be counted for that gene), and non-unique *none* (the read is counted as ambiguous if it aligns to overlapping genes, and as not unique if it aligns to multiple genes).
### Differential gene expresssion
Reads were filtered and Globin reads removed, then they were normalised and the dispersion calculated with the glmmSeq_prep_withchicks.R script. Then, different models were run with the glmmSeq_models.R script. The differentially expressed genes from these models were extracted and further analyses with a post hoc analysis wit the glmmPostHoc.R script. In there is also the code for the figures. Some additional age group related analyses are in PH_ageefftable.R.
### Cloacal temperature analysis
The effects of treatment, age, sex and size on cloacal temperature were analysed using chamber_exp_clean.R
### Gene overrepresentation analysis
Can be found in GSEA_glmmSeq_final.R, with functional annotations generated with eggnog.sh
