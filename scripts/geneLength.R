#!/bin/Rscript
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=02:00:00

# I need the gene length to calculate transcripts per million (tpm)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gtf"
FASTAfile = "/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed/Struthio_camelus_HiC.fasta"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="Struthio_camelus_HiC", asRangedData=F, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file="GC_lengths.tsv", sep="\t")