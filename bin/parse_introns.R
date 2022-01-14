#----------------------------------------
# Create a spliced+unspliced transcriptome from genome and GTF files using eisaR,
# see https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
#----------------------------------------

#/ parse arguments from command line:
args=commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
    stop("[invalid params]", "\n",
         "parse_introns.R <genome> <gtf> <gene_name> <gene_id> <gene_type> <chrM> <rrna>")
         
}

params <- list()
params$genome <- args[1]
params$gtf <- args[2]
params$gene_name <- args[3]
params$gene_id <- args[4]
params$gene_type <- args[5]
params$readlength <- 150 # just hardcode it to be consistent, and most sequencing is anyway 2x150
params$chrM <- args[6]
params$rrna <- args[7]

#/ load packages:
suppressMessages({
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(rtracklayer)  
})

#/ read GTF as GRanges:
gtf.gr <- import(params$gtf)

#/ Verify gene_name, gene_id and gene_type from params are part of that GTF:
checkExistCol <- function(x,  # gtf
                          y){ # colname
  if(!y %in% colnames(mcols(x))) {
    stop(
      "\n",
      "-----------------------------------------------------------","\n",
      "[Error] ", y, " is not a column name in the GTF!","\n",
      "-----------------------------------------------------------","\n")
  }
}
checkExistCol(gtf.gr, params$gene_name)
checkExistCol(gtf.gr, params$gene_type)
checkExistCol(gtf.gr, params$gene_id)

#/ check if chrM exists as chromosome and rRNA as gene_type:
if(!params$chrM %in% levels(seqnames(gtf.gr))) 
  stop("value of <chrM> is not a chromosome in the GTF!")

if(params$rrna %in% unique(gtf.gr$gene_type)) 
  stop("value of <rrna> is not a gene type in the GTF!")

#/ extract introns from GTF based on the given read length:
grl <- eisaR::getFeatureRanges(
  gtf=params$gtf, featureType=c("spliced", "intron"), intronType="separate", 
  flankLength=as.numeric(params$readlength)-1, joinOverlappingIntrons=FALSE,verbose=FALSE)

#/ load the genome and make sure name is not a whitespace-delimited mess:
genome <- Biostrings::readDNAStringSet(params$genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)

#/ get transcripts from genome
seqs <- GenomicFeatures::extractTranscriptSeqs(x=genome, transcripts=grl)

#/ write expanded transcriptome (exon+intron)
Biostrings::writeXStringSet(
  seqs, filepath="annotation.expanded.fa.gz", compress=TRUE)

#/ write expanded gtf 
gz <- gzfile("annotation.expanded.gtf.gz", "w")
eisaR::exportToGtf(grl, filepath=gz); close(gz)

#/ write exon/intron mappings
gz <- gzfile("annotation.expanded.features.tsv.gz", "w")
write.table(metadata(grl)$corrgene, file=gz, row.names=FALSE, 
            col.names=TRUE, quote=FALSE, sep="\t"); close(gz)

#/ write transcriptome to gene mapping file
df <- eisaR::getTx2Gene(grl, filepath="annotation.expanded.tx2gene.tsv")

#/ Parse mitochondrial and rRNA genes:
write.table(x=sort(unique(gtf.gr[gtf.gr$gene_type == "rRNA",]$gene_id)),
            file="annotation.rRNA.txt", sep="\n", col.names=FALSE, 
            row.names=FALSE, quote=FALSE)

write.table(x=sort(unique(gtf.gr[seqnames(gtf.gr)==params$chrM,]$gene_id)),
            file="annotation.mtRNA.txt", 
            sep="\n", col.names=FALSE, 
            row.names=FALSE, quote=FALSE)

#/ A data.frame with gene_id, gene_name and gene_type:
df_gene <- data.frame(gtf.gr[gtf.gr$type=="gene"])
write.table(x=df_gene[,c(params$gene_id, params$gene_name, params$gene_type)],
            file="annotation.gene2type.txt", sep="\t", 
            col.names=TRUE, row.names=FALSE, quote=FALSE)

