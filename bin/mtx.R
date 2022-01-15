#----------------------------------------
# Read quantifications from Alevin, split into exonic and intronic parts and save as mtx.gz
#----------------------------------------

#/ parse arguments from command line:
args=commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("[invalid params]", "\n",
         "mtx.R <basename> <alevin_quant> <expanded_features> <gene2type>")
}

params <- list()
params$basename <- args[1]
params$alevin_quant <- args[2]
params$expanded_features <- args[3]
params$gene2type <- args[4]

#/ load packages:
suppressMessages({
  library(fishpond)
  library(Matrix)
  library(rtracklayer)
  library(tximeta)
  library(SummarizedExperiment)
})

#/ Read the table mapping intronic to exonic features:
cg <- read.delim(params$expanded_features, header=TRUE, as.is=TRUE)
colnames(cg)[colnames(cg) == "intron"] <- "unspliced"

#/ Read quantifications with tximeta:
se <- tximeta::tximeta(
  coldata = data.frame(names=params$basename, 
                       files=paste0(params$alevin_quant,"/alevin/quants_mat.gz"),
                       stringsAsFactors=TRUE),
                       type="alevin", skipMeta=TRUE, dropInfReps=TRUE)

#/ split into spliced- and unspliced:
se <- tximeta::splitSE(se, cg, assayName = "counts")
assayNames(se)[1] <- "counts"

#/ add the gene_id, gene_name and gene_type data.frame to the rowData:  
gene2type <- read.delim(params$gene2type, header=TRUE, as.is=TRUE)

m <- merge(x=data.frame(gene_id=rownames(se)), 
           y=data.frame(gene2type),
           by.x="gene_id", by.y="gene_id", sort=FALSE, all.x=TRUE, all.y=FALSE)
  
rowData(se)  <- DataFrame(m)

#/ Save to disk as mtx and compress:
spl    <- paste(params$basename, "spliced.mtx", sep="_")
unspl <-  paste(params$basename, "unspliced.mtx", sep="_")

##/ Use the Matrix package for the conversion of the sparse matrix to mtx:
invisible(writeMM(obj=assay(se, "counts"), file=spl))
invisible(writeMM(obj=assay(se, "unspliced"), file=unspl))

#/ run gzip:
system(command = paste("gzip", spl))
system(command = paste("gzip", unspl))

#/ Save col- and rowdata and compress:
file.colnames <- gzfile(paste(params$basename, "barcodes.tsv.gz", sep="_"), "w")
write.table(x=data.frame(barcodes=colnames(se)),
            file=file.colnames, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
close(file.colnames)

file.rowdata <- gzfile(paste(params$basename, "features.tsv.gz", sep="_"), "w")
write.table(x=data.frame(rowData(se)), file=file.rowdata,
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
close(file.rowdata)