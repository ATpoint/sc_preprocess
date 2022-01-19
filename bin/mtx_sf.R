#----------------------------------------
# Read quantifications from Alevin, split into exonic and intronic parts and save as mtx.gz
#----------------------------------------

#/ parse arguments from command line:
args=commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("[invalid params]", "\n",
       "mtx.R <basename> <alevin_quant> <suffix>")
}

params <- list()
params$basename <- args[1]
params$alevin_quant <- args[2]
params$suffix <- args[3]

#/ load packages:
suppressWarnings(suppressMessages({
  library(fishpond)
  library(Matrix)
  library(tximeta)
  library(SummarizedExperiment)
}))

#/ Read quantifications with tximeta:
se <- tximeta::tximeta(
  coldata = data.frame(names=params$basename, 
                       files=paste0(params$alevin_quant,"/alevin/quants_mat.gz"),
                       stringsAsFactors=TRUE),
  type="alevin", skipMeta=TRUE, dropInfReps=TRUE)

#/ Save to disk as mtx and compress:
spl <- paste0(params$basename, params$suffix, ".mtx")

##/ Use the Matrix package for the conversion of the sparse matrix to mtx:
invisible(writeMM(obj=assay(se), file=spl))

#/ run gzip:
system(command = paste("gzip", spl))

#/ Save col- and rowdata and compress:
file.colnames <- gzfile(paste0(params$basename, params$suffix, "_barcodes.tsv.gz"), "w")
write.table(x=data.frame(barcodes=colnames(se)), file=file.colnames, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
close(file.colnames)

file.rowdata <- gzfile(paste0(params$basename, params$suffix, "_features.tsv.gz"), "w")
write.table(x=data.frame(features=rownames(se)), file=file.rowdata,
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
close(file.rowdata)
