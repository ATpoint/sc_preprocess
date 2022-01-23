#----------------------------------------
# Read Alevin quantifications and split into spliced and unspliced counts.
# If feature barcodes are present read also these and optionally perform barcode translation.
# Output the RNA experiment as ".mtx.gz" if no FBs present, or "_rna.mtx.gz" if FBs present.
# Output FBs as "_fb.mtx.gz"
# If FBs present then output ".mtx.gz" which only contain CBs present in both RNA and FBs.
# In any case the final "basename.mtx.gz" is always the one for downstream analysis while
# "basename_rna/fb.mtx.gz" is only returned to control the quality in case the matched output
# is super bad and one wants to check whether the crappiness comes from the RNA or FB experiment.
# AlevinQC is (in a separate process) returned for both RNA and FB separately.

#---------------------------------------------------------------------------------------------------------
# Packages and options from CL
#---------------------------------------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(fishpond)
  library(Matrix)
  library(optparse)
  library(rtracklayer)
  library(tximeta)
  library(SummarizedExperiment)
}))

o_in <- 
  parse_args(OptionParser(option_list=list( 
    make_option(c("--sampleid"), action="store", default="NULL"),
    make_option(c("--alevin_rna"), action="store", default="NULL"),
    make_option(c("--alevin_fb"), action="store", default="NULL"),
    make_option(c("--expanded_features"), action="store", default="NULL"),
    make_option(c("--gene2type"), action="store", default="NULL"),
    make_option(c("--translate_table"), action="store", default="NULL"),
    make_option(c("--do_translate"), action="store_true", default=FALSE)
  )))
opts <- sapply(o_in, function(x) if(x %in% c("NULL", "null")) NULL else x, simplify=FALSE)

#/ when no feature barcodes are present for that sampleid then the incoming nextflow channel will be empty,
#/ so set it to NULL
if(length(opts$alevin_fb)==0) opts$alevin_fb <- NULL

save.image("env.rdata")

#---------------------------------------------------------------------------------------------------------
# Read RNA experiment -- split into spliced and unspliced
#---------------------------------------------------------------------------------------------------------

#/ from here it is mainly https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
#/ starting from what they call "step 4"
cg <- read.delim(opts$expanded_features, header=TRUE, as.is=TRUE)
colnames(cg)[colnames(cg) == "intron"] <- "unspliced"

#/ Read quantifications with tximeta:
se_rna <- tximeta::tximeta(
  coldata = data.frame(names=opts$sampleid, 
                       files=paste0(opts$alevin_rna,"/alevin/quants_mat.gz"),
                       stringsAsFactors=TRUE),
                       type="alevin", skipMeta=TRUE, dropInfReps=TRUE)

#/ split into spliced- and unspliced:
se_rna <- tximeta::splitSE(se_rna, cg, assayName = "counts")
assayNames(se_rna)[1] <- "counts"

#/ add the gene_id, gene_name and gene_type data.frame to the rowData:  
gene2type <- read.delim(opts$gene2type, header=TRUE, as.is=TRUE)

m <- merge(x=data.frame(gene_id=rownames(se_rna)), 
           y=data.frame(gene2type),
           by.x="gene_id", by.y="gene_id", sort=FALSE, all.x=TRUE, all.y=FALSE)
  
rowData(se_rna)  <- DataFrame(m)

#---------------------------------------------------------------------------------------------------------
# Optionally: Read feature barcode (FB) experiment and perform (again optionally) barcode translation
#---------------------------------------------------------------------------------------------------------

if(!is.null(opts$alevin_fb)){
  
  se_fb <- tximeta::tximeta(
    coldata = data.frame(names=opts$sampleid, 
                         files=paste0(opts$alevin_fb,"/alevin/quants_mat.gz"),
                         stringsAsFactors=TRUE),
    type="alevin", skipMeta=TRUE, dropInfReps=TRUE)
  
  #/ optional barcode translation with the 10x 3M translation table
  if(opts$do_translate){
    
    map <- read.table(opts$translate_table, sep="\t", header=FALSE) # read translation table
    fb2rna <- map$V1                                # the sequences are rna
    names(fb2rna) <- map$V2                         # the names are fb
    subset_colnames <- colnames(se_fb)[colnames(se_fb) %in% names(fb2rna)] 
    se_fb <- se_fb[,subset_colnames]                # subset to colnames that are in the 3M list
    colnames(se_fb) <- fb2rna[colnames(se_fb)]      # if fb names are matched return rna CBs
    rm(map, fb2rna); invisible(gc(verbose=FALSE))
    
  }
  
  #/ match se_rna and se_fb: take intersect and make roworder the same
  intersected <- intersect(colnames(se_fb), colnames(se_rna))
  if(length(intersected)==0) {
    warning(paste("No overlap between RNA and feature barcode experiment for sample:", opts$sampleid)) 
  }
  
  ncells=data.frame(rna=ncol(se_rna), fb=ncol(se_fb), overlap=length(intersected))
  se_fb  <- se_fb[,intersected]
  se_rna <- se_rna[,intersected]
  
  #/ not sure whether SummarizedExperiments do that internally...ensure matching columns
  if(!all(colnames(se_rna)==colnames(se_fb))){
    matched <- as.numeric(na.omit(match(colnames(se_rna), colnames(se_fb))))
    se_fb <- se_fb[,matched]
  }
  
  assay_spliced   <- rbind(assay(se_rna, "counts"), assay(se_fb))
  assay_unspliced <- rbind(assay(se_rna, "unspliced"), assay(se_fb))
  barcodes        <- data.frame(barcodes=colnames(se_rna))
  features        <- data.frame(
                       rbind(rowData(se_rna),
                       DataFrame(gene_id=rownames(se_fb), gene_name=rownames(se_fb), gene_type="feature_barcode")))
  rm(se_fb)
  
} else {
  
  assay_spliced   <- assay(se_rna, "counts")
  assay_unspliced <- assay(se_rna, "unspliced")
  barcodes        <- data.frame(barcodes=colnames(se_rna))
  features        <- data.frame(rowData(se_rna))
  ncells=data.frame(rna=ncol(se_rna), fb=NA, overlap=NA)
  
}

rm(se_rna)
invisible(gc(verbose=FALSE))

#---------------------------------------------------------------------------------------------------------
# Save expression matrix, features and barcodes as mtx.gz or tsv.gz
#---------------------------------------------------------------------------------------------------------

#/ Use writeMM from the Matrix package to save expression matrix as mtx
spl    <- paste0(opts$sampleid, "_spliced.mtx")
unspl  <- paste0(opts$sampleid, "_unspliced.mtx")

invisible(writeMM(obj=assay_spliced, file=spl))
invisible(writeMM(obj=assay_unspliced, file=unspl))

system(command=paste("gzip -f", spl))
system(command=paste("gzip -f", unspl))

#/ Save barcodes and features
file.colnames <- gzfile(paste0(opts$sampleid, "_barcodes.tsv.gz"), "w")
write.table(x=barcodes, file=file.colnames, col.names=FALSE, row.names=FALSE, quote=FALSE)
close(file.colnames)

file.rowdata <- gzfile(paste0(opts$sampleid, "_features.tsv.gz"), "w")
write.table(x=features, file=file.rowdata, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
close(file.rowdata)

#/ write final number of cells:
write.table(cbind(data.frame(sample=opts$sampleid), ncells),
            file=paste0(opts$sampleid, "_ncells.txt"), 
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
