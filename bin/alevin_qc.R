#--------------------------------------------------
# QC reports using alevinQC
#--------------------------------------------------

suppressMessages({
  library(alevinQC)
})

args=commandArgs(trailingOnly=TRUE)
if (length(args)!= 3) {
  stop("[invalid params]", "\n",
       "Usage: alevin_qc.R <basename> <paths_to_alevins> <suffix>")
}

if(args[3]=="empty") {
  suffix <- ""
} else suffix <- args[3]
sampleid <- paste0(args[1], suffix)

alevinQCReport(baseDir=args[2], sampleId=sampleid, outputFile=paste0(sampleid, ".html"),
               outputDir=getwd(), outputFormat="html_document", quiet=TRUE)