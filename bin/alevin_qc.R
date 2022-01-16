#--------------------------------------------------
# QC reports using alevinQC
#--------------------------------------------------

suppressMessages({
  library(alevinQC)
})

args=commandArgs(trailingOnly=TRUE)
if (length(args)!= 2) {
  stop("[invalid params]", "\n",
       "Usage: alevin_qc.R <basename> <paths_to_alevins>")
}

alevinQCReport(baseDir=args[2], sampleId=args[1], outputFile=paste0(args[1], ".html"),
               outputDir=getwd(), outputFormat="html_document", quiet=TRUE)