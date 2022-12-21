#------------------------------------------------------------
# Run alevinQC -- that's it. Nothing more.
#------------------------------------------------------------

args=commandArgs(trailingOnly=TRUE)
sampleid <- args[1]
alevinQC::alevinQCReport(baseDir=args[2], sampleId=sampleid, outputFile=paste0(sampleid, ".html"), outputDir=getwd(), outputFormat="html_document", quiet=TRUE)