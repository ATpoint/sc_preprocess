#--------------------------------------------------
# SAMPLESHEET VALIDATION
#--------------------------------------------------

options(stringsAsFactors=FALSE)
list.error <- list()

#/ parse arguments from command line:
args=commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("[invalid params]", "\n",
       "Usage: validate_samplesheet.R <samplesheet>")
  
}

#/ read sample sheet
s <- suppressWarnings(read.delim(args[1], sep=",", header=TRUE))
if(ncol(s)<5) {
  list.error$colerror <- 
    paste("[SAMPLESHEET ERROR]",
          "samplesheet must have 5 columns!", 
          sep="\n")
}
s<-s[,1:5]

#/ validate colnames
if(!all(c("sample", "r1", "r2", "libtype", "is_fb") %in% colnames(s))){
  list.error$namerror <- 
    paste("[SAMPLESHEET ERROR]",
          "samplesheet must have column names <sample> <r1> <r2> <libtype> <is_fb>",
          sep="\n")
}

#/ validate that $1-3 are filled with non-empty content
s_empty <- s[,1:5]
s_empty[is.na(s_empty)] <- ""
s_rowsums <- rowSums(s_empty=="")
s_rowsums_is_empty <- which(s_rowsums>0)

if(length(s_rowsums_is_empty) > 0){
  list.error$emptyerror <-
    paste("[SAMPLESHEET ERROR]",
          "samplesheet has empty entries in the first three columns of these rows:",
          paste(s_rowsums_is_empty, collapse=","), sep="\n")
}

#/ validate that $4 is yes/no or TRUE/FALSE, or empty
if(!all(s$is_fb %in% c("", "true", "TRUE", "false", "FALSE", NA))){
  list.error$trueerror <-
    paste("[SAMPLESHEET ERROR",
          "Column 4 of the samplesheet must be either TRUE/FALSE or true/false, or empty",
          sep="\n")
}


#/ validate that there is not a SF library without matched non-SF:
check_matched <- lapply(unique(s$sample), function(x){
  
  r <- NULL
  sx <- s[s$sample==x,]
  sx_name <- sx[sx$is_fb=="true",]$sample
  if(length(sx_name)!=0){
    if(!"false" %in% sx[sx$sample==sx_name,]$is_fb) r <- sx_name
  }
  return(r)
})

not_matched <- setdiff(unlist(check_matched), NULL)
if(length(not_matched) > 0){
  list.error$fberror <-
    paste("[SAMPLESHEET ERROR]",
          "These samples only have feature barcode libraries and no matched transcriptomic libraries:",
          paste(not_matched, collapse=", "),
          sep="\n")

}

#/ validate that the r1 and r2 entries are unique:
s_tab <- table(c(s$r1, s$r2))
s_tab_nonunique <- s_tab[s_tab>1]
if(length(s_tab_nonunique)>0){
  list.error$nonuniqueerror <- 
    paste("[SAMPLESHEET ERROR]",
          "These r1/r2 entries appear more than once:",
          paste(names(s_tab_nonunique), collapse="\n"),
          sep="\n")
}

#/ exit if entries in list.error
if(length(list.error)>0){
  
  for(i in list.error) {
    message("")
    message(i)
    message("")
  }
  
  quit(status=1, save="no")
  
}
  
#/ ensure consistent naming
s$is_fb[s$is_fb %in% c("TRUE", "true", "yes")] <- "true"
s$is_fb[s$is_fb %in% c("", "false", "no")] <- "false"

#/ write validated sheet:
write.table(s, file="samplesheet_validated.txt", quote=FALSE, 
            col.names=TRUE, row.names=FALSE, sep=",")
       