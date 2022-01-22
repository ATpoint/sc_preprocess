#--------------------------------------------------
# SAMPLESHEET VALIDATION
#--------------------------------------------------

options(stringsAsFactors=FALSE)
list.error <- list()

#/ parse arguments from command line:
args=commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("[invalid params]", "\n",
       "Usage: validate_samplesheet.R <samplesheet> <basedir> <launchdir> <projectdir>")
  
}

# The *dirs are the absolute paths of the nf variables $baseDir, $launchDir and $projectDir.
# If the samplesheet contains any of these three as strings then this script will replace them
# with the absolute path. This is currently necessary because the samplesheet content (R1/R2 paths)
# is parsed as strings so nf variables such as $baseDir do not get expanded.

#/ read sample sheet
s <- suppressWarnings(read.delim(args[1], sep=",", header=TRUE))
if(ncol(s)<4) {
  list.error$colerror <- 
    paste("[SAMPLESHEET ERROR]",
          "samplesheet must have 4 columns!", 
          sep="\n")
}
s<-s[,1:4]

#/ validate colnames
if(!all(c("sample_id", "R1", "R2", "is_fb") %in% colnames(s))){
  list.error$namerror <- 
    paste("[SAMPLESHEET ERROR]",
          "samplesheet must have column names <sample_id> <R1> <R2> <is_fb>",
          sep="\n")
}

#/ validate that $1-3 are filled with non-empty content
s_empty <- s[,1:3]
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
if(!all(s[,4] %in% c("", "true", "TRUE", "false", "FALSE", NA))){
  list.error$trueerror <-
    paste("[SAMPLESHEET ERROR",
          "Column 4 of the samplesheet must be either TRUE/FALSE or true/false, or empty",
          sep="\n")
}

#/ if no or FALSE/false or empty set to "false"
s[,4][is.na(s[,4])] <- "false"

#/ validate that there is not a SF library without matched non-SF:
check_matched <- lapply(unique(s$sample_id), function(x){
  
  r <- NULL
  sx <- s[s$sample_id==x,]
  sx_name <- sx[sx$is_fb=="true",]$sample_id
  if(length(sx_name)!=0){
    if(!"false" %in% sx[sx$sample_id==sx_name,]$is_fb) r <- sx_name
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

#/ validate that the R1 and R2 entries are unique:
s_tab <- table(c(s$R1, s$R2))
s_tab_nonunique <- s_tab[s_tab>1]
if(length(s_tab_nonunique)>0){
  list.error$nonuniqueerror <- 
    paste("[SAMPLESHEET ERROR]",
          "These R1/R2 entries appear more than once:",
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
  
#/ ensure consistent naming in $4
s$is_fb[s$is_fb %in% c("TRUE", "true", "yes")] <- "true"
s$is_fb[s$is_fb %in% c("", "false", "no")] <- "false"

#/ replace $baseDir, $launchDir or $projectDir with the respective paths:
s[,2] <- gsub("\\$baseDir|\\$\\{baseDir\\}", as.character(args[2]), s[,2])
s[,3] <- gsub("\\$baseDir|\\$\\{baseDir\\}", as.character(args[2]), s[,3])

s[,2] <- gsub("\\$launchDir|\\$\\{launchDir\\}", as.character(args[3]), s[,2])
s[,3] <- gsub("\\$launchDir|\\$\\{launchDir\\}", as.character(args[3]), s[,3])

s[,2] <- gsub("\\$projectDir|\\$\\{baseDir\\}", as.character(args[4]), s[,2])
s[,3] <- gsub("\\$projectDir|\\$\\{projectDir\\}", as.character(args[4]), s[,3])

#/ write validated sheet:
write.table(s, file="samplesheet_validated.txt", quote=FALSE, 
            col.names=TRUE, row.names=FALSE, sep=",")
       