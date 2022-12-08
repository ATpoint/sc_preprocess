# Simple sort operation to get the command line collections into order
# based on the name of the first column

args = commandArgs(trailingOnly=TRUE)
indata <-args[1]

d <- read.delim(indata)

dx <- which(d$X!="")

df <- data.frame(matrix(ncol=2, nrow=0))
colnames(df) <- c("d1", "d2")
for(i in 1:length(dx)){
  
  d1 <- dx[i]
  d2 <- if(i==length(dx)) nrow(d) else dx[i+1]-1
  df <- rbind(df, data.frame(d1=d1, d2=d2))
  
}

df$first_column <- d[,1][dx]

dd <- df[order(df$first_column),]

final <- d[unlist(lapply(1:nrow(dd), function(x) dd[x,1]:dd[x,2])),]

write.table(final, file="command_lines.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
