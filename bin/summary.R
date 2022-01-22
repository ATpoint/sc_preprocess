#------------------------------------------------------------
# Number of detected cells after preprocessing per sample
# as table and barplot
#------------------------------------------------------------

f <- 
  lapply(list.files(getwd(), pattern="_ncells.txt"), function(x){
  read.delim(x, header=FALSE)
}) 

#/ summary table:
d <- do.call(rbind, f)
colnames(d) <- c("sample", "ncells_rna", "ncells_fb", "ncells_intersect")

#/ plot
library(ggplot2)
library(reshape2)
r <- reshape2::melt(d, id.vars="sample")
r$variable <- gsub("ncells_intersect", "Intersect (RNA -- Feature Barcodes)",
                   gsub("ncells_fb", "Feature Barcode Experiment",
                      gsub("ncells_rna", "RNA experiment", r$variable)))
r$variable <- factor(r$variable, levels=unique(r$variable))

if(all(is.na(r[r$variable=="Feature Barcode Experiment",]$value))){
  r <- r[r$variable=="RNA experiment",]
  r$variable <- droplevels(r$variable)
}

cols <- c("#E69F00", "#56B4E9", "#009E73")
cols <- cols[1:length(unique(r$variable))]

theme_set(theme_bw(base_size=10))
gg <- 
  ggplot(r, aes(x=sample, y=value, fill=variable)) +
    geom_bar(position="dodge", stat="identity") +
    xlab(element_blank()) + ylab("number of detected cells") +
    theme(legend.title=element_blank(), legend.position="top", legend.justification="left",
          axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
    scale_fill_manual(values=cols) +
    guides(fill=guide_legend(nrow=3), x=guide_axis(angle=45)) +
    facet_wrap(~sample, ncol=4, scales="free_x")

if(length(cols)==1) gg <- gg + theme(legend.position="none")

write.table(d, "summary_detected_cells.txt", sep="\t", col.names=TRUE, 
            row.names=FALSE, quote=FALSE)

pdf("summary_detected_cells.pdf"); print(gg); dev.off()