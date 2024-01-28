suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))

option_list = list(
  make_option(c("-f", "--infile1"), type="character", default=NULL, 
              help="input file 1 [file]", metavar="character"),
  make_option(c("-g", "--infile2"), type="character", default=NULL, 
              help="input file 2 [file]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile1) || is.null(opt$infile2) ||is.null(opt$outfile) ){
  print_help(opt_parser);
  stop("Please provide -f infile1 -g infile2 and -o outfile option", call.=FALSE);
}
inFile1 = opt$infile1
inFile2 = opt$infile2
outFile = opt$outfile

print(inFile1)
print(inFile2)
print(outFile)



pseudoU_anno_genetype_num<- read.table(inFile1, header=F,sep="\t")
# pseudoU_anno_genetype_num<-pseudoU_anno_genetype_num[order(pseudoU_anno_genetype_num[,2],decreasing = TRUE),]
percent<-round(100*pseudoU_anno_genetype_num$V2/sum(pseudoU_anno_genetype_num$V2),2)
percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$V2,')', sep = "")
pseudoU_anno_genetype_num$V1 <- paste(pseudoU_anno_genetype_num$V1, percent, sep = '')
pseudoU_anno_genetype_num$V1 <- factor(pseudoU_anno_genetype_num$V1, levels = pseudoU_anno_genetype_num$V1)
mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = V2, fill = V1)) + 
		geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
		theme(axis.text = element_blank()) + 
	  	theme(axis.ticks = element_blank()) + 
	  # scale_fill_discrete(labels = pie_data)+
	  	theme(panel.grid = element_blank(), 
	  	panel.background=element_blank(),
	  	axis.text = element_blank(), 
	  	axis.title = element_blank(), 
	  	axis.ticks = element_blank()) + 
	  	guides(fill = guide_legend(title = "gene biotype"))+
	  	scale_fill_manual(values = mycol)

pdf(paste(outFile,"_gene_biotype_piechart.pdf",sep=""))
pie_plot
dev.off()

pseudoU_anno_genetype_num<- read.table(inFile2, header=F,sep="\t")
# pseudoU_anno_genetype_num<-pseudoU_anno_genetype_num[order(pseudoU_anno_genetype_num[,2],decreasing = TRUE),]
percent<-round(100*pseudoU_anno_genetype_num$V2/sum(pseudoU_anno_genetype_num$V2),2)
percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$V2,')', sep = "")
pseudoU_anno_genetype_num$V1 <- paste(pseudoU_anno_genetype_num$V1, percent, sep = '')
pseudoU_anno_genetype_num$V1 <- factor(pseudoU_anno_genetype_num$V1, levels = pseudoU_anno_genetype_num$V1)
mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = V2, fill = V1)) + 
		geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
		theme(axis.text = element_blank()) + 
	  	theme(axis.ticks = element_blank()) + 
	  # scale_fill_discrete(labels = pie_data)+
	  	theme(panel.grid = element_blank(), 
	  	panel.background=element_blank(),
	  	axis.text = element_blank(), 
	  	axis.title = element_blank(), 
	  	axis.ticks = element_blank()) + 
	  	guides(fill = guide_legend(title = "gene feature"))+
	  	scale_fill_manual(values = mycol)

pdf(paste(outFile,"_gene_feature_piechart.pdf",sep=""))
pie_plot
dev.off()