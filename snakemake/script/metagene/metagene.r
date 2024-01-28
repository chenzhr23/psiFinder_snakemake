# reference https://github.com/olarerin/metaPlotR
suppressMessages(library("ggplot2"))
suppressMessages(library("optparse"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("scales"))
suppressMessages(library("dplyr"))

option_list = list(
  make_option(c("-f", "--metagenefile"), type="character", default=NULL, 
              help="distance of single sites [file]", metavar="character"),
  make_option(c("-o", "--outfile_prefix"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$metagenefile)|| is.null(opt$outfile_prefix) ){
  print_help(opt_parser);
  stop("Please provide -f metagenefile and -o outfile_prefix option", call.=FALSE);
}

Metagenefile = opt$metagenefile
outFile_prefix = opt$outfile_prefix

print(Metagenefile)
print(outFile_prefix)

pseudoU.dist <- read.delim(Metagenefile, header = T)
# Determine longest length transcript for each gene
trx_len <- pseudoU.dist$utr5_size + pseudoU.dist$cds_size + pseudoU.dist$utr3_size 
temp <- data.frame(paste(pseudoU.dist$chr,pseudoU.dist$coord,sep="_"), pseudoU.dist$refseqID, trx_len)
colnames(temp) <- c("coord", "gid", "trx_len") 
temp.df <- temp[order(temp$coord,  temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$coord),]

# m6a data to one transcript per gene (longest)
pseudoU.dist <- pseudoU.dist[pseudoU.dist$refseqID %in% temp.df$gid,]
pseudoU.dist$metagene_feature<-case_when(
  0 <= pseudoU.dist$rel_location & pseudoU.dist$rel_location < 1 ~ "5'UTR",
  1<= pseudoU.dist$rel_location & pseudoU.dist$rel_location < 2 ~ "CDS",
  2<= pseudoU.dist$rel_location & pseudoU.dist$rel_location <= 3 ~ "3'UTR"
)
table(pseudoU.dist$metagene_feature)
write.table(pseudoU.dist,paste(outFile_prefix,"_pseudoU.dist.uniq.txt",sep=""),row.names=F,quote=F)


####standard#####
metagene_pseudoU<- ggplot(pseudoU.dist,aes(x=rel_location))+
geom_density(alpha=0.8, color = "black",size=0.6,fill = "lightblue")+
geom_vline(xintercept = 1:2, col = brewer.pal(3, "Set1")[1],linetype="dashed")+
theme_classic()+
theme(legend.position="top",
  legend.title=element_blank(),
  panel.background = element_blank(),
  axis.title.x=element_blank())+
scale_x_continuous(limits = c(0, 3),expand=c(0,0))+
scale_y_continuous(expand=c(0,0))+
theme(axis.text.x = element_text(face="bold",size=14,hjust=1.8))+
theme(plot.title = element_text(hjust = 0.5),plot.margin=unit(c(2,2,2,2),units="cm"))
# +scale_x_discrete(limits = c('5_UTR','CDS  ','3_UTR'))
print(paste(outFile_prefix,"_metagene_pseudoU_standard.pdf",sep=""))
pdf(paste(outFile_prefix,"_metagene_pseudoU_standard.pdf",sep=""))
metagene_pseudoU
dev.off()


####normalize by region length#####
utr5.SF <- median(pseudoU.dist$utr5_size, na.rm = T)/median(pseudoU.dist$cds_size, na.rm = T)
utr3.SF <- median(pseudoU.dist$utr3_size, na.rm = T)/median(pseudoU.dist$cds_size, na.rm = T)

# assign the regions to new dataframes
utr5.pseudoU.dist <- pseudoU.dist[pseudoU.dist$rel_location < 1, ]
cds.pseudoU.dist <- pseudoU.dist [pseudoU.dist$rel_location < 2 & pseudoU.dist$rel_location >= 1, ]
utr3.pseudoU.dist <- pseudoU.dist[pseudoU.dist$rel_location >= 2, ]


# rescale 5'UTR and 3'UTR
utr5.pseudoU.dist$rel_location <- rescale(utr5.pseudoU.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.pseudoU.dist$rel_location <- rescale(utr3.pseudoU.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
pseudoU.metagene.coord <- data.frame(norm_value=c(utr5.pseudoU.dist$rel_location, cds.pseudoU.dist$rel_location, utr3.pseudoU.dist$rel_location),metagene_feature=c(rep("5'UTR",length(utr5.pseudoU.dist$rel_location)),rep("CDS",length(cds.pseudoU.dist$rel_location)),rep("3'UTR",length(utr3.pseudoU.dist$rel_location))))
pseudoU.metagene.coord<-arrange(pseudoU.metagene.coord,norm_value)

metagene_pseudoU<- ggplot(pseudoU.metagene.coord,aes(x=norm_value))+
geom_density(alpha=0.8, color = "black",size=0.6,fill = "lightblue")+
geom_density()+
theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))+
theme(axis.title.x=element_blank()) +
# scale_x_discrete(limits = c('5_UTR','CDS','3_UTR'))+
theme(axis.text.x = element_text(face="bold",size=8),axis.text.y = element_text(face="bold",size=8))+
scale_x_continuous(expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0)) +
theme(plot.margin=unit(c(2,2,2,2),units="cm"),legend.background = element_rect(colour = 'grey', fill = 'white', linetype='dashed'))+
geom_vline(xintercept = 1:2, col = "grey",size=0.7,linetype="dashed")+
geom_area(
    aes(x = stage(norm_value, after_stat = oob_censor(x, c(0, 1))),
    fill="UTR5"),
    stat = "density"
  )+
geom_area(
    aes(x = stage(norm_value, after_stat = oob_censor(x, c(1, 2))),
    fill="CDS"),
    stat = "density"
    # fill=brewer.pal(6,"Accent")[5],
  )+
geom_area(
    aes(x = stage(norm_value, after_stat = oob_censor(x, c(2, max(pseudoU.metagene.coord$norm_value)))),
    fill="UTR3"),
    stat = "density"
  )+
scale_fill_manual(values=c(UTR5=brewer.pal(6,"Set3")[4],CDS=brewer.pal(6,"Set3")[5],UTR3=brewer.pal(6,"Set3")[6]))


print(paste(outFile_prefix,"_metagene_pseudoU_norm_length.pdf",sep=""))
pdf(paste(outFile_prefix,"_metagene_pseudoU_norm_length.pdf",sep=""))
metagene_pseudoU
dev.off()
