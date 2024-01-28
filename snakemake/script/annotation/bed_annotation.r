#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("bedr"))


option_list = list(
  make_option(c("-f", "--infile1"), type="character", default=NULL, 
              help="input file 1 [file]", metavar="character"),
  make_option(c("-g", "--infile2"), type="character", default=NULL, 
              help="input file 2 [file]", metavar="character"),
  make_option(c("-s", "--infile3"), type="character", default=NULL, 
              help="input file 3 [file]", metavar="character"),
  make_option(c("-e", "--infile4"), type="character", default=NULL, 
              help="input file 4 [file]", metavar="character"),
  make_option(c("-i", "--infile5"), type="character", default=NULL, 
              help="input file 5 [file]", metavar="character"),
  make_option(c("-j", "--infile6"), type="character", default=NULL, 
              help="input file 6 [file]", metavar="character"),
  make_option(c("-k", "--infile7"), type="character", default=NULL, 
              help="input file 7 [file]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile1) || is.null(opt$infile2) ||is.null(opt$infile3) ||is.null(opt$infile4) ||is.null(opt$infile5) ||is.null(opt$infile6) ||is.null(opt$infile7) ||is.null(opt$outfile) ){
  print_help(opt_parser);
  stop("Please provide -f infile1 -g infile2 and -s infile3, -e infile4, -i infile5 , -j infile6, -k infile7, and -o outfile option", call.=FALSE);
}
anno.biotype.bed.file = opt$infile1
input.bed.file = opt$infile2
hg38.psiU.SingleSites.bed.file = opt$infile3
human.hg38.Pseudo.result.col29.xlsx.file = opt$infile4
genome_fasta = opt$infile5
snoRNA_fasta_rRNA = opt$infile6
snoRNA_fasta_snRNA = opt$infile7
outFile = opt$outfile

print(anno.biotype.bed.file)
print(input.bed.file)
print(hg38.psiU.SingleSites.bed.file)
print(human.hg38.Pseudo.result.col29.xlsx.file)
print(genome_fasta)
print(snoRNA_fasta_rRNA)
print(snoRNA_fasta_snRNA)
print(outFile)

anno.biotype.bed<-read.table(anno.biotype.bed.file,sep="\t")
anno.biotype.bed_uniqid<-paste(anno.biotype.bed$V1,anno.biotype.bed$V2,anno.biotype.bed$V3,anno.biotype.bed$V6,sep="_")
input.bed<-read.table(input.bed.file,sep="\t")
input.bed_uniqid<-paste(input.bed$V1,input.bed$V2,input.bed$V3,input.bed$V6,sep="_")
probe<-which(input.bed_uniqid %in% anno.biotype.bed_uniqid)
# probe<-which(input.bed$V4 %in% anno.biotype.bed$V4)
input.bed<-input.bed[probe,]
input.bed$input.bed_uniqid<-input.bed_uniqid
anno.biotype.bed$anno.biotype.bed_uniqid<-anno.biotype.bed_uniqid
# add_seq <- anno.biotype.bed %>% left_join(input.bed,by=c("V4"="V4"))
add_seq <- anno.biotype.bed %>% left_join(input.bed,by=c("anno.biotype.bed_uniqid"="input.bed_uniqid"))
add_seq<-add_seq %>% select(-anno.biotype.bed_uniqid)
seq_index <- which(as.data.frame(unlist(apply(add_seq,2,function(x){unique(str_detect(x[1],"^[AGCT].*[AGCT]$")&str_length(x)>1)})))=="TRUE")


extendSeq_split<-as.data.frame(str_split_fixed(add_seq[,seq_index], "Y", 2))
colnames(extendSeq_split)<-c("extendSeq_bef","extendSeq_aft")
add_seq<-cbind(add_seq,extendSeq_split)
seq_group<-data.frame(uniq_seq=unique(add_seq$extendSeq_aft))
seq_group$seq_id<-1:length(seq_group$uniq_seq)
add_seq <- add_seq %>% left_join(seq_group,by=c("extendSeq_aft"="uniq_seq"))


# priority<-c("rRNA","tRNA","snRNA","snoRNA","miRNA","mRNA","lncRNA","pseudogene","miscRNA","ribozyme","scRNA","srpRNA","rmsk","TEC","intergenic")
priority<-c("rRNA","tRNA","sncRNA","mRNA","lncRNA","circRNA","IG_gene","TR_gene","pseudogene","repeatMasker","intergenic")
priority_table<-data.frame(type=priority,priority_rank=c(seq(length(priority))))
add_seq <- add_seq %>% left_join(priority_table,by=c("V9.x"="type"))
add_seq$group_id<-paste("group",add_seq$seq_id,add_seq$priority_rank,sep="-")
write.table(add_seq,paste(outFile,"_add_seq_group.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
write.xlsx(add_seq,paste(outFile,"_add_seq_group.xlsx",sep=""), overwrite = TRUE)


add_seq_group_list<-split(add_seq,add_seq$seq_id)

chro_id<-c("chr21",paste("chr",1:20,sep=""),"chr22","chrM","chrX","chrY")
chro_priority<-data.frame(chro_id=chro_id,chro_rank=1:length(chro_id))
invisible(lapply(seq_along(add_seq_group_list),function(x){
  add_seq_group_list_tmp <- add_seq_group_list[[x]] %>% left_join( chro_priority,by=c("V1.x"="chro_id"))
  add_seq_group_list_tmp<-arrange(add_seq_group_list_tmp,priority_rank,chro_rank)
  add_seq_group_list_tmp<-add_seq_group_list_tmp %>% select(-chro_rank)
	add_seq_group_list[[x]]<<-head(add_seq_group_list_tmp,1)
	}))

add_seq_group_uniq<-do.call(rbind.data.frame, add_seq_group_list)


if(str_detect(input.bed.file,"roc_psi_prediction.bed")){
  colnames(add_seq_group_uniq)[1:54]<-c("chrom",
  "chromStart",
  "chromEnd",
  "name",
  "foldChange",
  "strand",
  "annotation",
  "gene_feature",
  "gene_biotype",
  "tolBaseNum",
  "tqueryCov",
  "tsampCov",
  "tupDist",
  "tdownDist",
  "chrom.y",
  "chromStart.y",
  "chromEnd.y",
  "name.y",
  "foldChange.y",
  "strand.y",
  "geneName",
  "geneStart",
  "geneEnd",
  "base",
  "treatPval",
  "ctrlPval",
  "minusPval",
  "treatStopNum",
  "treatStopRPM",
  "treatPreStopNum",
  "treatAfterStopNum",
  "treatReadthroughNum",
  "ctrlStopNum",
  "ctrlStopRPM",
  "ctrlPreStopNum",
  "ctrlAfterStopNum",
  "ctrlReadthroughNum",
  "stopRpmFC",
  "treatPreRpmFold",
  "ctrlPreRpmFold",
  "preRpmFoldRatio",
  "treatAfterRpmFold",
  "ctrlAfterRpmFold",
  "afterRpmFoldRatio",
  "treatStopRatio",
  "ctrlStopRatio",
  "stopRatioFC",
  "treatStopMeanNum",
  "treatStopMeanFold",
  "ctrlStopMeanNum",
  "ctrlStopMeanFold",
  "treatStopMeanFoldRatio",
  "extendSeq",
  "class")
add_seq_group_uniq$base<-"T"
add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
}else if(str_detect(input.bed.file,"svm_psi_prediction.bed")){
  colnames(add_seq_group_uniq)[1:59]<-c("chrom",
  "chromStart",
  "chromEnd",
  "name",
  "foldChange",
  "strand",
  "annotation",
  "gene_feature",
  "gene_biotype",
  "tolBaseNum",
  "tqueryCov",
  "tsampCov",
  "tupDist",
  "tdownDist",
  "chrom.y",
  "chromStart.y",
  "chromEnd.y",
  "name.y",
  "foldChange.y",
  "strand.y",
  "geneName",
  "geneStart",
  "geneEnd",
  "base",
  "treatPval",
  "ctrlPval",
  "minusPval",
  "treatStopNum",
  "treatStopRPM",
  "treatPreStopNum",
  "treatAfterStopNum",
  "treatReadthroughNum",
  "ctrlStopNum",
  "ctrlStopRPM",
  "ctrlPreStopNum",
  "ctrlAfterStopNum",
  "ctrlReadthroughNum",
  "stopRpmFC",
  "treatPreRpmFold",
  "ctrlPreRpmFold",
  "preRpmFoldRatio",
  "treatAfterRpmFold",
  "ctrlAfterRpmFold",
  "afterRpmFoldRatio",
  "treatStopRatio",
  "ctrlStopRatio",
  "stopRatioFC",
  "treatStopMeanNum",
  "treatStopMeanFold",
  "ctrlStopMeanNum",
  "ctrlStopMeanFold",
  "treatStopMeanFoldRatio",
  "extendSeq",
  "pred.decision.values",
  "svm_pred_desc_class",
  "psi_prob",
  "non_psi_prob",
  "pred",
  "svm_pred_prob_class")
add_seq_group_uniq$base<-"T"
add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")

}else if(str_detect(input.bed.file,"ann_psi_prediction.bed")){
  colnames(add_seq_group_uniq)[1:56]<-c("chrom",
  "chromStart",
  "chromEnd",
  "name",
  "foldChange",
  "strand",
  "annotation",
  "gene_feature",
  "gene_biotype",
  "tolBaseNum",
  "tqueryCov",
  "tsampCov",
  "tupDist",
  "tdownDist",
  "chrom.y",
  "chromStart.y",
  "chromEnd.y",
  "name.y",
  "foldChange.y",
  "strand.y",
  "geneName",
  "geneStart",
  "geneEnd",
  "base",
  "treatPval",
  "ctrlPval",
  "minusPval",
  "treatStopNum",
  "treatStopRPM",
  "treatPreStopNum",
  "treatAfterStopNum",
  "treatReadthroughNum",
  "ctrlStopNum",
  "ctrlStopRPM",
  "ctrlPreStopNum",
  "ctrlAfterStopNum",
  "ctrlReadthroughNum",
  "stopRpmFC",
  "treatPreRpmFold",
  "ctrlPreRpmFold",
  "preRpmFoldRatio",
  "treatAfterRpmFold",
  "ctrlAfterRpmFold",
  "afterRpmFoldRatio",
  "treatStopRatio",
  "ctrlStopRatio",
  "stopRatioFC",
  "treatStopMeanNum",
  "treatStopMeanFold",
  "ctrlStopMeanNum",
  "ctrlStopMeanFold",
  "treatStopMeanFoldRatio",
  "extendSeq",
  "psi_prob",
  "non_psi_prob",
  "pred")
add_seq_group_uniq$base<-"T"
add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")

}else{
  colnames(add_seq_group_uniq)[1:53]<-c("chrom",
  "chromStart",
  "chromEnd",
  "name",
  "foldChange",
  "strand",
  "annotation",
  "gene_feature",
  "gene_biotype",
  "tolBaseNum",
  "tqueryCov",
  "tsampCov",
  "tupDist",
  "tdownDist",
  "chrom.y",
  "chromStart.y",
  "chromEnd.y",
  "name.y",
  "foldChange.y",
  "strand.y",
  "geneName",
  "geneStart",
  "geneEnd",
  "base",
  "treatPval",
  "ctrlPval",
  "minusPval",
  "treatStopNum",
  "treatStopRPM",
  "treatPreStopNum",
  "treatAfterStopNum",
  "treatReadthroughNum",
  "ctrlStopNum",
  "ctrlStopRPM",
  "ctrlPreStopNum",
  "ctrlAfterStopNum",
  "ctrlReadthroughNum",
  "stopRpmFC",
  "treatPreRpmFold",
  "ctrlPreRpmFold",
  "preRpmFoldRatio",
  "treatAfterRpmFold",
  "ctrlAfterRpmFold",
  "afterRpmFoldRatio",
  "treatStopRatio",
  "ctrlStopRatio",
  "stopRatioFC",
  "treatStopMeanNum",
  "treatStopMeanFold",
  "ctrlStopMeanNum",
  "ctrlStopMeanFold",
  "treatStopMeanFoldRatio",
  "extendSeq")
add_seq_group_uniq$base<-"T"
add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")
}

#get extended fasta sequence
add_seq_group_uniq$index<-paste(add_seq_group_uniq$chrom,":",as.numeric(add_seq_group_uniq$chromStart-20),"-",as.numeric(add_seq_group_uniq$chromEnd+20),sep="")
chrM_probe<-which(add_seq_group_uniq$chrom=="chrM")
add_seq_group_uniq[chrM_probe,]$index<-paste(add_seq_group_uniq$chrom[chrM_probe],":",as.numeric(add_seq_group_uniq$chromStart[chrM_probe]-3),"-",as.numeric(add_seq_group_uniq$chromEnd[chrM_probe]+3),sep="")
sort_index<-bedr.sort.region(add_seq_group_uniq$index)
arrange_index<-match(sort_index,add_seq_group_uniq$index)
add_seq_group_uniq<-add_seq_group_uniq[arrange_index,]
get_fasta<-as.data.frame(get.fasta(add_seq_group_uniq$index,fasta=genome_fasta))
add_seq_group_uniq$extendSeq_20nt<-get_fasta$sequence
add_seq_group_uniq$extendSeq_20nt<-toupper(add_seq_group_uniq$extendSeq_20nt)

seq_rev <- function(char) {
  alphabets <- strsplit(char, split = "")[[1]]
  return(rev(alphabets))
}
seq_compl <- function(seq) {
  # Check if there's "T" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  cmplvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A", "N" = "N")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A", "N" = "N")
    }
  })
  return(paste(cmplvec, collapse = ""))
}

revcom <- function(input) {
  # Make sure the input is character and in upper case
  input <- as.character(input)
  input <- toupper(input)
  # Use regular expression to check if there's only legal bases
  # present in the sequence
  legal_char <- Reduce(`&`, grepl("^[A,T,C,G,U,N]*$", input))
  if (!legal_char) {
    stop("revcom() only applies to DNA/RNA sequences, and only A/T/C/G/U/N is allowed")
  }
  rev <- seq_rev(input)
  return(seq_compl(rev))
}

minus_strand<-which(add_seq_group_uniq$strand=="-")
invisible(lapply(minus_strand,function(x){
  add_seq_group_uniq$extendSeq_20nt[x]<<-revcom(add_seq_group_uniq$extendSeq_20nt[x])
}
))

if(dim(add_seq_group_uniq)[1]>0){
  pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq$gene_biotype)),desc(Freq))
  write.table(pseudoU_anno_genetype_num,paste(outFile,"_anno_gene_biotype_num.sort",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
  percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
  pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
  pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
  pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
      theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank()) + 
        theme(panel.grid = element_blank(), 
        panel.background=element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
        guides(fill = guide_legend(title = "gene biotype"))+
        scale_fill_manual(values = mycol)

  pdf(paste(outFile,"_gene_biotype_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())

  pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq$gene_feature)),desc(Freq))
  percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
  percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
  pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
  pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
  pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
      theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank()) + 
        theme(panel.grid = element_blank(), 
        panel.background=element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
        guides(fill = guide_legend(title = "gene feature"))+
        scale_fill_manual(values = mycol)

  pdf(paste(outFile,"_gene_feature_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())


  write.table(add_seq_group_uniq,paste(outFile,"_add_seq_group_uniq.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(add_seq_group_uniq,paste(outFile,"_add_seq_group_uniq.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  write.xlsx(add_seq_group_uniq,paste(outFile,"_add_seq_group_uniq.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(add_seq_group_uniq,paste(outFile,"_anno_group_redundance.xlsx",sep=""), overwrite = TRUE)

  #read human.hg38.Pseudo.result.col29.xlsx
  # cat("Detecting novel psi (human.hg38.Pseudo.result.col29.xlsx: all known pseudouridylation sites download from RMBase)\n")
  # human.hg38.Pseudo.result.col29.xlsx<-read.xlsx(human.hg38.Pseudo.result.col29.xlsx.file)#"human.hg38.Pseudo.result.col29.xlsx"
  # human.hg38.Pseudo.result.col29.xlsx.tab<-as.data.frame(table(add_seq_group_uniq$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id))
  # colnames(human.hg38.Pseudo.result.col29.xlsx.tab)<-c("known","hit")
  # print(human.hg38.Pseudo.result.col29.xlsx.tab, row.names = F)
  # novel_psi<-add_seq_group_uniq[add_seq_group_uniq$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id=="FALSE",]
  # common_psi<-add_seq_group_uniq[add_seq_group_uniq$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id=="TRUE",]
  # common_psi<-common_psi %>% left_join(human.hg38.Pseudo.result.col29.xlsx,by=c("uniq_id"="uniq_id"))
  cat("Detecting novel psi (human.hg38.Pseudo.result.col29.xlsx: all known pseudouridylation sites download from RMBase)\n")
  human.hg38.Pseudo.result.col29.xlsx<-read.xlsx(human.hg38.Pseudo.result.col29.xlsx.file)#"human.hg38.Pseudo.result.col29.xlsx"
  # human.hg38.Pseudo.result.col29.xlsx<-distinct(human.hg38.Pseudo.result.col29.xlsx, Seq, .keep_all = TRUE)
  human.hg38.Pseudo.result.col29.xlsx.tab<-as.data.frame(table(add_seq_group_uniq$extendSeq_20nt%in%human.hg38.Pseudo.result.col29.xlsx$Seq))
  colnames(human.hg38.Pseudo.result.col29.xlsx.tab)<-c("known","hit")
  print(human.hg38.Pseudo.result.col29.xlsx.tab, row.names = F)
  novel_psi<-add_seq_group_uniq[add_seq_group_uniq$extendSeq_20nt%in%human.hg38.Pseudo.result.col29.xlsx$Seq=="FALSE",]
  common_psi<-add_seq_group_uniq[add_seq_group_uniq$extendSeq_20nt%in%human.hg38.Pseudo.result.col29.xlsx$Seq=="TRUE",]
  # common_psi<-common_psi %>% left_join(human.hg38.Pseudo.result.col29.xlsx,by=c("extendSeq_20nt"="Seq"))
  probe<-match(common_psi$extendSeq_20nt,human.hg38.Pseudo.result.col29.xlsx$Seq)
  common_psi<-cbind(common_psi,human.hg38.Pseudo.result.col29.xlsx[probe,])
  write.xlsx(novel_psi,paste(outFile,"human.hg38.Pseudo.result.col29_novel.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(common_psi,paste(outFile,"human.hg38.Pseudo.result.col29_common.xlsx",sep=""), overwrite = TRUE)

  known_num<-arrange(as.data.frame(table(add_seq_group_uniq$extendSeq_20nt%in%human.hg38.Pseudo.result.col29.xlsx$Seq)),desc(Freq))
  known_num$Var1<-str_replace(known_num$Var1,"TRUE","Known")
  known_num$Var1<-str_replace(known_num$Var1,"FALSE","Novel")
  percent<-round(100*known_num$Freq/sum(known_num$Freq),2)
  percent <-paste('(',percent, "%", ", ", known_num$Freq,')', sep = "")
  known_num$Var1 <- paste(known_num$Var1, percent, sep = '')
  known_num$Var1 <- factor(known_num$Var1, levels = known_num$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
  pie_plot<-ggplot(data = known_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
      theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank()) + 
        theme(panel.grid = element_blank(), 
        panel.background=element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
        guides(fill = guide_legend(title = "Total Known/Novel psi"))+
        scale_fill_manual(values = mycol)

  pdf(paste(outFile,"_known_novel_psi_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())
}


#read hg38.psiU.SingleSites.bed/snoRNAbase.fa
add_seq_group_uniq_rRNA<-add_seq_group_uniq [add_seq_group_uniq$gene_biotype=="rRNA" ,]
if(dim(add_seq_group_uniq_rRNA)[1]>0){
  cat("Detecting novel rRNA psi (hg38.psiU.SingleSites.bed: all known pseudouridylation sites in rRNA)\n")
  hg38.psiU.SingleSites.bed<-read.table(hg38.psiU.SingleSites.bed.file,head=F)#"hg38.psiU.SingleSites.bed"
  colnames(hg38.psiU.SingleSites.bed)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  # hg38.psiU.SingleSites.bed$uniq_id<-paste(hg38.psiU.SingleSites.bed$chrom,hg38.psiU.SingleSites.bed$chromStart,hg38.psiU.SingleSites.bed$chromEnd,hg38.psiU.SingleSites.bed$strand,sep="_")
  hg38.psiU.SingleSites.bed$index<-paste(hg38.psiU.SingleSites.bed$chrom,":",as.numeric(hg38.psiU.SingleSites.bed$chromStart-20),"-",as.numeric(hg38.psiU.SingleSites.bed$chromEnd+20),sep="")
  sort_index<-bedr.sort.region(hg38.psiU.SingleSites.bed$index)
  arrange_index<-match(sort_index,hg38.psiU.SingleSites.bed$index)
  hg38.psiU.SingleSites.bed<-hg38.psiU.SingleSites.bed[arrange_index,]
  get_fasta<-as.data.frame(get.fasta(hg38.psiU.SingleSites.bed$index,fasta=genome_fasta))
  hg38.psiU.SingleSites.bed$extendSeq_20nt<-get_fasta$sequence
  hg38.psiU.SingleSites.bed$extendSeq_20nt<-toupper(hg38.psiU.SingleSites.bed$extendSeq_20nt)
  minus_strand<-which(hg38.psiU.SingleSites.bed$strand=="-")
  invisible(lapply(minus_strand,function(x){
    hg38.psiU.SingleSites.bed$extendSeq_20nt[x]<<-revcom(hg38.psiU.SingleSites.bed$extendSeq_20nt[x])
  }
  ))


  hg38.psiU.SingleSites.bed.tab<-as.data.frame(table(add_seq_group_uniq_rRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt))
  colnames(hg38.psiU.SingleSites.bed.tab)<-c("known","hit")
  print(hg38.psiU.SingleSites.bed.tab, row.names = F)
  
  novel_rRNA_psi<-add_seq_group_uniq_rRNA[add_seq_group_uniq_rRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt=="FALSE",]
  if(dim(novel_rRNA_psi)[1]==0){
    novel_rRNA_psi[nrow(novel_rRNA_psi) + 1,]<-rep(NA,length(colnames(novel_rRNA_psi)))
    novel_rRNA_psi$label<-"novel"
  }else{
    novel_rRNA_psi$label<-"novel"
  }

  common_rRNA_psi<-add_seq_group_uniq_rRNA[add_seq_group_uniq_rRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt=="TRUE",]
  if(dim(common_rRNA_psi)[1]==0){
    common_rRNA_psi[nrow(common_rRNA_psi) + 1,]<-rep(NA,length(colnames(common_rRNA_psi)))
    common_rRNA_psi$label<-"known"
  }else{
    common_rRNA_psi$label<-"known"
  }
  # common_rRNA_psi<-common_rRNA_psi %>% left_join(hg38.psiU.SingleSites.bed,by=c("extendSeq_20nt"="extendSeq_20nt"))
  # write.xlsx(novel_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_novel_rRNA.xlsx",sep=""), overwrite = TRUE)
  # write.xlsx(common_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_common_rRNA.xlsx",sep=""), overwrite = TRUE)

  snoRNAbase<-readLines(snoRNA_fasta_rRNA)
  snoRNAbase<-paste(snoRNAbase,collapse="")
  # snoRNAbase<-str_extract_all(snoRNAbase,">[0-9\\_]*S[A-Z]*")

  snoRNAbase_hit<-list()
  if(!is.na(novel_rRNA_psi$chrom)){
    invisible(lapply(seq_along(novel_rRNA_psi$extendSeq),function(x){
    extendSeq_reg<-str_replace_all(novel_rRNA_psi$extendSeq[x],"T","[YT]")
    if(grepl(extendSeq_reg,snoRNAbase)){
        snoRNAbase_hit[x]<<-x
      }else{
        snoRNAbase_hit[x]<<-NA
      }
    }))
  }else{
    print("No novel rRNA psi-sites")
  }
  

  snoRNAbase_hit_TRUE<-novel_rRNA_psi[which(!is.na(snoRNAbase_hit)),]
  if(!dim(snoRNAbase_hit_TRUE)[1]>0){
    print("No rRNA psi-sites overlap with snoRNAbase rRNA dataset")
    snoRNAbase_hit_TRUE<-NULL
    snoRNAbase_hit_FALSE<-NULL
  }else{
    snoRNAbase_hit_TRUE$label<-"known"
    snoRNAbase_hit_FALSE<-novel_rRNA_psi[which(is.na(snoRNAbase_hit)),]
    snoRNAbase_hit_FALSE$label<-"novel"
    novel_rRNA_psi<-snoRNAbase_hit_FALSE
  }

  common_rRNA_psi<-rbind(common_rRNA_psi,snoRNAbase_hit_TRUE)
  total_rRNA_psi<-rbind(novel_rRNA_psi,common_rRNA_psi)
  # total_rRNA_psi<-total_rRNA_psi %>% left_join(hg38.psiU.SingleSites.bed,by=c("extendSeq_20nt"="extendSeq_20nt"))
  probe<-match(total_rRNA_psi$extendSeq_20nt,hg38.psiU.SingleSites.bed$extendSeq_20nt)
  total_rRNA_psi<-cbind(total_rRNA_psi,hg38.psiU.SingleSites.bed[probe,])
  # total_rRNA_psi<-distinct(total_rRNA_psi,extendSeq_20nt,.keep_all = TRUE)
  write.xlsx(novel_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_novel_rRNA.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(common_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_common_rRNA.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(total_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_total_rRNA.xlsx",sep=""), overwrite = TRUE)
  # total_rRNA_psi<-distinct(total_rRNA_psi,extendSeq_20nt,.keep_all = TRUE)
  known_rRNA_num<-arrange(as.data.frame(table(total_rRNA_psi$label)),desc(Freq))
  percent<-round(100*known_rRNA_num$Freq/sum(known_rRNA_num$Freq),2)
  percent <-paste('(',percent, "%", ", ", known_rRNA_num$Freq,')', sep = "")
  known_rRNA_num$Var1 <- paste(known_rRNA_num$Var1, percent, sep = '')
  known_rRNA_num$Var1 <- factor(known_rRNA_num$Var1, levels = known_rRNA_num$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
  pie_plot<-ggplot(data = known_rRNA_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
      theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank()) + 
        theme(panel.grid = element_blank(), 
        panel.background=element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
        guides(fill = guide_legend(title = "rRNA Known/Novel psi"))+
        scale_fill_manual(values = mycol)

  pdf(paste(outFile,"_known_novel_rRNA_psi_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())
}

#read hg38.psiU.SingleSites.bed/snoRNAbase.fa
add_seq_group_uniq_snRNA<-add_seq_group_uniq [add_seq_group_uniq$gene_feature=="snRNA" ,]
if(dim(add_seq_group_uniq_snRNA)[1]>0){
  cat("Detecting novel snRNA psi (hg38.psiU.SingleSites.bed: all known pseudouridylation sites in snRNA)\n")
  hg38.psiU.SingleSites.bed<-read.table(hg38.psiU.SingleSites.bed.file,head=F)#"hg38.psiU.SingleSites.bed"
  colnames(hg38.psiU.SingleSites.bed)<-c("chrom","chromStart","chromEnd","snRNA_anno","score","strand")
  # hg38.psiU.SingleSites.bed$uniq_id<-paste(hg38.psiU.SingleSites.bed$chrom,hg38.psiU.SingleSites.bed$chromStart,hg38.psiU.SingleSites.bed$chromEnd,hg38.psiU.SingleSites.bed$strand,sep="_")
  hg38.psiU.SingleSites.bed$index<-paste(hg38.psiU.SingleSites.bed$chrom,":",as.numeric(hg38.psiU.SingleSites.bed$chromStart-20),"-",as.numeric(hg38.psiU.SingleSites.bed$chromEnd+20),sep="")
  sort_index<-bedr.sort.region(hg38.psiU.SingleSites.bed$index)
  arrange_index<-match(sort_index,hg38.psiU.SingleSites.bed$index)
  hg38.psiU.SingleSites.bed<-hg38.psiU.SingleSites.bed[arrange_index,]
  get_fasta<-as.data.frame(get.fasta(hg38.psiU.SingleSites.bed$index,fasta=genome_fasta))
  hg38.psiU.SingleSites.bed$extendSeq_20nt<-get_fasta$sequence
  hg38.psiU.SingleSites.bed$extendSeq_20nt<-toupper(hg38.psiU.SingleSites.bed$extendSeq_20nt)
  minus_strand<-which(hg38.psiU.SingleSites.bed$strand=="-")
  invisible(lapply(minus_strand,function(x){
    hg38.psiU.SingleSites.bed$extendSeq_20nt[x]<<-revcom(hg38.psiU.SingleSites.bed$extendSeq_20nt[x])
  }
  ))


  hg38.psiU.SingleSites.bed.tab<-as.data.frame(table(add_seq_group_uniq_snRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt))
  colnames(hg38.psiU.SingleSites.bed.tab)<-c("known","hit")
  print(hg38.psiU.SingleSites.bed.tab, row.names = F)
  
  novel_snRNA_psi<-add_seq_group_uniq_snRNA[add_seq_group_uniq_snRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt=="FALSE",]
  if(dim(novel_snRNA_psi)[1]==0){
    novel_snRNA_psi[nrow(novel_snRNA_psi) + 1,]<-rep(NA,length(colnames(novel_snRNA_psi)))
    novel_snRNA_psi$label<-"novel"
  }else{
    novel_snRNA_psi$label<-"novel"
  }

  common_snRNA_psi<-add_seq_group_uniq_snRNA[add_seq_group_uniq_snRNA$extendSeq_20nt%in%hg38.psiU.SingleSites.bed$extendSeq_20nt=="TRUE",]
  if(dim(common_snRNA_psi)[1]==0){
    common_snRNA_psi[nrow(common_snRNA_psi) + 1,]<-rep(NA,length(colnames(common_snRNA_psi)))
    common_snRNA_psi$label<-"known"
  }else{
    common_snRNA_psi$label<-"known"
  }
  
  # common_snRNA_psi<-common_snRNA_psi %>% left_join(hg38.psiU.SingleSites.bed,by=c("extendSeq_20nt"="extendSeq_20nt"))
  # write.xlsx(novel_snRNA_psi,paste(outFile,"hg38.psiU.SingleSites_novel_snRNA.xlsx",sep=""), overwrite = TRUE)
  # write.xlsx(common_snRNA_psi,paste(outFile,"hg38.psiU.SingleSites_common_snRNA.xlsx",sep=""), overwrite = TRUE)

  snoRNAbase<-readLines(snoRNA_fasta_snRNA)
  snoRNAbase<-paste(snoRNAbase,collapse="")
  # snoRNAbase<-str_extract_all(snoRNAbase,">[0-9\\_]*S[A-Z]*")

  snoRNAbase_hit<-list()
  if(!is.na(novel_snRNA_psi$chrom)){
    invisible(lapply(seq_along(novel_snRNA_psi$extendSeq),function(x){
    if(grepl(novel_snRNA_psi$extendSeq[x],snoRNAbase)){
        snoRNAbase_hit[x]<<-x
      }else{
        snoRNAbase_hit[x]<<-NA
      }
    }))
  }else{
      print("No novel snRNA psi-sites")
    }
  

  snoRNAbase_hit_TRUE<-novel_snRNA_psi[which(!is.na(snoRNAbase_hit)),]
  if(!dim(snoRNAbase_hit_TRUE)[1]>0){
    print("No snRNA psi-sites overlap with snoRNAbase snRNA dataset")
    snoRNAbase_hit_TRUE<-NULL
    snoRNAbase_hit_FALSE<-NULL
  }else{
    snoRNAbase_hit_TRUE$label<-"known"
    snoRNAbase_hit_FALSE<-novel_snRNA_psi[which(is.na(snoRNAbase_hit)),]
    snoRNAbase_hit_FALSE$label<-"novel"
    novel_snRNA_psi<-snoRNAbase_hit_FALSE
  }
  

  common_snRNA_psi<-rbind(common_snRNA_psi,snoRNAbase_hit_TRUE)
  total_snRNA_psi<-rbind(novel_snRNA_psi,common_snRNA_psi)
  # total_snRNA_psi<-total_snRNA_psi %>% left_join(hg38.psiU.SingleSites.bed,by=c("extendSeq_20nt"="extendSeq_20nt"))
  probe<-match(total_snRNA_psi$extendSeq_20nt,hg38.psiU.SingleSites.bed$extendSeq_20nt)
  total_snRNA_psi<-cbind(total_snRNA_psi,hg38.psiU.SingleSites.bed[probe,])
  # total_snRNA_psi<-distinct(total_snRNA_psi,extendSeq_20nt,.keep_all = TRUE)
  write.xlsx(novel_snRNA_psi,paste(outFile,"hg38.psiU.SingleSites_novel_snRNA.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(common_snRNA_psi,paste(outFile,"hg38.psiU.SingleSites_common_snRNA.xlsx",sep=""), overwrite = TRUE)
  write.xlsx(total_snRNA_psi,paste(outFile,"hg38.psiU.SingleSites_total_snRNA.xlsx",sep=""), overwrite = TRUE)
  # total_snRNA_psi<-distinct(total_snRNA_psi,extendSeq_20nt,.keep_all = TRUE)
  known_snRNA_num<-arrange(as.data.frame(table(total_snRNA_psi$label)),desc(Freq))
  percent<-round(100*known_snRNA_num$Freq/sum(known_snRNA_num$Freq),2)
  percent <-paste('(',percent, "%", ", ", known_snRNA_num$Freq,')', sep = "")
  known_snRNA_num$Var1 <- paste(known_snRNA_num$Var1, percent, sep = '')
  known_snRNA_num$Var1 <- factor(known_snRNA_num$Var1, levels = known_snRNA_num$Var1)
  mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
  pie_plot<-ggplot(data = known_snRNA_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + 
      theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank()) + 
        theme(panel.grid = element_blank(), 
        panel.background=element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank()) + 
        guides(fill = guide_legend(title = "snRNA Known/Novel psi"))+
        scale_fill_manual(values = mycol)

  pdf(paste(outFile,"_known_novel_snRNA_psi_piechart.pdf",sep=""))
  print(pie_plot)
  invisible(dev.off())
}