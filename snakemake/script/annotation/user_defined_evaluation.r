# library(data.table)
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("gridExtra"))


option_list = list(
  make_option(c("-f", "--rocfile"), type="character", default=NULL, 
              help="ROC of single sites [file]", metavar="character"),
  make_option(c("-p", "--positivefile"), type="character", default=NULL, 
              help="positive sites [file]", metavar="character"),
  make_option(c("-n", "--negativefile"), type="character", default=NULL, 
              help="negative sites [file]", metavar="character"),
  make_option(c("-o", "--outfile_prefix"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rocfile)|| is.null(opt$positivefile) || is.null(opt$negativefile) || is.null(opt$outfile_prefix) ){
  print_help(opt_parser);
  stop("Please provide -f rocfile, -p positivefile -n negativefile and -o outfile_prefix option", call.=FALSE);
}

ROCfile = opt$rocfile
positive = opt$positivefile
negative = opt$negativefile
outFile_prefix = opt$outfile_prefix

print(ROCfile)
print(positive)
print(negative)
print(outFile_prefix)


# read in files
# filelist_all<-c(ROCfile,"SRR1663493_GSM1554812_WT_Input_Rep1_Homo_sapiens_OTHER_versus_SRR1663494_GSM1554813_WT_Pulldown_Rep1_Homo_sapiens_OTHER_all.bed")
filelist_pos<-c(ROCfile,positive)
filelist_neg<-c(ROCfile,negative)
# datalist_all = lapply(filelist_all, function(x)fread(x, header=F)) 
datalist_pos = lapply(filelist_pos, function(x)read.table(x, header=F)) 
datalist_neg = lapply(filelist_neg, function(x)read.table(x, header=F)) 
# names(datalist_all)<-c("rRNA","all")
names(datalist_pos)<-c("rRNA","user_defined_positive")
names(datalist_neg)<-c("rRNA","user_defined_negative")

# add chrloc_id (i.e. chr1_4529052_4529080_4529080)
# lapply(seq_along(datalist_all), function(i){datalist_all[[i]]<<-within(datalist_all[[i]], chrloc_id <- paste(V1,V2,V3,V6,sep='_'))})
lapply(seq_along(datalist_pos), function(i){datalist_pos[[i]]<<-within(datalist_pos[[i]], chrloc_id <- paste(V1,V2,V3,V6,V7,sep='_'))})
lapply(seq_along(datalist_neg), function(i){datalist_neg[[i]]<<-within(datalist_neg[[i]], chrloc_id <- paste(V1,V2,V3,V6,V7,sep='_'))})

# # get common_all chrloc_id across datalist
# common_all <- Reduce(intersect, Map("[[", datalist_all, "chrloc_id"))
# # subset the datalist elements
# subset_list_all<-lapply(datalist_all, function(x) x[x$chrloc_id %in% common_all, ])
# library(dplyr)
# subset_list_arrange_all<-lapply(subset_list_all, function(x) arrange(x,by=chrloc_id))

# get common_pos chrloc_id across datalist
common_pos <- Reduce(intersect, Map("[[", datalist_pos, "chrloc_id"))
if(is.na(common_pos[1])){

	print("no common_pos id between rRNA dataset and User-defined result")

	TFP<-as.numeric(datalist_neg$rRNA$V38)
	TFP_tab<-table(TFP,deparse.level = 0)
	TFP_tab<-data.frame(rbind(TFP_tab))
	colnames(TFP_tab)<-c("FP","TP")
	TFP_tab[2]<-0
	attach(TFP_tab)
	precision<-TP/(TP+FP)
	
	}else{
		# subset the datalist elements
		subset_list_pos<-lapply(datalist_pos, function(x) x[x$chrloc_id %in% common_pos, ])
		subset_list_arrange_pos<-lapply(subset_list_pos, function(x) arrange(x,by=chrloc_id))
		TFP<-as.numeric(subset_list_arrange_pos$rRNA$V38==subset_list_arrange_pos$user_defined_positive$V38)
		TFP_tab<-table(TFP,deparse.level = 0)
		TFP_tab<-data.frame(rbind(TFP_tab))
		if(dim(TFP_tab)[2]==1){

		colnames(TFP_tab)<-"TP"
		TFP_tab$FP<-0
		attach(TFP_tab)
		precision<-TP/(TP+FP)

	}else{
		colnames(TFP_tab)<-c("FP","TP")
		attach(TFP_tab)
		precision<-TP/(TP+FP)
	}
		
	}




# get common_neg chrloc_id across datalist
common_neg <- Reduce(intersect, Map("[[", datalist_neg, "chrloc_id"))
if(is.na(common_neg[1])){

	print("no common_neg id between rRNA dataset and User-defined result")

	TFN<-as.numeric(datalist_neg$rRNA$V38)
	TFN_tab<-table(TFN,deparse.level = 0)
	TFN_tab<-data.frame(rbind(TFN_tab))
	colnames(TFN_tab)<-c("FN","TN")
	TFN_tab[1]<-0
	attach(TFN_tab)
	recall<-TP/(TP+FN)
	accuracy<-(TP+TN)/(TP+TN+FP+FN)
	F1_score<-(2*recall*precision)/(recall+precision)

	}else{

	# subset the datalist elements
	subset_list_neg<-lapply(datalist_neg, function(x) x[x$chrloc_id %in% common_neg, ])
	subset_list_arrange_neg<-lapply(subset_list_neg, function(x) arrange(x,by=chrloc_id))
	TFN<-as.numeric(subset_list_arrange_neg$rRNA$V38==subset_list_arrange_neg$user_defined_negative$V38)
	TFN_tab<-table(TFN,deparse.level = 0)
	TFN_tab<-data.frame(rbind(TFN_tab))
	
	if(dim(TFN_tab)[2]==1){

		colnames(TFN_tab)<-"TN"
		TFN_tab$FN<-0
		attach(TFN_tab)
		recall<-TP/(TP+FN)
		accuracy<-(TP+TN)/(TP+TN+FP+FN)
		F1_score<-(2*recall*precision)/(recall+precision)

	}else{

		colnames(TFN_tab)<-c("FN","TN")
		attach(TFN_tab)
		recall<-TP/(TP+FN)
		accuracy<-(TP+TN)/(TP+TN+FP+FN)
		F1_score<-(2*recall*precision)/(recall+precision)

	}

	}

evaluation<-cbind(precision,recall,accuracy,F1_score)
# colnames(evaluation)<-c("precision","recall","accuracy","F1_score")
print("Output User-defined evaluation indicators")
evaluation

#show afterRpmFoldRatio_eval_t_df as pdf table
tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = blues9[1:4], col=NA),
            fg_params=list(fontface=3)),
  colhead=list(fg_params=list(col="navyblue", fontface=4L)),
  rowhead=list(fg_params=list(col="orange", fontface=3L)))

pdf(paste(outFile_prefix, '_user_defined_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
grid.arrange(
  tableGrob(best_eval_t_df, theme=tt3),
  nrow=1)
dev.off() # Close the file



write.table(evaluation,paste(outFile_prefix, '_user_defined_evaluation.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)