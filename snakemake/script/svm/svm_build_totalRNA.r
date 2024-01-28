#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("e1071"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("cowplot"))
suppressMessages(library("pROC"))
suppressMessages(library("mccr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("reshape2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("openxlsx"))
suppressMessages(library("caTools"))
suppressMessages(library("factoextra"))
suppressMessages(library("ggpol"))


options(warn=-1)

option_list = list(
  make_option(c("-f", "--svmfile"), type="character", default=NULL, 
              help="ROC of single sites [file]", metavar="character"),
  make_option(c("-k", "--filtfile"), type="character", default=NULL, 
              help="filted file of single sites [file]", metavar="character"),
  make_option(c("-r", "--rRNAfile"), type="character", default=NULL, 
              help="rRNA of single sites [file]", metavar="character"),
  make_option(c("-s", "--rRNAfile2"), type="character", default=NULL, 
              help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
  make_option(c("-o", "--outfile_prefix"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$svmfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
  print_help(opt_parser);
  stop("Please provide -f svmfile, -k filtfile, -r rRNAfile, -s rRNAfile2, and -o outfile_prefix option", call.=FALSE);
}

svmfile = opt$svmfile
filtfile = opt$filtfile
rRNAfile = opt$rRNAfile
rRNAfile2 = opt$rRNAfile2
outfile_prefix = opt$outfile_prefix

print(svmfile)
print(filtfile)
print(rRNAfile)
print(rRNAfile2)
print(outfile_prefix)


#load data
cat("\n\n=====================Load data to perform classification model (factor as input)=====================\n")
SVM_data<-read.table(svmfile)
colnames(SVM_data)<-c("chrom",#1
  "chromStart",#2
  "chromEnd",#3
  "name",#4
  "foldChange",#5
  "strand",#6
  "geneName",#7
  "geneStart",#8
  "geneEnd",#9
  "base",#10
  "treatPval",#11
  "ctrlPval",#12
  "minusPval",#13
  "treatStopNum",#14
  "treatStopRPM",#15
  "treatPreStopNum",#16
  "treatAfterStopNum",#17
  "treatReadthroughNum",#18
  "ctrlStopNum",#19
  "ctrlStopRPM",#20
  "ctrlPreStopNum",#21
  "ctrlAfterStopNum",#22
  "ctrlReadthroughNum",#23
  "stopRpmFC",#24
  "treatPreRpmFold",#25 *
  "ctrlPreRpmFold",#26
  "preRpmFoldRatio",#27 *
  "treatAfterRpmFold",#28 *
  "ctrlAfterRpmFold",#29
  "afterRpmFoldRatio",#30 *
  "treatStopRatio",#31 *
  "ctrlStopRatio",#32
  "stopRatioFC",#33 *
  "treatStopMeanNum",#34
  "treatStopMeanFold",#35
  "ctrlStopMeanNum",#36
  "ctrlStopMeanFold",#37
  "treatStopMeanFoldRatio",#38
  "extendSeq",#39
  "class")#40
# SVM_data$base<-"T"
SVM_data$base<-str_replace_all(SVM_data$base, "TRUE", "T")

#rRNA-psi-non-psi visualization
SVM_data_sel<-SVM_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
SVM_data_melt<-melt(SVM_data_sel,id.vars = "class")
SVM_data_melt$class<-str_replace(as.character(SVM_data_melt$class),"0","non-psi")
SVM_data_melt$class<-str_replace(as.character(SVM_data_melt$class),"1","psi")

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

my_comparisons <- list( c("non-psi", "psi") )

#treatPreRpmFold
treatPreRpmFold<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
labs(x="group", y = "log2(treatPreRpmFold)")
treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

#treatAfterRpmFold
treatAfterRpmFold<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
labs(x="group", y = "log2(treatAfterRpmFold)")
treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

#preRpmFoldRatio
preRpmFoldRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
labs(x="group", y = "log2(preRpmFoldRatio)")
preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
# theme(text=element_text(size=8), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),legend.position = "none") + font("xy.text", size = 8)+
theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

#afterRpmFoldRatio
afterRpmFoldRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
  labs(x="group", y = "log2(afterRpmFoldRatio)")
afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

#treatStopRatio
treatStopRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
  labs(x="group", y = "log2(treatStopRatio)")
treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

stopRatioFC<-SVM_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
  labs(x="group", y = "log2(stopRatioFC)")
stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))


pdf(paste(outfile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
plot_grid(
  treatPreRpmFold_bp,
  preRpmFoldRatio_bp, 
  treatAfterRpmFold_bp,
  afterRpmFoldRatio_bp, 
  treatStopRatio_bp, 
  stopRatioFC_bp,
  align = "hv",
  labels = c('A','B','C','D','E','F'),ncol=2,nrow=3)
invisible(dev.off())


SVM_data_sel<-SVM_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
rownames(SVM_data_sel)<-SVM_data$name

SVM_data_sel$class<-str_replace_all(SVM_data_sel$class, "0", "non-psi")
SVM_data_sel$class<-str_replace_all(SVM_data_sel$class, "1", "psi")

#Encoding the target feature as factor
SVM_data_sel$class<-as.factor(SVM_data_sel$class)
cat("total classification: ")
table(SVM_data_sel$class)

variables.of.rRNA.psi <- t(SVM_data_sel %>% select(-class))
pca <- prcomp(variables.of.rRNA.psi, scale = TRUE)
fviz_eig_plot<-fviz_eig(pca,addlabels = T)
ggsave(paste(outfile_prefix,"_fviz_eig.pdf",sep=""),fviz_eig_plot)

group.list <- c("treatPreRpmFold","preRpmFoldRatio","treatAfterRpmFold","afterRpmFoldRatio","treatStopRatio","stopRatioFC")
group.list <- as.factor(group.list)
fviz_pca_ind_plot<-fviz_pca_ind(pca, geom.ind = c("point"), col.ind = group.list,
             palette = c(brewer.pal(9,"Set1")[1],brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[3],brewer.pal(9,"Set1")[4],brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[7]),
             legend.title = "Groups")
ggsave(paste(outfile_prefix,"_fviz_eig.pdf",sep=""),fviz_pca_ind_plot)



# Splitting the dataset into the Training set and Test set
set.seed(123)
split = sample.split(SVM_data_sel$class, SplitRatio = 0.7)
training_set = subset(SVM_data_sel, split == TRUE)
test_set = subset(SVM_data_sel, split == FALSE)
training_set_mean<-apply(training_set %>% select(-class),2,mean)# training_set_sd<-apply(training_set %>% select(-class),2,sd)
training_set_sd<-c(2,2,4,2,0.05,10)
training_set[-7] = scale(training_set[-7],training_set_mean,training_set_sd)# test_set[-7] = scale(test_set[-7],center=training_set_mean,scale=training_set_sd)
training_set_origin = subset(SVM_data, split == TRUE)
test_set_origin = subset(SVM_data, split == FALSE)
cat("\n","training set classification using sample.split: ")
table(training_set$class)
summary(training_set)#mean equals to mymodel[["x.scale"]][["scaled:center"]]; sd equals to mymodel[["x.scale"]][["scaled:scale"]] 
cat("\n","test set classification using sample.split: ")
table(test_set$class)
summary(test_set)

#evaluate contributation of each variables
cat("\n\n=====================Evaluate contributation for each variables=====================\n")
fit2 <- svm(class ~ ., data = SVM_data_sel)
w <- t(fit2$coefs) %*% fit2$SV                 # weight vectors
w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
w <- sort(w, decreasing = T)
print(w)


cat("\n\n=====================Get best model using tune (for modeling)=====================\n")

set.seed(100)
tmodel = tune(svm,
              class~.,
              data=training_set,
              type="C-classification",
              kernel="radial",
              ranges=list(cost=10^(-1:2), 
              gamma=c(.5,1,2)),
              probability=TRUE,
              scale=FALSE
)

pdf(paste(outfile_prefix, '_tmodel_plot.pdf', sep=""))
# plot(tmodel)
# plot(tmodel, transform.x = log2, transform.y = log2)
# plot(tmodel, type = "perspective", theta = 120, phi = 45)
# plot(tmodel, type = "perspective", theta = 120, phi = 30)
# plot(tmodel, type = "contour", theta = 120, phi = 45)
plot(tmodel, type = "perspective", theta = 30, phi = -10)
invisible(dev.off())

summary(tmodel)
mymodel <- tmodel$best.model
mymodel$scaled<-as.logical(rep("TRUE",6))
mymodel[["x.scale"]][["scaled:center"]]<-training_set_mean
attr(mymodel[["x.scale"]][["scaled:center"]],"names")<-colnames(SVM_data_sel)[1:6]
mymodel[["x.scale"]][["scaled:scale"]]<-training_set_sd
attr(mymodel[["x.scale"]][["scaled:scale"]],"names")<-colnames(SVM_data_sel)[1:6]

summary(mymodel)
str(mymodel)
save(mymodel, file = paste(outfile_prefix, '_svm_model.RData', sep=""))#"my-svm.RData"
saveRDS(mymodel, file = paste(outfile_prefix, '_svm_model.rds', sep=""))#"my-svm.rds"
write.svm(mymodel,svm.file = paste(outfile_prefix, '_svm_model.svm', sep=""),scale.file = paste(outfile_prefix, '_svm_model.scale', sep=""), yscale.file = paste(outfile_prefix, '_svm_model.yscale', sep=""))
pred <- predict(mymodel,test_set,probability=TRUE, decision.values=TRUE)

#Visualize the prediction effect of the model 
plot_confusion_matrix <- ggplot() +
geom_confmat(aes(x = test_set$class, y = pred),
                        normalize = TRUE, text.perc = TRUE) +
  labs(x = "Reference", y = "Prediction") +
  scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") + 
  theme_bw() +
  theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

pdf(paste(outfile_prefix, '_svm_best_test_confusion_matrix.pdf', sep = ""))
plot_confusion_matrix
invisible(dev.off())

#get attr
pred.decision.values<-as.vector(attr(pred, "decision.values"))
pred.probabilities<-attr(pred, "probabilities")
pred.probabilities<-as.data.frame(pred.probabilities)

pdf(paste(outfile_prefix, '_svm_roc_test_data_plot.pdf', sep=""))
svm_roc_test<-roc(main="Test Data ROC",test_set$class,pred.probabilities$psi,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c("non-psi","psi"), direction='<',auc=T, ci=T)
invisible(dev.off())

#add attr to SVM_test_data
SVM_test_data<-test_set_origin
SVM_test_data$pred.decision.values<-pred.decision.values
coords(svm_roc_test, "best", ret = "all", transpose = TRUE)
SVM_test_data$svm_test_prob_class<-ifelse(pred.probabilities$psi>coords(svm_roc_test, "best", ret = "all", transpose = TRUE)[1],"psi","non-psi")
SVM_test_data<-cbind(SVM_test_data,pred.probabilities,pred)
write.xlsx(SVM_test_data,paste(outfile_prefix, '_svm_test_data.xlsx', sep=""),overwrite = TRUE)

#output model info
cat("\n\n=====================tune best model confusion matrix (for modeling)=====================\n")
tab <- table(Predicted = pred,Actual = test_set$class)
tab
cat("tune best model error rate (for modeling): ",1-sum(diag(tab))/sum(tab),"\n")
cat("tune best model correct rate (for modeling): ",sum(diag(tab))/sum(tab),"\n")

#calculate MCC
actual <- test_set$class
actual <- gsub("non-psi", "0", actual)
actual <- gsub("psi", "1", actual)
pred <- gsub("non-psi", "0", pred)
pred <- gsub("psi", "1", pred)
SVM_MCC <- mccr(actual,pred)

#calculate evaluation indicators
cat("\n\n=====================Calculate evaluation indicators=====================\n")
confusion_matrix<-as.data.frame(tab)
SVM_TP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="psi",]$Freq#True Positives (TP) 
SVM_FP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="non-psi",]$Freq#False Positives (FP)
SVM_TN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="non-psi",]$Freq#True Negatives (TN)
SVM_FN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="psi",]$Freq#False Negatives (FN)
SVM_TPR <- SVM_TP / (SVM_TP + SVM_FN)#sensitivity (true positive rate, TPR)
SVM_TNR <- SVM_TN / (SVM_TN + SVM_FP)#specifity (selectivity or true negative rate, TNR)
SVM_FPR <- 1-SVM_TNR#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
SVM_FNR <- 1-SVM_TPR#False Negative Rate, FNR)
SVM_Prec <- SVM_TP / (SVM_TP + SVM_FP)#Precision
SVM_Recall <- SVM_TP / (SVM_TP + SVM_FN)#Recall
SVM_ACC <- (SVM_TP + SVM_TN) / (SVM_TP + SVM_TN + SVM_FP + SVM_FN)#accuracy
SVM_F1_score <- (2*SVM_Recall*SVM_Prec) / (SVM_Recall + SVM_Prec)#F1_score
# SVM_MCC <- (SVM_TP * SVM_TN - SVM_FP * SVM_FN) / sqrt((SVM_TP + SVM_FP) * (SVM_TP + SVM_FN) * (SVM_TN + SVM_FP) * (SVM_TN + SVM_FN))
eval<-cbind(SVM_TP,SVM_FP,SVM_TN,SVM_FN,SVM_TPR,SVM_TNR,SVM_FPR,SVM_FNR,SVM_Prec,SVM_Recall,SVM_ACC,SVM_F1_score,SVM_MCC)
eval<-round(eval,3)
eval
write.table(eval,paste(outfile_prefix, '_svm_eval.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

#show svm evaluation as pdf table
svm_eval_t_df<-as.data.frame(t(as.data.frame(eval)))
colnames(svm_eval_t_df)<-"value_or_percentage"
tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = blues9[1:4], col=NA),
            fg_params=list(fontface=3)),
  colhead=list(fg_params=list(col="navyblue", fontface=4L)),
  rowhead=list(fg_params=list(col="orange", fontface=3L)))

pdf(paste(outfile_prefix, '_svm_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
grid.arrange(
  tableGrob(svm_eval_t_df, theme=tt3),
  nrow=1)
invisible(dev.off()) # Close the file

#filt by svm best model
cat("\n\n=====================Filt by svm best model=====================\n")
cat("below is your input data ready to be predicted...\n")
# get prediction
to_pred<-read.table(filtfile)
colnames(to_pred)<-c("chrom",#1
  "chromStart",#2
  "chromEnd",#3
  "name",#4
  "foldChange",#5
  "strand",#6
  "geneName",#7
  "geneStart",#8
  "geneEnd",#9
  "base",#10
  "treatPval",#11
  "ctrlPval",#12
  "minusPval",#13
  "treatStopNum",#14
  "treatStopRPM",#15
  "treatPreStopNum",#16
  "treatAfterStopNum",#17
  "treatReadthroughNum",#18
  "ctrlStopNum",#19
  "ctrlStopRPM",#20
  "ctrlPreStopNum",#21
  "ctrlAfterStopNum",#22
  "ctrlReadthroughNum",#23
  "stopRpmFC",#24
  "treatPreRpmFold",#25
  "ctrlPreRpmFold",#26
  "preRpmFoldRatio",#27
  "treatAfterRpmFold",#28
  "ctrlAfterRpmFold",#29
  "afterRpmFoldRatio",#30
  "treatStopRatio",#31
  "ctrlStopRatio",#32
  "stopRatioFC",#33
  "treatStopMeanNum",#34
  "treatStopMeanFold",#35
  "ctrlStopMeanNum",#36
  "ctrlStopMeanFold",#37
  "treatStopMeanFoldRatio",#38
  "extendSeq")#39
# to_pred$base<-"T"
to_pred$base<-str_replace_all(to_pred$base, "TRUE", "T")
to_pred_var<-to_pred %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC)
str(to_pred_var)

#get prediction
pred <- predict(mymodel,to_pred_var,decision.values=TRUE,probability=TRUE)

#get attr
pred.decision.values<-as.vector(attr(pred, "decision.values"))
pred.probabilities<-attr(pred, "probabilities")
pred.probabilities<-as.data.frame(pred.probabilities)

#add attr to SVM_pred_data
SVM_pred_data<-to_pred
SVM_pred_data$pred.decision.values<-pred.decision.values
SVM_pred_data$svm_pred_desc_class<-ifelse(SVM_pred_data$pred.decision.values>0,"psi","non-psi")
SVM_pred_data<-cbind(SVM_pred_data,pred.probabilities,pred)
SVM_pred_data$svm_pred_prob_class<-ifelse(SVM_pred_data$psi>=0.5,"psi","non-psi")
table(SVM_pred_data$svm_pred_prob_class)
write.xlsx(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.xlsx', sep=""),overwrite = TRUE)
write.table(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
write.table(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

#read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

#read hg38.psiU.SingleSites.bed
evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

#known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites
SVM_data$rRNA_uniq_id<-paste(SVM_data$chrom,SVM_data$chromStart,SVM_data$chromEnd,SVM_data$strand,sep="_")
SVM_data_evidence<-SVM_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(SVM_data_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv",sep=""))
SVM_data_no_evidence<-evidence %>% left_join(SVM_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
write.csv(SVM_data_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv",sep=""))
recovery<-paste(round(length(unique(SVM_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
cat("rtsSeeker recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

#known_data miss/hit hg38.psiU.SingleSites.bed
SVM_data$rRNA_uniq_id<-paste(SVM_data$chrom,SVM_data$chromStart,SVM_data$chromEnd,SVM_data$strand,sep="_")
SVM_data_evidence2<-SVM_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(SVM_data_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_hit.csv",sep=""))
SVM_data_no_evidence2<-evidence2 %>% left_join(SVM_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
write.csv(SVM_data_no_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_miss.csv",sep=""))
recovery<-paste(round(length(unique(SVM_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
cat("rtsSeeker recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

#final_pred miss/hit
final_pred<-SVM_pred_data[SVM_pred_data$svm_pred_prob_class=="psi",]
final_pred<-final_pred[final_pred$foldChange>2,]
final_pred<-final_pred %>% arrange(desc(pred.decision.values))
write.table(final_pred,paste(outfile_prefix, '_svm_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
write.table(final_pred,paste(outfile_prefix, '_svm_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_hit.csv",sep=""))
final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_miss.csv",sep=""))
recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
cat("rtsSeeker+SVM recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")