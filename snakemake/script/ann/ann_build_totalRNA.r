#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library("optparse"))
suppressMessages(library("devtools"))
suppressMessages(library("caTools"))
suppressMessages(library("neuralnet"))
suppressMessages(library("NeuralNetTools"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("gridExtra"))
suppressMessages(library("cowplot"))
suppressMessages(library("pROC"))
suppressMessages(library("mccr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpol"))
suppressMessages(library("ggpubr"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("openxlsx"))
suppressMessages(library("reshape2"))
suppressMessages(library("factoextra"))

options(warn=-1)

option_list = list(
  make_option(c("-f", "--annfile"), type="character", default=NULL, 
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

if (is.null(opt$annfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
  print_help(opt_parser);
  stop("Please provide -f annfile, -k filtfile, -r rRNAfile, -s rRNAfile2, and -o outfile_prefix option", call.=FALSE);
}

annfile = opt$annfile
filtfile = opt$filtfile
rRNAfile = opt$rRNAfile
rRNAfile2 = opt$rRNAfile2
outfile_prefix = opt$outfile_prefix

print(annfile)
print(filtfile)
print(rRNAfile)
print(rRNAfile2)
print(outfile_prefix)


#load data
cat("\n\n=====================Load data to perform classification model (factor as input)=====================\n")
ANN_data<-read.table(annfile)
colnames(ANN_data)<-c("chrom",#1
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
# ANN_data$base<-"T"
ANN_data$base<-str_replace_all(ANN_data$base, "TRUE", "T")

#rRNA-psi-non_psi visualization
ANN_data_sel<-ANN_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
ANN_data_melt<-melt(ANN_data_sel,id.vars = "class")
ANN_data_melt$class<-str_replace(as.character(ANN_data_melt$class),"0","non_psi")
ANN_data_melt$class<-str_replace(as.character(ANN_data_melt$class),"1","psi")

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

my_comparisons <- list( c("non_psi", "psi") )

#treatPreRpmFold
treatPreRpmFold<-ANN_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
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
treatAfterRpmFold<-ANN_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
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
preRpmFoldRatio<-ANN_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_violin(trim=FALSE,alpha=0.8)+
  stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
labs(x="group", y = "log2(preRpmFoldRatio)")
preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
scale_fill_manual(values=c("#129a92","#1e95d4"))

#afterRpmFoldRatio
afterRpmFoldRatio<-ANN_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
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
treatStopRatio<-ANN_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
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

stopRatioFC<-ANN_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
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

ANN_data_sel<-ANN_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
rownames(ANN_data_sel)<-ANN_data$name
ANN_data_sel$class<-str_replace_all(ANN_data_sel$class, "0", "non_psi")
ANN_data_sel$class<-str_replace_all(ANN_data_sel$class, "1", "psi")
#Encoding the target feature as factor
ANN_data_sel$class<-as.factor(ANN_data_sel$class)
cat("total classification: ")
table(ANN_data_sel$class)

variables.of.rRNA.psi <- t(ANN_data_sel %>% select(-class))
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
ANN_data_sel$psi <- ANN_data_sel$class == "psi"
ANN_data_sel$non_psi <- ANN_data_sel$class == "non_psi"
split = sample.split(ANN_data_sel$class, SplitRatio = 0.7)
training_set = subset(ANN_data_sel, split == TRUE)
test_set = subset(ANN_data_sel, split == FALSE)

# feature scaling
training_set_mean <- apply(training_set %>% select(-class, -psi, -non_psi), 2, mean)
training_set_sd <- apply(training_set %>% select(-class, -psi, -non_psi), 2, sd)
training_set[, c(-7, -8, -9)] = scale(training_set[, c(-7, -8, -9)], center = training_set_mean, scale = training_set_sd)
test_set[, c(-7, -8, -9)] = scale(test_set[, c(-7, -8, -9)], center = training_set_mean, scale = training_set_sd)

training_set_origin = subset(ANN_data, split == TRUE)
test_set_origin = subset(ANN_data, split == FALSE)
cat("\n", "training set classification using sample.split: ")
table(training_set$class)
summary(training_set)
cat("\n", "test set classification using sample.split: ")
table(test_set$class)
summary(test_set)

#Network Aplication
# ANN_data_sel.net <- neuralnet(psi + non_psi ~
#                       treatPreRpmFold + preRpmFoldRatio + treatAfterRpmFold + afterRpmFoldRatio + treatStopRatio + stopRatioFC,
#                       data = training_set, hidden = c(12, 12), rep = 5, err.fct = "ce",
#                       linear.output = FALSE, lifesign = "minimal", stepmax = 1000000,
#                       threshold = 0.001)

ANN_data_sel.net <- neuralnet(psi + non_psi ~ treatPreRpmFold + preRpmFoldRatio + treatAfterRpmFold + afterRpmFoldRatio + treatStopRatio + stopRatioFC,data = training_set, hidden = c(12, 12), rep = 10)

#Get weights for a neural network in an organized list by extracting values from a neural network object
neuralweights(ANN_data_sel.net)

#import the function from Github
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
#plot each model
pdf(paste(outfile_prefix, '_ann_nnet_train_data_plot.pdf', sep = ""),width=16,height=8)
plot.nnet(ANN_data_sel.net)
invisible(dev.off())

#Visualize the network structure of the fully connected model
pdf(paste(outfile_prefix, '_ann_best_train_data_plot.pdf', sep = ""))
plot(ANN_data_sel.net, rep = "best")
invisible(dev.off())

#Visualize the prediction effect of the model 
mlppre <- predict(ANN_data_sel.net, test_set)
colnames(mlppre)<-c("psi","non_psi")
mlpprelab <- apply(mlppre, 1, which.max)
mlpprelab<-sub("1","psi",as.character(mlpprelab))
mlpprelab<-sub("2","non_psi",as.character(mlpprelab))
plot_confusion_matrix <- ggplot() +
geom_confmat(aes(x = test_set$class, y = mlpprelab),
                        normalize = TRUE, text.perc = TRUE) +
  labs(x = "Reference", y = "Prediction") +
  scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") + 
  theme_bw() +
  theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

pdf(paste(outfile_prefix, '_ann_best_test_confusion_matrix.pdf', sep = ""))
plot_confusion_matrix
invisible(dev.off())

# Predicting Result
ANN_data_sel.prediction <- neuralnet::compute(ANN_data_sel.net, test_set[-7:-9])
idx <- apply(ANN_data_sel.prediction$net.result, 1, which.max)
net.result <- as.data.frame(ANN_data_sel.prediction$net.result)
colnames(net.result) <- c('psi', 'non_psi')
pred <- c('psi', 'non_psi')[idx]
table(pred)

pdf(paste(outfile_prefix, '_ann_roc_test_data_plot.pdf', sep=""))
ANN_roc_test <- roc(main = "Test Data ROC", test_set$class, net.result$psi, smooth = FALSE, print.auc = TRUE, col = "#e41a1c", plot = TRUE, print.thres = "best", print.thres.best.method = "youden", levels = c("non_psi", "psi"), direction = '<', auc = T, ci = T)
invisible(dev.off())
cat("best roc threshold: ")
coords(ANN_roc_test, "best", ret = "all", transpose = TRUE)[1]

summary(ANN_data_sel.net)
str(ANN_data_sel.net)
save(ANN_data_sel.net, ANN_roc_test, training_set_mean, training_set_sd, file = paste(outfile_prefix, '_ann_model.RData', sep = "")) #"my-nn.RData"
# saveRDS(mymodel, file = paste(outfile_prefix, '_ann_model.rds', sep=""))#"my-nn.rds"

#add attr to ANN_test_data
ANN_test_data <- test_set_origin
# ANN_test_data$ANN_test_prob_class <- ifelse(net.result$psi >= coords(ANN_roc_test, "best", ret = "all", transpose = TRUE)[1], "psi", "non_psi")
ANN_test_data <- cbind(ANN_test_data, net.result, pred)
write.xlsx(ANN_test_data, paste(outfile_prefix, '_ann_test_data.xlsx', sep = ""), overwrite = TRUE)

#output model info
cat("\n\n=====================tune best model confusion matrix (for modeling)=====================\n")
tab <- table(Predicted = pred, Actual = test_set$class)
tab
cat("tune best model error rate (for modeling): ", 1 - sum(diag(tab)) / sum(tab), "\n")
cat("tune best model correct rate (for modeling): ", sum(diag(tab)) / sum(tab), "\n")

#calculate MCC
actual <- test_set$class
actual <- gsub("non_psi", "0", actual)
actual <- gsub("psi", "1", actual)
pred <- gsub("non_psi", "0", pred)
pred <- gsub("psi", "1", pred)
ANN_MCC <- mccr(actual,pred)

#calculate evaluation indicators
cat("\n\n=====================Calculate evaluation indicators=====================\n")
confusion_matrix <- as.data.frame(tab)
ANN_TP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "psi",]$Freq #True Positives (TP) 
ANN_FP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "non_psi",]$Freq #False Positives (FP)
ANN_TN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "non_psi",]$Freq #True Negatives (TN)
ANN_FN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "psi",]$Freq #False Negatives (FN)
ANN_TPR <- ANN_TP / (ANN_TP + ANN_FN) #sensitivity (true positive rate, TPR)
ANN_TNR <- ANN_TN / (ANN_TN + ANN_FP) #specifity (selectivity or true negative rate, TNR)
ANN_FPR <- 1 - ANN_TNR #False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
ANN_FNR <- 1 - ANN_TPR #False Negative Rate, FNR)
ANN_Prec <- ANN_TP / (ANN_TP + ANN_FP) #Precision
ANN_Recall <- ANN_TP / (ANN_TP + ANN_FN) #Recall
ANN_ACC <- (ANN_TP + ANN_TN) / (ANN_TP + ANN_TN + ANN_FP + ANN_FN) #accuracy
ANN_F1_score <- (2 * ANN_Recall * ANN_Prec) / (ANN_Recall + ANN_Prec) #F1_score
# ANN_MCC <- (ANN_TP * ANN_TN - ANN_FP * ANN_FN) / sqrt((ANN_TP + ANN_FP) * (ANN_TP + ANN_FN) * (ANN_TN + ANN_FP) * (ANN_TN + ANN_FN))
eval <- cbind(ANN_TP, ANN_FP, ANN_TN, ANN_FN, ANN_TPR, ANN_TNR, ANN_FPR, ANN_FNR, ANN_Prec, ANN_Recall, ANN_ACC, ANN_F1_score,ANN_MCC)
eval <- round(eval, 3)
eval
write.table(eval, paste(outfile_prefix, '_ann_eval.txt', sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

#show nn evaluation as pdf table
ANN_eval_t_df <- as.data.frame(t(as.data.frame(eval)))
colnames(ANN_eval_t_df) <- "value_or_percentage"
tt3 <- ttheme_minimal(
  core = list(bg_params = list(fill = blues9[1:5], col = NA),
            fg_params = list(fontface = 3)),
  colhead = list(fg_params = list(col = "navyblue", fontface = 4L)),
  rowhead = list(fg_params = list(col = "orange", fontface = 3L)))

pdf(paste(outfile_prefix, '_ann_evaluation.pdf', sep = ""), width = 7, height = 7) # Open a new pdf file
grid.arrange(
  tableGrob(ANN_eval_t_df, theme = tt3),
  nrow = 1)
invisible(dev.off()) # Close the file


#filt by nn best model
cat("\n\n=====================Filt by nn best model=====================\n")
cat("below is your input data ready to be predicted...\n")
# get prediction
to_pred <- read.table(filtfile)
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
to_pred$base<-str_replace_all(to_pred$base, "TRUE", "T")

to_pred_var <- to_pred %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC)
to_pred_var = scale(to_pred_var, center = training_set_mean, scale = training_set_sd)
str(to_pred_var)

# Predicting Result
ANN_data_sel.prediction <- neuralnet::compute(ANN_data_sel.net, to_pred_var)
idx <- apply(ANN_data_sel.prediction$net.result, 1, which.max)
net.result <- as.data.frame(ANN_data_sel.prediction$net.result)
colnames(net.result) <- c('psi', 'non_psi')
pred <- c('psi', 'non_psi')[idx]
table(pred)


#add attr to ANN_pred_data
ANN_pred_data <- to_pred
ANN_pred_data <- cbind(ANN_pred_data, net.result, pred)
write.xlsx(ANN_pred_data, paste(outfile_prefix, '_ann_total_prediction.xlsx', sep = ""), overwrite = TRUE)

#read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

#read hg38.psiU.SingleSites.bed
evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

#known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites
ANN_data$rRNA_uniq_id<-paste(ANN_data$chrom,ANN_data$chromStart,ANN_data$chromEnd,ANN_data$strand,sep="_")
ANN_data_evidence<-ANN_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(ANN_data_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_hit.csv",sep=""))
ANN_data_no_evidence<-evidence %>% left_join(ANN_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
write.csv(ANN_data_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_miss.csv",sep=""))
recovery<-paste(round(length(unique(ANN_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
cat("psiFinder recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

#known_data miss/hit hg38.psiU.SingleSites.bed
ANN_data$rRNA_uniq_id<-paste(ANN_data$chrom,ANN_data$chromStart,ANN_data$chromEnd,ANN_data$strand,sep="_")
ANN_data_evidence2<-ANN_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(ANN_data_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_hit.csv",sep=""))
ANN_data_no_evidence2<-evidence2 %>% left_join(ANN_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
write.csv(ANN_data_no_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_miss.csv",sep=""))
recovery<-paste(round(length(unique(ANN_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
cat("psiFinder recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

#final_pred miss/hit
final_pred<-ANN_pred_data[ANN_pred_data$pred=="psi",]
final_pred<-final_pred[final_pred$foldChange>2,]
write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
write.table(final_pred,paste(outfile_prefix, '_ann_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_hit.csv",sep=""))
final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_miss.csv",sep=""))
recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
cat("psiFinder+ANN recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")