# this script is to count number of cell lines that meet the conditions of cutoff for miRNA expression level
# 2023/10/30 made

# import CCLE miRNA expression data
# this data is located at ""
setwd("C:/Rdata/CCLE_data")
miRNA.exp <-read.table("CCLE_miRNA_20181103.gct.txt",sep="\t",header = T,stringsAsFactors = F)
miRNA.exp <-miRNA.exp[,-1]

# import CCLE miRNA list
# this list is located at ""
miRNA.list <-read.table("table_of_CCLE_miRNAs_correspond_to_miRBase_v20.txt",sep="\t",header = T,stringsAsFactors = F)
miRNA <-unique(miRNA.list[,3])

# remove virus miRNAs and dead miRNAs
m <-match(miRNA,miRNA.exp[,1])
miRNA.exp <-miRNA.exp[m,]

miRNA.exp <-miRNA.exp[,-1]

# log2 miRNA expression
miRNA.exp <-as.data.frame(apply(miRNA.exp, 2, log2))

# make empty list
all.miRNA.exp.list <-as.data.frame(matrix(nrow=637,ncol=11))
colnames(all.miRNA.exp.list) <-c("miRNA", "log2(exp)>=0", "log2(exp)>=1", "log2(exp)>=2", "log2(exp)>=3", "log2(exp)>=4",
                              "log2(exp)>=5","log2(exp)>=6","log2(exp)>=7","log2(exp)>=8","log2(exp)>=9")

# count number of cell lines that meet the conditions
for (i in 1:nrow(miRNA.exp)) {
  all.miRNA.exp.list[i,1] <-miRNA[i]
  all.miRNA.exp.list[i,2] <-length(which(miRNA.exp[i,]>0))
  all.miRNA.exp.list[i,3] <-length(which(miRNA.exp[i,]>=1))
  all.miRNA.exp.list[i,4] <-length(which(miRNA.exp[i,]>=2))
  all.miRNA.exp.list[i,5] <-length(which(miRNA.exp[i,]>=3))
  all.miRNA.exp.list[i,6] <-length(which(miRNA.exp[i,]>=4))
  all.miRNA.exp.list[i,7] <-length(which(miRNA.exp[i,]>=5))
  all.miRNA.exp.list[i,8] <-length(which(miRNA.exp[i,]>=6))
  all.miRNA.exp.list[i,9] <-length(which(miRNA.exp[i,]>=7))
  all.miRNA.exp.list[i,10] <-length(which(miRNA.exp[i,]>=8))
  all.miRNA.exp.list[i,11] <-length(which(miRNA.exp[i,]>=9))
  
}

# calculate mean miRNA expression level
miRNA.exp.mean <-as.data.frame(apply(miRNA.exp, 1, mean))
mean <-miRNA.exp.mean[,1]
all.miRNA.exp.list[,12] <-mean
colnames(all.miRNA.exp.list)[12] <-"mean expression level"

# output
setwd("C:/Rdata/20231019_CCLE_miRNA_cutoff")
write.table(all.miRNA.exp.list,"CCLE_number_of_cell_lines_meet_the_conditions.txt",sep="\t",row.names = F,quote = F)
write.table(all.miRNA.exp.list[all.miRNA.exp.list[,9]>=3,1],"list_of_CCLE_miRNA_meet_the_condition_log2_7_ÅÜ3.txt",sep="\t",row.names = F,col.names = F,quote = F)

# remain only hemato cell lines
h <-grep("HAEMATO",colnames(miRNA.exp))
hemato.miRNA.exp <-miRNA.exp[,h]

# make empty list
hemato.miRNA.exp.list <-as.data.frame(matrix(nrow=637,ncol=11))
colnames(hemato.miRNA.exp.list) <-c("miRNA", "log2(exp)>=0", "log2(exp)>=1", "log2(exp)>=2", "log2(exp)>=3", "log2(exp)>=4",
                              "log2(exp)>=5","log2(exp)>=6","log2(exp)>=7","log2(exp)>=8","log2(exp)>=9")

# count number of cell lines that meet the conditions
for (i in 1:nrow(hemato.miRNA.exp)) {
  hemato.miRNA.exp.list[i,1] <-miRNA[i]
  hemato.miRNA.exp.list[i,2] <-length(which(hemato.miRNA.exp[i,]>0))
  hemato.miRNA.exp.list[i,3] <-length(which(hemato.miRNA.exp[i,]>=1))
  hemato.miRNA.exp.list[i,4] <-length(which(hemato.miRNA.exp[i,]>=2))
  hemato.miRNA.exp.list[i,5] <-length(which(hemato.miRNA.exp[i,]>=3))
  hemato.miRNA.exp.list[i,6] <-length(which(hemato.miRNA.exp[i,]>=4))
  hemato.miRNA.exp.list[i,7] <-length(which(hemato.miRNA.exp[i,]>=5))
  hemato.miRNA.exp.list[i,8] <-length(which(hemato.miRNA.exp[i,]>=6))
  hemato.miRNA.exp.list[i,9] <-length(which(hemato.miRNA.exp[i,]>=7))
  hemato.miRNA.exp.list[i,10] <-length(which(hemato.miRNA.exp[i,]>=8))
  hemato.miRNA.exp.list[i,11] <-length(which(hemato.miRNA.exp[i,]>=9))

}

# calculate mean miRNA expression level
hemato.miRNA.exp.mean <-as.data.frame(apply(hemato.miRNA.exp, 1, mean))
hemato.mean <-hemato.miRNA.exp.mean[,1]
hemato.miRNA.exp.list[,12] <-hemato.mean
colnames(hemato.miRNA.exp.list)[12] <-"mean expression level"

# output
write.table(hemato.miRNA.exp.list,"CCLE_numer_of_hemato_cell_lines_meet_the_conditions.txt",sep="\t",row.names = F,quote = F)
write.table(hemato.miRNA.exp.list[hemato.miRNA.exp.list[,8]>=3,1],"list_of_CCLE_hemato_miRNA_meet_the_condition_log2_6_ÅÜ3.txt",sep="\t",row.names = F,col.names = F,quote = F)

# remain only liver cell lines
l <-grep("LIVER",colnames(miRNA.exp))
liver.miRNA.exp <-miRNA.exp[,l]

# make empty list
liver.miRNA.exp.list <-as.data.frame(matrix(nrow=637,ncol=11))
colnames(liver.miRNA.exp.list) <-c("miRNA", "log2(exp)>=0", "log2(exp)>=1", "log2(exp)>=2", "log2(exp)>=3", "log2(exp)>=4",
                                    "log2(exp)>=5","log2(exp)>=6","log2(exp)>=7","log2(exp)>=8","log2(exp)>=9")

# count number of cell lines that meet the conditions
for (i in 1:nrow(liver.miRNA.exp)) {
  liver.miRNA.exp.list[i,1] <-miRNA[i]
  liver.miRNA.exp.list[i,2] <-length(which(liver.miRNA.exp[i,]>0))
  liver.miRNA.exp.list[i,3] <-length(which(liver.miRNA.exp[i,]>=1))
  liver.miRNA.exp.list[i,4] <-length(which(liver.miRNA.exp[i,]>=2))
  liver.miRNA.exp.list[i,5] <-length(which(liver.miRNA.exp[i,]>=3))
  liver.miRNA.exp.list[i,6] <-length(which(liver.miRNA.exp[i,]>=4))
  liver.miRNA.exp.list[i,7] <-length(which(liver.miRNA.exp[i,]>=5))
  liver.miRNA.exp.list[i,8] <-length(which(liver.miRNA.exp[i,]>=6))
  liver.miRNA.exp.list[i,9] <-length(which(liver.miRNA.exp[i,]>=7))
  liver.miRNA.exp.list[i,10] <-length(which(liver.miRNA.exp[i,]>=8))
  liver.miRNA.exp.list[i,11] <-length(which(liver.miRNA.exp[i,]>=9))
  
}

# calculate mean miRNA expression level
liver.miRNA.exp.mean <-as.data.frame(apply(liver.miRNA.exp, 1, mean))
liver.mean <-liver.miRNA.exp.mean[,1]
liver.miRNA.exp.list[,12] <-liver.mean
colnames(liver.miRNA.exp.list)[12] <-"mean expression level"

# output
write.table(liver.miRNA.exp.list,"CCLE_numer_of_liver_cell_lines_meet_the_conditions.txt",sep="\t",row.names = F,quote = F)
write.table(liver.miRNA.exp.list[liver.miRNA.exp.list[,7]>=3,1],"list_of_CCLE_liver_miRNA_meet_the_condition_log2_5_ÅÜ3.txt",sep="\t",row.names = F,col.names = F,quote = F)
