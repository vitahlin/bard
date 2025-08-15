######课程来源：生信私学
######零基础学生信，就到生信私学
######保姆式教学，一对一指导
######生信私学，让您出类拔萃
######关注微信公众号：biofsci 回复 G2 获得脚本
######脚本经常更新，请和我们联系更换

rm(list = ls()) #清除环境变量
options(stringsAsFactors = F)

#安装包
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(data.table,dplyr,tidyverse,GEOquery)


# 第一步：生存fasta文件
# 读取包soft文件
filename = "GPL14613_family.soft"
#biof <- fread(filename, header=T,skip = "ID")
biof <- Table(getGEO(filename=filename))

#选择两列，一列是探针名称，一列是序列
colnames(biof)

#选择ID和Sequence两列,注意单词，有的是小写字母，需要修改
gpl_seq <- biof %>% select(ID,Sequence)
colnames(gpl_seq)

#过滤掉没有序列的行
gpl_seq <- gpl_seq %>% filter(nchar(Sequence)!=0) 

#整理成fasta文件
gpl_fasta <- paste0('>',gpl_seq$ID,'\n', gpl_seq$Sequence)
head(gpl_fasta)
write.table(gpl_fasta,'GPL.fasta', quote = F, row.names = F, col.names = F)

#第二步：比对序列




#第三步 整理文件

## 读入数据
probe2ID <- data.table::fread("gene_map.txt",data.table = F)
colnames(probe2ID)

## 整理数据，留下探针和转化id
probe2ID <- probe2ID %>%
  select(probe_id,trans_id) %>% 
  separate(trans_id,into = c("Ensembl",
                             "drop1","drop2","drop3",
                             "trans_Symble","gene_Symble","drop4","trans_biotype"),sep = "\\|") %>% 
  select(probe_id,Ensembl,trans_Symble,gene_Symble,trans_biotype)

write.table(probe2ID,'GPL.txt', quote = F, row.names = F, col.names = T,sep = "\t")

#第四步 表达矩阵转化为genesymbol

mydata <- fread("GSE84839_series_matrix.txt.gz", header=T,skip = "ID_REF")

#mydata <- fread("GSE84839_series_matrix.txt.gz", header=T,skip = 69)
#mydata <- fread("GSE84839_series_matrix.txt", header=T,skip = "ID_REF")
colnames(mydata)[1] <- "probe_id"

#匹配id
newdata <- inner_join(probe2ID,mydata,by="probe_id")

colnames(newdata)



newdata_mRNA <- newdata %>% 
  filter(trans_biotype=="protein_coding") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(mean)       #多行取平均值

newdata_mRNA2 <- newdata %>% 
  filter(trans_biotype=="protein_coding") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(max)        #多行取最大值



newdata_lncRNA <- newdata %>% 
  filter(trans_biotype=="lncRNA") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(mean)       #多行取平均值




newdata_lncRNA2 <- newdata %>% 
  filter(trans_biotype=="lncRNA") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(max)        #多行取最大值



newdata_miRNA <- newdata %>% 
  filter(trans_biotype=="miRNA") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(mean)       #多行取平均值


newdata_miRNA2 <- newdata %>% 
  filter(trans_biotype=="miRNA") %>%
  select(-c(1,2,3,5)) %>%
  group_by(gene_Symble) %>%
  summarise_all(max)        #多行取最大值



write.table(newdata_lncRNA,'GSE_lncRNA.txt', quote = F, row.names = F, col.names = T,sep = "\t")
write.table(newdata_mRNA,'GSE_mRNA.txt', quote = F, row.names = F, col.names = T,sep = "\t")
write.table(newdata_miRNA,'GSE_miRNA.txt', quote = F, row.names = F, col.names = T,sep = "\t")


######课程来源：生信私学
######零基础学生信，就到生信私学
######保姆式教学，一对一指导
######生信私学，让您出类拔萃
######关注微信公众号：biofsci 回复 G2 获得脚本
######脚本经常更新，请和我们联系更换