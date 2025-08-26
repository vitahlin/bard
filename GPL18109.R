# 清除环境变量
rm(list = ls())
options(stringsAsFactors = F)

#安装包
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
if (!requireNamespace("writexl", quietly = TRUE))
  install.packages("writexl")
library(pacman)
library(writexl)
p_load(data.table,dplyr,tidyverse,GEOquery)


# 第一步：生存fasta文件
# 读取包soft文件

# 指定soft文件，和当前代码的文件目录一致
filename = "GPL18109_family.soft"
#biof <- fread(filename, header=T,skip = "ID")
biof <- Table(getGEO(filename=filename)) 

# 打印soft文件的全部列数据
colnames(biof)

#选择ID和SEQUENCE两列,注意单词大小写
gpl_seq <- biof %>% select(ID,SEQUENCE)
colnames(gpl_seq)

#过滤掉没有序列SEQUENCE的行
gpl_seq <- gpl_seq %>% filter(nchar(SEQUENCE)!=0) 

#整理成fasta文件，注意列名
gpl_fasta <- paste0('>',gpl_seq$ID,'\n', gpl_seq$SEQUENCE)
head(gpl_fasta)
write.table(gpl_fasta,'GPL.fasta', quote = F, row.names = F, col.names = F)



#第二步：比对序列


## 读入geo_map数据
probe2ID <- data.table::fread("gene_map.txt",data.table = F)
colnames(probe2ID)
View(probe2ID)

## 整理数据，留下探针和转化id，并且把trans_id按｜切割
probe2ID <- probe2ID %>%
  select(probe_id,trans_id) %>% 
  separate(trans_id,into = c("Ensembl",
                             "drop1","drop2","drop3",
                             "trans_Symble","gene_Symble","drop4","trans_biotype"),sep = "\\|") %>% 
  select(probe_id,Ensembl,trans_Symble,gene_Symble,trans_biotype)
# 把上面的数据输出为GPL.txt
write.table(probe2ID,'GPL.txt', quote = F, row.names = F, col.names = T,sep = "\t")

#第四步 表达矩阵转化为genesymbol

mydata <- fread("GSE53625_series_matrix.txt.gz", header=T,skip = "ID_REF")
view(mydata)
#mydata <- fread("GSE84839_series_matrix.txt.gz", header=T,skip = 69)
#mydata <- fread("GSE84839_series_matrix.txt", header=T,skip = "ID_REF")

# 把第一列的ID_RED修改命名为probe_id
colnames(mydata)[1] <- "probe_id"

# 根据probe_id进行匹配关联
newdata <- inner_join(probe2ID,mydata,by="probe_id")
view(newdata)

library(writexl)
write_xlsx(newdata, path = "newdata.xlsx")


# 这里就是最终的数据
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