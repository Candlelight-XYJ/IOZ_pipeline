options(stringsAsFactors = F)
Ep_GO <-read.csv("F:\\学习文件存放地址\\IOZ_实验室\\zhangbin\\topGO\\第二次调整参数Ep_GO.csv")

#dropna <- na.omit(Ep_GO)

filter_GO <- Ep_GO[Ep_GO$V2!="test",]

topGO_input <- filter_GO[!duplicated(filter_GO$target),]


### 一行拆成多行

library(clusterProfiler)
tianniu_GO <- read.csv("F:\\学习文件存放地址\\IOZ_实验室\\zhangbin\\topGO\\天牛GO注释gofenlei.csv",stringsAsFactors=F,header=F)

res <- data.frame()
for(i in 1:nrow(tianniu_GO)){
  print(i)
  tmp <- as.vector(tianniu_GO[i,])
  nanum<-0
  
  for(j in 1:length(tmp)){
    
    if(nchar(tmp[j])==0){
      nanum <- j
      break}
  }
  compId <- as.vector(unlist(rep(tmp[1],times=(j-2))))
  GO_tmp <- as.vector(unlist(tmp[1,2:(j-1)]))
  tmp_res <- data.frame(compId=compId,GoId=as.vector(GO_tmp))
  res <-rbind(res,tmp_res)
}


  

### GO
library(clusterProfiler)

tianniu_func <- read.csv("F:\\学习文件存放地址\\IOZ_实验室\\zhangbin\\topGO\\天牛所有基因功能注释.csv",stringsAsFactors=F)
GO_anno <- merge(res,tianniu_func,by.x = "compId",by.y = "Trinity.assembly",all.x = T)
write.csv(GO_anno,"F:\\学习文件存放地址\\IOZ_实验室\\天牛\\天牛注释大全含GO和Kegg.csv",row.names=F)

###准备差异gene-list
x <- read.csv("F:\\学习文件存放地址\\IOZ_实验室\\zhangbin\\topGO\\第二次调整参数Ep_GO.csv",stringsAsFactors=F)
gene <- x$target
gene <- gene[!duplicated(gene)]

###准备go2geneid的文件
term2gene <- data.frame(GOID=GO_anno$GoId,compId=GO_anno$compId)
###准备go2description的文件
term2name <- data.frame(GOID=GO_anno$GoId,compname=GO_anno$Seq..descriptions)
###使用enricher完成富集分析
fout <- enricher(gene,TERM2GENE = term2gene,TERM2NAME = term2name)
###导出结果
write.csv(fout,"F:\\学习文件存放地址\\IOZ_实验室\\zhangbin\\topGO\\Ep_GO富集结果.csv")
###绘制GO富集分析的结果图
barplot(fout)

### KEGG

kegg2gene <- data.frame(GOID=GO_anno$K.item,compId=GO_anno$compId)
kegg2gene <- kegg2gene[kegg2gene$GOID!="-",]
kegg2gene <- na.omit(kegg2gene)

kegg2name <- data.frame(GOID=GO_anno$K.item,compname=GO_anno$Seq..descriptions)
kegg2name <- na.omit(kegg2name)


fout <- enricher(gene,TERM2GENE = kegg2gene,TERM2NAME = kegg2name)
