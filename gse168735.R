library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)

x <- getGEO(filename = "GPL29838_family.soft.gz")
y <- x@dataTable@table
colnames(y)
#> [1] "ID"     "GENEID"
ids <- y


e2s = bitr(ids$GENEID,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
#> Warning in bitr(ids$GENEID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb =
#> "org.Hs.eg.db"): 0.3% of input gene IDs are fail to map...
ids = merge(ids,e2s,by.x = "GENEID",by.y = "ENTREZID")
ids = ids[,2:3]
colnames(ids) = c("probe_id","symbol")
head(ids)
#>   probe_id   symbol
#> 1     1_at     A1BG
#> 2     2_at      A2M
#> 3     9_at     NAT1
#> 4    10_at     NAT2
#> 5    12_at SERPINA3
#> 6    13_at    AADAC

write.csv(ids, "GPL29838_probe2gene.csv", row.names = FALSE)#保存结果到本地
