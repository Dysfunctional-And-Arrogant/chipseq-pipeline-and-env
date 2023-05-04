library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library("dplyr")

#导入数据，list.files第一个参数直接用当前目录下的文件夹名，不需要./
samplefiles <- list.files("narrowPeak", pattern= ".narrowPeak", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c ("ECC-1","GM12878","MCF-7","MESSA","SK-N-SH","MESSA_subtract")
#导入hg38的knownGene数据库
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#assign peaks to genes
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)

#画图
plotAnnoPie(peakAnnoList)
plotAnnoPie(peakAnnoList$`ECC-1`)
plotAnnoPie(peakAnnoList$GM12878)
plotAnnoPie(peakAnnoList$`MCF-7`)
plotAnnoPie(peakAnnoList$MESSA)
plotAnnoPie(peakAnnoList$`SK-N-SH`)
plotAnnoPie(peakAnnoList$MESSA_subtract)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

#map EntrezIDs to gene symbols
ECC1_annot <- data.frame(peakAnnoList[["ECC-1"]]@anno)
GM12878_annot <- data.frame(peakAnnoList[["GM12878"]]@anno)
MCF7_annot <- data.frame(peakAnnoList[["MCF-7"]]@anno)
MESSA_annot <- data.frame(peakAnnoList[["MESSA"]]@anno)
SKNSH_annot <- data.frame(peakAnnoList[["SK-N-SH"]]@anno)
subtract_annot <- data.frame(peakAnnoList[["MESSA_subtract"]]@anno)

entrez_ECC1 <- ECC1_annot$geneId
entrez_GM12878 <- GM12878_annot$geneId
entrez_MCF7 <- MCF7_annot$geneId
entrez_MESSA <- MESSA_annot$geneId
entrez_SKNSH <- SKNSH_annot$geneId
entrez_subtract <- subtract_annot$geneId

annotations_edb_ECC1 <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_ECC1, columns = c("GENENAME"), keytype = "ENTREZID")
annotations_edb_GM12878 <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_GM12878, columns = c("GENENAME"), keytype = "ENTREZID")
annotations_edb_MCF7 <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_MCF7, columns = c("GENENAME"), keytype = "ENTREZID")
annotations_edb_MESSA <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_MESSA, columns = c("GENENAME"), keytype = "ENTREZID")
annotations_edb_SKNSH <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_SKNSH, columns = c("GENENAME"), keytype = "ENTREZID")
annotations_edb_subtract <- AnnotationDbi::select(EnsDb.Hsapiens.v75, keys = entrez_subtract, columns = c("GENENAME"), keytype = "ENTREZID")

annotations_edb_ECC1$ENTREZID <- as.character(annotations_edb_ECC1$ENTREZID)
annotations_edb_GM12878$ENTREZID <- as.character(annotations_edb_GM12878$ENTREZID)
annotations_edb_MCF7$ENTREZID <- as.character(annotations_edb_MCF7$ENTREZID)
annotations_edb_MESSA$ENTREZID <- as.character(annotations_edb_MESSA$ENTREZID)
annotations_edb_SKNSH$ENTREZID <- as.character(annotations_edb_SKNSH$ENTREZID)
annotations_edb_subtract$ENTREZID <- as.character(annotations_edb_subtract$ENTREZID)

#在进行下一步之前要创建peakAnnotationResults文件夹和对应的txt文件
ECC1_annot %>% 
  left_join(annotations_edb_ECC1, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/ECC1_peak_annotation.txt", sep="\t", quote=F, row.names=F)

GM12878_annot %>% 
  left_join(annotations_edb_GM12878, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/GM12878_peak_annotation.txt", sep="\t", quote=F, row.names=F)

MCF7_annot %>% 
  left_join(annotations_edb_MCF7, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/MCF7_peak_annotation.txt", sep="\t", quote=F, row.names=F)

MESSA_annot %>% 
  left_join(annotations_edb_MESSA, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/MESSA_peak_annotation.txt", sep="\t", quote=F, row.names=F)

SKNSH_annot %>% 
  left_join(annotations_edb_SKNSH, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/SKNSH_peak_annotation.txt", sep="\t", quote=F, row.names=F)

subtract_annot %>% 
  left_join(annotations_edb_subtract, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="peakAnnotationResults/subtract_peak_annotation.txt", sep="\t", quote=F, row.names=F)

#将基因名字输出到csv文件中，以便绘制Venn图
write.table(annotations_edb_ECC1, file = "peakAnnotationResults/ECC1_geneName.csv", sep = ",", quote = F, row.names = F)
write.table(annotations_edb_GM12878, file = "peakAnnotationResults/GM12878_geneName.csv", sep = ",", quote = F, row.names = F)
write.table(annotations_edb_MCF7, file = "peakAnnotationResults/MCF7_geneName.csv", sep = ",", quote = F, row.names = F)
write.table(annotations_edb_MESSA, file = "peakAnnotationResults/MESSA_geneName.csv", sep = ",", quote = F, row.names = F)
write.table(annotations_edb_SKNSH, file = "peakAnnotationResults/SKNSH_geneName.csv", sep = ",", quote = F, row.names = F)
write.table(annotations_edb_SKNSH, file = "peakAnnotationResults/subtract_geneName.csv", sep = ",", quote = F, row.names = F)

#GO富集分析,以MESSA为例，其余相同
ego_MESSA <- enrichGO(gene = entrez_MESSA, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_MESSA <- data.frame(ego_MESSA)
dotplot(ego_MESSA, showCategory = 10)

ego_subtract <- enrichGO(gene = entrez_subtract, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_subtract <- data.frame(ego_subtract)
dotplot(ego_subtract, showCategory = 10)

ego_ECC1 <- enrichGO(gene = entrez_ECC1, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_ECC1 <- data.frame(ego_ECC1)
dotplot(ego_ECC1, showCategory = 10)

ego_GM12878 <- enrichGO(gene = entrez_GM12878, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_GM12878 <- data.frame(ego_GM12878)
dotplot(ego_GM12878, showCategory = 10)

ego_MCF7 <- enrichGO(gene = entrez_MCF7, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_MCF7 <- data.frame(ego_MCF7)
dotplot(ego_MCF7, showCategory = 10)

ego_SKNSH <- enrichGO(gene = entrez_SKNSH, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
cluster_summary_SKNSH <- data.frame(ego_SKNSH)
dotplot(ego_SKNSH, showCategory = 10)

genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
compGO <- compareCluster(geneClusters = genes, fun = "enrichGO", pAdjustMethod = "BH", qvalueCutoff = 0.05, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
dotplot(compGO, showCategory = 3, title = "GO Pathway Enrichment Analysis")

#KEGG通路富集分析
ekegg_MESSA <- enrichKEGG(gene = entrez_MESSA,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
dotplot(ekegg_MESSA, showCategory = 6)

ekegg_MESSA_subtract <- enrichKEGG(gene = entrez_subtract, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_MESSA_subtract, showCategory = 6)

ekegg_ECC1 <- enrichKEGG(gene = entrez_ECC1, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_ECC1, showCategory = 6)

ekegg_GM12878 <- enrichKEGG(gene = entrez_GM12878, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_GM12878, showCategory = 6)

ekegg_SKNSH <- enrichKEGG(gene = entrez_SKNSH, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_SKNSH, showCategory = 6)

ekegg_MCF7 <- enrichKEGG(gene = entrez_MCF7, organism = 'hsa', pvalueCutoff = 0.05)
dotplot(ekegg_MCF7, showCategory = 6)

compKEGG <- compareCluster(geneClusters = genes, fun = 'enrichKEGG', organism = 'hsa', pvalueCutoff = 0.05)
dotplot(compKEGG, showCategory = 4, title = "KEGG Pathway Enrichment Analysis")
