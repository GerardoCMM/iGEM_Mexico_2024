# Differential expression analysis of Piper nigrum tissues

setwd("/home/gerardo/Documentos/Gera/SyntheticBiobots/2022/Metods/RNAseq") # Setting working directory

library(tximport) # Importing libraries

library(pheatmap)

library(DESeq2)

library(EnhancedVolcano)

library(limma)

library(RColorBrewer)

library(factoextra)

library(ggvenn)

library(webr)

library(dplyr)

quantdata <- c('./ERR4099794.tabular', # Reading count data output of kallisto alignment
               './ERR4099795.tabular',
               './ERR4099796.tabular',
               './ERR4099797.tabular',
               './ERR4099798.tabular',
               './ERR4099799.tabular',
               './ERR4099800.tabular',
               './ERR4099801.tabular',
               './ERR4099802.tabular',
               './ERR4099803.tabular',
               './ERR4099804.tabular',
               './ERR4099805.tabular')

txi_output <- tximport(quantdata, type = "kallisto",txOut=T)

sampleTable <- data.frame(condition = factor(c(rep("fruit_20d",3),rep("fruit_40d",3),rep("leaf",3),rep("panicle",3))))

coldata <- sampleTable

coldata$condition <- relevel(coldata$condition,ref="fruit_40d")

dds <- DESeqDataSetFromTximport(txi = txi_output,
                                colData = coldata,
                                design= ~ condition)

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

dds <- DESeq(dds)

r_crudos <- as.data.frame(results(dds))

vsd <- vst(object=dds,blind=T)

plotPCA(vsd)+geom_point(size=5)+ggtitle("PCA over samples","Counts as variables")+my_theme

##########################################################################

df <- as.data.frame(colData(dds)[,"condition"])

muestras <- SraRunInfo$LibraryName

rownames(df) <- muestras

colnames(df)<- c("Sample")

f <- factor(c(rep("fruit",6),rep("not_fruit",6)))

df$Tissue <- f

#############################################################################

f40vsf20 <- results(dds,c("condition","fruit_40d","fruit_20d"))

f40vspan <- results(dds,c("condition","fruit_40d","panicle"))

f40vsleaf <- results(dds,c("condition","fruit_40d","leaf"))

f40vsf20 <- as.data.frame(f40vsf20)

f40vspan <- as.data.frame(f40vspan)

f40vsleaf <- as.data.frame(f40vsleaf)

f40vsf20 <- f40vsf20[!is.na(f40vsf20$padj), ]

f40vspan <- f40vspan[!is.na(f40vspan$padj), ]

f40vsleaf <- f40vsleaf[!is.na(f40vsleaf$padj), ]

f40vsf20 <- f40vsf20[ , c("log2FoldChange", "padj")]

f40vspan <- f40vspan[ , c("log2FoldChange", "padj")]

f40vsleaf <- f40vsleaf[ , c("log2FoldChange", "padj")]

colnames(f40vsf20) <- paste0("f40vsf20", "_", colnames(f40vsf20))

colnames(f40vspan) <- paste0("f40vspan", "_", colnames(f40vspan))

colnames(f40vsleaf) <- paste0("f40vsleaf", "_", colnames(f40vsleaf))

fullTable <- merge(f40vsf20, f40vspan, by =0)

rownames(fullTable) <- fullTable$Row.names

fullTable <- fullTable[,2:5]

fullTable <- merge(fullTable, f40vsleaf, by =0)

difInf40 <- abs(fullTable$f40vsf20_log2FoldChange) > 1 & fullTable$f40vsf20_padj < 0.05 & 
  abs(fullTable$f40vspan_log2FoldChange) > 1 & fullTable$f40vspan_padj < 0.05 & 
  abs(fullTable$f40vsleaf_log2FoldChange) > 1 & fullTable$f40vsleaf_padj < 0.05 

difInf40Table <- fullTable[difInf40, ]

forVennDif <- data.frame(Difvsf20 = abs(fullTable$f40vsf20_log2FoldChange) > 1 & fullTable$f40vsf20_padj < 0.05,
                         DifvsPan = abs(fullTable$f40vspan_log2FoldChange) > 1 & fullTable$f40vspan_padj < 0.05,
                         DifvsLeaf = abs(fullTable$f40vsleaf_log2FoldChange) > 1 & fullTable$f40vsleaf_padj < 0.05)

vennDiagram(forVennDif,circle.col = c("red","blue","green"))

upInf40 <- fullTable$f40vsf20_log2FoldChange > 1 & fullTable$f40vsf20_padj < 0.05 & 
  fullTable$f40vspan_log2FoldChange > 1 & fullTable$f40vspan_padj < 0.05 & 
  fullTable$f40vsleaf_log2FoldChange > 1 & fullTable$f40vsleaf_padj < 0.05 

upInf40Table <- fullTable[upInf40, ]

forVennUp <- data.frame(Upvsf20 = fullTable$f40vsf20_log2FoldChange > 1 & fullTable$f40vsf20_padj < 0.05,
                        UpvsPan = fullTable$f40vspan_log2FoldChange > 1 & fullTable$f40vspan_padj < 0.05,
                        UpvsLeaf = fullTable$f40vsleaf_log2FoldChange > 1 & fullTable$f40vsleaf_padj < 0.05)

vennDiagram(forVennUp,circle.col = c("red","blue","green"))

########################################################################################################

dds2 <- DESeq(dds, test = "LRT", reduced = ~1) 

acrossGroups <- results(dds2)

acrossGroups <- acrossGroups[order(acrossGroups$padj), ]

acrossGroupsTable <- as.data.frame(acrossGroups)[(acrossGroups$padj < 0.01 & !is.na(acrossGroups$padj) & abs(acrossGroups$log2FoldChange)>1),]

sigChanges <- rownames(acrossGroups)[acrossGroups$padj < 0.01 & !is.na(acrossGroups$padj) & abs(acrossGroups$log2FoldChange)>1]

x <-rownames(assay(vsd))%in%sigChanges[1:3000]

counts <- as.data.frame(subset(assay(vsd),x))

colnames(counts)<- muestras

sampleCor <- cor(counts)

sampleDists <- as.dist(1 - cor(counts))

sampleDistMatrix <- as.matrix(sampleDists)

blueColours <- brewer.pal(9, "Blues")

colors <- colorRampPalette(rev(blueColours))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,cex=1.15,color=colors)

pheatmap(counts,treeheight_row = 0,treeheight_col = 25, show_rownames=F,
         cluster_cols=T,scale="row",annotation_col=df)

fviz_nbclust(counts, FUNcluster = kmeans, method = "wss", k.max = 10,linecolor = "black")

k <- pheatmap(counts,treeheight_row = 0,treeheight_col = 25, show_rownames=F,
              cluster_cols=T,scale="row",annotation_col=df,kmeans_k = 6)

clusterDF <- as.data.frame(factor(k$kmeans$cluster))

colnames(clusterDF) <- "Cluster"

OrderByCluster <- counts[order(clusterDF$Cluster), ]

pheatmap(OrderByCluster, scale = "row", annotation_row = clusterDF,annotation_col = df, show_rownames = FALSE,
         cluster_rows = FALSE)

#############################################################################################################

x <-rownames(assay(vsd))%in%sigChanges

counts2 <- as.data.frame(subset(assay(vsd),x))

difInf40Counts <- counts2[difInf40Table$Row.names, ]

upInf40Counts <- counts2[upInf40Table$Row.names, ]

colnames(difInf40Counts)<-colnames(counts)

difInf40Counts <- na.exclude(difInf40Counts)

sampleCor2 <- cor(difInf40Counts)

sampleDists2 <- as.dist(1 - cor(difInf40Counts))

sampleDistMatrix2 <- as.matrix(sampleDists2)

pheatmap(sampleDistMatrix2, clustering_distance_rows = sampleDists2, clustering_distance_cols = sampleDists2,cex=1.15,color=colors)

pheatmap(difInf40Counts,treeheight_row = 0,treeheight_col = 25, show_rownames=F,
         cluster_cols=T,scale="row",annotation_col=df)

fviz_nbclust(difInf40Counts, FUNcluster = kmeans, method = "wss", k.max = 10,linecolor = "black")

k2 <- pheatmap(difInf40Counts,treeheight_row = 0,treeheight_col = 25, show_rownames=F,
               cluster_cols=T,scale="row",annotation_col=df,kmeans_k =8)

clusterDF2 <- as.data.frame(factor(k2$kmeans$cluster))

colnames(clusterDF2) <- "Cluster"

OrderByCluster2 <- difInf40Counts[order(clusterDF2$Cluster), ]

pheatmap(OrderByCluster2, scale = "row", annotation_row = clusterDF2, annotation_col = df, show_rownames = FALSE,
         cluster_rows = FALSE)

genes_dif <- abs(fullTable$f40vsf20_log2FoldChange) > 1 & fullTable$f40vsf20_padj < 0.05 | 
  abs(fullTable$f40vspan_log2FoldChange) > 1 & fullTable$f40vspan_padj < 0.05 | 
  abs(fullTable$f40vsleaf_log2FoldChange) > 1 & fullTable$f40vsleaf_padj < 0.05 

genes_dif_Table <- fullTable[genes_dif, ]

genes_up <- fullTable$f40vsf20_log2FoldChange > 1 & fullTable$f40vsf20_padj < 0.05 | 
  fullTable$f40vspan_log2FoldChange > 1 & fullTable$f40vspan_padj < 0.05 | 
  fullTable$f40vsleaf_log2FoldChange > 1 & fullTable$f40vsleaf_padj < 0.05 

genes_up_Table <- fullTable[genes_up, ]

#########################################################################################################

difForVolcano <- rownames(acrossGroupsTable)[rownames(acrossGroupsTable)%in%rownames(difInf40Counts)][1:50]

upForVolcano <- rownames(acrossGroupsTable)[rownames(acrossGroupsTable)%in%upInf40Table$Row.names][1:25]

forVolcano <- as.data.frame(acrossGroups)

EnhancedVolcano(forVolcano,lab = rownames(forVolcano),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.01)

EnhancedVolcano(forVolcano,selectLab = difForVolcano, lab = rownames(forVolcano),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.01,drawConnectors = T,maxoverlapsConnectors = 50,max.overlaps = 50)

EnhancedVolcano(forVolcano,selectLab = upForVolcano, lab = rownames(forVolcano),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.01,drawConnectors = T,maxoverlapsConnectors = 50,max.overlaps = 50)

######## Not by LTR

EnhancedVolcano(r_crudos,lab = rownames(r_crudos),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.05)

EnhancedVolcano(r_crudos,selectLab = difForVolcano, lab = rownames(r_crudos),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.01,drawConnectors = T,maxoverlapsConnectors = 50,max.overlaps = 50)

EnhancedVolcano(r_crudos,selectLab = upForVolcano, lab = rownames(r_crudos),x = 'log2FoldChange',
                y = 'padj',FCcutoff = 1, pCutoff = 0.01,drawConnectors = T,maxoverlapsConnectors = 50,max.overlaps = 50)

##########################################################################################################################

venn_up <- list("Up vs fruit_20d" = fullTable$Row.names[fullTable$f40vsf20_log2FoldChange > 1 & fullTable$f40vsf20_padj < 0.05],
                "Up vs Panicle" = fullTable$Row.names[fullTable$f40vspan_log2FoldChange > 1 & fullTable$f40vspan_padj < 0.05],
                "Up vs Leaf" = fullTable$Row.names[fullTable$f40vsleaf_log2FoldChange > 1 & fullTable$f40vsleaf_padj < 0.05])
ggvenn(
  venn_up, 
  fill_color = c("#FF0000", "#00FF00", "#0000FF"),
  stroke_size = 0.5, set_name_size = 10,fill_alpha = 0.3,text_size = 7
)

venn_dif <- list("DE vs fruit_20d" = fullTable$Row.names[abs(fullTable$f40vsf20_log2FoldChange) > 1 & fullTable$f40vsf20_padj < 0.05],
                 "DE vs Panicle" = fullTable$Row.names[abs(fullTable$f40vspan_log2FoldChange) > 1 & fullTable$f40vspan_padj < 0.05],
                 "DE vs Leaf" = fullTable$Row.names[abs(fullTable$f40vsleaf_log2FoldChange) > 1 & fullTable$f40vsleaf_padj < 0.05])
ggvenn(
  venn_dif, 
  fill_color = c("#FF0000", "#00FF00", "#0000FF"),
  stroke_size = 0.5, set_name_size = 10,fill_alpha = 0.3,text_size = 7
)


######  With LRT

venn_up <- list("Up vs fruit_20d" = fullTable$Row.names[fullTable$f40vsf20_log2FoldChange > 1 & fullTable$f40vsf20_padj < 0.05 & fullTable$Row.names%in%sigChanges],
                "Up vs Panicle" = fullTable$Row.names[fullTable$f40vspan_log2FoldChange > 1 & fullTable$f40vspan_padj < 0.05 & fullTable$Row.names%in%sigChanges],
                "Up vs Leaf" = fullTable$Row.names[fullTable$f40vsleaf_log2FoldChange > 1 & fullTable$f40vsleaf_padj < 0.05 & fullTable$Row.names%in%sigChanges])
ggvenn(
  venn_up, 
  fill_color = c("#FF0000", "#00FF00", "#0000FF"),
  stroke_size = 0.5, set_name_size = 10,fill_alpha = 0.3,text_size = 7
)

venn_dif <- list("DE vs fruit_20d" = fullTable$Row.names[abs(fullTable$f40vsf20_log2FoldChange) > 1 & fullTable$f40vsf20_padj < 0.05 & fullTable$Row.names%in%sigChanges],
                 "DE vs Panicle" = fullTable$Row.names[abs(fullTable$f40vspan_log2FoldChange) > 1 & fullTable$f40vspan_padj < 0.05 & fullTable$Row.names%in%sigChanges],
                 "DE vs Leaf" = fullTable$Row.names[abs(fullTable$f40vsleaf_log2FoldChange) > 1 & fullTable$f40vsleaf_padj < 0.05 & fullTable$Row.names%in%sigChanges])
ggvenn(
  venn_dif, 
  fill_color = c("#FF0000", "#00FF00", "#0000FF"),
  stroke_size = 0.5, set_name_size = 10,fill_alpha = 0.3,text_size = 7
)

#######################################################################################################################

Genes_Anot <- read.delim("./Genes_Anot.tab")

upAnot<-Genes_Anot[which(Genes_Anot$Gene.ID%in%upInf40Table$Row.names),c(1,2,35,36,37,39,40,41)]

rownames(upAnot)<-upAnot$Gene.ID

upAnot <- upAnot[2:7]

upAnot <- merge(upAnot, acrossGroupsTable, by =0)

rownames(upAnot)<-upAnot$Row.names

upAnot <- upAnot[2:13]

write.csv(upAnot,"upAnot.csv")

clusters_3000 <- Genes_Anot[which(Genes_Anot$Gene.ID%in%rownames(counts)),c(1,2,35,36,37,39,40,41)]

rownames(clusters_3000) <- clusters_3000$Gene.ID

clusters_3000 <- clusters_3000[2:7]

clusters_3000 <- merge(clusters_3000, clusterDF, by =0)

rownames(clusters_3000) <- clusters_3000$Row.names

clusters_3000 <- clusters_3000[,2:8]

write.csv(clusters_3000,"Clusters.csv")

clusters_dif40 <- Genes_Anot[which(Genes_Anot$Gene.ID%in%difInf40Table$Row.names),c(1,2,35,36,37,39,40,41)]

rownames(clusters_dif40) <- clusters_dif40$Gene.ID

clusters_dif40 <- clusters_dif40[2:7]

clusters_dif40 <- merge(clusters_dif40, clusterDF2, by =0)

rownames(clusters_dif40) <- clusters_dif40$Row.names

clusters_dif40 <- clusters_dif40[,2:8]

write.csv(clusters_dif40,"Clusters_dif40.csv")

#########################################################################################################

Clusters_2 <- read.csv("./Clusters_2.csv", row.names=1)

Clusters_dif40_2 <- read.csv("./Clusters_dif40_2.csv", row.names=1)

Clusters_dif40_2$Freq <- rep(1,479)

PD = Clusters_dif40_2 %>% group_by(Cluster, GO.Molecular.Functions) %>% summarise(n = sum(Freq))

PD <- PD[-c(18,44,86,126,157,178),]

PieDonut(PD, aes(Cluster, GO.Molecular.Functions, count=n), title = "GO Molecular function of diferentially genes")

Clusters_2$Freq <- rep(1,2826)

PD = Clusters_2 %>% group_by(Cluster, GO.Molecular.Functions) %>% summarise(n = sum(Freq))

PD <- PD[-c(37,129,230,328,413,487),]

PieDonut(PD, aes(Cluster, GO.Molecular.Functions, count=n), title = "GO Molecular function of diferentially expressed genes")

#####################################################################################################

forty <- rownames(upInf40Counts)

forty <- as.data.frame(forty)

write.csv(forty, "upIn40.csv",row.names = F)

difforty <- rownames(difInf40Counts)

difforty <- as.data.frame(difforty)

write.csv(difforty,"difIn40.csv",row.names = F)

dif3000 <- sigChanges[1:3000]

dif3000 <- as.data.frame(dif3000)

write.csv(dif3000,"dif3000.csv",row.names = F)

#####################################################################################################

#For genes overexpressed vs panicle and leaf but nof f20

f20vspan <- results(dds,c("condition","fruit_20d","panicle"))

f20vsleaf <- results(dds,c("condition","fruit_20d","leaf"))

f20vspan <- as.data.frame(f20vspan)

f20vsleaf <- as.data.frame(f20vsleaf)

f20vspan <- f20vspan[!is.na(f20vspan$padj), ]

f20vsleaf <- f20vsleaf[!is.na(f20vsleaf$padj), ]

f20vspan <- f20vspan[ , c("log2FoldChange", "padj")]

f20vsleaf <- f20vsleaf[ , c("log2FoldChange", "padj")]

colnames(f20vspan) <- paste0("f20vspan", "_", colnames(f20vspan))

colnames(f20vsleaf) <- paste0("f20vsleaf", "_", colnames(f20vsleaf))

fullTable2 <- merge(f40vsleaf, f40vspan, by =0)

rownames(fullTable2) <- fullTable2$Row.names

fullTable2 <- fullTable2[,2:5]

fullTable2 <- merge(fullTable2, f20vsleaf, by =0)

rownames(fullTable2) <- fullTable2$Row.names

fullTable2 <- fullTable2[,2:7]

fullTable2 <- merge(fullTable2, f20vspan, by =0)

upInF <-  fullTable2$f20vspan_log2FoldChange > 1 & fullTable2$f20vspan_padj < 0.05 & 
  fullTable2$f20vsleaf_log2FoldChange > 1 & fullTable2$f20vsleaf_padj < 0.05 &
  fullTable2$f40vspan_log2FoldChange > 1 & fullTable2$f40vspan_padj < 0.05 & 
  fullTable2$f40vsleaf_log2FoldChange > 1 & fullTable2$f40vsleaf_padj < 0.05 

upInF <- fullTable2[upInF, ]

upInF<-Genes_Anot[which(Genes_Anot$Gene.ID%in%upInF$Row.names),c(1,2,35,36,37,39,40,41)]

rownames(upInF)<-upInF$Gene.ID

upInF <- upInF[2:7]

upInF <- merge(upInF, acrossGroupsTable, by =0)

rownames(upInF)<-upInF$Row.names

upInF <- upInF[2:13]

write.csv(upInF,"upInF.csv")

