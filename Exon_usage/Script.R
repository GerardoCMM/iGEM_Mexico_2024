setwd("../Exon_usage") # Moving to Exon_usage directory

library("DEXSeq") # Loading DEXSeq library

countFiles <- c('./ERR4099794.tabular',
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
               './ERR4099805.tabular') # Files with counts for every exon of every gene

#flattenedFile = "genome_flattened.gff"

flattenedFile = "./genome_flattened_2.gff" # This file includes all the know transcripts for the proteins previously know

sampleTable <- data.frame(condition = factor(c(rep("fruit_20d",3),rep("fruit_40d",3),rep("leaf",3),rep("panicle",3))))

sampleTable$condition <- relevel(sampleTable$condition,ref="fruit_40d") # Using fruit at 40 days stage as reference

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile ) # Loading files and experimental design

extra = c("Pn8.2617", "Pn2.84", "Pn1.1317", "Pn3.4770", "Pn16.1237",
          "Pn4.3222", "Pn2.2377", "Pn12.1813","Pn7.1626", "Pn16.1198")

genesForSubset = c(rownames(difInf40Counts), extra) # Using only differentially expressed genes and genes of importance for the pathway

dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

dxd = estimateSizeFactors( dxd ) # Normalisation

dxd = estimateDispersions( dxd ) # Dispersion estimation

plotDispEsts( dxd, cex = 1) # Plot dispersions

dxd = testForDEU( dxd ) # Test for diferential exon usage

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition") # Exon usage fold changes

dxr1 = DEXSeqResults( dxd ) # Obtaining results

dfdxr1 = as.data.frame(dxr1) # Conversion to data.frame 

table ( dxr1$padj < 0.1 ) # How many exonic regions are significant

table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) ) # How many genes are significant

numbOfGenes <- sum( perGeneQValue(dxr1) < 0.1) # Number of genes with at least one differentially used exon (False discovery rate == 0.05)

plotMA( dxr1, cex = 1, colNonSig = "cyan3") # Plotting mean expression of exonic part, log 2 fold change and significance (FDR = 0.1)



### Visualization

plotDEXSeq( dxr1, "Pn7.1626", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 1 Expression

plotDEXSeq( dxr1, "Pn7.1626", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 2 Expression with transcripts

plotDEXSeq( dxr1, "Pn7.1626", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 3 Normalized counts

plotDEXSeq( dxr1, "Pn7.1626", splicing = TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 4 Exon usage

plotDEXSeq( dxr1, "Pn7.1626", splicing = TRUE, displayTranscripts=TRUE, norCounts=FALSE,expression=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 5 Exon usage with transcripts

plotDEXSeq( dxr1, "Pn7.1626", splicing = TRUE, norCounts=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 6 Expression and exon usage

plotDEXSeq( dxr1, "Pn7.1626", splicing = TRUE, displayTranscripts=TRUE, norCounts=FALSE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ) # 7 Expression and exon usage with transcripts


### Change the name of the gene in the visialuzation block to see results for other genes
