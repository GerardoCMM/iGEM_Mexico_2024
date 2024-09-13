# Analisis enzimas candidatas

# Librerias

library(bio3d)
library(ggplot2)
library(pheatmap)
library(WGCNA)
library(igraph)
library(dplyr)
library(stringr)

# Cargar archivos

sample = "6"
Thres = 0.25

setwd(paste("~/Documents/Alebricks/iGEM_Mexico_2024/Enzyme_alternatives",
            sample,"/Selection",sep = ""))

data <- read.csv(paste("E",sample,"_candidates.csv",sep=""),
                 stringsAsFactors=TRUE)

rownames(data) = data$Acceso.NCBI

data_fasta = read.fasta(paste("E",sample,
                              "_candidates_aa_aln.fasta",sep=""))

names = data_fasta$id

for(i in 1:length(names)){
  if(grepl("|",names[i], fixed=TRUE)){
    names[i] = strsplit(names[i],split="\\|")[[1]][2]
  }
}

data_fasta$id = names

# Visualizacion inicial

library(ggplot2)

ggplot(data,aes(x=Homologia)) +
  geom_histogram(aes(y=..density..),position="identity", alpha=0.5,
                 fill="blue",col="black") +
  geom_density(alpha=0.6) +
  ylab("Densidad") +
  theme_bw() 

# Identity Matrix

ident_mat = seqidentity(data_fasta)

# heatmap para ver como se ven

pheatmap(ident_mat, treeheight_row = 0, treeheight_col = 0,
         show_rownames = F, show_colnames = F)

# Clusters Hierarchical clustering

adj = ident_mat

diss=1-adj

hier_adj=hclust(as.dist(diss), method="average")

colorStatic_adj=as.character(
  cutreeStaticColor(hier_adj,
                    cutHeight=Thres, minSize=1))

plotDendroAndColors(dendro = hier_adj,
                    colors=data.frame(colorStatic_adj),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Dendrogram and module colors")


# Clusters CdHit

clusters_cdhit <- "clusters.clstr"

conn <- file(clusters_cdhit,open="r")

linn <-readLines(conn)

data$cluster_cd_hit = rep(NA,length(data$Descripccion))

for (i in 1:length(linn)){
  if(!grepl(",",linn[i], fixed=TRUE)){
    cur_cluster = strsplit(linn[i]," ")[[1]][2]
  }else{
    temp_id = str_extract(linn[i],regex(">.+\\.\\.",dotall=T))
    temp_id = substring(temp_id,2,nchar(temp_id)-3)
    if(grepl("|",temp_id, fixed=TRUE)){
      temp_id = str_extract(temp_id,regex("\\|.+\\|",dotall=T))
      temp_id = substring(temp_id,2,nchar(temp_id)-1)
    }
    data[temp_id,"cluster_cd_hit"] = cur_cluster
  }
}

close(conn)

data$cluster_cd_hit = as.factor(data$cluster_cd_hit)

color_cd_hit = data[rownames(diss),"cluster_cd_hit"]

levels(color_cd_hit) = levels(as.factor(colorStatic_adj))

names(color_cd_hit) = rownames(diss)

plotDendroAndColors(dendro = hier_adj,
                    colors=data.frame(color_cd_hit),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Dendrogram and cluster colors")


plotDendroAndColors(dendro = hier_adj,
                    colors=data.frame(colorStatic_adj,color_cd_hit),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Dendrogram, module and cluster colors")




# Network

start = c()

end = c()

edge_width = c()

edge_weight = c()

adj_4net = as.data.frame(adj)

adj_4net = adj_4net[names%in%row.names(data),names%in%row.names(data)]

for(i in 1:length(adj_4net[,1])){
  for(j in i:length(adj_4net[,1])){
    if(i!= j & adj_4net[i,j]>(1-Thres)){
      start = c(start,rownames(adj_4net[i,]))
      end = c(end,rownames(adj_4net[j,]))
      edge_width = c(edge_width,adj_4net[i,j])
      edge_weight = c(edge_weight,adj_4net[i,j])
    }
  }
}

links = cbind(start,end)

nodes = rownames(adj_4net)

color_net = colorStatic_adj[names%in%row.names(data)]

size_net = data$Homologia[match(nodes,data$Acceso.NCBI)]

edge_color = c()

for(i in 1:length(start)){
  if(color_net[which(nodes==start[i])]==color_net[which(nodes==end[i])]){
    edge_color = c(edge_color,"black")
  }else{
    edge_color = c(edge_color,"red")
  }
}

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)

V(net)$color <- color_net

V(net)$size <- (size_net^3)*15

E(net)$width <- (edge_width^3)*3

E(net)$color = adjustcolor(edge_color,alpha.f = 0.5)

color_cdhit_sel = data$Cluster[match(nodes,data$Acceso.NCBI)]

color_hom_sel = data$Top.20.homología[match(nodes,data$Acceso.NCBI)]

color_cdhit_sel[is.na(color_cdhit_sel)] = "black"

color_sel = adjustcolor(color_cdhit_sel,alpha.f = 0.5)

color_sel[color_cdhit_sel==1] = "darkorange"

color_sel[color_hom_sel==1] = "darkblue"

color_cd_hit_net = as.character(color_cd_hit[nodes])

color_cd_hit_net = as.factor(color_cd_hit_net)

#Plot Network and layout

par(mar=c(0,0,0,0),mfrow=c(1,1))

l = layout_components(net)

l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

plot(net, edge.arrow.size=0,vertex.label=NA,
     rescale=F,layout=l)

plot(net, edge.arrow.size=0,vertex.label=NA,
     rescale=F,layout=l,
     vertex.color=color_cd_hit_net)

plot(net, edge.arrow.size=0,vertex.label=NA,
     rescale=F,layout=l,
     vertex.color=color_sel)








