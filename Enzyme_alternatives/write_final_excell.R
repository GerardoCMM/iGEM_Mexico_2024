library(ape)

E6 <- read.csv("~/Documents/Alebricks/Ginkgo/6/Selection/E6_candidates.csv", stringsAsFactors=TRUE)

E6_fasta = read.FASTA("/home/gerardo/Documents/Alebricks/Ginkgo/6/Selection/E6_selected.fasta")

E6$DNA.Length = rep(0,length(E6$Descripccion))

for(i in 2:length(E6$Descripccion)){
  if(E6$Seleccionada[i] == 1){
    for(j in 1:length(names(E6_fasta))){
      if(grepl(E6$Acceso.NCBI[i],names(E6_fasta)[j], fixed=TRUE)){
        E6$DNA.Length[i] = length(E6_fasta[[j]])
      }
    }
  }
}

write.csv(E6,
          "~/Documents/Alebricks/Ginkgo/6/Selection/E6_selection.csv")
