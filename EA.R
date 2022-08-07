
#install.packages('devtools')
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
#install.packages("TCGAbiolinks")
library(TCGAbiolinks)

'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")'

#BiocManager::install("EDASeq")
library(EDASeq)



df = read.csv("C:/Users/USER/OneDrive - North Dakota University System/COBRE/project_COBRE/data/genes_long_format/sure/27k_All.csv", header=TRUE)
print(length(df$Gene_Symbol))
Genelist_two <- unique(df$Gene_Symbol)
print(length(Genelist_two))
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist_two)
#print(ansEA)
#system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist_two))
par(mar=c(1,1,1,1))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        
                        GOBPTab = ansEA$ResBP,
                        
                        GOCCTab = ansEA$ResCC,
                        
                        GOMFTab = ansEA$ResMF,
                        
                        PathTab = ansEA$ResPat,
                        
                        nRGTab = Genelist_two,
                        
                        nBar = 10,
                        
                        filename=NULL)

