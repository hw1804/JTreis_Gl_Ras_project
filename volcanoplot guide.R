#enhanced volcano package is an easy to use, established package for generating volcano plot in R. Installation guide is here: https://github.com/kevinblighe/EnhancedVolcano
#Install ggrepel if you haven't installed.     install.packages("ggrepel")
#install ggplot2 and readxl packages with the insructions here: https://ggplot2.tidyverse.org/

#Load these packages if you are starting a new R session or never loaded them before.
library(ggrepel)
library(EnhancedVolcano)
library(ggplot2)
library(readxl)
#load the Differential expression dataset 
volcanoallglras <- read_excel("~/Desktop/volcanoallglras.xlsx")

#Generate the volcano plot with your genes of interest with selectLab function
EnhancedVolcano(volcanoallglras,
                lab = as.character(volcanoallglras$glvsglrasvolcano),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('pph13','trp','Syt4',
                              'scrt','CG12605','lz','futsch'),
                xlab = bquote('log2 fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                title = 'glRas vs gl volcano plot',
                colConnectors = 'black')
