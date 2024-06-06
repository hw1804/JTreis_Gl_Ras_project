#The scCustomize is good at ploting gene feature plot with the same scale on mutant and wild type
#viridis should be already installed if you did the monocle 3 trajectory inference code
#Installation and loading commands are here:
install.packages("scCustomize")
library(scCustomize)
#install.packages("viridis")
library(viridis)

#Also load seurat if your R section was brand new
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Seurat)


#First adjust the sequence of wt and gl mutant so that wt appear on the left of the future graphs 
wppgl@meta.data$condition <- factor(x = wppgl@meta.data$condition, levels = c("wt", "gl mutant"))
#Use defaultassay RNA to plot gene expression levels so it is not modified due to integration processes
DefaultAssay(wppgl)<-"RNA"

#You can subset only PRs for feature plot, or you can use the full wppgl file for a full map
#Here is how you subset only PRs in seurat
wppgl_PR <-subset(x=wppgl, idents = c("R8", "R2/5", "R3/4", "R1/6", "R7", "diff.PR", "Furrow"), invert = FALSE)

#Here is how to generate feature plot using scCustom package
#first load the following colormap if you haven't done previously 
cpal1 <-c("#CFCFCF","#4662D7FF" ,"#36AAf9FF", "#1AE4B6FF", "#72FE5EFF" ,"#C7EF34FF",
          "#FABA39FF" ,"#F66B19FF" ,"#CB2A04FF" ,"#7A0403FF")
#Then plot the graph. The text sizes are used for plot zoom in the Zoom section of the plots in R. If you don't know how to zoom, you may need to decrease text sizes.
FeaturePlot_scCustom(seurat_object =wppgl, features = "sens", split.by = "condition", colors_use = cpal1, pt.size = 0.8) & theme(axis.text = element_text(size = 20),axis.title.x.top = element_text(size = 50), axis.title = element_text(size = 20), plot.title = element_text(size = 20),legend.title = element_text(size = 25),legend.text = element_text(size = 15), text = element_text(size = 60, face = "bold"))
#We can try another gene, such as CAP, in R7
FeaturePlot_scCustom(seurat_object =wppgl, features = "CAP", split.by = "condition", colors_use = cpal1, pt.size = 0.8) & theme(axis.text = element_text(size = 20),axis.title.x.top = element_text(size = 50), axis.title = element_text(size = 20), plot.title = element_text(size = 20),legend.title = element_text(size = 25),legend.text = element_text(size = 15), text = element_text(size = 60, face = "bold"))
#You can swap CAP in the code with another gene you want to plot. For example CG34377 in all PRs
FeaturePlot_scCustom(seurat_object =wppgl, features = "CG34377", split.by = "condition", colors_use = cpal1, pt.size = 0.8) & theme(axis.text = element_text(size = 20),axis.title.x.top = element_text(size = 50), axis.title = element_text(size = 20), plot.title = element_text(size = 20),legend.title = element_text(size = 25),legend.text = element_text(size = 15), text = element_text(size = 60, face = "bold"))

#Here is how to plot vln plot. You can let label appear in desired order by specifying it in scale_x_discretes limits. You can change the gene by change the feature, or you can change the cell type by change the idents.
#If you want to plot genes in another cell, don't forget change the scale x discretes too.
VlnPlot(wppgl, features = c("CG34377"), split.by = "condition", idents = c("R8", "R2/5", "R3/4", "R1/6", "R7")) + scale_x_discrete(limits = c("R8", "R2/5", "R3/4", "R1/6", "R7")) & theme(text=element_text(size = 25, face = "bold") ,axis.text = element_text(size = 30, face = "bold"))

#Another example of gene B-H1 in only one cell cluster, R1/6
VlnPlot(wppgl, features = c("B-H1"), split.by = "condition", idents = c("R1/6")) + scale_x_discrete(limits = c("R1/6")) & theme(text=element_text(size = 25, face = "bold") ,axis.text = element_text(size = 30, face = "bold"))

