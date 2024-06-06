# use install.packages("data.table") to install data.table prior to library commands down below if you never installed on your R
# magrittr comes with tidyverse. If you ran previous code you will have tidyverse. if not, install.packages("tidyverse")

# load these packages using the library command
library(data.table)
library(magrittr)

# Check the default assay type of wppgl to make sure it is integrated not RNA, which is the raw data
DefaultAssay(wppgl)

# You can get cell counts per cluster when you extract meta data
# extract meta data
extractedmetadata <- wppgl@meta.data %>% as.data.table

# count the number of cells per unique combinations of "condition" and "integrated_snn_res.0.55". The resolution is set at 0.55 as explained previously
# in the seurat data integration and marker calling file.

# The numbers that would be displayed by the following command was used in making figure 1B. Including gl mutant and wt conditions and in each cluster
extractedmetadata[, .N, by = c("condition", "integrated_snn_res.0.55")]


