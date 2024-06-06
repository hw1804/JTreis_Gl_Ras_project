# JTreis_Gl_Ras_project
Dear viewers,

The code we provided here is not newly invented, nor we developed new software to run the code. They are only for documentation purpose to show how we analyzed our data from a less processed form using all widely-used and well-maintained public packages and software.

The software we run code is called R, and we tested the code in the following version: RStudio 2022.07.1+554 "Spotted Wakerobin" Release (7872775ebddc40635780ca1ed238934c3345c5de, 2022-07-22) for macOS Mozilla/5.0 (Macintosh; Intel Mac OS X 14_3_0) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.10 Chrome/69.0.3497.128 Safari/537.36.

Please note that some R packages that we use have less support on the latest R version 4.3.3, as some developers may have limited time and effort to update their packages, thus more errors may emerge during package installations. Different devices may have different interactions with error, and stack overflow is a recommended debugging resource that compiles a big number of errors encountered by previous users.

Package installation instructions are written in details in the code. There is the .R format ready for you to run, or the marked pdf file for quick and easy visualization or printing. Any sentence that starts with # is a comment instead of a code execution. They should appear to be green by default in R as well as in marked PDF.  The package installation processes should take less than 1 min per package if you run R using your own computer power. However, the package ìSeurat wrappersî may ask for your own GitHub tokens and monocle 3 may take a while because some of their dependencies may be out of date. Please contact the package authors and submit a ticket if you encounter difficulties, as their support is still prompt. 

The runtime of each line is usually below 15s, however we indicated the more time-consuming commands in the comments within the code.

We recommend that you run our code line by line using the run button in R. It is more controllable and easy to stop. The code is chronologically organized, however, some intermediate output graphs may be overridden to the final graph if you run all code all at once.

The R code deposited here is mainly for scRNA-Seq downstream analysis since it is less commonly written compared to the other code. The other analyses include upstream RNA-Seq, scRNA-Seq QC, Dam-ID are very much streamlined, documented in the method section and used very popular and published software with well documented manuals and are well tested. The downstream Dam-ID analysis is mainly done in the excel using conventional tools such as filters, ìvlookupî function, and conditional formatting. No excel macro was used. The data that generated the graph are in the supplementary files. We also have deposited our raw data and intermediate data in NCBI GEO GSE256221.

We currently deposited all these code to GitHub for public view to finalize our publication.

Please contact us if you encounter any problems or have new inquiries. We will try our best to get back to you quickly. 



 



