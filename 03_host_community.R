## Part III: Fish Community Composition between Reefs
#Here again, we will remove all the object to clean the memory 
knitr::opts_chunk$set(eval = FALSE)
remove(list = ls())
cran_packages   <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny", "microbiome",
                     "tidyverse", "miniUI", "caret", "pls", "e1071", "ggplot2", 
                     "randomForest","entropart", "vegan", "plyr", "dplyr", "here",
                     "ggrepel", "nlme", "R.utils", "gridExtra","grid", "googledrive",
                     "googlesheets", "phangorn", "devtools", "rmarkdown", "sys",
                     "reshape2", "devtools", "PMA","structSSI","ade4", "ape",
                     "Biostrings", "igraph", "ggnetwork", "intergraph", "ips",
                     "scales", "kableExtra", "pgirmess", "treemap", "knitr","kableExtra",
                     "rstudioapi" ,"data.table","DT","pander","formatR","grDevices","svgPanZoom",
                     "RCurl","plotly","pairwiseAdonis", "stringr")
github_packages <- c("jfukuyama/phyloseqGraphTest")
bioc_packages   <- c("phyloseq", "genefilter", "impute", "dada2", "DECIPHER")
# Install CRAN packages (if not already installed)
#Some packages would be not availbale for your R version
inst <- cran_packages %in% installed.packages()
if (any(! inst)) {
  install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
# 
inst <- github_packages %in% installed.packages()
if (any(! inst)) {
  devtools::install_github(github_packages[!inst]) }

# Load libraries
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
sessionInfo()
set.seed(1000)

# Set wd------
knitr::opts_knit$set(root.dir = getwd())
path = getwd()
# This will setwd to wherever the .Rmd file is opened.
dir_sample_selection <- paste0(path,"/analyses/01_select_samples/")
dir_seq_processing   <- paste0(path,"/analyses/02_process_sequences/")
dir_taxa_assign      <- paste0(path,"/analyses/03_assign_taxonomy/")
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/")
dir_primers   <- paste0(path,"/dir_data_source/primers_sequences/")
dir_refdb   <- paste0(path,"/dir_data_source/reference_databases/")
dir_fastq_source <- paste0(path,"/dir_data_source/sequences/")

#In this is part, we will focus on our sampling dataset and analyze the distribution of the hosts between the IRs and HRs.
envdata <- read.csv2(paste0(dir_refdb, "env.csv"))
colnames(envdata)[1] = "ID"
load(paste0(dir_taxa_assign, "sey_gut.RData"))
load(paste0(dir_taxa_assign, "sey_final.RData"))

#We will first write the table of the species and their corresponding diet.
diet_sp <- as.data.frame(table(envdata$diet4, envdata$tax1))[which(as.data.frame(table(envdata$diet4, envdata$tax1))[,3]> 0),][,c(1,2)]
sp_reef <- table(envdata$tax1, envdata$geomorpho)
sp_reef <- cbind(as.data.frame(sp_reef[,1]), as.data.frame(sp_reef[,2]))
table_diet_sp <- cbind(diet_sp[,2],diet_sp[,1],sp_reef[,c(1,2)])

table_to_print <- rbind(table_diet_sp[-c(26,37),], table_diet_sp[c(26,37),])
datatable(table_to_print, rownames = F, width = "100%",
          colnames = c("Species", "Diet", "HR","IR"),
          caption = htmltools::tags$caption(style = "caption-side: 
                                            bottom; text-align: left;", 
                                            "Table: ", 
                                            htmltools::em("Sampling table of the species and diet between reefs.")), 
          extensions = "Buttons", 
          options = list(columnDefs = 
                           list(list(className = "dt-left", targets = 0)), 
                         dom = "Blfrtip", pageLength = 5, 
                         lengthMenu = c(5, 10, 25, 50), 
                         buttons = c("csv", "copy"), 
                         scrollX = TRUE, scrollCollapse = TRUE))


##### PCoA on fish community------
#The distribution of the sampling set is a first result showing the influence of the shift on the fish communities. How are they distributed?
samp_data <- envdata[envdata$type=="gut",]
fam_tab <- table(samp_data$site,samp_data$family)
# Transform to log
fam.log <- log1p(fam_tab)  # Equivalent: log(fam_tab + 1)
# Principal coordinate analysis and simple ordination plot
fam.D <- vegdist(fam.log, "bray")
res <- pcoa(fam.D)
#res$values
biplot(res, fam.log)
#round(res$values$Relative_eig[1]*100, 1) # 57.8 %
#round(res$values$Relative_eig[2]*100, 1) # 26 %

site1 <- c("C1","C2","C3","C4", "M1","M2","M3")
site2 <- c(rep("coral", 4), rep("macroalgal",3))
site_data <- cbind(site1,site2)
colnames(site_data) <- c("site", "geomorpho")
site_data <- as.data.frame(site_data)

adonis(fam.D ~ geomorpho, data = site_data)
beta_reef <- betadisper(fam.D, site_data$geomorpho)
permutest(beta_reef)

#Now, we will do the same on trophic srtucture with the diet

samp_data <- envdata[envdata$type=="gut",]
diet_tab <- table(samp_data$site,samp_data$diet4)
# Transform to log
diet.log <- log1p(diet_tab)  # Equivalent: log(diet_tab + 1)
# Principal coordinate analysis and simple ordination plot
diet.D <- vegdist(diet.log, "bray")
res <- pcoa(diet.D)
#res$values
par(mfrow=c(1,2))
biplot(res, diet.log)
percent(res$values$Relative_eig[1])
percent(res$values$Relative_eig[2])

site1 <- c("C1","C2","C3","C4", "M1","M2","M3")
site2 <- c(rep("coral", 4), rep("macroalgal",3))
site_data <- cbind(site1,site2)
colnames(site_data) <- c("site", "geomorpho")
site_data <- as.data.frame(site_data)

adonis(diet.D ~ geomorpho, data = site_data)
beta_reef <- betadisper(diet.D, site_data$geomorpho)
permutest(beta_reef)

#Because we will focus our sutdy on herbivores and invertivores, we need to know if families are equally distributed.
samp_data_inv <- samp_data[samp_data$diet4 == "Mobile invertebrate",]
fam_tab <- table(samp_data_inv$site,samp_data_inv$family)
# Transform to log
fam.log <- log1p(fam_tab)  # Equivalent: log(fam_tab + 1)
# Principal coordinate analysis and simple ordination plot
fam.D <- vegdist(fam.log, "bray")
res <- pcoa(fam.D)
#res$values
biplot(res, fam.log)
percent(res$values$Relative_eig[1])
percent(res$values$Relative_eig[2])

site1 <- c("C1","C2","C3","C4", "M1","M2","M3")
site2 <- c(rep("coral", 4), rep("macroalgal",3))
site_data <- cbind(site1,site2)
colnames(site_data) <- c("site", "geomorpho")
site_data <- as.data.frame(site_data)

adonis(fam.D ~ geomorpho, data = site_data)
beta_reef <- betadisper(fam.D, site_data$geomorpho)
permutest(beta_reef)

