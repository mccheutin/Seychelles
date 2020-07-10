## Part II: Phyloseq process ------
#Now you can create your phyloseq object with your filtered and assigned sequences. Before running the phyloseq process, you need to clean your worspace in order to free the r stack memory
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

#Set wd -----
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

#### Phyloseq object creation-----
#We can create our phyloseq object from the seqtab.nochim file created previously in the 01_dada.
load(paste0(dir_taxa_assign,"dada2_files.rds"))
load(paste0(dir_taxa_assign,"seqtab.nochim_515F-Y.926R.RData"))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),tax_table(taxaRC)) 
# You can use tax_table(taxaSp) if you need the assignment until species level


#Now, we will modify our table with ASVs names 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
save(ps, taxaRC,seqtab.nochim, seqtab, file=paste0(dir_taxa_assign ,"ps_515F-Y.926R.RData"))
# First results 
ntaxa(ps)
nsamples(ps)
sample_names(ps)[1:5]
rank_names(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:6]

#### Presentation of the envdata -----
envdata <- read.csv2(paste0(dir_refdb, "env.csv"))
colnames(envdata)[1] = "ID"
datatable(envdata[-c(3,6,10, 12)], rownames = FALSE, width = "100%",
          colnames = c("ID", "Species","Species number" ,"Compartment","Order","Family","Genus","Diet abreviation","Diet High scale","Diet Low scale","Trophic position", "Mass (g)","Total Length (cm)","Gut mass (g)","Sex","Site","Reef condition","GPS latitude","GPS longitude", "Island", "Substrat"), 
          caption = htmltools::tags$caption(style = "caption-side: 
                                            bottom; text-align: left;", 
                                            "Table: ", 
                                            htmltools::em("Sample presentation.")), 
          extensions = "Buttons", 
          options = list(columnDefs = 
                           list(list(className = "dt-left", targets = 0)), 
                         dom = "Blfrtip", pageLength = 5, 
                         lengthMenu = c(5, 10, 25, 50), 
                         buttons = c("csv", "copy"), 
                         scrollX = TRUE, scrollCollapse = TRUE))

#It's just an parenthesis for ordering envdata depending of your phyloseq objects. It is a step a bit borring but necessary if your env file is not synchro with your phyloseq object, meaning that the samples names in your env file are not the same that your sample names in the phyloseq object. 
#First, we will order the names
load(paste0(dir_taxa_assign, "dada2_files.RData"))
nochim_names <- rownames(seqtab.nochim)
env_names <- as.character(envdata$ID)

#Now, the aim is to have exactly the same names in the phyloseq object and in the env table. 
lecture <- cbind(sort(nochim_names),sort(env_names))
identical(lecture[,1], lecture[,2]) 

#If it's True , you're names are the same, directly pass to the merging step; If False, you have to replace the correct names. Also, you have to check your env file with the names of your fish species, genus or family which could be incorrect. 
levels(factor(envdata$family))
levels(factor(envdata$species))

#### Merging data ------
#Now that the seqtab and the envdata are corresponding, we will merge them into phyloseq object.
load(paste0(dir_taxa_assign,"ps_515F-Y.926R.RData"))
envdata <- read.csv2(paste0(dir_refdb, "env.csv"))
colnames(envdata)[1] = "ID"
DAT <- sample_data(envdata)
DAT
sort(sample_names(DAT)) == sort(sample_names(ps)) # must be true to be merged

ps1<- merge_phyloseq(ps, DAT)
ps1
save(ps1, envdata,file=paste0(dir_taxa_assign ,"ps1_515F-Y.926R.RData"))

#Plot the sample minimum ASV -----
min(rowSums(ps1@otu_table@.Data))
readsumsdf = data.frame(nreads = sort(taxa_sums(ps1),TRUE), sorted = 1:ntaxa(ps1), type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps1), TRUE), sorted = 1:nsamples(ps1), type = "Samples"))

title = "Total number of reads"
pdf(file = paste0(dir_quality_plots, "total_number_of_reads.pdf"), he = 7, wi = 7)
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()

#### Phyloseq creation-------
#We will filter our ps1 to only keep the Prokaryotes
load(paste0(dir_taxa_assign, "ps1_515F-Y.926R.RData"))
ps_sey_proka <- subset_taxa(ps_sey , Kingdom %in% c("Archaea", "Bacteria"))
save(ps_sey_proka, file = paste0(dir_taxa_assign, "ps_sey_proka.RData"))

#Then, we will obtain the phylogenetic distances between ASVs
#The tree was obtained after aligning the sequences on mothur and was names "Sey.tree" and was load in the *\dir_taxa_assign* folder. The assignment was made against the ARB database available on SILVA database reference (v132) and proceeded on [**Mothur**](https://aem.asm.org/content/75/23/7537) platform.
sey_tree <- read.tree(paste(dir_taxa_assign, "Sey.tree"))
sey_tree2 <- root(sey_tree, "ASV40619", resolve.root = T)
sey_tree2 <- drop.tip(sey_tree2,"ASV40619" ) # To delete the outgroup ASV

library(picante)
cal1<-makeChronosCalib(sey_tree2, node = "root", age.min = 1, age.max = 1, interactive = FALSE, soft.bounds = FALSE) #calibration for ultrametric branchs.
sey_chronogramme<-chronos(sey_tree2, lambda=0, model = "discrete", cal=cal1, quiet = FALSE, control=chronos.control(nb.rate.cat=1))
save(sey_chronogramme ,file = paste0(dir_taxa_assign ,"sey_chronogramme.Rdata"))

ps_sey_tree <- merge_phyloseq(ps_sey_proka, sey_tree2)
save(ps_sey_tree, file = paste0(dir_taxa_assign, "ps_sey_tree.RData"))

 
#Now that the phyloseq is only made by Prokaryotes which phylogenetic distances are associated, we will finally filter all organels (Chloroplast and Mitochondria) and potential contaminants of the exctraction kits and steps of amplification as mentionned in [Salter et al., 2014](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z).

#### Clean the phyloseq object ------
#Filter the contaminants and organels
cont_list <- read.csv("/Users/marie-charlottecheutin/Drive/Thesis/Seychelles/WP5.3_data/data_sources/reference_databases/list_potential_contaminants.csv")
contaminant <- cont_list$Genus
sey_final <- subset_taxa(ps_sey_tree, Order != "Chloroplast")
sey_final <- subset_taxa(sey_final , Family != "Mitochondria")
sey_contaminant <- subset_taxa(sey_final , Genus %in% contaminant)
save(sey_contaminant , file = paste0(dir_taxa_assign, "sey_contaminant.Rdata"))
sey_final <- subset_taxa(sey_final , !Genus %in% contaminant)
save(sey_final, file = paste0(dir_taxa_assign,"sey_final.RData"))
seyrff_final  <- prune_samples(sample_sums(sey_final) >= min(sample_sums(sey_final)) , sey_final)
seyrff_final <- rarefy_even_depth(seyrff_final, sample.size = min(sample_sums(sey_final)))
save(seyrff_final, file = paste0(dir_taxa_assign, "seyrff_final.RData"))

#Now we have the final phyloseq for all the samples, we can subset the phyloseq object for the gut, the algae and the turf and rarefy them. To be sure that the minimal abundance of the data set is enough to capture all the diversity, we use the function **ggrare()** 
ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {

  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)

  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }

  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }

  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))

  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")

  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                       size = 4, hjust = 0)
  }

  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

And now use it and create the phyloseq objects

#### Subset the phyloseq objects------
load(paste0(dir_taxa_assign, "sey_final.RData"))
# Fish ps ----
sey_gut <- subset_samples(sey_final, type == "gut")
sey_gut <- prune_taxa(names(which(colSums(sey_gut@otu_table)>0)), sey_gut)
save(sey_gut, file = paste0(dir_taxa_assign, "sey_gut.RData"))
sort(sample_sums(sey_gut))
set.seed(10000)
p_gut <- ggrare(sey_gut, step = 500, color = "geomorpho", label = "tax1", se = FALSE)
seyrff_gut <- prune_samples(sample_sums(sey_gut) >= min(sample_sums(sey_gut)) , sey_gut)
seyrff_gut <- rarefy_even_depth(seyrff_gut, sample.size = min(sample_sums(sey_gut)))
save(seyrff_gut, file = paste0(dir_taxa_assign, "seyrff_gut.RData"))
# Algae ps -----
sey_algae <- subset_samples(sey_final, tax1 %in% c("macroalgae", "turf"))
sey_algae <- prune_taxa(names(which(colSums(sey_algae@otu_table)>0)), sey_algae)
save(sey_algae, file = paste0(dir_taxa_assign, "sey_algae.RData"))
set.seed(10000)
p_algae <- ggrare(sey_algae, step = 500, color = "geomorpho", label = "tax1", se = FALSE)
seyrff_algae <- prune_samples(sample_sums(sey_algae) >= min(sample_sums(sey_algae)) , sey_algae)
seyrff_algae <- rarefy_even_depth(seyrff_algae, sample.size = min(sample_sums(sey_algae)))
save(seyrff_algae, file = paste0(dir_taxa_assign, "seyrff_algae.RData"))
