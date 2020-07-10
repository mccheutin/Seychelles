## Part I: DADA2 process----
#### Preparation of your sequence----

#First, you need to load the original R1 R2 reads in the **/dir_fastq_source** folder. Then, download the [SILVA database references](https://zenodo.org/record/1172783#.XpSB25MzZQI) into the **/dir_ref_db** folder. Concerning the primers **515F-Y** & **926R** used for this publication, their sequences are written right here. You can also write a .txt files for each primer you used and load it in the **/dir_primers** folder to call them when you need it.

#r primers description
primer_515FY <- "GTGYCAGCMGCCGCGGTAA"
primer_926R <- "CCGYCAATTYMTTTRAGTTT"


#We need to read the folders which contain our sequences,primers sequences, environmental dataset and references.
library(stringr)
nms_seq_runs  <- list.files(dir_fastq_source)
paths_seq_runs <- list.files(dir_fastq_source, full.names = TRUE) %>% setNames(nms_seq_runs)

nms_refdb   <- list.files(dir_refdb)
paths_refdb <- list.files(dir_refdb, full.names = TRUE) %>% setNames(nms_refdb)

nms_primers   <- list.files(dir_primers)
paths_primers <- list.files(dir_primers, full.names = TRUE) %>% setNames(nms_primers)


#You have R1 and R2 reads with the corresponding sequencing names which can be very ugly. First, get the list of your samples and exctract the sample names: 
fns <- sort(list.files(dir_sample_selection, full.names = TRUE))
fns <- fns[str_detect(basename(fns), ".fastq")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]
if(length(fns_R1) != length(fns_R2)) stop("Forward and reverse files do not match.")
sample_names <- sapply(strsplit(basename(fns_R1), "_"), `[`, 1)
sample_names


#Once sample names cut, you can remove the primer sequence from your data by "counting" the number of nucleotides of each primer. Because nucleotides maybe me flowing, it not advised to remove the sequence of the primer by its nucleotide composition.
library(readr)
primer_set_fwd <- read_lines(paste0(dir_primers, "primer_", "515F-Y" , ".txt"))
primer_set_rev <- read_lines(paste0(dir_primers, "primer_", "926R", ".txt"))
primer_length_fwd <- str_length(primer_set_fwd[2])
primer_length_rev <- str_length(primer_set_rev[1])

# You can also use the sequence of the primers description.


#### Quality profiles and filter/trim the sequences -----
#Now you can plot the quality profiles of your reads (which is usely better on the R1) in order to truncate the sequences before quality drastically decrease. Then, you can filter your sequences and trim the primer length. 
dir_quality_plots <- paste0(dir_seq_processing, "/quality_pdf/plots/")
dir_create(dir_quality_plots , recursive = T)
qual_R1 <- plotQualityProfile(fns_R1[1])
qual_R2 <- plotQualityProfile(fns_R2[1])

ggsave(qual_R1, file = paste0(dir_quality_plots , "quality_profile_R1.png"))
ggsave(qual_R2, file = paste0(dir_quality_plots, "quality_profile_R2.png"))

filt_R1 <- str_c(dir_filtered, sample_names, "_R1_filt.fastq")
filt_R2 <- str_c(dir_filtered, sample_names, "_R2_filt.fastq")
names(filt_R1) <- sample_names
names(filt_R2) <- sample_names
set.seed(1000)

out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(primer_length_fwd,primer_length_rev) , rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

sample_names <- sapply(strsplit(basename(filt_R1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample_namesR <- sapply(strsplit(basename(filt_R2), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample_names, sample_namesR)) stop("Forward and reverse files do not match.")
names(filt_R1) <- sample_names
names(filt_R2) <- sample_names
set.seed(1000)

#Finally, you can learn the error rates for both reads and save them. 
### Learn the error rate ----
errF <- learnErrors(filt_R1,nbases=1e8, multithread=TRUE)
errR <- learnErrors(filt_R2, nbases=1e8, multithread=TRUE)
plotErrors_F <- plotErrors(errF, nominalQ=TRUE)
plotErrors_R <- plotErrors(errR, nominalQ=TRUE)
ggsave(plotErrors_F, file = paste0(dir_quality_plots , "plotErrors_F.png"))
ggsave(plotErrors_R, file = paste0(dir_quality_plots, "plotErrors_R.png"))


# Ps:  If you have big data, please proceed directly at the [Part Ib: DADA2 process for BIG DATA](# Part Ib: DADA2 process for BIG DATA).

#### Dereplication and merging sequences-----
#Now, we dereplicate the filtered reads in order to infer them in dada file to merge.
#Dereplication step
derepFs <- derepFastq(filt_R1, verbose=TRUE)
derepRs <- derepFastq(filt_R2, verbose=TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names
# Inference step
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = T)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = T)
dadaFs[[1]]
dadaRs[[1]]
# Merging step
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#Now, we can construct our sequence table "seqtab" and remove chimeras. It is recommandable to track the reads through the pipeline process
seqtab <- makeSequenceTable(mergers)
row.names(seqtab) = sample_names
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, paste0(dir_seq_processing,"seqtab.nochim.rds"))

#### Track reads through the pipeline ----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_names
head(track)
save(track, file = "track_515F-Y_926R.RData")
#### ------------------------------------------------------------------------------------------------------------------------------------------------------
### Part Ib: DADA2 process for BIG DATA
#These chunks are highly recommendable if you have a huge mount of sequences !

mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names
for(sam in sample_names) {
  cat("Processing:", sam, "\n")
  derepFs <- derepFastq(filt_R1[[sam]])
  ddF <- dada(derepFs, err=errF, multithread=TRUE)
  derepRs <- derepFastq(filt_R2[[sam]])
  ddR <- dada(derepRs, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepFs, ddR, derepRs)
  mergers[[sam]] <- merger
}
rm(derepFs); rm(derepRs)

#Now the construction of the sequence table 
seqtab <- makeSequenceTable(mergers)
row.names(seqtab) = sample_names
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(dir_seq_processing,"seqtab.rds"))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "tabled")
rownames(track) <- sample_names
head(track)
save(track, file = "table_track_A515F-Y_926R.RData")


#Remove chimeras -----
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, paste0(dir_seq_processing,"seqtab.nochim.rds"))

########-----------------------------------------------------------------------------------------------------------------------------------------
#### Taxonomy assignment-----

#You have already load the SILVA ref files, then you have to assign from them the sequence tables
load(dir_seq_processing , "seqtab.nochim.rds")
path_reference_db <- paste0(dir_refdb, "silva_nr_v132_train_set.fa.gz")
path_species_db   <- paste0(dir_refdb, "silva_species_assignment_v132.fa.gz")
taxaRC <- assignTaxonomy(seqtab.nochim, path_reference_db, tryRC=TRUE)
taxaSp <- addSpecies(taxaRC, path_species_db) # This step can be very long
taxa.print <- taxaSp # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
saveRDS(seqtab.nochim,taxaRC, taxaSp, file= paste0(dir_seq_processing,"dada2_files.rds"))

