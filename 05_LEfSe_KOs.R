# Process the T4F algorithm on the core bacteriome of enteric reef fish----
dir.create(paste0(dir_data_cleaning),"/Fish/T4F/LEfSe/All", recursive = T)
setwd(paste0(dir_data_cleaning),"/Fish/T4F/LEfSe/All")

otu <-as.data.frame(otu_table(sey_gut_core))
write.table(t(otu), sep = '\t', 'otu_table.txt', quote = F)

Biostrings::writeXStringSet(sey_gut_core@refseq, file = "core.fasta")

runRefBlast(path_to_otus = "core.fasta", path_to_reference_data = "../Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "core_Ref99NR", database_mode = "Ref99NR", use_force = T, num_threads = 1)

makeFunctionalPrediction(path_to_otu_table = "otu_table.txt", path_to_reference_data = "../Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "core_Ref99NR", database_mode = "Ref99NR", normalize_by_copy_number = TRUE, min_identity_to_reference = round(min(read.table("./core_Ref99NR/ref_blast.txt")[,3]),0), normalize_pathways = T)

calculateFunctionalRedundancy(path_to_otu_table = "otu_table.txt",
                              path_to_reference_data = "../Tax4Fun2_ReferenceData_v2",
                              path_to_temp_folder = "core_Ref99NR",
                              database_mode = "Ref99NR",
                              min_identity_to_reference = 0.97)


prediction_core <- as.data.frame(read_delim("core_Ref99NR/functional_prediction.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
colnames(prediction_core)[1] <- "X.OTU_IDs"
write.table(prediction_core , file = "prediction_core.txt", sep = "\t", quote = F)
prediction_core <- read.table(file = "prediction_core.txt", sep = "\t" , header = T)

prediction_core2 <- left_join(prediction_core, ko_pathway_sort, by = "X.OTU_IDs")
prediction_core2$kegg[which(is.na(prediction_core$kegg) == T)] <- "None"
# How many KO redudant?
length(which(table(prediction_core$X.OTU_IDs)>1, T))
redundant_ko <- rownames(which(table(prediction_core$X.OTU_IDs)>1, T))

prediction_core2 <- prediction_core[-which(prediction_core$X.OTU_IDs %in% redundant_ko), ]
prediction_core3 <- prediction_core2[-which(prediction_core2$kegg == "None"),]
prediction_core4 <- prediction_core3[-which(prediction_core3$kegg == "hypothetical protein"),]

db_core <- as.data.frame(sample_data(sey_gut_core))
rownames(db_core) <- str_replace(rownames(db_core), "-", ".")

######## Which are the funtions lost or gained with the macroalgal transition system? ---- ##########
# ____ LEFSE on the KO of all samples with CCR vs MSR----
lefse <- data.frame(db_core[,19], stringsAsFactors = FALSE)
levels(lefse$geomorpho) <- c("CCR","MSR")
ko_core = prediction_core4
rownames(ko_core) <- ko_core$X.OTU_IDs
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core[1:99,])
taxo_table_core <- prediction_core4$X.OTU_IDs
tax_name_core <- taxo_table_core
tax_name_core <- c("Reef",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_all_reef.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size_reef <- read.delim("LDA_Effect_Size_reef", header=FALSE)
LDA_MSR<- subset(LDA_Effect_Size_reef , LDA_Effect_Size_reef$V3 == "MSR")
LDA_CCR<- subset(LDA_Effect_Size_reef , LDA_Effect_Size_reef$V3 == "CCR")
LDA_all_reef <- rbind(LDA_MSR,LDA_CCR)
colnames(LDA_all_reef) <- c("X.OTU_IDs" , "LDA_res" , "Reef" , "LDA_res2" , "sig")
write.table(LDA_all_reef, file = "LDA_all_reef.txt", sep = "\t", quote = T)

prediction_lda_reef <- left_join(LDA_all_reef[,c(1,3)], prediction_core4, by = "X.OTU_IDs")

CCR_others_id <- rownames(db[db$geomorpho == "coral" & !db$family %in% c("Scaridae", "Siganidae")])
Other_C <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% CCR_others_id])
MSR_others_id <- rownames(db[db$geomorpho == "macroalgal" & !db$family %in% c("Scaridae", "Siganidae")])
Other_M <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% MSR_others_id])

CCR_sc_id <- rownames(db[db$geomorpho == "coral" & db$family == "Scaridae"])
sc_C <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% CCR_sc_id])
MSR_sc_id <- rownames(db[db$geomorpho == "macroalgal" & db$family == "Scaridae"])
sc_M <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% MSR_sc_id])

CCR_sig_id <- rownames(db[db$geomorpho == "coral" & db$family == "Siganidae"])
sig_C <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% CCR_sig_id])
MSR_sig_id <- rownames(db[db$geomorpho == "macroalgal" & db$family == "Siganidae"])
sig_M <- rowSums(prediction_lda_reef[colnames(prediction_lda_reef) %in% MSR_sig_id])

lda_ab_reef <- as.data.frame(cbind(prediction_lda_reef$X.OTU_IDs, Other_C, Other_M, sc_C, sc_M, sig_C, sig_M))
colnames(lda_ab_reef)[1] <- "X.OTU_IDs"
functions <-prediction_core[,c(1,101,102,103, 104, 105)]
lda_ab_fun_reef <- left_join(as.data.frame(lda_ab_reef), functions, by = "X.OTU_IDs")

write.table(lda_ab_fun_reef, file = "lda_ab_fun_reef.txt" , sep ="\t", quote = T)

# Change the table in order to produce a barplot
barplot_lda_reef <- read_delim("barplot_reef.txt", "\t", escape_double = FALSE, col_types = cols(Abundance = col_number()), trim_ws = TRUE)

barplot <- ggplot(barplot_lda_reef, aes(x = Reef, y = Abundance, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("KO Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  xlab("Abundance")+
  scale_fill_manual(values = c("black", "darkred", "darkgreen")) +
  theme_bw() +
  theme(axis.line = element_line(size = 0),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, family = "serif"),
        axis.text.y = element_text(family = "serif", size = 16),
        plot.title = element_text(family = "serif", size = 15),
        legend.text = element_text(size = 0,  family = "serif"),
        legend.title = element_text(size = 0, family = "serif"),
        strip.text.x = element_text(size=8, family = "serif")) +
  facet_wrap(~ C, scale ="free")
barplot

pdf(file ="barplot_C.pdf", he = 10 , wi = 15)
barplot
dev.off()


######## Wich are those in the Scaridae and Siganidae?---- ##########
# ____ LEFSE on the KO of all samples with Scaridae, Siganidae and the others.----
lefse <- data.frame(db_core[,7], stringsAsFactors = FALSE)
levels(lefse$family) <- c(rep("Other",6), "Scaridae", "Other", "Siganidae")
ko_core = prediction_core4
rownames(ko_core) <- ko_core$X.OTU_IDs
ko_core <- ko_core[,-1]
ko_core <- t(ko_core)

lefse <- cbind(lefse, ko_core[1:99,])
taxo_table_core <- prediction_core4$X.OTU_IDs
tax_name_core <- taxo_table_core
tax_name_core <- c("fam",tax_name_core)
lefse2 <- cbind(tax_name_core,t(lefse))
write.table(lefse2, file="lefse_all2.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
LDA_Effect_Size <- read.delim("LDA_Effect_Size", header=FALSE)
LDA_sig_all<- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Siganidae")
LDA_sc_all<- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Scaridae")
LDA_all <- rbind(LDA_sig_all,LDA_sc_all)
colnames(LDA_all) <- c("X.OTU_IDs" , "LDA_res" , "Fam" , "LDA_res2" , "sig")
write.table(LDA_all, file = "LDA_all.txt", sep = "\t", quote = T)

prediction_lda <- left_join(LDA_all[,c(1,3)], prediction_core4, by = "X.OTU_IDs")

CCR_others_id <- rownames(db[db$geomorpho == "coral" & !db$family %in% c("Scaridae", "Siganidae")])
Other_C <- rowSums(prediction_lda[colnames(prediction_lda) %in% CCR_others_id])
MSR_others_id <- rownames(db[db$geomorpho == "macroalgal" & !db$family %in% c("Scaridae", "Siganidae")])
Other_M <- rowSums(prediction_lda[colnames(prediction_lda) %in% MSR_others_id])

CCR_sc_id <- rownames(db[db$geomorpho == "coral" & db$family == "Scaridae"])
sc_C <- rowSums(prediction_lda[colnames(prediction_lda) %in% CCR_sc_id])
MSR_sc_id <- rownames(db[db$geomorpho == "macroalgal" & db$family == "Scaridae"])
sc_M <- rowSums(prediction_lda[colnames(prediction_lda) %in% MSR_sc_id])

CCR_sig_id <- rownames(db[db$geomorpho == "coral" & db$family == "Siganidae"])
sig_C <- rowSums(prediction_lda[colnames(prediction_lda) %in% CCR_sig_id])
MSR_sig_id <- rownames(db[db$geomorpho == "macroalgal" & db$family == "Siganidae"])
sig_M <- rowSums(prediction_lda[colnames(prediction_lda) %in% MSR_sig_id])

lda_ab <- as.data.frame(cbind(prediction_lda$X.OTU_IDs, Other_C, Other_M, sc_C, sc_M, sig_C, sig_M))
colnames(lda_ab)[1] <- "X.OTU_IDs"
functions <-prediction_core[,c(1,101,102,103, 104, 105)]
lda_ab_fun <- left_join(as.data.frame(lda_ab), functions, by = "X.OTU_IDs")

write.table(lda_ab_fun, file = "lda_ab_fun3.txt" , sep ="\t", quote = T)

pathway_core <- as.data.frame(read_delim("core_Ref99NR/pathway_prediction.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
colnames(pathway_core)[1] <- "X.OTU_IDs"
pathway_core2 <- left_join(pathway_core, ko_pathway_sort, by = "X.OTU_IDs")
write.table(pathway_core2, file = "pathway_prediction.txt" , sep ="\t", quote = T)

# Change the table in order to produce a barplot
barplot_lda_kegg <- read_delim("barplot.txt", "\t", escape_double = FALSE, col_types = cols(Abundance = col_number()), trim_ws = TRUE)

# Which are those in common which could be those explained by the loss of Scaridae and apparition of Siganidae? ----
reef_discriminant <- levels(factor(barplot_lda_reef$X.OTU_IDs))
fish_discriminant <- levels(factor(barplot_lda_kegg$X.OTU_IDs))

barplot_ko <- barplot_lda_reef[barplot_lda_reef$X.OTU_IDs %in% fish_discriminant,]
discriminant <- levels(factor(barplot_ko$X.OTU_IDs))
final_prediction <- prediction_lda[prediction_lda$X.OTU_IDs %in% discriminant,]
# How many  KOs discriminated by the shift are mainly cared by Scaridae's bacteriome or Siganidae's bacteriome?
length(levels(factor(final_prediction[final_prediction$Fam == "Siganidae", ]$X.OTU_IDs)))


## Print the Barplots -----
C_level_selected <- read_delim("barplot_families.txt", "\t", escape_double = FALSE, col_types = cols(Abundance = col_number()), trim_ws = TRUE)

barplot <- ggplot(C_level_selected, aes(x = `Reef condition`, y = Abundance, fill = Family)) +
  geom_bar(stat="identity") +
  ylab("KO Abundance by metabolic cathegory") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  scale_fill_manual(values = c("black", "darkred", "darkgreen")) +
  theme_bw() +
  theme(axis.line = element_line(size = 0),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(size = 8, family = "serif"),
        axis.text.y = element_text(family = "serif", size = 16),
        plot.title = element_text(family = "serif", size = 15),
        legend.text = element_text(size = 12,  family = "serif"),
        legend.title = element_text(size = 0, family = "serif"),
        strip.text.x = element_text(size=8, family = "serif")) +
  facet_wrap(~ A, scale ="free")
barplot

pdf(file ="barplot_C.pdf", he = 10 , wi = 10)
barplot
dev.off()




