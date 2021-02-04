## ** LEFsE on scaridae and siganidae ** ------
# ____Bacteria LEfSE file creation ----
dir.create("./LEfSe")
setwd("./LEfSe")

db <- as.data.frame(sample_data(sey_gut_core))
write.table(db, file= "db.txt", sep ="\t", quote = F)
db <- read_delim("db.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
db <-as.data.frame(db)
rownames(db) <- db$X1
db <- db[,-1]

# Between Siganidae and other algal covered reefs inhabitants -----
dir.create("./Sig")
setwd("./Sig")

ps.M <- subset_samples(sey_gut_core, geomorpho == "macroalgal")
db_M <- as.data.frame(sample_data(ps.M))
write.table(db_M, file= "db_M.txt", sep ="\t", quote = F)
db_M <- read_delim("db_M.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
db_M <-as.data.frame(db_M)
rownames(db_M) <- db_M$X1
db_M <- db_M[,-1]

db.genus <- tax_glom(ps.M, "Genus", NArm=T)
taxo = data.frame(tax_table(db.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Sig",tax_name)
save(taxo, tax_name2, file="sig.genus_taxo.RData")

algae_lefse <- data.frame(db_M[,7], stringsAsFactors = FALSE)
otu = data.frame(otu_table(db.genus))

# On family
lefse <- cbind(algae_lefse, otu)
lefse2 <- cbind(tax_name2,t(lefse))
write.table(lefse2, file="lefse_algae_sey.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
#
library(plyr)
LDA_Effect_Size <- read.delim("LDA_results", header=FALSE)
LDA_sig <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Siganidae")
LDA_others <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Others")
LDA_M <-  rbind(LDA_sig,LDA_others)
LDA_phylum <- LDA_M$V1
names(LDA_phylum) <- as.factor(".Phylum.Class.Order.Family.Genus") #paste the names(LDA_phylum) in excel
write.table(LDA_phylum, file="LDA_phylum.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
library(readr)
LDA_phylum <- read_delim("LDA_phylum.txt",  ".", escape_double = FALSE, col_types = cols(X1 = col_skip()),  trim_ws = TRUE)

LDA_phylum_M <- cbind(LDA_phylum , LDA_M[,2:5])
colnames(LDA_phylum_M) <- c(names(LDA_phylum)[1:5] , "LDA_res" , "Family" , "LDA_res2" , "sig")
write.table(LDA_phylum_M, file= "LDA_phylum_M.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

# Select your significativ value and taxa
LEFSE_sig <- subset(LDA_phylum_M, LDA_phylum_M$LDA_res2 >= 3)
LEFSE_sig_genus <- subset(LEFSE_sig, LEFSE_sig$Genus %in% as.character(na.exclude(LEFSE_sig$Genus)))
LEFSE_sig_genus_vec <- paste(LEFSE_sig_genus$Phylum, LEFSE_sig_genus$Class, LEFSE_sig_genus$Order, LEFSE_sig_genus$Family, LEFSE_sig_genus$Genus, sep = "|")

biomarkers <- paste(LEFSE_sig_genus$Order, LEFSE_sig_genus$Genus, sep = "|")
save(biomarkers, file = "biomarkers.RData")
library(dplyr)

physeq_phylum <- db.genus %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%
  filter(Abundance > 0) %>% # Melt to long format
  arrange(Genus)

physeq_phylum$Abundance <- physeq_phylum$Abundance/nrow(sample_data(db.genus))

physeq_phylum$tax <- paste(physeq_phylum$Phylum, physeq_phylum$Class, physeq_phylum$Order, physeq_phylum$Family, physeq_phylum$Genus, sep = "|")
save(physeq_phylum, file= "physeq_phylum_M.RData")
#write.table(physeq_phylum ,file = paste0(dir_data_cleaning, "physeq_phylum.txt"), sep="\t",quote = F)

M_Pol_sey <- physeq_phylum
levels(M_Pol_sey$family) <- c("Others", "Others", "Others", "Others","Others","Others", "Others", "Siganidae")
M_Pol_sey$tax <-  paste(M_Pol_sey$Order, M_Pol_sey$Genus, sep = "|")
length(physeq_phylum$tax) #482

df_sey <- M_Pol_sey[,c("family","tax3","tax","Abundance")]
colnames(df_sey) <- c("family","item","score","value")
save(df_sey, file =  "df_sey.M.RData")
library(plyr)
print_biomarkers <- levels(factor(df_sey[df_sey$score %in% biomarkers,]$score)) # All biomarkers present in the sey.diet_genus RData file
others <- levels(factor(df_sey[!df_sey$score %in% biomarkers,]$score)) # All taxa which are not biomarkers
df_sey[!df_sey$score %in% biomarkers,]$score <- "Z-Other" # Call them other (Z to figure at the end of the list)

# Some biomarkers are in common so Output the frequences.
freq_sig <- table(df_sey[df_sey$family == "Siganidae",]$score)/length(df_sey[df_sey$family == "Siganidae",]$score)*100
freq_sig <- freq_sig[!names(freq_sig) == "Z-Other"]
freq_others <- table(df_sey[df_sey$family == "Others",]$score)/length(df_sey[df_sey$family == "Others",]$score)*100
freq_others <- freq_others[!names(freq_others) == "Z-Other"]
# Keep the common biomarkers of the algae and gut and the respective frequencies in each compartment
biom_common <- names(freq_sig[which(names(freq_sig) %in% names(freq_others))])
freq_others_common <- as.data.frame(freq_others[names(freq_others) %in% biom_common])
freq_sig_common <- as.data.frame(freq_sig[names(freq_sig) %in% biom_common])
freq_common <- cbind(freq_others_common, freq_sig_common$Freq)
colnames(freq_common) <- c("Biomarkers", "Freq_others" ,"Freq_sig")

# Class the biomarkers in the compartment where the value is the higher
biom_others_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_others > freq_common$Freq_sig, T)]
biom_sig_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_sig > freq_common$Freq_others, T)]

# See which are the biomarkers for the gut and which are the biomarkers for the algae
biom_others <- levels(factor(df_sey[df_sey$family == "Others",]$score))
biom_others <- biom_others[!biom_others == "Z-Other"]
biom_others <- biom_others[!biom_others %in% biom_sig_to_keep]
biom_sig <- levels(factor(df_sey[df_sey$family == "Siganidae",]$score))
biom_sig <- biom_sig[!biom_sig == "Z-Other"]
biom_sig <- biom_sig[!biom_sig %in% biom_others_to_keep]

# Color choice
sig_colors <- c("darkgreen" , "darkkhaki", "darkolivegreen",
                "darkolivegreen2","forestgreen","chartreuse",
                "aquamarine","aquamarine3", "darkcyan",
                "darkseagreen", "yellowgreen", "darkslategray",
                "gold3")

others_colors <- c("brown", "brown1", "burlywood4",
                   "chocolate","chocolate4", "coral2",
                   "darksalmon")
polar_col <- c(others_colors , sig_colors , 'Z_Other'="black")
# Add a prefix G (for gut) or A (for algae) in their respectiv biomarkers. In this way, they will be listed following
# their respectiv compartment and no their names (easier for the add of the color in the ggplot)
df_sey[df_sey$score %in% biom_others,]$score <- paste("O",df_sey[df_sey$score %in% biom_others,]$score, sep= "-")
df_sey[df_sey$score %in% biom_sig,]$score <- paste("Sig",df_sey[df_sey$score %in% biom_sig,]$score, sep= "-")

library(ggplot2)

pdf("polar_M.res3.pdf", he= 7, wi= 20)
p <- polarHistogram(df_sey , familyLabel=T, innerRadius = 0.2, spaceFamily =7, circleProportion = 0.85)
p + ggplot2::scale_fill_manual(values= polar_col)
dev.off()

# Between scaridae and other coral inhabitants -----
dir.create("./Sig")
setwd("./Sig")

ps.C <- subset_samples(sey_gut_core, geomorpho == "coral")
db_C <- as.data.frame(sample_data(ps.C))
write.table(db.C, file= "db.C.txt", sep ="\t", quote = F)
db_C <- read_delim("db.C.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
db_C <-as.data.frame(db_C)
rownames(db_C) <- db_C$X1
db_C <- db_C[,-1]

db.genus <- tax_glom(ps.C, "Genus", NArm=T)
taxo = data.frame(tax_table(db.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Scar",tax_name)
save(taxo, tax_name2, file="scar.genus_taxo.RData")

coral_lefse <- data.frame(db_C[,7], stringsAsFactors = FALSE)
otu = data.frame(otu_table(db.genus))

# On family
lefse <- cbind(coral_lefse, otu)
lefse2 <- cbind(tax_name2,t(lefse))
write.table(lefse2, file="lefse_coral_sey.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Run your file in LEfSE software on Galaxy, download the results, move it on your active folder "dir_data_cleaning" and rename the file "LDA_results".
#
library(plyr)
LDA_Effect_Size <- read.delim("LDA_results", header=FALSE)
LDA_sc <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Scaridae")
LDA_others <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Others")
LDA_C <-  rbind(LDA_sc,LDA_others)
LDA_phylum <- LDA_C$V1
names(LDA_phylum) <- as.factor(".Phylum.Class.Order.Family.Genus") #paste the names(LDA_phylum) in excel
write.table(LDA_phylum, file="LDA_phylum.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
library(readr)
LDA_phylum <- read_delim("LDA_phylum.txt",  ".", escape_double = FALSE, col_types = cols(X1 = col_skip()),  trim_ws = TRUE)

LDA_phylum_C <- cbind(LDA_phylum , LDA_C[,2:5])
colnames(LDA_phylum_C) <- c(names(LDA_phylum)[1:5] , "LDA_res" , "Family" , "LDA_res2" , "sig")
write.table(LDA_phylum_C, file= "LDA_phylum_C.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

# Select your significativ value and taxa
LEFSE_sig <- subset(LDA_phylum_C, LDA_phylum_C$LDA_res2 >= 3)
LEFSE_sig_genus <- subset(LEFSE_sig, LEFSE_sig$Genus %in% as.character(na.exclude(LEFSE_sig$Genus)))
LEFSE_sig_genus_vec <- paste(LEFSE_sig_genus$Phylum, LEFSE_sig_genus$Class, LEFSE_sig_genus$Order, LEFSE_sig_genus$Family, LEFSE_sig_genus$Genus, sep = "|")

biomarkers <- paste(LEFSE_sig_genus$Order, LEFSE_sig_genus$Genus, sep = "|")
save(biomarkers, file = paste0(dir_data_cleaning, "biomarkers.RData"))
library(dplyr)

physeq_phylum <- db.genus %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%
  filter(Abundance > 0) %>% # Melt to long format
  arrange(Genus)

physeq_phylum$Abundance <- physeq_phylum$Abundance/nrow(sample_data(db.genus))

physeq_phylum$tax <- paste(physeq_phylum$Phylum, physeq_phylum$Class, physeq_phylum$Order, physeq_phylum$Family, physeq_phylum$Genus, sep = "|")
save(physeq_phylum, file= "physeq_phylum_C.RData")
#write.table(physeq_phylum ,file = paste0(dir_data_cleaning, "physeq_phylum.txt"), sep="\t",quote = F)

C_Pol_sey <- physeq_phylum
levels(C_Pol_sey$family) <- c("Others", "Others", "Others", "Others","Scaridae","Others","Others")
C_Pol_sey$tax <-  paste(C_Pol_sey$Order, C_Pol_sey$Genus, sep = "|")
length(physeq_phylum$tax) #463


df_sey <- C_Pol_sey[,c("family","tax3","tax","Abundance")]
colnames(df_sey) <- c("family","item","score","value")
save(df_sey, file =  "df_sey.C.RData")
library(plyr)
print_biomarkers <- levels(factor(df_sey[df_sey$score %in% biomarkers,]$score)) # All biomarkers present in the sey.diet_genus RData file
others <- levels(factor(df_sey[!df_sey$score %in% biomarkers,]$score)) # All taxa which are not biomarkers
df_sey[!df_sey$score %in% biomarkers,]$score <- "Z-Other" # Call them other (Z to figure at the end of the list)

# Some biomarkers are in common so Output the frequences.
freq_sc <- table(df_sey[df_sey$family == "Scaridae",]$score)/length(df_sey[df_sey$family == "Scaridae",]$score)*100
freq_sc <- freq_sc[!names(freq_sc) == "Z-Other"]
freq_others <- table(df_sey[df_sey$family == "Others",]$score)/length(df_sey[df_sey$family == "Others",]$score)*100
freq_others <- freq_others[!names(freq_others) == "Z-Other"]
# Keep the common biomarkers of the algae and gut and the respective frequencies in each compartment
biom_common <- names(freq_sc[which(names(freq_sc) %in% names(freq_others))])
freq_others_common <- as.data.frame(freq_others[names(freq_others) %in% biom_common])
freq_sc_common <- as.data.frame(freq_sc[names(freq_sc) %in% biom_common])
freq_common <- cbind(freq_others_common, freq_sc_common$Freq)
colnames(freq_common) <- c("Biomarkers", "Freq_others" ,"Freq_sc")

# Class the biomarkers in the compartment where the value is the higher
biom_others_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_others > freq_common$Freq_sc, T)]
biom_sc_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_sc > freq_common$Freq_others, T)]

# See which are the biomarkers for the gut and which are the biomarkers for the algae
biom_others <- levels(factor(df_sey[df_sey$family == "Others",]$score))
biom_others <- biom_others[!biom_others == "Z-Other"]
biom_others <- biom_others[!biom_others %in% biom_sc_to_keep]
biom_sc <- levels(factor(df_sey[df_sey$family == "Scaridae",]$score))
biom_sc <- biom_sc[!biom_sc == "Z-Other"]
biom_sc <- biom_sc[!biom_sc %in% biom_others_to_keep]

# Color choice
sc_colors <- c("brown", "brown1", "burlywood4",
               "chocolate","chocolate4", "coral2")

others_colors <- c("darkgreen" , "darkkhaki")
polar_col <- c(others_colors , sc_colors , 'Z_Other'="black")
# Add a prefix G (for gut) or A (for algae) in their respectiv biomarkers. In this way, they will be listed following
# their respectiv compartment and no their names (easier for the add of the color in the ggplot)
df_sey[df_sey$score %in% biom_others,]$score <- paste("O",df_sey[df_sey$score %in% biom_others,]$score, sep= "-")
df_sey[df_sey$score %in% biom_sc,]$score <- paste("Sc",df_sey[df_sey$score %in% biom_sc,]$score, sep= "-")

library(ggplot2)

pdf("polar_C.res3.pdf", he= 7, wi= 20)
p <- polarHistogram(df_sey , familyLabel=T, innerRadius = 0.2, spaceFamily =7, circleProportion = 0.85)
p + ggplot2::scale_fill_manual(values= polar_col)
dev.off()

# Proportions biomarkers -----
# __ For Scaridae
setwd("../Sc")

gen_biom_sc <- c("Anaeroplasma", "Odoribacter", "Cellulosilyticum", "Lachnoclostridium", "Cetobacterium","Fusobacterium")
ps.biom.sc <- subset_taxa(sey_gut_core, Genus %in% gen_biom_sc)
biom_sc_ab <- ps.biom.sc %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0) %>%                         # Filter out low abundance taxa
  arrange(Genus)

save(ps.biom.sc, file= "./Sc/ps.biom.sc.RData")
save(biom_sc_ab, file= "./Sc/biom_sc_ab.RData")

levels(biom_sc_ab$family) <- c(rep("Others", 6), "Scaridae","Others","Siganidae")
levels(biom_sc_ab$geomorpho) <- c("C" , "M")
# __ For Siganidae
setwd("../Sig")

gen_biom_sig <- c("Labilibacter","Rikenella","Rikenellaceae_RC9_gut_group","Brachyspira", "Anaerofilum","Epulopiscium","Mucispirillum","Desulfovibrio",
                  "Erysipelatoclostridium","Fusobacterium","Ureaplasma","Treponema_2","Akkermansia")
ps.biom.sig <- subset_taxa(sey_gut_core, Genus %in% gen_biom_sig)
biom_sig_ab <- ps.biom.sig %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0) %>%                         # Filter out low abundance taxa
  arrange(Genus)

save(ps.biom.sig, file= "./Per conditions/Sig/ps.biom.sig.RData")
save(biom_sig_ab, file= "./Per conditions/Sig/biom_sig_ab.RData")

levels(biom_sig_ab$family) <- c(rep("Others", 6),"Scaridae","Others","Siganidae")
levels(biom_sig_ab$geomorpho) <- c("C" , "M")

# __ Compute all
setwd("../.")

all_ab <- rbind(biom_sig_ab, biom_sc_ab)
pdf(file = "abundance of biom.pdf", he=10 , wi = 15)
ggplot(all_ab, aes(x = geomorpho, y = Abundance, fill = family)) +
  geom_bar(stat="identity") +
  ylab("Abundance of each biomarkers") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Contribution of the biomarkers in the reefs") +
  scale_fill_manual(values = c("black","darkred", "darkgreen")) +
  theme_bw() +
  theme(axis.line = element_line(size = 0),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(family = "serif", size = 13),
        plot.title = element_text(family = "serif", size = 15),
        legend.text = element_text(size = 14,  family = "serif"),
        legend.title = element_text(size = 14, family = "serif"),
        strip.text.x = element_text(size=12, family = "serif")) +
  facet_wrap(~ Genus, scales = "free_y")
dev.off()
