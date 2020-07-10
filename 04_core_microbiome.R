## Part IV: Determination and description of the enteric and algal core microbiome-----
### Colors-----
#Microbial community of fish is very high and it is very difficult to discern less than 10 taxa. That's why, we kept the same color for the same taxa, in order to be able to compare the graphs between compartemtn, diet or fish families. Of course, if your are only interested in some taxa, better choose few colors. For the bigger table (*i.e* tax table at Order level), **distinctColorPalette()** could be very useful but be awared that colors are not repeated.  Here is the code : 
phylum_colors <- c('Acidobacteria'='lavenderblush4',
                   'Actinobacteria'='darkblue',
                   'Armatimonadetes'='cadetblue3',
                   'Bacteroidetes'='cornflowerblue',
                   'Calditrichaeota'='azure3',
                   'Chloroflexi'='#DCE1D2',
                   'Cyanobacteria'='#DE6554',
                   'Dadabacteria'='brown2',
                   'Deferribacteres'='darkslategray1',
                   'Deinococcus-Thermus'='salmon4',
                   'Dependentiae'='sandybrown',
                   'Epsilonbacteraeota'='darkslateblue',
                   'Elusimicrobia'='plum1',
                   'Euryarchaeota'='hotpink4',
                   'Firmicutes'='brown4',
                   'Fusobacteria'='orange',
                   'Gemmatimonadetes'='darkolivegreen',
                   'Kiritimatiellaeota'='darkkhaki',
                   'Marinimicrobia_(SAR406_clade)'='darkgoldenrod4',
                   'Latescibacteria'='darkseagreen2',
                   'Lentisphaerae'='darkseagreen4',
                   'Patescibacteria'='darkturquoise',
                   'Planctomycetes'='darkslategray',
                   'Proteobacteria'='aquamarine4',
                   'Spirochaetes'='darkolivegreen3',
                   'Tenericutes'='#CA8EA7',
                   'Thaumarchaeota'='gold3',
                   'Verrucomicrobia'='darkgreen',
                   'WPS-2'='thistle2',
                   'Other'= 'black',
                   'Z-Other' = 'black')

class_colors <- c("Acidimicrobiia"="darksalmon",
                  "Acidobacteriia"="lavenderblush4",
                  "Actinobacteria"="darkblue",
                  "Alphaproteobacteria"="lightseagreen",
                  "Babeliae"="peachpuff",
                  "Anaerolineae"="tomato2",
                  "Bacilli"="brown4",
                  "Bacteroidia"="cornflowerblue",
                  "Brachyspirae"="darkolivegreen2",
                  "Campylobacteria"="royalblue",
                  "Clostridia"="orange3",
                  "Coriobacteriia"="deepskyblue4",
                  "Deferribacteres"="darkslategray1",
                  "Deinococci"="skyblue3",
                  "Deltaproteobacteria"="skyblue4",
                  "Erysipelotrichia"="yellow",
                  "Fusobacteriia"="orange",
                  "Fimbriimonadia" = "darkseagreen",
                  "Gammaproteobacteria"="aquamarine4",
                  "Kiritimatiellae"="darkgray",
                  "Lentisphaeria"="darkseagreen4",
                  "Microgenomatia"="seashell3",
                  "Mollicutes"="#CA8EA7",
                  "Negativicutes"="palevioletred4",
                  "Nitrososphaeria"="sandybrown",
                  "Oxyphotobacteria"="#DE6554",
                  "Phycisphaerae"="rosybrown4",
                  "Planctomycetacia"="darkslategray",
                  "Rhodothermia"="cornsilk3",
                  "Spirochaetia"="darkolivegreen3",
                  "Thermoanaerobaculia"="purple",
                  "Thermoleophilia"="deeppink4",
                  "Verrucomicrobiae"="darkgreen")

#library(randomcoloR)
#n <- length(levels(core_physeq_order$Order))
#order_colors <- distinctColorPalette(n)

####--------------------------------------------------------------------------------------
#Let's analyze our data previously proceeded. In order to descrimine which ASVs are rare (transients) and which would be considered as permanent (core), we use a dispersion index for each ASV and compare it at the Poisson distribution to determine which ASVs are signifcantly randomly distributed (transiant) and which are not (core) following the method of [Fillol et al. (2016)](https://www.nature.com/articles/ismej2015143).
### Part IVa: Enteric core microbiome------

#We create a new folder for the enteric microbiome of reef fish and determine the core.
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/Fish/Core/")
dir.create(dir_data_cleaning, recursive = T)
load(paste0(dir_taxa_assign, "sey_gut.RData"))
load(paste0(dir_taxa_assign, "seyrff_gut.RData"))

#### Core determination------
library(labdsv)
otu_table <- seyrff_gut@otu_table@.Data
abuoccplot_otu <- abuocc(otu_table)
#sub_objects of abuocc objects
str(abuoccplot_otu)
# transform spc.plt vector into table in order to calculate specific richness
richness_otu <- data.frame(abuoccplot_otu$spc.plt)
# occurence of each OTU
otu_occurence <- data.frame(abuoccplot_otu$plt.spc)
mean.abun_otu <- colSums(otu_table)/otu_occurence
square_otu <- otu_table^2
ss_otu <- data.frame(colSums(square_otu))
# Variance calculation
variance_otu=ss_otu/otu_occurence-mean.abun_otu^2
disp_otu <- (variance_otu/mean.abun_otu)*otu_occurence
# IC calculation for Poisson distribution using Chi square distribution (value and formula within Zar p574)
library(epitools)
poisic_otu = pois.exact(otu_occurence, conf.level = 0.95)
gut_dstat_otu <- cbind(mean.abun_otu, disp_otu, otu_occurence, poisic_otu)
names(gut_dstat_otu) <- c("average","disp", "occurence", "x", "pt", "rate", "lower", "upper", "prob")
save(gut_dstat_otu, file=paste0(dir_data_cleaning, "gut_sey_dstat_asv.RData"))
# Selection of core ASVs
gut_sey_core_otu <- gut_dstat_otu[gut_dstat_otu$disp > gut_dstat_otu$upper,]
gut_sey_core_otu <- na.exclude(gut_sey_core_otu)
gut_tax <- data.frame(seyrff_gut@tax_table@.Data)
gut_sey_core_otu$tax <- gut_tax[rownames(gut_tax) %in% row.names(gut_sey_core_otu),6]
gut_sey_core_otu$phylum <- gut_tax[rownames(gut_tax) %in% row.names(gut_sey_core_otu),2]
save(gut_sey_core_otu, file = paste0(dir_data_cleaning,"gut_sey_core_otu.Rdata"))

sey_gut_core <-prune_taxa(rownames(gut_sey_core_otu), seyrff_gut)
save(sey_gut_core, file = paste0(dir_data_cleaning, "sey_gut_core.RData"))


#To blast the reads of the core, we need to generate a fasta file with the function `writeXStringSet()`.
gut_names <- colnames(sey_gut_core@otu_table)
gut_tree <- subset_taxa(sey_gut_core, rownames(sey_gut_core@tax_table) %in% gut_names)
gut_tree
Biostrings::writeXStringSet(sey_gut_core@refseq, file = paste0(dir_data_cleaning,"sey_gut_core.fasta"))


#### Composition of the core-------
library(dplyr)
gut_core_order <- sey_gut_core %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Order)      

save(gut_core_order , file = paste0(dir_data_cleaning, "gut_core_order.RData"))

#Plot a treemap of the core thanks to the `treemap()` function for the main contributor of the core.
gut_core_order$Phylum = as.character(gut_core_order$Phylum) # Avoid error message with factor for next step
sort(table(factor(gut_core_order$Phylum)), T)
gut_core_order[!gut_core_order$Phylum %in% c("Proteobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Planctomycetes", "Spirochaetes", "Verrucomicrobia","Fusobacteria","Tenericutes"),which(names(gut_core_order) == "Phylum", T)] <- "Other" #Change phylum to other for those not included in the list

group <-  gut_core_order$Phylum
subgroup <- gut_core_order$Order
value <- gut_core_order$Abundance

gut_core_treemap_data=data.frame(group,subgroup,value)

library(treemap)
gut_core_treemap <- treemap(gut_core_treemap_data,
                            index=c("group","subgroup"), vSize = "value", type = "index",
                            fontcolor.labels=c("white","black"),
                            fontsize.labels=c(12),bg.labels=c("transparent"),
                            fontface.labels=c(2,3),
                            border.col=c("black","white"), border.lwds=c(4,2), 
                            align.labels=list(c("center", "center"),c("left", "bottom")),
                            title="Seychelles Enteric Core Treemap",fontsize.title=12,
                            fontfamily.title ="serif")


#and visualize the relative contribution of the different taxa (with the `fill=` argument for the level of the taxonomy you want represent) for the different reef condition or site with the argument `x=`.

#*Ex* : The composition of the enteric core microbiome between reef fish families at phylum level. 
gut_core_order$Phylum <- str_replace_all(gut_core_order$Phylum,c("Other" = "Z-Other"))

ggplot(gut_core_order, aes(x = family , y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Relative Abundance (Order > 2%)") +
  xlab("Family")+
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum Composition of Core Enteric Microbiome") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 12, family = "serif", angle = 45,vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 10,  family = "serif"), 
        legend.title = element_text(size = 12, family = "serif"))
#### Core ratio-------
#To have a representation of the size of the core in the entire community, we can calculate the ratio of the sequences. This can give us informations on the inter individual dispersion : the higher the intra-species dispersion, the lower the size of the core and the larger the variable bacterial community of the compartment or species. 
ratio_gut <- sum(sample_sums(sey_gut_core)) / sum(sample_sums(seyrff_gut))
percent(ratio_gut)

#### Itol----------
#In order to represent the taxa in the whole microbial description via Blast, we can then represent the source of our sequence (where they have already been described) in an interactive tree in [Itol](https://itol.embl.de/) (Interactive Tree of Life). To do that, we need to transform our phylogenetic tree from the phyloseq object (obtained after an alignement in ARB here) in an xml output thanks to the `write_xml()` function. 
devtools::install_github("USCBiostats/rphyloxml")
library(ape)
library(rphyloxml)
set.seed(12)
gut_tree <- sey_gut_core@phy_tree
gut_phyloxml <- write_phyloxml(gut_tree)
cat(as.character(gut_phyloxml))
xml2::write_xml(gut_phyloxml, "gut_phyloxml.xml")

#Then, in order to add supplementary envdata (other than Blast, as a gradient, host species, geography etc....) you have to use a special template .txt for Itol. The function `create_itol_files(env.xlsx)` from the script [table2itol.R](https://github.com/mgoeker/table2itol/blob/master/table2itol.R) is very useful to made automatically the templates from your env file. You are free to change the colors and shape as you want next. It faster than made a manually entry in the Itol website. 

###--------------------------------------------------------------------------------------------------------------------------------------
### Part IVb: Algal core microbiome---------
#We create a new folder for the algal microbiome and determine the core.
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/Algae/Core/")
dir.create(dir_data_cleaning, recursive = T)
load(paste0(dir_taxa_assign, "sey_algae.RData"))
load(paste0(dir_taxa_assign, "seyrff_algae.RData"))

#### Core determination-----
library(labdsv)
otu_table <- seyrff_algae@otu_table@.Data
abuoccplot_otu <- abuocc(otu_table)
#sub_objects of abuocc objects
str(abuoccplot_otu)
# transform spc.plt vector into table in order to calculate specific richness
richness_otu <- data.frame(abuoccplot_otu$spc.plt)
# occurence of each OTU
otu_occurence <- data.frame(abuoccplot_otu$plt.spc)
mean.abun_otu <- colSums(otu_table)/otu_occurence
square_otu <- otu_table^2
ss_otu <- data.frame(colSums(square_otu))
# Variance calculation
variance_otu=ss_otu/otu_occurence-mean.abun_otu^2
disp_otu <- (variance_otu/mean.abun_otu)*otu_occurence
# IC calculation for Poisson distribution using Chi square distribution (value and formula within Zar p574)
library(epitools)
poisic_otu = pois.exact(otu_occurence, conf.level = 0.95)
algae_dstat_otu <- cbind(mean.abun_otu, disp_otu, otu_occurence, poisic_otu)
names(algae_dstat_otu) <- c("average","disp", "occurence", "x", "pt", "rate", "lower", "upper", "prob")
save(algae_dstat_otu, file=paste0(dir_data_cleaning, "algae_sey_dstat_asv.RData"))
# Selection of core ASVs
algae_sey_core_otu <- algae_dstat_otu[algae_dstat_otu$disp > algae_dstat_otu$upper,]
algae_sey_core_otu <- na.exclude(algae_sey_core_otu)
algae_tax <- data.frame(seyrff_algae@tax_table@.Data)
algae_sey_core_otu$tax <- algae_tax[rownames(algae_tax) %in% row.names(algae_sey_core_otu),6]
algae_sey_core_otu$phylum <- algae_tax[rownames(algae_tax) %in% row.names(algae_sey_core_otu),2]
save(algae_sey_core_otu, file = paste0(dir_data_cleaning,"algae_sey_core_otu.Rdata"))

sey_algae_core <-prune_taxa(rownames(algae_sey_core_otu), seyrff_algae)
save(sey_algae_core, file = paste0(dir_data_cleaning, "sey_algae_core.RData"))

#To blast the reads of the core, we need to generate a fasta file with the function `writeXStringSet()`.
algae_names <- colnames(sey_algae_core@otu_table)
algae_tree <- subset_taxa(sey_algae_core, rownames(sey_algae_core@tax_table) %in% algae_names)
algae_tree
Biostrings::writeXStringSet(sey_algae_core@refseq, file = paste0(dir_data_cleaning,"sey_algae_core.fasta"))

#### Composition of the core------
library(dplyr)
algae_core_order <- sey_algae_core %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Order)      

save(algae_core_order , file = paste0(dir_data_cleaning, "algae_core_order.RData"))

#Plot a treemap of the core thanks to the `treemap()` function for the main contributor of the core. 
algae_core_order$Phylum = as.character(algae_core_order$Phylum) # Avoid error message with factor for next step
sort(table(factor(algae_core_order$Phylum)), T)
algae_core_order[!algae_core_order$Phylum %in% 
                   c("Proteobacteria", "Bacteroidetes","Cyanobacteria","Verrucomicrobia","Planctomycetes"),which(names(algae_core_order) == "Phylum", T)] <- "Other" #Change phylum to other for those not included in the list

group <-  algae_core_order$Phylum
subgroup <- algae_core_order$Order
value <- algae_core_order$Abundance

algae_core_treemap_data=data.frame(group,subgroup,value)

library(treemap)
algae_core_treemap <- treemap(algae_core_treemap_data,
                            index=c("group","subgroup"), vSize = "value", type = "index",
                            fontcolor.labels=c("white","black"),
                            fontsize.labels=c(12),bg.labels=c("transparent"),
                            fontface.labels=c(2,3),
                            border.col=c("black","white"), border.lwds=c(4,2), 
                            align.labels=list(c("center", "center"),c("left", "bottom")),
                            title="Seychelles Algae Core Treemap",fontsize.title=12,
                            fontfamily.title ="serif")

#and visualize the relative contribution of the different taxa (with the `fill=` argument for the level of the taxonomy you want represent) for the different reef condition or site with the argument `x=`.
#*Ex* : The composition of the core algal microbiome between Mahe and Praslin at phylum level. 
algae_core_order$Phylum <- str_replace_all(algae_core_order$Phylum,c("Other" = "Z-Other"))
ggplot(algae_core_order, aes(x = site, y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Relative Abundance (Order > 2%)") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum Composition of Algae between sites") +
  xlab("Site")+
  scale_fill_manual(values = phylum_colors) +
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
        legend.text = element_text(size = 10,  family = "serif"), 
        legend.title = element_text(size = 12, family = "serif"))

#### Core ratio------
#To have a representation of the size of the core in the entire community, we can calculate the ratio of the sequences. This can give us informations on the inter individual dispersion : the higher the intra-species dispersion, the lower the size of the core and the larger the variable bacterial community of the compartment or species. 
ratio_algae <- sum(sample_sums(sey_algae_core)) / sum(sample_sums(seyrff_algae))
percent(ratio_algae)

####-----------------------------------------------------------------------------------------------------------------
### Part IVc: Global core and comparison------------
#We create a new folder for the total core microbiome and compare both compatment in composition (LefSe) and in diversity (alpha and beta).
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/Total Core/")
dir.create(dir_data_cleaning, recursive = T)
load(paste0(dir_taxa_assign, "sey_final.RData"))
load(paste0(dir_taxa_assign, "seyrff_final.RData"))

#We use the otu tables of the enteric core and algal core to keep the ASVs names and subset them in the initial phyloseq object `sey_final`.
core_ASV <- unique(c(row.names(gut_sey_core_otu), row.names(algae_sey_core_otu)))
global_core <- subset_taxa(sey_final, taxa_names(sey_final) %in% core_ASV)
save(global_core, file = paste0(dir_data_cleaning, "global_core.RData"))

core_physeq_order <- global_core %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)
save(core_physeq_order , file = paste0(dir_data_cleaning, "core_physeq_order.RData"))

#Now we can visualize easily by compartments the bacterial composition with the chosen taxa levels (e.g. Phylum here). 
core_physeq_order$lineage <- str_replace_all(core_physeq_order$lineage, c("plantae" = "Algae", "vertebrate"="Gut"))

ggplot(core_physeq_order, aes(x = lineage , y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  scale_fill_manual(values = phylum_colors) +
  ylab("Relative Abundance (Phylum > 2%)") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Core Phylum Composition between Gut and Algae") +
  theme_bw() +
  xlab("Family")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank())  + 
  theme(axis.line = element_line(size = 0), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 12, family = "serif", angle = 0), 
        axis.text.y = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 10,  family = "serif"), 


#### Alpha diversity comparisons---------------------------

#We compare the taxonomic diversity (observed richness and exp(shannon index) for the effective number of species (ENS)) between the normalized tables (478 ASVs normalized) of the enteric microbiome and epiphytes bacteria. Note that 478 was chosen because Lethrinus samples were very poor but essential for the comparisons between health status of the reef. Following the richness curves made higher during the [Phyloseq process](#Part II: Phyloseq process)
box1 = plot_richness(seyrff_final , measures = c("Observed","Shannon") , color = "lineage")
box1$data[box1$data$variable == "Shannon",]$value = exp(box1$data[box1$data$variable == "Shannon",]$value)
levels(box1$data$variable)<- c("Richness observed", "Shannon (ENS)")
box1$data$lineage <- str_replace_all(box1$data$lineage, c("plantae" = "Algae", "vertebrate"="Gut"))

library(ggsignif)

palette = c("darkgreen" , "darkorange")

sey_final_alpha_div =  ggplot(box1$data, aes(x = lineage , y = value, color = lineage)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_colour_manual(values= palette) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=18, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 18, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=12, family = "serif"),
        legend.position="none") + 
  facet_wrap( ~ variable, nrow=1, ncol=5, scales = "free") +
  geom_signif(comparisons = list(c("Algae", "Gut")), 
              map_signif_level = TRUE, textsize=6, color="black", family = "serif", vjust = 0.5)
sey_final_alpha_div

library(pgirmess)
MW.richness = wilcox.test(box1$data$value[box1$data$variable =="Richness observed"] ~ box1$data$lineage[box1$data$variable=="Richness observed"], method = "bonferroni")
MW.richness = cbind(MW.richness$statistic , MW.richness$p.value)
MW.shannon = wilcox.test(box1$data$value[box1$data$variable=="Shannon (ENS)"] ~ box1$data$lineage[box1$data$variable=="Shannon (ENS)"], method = "bonferroni")
MW.shannon = cbind(MW.shannon$statistic , MW.shannon$p.value)
# Plot the results  
MW.test <- rbind(MW.richness, MW.shannon)
rownames(MW.test) = c("Observed richness","Shannon (ENS)")
colnames(MW.test) = c('W stat' , 'p-value')
MW.test

#### Beta diversity comparisons-----------
#Even we work on the core of normalized tables, we calculate the beta diversity ont the relative ASVs counts with the Bray-curtis distances with `vegdist` function from [vegan package (Oksanen et al., 2019)](https://CRAN.R-project.org/package=vegan).
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/Total Core/")
load(paste0(dir_data_cleaning, "global_core.RData"))
global_core_rel <- transform_sample_counts(global_core, function(x) x / sum(x) )
save(global_core_rel, file = paste0(dir_data_cleaning, "global_core_rel.RData"))
# PCOA coda & Permanova
library(vegan)
otu <- vegdist(global_core_rel@otu_table, method = "bray")
save(otu, file = paste0(dir_data_cleaning, "beta_matrices.RData"))
pcoa.sub <- pcoa(otu)
pcoa_coord <- pcoa.sub$vectors[,1:3]
save(pcoa.sub,pcoa_coord ,file = paste0(dir_data_cleaning, "pcoa.values.RData"))

# Contruction of the table for graphic 
library(stringr)
samp_data <- data.frame(sample_data(global_core_rel))
names(samp_data)[6] = "Type"
hull <- cbind(pcoa_coord, samp_data)
hull$Type <- str_replace_all(hull$Type, c("plantae" = "Algae", "vertebrate"="Gut"))

# What is the percentage of the explicative variance? 
paste("Axis 1 :",percent(pcoa.sub$values$Relative_eig[1])) # 8.6 %
paste("Axis 2 :",percent(pcoa.sub$values$Relative_eig[2])) # 6.6 %
paste("Axis 3 :",percent(pcoa.sub$values$Relative_eig[3])) # 5.1 %
paste("Axis 4 :",percent(pcoa.sub$values$Relative_eig[4])) # 4.9 %
palette = c("darkgreen", "darkorange")

# Plot
pcoa.bray <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = palette) +
  geom_point(data = hull, aes(x=Axis.1, y=Axis.2, color = Type), alpha = 0.7, size = 5, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=20, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=20, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 20, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"))+
  labs(colour = "Type", fill = "Type")

pcoa.bray

#We then compare the beta diversity with the `adonis` function and test the dispersion between samples with `betadisper`. 
library(pairwiseAdonis)
adonis(otu ~ Type, data = samp_data) # between gut, macroalgae and turf
pairwise.adonis(otu , samp_data$Type)

beta <- betadisper(otu, samp_data$Type)
permutest(beta)


#### LEfSe on compartments-------
#In order to detect biomarkers proper to each compartment, we proceed to a [Linear discriminant analysis on effect size (LEfSe) available on Galaxy](https://galaxyproject.org/learn/visualization/custom/lefse/). Here is the code to prepare the data input to proceed in Galaxy hub. 
dir_data_cleaning <- paste0(path, "/analyses/04_data_cleaning/Total Core/LDA/")
dir.create(dir_data_cleaning, recursive = T)

sey_compartment.genus <- tax_glom(global_core, "Genus", NArm = TRUE) #merges species that have the same Genus + filtre Genres NA
save(sey_compartment.genus, file = paste0(dir_data_cleaning, "sey_compartment.genus.RData")) 
taxo = data.frame(tax_table(sey_compartment.genus))
tax_name <- paste(taxo$Domain, taxo$Phylum, taxo$Class, taxo$Order, taxo$Family, taxo$Genus,  sep="|")
tax_name2 <- c("Type",tax_name)
save(taxo, tax_name2, file=paste0(dir_data_cleaning , "sey_compartment.genus_taxo.RData"))
samp_data <- data.frame(sample_data(sey_compartment.genus))
levels(samp_data$lineage) <- c("Algae", "Gut")
type_lefse <- data.frame(samp_data[,6], stringsAsFactors = FALSE)
otu = data.frame(otu_table(sey_compartment.genus))

lefse <- cbind(type_lefse, otu)
lefse2 <- cbind(tax_name2,t(lefse))
write.table(lefse2, file=paste0(dir_data_cleaning, "lefse_type_sey.txt"), 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#Proceed the "lefse_type_sey.txt" in a LEfse Galaxy and retrieve data "LDA Effect Size" with the associated "Plot LEfSe Results" in the LDA folder.You have to rename the former file "LDA_results" and run the following code. 
library(plyr)
LDA_Effect_Size <- read.delim(paste0(dir_data_cleaning, "LDA_results"), header=FALSE)
LDA_gut <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Gut")
LDA_algae <- subset(LDA_Effect_Size , LDA_Effect_Size$V3 == "Algae")
LDA_type <- rbind(LDA_gut,LDA_algae)
LDA_phylum <- LDA_type$V1
names(LDA_phylum) <- ".Phylum.Class.Order.Family.Genus"
write.table(LDA_phylum, file=paste0(dir_data_cleaning, "LDA_phylum.txt"), sep="\t", row.names=F, col.names=F, quote=FALSE)
#paste the names(LDA_phylum) in excel then read it again.

library(readr)
LDA_phylum <- read_delim(paste0(dir_data_cleaning, "LDA_phylum.txt"),".", escape_double = FALSE, col_types = cols(X1 = col_skip()),trim_ws = TRUE)

#It will transform it in tabble and convert "." by separator 
LDA_phylum_type <- cbind(LDA_phylum , LDA_type[,2:5])
colnames(LDA_phylum_type) <- c(names(LDA_phylum)[1:5] , "LDA_res" , "Type" , "LDA_res2" , "sig")
write.table(LDA_phylum_type, file=paste0(dir_data_cleaning, "LDA_phylum_type.txt"), sep="\t", row.names=F, col.names=T, quote=FALSE)
LDA_phylum_type <- read_delim(paste0(dir_data_cleaning, "LDA_phylum_type.txt"),"\t", escape_double = FALSE, trim_ws = TRUE)
# Select your significativ value and taxa
LEFSE_sig <- subset(LDA_phylum_type, LDA_phylum_type$LDA_res2 >= 3)
LEFSE_sig_genus <- subset(LEFSE_sig, LEFSE_sig$Genus %in% as.character(na.exclude(LEFSE_sig$Genus)))
table(LEFSE_sig_genus$Type)
LEFSE_sig_genus_vec <- paste(LEFSE_sig_genus$Phylum, LEFSE_sig_genus$Class, LEFSE_sig_genus$Order, LEFSE_sig_genus$Family, LEFSE_sig_genus$Genus, sep = "|")

biomarkers <- paste(LEFSE_sig_genus$Order, LEFSE_sig_genus$Genus, sep = "|")
save(biomarkers, file = paste0(dir_data_cleaning, "biomarkers.RData"))

#sey_samp <- sey_compartment.genus@sam_data
#write.table(sey_samp, file=paste0(dir_data_cleaning, "sey_samp.txt"), sep="\t", row.names=F, col.names=T, quote=FALSE) # You need to numerate the tax3 name for the turf and sargassum samples and re-read it. 
sey_samp <- read_table2(paste0(dir_data_cleaning, "sey_samp.txt"))
sample_data(sey_compartment.genus)$tax3 <- sey_samp$tax3
levels(sample_data(sey_compartment.genus)$lineage) <- c("Algae", "Gut")

physeq_phylum <- sey_compartment.genus %>%    
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%  
  filter(Abundance > 0) %>% # Melt to long format
  arrange(Genus)

physeq_phylum$Abundance <- physeq_phylum$Abundance/nrow(sample_data(sey_compartment.genus))

physeq_phylum$tax <- paste(physeq_phylum$Phylum, physeq_phylum$Class, physeq_phylum$Order, physeq_phylum$Family, physeq_phylum$Genus, sep = "|")
save(physeq_phylum, file= paste0(dir_data_cleaning, "physeq_phylum.RData"))
write.table(physeq_phylum ,file = paste0(dir_data_cleaning, "physeq_phylum.txt"), sep="\t",quote = F)

Diet_Pol_sey <- physeq_phylum
Diet_Pol_sey$tax <-  paste(Diet_Pol_sey$Order, Diet_Pol_sey$Genus, sep = "|")
length(physeq_phylum$tax)

df_sey <- Diet_Pol_sey[,c("lineage","Sample","tax","Abundance")]
colnames(df_sey) <- c("family","item","score","value")
save(df_sey, file = paste0(dir_data_cleaning, "df_sey.RData"))


#Let's organize our data to plot corresponding colours depending on the  compartment. 
#Polarhisogram -------
load(paste0(dir_data_cleaning , "df_sey.RData"))
load(paste0(dir_data_cleaning , "biomarkers.RData"))

library(plyr)
print_biomarkers <- levels(factor(df_sey[df_sey$score %in% biomarkers,]$score)) # All biomarkers present in the sey.compartment_genus RData file
others <- levels(factor(df_sey[!df_sey$score %in% biomarkers,]$score)) # All taxa which are not biomarkers
df_sey[!df_sey$score %in% biomarkers,]$score <- "Z-Other" # Call them other (Z to figure at the end of the list)

# Some biomarkers are in common so Output the frequences.
freq_gut <- table(df_sey[df_sey$family == "Gut",]$score)/length(df_sey[df_sey$family == "Gut",]$score)*100 
freq_gut <- freq_gut[!names(freq_gut) == "Z-Other"]
freq_algae <- table(df_sey[df_sey$family == "Algae",]$score)/length(df_sey[df_sey$family == "Algae",]$score)*100 
freq_algae <- freq_algae[!names(freq_algae) == "Z-Other"]

# Keep the common biomarkers of the algae and gut and the respective frequencies in each compartment
biom_common <- names(freq_algae[which(names(freq_algae) %in% names(freq_gut))])
freq_gut_common <- as.data.frame(freq_gut[names(freq_gut) %in% biom_common])
freq_algae_common <- as.data.frame(freq_algae[names(freq_algae) %in% biom_common])
freq_common <- cbind(freq_gut_common, freq_algae_common$Freq)
colnames(freq_common) <- c("Biomarkers", "Freq_gut" ,"Freq_algae")
# Class the biomarkers in the compartment where the value is the higher
biom_gut_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_gut > freq_common$Freq_algae, T)]
biom_algae_to_keep <- freq_common$Biomarkers[which(freq_common$Freq_algae > freq_common$Freq_gut, T)]
# See which are the biomarkers for the gut and which are the biomarkers for the algae
biom_gut <- levels(factor(df_sey[df_sey$family == "Gut",]$score))
biom_gut <- biom_gut[!biom_gut == "Z-Other"]
biom_gut <- biom_gut[!biom_gut %in% biom_algae_to_keep]
biom_algae <- levels(factor(df_sey[df_sey$family == "Algae",]$score))
biom_algae <- biom_algae[!biom_algae == "Z-Other"]
biom_algae <- biom_algae[!biom_algae %in% biom_gut_to_keep]

# Color choice
algae_colors <- c("darkgreen" , "darkkhaki", "darkolivegreen",
                  "darkolivegreen2","forestgreen","chartreuse",
                  "aquamarine","aquamarine3", "darkcyan",
                  "darkseagreen", "yellowgreen", "darkslategray",
                  "gold3", "gold4", "goldenrod",
                  "gold","palegreen4","chartreuse3",
                  "palegreen3", "seagreen3","palegoldenrod","yellow")
pie(rep(1, 22), col= algae_colors)

gut_colors <- c("brown", "brown1", "burlywood4",
                "chocolate","chocolate4", "coral2",
                "darksalmon","darkred","darkorchid4",
                "darkmagenta", "blueviolet","darkblue",
                "blue3", "deeppink4","deeppink")
pie(rep(1, 15), col= gut_colors)
polar_col <- c(algae_colors , gut_colors , 'Z_Other'="black")
# Add a prefix G (for gut) or A (for algae) in their respectiv biomarkers. In this way, they will be listed following 
# their respectiv compartment and no their names (easier for the add of the color in the ggplot)
df_sey[df_sey$score %in% biom_gut,]$score <- paste("G",df_sey[df_sey$score %in% biom_gut,]$score, sep= "-")
df_sey[df_sey$score %in% biom_algae,]$score <- paste("A",df_sey[df_sey$score %in% biom_algae,]$score, sep= "-")


#We use the script of [Ladroue et al.,(2012)](http://dx.doi.org/10.1038/ng.1073.) to draw the histogram with the `polarHistogram` function available [here](https://github.com/burakaydin/materyaller/blob/master/mugla_ejeR/polarHistogram.R).

#Print the polarhistogramm
source(paste0(dir_refdb, "polarHistogram.R"))
p <- polarHistogram(df_sey , familyLabel=T, innerRadius = 0.2, spaceFamily =7, circleProportion = 0.85)
p + ggplot2::scale_fill_manual(values= polar_col) +
  theme(legend.text = element_text(family = "serif", size = 3)) +
  theme(legend.position= "none") # to hide the legend because of the size

