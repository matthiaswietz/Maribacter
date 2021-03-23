##########################################################
 ### Wolter et al; Maribacter 62-1
 ### Front. Microbiol. | doi: 10.3389/fmicb.2021.628055 
##########################################################

# This Rscript evaluates ANI, CAZyme and Growth data
# See workflow at https://github.com/matthiaswietz/Maribacter 

library(reshape2)
library(dplyr)
library(stringr)
library(splitstackshape)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(PNWColors)


####################################################

## ANI MATRIX ##

# Import full matrix 
ANI <- read.table(
  "ANI_matrix.txt", sep="\t", row.names=1,
  header=T,  stringsAsFactors=F) 

# convert
ANI <- as.matrix(as.matrix(ANI)) 

##########################

# Get lower matrix triangle 
get_lower_tri<-function(ANI){
  ANI[upper.tri(ANI)] <- NA
  return(ANI)}

# Get upper matrix triangle 
get_upper_tri <- function(ANI){
  ANI[lower.tri(ANI)]<- NA
  return(ANI)}

##########################

# reorder to match tree

ANI <- ANI[, c(
  "Muricauda_alvinocaridis_SCR12",
  "Arenibacter_aquaticus_GUO",
  "Pricia_antarctica_DSM23421",
  "Pseudozobellia_thermophila_DSM19858",
  "Zobellia_amurskyensis_KMM3526",
  "Zobellia_uliginosa_MAR.2009.138",
  "Zobellia_uliginosa_DSM2061",
  "Zobellia_galactanivorans_DsiJ",
  "Maribacter_thermophilus_HT7.2",
  "Maribacter_MAR.2009.72",
  "Maribacter_sedimenticola_DSM19840",
  "Maribacter_ulvicola_DSM15366",
  "Maribacter_forsetii_DSM18668",
  "Maribacter_caenipelagi_CECT8455",
  "Maribacter_aquivivus_DSM16478",
  "Maribacter_stanieri_DSM19891",
  "Maribacter_Hel_I_7",
  "Maribacter_litoralis_SDRB.Phe2",
  "Maribacter_dokdonensis_MAR.2009.60",
  "Maribacter_dokdonensis_MAR.2009.71",
  "Maribacter_2014MBL_MicDiv",
  "Maribacter_621",
  "Maribacter_dokdonensis_DSW.8")]

ANI <- ANI[c(
  "Muricauda_alvinocaridis_SCR12",
  "Arenibacter_aquaticus_GUO",
  "Pricia_antarctica_DSM23421",
  "Pseudozobellia_thermophila_DSM19858",
  "Zobellia_amurskyensis_KMM3526",
  "Zobellia_uliginosa_MAR.2009.138",
  "Zobellia_uliginosa_DSM2061",
  "Zobellia_galactanivorans_DsiJ",
  "Maribacter_thermophilus_HT7.2",
  "Maribacter_MAR.2009.72",
  "Maribacter_sedimenticola_DSM19840",
  "Maribacter_ulvicola_DSM15366",
  "Maribacter_forsetii_DSM18668",
  "Maribacter_caenipelagi_CECT8455",
  "Maribacter_aquivivus_DSM16478",
  "Maribacter_stanieri_DSM19891",
  "Maribacter_Hel_I_7",
  "Maribacter_litoralis_SDRB.Phe2",
  "Maribacter_dokdonensis_MAR.2009.60",
  "Maribacter_dokdonensis_MAR.2009.71",
  "Maribacter_2014MBL_MicDiv",
  "Maribacter_621",
  "Maribacter_dokdonensis_DSW.8"),]

# Melt the matrix
ANI.matrix <- get_lower_tri(ANI)
ANI.matrix <- melt(ANI.matrix, na.rm=T)

##########################

ggplot(data = ANI.matrix, 
  aes(Var2, Var1, fill=value)) +
geom_tile(color = "white") +
scale_fill_gradientn(
  colors = pnw_palette("Sailboat", 3)) +
theme_bw() + 
theme(axis.text.x = element_text(
    angle=90, size=10, vjust=0.5, hjust=1))+
coord_fixed() 

# remove temp-data
rm(get_lower_tri, get_upper_tri)


####################################################

## CAZYmes -- dbCan ##

# faa file of each strain annotated using dbCAN2 
# hmm output saved as txtfile
# change Gene ID = locus_tag / E Value = evalue
# in Excel: eliminate duplicate assignments of locus-tags
# Combine all in one txtfile "dbcanOutput"

dbcan <- read.table(
  "dbcanOutput.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors=F)

# Filter & subset 
dbcan = filter(
  dbcan, Coverage >= 0.65 & evalue < 1e-15)
dbcan$HMMProfile = gsub(
  dbcan$HMMProfile, pattern="_.*|.hmm", replacement = "")

# Merge with full strain info
StrainInfo = read.csv(
  "AllStrains_anno.csv", header=T)

dbcan = merge(
  dbcan, StrainInfo, by=c("protein_id","Strain"), all.x=T)

# Export full, cleaned results 
write.table(
  dbcan, file = "../dbcan/dbcanFinal.txt", 
  sep ="\t", row.names=F)

############################

# Modify for count summary
dbcan$CAZyme = dbcan$HMMProfile
dbcan$CAZyme =  gsub('[[:digit:]]+', '', dbcan$CAZyme)

dbcanCount <- plyr::count(dbcan, c("Strain","CAZyme")) 

# Order according to genome phylogeny
dbcanCount$Strain <- factor(dbcanCount$Strain, levels = c(
  "Maribacter dokdonensis DSW-8",
  "Maribacter 62-1",
  "Maribacter_2014MBL_MicDiv",
  "Maribacter dokdonensis MAR_2009_71",
  "Maribacter dokdonensis MAR_2009_60",
  "Maribacter litoralis SDRB-Phe2",
  "Maribacter Hel_I_7",
  "Maribacter stanieri DSM19891",
  "Maribacter forsetii DSM18668",
  "Maribacter ulvicola DSM15366",
  "Maribacter aquivivus DSM16478",
  "Maribacter caenipelagi CECT8455",
  "Maribacter sedimenticola DSM19840",
  "Maribacter sp. MAR_2009_72",
  "Maribacter thermophilus HT7-2",
  "Zobellia galactanivorans DsiJ",
  "Zobellia uliginosa DSM2061",
  "Zobellia uliginosa MAR_2009_138",
  "Zobellia amurskyensis KMM3526",
  "Pseudozobellia thermophila DSM19858",
  "Pricia antarctica DSM23421",
  "Arenibacter aquaticus GUO",
  "Muricauda alvinocaridis SCR12"))

# Reshape 
dbcanCount <- reshape2::dcast(
  dbcanCount, Strain ~ CAZyme)

# Export
write.table(
  dbcanCount, file = "dbcanCount.txt", 
  sep ="\t", row.names=F)

############################

# Count for heatmap
dbcan = distinct(
  dbcan, HMMProfile, locus_tag, Strain) %>%
  group_by(HMMProfile, Strain) %>%
  filter(str_detect(HMMProfile, "PL|GH|CBM")) %>%
  tally()

# Order according to genome phylogeny
dbcan$Strain <- factor(dbcan$Strain, levels = c(
  "Maribacter dokdonensis DSW-8",
  "Maribacter 62-1",
  "Maribacter_2014MBL_MicDiv",
  "Maribacter dokdonensis MAR_2009_71",
  "Maribacter dokdonensis MAR_2009_60",
  "Maribacter litoralis SDRB-Phe2",
  "Maribacter Hel_I_7",
  "Maribacter stanieri DSM19891",
  "Maribacter forsetii DSM18668",
  "Maribacter ulvicola DSM15366",
  "Maribacter aquivivus DSM16478",
  "Maribacter caenipelagi CECT8455",
  "Maribacter sedimenticola DSM19840",
  "Maribacter sp. MAR_2009_72",
  "Maribacter thermophilus HT7-2",
  "Zobellia galactanivorans DsiJ",
  "Zobellia uliginosa DSM2061",
  "Zobellia uliginosa MAR_2009_138",
  "Zobellia amurskyensis KMM3526",
  "Pseudozobellia thermophila DSM19858",
  "Pricia antarctica DSM23421",
  "Arenibacter aquaticus GUO",
  "Muricauda alvinocaridis SCR12"))

# Subset to specific CAZymes; selected from test heatmap
dbcan = filter(dbcan, HMMProfile %in% c(
  "GH105","GH110","GH130","GH28",
  "PL1", "PL12","PL17","PL29","PL10",
  "PL33","PL38","PL40","PL6","PL7"))

# Reorder CAZymes
dbcan$HMMProfile <- factor(dbcan$HMMProfile, levels = c(
  "GH110","GH85","GH28","GH105","GH130",
  "PL1","PL10","PL29","PL12","PL33","PL40",
  "PL6","PL7","PL17","PL38"))

# Reshape 
dbcan <- reshape2::dcast(
  dbcan, Strain ~ HMMProfile)
rownames(dbcan) = dbcan$Strain
dbcan$Strain = NULL

pheatmap(
  dbcan, 
  cellwidth = 15, 
  cellheight = 15,
  show_rownames = T, 
  show_colnames = T, 
  fontsize = 11,
  na_col = "white",
  cluster_rows = F, 
  cluster_cols = F,
  col = colorRampPalette(c(
    "#AEDBE7","#93CFE0","#43ABC9",
    "#1287a8","#3c6478","#0c374d",
    "#093145","gray11"))(80)) 

####################################################

## CAZYmes -- comparison with 76-1 ##

dbcanComp <- read.table(
  "dbcanComp.txt",
  h = T, 
  sep = "\t",                                               
  stringsAsFactors=F)

# Filter & subset 
dbcanComp = filter(
  dbcanComp, Coverage >= 0.8 & evalue < 1e-15)
dbcanComp$HMMProfile = gsub(
  dbcanComp$HMMProfile, pattern="_.*|.hmm", replacement = "")

# Import amino acid seqs
dbcanAA <- read.table(
  "dbcanComp.faa",
  h = T, 
  sep = "\t", 
  stringsAsFactors=F)

# Merge
dbcanComp <- merge(
  dbcanAA, dbcanComp, by="locus_tag", all.x=F)

# Export
write.table(
  dbcanComp, file = "dbcanComp_aa.txt", 
  sep ="\t", row.names=T)
  
############################

# Count for heatmap
dbcanComp = distinct(
  dbcanComp, HMMProfile, locus_tag, Strain) %>%
  group_by(HMMProfile, Strain) %>%
  filter(str_detect(HMMProfile, "PL|GH|CBM")) %>%
  tally()

# Subset to specific CAZymes; selected from test heatmap
dbcanComp = filter(dbcanComp, HMMProfile %in% c(
  "PL6","PL7","PL12","PL18",
  "PL22","PL24","PL25","PL29","PL33","PL40","GH10","GH95","GH141",
  "GH144", "CBM16","CBM32","CBM41","GH103","GH105",
   "GH106","GH109","GH113","GH39"))

# Reorder CAZymes
dbcanComp$HMMProfile <- factor(dbcanComp$HMMProfile, levels = c(
  "PL6","PL7","PL12","PL18","PL22","PL24","PL25",
  "PL29","PL33","PL40","GH10","GH95","GH109","GH113",
  "GH141","GH144","CBM16","CBM32","CBM41","GH103","GH105",
  "GH106","GH39"))

# Reshape 
dbcanComp <- reshape2::dcast(
  dbcanComp, HMMProfile ~ Strain)
rownames(dbcanComp) = dbcanComp$HMMProfile
dbcanComp$HMMProfile = NULL

pheatmap(
  dbcanComp, 
  cellwidth = 15, 
  cellheight = 15,
  show_rownames = T, 
  show_colnames = T, 
  fontsize = 11,
  na_col = "white",
  cluster_rows = F, 
  cluster_cols = F,
  col = rev(pnw_palette("Moth",4)))

######################################################

## GROWTH DATA ##

# Define summary function
summarySE <- function (
  data=NULL, measurevar, 
  groupvars=NULL, na.rm=T, conf.interval=0.95) {
  library(data.table)
  data <- data.table(data)
length2 <- function(x, na.rm=F) {
    if (na.rm) 
      sum(!is.na(x))
    else length(x)}
datac <- data[, .(
  lapply(.SD, length2, na.rm = na.rm), 
  lapply(.SD, mean, na.rm = na.rm),
  lapply(.SD, sd, na.rm = na.rm)),
  by = groupvars, .SDcols = measurevar]
names(datac) <- c(groupvars, "N", measurevar, "sd")
setkeyv(datac, groupvars)
datac[, se := unlist(sd) / sqrt(unlist(N))] 
ciMult <- qt(conf.interval / 2 + 0.5, unlist(datac$N) - 1)
datac[, ci := se * ciMult]
datac <- data.frame(datac)
return(datac)}

# Import growth / HPLC data
Growth <- read.table(
  "Growth.txt", sep="\t", header=T,
  stringsAsFactors=F, check.names=F)

# Melt + split 
Growth <- melt(Growth, id.vars = "hours") 
Growth = cSplit(Growth, "variable", "_")

# Rename columns
colnames(Growth)[colnames(Growth)=="variable_1"] <- "datatype"
colnames(Growth)[colnames(Growth)=="variable_2"] <- "substrate"
colnames(Growth)[colnames(Growth)=="variable_3"] <- "replicate"

# Add ID-column for later plotting
Growth = Growth %>%
  mutate(type = case_when(
    substrate=="Man"|substrate=="Gul" ~ "Alginate",
    substrate=="Glc" ~ "Glucose"),
    TRUE ~ "none") 
Growth$type = factor(
  Growth$type, levels=c("Alginate","Glucose"))
Growth$datatype = factor(
  Growth$datatype, levels=c("OD","HPLC"))

# Calculate mean + SD
Growth <- summarySE(
  Growth, measurevar="value", groupvars=c(
  "hours","datatype","substrate","type"), na.rm=T)

# Reformat
Growth$condition <- paste(Growth$datatype, Growth$substrate, sep="_")
Growth$value = as.numeric(Growth$value) 

##############################

plot1 <- ggplot(subset(
  Growth, condition %in% c("OD_Alg","HPLC_Man","HPLC_Gul")),
  aes(x=hours, y=value, group=condition)) +
geom_line(aes(color=condition),size=1.3) +
geom_point(aes(color=condition),size=3) +
geom_errorbar(
  aes(ymin=value-se, ymax=value+se, color=condition),  
  width=.3) +
scale_color_manual(values=c(
  "OD_Alg"="darkturquoise","HPLC_Man"="darkgreen",
  "HPLC_Gul"="aquamarine3")) +
facet_wrap(.~datatype, ncol=1, scales="free") +
theme_bw() +
theme(strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size=10),
      legend.position = "top")

plot2 <- ggplot(subset(
  Growth, condition %in% c("OD_Glc","HPLC_Glc")),
  aes(x=hours, y=value, group=condition)) +
geom_line(aes(color=condition),size=1.3) +
geom_point(aes(color=condition),size=3) +
geom_errorbar(
  aes(ymin=value-se, ymax=value+se, color=condition),  
  width=.3) +
scale_color_manual(values=c(
  "OD_Glc"="brown", 
  "HPLC_Glc"="darkorange")) +
facet_wrap(.~datatype, ncol=1, scales="free") +
theme_bw() +
theme(strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size=10),
      legend.position = "top")

plot_grid(
  plot1, 
  plot2, 
  ncol=2,
  rel_heights = c(0.9,1.1),
  rel_widths = c(1.1,1),
  align="h", 
  axis = "tblr")

######################################################

## GROWTH -- comparison with 76-1 ##

# Import growth / HPLC data
GrowthComp <- read.table(
  "GrowthComp.txt", sep="\t", header=T, 
  stringsAsFactors=F, check.names=F) 

# Melt + split 
GrowthComp <- melt(GrowthComp, id.vars = "hours") 
GrowthComp = cSplit(GrowthComp, "variable", "_")

# Rename columns
colnames(GrowthComp)[colnames(GrowthComp)=="variable_1"] <- "datatype"
colnames(GrowthComp)[colnames(GrowthComp)=="variable_2"] <- "strain"
colnames(GrowthComp)[colnames(GrowthComp)=="variable_3"] <- "replicate"

# Calculate mean + SD
GrowthComp <- summarySE(
  GrowthComp, measurevar="value", groupvars=c(
    "hours","datatype","strain"), na.rm=T)
GrowthComp$value = as.numeric(GrowthComp$value) 

ggplot(GrowthComp,
  aes(x=hours, y=value, group=strain)) +
geom_line(aes(color=strain),size=1.3) +
geom_point(aes(color=strain),size=3) +
geom_errorbar(
  aes(ymin=value-se, ymax=value+se, color=strain),  
  width=.3) +
scale_color_manual(values=c(
  "Alt761"="gray22",
  "Mar621"="darkturquoise")) +
theme_bw() +
theme(strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid.minor = element_blank())

######################################################

sessionInfo()
save.image("Maribacter_62-1.Rdata")
