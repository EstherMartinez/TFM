# ClinVar 

clinvar <- read.csv('./data/clinvar.csv', header = TRUE, sep = '\t')
clinvar$Name<- as.character(clinvar$Name)

protein <- unlist(strsplit(clinvar$Name, " "))
clinvar$Aminoacid_change <-protein[-seq(1,length(protein),2)] 
cDNA <- protein[seq(1,length(protein),2)]
cDNA <- unlist(strsplit(cDNA, ":"))
clinvar$cDNA_change <- cDNA[-seq(1,length(cDNA),2)]
grup_1_clinvar <- clinvar %>%
  filter(clinvar$Gene.s.=="PTPN11")
grup_2_clinvar <- clinvar %>%
  filter(clinvar$Gene.s.=="BRAF" | clinvar$Gene.s.=="RAF1")
grup_3_clinvar <- clinvar %>%
  filter(clinvar$Gene.s.=="HRAS" | clinvar$Gene.s.=="KRAS" | clinvar$Gene.s.=="NRAS")
grup_4_clinvar <- clinvar %>%
  filter(clinvar$Gene.s.=="MAP2K1" | clinvar$Gene.s.== "MAP2K2")
grup_5_clinvar <- clinvar %>%
  filter(clinvar$Gene.s.=="SOS1" | clinvar$Gene.s.=="SOS2")

# Gelb_2018 

estudis_funcionals_Gelb2018 <- read.csv('./data/Gelb2018_Estudis_funcionals.csv', header = TRUE, sep = '\t')

# GnomAD_v3.1

GnomAD <- read.csv('./data/gnomAD_v3.1.csv', header = TRUE, sep = ',')
