# Criteris ACMG_AMP RASopaties
Criteris_ACMG_AMP <- function(gen,variant_cDNA,variant_prot,de_novo,confirmats,estudis_funcionals,
                              casos_independents,individus_control,meiosis,computational_evidence_FM,trans_cis,
                              computational_evidence_FI,computational_evidence_ALT){
  
criteris <- data.frame(matrix(data = NA, nrow = 1, ncol = 7))
colnames(criteris) <- c('standalone_P','strong_P','moderate_P','supporting_P','standalone_B','strong_B','supporting_B')
rownames(criteris) <- c('RASopaties_suma')
  
# PS1_strong_P Base de dades clinvar

if (gen == "BRAF" | gen == "RAF1") {
  PS1_strong_P <- sum(ifelse(variant_prot %in% grup_2_clinvar$Aminoacid_change,1,0))
}
if (gen == "HRAS" | gen == "KRAS" | gen == "NRAS") {
  PS1_strong_P <- sum(ifelse(variant_prot %in% grup_3_clinvar$Aminoacid_change,1,0))
}
if (gen == "MAP2K1" | gen == "MAP2K2") {
  PS1_strong_P <- sum(ifelse(variant_prot %in% grup_4_clinvar$Aminoacid_change,1,0))
}
if (gen == "SOS1" | gen == "SOS2") {
  PS1_strong_P <- sum(ifelse(variant_prot %in% grup_5_clinvar$Aminoacid_change,1,0))
}

# PS2 / PM6 Pregunta a l'usuari si la variant es presenta de novo en un pacient amb fenotip i si hi ha confirmació parental

PS2_standalone_P <- sum(ifelse(de_novo >= 2 & confirmats >=2 | de_novo == 3 & confirmats >= 1, 1, 0))
PS2_strong_P <- sum(ifelse(de_novo == 1 & confirmats >=1, 1, 0))
PS2_supporting_P <- sum(ifelse(de_novo == 1 & confirmats == 0, 1, 0))
PM6_standalone_P <- sum(ifelse(de_novo >=4 & confirmats >= 0, 1, 0))
PM6_strong_P <- sum(ifelse(de_novo == 2 & confirmats == 0 | de_novo == 3 & confirmats == 0, 1, 0))
PM6_supporting_P <- sum(ifelse(de_novo == 1 & confirmats == 0, 1, 0))

# PS3 / BS3 Es pregunta a l'usuari si hi ha estudis funcionals in vitro o in vivo describint si la variant té efectes "malignes" sobre el gen
# També consulta la informació obtinguda de la literatura en concret la supplementary table S6 de Gelb et al 2018

if (gen == "MAP2K2" | gen == "BRAF" | gen == "HRAS" | gen == "KRAS" | gen == "PTPN11" | gen == "RAF1" | gen == "SHOC2" | gen == "SOS1") {
  protein_query <- str_sub(variant_prot, start=2, end=-2)
  estudis_funcionals_literatura <- sum(ifelse(protein_query %in% estudis_funcionals_Gelb2018$protein_change,1,0))
}

PS3_strong_P <- sum(ifelse(estudis_funcionals > 0, 1, 0),estudis_funcionals_literatura)
BS3_strong_B <- sum(ifelse(PS3_strong_P == 0, 1, 0))

# PS4 / BS2 Es pregunta per la prevalença de la variant

PS4_strong_P <- sum(ifelse(casos_independents >= 5, 1, 0))
PS4_moderate_P <- sum(ifelse(casos_independents == 3 | casos_independents == 4, 1, 0))
PS4_supporting_P <- sum(ifelse(casos_independents == 1 | casos_independents == 2, 1, 0))
BS2_strong_B <- sum(ifelse(individus_control >= 3, 1, 0))

# PM1 Mira si la variant es troba en un hotspot i/o en un domini funcional crític sense ser benigne

proteina <- unlist(strsplit(variant_prot,"p.")) 
posicio_aa <- as.numeric(str_sub(proteina[2],4,5))  

if (gen == "BRAF" | gen == "RAF1") {
  # Dominis funcionals dels gens del grup 2
  BRAF_df <- c(157:200,201:261)
  RAF1_df <- c(58:101,105:185)
  # Regions hotspot  
  BRAF_hotspot <-c (238:286,439:474,531,459:474,594:627)
  RAF1_hotspot <- c(251:266)
  dominis <- sum(ifelse(posicio_aa == BRAF_df, 1, 0), ifelse(posicio_aa == RAF1_df, 1, 0), ifelse(posicio_aa == BRAF_hotspot, 1, 0), ifelse(posicio_aa == RAF1_hotspot, 1, 0))
}

if (gen == "HRAS" | gen == "KRAS" | gen == "NRAS") {
  # Dominis funcionals dels gens del grup 3
  HRAS_df <- c(32:40,10:17,57:61,116:119)
  KRAS_df <- c(32:38,10:18,59:60,116:119)
  NRAS_df <- c(0:0,10:18,57:61,116:119)  
  # Regions hotspot 
  HRAS_hotspot <- c(12,13,14,58,59,60,61,62,63)
  KRAS_hotspot <- c(12,13,14,58,59,60,61,62,63)
  dominis <- sum(ifelse(posicio_aa == HRAS_df, 1, 0), ifelse(posicio_aa == KRAS_df, 1, 0), ifelse(posicio_aa == NRAS_df, 1, 0), ifelse(posicio_aa == HRAS_hotspot, 1, 0), ifelse(posicio_aa == KRAS_hotspot, 1, 0))
}

if (gen == "MAP2K2" | gen == "MAP2K1") {
  # Dominis funcionals dels gens del grup 4
  MAP2K1_df <- c(32:44,44:51,74:82,143:146,192:195,208:233,262:307,362:396)
  MAP2K2_df <- c(36:48,48:55,78:96,147:150,196:199,212:237,266:311,370:400)
  # Regions hotspot  
  MAP2K1_hotspot <- c(43:61,124:134)
  MAP2K2_hotspot <- c(47:65,128:138)
  dominis <- sum(ifelse(posicio_aa == MAP2K1_df, 1, 0), ifelse(posicio_aa == MAP2K2_df, 1, 0), ifelse(posicio_aa == MAP2K1_hotspot, 1, 0), ifelse(posicio_aa == MAP2K2_hotspot, 1, 0))
}

if (gen == "SOS1" | gen == "SOS2") {
  # Dominis funcionals dels gens del grup 5
  SOS1_df <- c(200:390,444:548,597:741,780:1019)
  SOS2_df <- c(198:388,442:546,595:739,778:1017)
  # Regions hotspot  
  SOS1_hotspot <- c(269,552,420:500)
  dominis <- sum(ifelse(posicio_aa == SOS1_df, 1, 0), ifelse(posicio_aa == SOS2_df, 1, 0), ifelse(posicio_aa == SOS1_hotspot, 1, 0))
}

if (gen == "PTPN11") {
  PTPN11_hotspot <- c(308,4,7,8,9,58:63,69:77,247,251,255,256,258,261,265,278:281,284)
  dominis <- sum(ifelse(posicio_aa == PTPN11_hotspot, 1, 0))
}
if (gen == "SHOC2") {
  SHOC2_hotspot <- c(2)
  dominis <- sum(ifelse(posicio_aa == SHOC2_hotspot, 1, 0))  
}
if (gen == "SPRED1") {
  SPRED1_hotspot <- c(6:123, 233-285, 334:442)
  dominis <- sum(ifelse(posicio_aa == SPRED1_hotspot, 1, 0))    
}
if (gen == "NF1") {
  NF1_hotspot <- c(543:909,1095:1176, 1170:1218, 1198:1530, 1510:1570, 1560:1696,1713:1816, 2619:2719)
  dominis <- sum(ifelse(posicio_aa == NF1_hotspot, 1, 0))
}

PM1_supporting_P <- ifelse(dominis > 0, 1, 0)

# PM2 / BA1 / BS1 Les dades de freqüència al·lèlica les consulta de GnomAD_v3.1

AF <- GnomAD[GnomAD$Gene == gen & GnomAD$Transcript.Consequence == variant_cDNA,]$Allele.Frequency
PM2_supporting_P <- sum(ifelse(AF == 0 | AF < 0.00025, 1, 0))
BA1_standalone_B <- sum(ifelse(AF >= 0.0005, 1, 0))
BS1_strong_B <- sum(ifelse(AF >= 0.00025 & AF < 0.0005, 1, 0))

# PM4 / BP3 Només fem variants missense 
PM4_supporting_P <- 0
BP3_supporting_B <- 0

# PM5 Base de dades ClinVar

residu_aa <- str_sub(proteina[2],1,5)

PM5_strong_P <- sum(ifelse(length(grep(residu_aa, clinvar$Aminoacid_change)) >= 2, 1, 0))

PM5_supporting_P <- sum(ifelse(PM5_strong_P == 0 & length(grep(residu_aa, clinvar$Aminoacid_change)) > 0, 1, 0))

# PP1 / BS4 Es pregunta a l'usuari sobre la cosegregació de la variant amb el fenotip

PP1_strong_P <- sum(ifelse( meiosis >= 7, 1, 0))
PP1_moderate_P <- sum(ifelse( meiosis == 5 | meiosis == 6, 1, 0))
PP1_supporting_P <- sum(ifelse( meiosis == 3 | meiosis == 4, 1, 0))
BS4_strong_B <- sum(ifelse(meiosis == 0, 1, 0))

# PP2 / BP1_Només per SPRED1 i NF1 només fem variants missense 

PP2_supporting_P <- ifelse(gen == 'NF1' | gen == 'SPRED1', 0, 1)
BP1_supporting_B <- ifelse(gen == 'NF1' | gen == 'SPRED1', 1, 0)

# PP3 Es pregunta a l'usuari si hi ha evidències computacionals que apuntin cap a un efecte perjudicial en el gen

PP3_supporting_P <- sum(ifelse(computational_evidence_FM > 0, 1, 0))

# BP2 Es pregunta a l'usuari si la variant es troba en trans amb una variant patogènica d'un gen amb penetrancia completa dominant 
# o si es troba en cis amb una variant patogènica amb qualsevol tipus d'herencia

BP2_supporting_B <- sum(ifelse(trans_cis > 0, 1, 0))

# BP4 Es pregunta a l'usuari si hi ha evidències computacionals que apuntin cap a un efecte que no compromet la funció de la proteïna

BP4_supporting_B <- sum(ifelse(computational_evidence_FI > 0, 1, 0))

# BP5 Es pregunta a l'usuari si hi ha evidències computacionals d'alteracions genètiques alternatives

BP5_supporting_B <- ifelse(computational_evidence_ALT > 0, 1, 0)

criteris$standalone_P <- PS2_standalone_P + PM6_standalone_P
criteris$strong_P <- PS1_strong_P + PS2_strong_P + PS3_strong_P + PS4_strong_P + PM5_strong_P + PM6_strong_P + PP1_strong_P
criteris$moderate_P <- PS4_moderate_P + PP1_moderate_P
criteris$supporting_P <- PS2_supporting_P + PS4_supporting_P + PM1_supporting_P + PM2_supporting_P + PM4_supporting_P +
  PM5_supporting_P + PM6_supporting_P + PP1_supporting_P + PP2_supporting_P + PP3_supporting_P
criteris$standalone_B <- BA1_standalone_B
criteris$strong_B <- BS1_strong_B + BS2_strong_B + BS3_strong_B + BS4_strong_B
criteris$supporting_B <- BP1_supporting_B + BP2_supporting_B + BP3_supporting_B + BP4_supporting_B + BP5_supporting_B

return(criteris)
}


# Classificació ACMG-AMP-RASopaties

Classificacio_ACMG_AMP <- function(criteris){
  
ACMG_AMP_Clas <- data.frame(matrix(data = NA, nrow = 1, ncol = 5))
colnames(ACMG_AMP_Clas) <- c('Pathogenic', 'Likely_Pathogenic', 'Benign', 'Likely_Benign', 'VUS')

ACMG_AMP_Clas$Pathogenic[criteris$standalone_P >= 1 & criteris$strong_P >=1 |
                           criteris$standalone_P >= 1 & criteris$moderate_P >= 2 |
                           criteris$standalone_P >= 1 & criteris$moderate_P == 1 & criteris$supporting_P == 1 |
                           criteris$standalone_P >= 1 & criteris$supporting_P >= 2 |
                           criteris$strong_P >= 2 |
                           criteris$strong_P == 1 & criteris$moderate_P >= 3 |
                           criteris$strong_P == 1 & criteris$moderate_P == 2 & criteris$supporting_P >= 2 |
                           criteris$strong_P == 1 & criteris$moderate_P == 1 & criteris$supporting_P >= 4] <- 1

ACMG_AMP_Clas$Likely_Pathogenic[criteris$standalone_P == 1 & criteris$moderate_P == 1 |
                                  criteris$strong_P == 1 & criteris$moderate_P == 1 |
                                  criteris$strong_P == 1 & criteris$moderate_P == 2 |
                                  criteris$strong_P == 1 & criteris$supporting_P >= 2 |
                                  criteris$moderate_P >= 3 |
                                  criteris$moderate_P == 2 & criteris$supporting_P >= 2 |
                                  criteris$moderate_P == 1 & criteris$supporting_P >= 4] <- 1
ACMG_AMP_Clas$Benign[criteris$standalone_B == 1 |
                       criteris$strong_B >= 2] <- 1
ACMG_AMP_Clas$Likely_Benign[criteris$strong_B == 1 & criteris$supporting_B == 1 |
                              criteris$supporting_B >= 2] <- 1
ACMG_AMP_Clas$VUS[ACMG_AMP_Clas$Pathogenic == 0 &
                    ACMG_AMP_Clas$Likely_Pathogenic == 0 &
                    ACMG_AMP_Clas$Benign == 0 &
                    ACMG_AMP_Clas$Likely_Benign == 0 |
                    ACMG_AMP_Clas$Pathogenic == 1 & ACMG_AMP_Clas$Benign == 1 |
                    ACMG_AMP_Clas$Pathogenic == 1 & ACMG_AMP_Clas$Likely_Benign == 1 |
                    ACMG_AMP_Clas$Likely_Pathogenic == 1 & ACMG_AMP_Clas$Benign == 1 |
                    ACMG_AMP_Clas$Likely_Pathogenic == 1 & ACMG_AMP_Clas$Likely_Benign ==1] <- 1

return(ACMG_AMP_Clas) 
}
classificacio <- Classificacio_ACMG_AMP(criteris)
classificacio