# Script classificació de la variant c.171T>A (p.Phe57Leu) del gen MAK2K2
# L'usuari introdueix el nom del gen
print("Introdueix el nom del gen")
gen <- as.character(scan(n=1, what="character"))

# L'usuari introdueix la variant en format cDNA
print("Introdueix la variant en format cDNA")
variant_cDNA <- as.character(scan(n=1, what="character"))
# L'usari introdueix la variant en format proteïna
print("Introdueix la variant en format proteïna")
variant_prot <- as.character(scan(n=1, what="character"))

# Instal·lació de paquets
library(stringr)
# Criteris RASopaties

# PS1_strong_P (script) de moment només mira pel grup 4 base de dades clinvar

clinvar <- read.csv('../clinvar_result_MAP2K1_MAP2K2.csv', header = TRUE, sep = ';')

PS1_strong_P <- sum(ifelse(clinvar$Aminoacid_change == variant_prot, 1, 0))

# PS2 / PM6 (usuari)

# Es preguntarà a l'usuari si la variant es presenta de novo en un pacient amb fenotip i si hi ha confirmació parental

print("Quants casos de novo s'han trobat?")
de_novo <- as.integer(scan(n=1, what="integer"))
print("Quants familiars hi ha confirmats?")
confirmats <- as.integer((scan(n=1, what="integer")))

PS2_standalone_P <- sum(ifelse(de_novo >= 2 & confirmats >=2 | de_novo == 3 & confirmats >= 1, 1, 0))

PS2_strong_P <- sum(ifelse(de_novo == 1 & confirmats >=1, 1, 0))

PS2_supporting_P <- sum(ifelse(de_novo == 1 & confirmats == 0, 1, 0))

PM6_standalone_P <- sum(ifelse(de_novo >=4 & confirmats >= 0, 1, 0))

PM6_strong_P <- sum(ifelse(de_novo == 2 & confirmats == 0 | de_novo == 3 & confirmats == 0, 1, 0))

PM6_supporting_P <- sum(ifelse(de_novo == 1 & confirmats == 0, 1, 0))

# PS3 / BS3 (usuari)

# Es preguntarà a l'usuari si hi ha estudis funcionals in vitro o in vivo describint si la variant té efectes "malignes" sobre el gen
# Més endavant es pot consultar la supplementary table S6 de Gelb et al 2018

print("Hi ha estudis funcionals in vitro o in vivo, quants?")
estudis_funcionals <- as.integer((scan(n=1, what="integer")))
PS3_strong_P <- sum(ifelse(estudis_funcionals > 0, 1, 0))
BS3_strong_B <- sum(ifelse(PS3_strong_P == 0, 1, 0))

# PS4 / BS2 (usuari?)

print("La prevalencia de la variant en individus afectats és significativament elevada en comparació a individus controls?")
casos_independents <- as.integer((scan(n=1, what = "integer")))
PS4_strong_P <- sum(ifelse(casos_independents >= 5, 1, 0))
PS4_moderate_P <- sum(ifelse(casos_independents == 3 | casos_independents == 4, 1, 0))
PS4_supporting_P <- sum(ifelse(casos_independents == 1 | casos_independents == 2, 1, 0))
individus_control <- as.integer((scan(n=1, what = "integer")))
BS2_strong_B <- sum(ifelse(individus_control >= 3, 1, 0))

# PM1 (script)

# Per començar es mira si la variant es troba en un hotspot i/o en un domini funcional crític sense ser benigne del grup 4 (MAP2K1-MAP2K2)
# Dominis funcionals dels gens del grup 4
MAP2K1_df <- c(32:44,44:51,74:82,143:146,192:195,208:233,262:307,362:396)
MAP2K2_df <- c(36:48,48:55,78:96,147:150,196:199,212:237,266:311,370:400)

# Regions hotspot dels gens del grup 4

MAP2K1_hotspot <- c(43:61,124:134)
MAP2K2_hotspot <- c(47:65,128:138)

proteina <- unlist(strsplit(variant_prot,"p.")) 
posicio_aa <- as.numeric(str_sub(proteina[2],4,5))  

dominis <- sum(ifelse(posicio_aa == MAP2K1_df, 1, 0), ifelse(posicio_aa == MAP2K2_df, 1, 0), ifelse(posicio_aa == MAP2K1_hotspot, 1, 0), ifelse(posicio_aa == MAP2K2_hotspot, 1, 0))
PM1_supporting_P <- ifelse(dominis > 0, 1, 0)

# PM2 / BA1 / BS1 (script) Per ara preguntar a l'usuari les dades de freqüència al·lèlica, més endavant consultar: Exome Seq Project, 1000 Genomes, ExAC

print("Quina és la freqüència al·lèlica de la variant en la població?")

frequencia_alel <- as.numeric(scan(n=1, what="numeric"))

PM2_supporting_P <- sum(ifelse(frequencia_alel == 0, 1, 0))
BA1_standalone_B <- sum(ifelse(frequencia_alel >= 0.0005, 1, 0))
BS1_strong_B <- sum(ifelse(frequencia_alel >= 0.00025, 1, 0))

# PM4 (script) Per ara preguntar a l'usuari si el canvi en la proteïna es deu a una deleció o inserció in-frame en una regió no repetitiva o stop-loss variants
print("Deleció o inserció in-frame en una regió no repetitiva o variant stop-loss?")
prot_length <- as.integer(scan(n=1, what = "integer")) # Si és que sí posar un 1, si és que no posar un 0
PM4_supporting_P <- sum(ifelse(prot_length == 1, 1, 0))

# PM5 (script)
# Es busca un canvi missense en un residu aa on un altre canvi missense s'ha descrit com a patogènic (per ara consulta la base de dades ClinVar)

residu_aa <- str_sub(proteina[2],1,5)

PM5_strong_P <- sum(ifelse(length(grep(residu_aa, clinvar$Aminoacid_change)) >= 2, 1, 0))

PM5_supporting_P <- sum(ifelse(PM5_strong_P == 0 & length(grep(residu_aa, clinvar$Aminoacid_change)) > 0, 1, 0))


# PP1 (usuari)

# Es preguntarà a l'usuari sobre la cosegregació de la variant amb el fenotip

print("Hi ha cosegregació de la variant amb el fenotip? Si és que sí indica-ho amb un 1, si és que no amb un 0. Indica el nombre de meiosis")
meiosis <- as.integer((scan(n=2, what="integer")))
PP1_strong_P <- sum(ifelse(meiosis[1] == 1 & meiosis[2] >= 7, 1, 0))
PP1_moderate_P <- sum(ifelse(meiosis[1] == 1 & meiosis[2] == 5 | meiosis == 6, 1, 0))
PP1_supporting_P <- sum(ifelse(meiosis[1] == 1 & meiosis[2] == 3 | meiosis == 4, 1, 0))

BS4_strong_B <- sum(ifelse(meiosis[1] == 0 & meiosis[2] >= 1, 1, 0))

# PP2 (usuari?)

# Per ara es pregunta a l'usuari si la variant és una missense comú en la malaltia
print("És una variant missense comú per la malaltia?")
low_rate_B <- as.integer((scan(n=1, what="integer"))) # Si és que sí posa 1, si és que no 0

PP2_supporting_P <- sum(ifelse(low_rate_B == 1, 1, 0))
# PP3 / BP4 (usuari?)

# Per ara es pregunta a l'usuari si hi ha evidències computacionals que apuntin cap a un efecte perjudicial en el gen
print("Hi ha evidències computacionals que descriguin efectes en el gen provocats per la variant, quantes?")
computational_evidence <- as.integer((scan(n=1, what="integer")))
PP3_supporting_P <- sum(ifelse(computational_evidence > 0, 1, 0))

BP4_supporting_B <- sum(ifelse(PP3_supporting_P == 0, 1, 0))

# BP7?
BP7 <- 0
# BP1 (script)
# Per ara es pregunta a l'usuari si la variant és una missense que afecta a gens en els quals variants truncants provoquen la malaltia
print("És una variant missense en un gen on truncaments provoquen la malaltia?")
GOF <- as.integer((scan(n=1, what = "integer"))) # Si és que sí posa 1, sinó 0

BP1_supporting_B <- sum(ifelse(GOF > 0, 1, 0))

# BP2 (script)
# Per ara es pregunta a l'usuari si la variant es troba en trans amb una variant patogènica d'un gen amb penetrancia completa dominant 
# o si es troba en cis amb una variant patogènica amb qualsevol tipus d'herencia
print("Variant en trans o cis?")
trans_cis <- as.integer((scan(n=1, what = "integer"))) # Es posa 1 si és en trans o cis i 0 si no es troba en cap d'aquestes

BP2_supporting_B <- sum(ifelse(trans_cis > 0, 1, 0))

# BP3 (script)
# De moment només fem variants missense pels gens que no són NF1 o SPRED1
BP3_supporting_B <- 0

# BP5 (script)
# variant trobada en casos amb bases moleculars alternatives per la malaltia
BP5_supporting_B <- 0

# BP7 (script)
# synonymous variant? es pot definir junt a BP4



# Criteris ACMG_AMP RASopaties

criteris <- data.frame(matrix(data = NA, nrow = 1, ncol = 7))
colnames(criteris) <- c('standalone_P','strong_P','moderate_P','supporting_P','standalone_B','strong_B','supporting_B')
rownames(criteris) <- c('RASopaties_suma')

criteris$standalone_P <- PS2_standalone_P + PM6_standalone_P
criteris$strong_P <- PS1_strong_P + PS2_strong_P + PS3_strong_P + PS4_strong_P + PM5_strong_P + PM6_strong_P + PP1_strong_P
criteris$moderate_P <- PS4_moderate_P + PP1_moderate_P
criteris$supporting_P <- PS2_supporting_P + PS4_supporting_P + PM1_supporting_P + PM2_supporting_P + PM4_supporting_P +
  PM5_supporting_P + PM6_supporting_P + PP1_supporting_P + PP2_supporting_P + PP3_supporting_P
criteris$standalone_B <- BA1_standalone_B
criteris$strong_B <- BS1_strong_B + BS2_strong_B + BS3_strong_B + BS4_strong_B
criteris$supporting_B <- BP1_supporting_B + BP2_supporting_B + BP3_supporting_B + BP4_supporting_B + BP7 + BP5_supporting_B


# Classificació ACMG-AMP-RASopaties

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

  