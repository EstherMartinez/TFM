
# Instal·lació de paquets o llibreries
libraries <- c("stringr", "dplyr")
check.libraries <- is.element(libraries, installed.packages()[, 1])==FALSE
libraries.to.install <- libraries[check.libraries]
if (length(libraries.to.install!=0)) {
  install.packages(libraries.to.install)
}
library(stringr)
library(dplyr)

# Bases de dades

source(file = "./scripts/02_bases de dades.R")

# Criteris RASopaties

source(file = "./scripts/04_funcions.R")

criteris <- Criteris_ACMG_AMP(gen='MAP2K2',variant_cDNA = 'c.171T>A',variant_prot = '(p.Phe57Leu)',
                              de_novo = 0,confirmats = 5,estudis_funcionals = 1,
                              casos_independents = 3,individus_control = 0,meiosis = 7,
                              computational_evidence_FM = 2,trans_cis = 0,
                              computational_evidence_FI = 0,computational_evidence_ALT = 0)
# Classificació RASopaties

classificacio <- Classificacio_ACMG_AMP(criteris)

