sink("Rscript.txt")

library(remotes)
library(dplyr)
library(TwoSampleMR)
library(data.table)



# ASAT --------------------------------------------------------------------


exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

#Fasting insulin
exposure_dat2 <- read_exposure_data(
  filename = "data/FI_GWAS/FI_combined_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)

snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")

exposure_dat2 <- read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/FI_GWAS/FI_combined_1000G_density_formatted_21-03-29.txt",
  sep = "\t",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "a2",
  other_allele_col = "a1",
  samplesize_col = "n"
)

dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006464.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_FI_EC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_FI_EC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_FI","se_FI","beta_EC","se_EC")

ASAT_FI_EC$SNP<-dat1$SNP
ASAT_FI_EC$effect_alllele<-dat1$effect_allele.exposure
ASAT_FI_EC$other_allele<-dat1$other_allele.exposure
ASAT_FI_EC$beta_ASAT<-dat1$beta.exposure
ASAT_FI_EC$se_ASAT<-dat1$se.exposure
ASAT_FI_EC$beta_FI<-dat1$beta.outcome
ASAT_FI_EC$se_FI<-dat1$se.outcome
ASAT_FI_EC$beta_EC<-dat2$beta.outcome
ASAT_FI_EC$se_EC<-dat2$se.outcome

fwrite(ASAT_FI_EC,"data/LiverFatMR/MVMR_files/ASAT_FI_EC.csv")

#Endometrioid endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006465.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_FI_EEC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_FI_EEC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_FI","se_FI","beta_EEC","se_EEC")

ASAT_FI_EEC$SNP<-dat1$SNP
ASAT_FI_EEC$effect_alllele<-dat1$effect_allele.exposure
ASAT_FI_EEC$other_allele<-dat1$other_allele.exposure
ASAT_FI_EEC$beta_ASAT<-dat1$beta.exposure
ASAT_FI_EEC$se_ASAT<-dat1$se.exposure
ASAT_FI_EEC$beta_FI<-dat1$beta.outcome
ASAT_FI_EEC$se_FI<-dat1$se.outcome
ASAT_FI_EEC$beta_EEC<-dat2$beta.outcome
ASAT_FI_EEC$se_EEC<-dat2$se.outcome

fwrite(ASAT_FI_EEC,"data/LiverFatMR/MVMR_files/ASAT_FI_EEC.csv")


#SHBG
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/IEUdownMVMR/ieu-b-4870_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)


snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)

exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <-read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/IEUdownMVMR/ieu-b-4870.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006464.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_SHBG_EC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_SHBG_EC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_SHBG","se_SHBG","beta_EC","se_EC")

ASAT_SHBG_EC$SNP<-dat1$SNP
ASAT_SHBG_EC$effect_alllele<-dat1$effect_allele.exposure
ASAT_SHBG_EC$other_allele<-dat1$other_allele.exposure
ASAT_SHBG_EC$beta_ASAT<-dat1$beta.exposure
ASAT_SHBG_EC$se_ASAT<-dat1$se.exposure
ASAT_SHBG_EC$beta_SHBG<-dat1$beta.outcome
ASAT_SHBG_EC$se_SHBG<-dat1$se.outcome
ASAT_SHBG_EC$beta_EC<-dat2$beta.outcome
ASAT_SHBG_EC$se_EC<-dat2$se.outcome

fwrite(ASAT_SHBG_EC,"data/LiverFatMR/MVMR_files/ASAT_SHBG_EC.csv")

#Endometrioid endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006464.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_SHBG_EEC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_SHBG_EEC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_SHBG","se_SHBG","beta_EEC","se_EEC")

ASAT_SHBG_EEC$SNP<-dat1$SNP
ASAT_SHBG_EEC$effect_alllele<-dat1$effect_allele.exposure
ASAT_SHBG_EEC$other_allele<-dat1$other_allele.exposure
ASAT_SHBG_EEC$beta_ASAT<-dat1$beta.exposure
ASAT_SHBG_EEC$se_ASAT<-dat1$se.exposure
ASAT_SHBG_EEC$beta_SHBG<-dat1$beta.outcome
ASAT_SHBG_EEC$se_SHBG<-dat1$se.outcome
ASAT_SHBG_EEC$beta_EEC<-dat2$beta.outcome
ASAT_SHBG_EEC$se_EEC<-dat2$se.outcome

fwrite(ASAT_SHBG_EEC,"data/LiverFatMR/MVMR_files/ASAT_SHBG_EEC.csv")

#Bioavailable testosterone (female)
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/IEUdownMVMR/ieu-b-4869_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <-read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/IEUdownMVMR/ieu-b-4869.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006464.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_BT_EC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_BT_EC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_TT","se_TT","beta_OC","se_OC")

ASAT_BT_EC$SNP<-dat1$SNP
ASAT_BT_EC$effect_alllele<-dat1$effect_allele.exposure
ASAT_BT_EC$other_allele<-dat1$other_allele.exposure
ASAT_BT_EC$beta_ASAT<-dat1$beta.exposure
ASAT_BT_EC$se_ASAT<-dat1$se.exposure
ASAT_BT_EC$beta_TT<-dat1$beta.outcome
ASAT_BT_EC$se_TT<-dat1$se.outcome
ASAT_BT_EC$beta_OC<-dat2$beta.outcome
ASAT_BT_EC$se_OC<-dat2$se.outcome

fwrite(ASAT_BT_EC,"data/LiverFatMR/MVMR_files/ASAT_BT_EC.csv")

#Endometrioid endometrial cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006464.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_BT_EEC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_BT_EEC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_FI","se_FI","beta_EEC","se_EEC")

ASAT_BT_EEC$SNP<-dat1$SNP
ASAT_BT_EEC$effect_alllele<-dat1$effect_allele.exposure
ASAT_BT_EEC$other_allele<-dat1$other_allele.exposure
ASAT_BT_EEC$beta_ASAT<-dat1$beta.exposure
ASAT_BT_EEC$se_ASAT<-dat1$se.exposure
ASAT_BT_EEC$beta_FI<-dat1$beta.outcome
ASAT_BT_EEC$se_FI<-dat1$se.outcome
ASAT_BT_EEC$beta_EEC<-dat2$beta.outcome
ASAT_BT_EEC$se_EEC<-dat2$se.outcome

fwrite(ASAT_BT_EEC,"data/LiverFatMR/MVMR_files/ASAT_BT_EEC.csv")

#CXCL8
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/CXCL8/CXCL8_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)

snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <- read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/Ferkingstad_et_al_2021_proteins/CXCL8/3447_64_CXCL8_IL_8.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)

dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Ovarian cancer
outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ieu-a-1120.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_CXC_OC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_CXC_OC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_CXC","se_CXC","beta_OC","se_OC")

ASAT_CXC_OC$SNP<-dat1$SNP
ASAT_CXC_OC$effect_alllele<-dat1$effect_allele.exposure
ASAT_CXC_OC$other_allele<-dat1$other_allele.exposure
ASAT_CXC_OC$beta_ASAT<-dat1$beta.exposure
ASAT_CXC_OC$se_ASAT<-dat1$se.exposure
ASAT_CXC_OC$beta_CXC<-dat1$beta.outcome
ASAT_CXC_OC$se_CXC<-dat1$se.outcome
ASAT_CXC_OC$beta_OC<-dat2$beta.outcome
ASAT_CXC_OC$se_OC<-dat2$se.outcome

fwrite(ASAT_CXC_OC,"data/LiverFatMR/MVMR_files/ASAT_CXC_OC.csv")

#IGFBP
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/IGFBP_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)

snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <- read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/2771_35_IGFBP1_IGFBP_1.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)

dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Esophagus adenocarcinoma
outcome_dat <- read_outcome_data(
  snps = dat1$SNP,
  filename = "data/Obesity_related_cancer_GWAS/Esophagus_adenocarcinoma/GCST003739_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele"
)


dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_IGFBP_OAC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_IGFBP_OAC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_IGFBP","se_IGFBP","beta_OAC","se_OAC")

ASAT_IGFBP_OAC$SNP<-dat1$SNP
ASAT_IGFBP_OAC$effect_alllele<-dat1$effect_allele.exposure
ASAT_IGFBP_OAC$other_allele<-dat1$other_allele.exposure
ASAT_IGFBP_OAC$beta_ASAT<-dat1$beta.exposure
ASAT_IGFBP_OAC$se_ASAT<-dat1$se.exposure
ASAT_IGFBP_OAC$beta_IGFBP<-dat1$beta.outcome
ASAT_IGFBP_OAC$se_IGFBP<-dat1$se.outcome
ASAT_IGFBP_OAC$beta_OAC<-dat2$beta.outcome
ASAT_IGFBP_OAC$se_OAC<-dat2$se.outcome

fwrite(ASAT_IGFBP_OAC,"data/LiverFatMR/MVMR_files/ASAT_IGFBP_OAC.csv")

#HDL cholesterol
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/IEUdownMVMR/ieu-b-109_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <-read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/IEUdownMVMR/ieu-b-109.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)


dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

#Triple negative breast cancer
outcome_dat <- read_outcome_data(
  snps = dat1$SNP,
  filename = "data/Obesity_related_cancer_GWAS/Breast_cancer/annotated/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt",
  sep = ",",
  snp_col = "rsID",
  beta_col = "Triple_Neg_log_or_meta",
  se_col = "Triple_Neg_se_meta",
  effect_allele_col = "Effect.Meta",
  other_allele_col = "Baseline.Meta",
  eaf_col = "Freq.Gwas"
)

dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

ASAT_HDL_TNBC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(ASAT_HDL_TNBC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_HDL","se_HDL","beta_TNBC","se_TNBC")

ASAT_HDL_TNBC$SNP<-dat1$SNP
ASAT_HDL_TNBC$effect_alllele<-dat1$effect_allele.exposure
ASAT_HDL_TNBC$other_allele<-dat1$other_allele.exposure
ASAT_HDL_TNBC$beta_ASAT<-dat1$beta.exposure
ASAT_HDL_TNBC$se_ASAT<-dat1$se.exposure
ASAT_HDL_TNBC$beta_HDL<-dat1$beta.outcome
ASAT_HDL_TNBC$se_HDL<-dat1$se.outcome
ASAT_HDL_TNBC$beta_TNBC<-dat2$beta.outcome
ASAT_HDL_TNBC$se_TNBC<-dat2$se.outcome

fwrite(ASAT_HDL_TNBC,"data/LiverFatMR/MVMR_files/ASAT_HDL_TNBC.csv")


# VAT ---------------------------------------------------------------------


# GFAT --------------------------------------------------------------------


#Adiponectin
exposure_dat1 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat2 <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/Adiponectin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
snplist<-c(exposure_dat1$SNP,exposure_dat2$SNP)

exposure_dat1 <- read_outcome_data(
  snps = snplist,
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT.allSNPs",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat1<-data.frame(exposure_dat1)
exposure_dat1 <- format_data(exposure_dat1,
                             type="exposure",
                             snp_col="SNP",
                             beta_col = "beta.outcome",
                             se_col="se.outcome",
                             effect_allele_col = "effect_allele.outcome",
                             other_allele_col = "other_allele.outcome",
                             eaf_col = "eaf.outcome")


exposure_dat2 <- read_outcome_data(
  snps = exposure_dat1$SNP,
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/3554_24_ADIPOQ_Adiponectin.txt",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  eaf_col = "ImpMAF"
)
dat1 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = exposure_dat2,action=3)
dat1<-dat1[dat1$mr_keep=="TRUE",]

outcome_dat <-read_outcome_data(
  snps = dat1$SNP,
  filename = "data/IEUdownMVMR/ebi-a-GCST006466.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
dat2 <- harmonise_data(exposure_dat = exposure_dat1, outcome_dat = outcome_dat)
dat2<-dat2[dat2$mr_keep=="TRUE",]

dat1<-dat1[dat1$SNP %in% dat2$SNP,]

GFAT_AD_NEC<-data.frame(matrix(ncol=9,nrow=nrow(dat1)))
colnames(GFAT_AD_NEC)<-c("SNP","effect_alllele","other_allele","beta_ASAT","se_ASAT","beta_FI","se_FI","beta_EEC","se_EEC")

GFAT_AD_NEC$SNP<-dat1$SNP
GFAT_AD_NEC$effect_alllele<-dat1$effect_allele.exposure
GFAT_AD_NEC$other_allele<-dat1$other_allele.exposure
GFAT_AD_NEC$beta_ASAT<-dat1$beta.exposure
GFAT_AD_NEC$se_ASAT<-dat1$se.exposure
GFAT_AD_NEC$beta_FI<-dat1$beta.outcome
GFAT_AD_NEC$se_FI<-dat1$se.outcome
GFAT_AD_NEC$beta_EEC<-dat2$beta.outcome
GFAT_AD_NEC$se_EEC<-dat2$se.outcome

fwrite(GFAT_AD_NEC,"data/LiverFatMR/MVMR_files/GFAT_AD_NEC.csv")

sink()
