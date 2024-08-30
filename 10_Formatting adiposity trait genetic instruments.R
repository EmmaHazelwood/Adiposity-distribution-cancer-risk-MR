library(MRInstruments)
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(gwasvcf)
gwasvcf::set_bcftools()


# Measures of adiposity ---------------------------------------------------
#Read in instruments
exposure_dat <- extract_instruments(outcomes="ebi-a-GCST90016673", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure[exposure_dat$exposure=="Percent liver fat || id:ebi-a-GCST90016673"]<-"Liver fat"

#Remove any with an F stat below 10
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]

#Save result
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


exposure_dat <- extract_instruments(outcomes="ebi-a-GCST90016675", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure[exposure_dat$exposure=="Pancreas fat || id:ebi-a-GCST90016675"]<-"Pancreas fat"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"VAT"
exposure_dat$samplesize.exposure<-38965
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"ASAT"
exposure_dat$samplesize.exposure<-38965
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"GFAT"
exposure_dat$samplesize.exposure<-38965
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

# Molecular traits --------------------------------------------------------

#Total testosterone female
exposure_dat <- extract_instruments(outcomes="ieu-b-4864", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
exposure_dat$exposure<-"Total testosterone female"
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Total testosterone male
exposure_dat <- extract_instruments(outcomes="ieu-b-4865", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
exposure_dat$exposure<-"Total testosterone male"
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#SHBG female
exposure_dat <- extract_instruments(outcomes="ieu-b-4870", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
exposure_dat$exposure<-"SHBG female"
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#SHBG male
exposure_dat <- extract_instruments(outcomes="ieu-b-4871", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
exposure_dat$exposure<-"SHBG male"
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Fasting insulin
exposure_dat <- read_exposure_data(
  filename = "data/FI_GWAS/FI_combined_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$samplesize.exposure<-98210
exposure_dat$exposure<-"Fasting insulin"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#IGF-1
exposure_dat <- read_exposure_data(
  filename = "data/Sinnott_Armstrong_et_al_2019/IGF_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$samplesize.exposure<-317114
exposure_dat$exposure<-"IGF-1"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


#IGFBP1
exposure_dat<- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGFBP1/IGFBP_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"IGFBP1"
exposure_dat$samplesize.exposure<-35560
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


#Resistin
exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/RETN/Resistin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"Resistin"
exposure_dat$samplesize.exposure<-35560
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Leptin
exposure_dat <- extract_instruments(outcomes="ebi-a-GCST90007310", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"Leptin"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Adiponectin
exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/ADIPOQ/Adiponectin_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"Adiponectin"
exposure_dat$samplesize.exposure<-35560
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#CXCL8
exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/CXCL8/CXCL8_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"CXCL8"
exposure_dat$samplesize.exposure<-35560
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#HDL cholesterol
exposure_dat <- extract_instruments(outcomes="ieu-b-109", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"HDL cholesterol"
exposure_dat$samplesize.exposure<-403943
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Triglycerides
exposure_dat <- extract_instruments(outcomes="ieu-b-111", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"Triglycerides"
exposure_dat$samplesize.exposure<-441016
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Bioavailable testosterone female
exposure_dat <- extract_instruments(outcomes="ieu-b-4869", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"Bioavailable testosterone female"
exposure_dat$samplesize.exposure<-441016
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Bioavailable testosterone male
exposure_dat <- extract_instruments(outcomes="ieu-b-4868", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"Bioavailable testosterone male"
exposure_dat$samplesize.exposure<-441016
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#BMI
exposure_dat <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"BMI"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

# Sex-specific ------------------------------------------------------------
exposure_dat<- extract_instruments(outcomes="ieu-a-974", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"Female BMI"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT_Female.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Female ASAT"
exposure_dat$samplesize.exposure<-19872
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT_Female.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Female VAT"
exposure_dat$samplesize.exposure<-19872
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT_Female.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Female GFAT"
exposure_dat$samplesize.exposure<-19872
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Fasting insulin
exposure_dat <- read_exposure_data(
  filename = "data/FI_GWAS/FI_female_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"Female fasting insulin"
exposure_dat$samplesize.exposure<-50404
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#IGF1
exposure_dat <- read_exposure_data(
  filename = "data/Neale_lab_2018/IGF1/female_IGF1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"Female IGF-1"
exposure_dat$samplesize.exposure<-194174
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#HDL cholesterol
exposure_dat <- read_exposure_data(
  filename = "data/Neale_lab_2018/HDL_cholesterol/female_HDL_cholesterol_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$samplesize.exposure<-194174
exposure_dat$exposure<-"Female HDL cholesterol"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Triglycerides
exposure_dat <- read_exposure_data(
  filename = "data/Neale_lab_2018/Triglycerides/female_Triglycerides_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"Female triglycerides"
exposure_dat$samplesize.exposure<-194174
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Male measures of adiposity
exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT_Male.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Male ASAT"
exposure_dat$samplesize.exposure<-19093
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT_Male.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Male VAT"
exposure_dat$samplesize.exposure<-19093
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT_Male.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Male GFAT"
exposure_dat$samplesize.exposure<-19093
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- extract_instruments(outcomes=c("ieu-a-2"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"BMI No UK Biobank"
exposure_dat$F<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
exposure_dat<-exposure_dat[exposure_dat$F>=10,]
fwrite(exposure_dat,paste("data/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))






