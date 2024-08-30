#sink("Rscript.txt")
library(MRInstruments)
library(TwoSampleMR)
library(data.table)
library(dplyr)


# Set up functions --------------------------------------------------------

r2 <- function(beta, se, af, n) {
  num <- 2*(beta^2)*af*(1-af)
  den <- 2*(beta^2)*af*(1-af) + (se^2)*2*n*af*(1-af)
  pve <- num/den
  return(pve)
}

f_and_r2 <- function(exposure_dat) {
  F_all<-(exposure_dat$beta.exposure/exposure_dat$se.exposure)^2
  F<-mean(F_all)
  r2a<-r2(beta=exposure_dat$beta.exposure,exposure_dat$se.exposure,af=exposure_dat$eaf.exposure,n=exposure_dat$samplesize.exposure)
  r=sum(r2a)
  res<<-c(exposure_dat$exposure[1],F,r)
  
  paste(c("Minimum F stat for ", exposure_dat$exposure[1],"is ",min(F_all)),sep=)
}

results<-data.frame()

# Measures of adiposity ---------------------------------------------------

exposure_dat <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Liver fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$exposure<-"Liver fat"
exposure_dat$samplesize.exposure<-32858
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
f_and_r2(exposure_dat)
results<-rbind(results,res)
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/LiverFatMR/Instruments/Pancreas fat_instruments.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure[exposure_dat$exposure=="Pancreas fat || id:ebi-a-GCST90016675"]<-"Pancreas fat"
exposure_dat$samplesize.exposure<-25617
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"VAT"
exposure_dat$samplesize.exposure<-38965
f_and_r2(exposure_dat)
results<-rbind(results,res) 
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"ASAT"
exposure_dat$samplesize.exposure<-38965
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"GFAT"
exposure_dat$samplesize.exposure<-38965
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

# Molecular traits --------------------------------------------------------

#SHBG female
exposure_dat <- extract_instruments(outcomes=c("ieu-b-4870"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-214989
exposure_dat$exposure<-"SHBG (female)"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#SHBG male
exposure_dat <- extract_instruments(outcomes=c("ieu-b-4871"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-185221
exposure_dat$exposure<-"SHBG (male)"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Bioavailable testosterone female
exposure_dat <- extract_instruments(outcomes=c("ieu-b-4869"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-180386
exposure_dat$exposure<-"Bioavailable testosterone (female)"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Bioavailable testosterone male
exposure_dat <- extract_instruments(outcomes=c("ieu-b-4868"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-184205
exposure_dat$exposure<-"Bioavailable testosterone (male)"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Total testosterone male
exposure_dat <- extract_instruments(outcomes=c("ieu-b-4865"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-199569
exposure_dat$exposure<-"Total testosterone (male)"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-98210
exposure_dat$exposure<-"Fasting insulin"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-317114
exposure_dat$exposure<-"IGF-1"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Resistin"
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Leptin
exposure_dat <- extract_instruments(outcomes="ebi-a-GCST90007310", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Leptin"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Adiponectin"
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#PAI1
exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/PAI1/PAI1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"PAI1"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#FASN
exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/FASN/FASN_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure"
)
exposure_dat$exposure<-"FASN"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-35560
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#HDL cholesterol
exposure_dat <- extract_instruments(outcomes="ieu-b-109", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$exposure<-"HDL cholesterol"
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-403943
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Triglycerides
exposure_dat <- extract_instruments(outcomes="ieu-b-111", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Triglycerides"
exposure_dat$samplesize.exposure<-441016
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#Bioavailable testosterone
exposure_dat <- extract_instruments(outcomes="ebi-a-GCST90012104", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Bioavailable testosterone"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

#BMI
exposure_dat <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"BMI"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


# Sex-specific ------------------------------------------------------------
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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-19872
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-19872
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-19872
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


exposure_dat<- extract_instruments(outcomes="ieu-a-974", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Female BMI"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat<- extract_instruments(outcomes="ieu-a-785", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Male BMI"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Female fasting insulin"
exposure_dat$samplesize.exposure<-50404
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-194174
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-194174
exposure_dat$exposure<-"Female HDL cholesterol"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Female triglycerides"
exposure_dat$samplesize.exposure<-194174
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Male ASAT"
exposure_dat$samplesize.exposure<-19093
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$samplesize.exposure<-19093
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

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
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"Male GFAT"
exposure_dat$samplesize.exposure<-19093
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


#Male molecular traits


# Sample overlap ----------------------------------------------------------

exposure_dat <- extract_instruments(outcomes=c("ieu-a-2"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"BMI No UK Biobank"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- extract_instruments(outcomes=c("ebi-a-GCST002223"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"HDL cholesterol No UK Biobank"
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))

exposure_dat <- read_exposure_data(
  filename = "data/Ferkingstad_et_al_2021_proteins/IGF1/IGF_1_Instrument.txt",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat$pval.exposure[exposure_dat$pval.exposure==0]<-1*10^-200
exposure_dat$exposure<-"IGF1 No UK Biobank"
exposure_dat$samplesize.exposure<-38965
f_and_r2(exposure_dat)
results<-rbind(results,res)
if ("eaf.exposure" %in% colnames(exposure_dat)==FALSE){
  exposure_dat$eaf.exposure<-"Missing/unclear"
}
exposure_dat2<-select(exposure_dat,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure)
colnames(exposure_dat2)<-c("SNP","Effect_allele","Other_allele","EAF","Beta","Standard_error","P_value")
fwrite(exposure_dat2,paste("results/LiverFatMR/Instruments/",exposure_dat$exposure[1],"_instruments.csv",sep=""))


fwrite(results,"results/LiverFatMR/Adiposity_measures_instruments_F_stat.csv")

sink()

