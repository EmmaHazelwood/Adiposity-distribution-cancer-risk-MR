
library(remotes)
library(dplyr)
library(TwoSampleMR)
library(data.table)


# Get data ready ----------------------------------------------------------


## Sample overlap - HDL ----------------------------------------------------

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

fwrite(ASAT_HDL_TNBC,"data/LiverFatMR/MVMR_sensitivity analyses/ASAT_HDL_nooverlap_TNBC.csv")


# Do MVMR -----------------------------------------------------------------
files<-list.files("data/LiverFatMR/MVMR_sensitivity analyses/")

for (a in files){
  print(a)
  print(which(files==a))
  mvmrdata <- fread(paste("data/LiverFatMR/MVMR_sensitivity analyses/",a,sep=""))
  
  a<-gsub("\\..*","",a)
  
  
  F.data<-format_mvmr(BXGs=mvmrdata[,c(4,6)],
                      BYG=mvmrdata[,8],
                      seBXGs=mvmrdata[,c(5,7)],
                      seBYG=mvmrdata[,9],
                      RSID=mvmrdata[,1])
  
  sres<-strength_mvmr(r_input=F.data,gencov=0)
  
  write.table(sres, file = paste("results/LiverFatMR/MVMR_sensitivity analyses/",a,"_F_stats.csv",sep=""), sep = ",", col.names = NA,
              qmethod = "double")
  
  pres<-pleiotropy_mvmr(r_input=F.data,gencov=0)
  write.table(pres, file = paste("results/LiverFatMR/MVMR_sensitivity analyses/",a,"_Q_stats.csv",sep=""), sep = ",", col.names = NA,
              qmethod = "double")
  
  res<-ivw_mvmr(r_input=F.data)
  write.table(res, file = paste("results/LiverFatMR/MVMR_sensitivity analyses/",a,"_MVMR_results.csv",sep=""), sep = ",", col.names = NA,
              qmethod = "double")
  
}

results2<-data.frame(matrix(ncol=5))
colnames(results2)<-c("analysis","percent","LCI","UCI","p")

#Sobel's SE
#indirect effect = Effect of adiposity measure on molecular trait (a) multiplied by effect of the molecular trait on the outcome (adjusted for measure of adiposity) (b)
#Sobel's SE = sqrt(a^2*aSE^2+b^2*bSE^2)

adiposity_molecular <- fread("results/LiverFatMR/Adiposity_measures_molecular_traits_results_new_proxies.csv")

#ASAT_HDL_TNBC
MVMR <- read.csv("results/LiverFatMR/MVMR_sensitivity analyses/ASAT_HDL_nooverlap_TNBC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="HDL cholesterol || id:ieu-b-109",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#Update with no overlap betas and se
adi_mol$beta<--0.109748713
adi_mol$se<-0.108967542

#a
ASAT_HDL <- adi_mol[1,7]
ASAT_HDL_se <-adi_mol[1,8]
#b
HDL_TNBC <- MVMR[2,2]
HDL_TNBC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_HDL^2)*(HDL_TNBC_se^2))+((HDL_TNBC^2)*(ASAT_HDL_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Triple negative or basal-like breast cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_HDL*HDL_TNBC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_HDL_TNBC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p
colnames(results2)<-c("analysis","percent","LCI","UCI","p")

results2[,2]<-unlist(results2[,2])
results2[,3]<-unlist(results2[,3])
results2[,4]<-unlist(results2[,4])

fwrite(results2,"results/LiverFatMR/MVMR_sensitivity analyses/ASAT_HDL_TNBC_mediation_results.csv")
