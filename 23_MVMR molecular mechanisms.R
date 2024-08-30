library(remotes)
library(MVMR)
library(dplyr)
library(TwoSampleMR)
library(data.table)


files<-list.files("data/LiverFatMR/MVMR_files/")

for (a in files){
  print(a)
  print(which(files==a))
mvmrdata <- fread(paste("data/LiverFatMR/MVMR_files/",a,sep=""))

a<-gsub("\\..*","",a)


F.data<-format_mvmr(BXGs=mvmrdata[,c(4,6)],
                    BYG=mvmrdata[,8],
                    seBXGs=mvmrdata[,c(5,7)],
                    seBYG=mvmrdata[,9],
                    RSID=mvmrdata[,1])

sres<-strength_mvmr(r_input=F.data,gencov=0)

write.table(sres, file = paste("results/LiverFatMR/",a,"_F_stats.csv",sep=""), sep = ",", col.names = NA,
            qmethod = "double")

pres<-pleiotropy_mvmr(r_input=F.data,gencov=0)
write.table(pres, file = paste("results/LiverFatMR/",a,"_Q_stats.csv",sep=""), sep = ",", col.names = NA,
            qmethod = "double")

res<-ivw_mvmr(r_input=F.data)
write.table(res, file = paste("results/LiverFatMR/",a,"_MVMR_results.csv",sep=""), sep = ",", col.names = NA,
            qmethod = "double")

}


# Calculating percent mediated and CIs ------------------------------------

results<-data.frame()
results2<-data.frame(matrix(ncol=5))
colnames(results2)<-c("analysis","percent","LCI","UCI","p")

#Sobel's SE
#indirect effect = Effect of adiposity measure on molecular trait (a) multiplied by effect of the molecular trait on the outcome (adjusted for measure of adiposity) (b)
#Sobel's SE = sqrt(a^2*aSE^2+b^2*bSE^2)

adiposity_molecular <- fread("results/LiverFatMR/Adiposity_measures_molecular_traits_results_new_proxies.csv")

#ASAT_FI_EC
MVMR <- read.csv("results/LiverFatMR/ASAT_FI_EC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Fasting insulin",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_FI <- adi_mol[1,7]
ASAT_FI_se <-adi_mol[1,8]
#b
FI_EC <- MVMR[2,2]
FI_EC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_FI^2)*(FI_EC_se^2))+((FI_EC^2)*(ASAT_FI_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_FI*FI_EC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_FI_EC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#ASAT_FI_EEC
MVMR <- read.csv("results/LiverFatMR/ASAT_FI_EEC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Fasting insulin",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_FI <- adi_mol[1,7]
ASAT_FI_se <-adi_mol[1,8]
#b
FI_EEC <- MVMR[2,2]
FI_EEC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_FI^2)*(FI_EEC_se^2))+((FI_EEC^2)*(ASAT_FI_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrioid endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_FI*FI_EEC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_FI_EEC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#ASAT_SHBG_EC
MVMR <- read.csv("results/LiverFatMR/ASAT_SHBG_EC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="Female ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_SHBG <- adi_mol[1,7]
ASAT_SHBG_se <-adi_mol[1,8]
#b
SHBG_EC <- MVMR[2,2]
SHBG_EC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_SHBG^2)*(SHBG_EC_se^2))+((SHBG_EC^2)*(ASAT_SHBG_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_SHBG*SHBG_EC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_SHBG_EC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#ASAT_SHBG_EEC
MVMR <- read.csv("results/LiverFatMR/ASAT_SHBG_EEC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="Female ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_SHBG <- adi_mol[1,7]
ASAT_SHBG_se <-adi_mol[1,8]
#b
SHBG_EEC <- MVMR[2,2]
SHBG_EEC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_SHBG^2)*(SHBG_EEC_se^2))+((SHBG_EEC^2)*(ASAT_SHBG_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrioid endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_SHBG*SHBG_EEC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_SHBG_EEC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)


#ASAT_CXC_OC
MVMR <- read.csv("results/LiverFatMR/ASAT_CXC_OC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="CXCL8",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_CXC <- adi_mol[1,7]
ASAT_CXC_se <-adi_mol[1,8]
#b
CXC_OC <- MVMR[2,2]
CXC_OC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_CXC^2)*(CXC_OC_se^2))+((CXC_OC^2)*(ASAT_CXC_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Ovarian cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_CXC*CXC_OC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_CXC_OC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)



#ASAT_HDL_TNBC
MVMR <- read.csv("results/LiverFatMR/ASAT_HDL_TNBC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="HDL cholesterol || id:ieu-b-109",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

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

results<-rbind(results,results2)

#ASAT_IGFBP_OAC
MVMR <- read.csv("results/LiverFatMR/ASAT_IGFBP_OAC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="IGFBP1",]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_IGFBP <- adi_mol[1,7]
ASAT_IGFBP_se <-adi_mol[1,8]
#b
IGFBP_OAC <- MVMR[2,2]
IGFBP_OAC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_IGFBP^2)*(IGFBP_OAC_se^2))+((IGFBP_OAC^2)*(ASAT_IGFBP_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Esophagus adenocarcinoma",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_IGFBP*IGFBP_OAC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_IGFBP_OAC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#ASAT_BT_EC
MVMR <- read.csv("results/LiverFatMR/ASAT_BT_EC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="Female ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Bioavailable Testosterone || id:ieu-b-4869"  ,]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_BT <- adi_mol[1,7]
ASAT_BT_se <-adi_mol[1,8]
#b
BT_EC <- MVMR[2,2]
BT_EC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_BT^2)*(BT_EC_se^2))+((BT_EC^2)*(ASAT_BT_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_BT*BT_EC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_BT_EC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#ASAT_BT_EEC
MVMR <- read.csv("results/LiverFatMR/ASAT_BT_EEC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="Female ASAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Bioavailable Testosterone || id:ieu-b-4869"  ,]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
ASAT_BT <- adi_mol[1,7]
ASAT_BT_se <-adi_mol[1,8]
#b
BT_EEC <- MVMR[2,2]
BT_EEC_se <- MVMR[2,3]

sigmay <- sqrt(((ASAT_BT^2)*(BT_EEC_se^2))+((BT_EEC^2)*(ASAT_BT_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrioid endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- ASAT_BT*BT_EEC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"ASAT_BT_EEC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

#GFAT_AD_NEC
MVMR <- read.csv("results/LiverFatMR/GFAT_AD_NEC_MVMR_results.csv")
adi_mol<-adiposity_molecular[adiposity_molecular$exposure=="GFAT",]
adi_mol<-adi_mol[adi_mol$outcome=="Adiponectin" ,]
adi_mol<-adi_mol[(adi_mol$method=="Inverse variance weighted" | adi_mol$method=="Wald ratio"),]

#a
GFAT_AD <- adi_mol[1,7]
GFAT_AD_se <-adi_mol[1,8]
#b
AD_NEC <- MVMR[2,2]
AD_NEC_se <- MVMR[2,3]

sigmay <- sqrt(((GFAT_AD^2)*(AD_NEC_se^2))+((AD_NEC^2)*(GFAT_AD_se^2)))

#total effect
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="GFAT",]
res<-res[res$outcome=="Non-endometrioid endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
total_effect <- res[1,7]

#se for % mediated
se_mediated <- sqrt((sigmay^2)/(total_effect^2))

#Calculating mediation percent
Indirect <- GFAT_AD*AD_NEC
mediated <- Indirect/total_effect
LCI_mediated <- mediated - 1.96*se_mediated
UCI_mediated <- mediated + 1.96*se_mediated
#p-value
z <- as.numeric(mediated/se_mediated)
p <- 2*pnorm(-abs(z))

results2$analysis<-"GFAT_AD_NEC"
results2$percent<-mediated[1,1]
results2$LCI<-LCI_mediated[1,1]
results2$UCI<-UCI_mediated[1,1]
results2$p<-p

results<-rbind(results,results2)

results[,2]<-unlist(results[,2])
results[,3]<-unlist(results[,3])
results[,4]<-unlist(results[,4])

fwrite(results,"results/LiverFatMR/MVMR_Results.csv")
