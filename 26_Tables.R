library(data.table)
library(dplyr)
library(readxl)




# Table 1 -----------------------------------------------------------------

tb_1<-read_excel("Tables.xlsx",sheet=1)
#Just do it manually


# Table 2 -----------------------------------------------------------------

tb_2<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
tb_2$`Odds ratio`<-exp(tb_2$b)
tb_2$LCI<-exp(tb_2$b-1.96*tb_2$se)
tb_2$UCI<-exp(tb_2$b+1.96*tb_2$se)
tb_2$`Odds ratio`<-format(round(tb_2$`Odds ratio`,digits=2),scientific=FALSE)
tb_2$LCI<-format(round(tb_2$LCI,digits=2),scientific=FALSE)
tb_2$UCI<-format(round(tb_2$UCI,digits=2),scientific=FALSE)

tb_2$pval<-formatC(tb_2$pval, format = "e", digits = 2)

tb_2<-tb_2[!(tb_2$method=="Weighted median" & tb_2$nsnp<10),]
tb_2<-tb_2[!(tb_2$method=="Weighted mode" & tb_2$nsnp<10),]


tb_2$`95% confidence interval`<-paste(tb_2$LCI," to ",tb_2$UCI,sep="")

tb_2<-dplyr::select(tb_2,exposure,outcome,method,nsnp,`Odds ratio`,`95% confidence interval`,pval)
colnames(tb_2)<-c("Exposure","Outcome","MR model","Number of SNPs","Odds ratio","95% confidence interval","P-value")
tb_2$Outcome<-gsub("uminal B-HER2","uminal B/HER2",tb_2$Outcome)
tb_2$Outcome<-gsub("Esophagus","Oesophageal",tb_2$Outcome)
tb_2<-distinct(tb_2)


# Table 3 -----------------------------------------------------------------

tb_3<-fread("results/LiverFatMR/Adiposity_measures_molecular_traits_results_new_proxies.csv")
tb_3$Beta<-tb_3$b
tb_3$LCI<-tb_3$b-1.96*tb_3$se
tb_3$UCI<-tb_3$b+1.96*tb_3$se
tb_3$Beta<-format(round(tb_3$Beta,digits=2),scientific=FALSE)
tb_3$LCI<-format(round(tb_3$LCI,digits=2),scientific=FALSE)
tb_3$UCI<-format(round(tb_3$UCI,digits=2),scientific=FALSE)

tb_3$pval<-formatC(tb_3$pval, format = "e", digits = 2)

tb_3<-tb_3[!(tb_3$method=="Weighted median" & tb_3$nsnp<10),]
tb_3<-tb_3[!(tb_3$method=="Weighted mode" & tb_3$nsnp<10),]
tb_3$`95% confidence interval`<-paste(tb_3$LCI," to ",tb_3$UCI,sep="")

tb_3<-distinct(tb_3)

tb_3<-dplyr::select(tb_3,exposure,outcome,method,nsnp,Beta,`95% confidence interval`,pval)
colnames(tb_3)<-c("Exposure","Outcome","MR model","Number of SNPs","Beta","95% confidence interval","P-value")
tb_3$Outcome[tb_3$Outcome=="Total Testosterone || id:ieu-b-4864"]<-"Total testosterone (female)"
tb_3$Outcome[tb_3$Outcome=="Total Testosterone || id:ieu-b-4865"]<-"Total testosterone (male)"
tb_3$Outcome[tb_3$Outcome=="Total Testosterone || id:ieu-b-4868"]<-"Bioavailable testosterone (male)"
tb_3$Outcome[tb_3$Outcome=="Total Testosterone || id:ieu-b-4869"]<-"Bioavailable testosterone (female)"
tb_3$Outcome[tb_3$Outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870"]<-"SHBG (male)"
tb_3$Outcome[tb_3$Outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4871"]<-"SHBG (female)"


tb_3$Outcome<-gsub("\\|.*","",tb_3$Outcome)
tb_3$Outcome<-gsub("levels","",tb_3$Outcome)
tb_3$Outcome<-gsub("triglycerides ","Triglycerides",tb_3$Outcome)
tb_3$Outcome<-gsub("CXCL8","CXCL-8",tb_3$Outcome)
tb_3$Outcome<-gsub("MCP1","MCP-1",tb_3$Outcome)
tb_3$Outcome<-gsub("PAI1","PAI-1",tb_3$Outcome)
tb_3$Outcome<-gsub("IGFBP1","IGFBP-1",tb_3$Outcome)
tb_3$Outcome<-gsub("IGF1","IGF-1",tb_3$Outcome)
tb_3$Outcome<-gsub("IGF2","IGF-2",tb_3$Outcome)
tb_3$Outcome<-gsub("IGFBP3","IGFBP-3",tb_3$Outcome)
tb_3$Outcome<-gsub("TNF-a","TNF-α",tb_3$Outcome)
tb_3$Outcome<-gsub("IFN-a","IFN-α",tb_3$Outcome)
tb_3$Outcome<-gsub("IFN-B","IFN-β",tb_3$Outcome)
tb_3$Outcome<-gsub("IL-1B","IL-1β",tb_3$Outcome)
tb_3$Exposure<-gsub("Male ","",tb_3$Exposure)
tb_3$Exposure<-gsub("Female ","",tb_3$Exposure)
tb_3$Exposure<-gsub("body mass index || id:ieu-b-40","BMI",tb_3$Exposure)


tb_3<-distinct(tb_3)


# Table 4 -----------------------------------------------------------------

tb_4<-fread("results/LiverFatMR/Molecular_traits_cancers_results.csv")
tb_4$`Odds ratio`<-exp(tb_4$b)
tb_4$LCI<-exp(tb_4$b-1.96*tb_4$se)
tb_4$UCI<-exp(tb_4$b+1.96*tb_4$se)
tb_4$`Odds ratio`<-format(round(tb_4$`Odds ratio`,digits=2),scientific=FALSE)
tb_4$LCI<-format(round(tb_4$LCI,digits=2),scientific=FALSE)
tb_4$UCI<-format(round(tb_4$UCI,digits=2),scientific=FALSE)

tb_4$pval<-formatC(tb_4$pval, format = "e", digits = 2)

tb_4<-tb_4[!(tb_4$method=="Weighted median" & tb_4$nsnp<10),]
tb_4<-tb_4[!(tb_4$method=="Weighted mode" & tb_4$nsnp<10),]


tb_4$`95% confidence interval`<-paste(tb_4$LCI," to ",tb_4$UCI,sep="")

tb_4$exposure[tb_4$exposure=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870"]<-"SHBG (female)"
tb_4$exposure[tb_4$exposure=="Sex hormone binding globulin (SHBG) || id:ieu-b-4871"]<-"SHBG (male)"
tb_4$exposure[tb_4$exposure=="triglycerides || id:ieu-b-111"]<-"Triglycerides"
tb_4$exposure[tb_4$exposure=="HDL cholesterol || id:ieu-b-109"]<-"HDL cholesterol"
tb_4$exposure[tb_4$exposure=="Total Testosterone || id:ieu-b-4864"]<-"Total testosterone (female)"
tb_4$exposure[tb_4$exposure=="Total Testosterone || id:ieu-b-4865"]<-"Total testosterone (male)"
tb_4$exposure[tb_4$exposure=="Bioavailable Testosterone || id:ieu-b-4868"]<-"Bioavailable testosterone (male)"
tb_4$exposure[tb_4$id.exposure=="ieu-b-4864"]<-"Total testosterone (female)"
tb_4$exposure[tb_4$id.exposure=="ieu-b-4868"]<-"Bioavailable testosterone (male)"
tb_4$exposure[tb_4$exposure=="Bioavailable Testosterone || id:ieu-b-4869"]<-"Bioavailable testosterone (female)"
tb_4$exposure[tb_4$id.exposure=="ieu-b-4869"]<-"Bioavailable testosterone (female)"
tb_4$exposure[tb_4$exposure=="IGFBP1"]<-"IGFBP-1"
tb_4$exposure[tb_4$exposure=="CXCL8"]<-"CXCL-8"
tb_4$exposure[tb_4$exposure=="PAI1"]<-"PAI-1"
tb_4$exposure[tb_4$exposure=="Adiponectin"]<-"              Adiponectin"
tb_4$outcome<-gsub("uminal B-HER2","uminal B/HER2",tb_4$outcome)

tb_4<-dplyr::select(tb_4,exposure,outcome,method,nsnp,`Odds ratio`,`95% confidence interval`,pval)

colnames(tb_4)<-c("Exposure","Outcome","MR model","Number of SNPs","Odds ratio","95% confidence interval","P-value")

tb_4$Outcome<-gsub("\\|.*","",tb_4$Outcome)
tb_4$Exposure<-gsub("\\|.*","",tb_4$Exposure)
tb_4$Exposure<-gsub("levels","",tb_4$Exposure)
tb_4$Exposure<-gsub("triglycerides ","Triglycerides",tb_4$Exposure)
tb_4$Outcome<-gsub("uminal B-HER2","uminal B/HER2",tb_4$Outcome)
tb_4$Exposure<-gsub("CXCL8","CXCL-8",tb_4$Exposure)
tb_4$Exposure<-gsub("MCP1","MCP-1",tb_4$Exposure)
tb_4$Exposure<-gsub("PAI1","PAI-1",tb_4$Exposure)
tb_4$Exposure<-gsub("IGFBP1","IGFBP-1",tb_4$Exposure)
tb_4$Exposure<-gsub("              Adiponectin","Adiponectin",tb_4$Exposure)
tb_4$Outcome<-gsub("Esophagus","Oesophageal",tb_4$Outcome)

tb_4<-distinct(tb_4)

#Check against ones we should have

x<-fread("Liver\ Fat\ MR/Analyses molecular MR.csv")

x$outcome<-gsub("\\|.*","",x$outcome)
x$exposure<-gsub("\\|.*","",x$exposure)
x$exposure<-gsub("levels","",x$exposure)
x$exposure<-gsub("triglycerides ","Triglycerides",x$exposure)
x$outcome<-gsub("uminal B-HER2","uminal B/HER2",x$outcome)
x$exposure<-gsub("circulating leptin","Leptin",x$exposure)
x$exposure<-gsub("IGF1","IGF-1",x$exposure)
x$exposure<-gsub("PAI1","PAI-1",x$exposure)
x$x<-"shouldbe"

tb_4x<-dplyr::select(tb_4,Exposure,Outcome)
tb_4x<-distinct(tb_4x)
tb_4x$x<-"is"

both<-merge(tb_4x,x,by.x=c("Exposure","Outcome"),by.y=c("exposure","outcome"),allow.cartesian=TRUE,all.x=T,all.y=T)


# Table 5 -----------------------------------------------------------------

tb_5<-fread("results/LiverFatMR/MVMR_Results.csv")
tb_5b<-fread("results/LiverFatMR/Adiposity_cancers_mediation_results.csv")

tb_5$`Unadjusted OR`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")
tb_5$`Unadjusted OR`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")

tb_5$`Unadjusted OR`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")

tb_5$`Unadjusted OR`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")

tb_5$`Unadjusted OR`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")

tb_5$`Unadjusted OR`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Unadjusted LCI`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"]- 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Unadjusted UCI`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Unadjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"]+ 1.96*tb_5b$`Unadjusted SE`[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Unadjusted 95% confidence interval`<-paste((round(tb_5$`Unadjusted LCI`,digits=2))," to ",(round(tb_5$`Unadjusted UCI`,digits=2)),sep="")

tb_5$`Adjusted OR`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Adjusted LCI`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Adjusted UCI`[tb_5$analysis=="ASAT_FI_EC" | tb_5$analysis=="ASAT_SHBG_EC" | tb_5$analysis=="ASAT_BT_EC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrial cancer"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")
tb_5$`Adjusted OR`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Adjusted LCI`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Adjusted UCI`[tb_5$analysis=="ASAT_FI_EEC" | tb_5$analysis=="ASAT_SHBG_EEC" | tb_5$analysis=="ASAT_BT_EEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Endometrioid endometrial cancer"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")

tb_5$`Adjusted OR`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Adjusted LCI`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Adjusted UCI`[tb_5$analysis=="ASAT_CXC_OC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Ovarian cancer"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")

tb_5$`Adjusted OR`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Adjusted LCI`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Adjusted UCI`[tb_5$analysis=="ASAT_HDL_TNBC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Triple negative or basal-like breast cancer"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")

tb_5$`Adjusted OR`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Adjusted LCI`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Adjusted UCI`[tb_5$analysis=="ASAT_IGFBP_OAC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="ASAT" & tb_5b$Outcome=="Oesophageal adenocarcinoma"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")

tb_5$`Adjusted OR`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Adjusted LCI`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"]- 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Adjusted UCI`[tb_5$analysis=="GFAT_AD_NEC"]<-exp(tb_5b$Adjusted[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"]+ 1.96*tb_5b$`Adjusted SE`[tb_5b$Exposure=="GFAT" & tb_5b$Outcome=="Non-endometrioid endometrial cancer"])
tb_5$`Adjusted 95% confidence interval`<-paste((round(tb_5$`Adjusted LCI`,digits=2))," to ",(round(tb_5$`Adjusted UCI`,digits=2)),sep="")

tb_5$analysis<-gsub("ASAT_FI_EC","Exposure 1: ASAT; Exposure 2: fasting insulin; Outcome: Endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_FI_EEC","Exposure 1: ASAT; Exposure 2: fasting insulin; Outcome: Endometrioid endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_SHBG_EC","Exposure 1: ASAT; Exposure 2: SHBG; Outcome: Endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_SHBG_EEC","Exposure 1: ASAT; Exposure 2: SHBG; Outcome: Endometrioid endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("GFAT_AD_NEC","Exposure 1: GFAT; Exposure 2: adiponectin; Outcome: Non-endometrioid endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_HDL_TNBC","Exposure 1: ASAT; Exposure 2: HDL cholesterol; Outcome: Triple negative or basal-like breast cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_CXC_OC","Exposure 1: ASAT; Exposure 2: CXCL-8; Outcome: Ovarian cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_IGFBP_OAC","Exposure 1: ASAT; Exposure 2: IGFBP-1; Outcome: Oesophageal adenocarcinoma",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_BT_EC","Exposure 1: ASAT; Exposure 2:Bioavailable testosterone; Outcome: Endometrial cancer",tb_5$analysis)
tb_5$analysis<-gsub("ASAT_BT_EEC","Exposure 1: ASAT; Exposure 2:Bioavailable testosterone; Outcome: Endometrioid endometrial cancer",tb_5$analysis)

tb_5$percent<-round(tb_5$percent*100,digits=2)
tb_5$`95% confidence interval`<-paste((round(tb_5$LCI*100,digits=2))," to ",(round(tb_5$UCI*100,digits=2)),sep="")

tb_5$p<-formatC(tb_5$p, format = "e", digits = 2)

tb_5<-dplyr::select(tb_5,analysis,`Unadjusted OR`,`Unadjusted 95% confidence interval`,`Adjusted OR`,`Adjusted 95% confidence interval`,percent,`95% confidence interval`,p)
colnames(tb_5)<-c("Analysis","Unadjusted OR","Unadjusted 95% confidence interval","Adjusted OR","Adjusted 95% confidence interval","Percent mediated","95% confidence interval for percent mediated","P-value for percent mediated")

# Write as excel file -----------------------------------------------------

library(openxlsx)

xl_lst <- list('Table 1' = tb_1, 'Table 2' = tb_2, 'Table 3' = tb_3, 'Table 4' = tb_4, 'Table 5' = tb_5)

write.xlsx(xl_lst, file = "Tables.xlsx")




# Supplementary tables ----------------------------------------------------


# Supplementary table 1 ---------------------------------------------------

st_1<-read_excel("Supplementary tables.xlsx")


# Supplementary table 2 ---------------------------------------------------

files<-list.files("results/LiverFatMR/Instruments/")
files<-paste("results/LiverFatMR/Instruments/",files,sep="")
st_2<-data.frame()
for (a in files){
  df<-fread(a)
  SNP<-a
  a<-sub("results/LiverFatMR/Instruments/","",a)
  a<-sub(".csv","",a)
  df$Exposure<-a
  st_2<-rbind(st_2,df)
}
st_2$Exposure[st_2$Exposure=="exposure"]<-"Pancreas fat"
st_2<-dplyr::select(st_2,Exposure,SNP,Effect_allele,Other_allele,EAF,Beta,Standard_error,P_value)


# Supplementary table 3 ---------------------------------------------------
cancerlist<-c("Thyroid cancer (FinnGen)","Multiple myeloma (FinnGen)" ,"Liver cancer (FinnGen)","Endometrial cancer" ,"Endometrioid endometrial cancer","Oesophageal adenocarcinoma","Proximal colon cancer","Low malignant potential ovarian cancer","Colon cancer","Colorectal cancer (overall)" ,"Non-endometrioid endometrial cancer","Triple negative or basal-like breast cancer","Breast cancer" ,"Luminal B/HER2-negative-like breast cancer","Luminal B-like breast cancer","HER2-enriched-like breast cancer","Invasive mucinous ovarian cancer","Ovarian cancer" ,"Distal colon cancer")
st_3<-fread("results/LiverFatMR/Repeating_no_sample_overlap.csv")
st_3$outcome<-gsub("uminal B-HER2","uminal B/HER2",st_3$outcome)
st_3$outcome<-gsub(" no UK Biobank"," (no UK Biobank)",st_3$outcome)

st_3<-st_3[!(st_3$method=="Weighted median" & st_3$nsnp<10),]
st_3<-st_3[!(st_3$method=="Weighted mode" & st_3$nsnp<10),]

st_3$LCI<-st_3$b-1.96*st_3$se
st_3$UCI<-st_3$b+1.96*st_3$se

st_3$b[st_3$outcome %in% cancerlist]<-exp(st_3$b[st_3$outcome %in% cancerlist])

st_3$`Odds ratio`<-st_3$b

st_3$LCI[st_3$outcome %in% cancerlist]<-exp(st_3$LCI[st_3$outcome %in% cancerlist])
st_3$UCI[st_3$outcome %in% cancerlist]<-exp(st_3$UCI[st_3$outcome %in% cancerlist])
st_3$`Odds ratio`<-format(round(st_3$`Odds ratio`,digits=2),scientific=FALSE)
st_3$LCI<-format(round(st_3$LCI,digits=2),scientific=FALSE)
st_3$UCI<-format(round(st_3$UCI,digits=2),scientific=FALSE)

st_3$pval<-formatC(st_3$pval, format = "e", digits = 2)

st_3$`95% confidence interval`<-paste(st_3$LCI," to ",st_3$UCI,sep="")

st_3<-dplyr::select(st_3,exposure,outcome,method,nsnp,`Odds ratio`,`95% confidence interval`,pval)
colnames(st_3)<-c("Exposure","Outcome","MR model","Number of SNPs","Odds ratio/Beta","95% confidence interval","P-value")
st_3<-distinct(st_3)

st_3$Outcome<-gsub("\\|.*","",st_3$Outcome)
st_3$Exposure<-gsub("\\|.*","",st_3$Exposure)


# Supplementary table 4 ---------------------------------------------------
cancerlist<-c("Endometrial cancer" ,"Endometrioid endometrial cancer","Non-endometrioid endometrial cancer","Ovarian cancer","High grade serous carcinoma ovarian cancer","Low grade serous carcinoma ovarian cancer","Invasive mucinous ovarian cancer","Endometrioid ovarian cancer","Clear cell ovarian cancer","Low malignant potential ovarian cancer","Breast cancer","Luminal A-like breast cancer" ,"Luminal B/HER2-negative-like breast cancer","Luminal B-like breast cancer" ,"HER2-enriched-like breast cancer","Triple negative or basal-like breast cancer" ,"Luminal B-HER2-negative-like breast cancer","Liver cancer","Meningioma" ,"Multiple myeloma" ,"Proximal colon cancer")
st_4<-fread("results/LiverFatMR/Sex_specific_results.csv")
st_4$LCI<-st_4$b-1.96*st_4$se
st_4$UCI<-st_4$b+1.96*st_4$se

st_4<-st_4[!(st_4$method=="Weighted median" & st_4$nsnp<10),]
st_4<-st_4[!(st_4$method=="Weighted mode" & st_4$nsnp<10),]


st_4$b[st_4$outcome %in% cancerlist]<-exp(st_4$b[st_4$outcome %in% cancerlist])

st_4$`Odds ratio`<-st_4$b

st_4$LCI[st_4$outcome %in% cancerlist]<-exp(st_4$LCI[st_4$outcome %in% cancerlist])
st_4$UCI[st_4$outcome %in% cancerlist]<-exp(st_4$UCI[st_4$outcome %in% cancerlist])
st_4$`Odds ratio`<-format(round(st_4$`Odds ratio`,digits=2),scientific=FALSE)
st_4$LCI<-format(round(st_4$LCI,digits=2),scientific=FALSE)
st_4$UCI<-format(round(st_4$UCI,digits=2),scientific=FALSE)

st_4$pval<-formatC(st_4$pval, format = "e", digits = 2)

st_4$`95% confidence interval`<-paste(st_4$LCI," to ",st_4$UCI,sep="")

st_4<-dplyr::select(st_4,exposure,outcome,method,nsnp,`Odds ratio`,`95% confidence interval`,pval)
colnames(st_4)<-c("Exposure","Outcome","MR model","Number of SNPs","Odds ratio/Beta","95% confidence interval","P-value")

st_4$Outcome<-gsub("\\|.*","",st_4$Outcome)
st_4$Exposure<-gsub("\\|.*","",st_4$Exposure)
st_4$Exposure<-gsub("IGF1","IGF-1",st_4$Exposure)
st_4$Outcome<-gsub("uminal B-HER2","uminal B/HER2",st_4$Outcome)


# Supplementary table 5 ---------------------------------------------------
st_6<-fread("results/LiverFatMR/BMI_obesity_cancers_results.csv")

st_6<-st_6[!(st_6$method=="Weighted median" & st_6$nsnp<10),]
st_6<-st_6[!(st_6$method=="Weighted mode" & st_6$nsnp<10),]


st_6$`Odds ratio`<-exp(st_6$b)
st_6$LCI<-exp(st_6$b-1.96*st_6$se)
st_6$UCI<-exp(st_6$b+1.96*st_6$se)
st_6$`Odds ratio`<-format(round(st_6$`Odds ratio`,digits=2),scientific=FALSE)
st_6$LCI<-format(round(st_6$LCI,digits=2),scientific=FALSE)
st_6$UCI<-format(round(st_6$UCI,digits=2),scientific=FALSE)

st_6$pval<-formatC(st_6$pval, format = "e", digits = 2)

st_6$`95% confidence interval`<-paste(st_6$LCI," to ",st_6$UCI,sep="")

st_6<-dplyr::select(st_6,exposure,outcome,method,nsnp,`Odds ratio`,`95% confidence interval`,pval)
colnames(st_6)<-c("Exposure","Outcome","MR model","Number of SNPs","Odds ratio","95% confidence interval","P-value")

st_6$Outcome<-gsub("\\|.*","",st_6$Outcome)
st_6$Exposure<-gsub("\\|.*","",st_6$Exposure)
st_6$Outcome<-gsub("Esophagus","Oesophageal",st_6$Outcome)

st_6$Exposure<-"Body mass index"

st_6<-distinct(st_6)

# Supplementary table 6 ---------------------------------------------------

files<-list.files("results/LiverFatMR/",pattern="_F_stats.csv")

st_7<-data.frame()

for (a in files){
  df<-fread(paste("results/LiverFatMR/",a,sep=""))
  df$analysis<-a
  st_7<-rbind(st_7,df)
}

st_7$analysis<-gsub("ASAT_FI_EC_F_stats.csv","Exposure 1: ASAT; Exposure 2: fasting insulin",st_7$analysis)
st_7$analysis<-gsub("ASAT_FI_EEC_F_stats.csv","Exposure 1: ASAT; Exposure 2: fasting insulin",st_7$analysis)
st_7$analysis<-gsub("ASAT_SHBG_EC_F_stats.csv","Exposure 1: ASAT; Exposure 2: SHBG",st_7$analysis)
st_7$analysis<-gsub("ASAT_SHBG_EEC_F_stats.csv","Exposure 1: ASAT; Exposure 2: SHBG",st_7$analysis)
st_7$analysis<-gsub("GFAT_AD_NEC_F_stats.csv","Exposure 1: GFAT; Exposure 2: adiponectin",st_7$analysis)
st_7$analysis<-gsub("ASAT_HDL_TNBC_F_stats.csv","Exposure 1: ASAT; Exposure 2: HDL cholesterol",st_7$analysis)
st_7$analysis<-gsub("ASAT_CXC_OC_F_stats.csv","Exposure 1: ASAT; Exposure 2: CXCL-8",st_7$analysis)
st_7$analysis<-gsub("ASAT_IGFBP_OAC_F_stats.csv","Exposure 1: ASAT; Exposure 2: IGFBP-1",st_7$analysis)
st_7$analysis<-gsub("ASAT_BT_EC_F_stats.csv","Exposure 1: ASAT; Exposure 2:Bioavailable testosterone",st_7$analysis)
st_7$analysis<-gsub("ASAT_BT_EEC_F_stats.csv","Exposure 1: ASAT; Exposure 2:Bioavailable testosterone",st_7$analysis)
st_7<-distinct(st_7)

st_7<-dplyr::select(st_7,analysis,exposure1,exposure2)
st_7$exposure1<-round(st_7$exposure1,digits=2)
st_7$exposure2<-round(st_7$exposure2,digits=2)

colnames(st_7)<-c("Analysis","Exposure 1 conditional F-statistic","Exposure 2 conditional F-statistic")

st_7<-distinct(st_7)


# Supplementary table 7 ---------------------------------------------------

st_8b<- fread("results/LiverFatMR/MVMR_sensitivity analyses/ASAT_HDL_nooverlap_TNBC_MVMR_results.csv")
st_8<-fread("results/LiverFatMR/MVMR_sensitivity analyses/ASAT_HDL_TNBC_mediation_results.csv")


st_8$`Unadjusted OR`<-tb_2$`Odds ratio`[tb_2$Exposure=="ASAT" & tb_2$Outcome=="Triple negative or basal-like breast cancer"]
st_8$`Unadjusted 95% confidence interval`<-tb_2$`95% confidence interval`[tb_2$Exposure=="ASAT" & tb_2$Outcome=="Triple negative or basal-like breast cancer"]

st_8$`Adjusted OR`<-exp(st_8b$Estimate[1])
st_8$`Adjusted LCI`<-exp(st_8b$Estimate[1]-st_8b$`Std. Error`[1]*1.96)
st_8$`Adjusted UCI`<-exp(st_8b$Estimate[1]+st_8b$`Std. Error`[1]*1.96)
st_8$`Adjusted 95% confidence interval`<-paste((round(st_8$`Adjusted LCI`,digits=2))," to ",(round(st_8$`Adjusted UCI`,digits=2)),sep="")


st_8$analysis<-gsub("ASAT_HDL_TNBC","Exposure 1: ASAT; Exposure 2: HDL cholesterol (no UK Biobank); Outcome: Triple negative or basal-like breast cancer",st_8$analysis)

st_8$percent<-round(st_8$percent*100,digits=2)
st_8$`95% confidence interval`<-paste((round(st_8$LCI*100,digits=2))," to ",(round(st_8$UCI*100,digits=2)),sep="")

st_8$p<-formatC(st_8$p, format = "e", digits = 2)

st_8<-dplyr::select(st_8,analysis,`Unadjusted OR`,`Unadjusted 95% confidence interval`,`Adjusted OR`,`Adjusted 95% confidence interval`,percent,`95% confidence interval`,p)
colnames(st_8)<-c("Analysis","Unadjusted OR","Unadjusted 95% confidence interval","Adjusted OR","Adjusted 95% confidence interval","Percent mediated","95% confidence interval for percent mediated","P-value for percent mediated")

# Supplementary table 8 ---------------------------------------------------
files<-list.files("results/LiverFatMR/MVMR_sensitivity analyses/",pattern="_F_stats.csv")

st_9<-data.frame()

for (a in files){
  df<-fread(paste("results/LiverFatMR/MVMR_sensitivity analyses/",a,sep=""))
  df$analysis<-a
  st_9<-rbind(st_9,df)
}

st_9$analysis<-gsub("ASAT_HDL_TNBC_F_stats.csv","Exposure 1: ASAT; Exposure 2: HDL cholesterol (no UK Biobank)",st_9$analysis)
st_9<-distinct(st_9)

st_9<-dplyr::select(st_9,analysis,exposure1,exposure2)
st_9$exposure1<-round(st_9$exposure1,digits=2)
st_9$exposure2<-round(st_9$exposure2,digits=2)

colnames(st_9)<-c("Analysis","Exposure 1 conditional F-statistic","Exposure 2 conditional F-statistic")

st_9<-distinct(st_9)


# Write supplementary tables ----------------------------------------------

xl_lst <- list('Supplementary table 1' = st_1, 'Supplementary table 2' = st_2, 'Supplementary table 3' = st_6, 'Supplementary table 4' = st_3, 'Supplementary table 5' = st_4, 'Supplementary table 6' = st_7, 'Supplementary table 7' = st_8, 'Supplementary table 8' = st_9)

write.xlsx(xl_lst, file = "Supplementary tables.xlsx")



