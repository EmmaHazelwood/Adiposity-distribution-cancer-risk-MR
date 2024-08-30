library(data.table)
library(dplyr)
library(TwoSampleMR)
library(LDlinkR)
library(tidyr)

#Liver cancer
meta1<-fread("data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_formatted_allSNPs.csv")

#Gallbladder
meta1<-fread("data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS_formatted_allSNPs.csv")

#Kidney
meta1<-fread("data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS_formatted_allSNPs.csv")

#Thyroid
meta1<-fread("data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS_formatted_allSNPs.csv")

#Multiple_myeloma
meta1<-fread("data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted_allSNPs.csv")

#Pancreas_cancer
meta1<-fread("data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_cancer_GWAS_formatted_allSNPs.csv")
