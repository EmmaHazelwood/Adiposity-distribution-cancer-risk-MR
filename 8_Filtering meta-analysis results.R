#Used equation from supplementary note of this paper - https://www.nature.com/articles/ng.3538#MOESM39 to estimate effects from Z scores

library(data.table)
library(dplyr)

#Liver cancer
meta1<-fread("data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Liver_cancer/Meta_analysis_Liver_cancer_GWAS_formatted.csv")

#Gallbladder
meta1<-fread("data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Gallbladder_cancer/Meta_analysis_Gallbladder_cancer_GWAS_formatted.csv")

#Kidney
meta1<-fread("data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Kidney_cancer/Meta_analysis_Kidney_cancer_GWAS_formatted.csv")

#Thyroid
meta1<-fread("data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Thyroid_cancer/Meta_analysis_Thyroid_cancer_GWAS_formatted.csv")

#Multiple_myeloma
meta1<-fread("data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Multiple_myeloma/Meta_analysis_Multiple_myeloma_cancer_GWAS_formatted.csv")

#Pancreas
meta1<-fread("data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_GWAS1.csv")
meta2<-fread("data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_GWAS_GC1.csv")
meta1<-select(meta1,-c("Effect","StdErr"))
meta2<-select(meta2,MarkerName,Effect,StdErr)
meta<-merge(meta1,meta2,by="MarkerName")
meta<-meta[meta$HetPVal>0.05,]
meta$beta<-meta$Effect
meta$se<-meta$StdErr
meta$Allele1<-toupper(meta$Allele1)
meta$Allele2<-toupper(meta$Allele2)
fwrite(meta,"data/Obesity_related_cancer_GWAS/Pancreas_cancer/Meta_analysis_Pancreas_cancer_GWAS_formatted.csv")
