#Kidney
library(data.table)

GWAS<-fread("data/Obesity_related_cancer_GWAS/Kidney_cancer/C_KIDNEY_NOTRENALPELVIS.gwas.imputed_v3.both_sexes.tsv")
var<-fread("data/Neale_lab_annotation_file_2018/variants.tsv")

gwas<-merge(GWAS,var)

u=1114/463010
  
gwas$logOR<-gwas$beta/(u * (1 - u))
gwas$logOR_se<-gwas$se/(u * (1 - u))

fwrite(gwas,"data/Obesity_related_cancer_GWAS/Kidney_cancer/C_KIDNEY_NOTRENALPELVIS.gwas.imputed_v3.both_sexes_formatted.tsv")

