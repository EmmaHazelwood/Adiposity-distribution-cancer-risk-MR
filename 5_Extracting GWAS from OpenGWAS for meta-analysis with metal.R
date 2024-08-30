#Run on linux server:
#wget https://objectstorage.us-ashburn-1.oraclecloud.com/n/idrvm4tkz2a8/b/OpenGWAS/o/finn-b/finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC/finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC.vcf.gz

#remotes::install_github("mrcieu/gwasvcf")
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)


remotes::install_github('mrcieu/genetics.binaRies')
set_plink()
set_bcftools()

vcf <- readVcf("data/Obesity_related_cancer_GWAS/Kidney_cancer/finngen/finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC.vcf.gz")

#vcf_to_granges(vcf) %>% dplyr::as_tibble()

gwas <- gwasglue::gwasvcf_to_TwoSampleMR(vcf, type="exposure")

library(data.table)
fwrite(gwas,"data/Obesity_related_cancer_GWAS/Kidney_cancer/finngen/finn-b-C3_KIDNEY_NOTRENALPELVIS_EXALLC.csv")
