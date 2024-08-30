library(data.table)
library(TwoSampleMR)
library(dplyr)

var<-fread("data/Neale_lab_annotation_file_2018/variants.tsv")

#HDL cholesterol
GWAS<-fread("data/Neale_lab_2018/HDL_cholesterol/30760_irnt.gwas.imputed_v3.female.varorder.tsv")
gwas<-merge(GWAS,var,by="variant")
fwrite(gwas,"data/Neale_lab_2018/HDL_cholesterol/female_HDL_cholesterol.tsv")

gwas<-data.frame(gwas)
exposure_dat<-format_data(
  gwas,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref"
)

exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)

exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")

fwrite(exposure_dat,"data/Neale_lab_2018/HDL_cholesterol/female_HDL_cholesterol_Instrument.txt")

#Triglycerides
GWAS<-fread("data/Neale_lab_2018/Triglycerides/30870_irnt.gwas.imputed_v3.female.varorder.tsv")
gwas<-merge(GWAS,var,by="variant")
fwrite(gwas,"data/Neale_lab_2018/Triglycerides/female_Triglycerides.tsv")
gwas<-data.frame(gwas)
exposure_dat<-format_data(
  gwas,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref"
)

exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)

exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")

fwrite(exposure_dat,"data/Neale_lab_2018/Triglycerides/female_Triglycerides_Instrument.txt")

#IGF-1
GWAS<-fread("data/Neale_lab_2018/IGF1/30770_irnt.gwas.imputed_v3.female.varorder.tsv")
gwas<-merge(GWAS,var,by="variant")
fwrite(gwas,"data/Neale_lab_2018/IGF1/female_IGF1.tsv")

gwas<-gwas[gwas$pval<5e-8,]

gwas<-gwas%>%tidyr::separate(variant,c("CHROM","POS"),sep=":")

gwas$CHROM<-as.numeric(gwas$CHROM)
gwas$POS<-as.numeric(gwas$POS)

gwas<-gwas[gwas$CHROM==12,]
gwas<-gwas[gwas$POS>=102289652,] #500kb of the gene coding window
gwas<-gwas[gwas$POS<=103374341,]

gwas<-data.frame(gwas)
exposure_dat<-format_data(
  gwas,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref"
)


exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)

exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")

fwrite(exposure_dat,"data/Neale_lab_2018/IGF1/female_IGF1_Instrument.txt")



# Fasting insulin ---------------------------------------------------------
out_dat<-fread("data/FI_GWAS/FI_female_1000G_density_formatted_21-03-29.txt",sep = "\t")

out_dat<-data.frame(out_dat)
exposure_dat<-format_data(
  out_dat,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "a2",
  other_allele_col = "a1"
)

exposure_dat <- subset(exposure_dat, pval.exposure<5e-8)

exposure_dat <- clump_data(exposure_dat,
                           clump_kb = 10000,
                           clump_r2 = 0.001,
                           clump_p1 = 0.00000005,
                           clump_p2 = 0.00000005,
                           pop = "EUR")

fwrite(exposure_dat,"data/FI_GWAS/FI_female_Instrument.txt")


print("fin")