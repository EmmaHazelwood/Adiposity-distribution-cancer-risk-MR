sink("formatting.txt")
library(dplyr)
library(TwoSampleMR)
library(data.table)
library(biomaRt)
library(tidyr)
library(ieugwasr)

#New approach - find position and chromosome of all exposure SNPs, and work out 100kb window around them
#Filter outcome data for only SNPs in these regions (in a loop)
#Use biomart to find rsID for these SNPs (and GWS SNPs)

#Get SNPs we need for exposures
exposure_dat1 <- extract_instruments(outcomes=c("ebi-a-GCST90016673","ebi-a-GCST90016675"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat2 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat2$samplesize.exposure<-NA
exposure_dat2$exposure<-"ASAT"

exposure_dat3 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat3$samplesize.exposure<-NA
exposure_dat3$exposure<-"VAT"

exposure_dat4 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat4$samplesize.exposure<-NA
exposure_dat4$exposure<-"GFAT"

exposure_dat2<-rbind(exposure_dat2,exposure_dat3,exposure_dat4)

#Get position and chromosome by rsID
ensembl <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                      host="grch37.ensembl.org", #Open GWAS and outcome also gr37
                      dataset="hsapiens_snp")

results<-data.frame()
for (j in 1:length(unique(exposure_dat2$SNP))){
  print(j)
  SNP<-unique(exposure_dat2$SNP)[j]
  res<-getBM(attributes = c('refsnp_id', 'chrom_start','chr_name'), 
             filters = 'snp_filter', 
             values = SNP, 
             mart = ensembl,
             useCache = FALSE)
  results<-rbind(results,res)
}

exposure_dat2<-merge(exposure_dat2,results,by.x="SNP",by.y="refsnp_id")
colnames(exposure_dat2)[colnames(exposure_dat2)=="chrom_start"]<-"pos.exposure"
colnames(exposure_dat2)[colnames(exposure_dat2)=="chr_name"]<-"chr.exposure"

exposure_dat2$chr.exposure[exposure_dat2$chr.exposure=="HSCHR6_MHC_COX"]<-6
exposure_dat2$chr.exposure[exposure_dat2$chr.exposure=="HSCHR6_MHC_QBL"]<-6

exposure_dat<-rbind(exposure_dat1,exposure_dat2)
SNPs<-dplyr::select(exposure_dat,SNP,chr.exposure,pos.exposure)

outcome_dat<-fread("data/Sinnott_Armstrong_et_al_2019/IGF_1.imp")
colnames(outcome_dat)[1]<-"CHROM"
SNPs$chr.exposure<-as.numeric(SNPs$chr.exposure)
SNPs$pos.exposure<-as.numeric(SNPs$pos.exposure)
outcome_dat<-merge(outcome_dat,SNPs,by.x=c("CHROM","POS"),by.y=c("chr.exposure","pos.exposure"))
outcome_dat_exposuresnps<-outcome_dat

outcome_dat<-data.frame(outcome_dat)
outcome_dat<-format_data(outcome_dat,
                         type="outcome",
                         snps=exposure_dat$SNP,
                         snp_col="SNP",
                         beta_col = "Effect",
                         se_col = "StdErr",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF")

#Harmonise to find missing/palindromic SNPs that might need proxies
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat<-dat[dat$mr_keep==TRUE,]

#Missing SNPs
snps<-setdiff(exposure_dat$SNP,dat$SNP)

exposure_dat<-exposure_dat[exposure_dat$SNP %in% snps,]

exposure_dat$minpos<-exposure_dat$pos.exposure-500
exposure_dat$maxpos<-exposure_dat$pos.exposure+500

outcome_dat<-fread("data/Sinnott_Armstrong_et_al_2019/IGF_1.imp")
colnames(outcome_dat)[1]<-"CHROM"
outcome_dat$`P-value`<-as.numeric(outcome_dat$`P-value`)
df<-outcome_dat[outcome_dat$`P-value`<0.00000005,]

#Get cis SNPs to limit numbers we're looking up
df<-df[df$CHROM==12,]
df<-df[df$POS>=102689652,]
df<-df[df$POS>=102974341,]

dat<-data.frame()
for (i in 1:nrow(exposure_dat)){
  print(i)
  dat2<-exposure_dat[i,]
  outcome_dat2<-outcome_dat[outcome_dat$CHROM==dat2$chr.exposure,]
  outcome_dat2<-outcome_dat2[outcome_dat2$POS<=dat2$maxpos,]
  outcome_dat2<-outcome_dat2[outcome_dat2$POS>=dat2$minpos,]
  dat<-rbind(dat,outcome_dat2)
  
}

dat<-rbind(dat,df)

df<-dplyr::select(dat,CHROM,POS)
df$POS2<-df$POS

## combine the positions in to a single vector
position <- apply(df, 1, paste, collapse = ":")

results<-data.frame()

for (a in 1:length(position)){
  print(a)
  SNP<-position[a]
  res<-getBM(attributes = c('refsnp_id', 'allele', 'chrom_start','chr_name'), 
             filters = 'chromosomal_region', 
             values = SNP, 
             mart = ensembl,
             useCache = FALSE)
  results<-rbind(results,res)
}

results<-results%>% separate(allele, c("allele1", "allele2"))

outcome_dat$chr_name<-as.numeric(outcome_dat$chr_name)
outcome_dat$chrom_start<-as.numeric(outcome_dat$chrom_start)
outcome_dat1<-merge(outcome_dat,results,by.x=c("CHROM","POS","REF"),by.y=c("chr_name","chrom_start","allele1"))
outcome_dat1<-outcome_dat1[,!"allele2"]
outcome_dat2<-merge(outcome_dat,results,by.x=c("CHROM","POS","REF"),by.y=c("chr_name","chrom_start","allele2"))
outcome_dat2<-outcome_dat2[,!"allele1"]
outcome_dat<-rbind(outcome_dat1,outcome_dat2)
outcome_dat$SNP<-outcome_dat$refsnp_id
outcome_dat<-outcome_dat[,-c(16,17,18)]

outcome_dat<-rbind(outcome_dat,outcome_dat_exposuresnps)

fwrite(outcome_dat,"data/Sinnott_Armstrong_et_al_2019/IGF_1_Instruments.txt")

#CRP
#Get SNPs we need for exposures
exposure_dat1 <- extract_instruments(outcomes=c("ebi-a-GCST90016673","ebi-a-GCST90016675"), p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)
exposure_dat2 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/ASAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat2$samplesize.exposure<-NA
exposure_dat2$exposure<-"ASAT"

exposure_dat3 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/VAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)
exposure_dat3$samplesize.exposure<-NA
exposure_dat3$exposure<-"VAT"

exposure_dat4 <- read_exposure_data(
  filename = "data/Measures_of_adiposity_Agrawal_2022/GFAT.uvinput",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure"
)

exposure_dat4$samplesize.exposure<-NA
exposure_dat4$exposure<-"GFAT"

exposure_dat2<-rbind(exposure_dat2,exposure_dat3,exposure_dat4)

#Get position and chromosome by rsID
ensembl <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
                      host="grch37.ensembl.org", #Open GWAS and outcome also gr37
                      dataset="hsapiens_snp")

results<-data.frame()
for (j in 1:length(unique(exposure_dat2$SNP))){
  print(j)
  SNP<-unique(exposure_dat2$SNP)[j]
  res<-getBM(attributes = c('refsnp_id', 'chrom_start','chr_name'), 
             filters = 'snp_filter', 
             values = SNP, 
             mart = ensembl,
             useCache = FALSE)
  results<-rbind(results,res)
}

exposure_dat2<-merge(exposure_dat2,results,by.x="SNP",by.y="refsnp_id")
colnames(exposure_dat2)[colnames(exposure_dat2)=="chrom_start"]<-"pos.exposure"
colnames(exposure_dat2)[colnames(exposure_dat2)=="chr_name"]<-"chr.exposure"

exposure_dat2$chr.exposure[exposure_dat2$chr.exposure=="HSCHR6_MHC_COX"]<-6
exposure_dat2$chr.exposure[exposure_dat2$chr.exposure=="HSCHR6_MHC_QBL"]<-6

exposure_dat<-rbind(exposure_dat1,exposure_dat2)
SNPs<-dplyr::select(exposure_dat,SNP,chr.exposure,pos.exposure)

outcome_dat<-fread("data/Sinnott_Armstrong_et_al_2019/C_reactive_protein.imp")
colnames(outcome_dat)[1]<-"CHROM"
SNPs$chr.exposure<-as.numeric(SNPs$chr.exposure)
SNPs$pos.exposure<-as.numeric(SNPs$pos.exposure)
outcome_dat<-merge(outcome_dat,SNPs,by.x=c("CHROM","POS"),by.y=c("chr.exposure","pos.exposure"))
outcome_dat_exposuresnps<-outcome_dat

outcome_dat<-data.frame(outcome_dat)
outcome_dat<-format_data(outcome_dat,
                         type="outcome",
                         snps=exposure_dat$SNP,
                         snp_col="SNP",
                         beta_col = "Effect",
                         se_col = "StdErr",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF")

#Harmonise to find missing/palindromic SNPs that might need proxies
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat<-dat[dat$mr_keep==TRUE,]

#Missing SNPs
snps<-setdiff(exposure_dat$SNP,dat$SNP)

exposure_dat<-exposure_dat[exposure_dat$SNP %in% snps,]

exposure_dat$minpos<-exposure_dat$pos.exposure-500
exposure_dat$maxpos<-exposure_dat$pos.exposure+500

outcome_dat<-fread("data/Sinnott_Armstrong_et_al_2019/C_reactive_protein.imp")
colnames(outcome_dat)[1]<-"CHROM"
outcome_dat$`P-value`<-as.numeric(outcome_dat$`P-value`)
df<-outcome_dat[outcome_dat$`P-value`<0.00000005,]

#Get cis SNPs to limit numbers we're looking up
df<-df[df$CHROM==12,]
df<-df[df$POS>=102689652,]
df<-df[df$POS>=102974341,]

dat<-data.frame()
for (i in 1:nrow(exposure_dat)){
  print(i)
  dat2<-exposure_dat[i,]
  outcome_dat2<-outcome_dat[outcome_dat$CHROM==dat2$chr.exposure,]
  outcome_dat2<-outcome_dat2[outcome_dat2$POS<=dat2$maxpos,]
  outcome_dat2<-outcome_dat2[outcome_dat2$POS>=dat2$minpos,]
  dat<-rbind(dat,outcome_dat2)
  
}

dat<-rbind(dat,df)

df<-dplyr::select(dat,CHROM,POS)
df$POS2<-df$POS

## combine the positions in to a single vector
position <- apply(df, 1, paste, collapse = ":")

results<-data.frame()

for (a in 1:length(position)){
  print(a)
  SNP<-position[a]
  res<-getBM(attributes = c('refsnp_id', 'allele', 'chrom_start','chr_name'), 
             filters = 'chromosomal_region', 
             values = SNP, 
             mart = ensembl,
             useCache = FALSE)
  results<-rbind(results,res)
}

results<-results%>% separate(allele, c("allele1", "allele2"))

outcome_dat$chr_name<-as.numeric(outcome_dat$chr_name)
outcome_dat$chrom_start<-as.numeric(outcome_dat$chrom_start)
outcome_dat1<-merge(outcome_dat,results,by.x=c("CHROM","POS","REF"),by.y=c("chr_name","chrom_start","allele1"))
outcome_dat1<-outcome_dat1[,!"allele2"]
outcome_dat2<-merge(outcome_dat,results,by.x=c("CHROM","POS","REF"),by.y=c("chr_name","chrom_start","allele2"))
outcome_dat2<-outcome_dat2[,!"allele1"]
outcome_dat<-rbind(outcome_dat1,outcome_dat2)
outcome_dat$SNP<-outcome_dat$refsnp_id
outcome_dat<-outcome_dat[,-c(16,17,18)]

outcome_dat<-rbind(outcome_dat,outcome_dat_exposuresnps)

fwrite(outcome_dat,"data/Sinnott_Armstrong_et_al_2019/C_reactive_protein_Instruments.txt")


print("fin")
