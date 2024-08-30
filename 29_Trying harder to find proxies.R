library(dplyr)
library(TwoSampleMR)
library(data.table)
library(LDlinkR)
library(tidyr)


#Read in 1000Genomes file for proxy SNPs
kg<-fread("data/1000GenomesReferenceFiles/EUR.bim")
kgsnps<-kg$V2

results<-data.frame(matrix(nrow=0,ncol=0))


# BMI -> CRP --------------------------------------------------------------
exposure_dat <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)

outfilepath<-"data/Sinnott_Armstrong_et_al_2019/C_reactive_protein_Instruments.txt"
print("CRP")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "ALT",
  other_allele_col = "REF"
)

#ask if proxy SNPs needed
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  #Find which SNPs needed
  new_outcome_dat<-data.frame(matrix(ncol=29,nrow=0))
  missing<-setdiff(exposure_dat$SNP,outcome_dat$SNP)
  
  #Limit missing SNPs to those in 1000 genomes (as need to be in order to look up LD structure)
  missing<-missing[missing %in% kgsnps]
  
  #Find proxy for each missing SNP
  for (j in missing){
    #Get proxies
    print(j)
    rm(proxy)
    
    #Find all proxy SNPs using 1000 genomes
    proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
    
    #Limit to those with R2 > 0.8
    proxy <- proxy %>% dplyr::filter(R2 > 0.8) # minimum R2 0.8    
    
    #Check if any of these SNPs are present in outcome data
    out<-fread(outfilepath)
    
    #If none of these SNPs are in outcome data, leave it there and go onto next SNP. If not...
    if(length(intersect(out$MarkerName,proxy$RS_Number))!=0){
      #Filter for those in outcome data
      outcome_dat_proxy <- read_outcome_data(
        snps = proxy$RS_Number,
        filename = outfilepath,
        sep = ",",
        snp_col = "SNP",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "ALT",
        other_allele_col = "REF"
      )
      
      proxy<-proxy[proxy$RS_Number %in% outcome_dat_proxy$SNP,]
      #subset for those with compatible alleles
      #For each potential proxy SNP
      for (a in proxy$RS_Number){
        #Limit to just that SNP
        proxy2<-proxy[proxy$RS_Number==a,]
        #Make the alleles into a dataframe
        Correlated_Alleles<-data.frame(proxy2$Correlated_Alleles)
        Correlated_Alleles$proxy2.Correlated_Alleles<- as.character(Correlated_Alleles$proxy2.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy2.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        #a1exp = original SNP, allele1, a1out = proxy SNP, allele1
        #a2exp = original SNP, allele2, a2out = proxy SNP, allele2
        
        #Make sure the alleles available in the outcome dataset are the same as those needed for the proxy
        new_outcome_dat<-outcome_dat_proxy[outcome_dat_proxy$SNP==proxy2$RS_Number,]
        original_alleles<-c(new_outcome_dat$effect_allele.outcome,new_outcome_dat$other_allele.outcome)
        new_alleles<-c(Correlated_Alleles$a1out,Correlated_Alleles$a2out)
        
        #If the alleles don't match up between proxy and those available in outcome data drop this proxy SNP
        if(length(intersect(original_alleles,new_alleles))!=2){
          proxy<-proxy[-which(proxy$RS_Number==a),]
        }}
      
      #Get the proxy SNP with highest R2 that remains and match up with original alleles
      proxy<-proxy[1,]
      Correlated_Alleles<-data.frame(proxy$Correlated_Alleles)
      Correlated_Alleles$proxy.Correlated_Alleles<- as.character(Correlated_Alleles$proxy.Correlated_Alleles)
      Correlated_Alleles<-Correlated_Alleles %>% separate(proxy.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
      Correlated_Alleles<-na.omit(Correlated_Alleles)
      
      new_outcome_dat<-outcome_dat_proxy[outcome_dat_proxy$SNP==proxy$RS_Number,]
      original_alleles<-c(new_outcome_dat$effect_allele.outcome,new_outcome_dat$other_allele.outcome)
      new_alleles<-c(Correlated_Alleles$a1out,Correlated_Alleles$a2out)
      
      #If alleles compatible
      if(length(intersect(original_alleles,new_alleles))==2){
        #Call the proxy SNP the RS number of the original SNP
        new_outcome_dat$SNP<-j
        
        #Check for both combinations (if a1 and a2 are different way round in proxy and original)
        if (new_outcome_dat$effect_allele.outcome==Correlated_Alleles$a1out){
          new_outcome_dat$effect_allele.outcome<-Correlated_Alleles$a1exp
          new_outcome_dat$other_allele.outcome<-Correlated_Alleles$a2exp
        }
        if (new_outcome_dat$other_allele.outcome==Correlated_Alleles$a1out){
          new_outcome_dat$other_allele.outcome<-Correlated_Alleles$a1exp
          new_outcome_dat$effect_allele.outcome<-Correlated_Alleles$a2exp
        }
        
        new_outcome_dat$id.outcome<-unique(outcome_dat$id.outcome)
        
        #Add the proxy SNP to the outcome data
        outcome_dat<-rbind(outcome_dat,new_outcome_dat)
      }
      
    }
  }
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-318271
dat <- steiger_filtering(dat)
dat$outcome <- "CRP"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)




# BMI -> IGF-1 ------------------------------------------------------------
exposure_dat <- extract_instruments(outcomes="ieu-b-40", p1 = 5e-08, clump = TRUE, p2 = 5e-08, r2 = 0.001, kb = 10000)

outfilepath<-"data/Sinnott_Armstrong_et_al_2019/IGF_1_Instruments.txt"
print("IGF1")
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = outfilepath,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "ALT",
  other_allele_col = "REF"
)
#ask if proxy SNPs needed
if (nrow(exposure_dat)!=nrow(outcome_dat)){
  #Find which SNPs needed
  new_outcome_dat<-data.frame(matrix(ncol=29,nrow=0))
  missing<-setdiff(exposure_dat$SNP,outcome_dat$SNP)
  
  #Limit missing SNPs to those in 1000 genomes (as need to be in order to look up LD structure)
  missing<-missing[missing %in% kgsnps]
  
  #Find proxy for each missing SNP
  for (j in missing){
    #Get proxies
    print(j)
    rm(proxy)
    
    #Find all proxy SNPs using 1000 genomes
    proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
    
    #Limit to those with R2 > 0.8
    proxy <- proxy %>% dplyr::filter(R2 > 0.8) # minimum R2 0.8    
    
    #Check if any of these SNPs are present in outcome data
    out<-fread(outfilepath)
    
    #If none of these SNPs are in outcome data, leave it there and go onto next SNP. If not...
    if(length(intersect(out$MarkerName,proxy$RS_Number))!=0){
      #Filter for those in outcome data
      outcome_dat_proxy <- read_outcome_data(
        snps = proxy$RS_Number,
        filename = outfilepath,
        sep = ",",
        snp_col = "SNP",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "ALT",
        other_allele_col = "REF"
      )
      
      proxy<-proxy[proxy$RS_Number %in% outcome_dat_proxy$SNP,]
      #subset for those with compatible alleles
      #For each potential proxy SNP
      for (a in proxy$RS_Number){
        #Limit to just that SNP
        proxy2<-proxy[proxy$RS_Number==a,]
        #Make the alleles into a dataframe
        Correlated_Alleles<-data.frame(proxy2$Correlated_Alleles)
        Correlated_Alleles$proxy2.Correlated_Alleles<- as.character(Correlated_Alleles$proxy2.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy2.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        #a1exp = original SNP, allele1, a1out = proxy SNP, allele1
        #a2exp = original SNP, allele2, a2out = proxy SNP, allele2
        
        #Make sure the alleles available in the outcome dataset are the same as those needed for the proxy
        new_outcome_dat<-outcome_dat_proxy[outcome_dat_proxy$SNP==proxy2$RS_Number,]
        original_alleles<-c(new_outcome_dat$effect_allele.outcome,new_outcome_dat$other_allele.outcome)
        new_alleles<-c(Correlated_Alleles$a1out,Correlated_Alleles$a2out)
        
        #If the alleles don't match up between proxy and those available in outcome data drop this proxy SNP
        if(length(intersect(original_alleles,new_alleles))!=2){
          proxy<-proxy[-which(proxy$RS_Number==a),]
        }}
      
      #Get the proxy SNP with highest R2 that remains and match up with original alleles
      proxy<-proxy[1,]
      Correlated_Alleles<-data.frame(proxy$Correlated_Alleles)
      Correlated_Alleles$proxy.Correlated_Alleles<- as.character(Correlated_Alleles$proxy.Correlated_Alleles)
      Correlated_Alleles<-Correlated_Alleles %>% separate(proxy.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
      Correlated_Alleles<-na.omit(Correlated_Alleles)
      
      new_outcome_dat<-outcome_dat_proxy[outcome_dat_proxy$SNP==proxy$RS_Number,]
      original_alleles<-c(new_outcome_dat$effect_allele.outcome,new_outcome_dat$other_allele.outcome)
      new_alleles<-c(Correlated_Alleles$a1out,Correlated_Alleles$a2out)
      
      #If alleles compatible
      if(length(intersect(original_alleles,new_alleles))==2){
        #Call the proxy SNP the RS number of the original SNP
        new_outcome_dat$SNP<-j
        
        #Check for both combinations (if a1 and a2 are different way round in proxy and original)
        if (new_outcome_dat$effect_allele.outcome==Correlated_Alleles$a1out){
          new_outcome_dat$effect_allele.outcome<-Correlated_Alleles$a1exp
          new_outcome_dat$other_allele.outcome<-Correlated_Alleles$a2exp
        }
        if (new_outcome_dat$other_allele.outcome==Correlated_Alleles$a1out){
          new_outcome_dat$other_allele.outcome<-Correlated_Alleles$a1exp
          new_outcome_dat$effect_allele.outcome<-Correlated_Alleles$a2exp
        }
        
        new_outcome_dat$id.outcome<-unique(outcome_dat$id.outcome)
        
        #Add the proxy SNP to the outcome data
        outcome_dat<-rbind(outcome_dat,new_outcome_dat)
      }
      
    }
  }
}

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat,action=3)
dat$samplesize.outcome<-317114
dat <- steiger_filtering(dat)
dat$outcome <- "IGF1"
res <- mr(dat, method_list=c("mr_ivw", "mr_weighted_mode", "mr_weighted_median","mr_wald_ratio"))
results<-rbind(results,res)




