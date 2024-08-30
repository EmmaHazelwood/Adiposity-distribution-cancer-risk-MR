library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(gwascat)

kgenomes<-fread("data/1000GenomesReferenceFiles/EUR.bim") #grch37

colnames(kgenomes) <- c("Chromosome","SNP","NA","Position","Allele1","Allele2")
kgenomes <- kgenomes[, c("Chromosome", "Position", "Allele1", "Allele2", "SNP")]

out_dat<-fread("data/Sinnott_Armstrong_et_al_2019/IGF_1.imp")
colnames(out_dat)[1]<-"CHROM"

out_dat<-out_dat[out_dat$CHROM==12,]
out_dat<-out_dat[out_dat$POS>=102289652,] #500mb of the gene coding window
out_dat<-out_dat[out_dat$POS<=103374341,]

out_dat$Pvalue<-as.numeric(out_dat$`P-value`)
out_dat<-out_dat[out_dat$Pvalue<=5*10^-8,]

rsid1<-merge(out_dat,kgenomes,by.x=c("CHROM","POS","REF","ALT"),by.y=c("Chromosome","Position","Allele1","Allele2"))
rsid2<-merge(out_dat,kgenomes,by.x=c("CHROM","POS","ALT","REF"),by.y=c("Chromosome","Position","Allele1","Allele2"))

rsid<-rbind(rsid1,rsid2)

fwrite(rsid,"data/Sinnott_Armstrong_et_al_2019/IGF_1_withrsID.imp")

out_dat<-fread("data/Sinnott_Armstrong_et_al_2019/C_reactive_protein.imp",sep = "\t")
colnames(out_dat)[1]<-"CHROM"

out_dat<-out_dat[out_dat$CHROM==1,]
out_dat<-out_dat[out_dat$POS>=159182079,] #500mb of the gene coding window
out_dat<-out_dat[out_dat$POS<=160184379,]

out_dat$Pvalue<-as.numeric(out_dat$`P-value`)
out_dat<-out_dat[out_dat$Pvalue<=5*10^-8,]

rsid1<-merge(out_dat,kgenomes,by.x=c("CHROM","POS","REF","ALT"),by.y=c("Chromosome","Position","Allele1","Allele2"))
rsid2<-merge(out_dat,kgenomes,by.x=c("CHROM","POS","ALT","REF"),by.y=c("Chromosome","Position","Allele1","Allele2"))

rsid<-rbind(rsid1,rsid2)

fwrite(rsid,"data/Sinnott_Armstrong_et_al_2019/C_reactive_protein_withrsID.imp")
