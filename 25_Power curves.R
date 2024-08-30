library(dplyr)
library(ggplot2)
library(readxl)

dat<-read_excel("Power calculations liver fat MRs.xlsx")

dat$Power<-as.numeric(dat$Power)

dat<-dat%>%arrange(desc(Power))
dat_vector<-dat[dat$OR==max(dat$OR),]
dat_vector<-dat_vector%>%arrange(desc(Power))
vector<-dat_vector$Outcome
vector<-unique(vector)
dat$Outcome<-factor(dat$Outcome,levels=vector)

dat$Exposure<-factor(dat$Exposure,levels=c("ASAT","VAT","GFAT","Pancreas fat","Liver fat"))

p<-ggplot(data=dat,aes(x=OR,y=Power,color=Outcome)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=80,linetype="dashed") +
  scale_x_continuous(breaks=c(1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))+
  scale_y_continuous(limits=c(0,100))+
  facet_wrap(~Exposure, ncol=2,scales = "free")
p

ggsave(filename="Liver Fat MR/Paper/Supplementary figure 1.png", plot=last_plot(),width = 220, height = 200, units = "mm",bg="white",dpi=1000)

p<-ggplot(data=dat,aes(x=OR,y=Power,color=Outcome)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept=80,linetype="dashed") +
  scale_x_continuous(breaks=c(1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2))+
  scale_y_continuous(limits=c(0,100))+
  facet_wrap(~Exposure, ncol=2,scales = "free") + 
  theme(axis.text.x = element_text(angle = 45)) +
  xlab("Odds ratio")
p

ggsave(filename="Liver Fat MR/Paper/Supplementary figure 1 wide.png", plot=last_plot(),width = 0.7*220, height = 0.7*240, units = "mm",bg="white",dpi=1000)
