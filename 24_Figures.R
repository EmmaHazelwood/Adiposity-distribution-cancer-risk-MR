library(dplyr)
library(ggforestplot)
library(ggplot2)
library(data.table)
library(qqman)
library(biomaRt)
library(tidyr)




# Figure 2 ----------------------------------------------------------------
fg_2a<-fread("LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
fg_2b<-fread("LiverFatMR/BMI_obesity_cancers_results.csv")
fg_2<-rbind(fg_2a,fg_2b,fill=TRUE)
fg_2$exposure[fg_2$exposure=="Percent liver fat || id:ebi-a-GCST90016673"]<-"Liver fat"
fg_2$exposure[fg_2$exposure=="Pancreas fat || id:ebi-a-GCST90016675"]<-"Pancreas fat"
fg_2$exposure[fg_2$exposure=="body mass index || id:ieu-b-40"]<-"BMI"
fg_2$exposure<-factor(fg_2$exposure,levels=c("ASAT","VAT","GFAT","Liver fat","Pancreas fat","BMI"))
fg_2$outcome<-gsub("uminal B-HER2","uminal B/HER2",fg_2$outcome)

fg_2$BP[fg_2$pval<(0.05/12)]<-0
fg_2$BP[fg_2$pval>=(0.05/12)]<-1
fg_2<-fg_2 %>%
group_by(exposure) %>%
arrange(pval, .by_group = TRUE)
fg_2$method<-factor(fg_2$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","MR Egger", "Wald ratio")))
fg_2$type<-"Subtype"
fg_2$type[fg_2$outcome=="Esophagus adenocarcinoma"|fg_2$outcome=="Colorectal cancer (overall)"|fg_2$outcome=="Liver cancer"|fg_2$outcome=="Gallbladder cancer"|fg_2$outcome=="Pancreas cancer"|fg_2$outcome=="Breast cancer"|fg_2$outcome=="Endometrial cancer"|fg_2$outcome=="Ovarian cancer"|fg_2$outcome=="Kidney (renal-cell) cancer"|fg_2$outcome=="Thyroid cancer"|fg_2$outcome=="Multiple myeloma"|fg_2$outcome=="Meningioma"]<-"Overall"

fg_2<-fg_2[!(fg_2$method=="Weighted median" & fg_2$nsnp<10),]
fg_2<-fg_2[!(fg_2$method=="Weighted mode" & fg_2$nsnp<10),]

fg_2$outcome<-gsub("Esophagus adenocarcinoma","Oesophageal adenocarcinoma",fg_2$outcome)
fg_2$outcome[fg_2$outcome=="Colorectal cancer (overall)"]<-"Colorectal cancer"
fg_2$method<-factor(fg_2$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","MR Egger")))

fg_2_all<-fg_2


#Overall
fg_2<-fg_2_all[fg_2_all$type=="Overall",]
fg_2<-ungroup(fg_2)

p1<- ggforestplot::forestplot(
df = fg_2,
name=exposure,
estimate = b,
pvalue = pval,
psignif = 0.05,
xlab = "Odds ratio (95% CI)",
logodds=TRUE,
se=se,
colour = method
)+
ggplot2::facet_wrap(vars(outcome)) +
scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))
#Rewrite code above but remove legend
p1<- ggforestplot::forestplot(
  df = fg_2,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome)) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "none", text=element_text(size=17))+
  coord_cartesian(xlim=c(0.05,15))

p1

p1_poster<- ggforestplot::forestplot(
  df = fg_2,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=2) +
  scale_color_manual(values=rev(c("#CA4275","#DFB845","#4A9A9A")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "bottom", text=element_text(size=17))+
  coord_cartesian(xlim=c(0.05,15))

p1_poster
ggplot2::ggsave(filename="Figure 2 poster.png", plot=ggplot2::last_plot(),width = 0.45*380, height = 0.44*580, units = "mm",bg="white",dpi=1000)
ggplot2::ggsave(filename="Figure 2 poster legend.png", plot=ggplot2::last_plot(),width = 0.7*380, height = 0.43*580, units = "mm",bg="white",dpi=1000)


p1

#Just IVW and certain cancers
fg_2_ivw<-fg_2[fg_2$method=="Inverse variance weighted",]
fg_2_ivw<-fg_2_ivw[fg_2_ivw$outcome=="Oesophageal adenocarcinoma"|fg_2_ivw$outcome=="Liver cancer"|fg_2_ivw$outcome=="Endometrial cancer"|fg_2_ivw$outcome=="Breast cancer",]

p1_ivw<- ggforestplot::forestplot(
  df = fg_2_ivw,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),nrow=1) +
  theme(legend.position = "none", text=element_text(size=17))

p1_ivw

ggplot2::ggsave(filename="Figure 2 wide ivw.png", plot=ggplot2::last_plot(),width = 0.9*0.5*0.8*1000, height = 0.4*0.8*0.7*480, units = "mm",bg="white",dpi=1000)


#subtypes
fg_2<-fg_2_all[fg_2_all$type=="Subtype",]
fg_2<-ungroup(fg_2)
fg_2$outcome<-factor(fg_2$outcome,levels=c("Luminal B/HER2-negative-like breast cancer","Triple negative or basal-like breast cancer","Luminal B-like breast cancer","Luminal A-like breast cancer", "HER2-enriched-like breast cancer",
 "Endometrioid ovarian cancer","High grade serous carcinoma ovarian cancer", "Clear cell ovarian cancer","Low malignant potential ovarian cancer","Invasive mucinous ovarian cancer","Low grade serous carcinoma ovarian cancer",
 "Endometrioid endometrial cancer","Non-endometrioid endometrial cancer",
 "Colon cancer","Distal colon cancer","Proximal colon cancer","Rectal cancer"))

p2<- ggforestplot::forestplot(
df = fg_2,
name=exposure,
estimate = b,
pvalue = pval,
psignif = 0.05,
xlab = "Odds ratio (95% CI)",
logodds=TRUE,
se=se,
colour = method
)+
ggplot2::facet_wrap(vars(outcome),ncol=4) +
scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "bottom", text=element_text(size=17))

p2

library(ggpubr)
fg_2<-ggarrange(p1,p2,ncol=1,labels=c("A","B"),heights=c(3,5))
fg_2

ggplot2::ggsave(filename="Figure 2.png", plot=ggplot2::last_plot(),width = 480, height = 520, units = "mm",bg="white",dpi=1000)
ggplot2::ggsave(filename="Figure 2 wide.png", plot=ggplot2::last_plot(),width = 0.9*0.5*0.8*1000, height = 0.8*0.7*480, units = "mm",bg="white",dpi=1000)




#Wider subtypes plot for presentations
fg_2<-fg_2_all[fg_2_all$type=="Subtype",]
fg_2<-ungroup(fg_2)
fg_2$outcome<-factor(fg_2$outcome,levels=c("Luminal B/HER2-negative-like breast cancer","Triple negative or basal-like breast cancer","Luminal B-like breast cancer","Luminal A-like breast cancer", "HER2-enriched-like breast cancer",
                                           "Endometrioid ovarian cancer","High grade serous carcinoma ovarian cancer", "Clear cell ovarian cancer","Low malignant potential ovarian cancer","Invasive mucinous ovarian cancer","Low grade serous carcinoma ovarian cancer",
                                           "Endometrioid endometrial cancer","Non-endometrioid endometrial cancer",
                                           "Colon cancer","Distal colon cancer","Proximal colon cancer","Rectal cancer"))


p2<- ggforestplot::forestplot(
  df = fg_2,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=7) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))

p2


ggplot2::ggsave(filename="Figure 2 wide just subtypes.png", plot=ggplot2::last_plot(),width = 0.7*1000, height = 0.7*480, units = "mm",bg="white",dpi=1000)



# Figure 3 ----------------------------------------------------------------
fg_3<-fread("LiverFatMR/Adiposity_measures_molecular_traits_results_new_proxies.csv")
fg_3<-fg_3[!fg_3$method=="MR Egger",]
fg_3$exposure[fg_3$exposure=="Percent liver fat || id:ebi-a-GCST90016673"]<-"Liver fat"
fg_3$exposure[fg_3$exposure=="Pancreas fat || id:ebi-a-GCST90016675"]<-"Pancreas fat"
fg_3$outcome[fg_3$outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870"]<-"Female SHBG"
fg_3$outcome[fg_3$outcome=="Sex hormone binding globulin (SHBG) || id:ieu-b-4871"]<-"Male SHBG"
fg_3$outcome[fg_3$outcome=="Bioavailable Testosterone || id:ieu-b-4869"]<-"Female bioavailable testosterone"
fg_3$outcome[fg_3$outcome=="Bioavailable Testosterone || id:ieu-b-4868"]<-"Male bioavailable testosterone"
fg_3$outcome[fg_3$outcome=="Total Testosterone || id:ieu-b-4864"]<-"Female total testosterone"
fg_3$outcome[fg_3$outcome=="Total Testosterone || id:ieu-b-4865"]<-"Male total testosterone"
fg_3$exposure[fg_3$exposure=="body mass index || id:ieu-b-40"]<-"BMI"
fg_3$exposure[fg_3$exposure=="Body mass index || id:ieu-a-974"]<-"BMI"
fg_3$exposure[fg_3$exposure=="Body mass index || id:ieu-a-785"]<-"BMI"
fg_3$exposure<-gsub("Female ","",fg_3$exposure)
fg_3$exposure<-gsub("Male ","",fg_3$exposure)
fg_3$exposure<-factor(fg_3$exposure,levels=c("ASAT","VAT","GFAT","Liver fat","Pancreas fat","BMI"))
fg_3$BP[fg_3$pval<(0.05/12)]<-0
fg_3$BP[fg_3$pval>=(0.05/12)]<-1
fg_3<-fg_3 %>%
  group_by(exposure) %>%
  arrange(pval, .by_group = TRUE)
fg_3$method<-factor(fg_3$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","MR Egger", "Wald ratio")))

fg_3$outcome<-gsub("\\|.*","",fg_3$outcome)
fg_3$outcome[fg_3$outcome=="triglycerides "]<-"Triglycerides"
fg_3$outcome[fg_3$outcome=="CXCL8"]<-"CXCL-8"
fg_3$outcome[fg_3$outcome=="MCP1"]<-"MCP-1"
fg_3$outcome[fg_3$outcome=="PAI1"]<-"PAI-1"
fg_3$outcome[fg_3$outcome=="IGFBP1"]<-"IGFBP-1"
fg_3$outcome[fg_3$outcome=="IGFBP3"]<-"IGFBP-3"
fg_3$outcome[fg_3$outcome=="IGF1"]<-"IGF-1"
fg_3$outcome[fg_3$outcome=="IGF2"]<-"IGF-2"
fg_3$outcome<-gsub("TNF-a","TNF-α",fg_3$outcome)
fg_3$outcome<-gsub("IFN-a","IFN-α",fg_3$outcome)
fg_3$outcome<-gsub("IFN-b","IFN-β",fg_3$outcome)
fg_3$outcome<-gsub("IL-1B","IL-1β",fg_3$outcome)
fg_3$outcome<-gsub(" levels","",fg_3$outcome)

fg_3<-fg_3[!(fg_3$method=="Weighted median" & fg_3$nsnp<10),]
fg_3<-fg_3[!(fg_3$method=="Weighted mode" & fg_3$nsnp<10),]
fg_3<-ungroup(fg_3)

fg_3$outcome<-factor(fg_3$outcome,levels=c("Female SHBG",                
                                           "Female bioavailable testosterone",
                                           "Female total testosterone",       
                                           "Male SHBG",                       
                                           "Male total testosterone",         
                                           "Male bioavailable testosterone",  
                                           "Fasting insulin","IGF-1","IGF-2","IGFBP-1","IGFBP-3",
                                           "MCP-1","CXCL-8",
                                           "IL-1β","IL-6","TNF-α","CRP","IFN-α","IFN-β","PAI-1",
                                           "Triglycerides","FASN","HDL cholesterol ",
                                           "Leptin","Visfatin" ,"Resistin","Adiponectin" ))

fg_3$method<-factor(fg_3$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","Wald ratio")))


fg_3<-fg_3[!is.na(fg_3$b),]

p1<- ggforestplot::forestplot(
  df = fg_3,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Beta (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=6) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "bottom", text=element_text(size=10))

p1

ggplot2::ggsave(filename="Figure 3.png", plot=ggplot2::last_plot(),width = 0.75*0.9*480, height = 0.75*0.5*520, units = "mm",bg="white",dpi=1000)

#Making wider for presentations
p1<- ggforestplot::forestplot(
  df = fg_3,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Beta (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=8) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))

p1

ggplot2::ggsave(filename="Figure 3 wide.png", plot=ggplot2::last_plot(),width = 0.58*850, height = 0.58*0.8*520, units = "mm",bg="white",dpi=1000)


# Figure 4 ----------------------------------------------------------------

fg_4<-fread("LiverFatMR/Molecular_traits_cancers_results.csv")

fg_4$method<-factor(fg_4$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","MR Egger","Wald ratio")))

fg_4$exposure[fg_4$exposure=="Sex hormone binding globulin (SHBG) || id:ieu-b-4870"]<-"SHBG (female)"
fg_4$exposure[fg_4$exposure=="Sex hormone binding globulin (SHBG) || id:ieu-b-4871"]<-"SHBG (male)"
fg_4$exposure[fg_4$exposure=="triglycerides || id:ieu-b-111"]<-"Triglycerides"
fg_4$exposure[fg_4$exposure=="HDL cholesterol || id:ieu-b-109"]<-"HDL cholesterol"
fg_4$exposure[fg_4$exposure=="Total Testosterone || id:ieu-b-4864"]<-"Total testosterone (female)"
fg_4$exposure[fg_4$exposure=="Total Testosterone || id:ieu-b-4865"]<-"Total testosterone (male)"
fg_4$exposure[fg_4$exposure=="Bioavailable Testosterone || id:ieu-b-4868"]<-"Bioavailable testosterone (male)"
fg_4$exposure[fg_4$id.exposure=="ieu-b-4864"]<-"Total testosterone (female)"
fg_4$exposure[fg_4$id.exposure=="ieu-b-4868"]<-"Bioavailable testosterone (male)"
fg_4$exposure[fg_4$exposure=="Bioavailable Testosterone || id:ieu-b-4869"]<-"Bioavailable testosterone (female)"
fg_4$exposure[fg_4$id.exposure=="ieu-b-4869"]<-"Bioavailable testosterone (female)"
fg_4$exposure[fg_4$exposure=="IGFBP1"]<-"IGFBP-1"
fg_4$exposure[fg_4$exposure=="CXCL8"]<-"CXCL-8"
fg_4$exposure[fg_4$exposure=="PAI1"]<-"PAI-1"
fg_4$exposure[fg_4$exposure=="Adiponectin"]<-"              Adiponectin"
fg_4$outcome<-gsub("uminal B-HER2","uminal B/HER2",fg_4$outcome)

fg_4_all<-fg_4

fg_4$exposure<-factor(fg_4$exposure,levels=c("Total testosterone (female)","Bioavailable testosterone (female)","SHBG (female)","Total testosterone (male)","Bioavailable testosterone (male)","SHBG (male)","Fasting insulin","IGF-1","IGFBP-1","Resistin","              Adiponectin","PAI-1","CXCL-8","HDL cholesterol","Triglycerides","FASN"))
fg_4<-fg_4[order(fg_4$exposure),]

fg_4<-fg_4[!fg_4$method=="MR Egger",]
fg_4<-fg_4[!(fg_4$method=="Weighted median" & fg_4$nsnp<10),]
fg_4<-fg_4[!(fg_4$method=="Weighted mode" & fg_4$nsnp<10),]
fg_4<-ungroup(fg_4)

fg_4$outcome<-gsub("Esophagus adenocarcinoma","Oesophageal adenocarcinoma",fg_4$outcome)

asat_list<-c("SHBG (female)","SHBG (male)","Fasting insulin","IGF-1","IGFBP-1","CXCL-8","IL-6","TNF-a","IFN-b","PAI-1","Triglycerides","HDL cholesterol","Resistin","              Adiponectin","Bioavailable testosterone (female)","Bioavailable testosterone (male)")
gfat_list<-c("SHBG (male)","IGF-1","IGFBP-1","HDL cholesterol","              Adiponectin")
liver_list<-c("MCP-1","CXCL-8","IL-1B","IL-6","TNF-a","IFN-a","IFN-b","PAI-1","FASN","              Adiponectin")
pancreas_list<-c("Triglycerides","Total testosterone (female)")

asat_list_cancers<-c("Endometrial cancer","Oesophageal adenocarcinoma","Liver cancer","Ovarian cancer","Luminal B/HER2-negative-like breast cancer","Triple negative or basal-like breast cancer","Luminal B-like breast cancer","Endometrioid ovarian cancer","Endometrioid endometrial cancer")
gfat_list_cancers<-c("Breast cancer","Meningioma","Triple negative or basal-like breast cancer","Invasive mucinous ovarian cancer","Non-endometrioid endometrial cancer")
liver_list_cancers<-c("Liver cancer")
pancreas_list_cancers<-c("Multiple myeloma","Endometrioid ovarian cancer","Proximal colon cancer")

fg_4a<-fg_4[fg_4$exposure %in% asat_list,]
fg_4a<-fg_4a[fg_4a$outcome %in% asat_list_cancers,]

fg_4a$exposure<-factor(fg_4a$exposure,levels=asat_list)
fg_4a$outcome<-factor(fg_4a$outcome,levels=asat_list_cancers)

point <- scales::format_format(scientific = FALSE,accuracy=0.1)

p1<- ggforestplot::forestplot(
  df = fg_4a,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=3) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  scale_x_continuous(trans = scales::log_trans(),breaks=c(0.05,0.1,0.3,1,3,10),labels=point)+
  theme(legend.position = "none", text=element_text(size=17)) +
  coord_cartesian(xlim=c(0.05,15)) 

p1


subset_data <- subset(fg_4a, outcome == "Liver cancer" & exposure == "Fasting insulin" & method == "Weighted mode")

p1 <- p1 +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 9, yend = 9 ),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#4472ED",
    data = subset_data
  ) 

p1 <- p1 +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 8.85, yend = 8.85 ),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#FFB800",
    data = subset_data
  ) 

subset_data <- subset(fg_4a, outcome == "Oesophageal adenocarcinoma" & exposure == "Fasting insulin" & method == "Weighted mode")

p1 <- p1 +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 9, yend = 9 ),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#4472ED",
    data = subset_data
  ) 


p1

p1 <- p1+ theme(plot.margin = unit(c(0, 0, 0, 1.4), "cm"))

p1


fg_4c<-fg_4[fg_4$exposure %in% gfat_list,]
fg_4c<-fg_4c[fg_4c$outcome %in% gfat_list_cancers,]

fg_4c$exposure<-factor(fg_4c$exposure,levels=gfat_list)
fg_4c$outcome<-factor(fg_4c$outcome,levels=gfat_list_cancers)


p3<- ggforestplot::forestplot(
  df = fg_4c,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=3) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "none", text=element_text(size=17))

p3

fg_4d<-fg_4[fg_4$exposure %in% liver_list,]
fg_4d<-fg_4d[fg_4d$outcome %in% liver_list_cancers,]

fg_4d$exposure<-factor(fg_4d$exposure,levels=liver_list)
fg_4d$outcome<-factor(fg_4d$outcome,levels=liver_list_cancers)



p4<- ggforestplot::forestplot(
  df = fg_4d,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=3,drop=FALSE) +
  scale_color_manual(values=rev(c("#D11D22","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "none", text=element_text(size=17))


p4

p4 <- p4 + theme(plot.margin = unit(c(0, 25, 0, 1.2), "cm"))



fg_4e<-fg_4[fg_4$exposure %in% pancreas_list,]
fg_4e<-fg_4e[fg_4e$outcome %in% pancreas_list_cancers,]

fg_4e$exposure<-factor(fg_4e$exposure,levels=pancreas_list)
fg_4e$outcome<-factor(fg_4e$outcome,levels=pancreas_list_cancers)


p5<- ggforestplot::forestplot(
  df = fg_4e,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=3) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "bottom", text=element_text(size=17))


p5

library(ggpubr)
labels=c("ASAT","GFAT","Liver fat","Pancreas fat")

fg_4<-ggarrange(p1,p3,p4,p5,ncol=1,labels=c("A","B","C","D"),heights=c(36,15,5.4,6.5),widths=c(1,1,1,1))
fg_4


ggsave(filename="Figure 4.png", plot=last_plot(),width = 420, height = 3.5*200, units = "mm",bg = "white")

#Making wider for presentations
p1<- ggforestplot::forestplot(
  df = fg_4a,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=6) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  scale_x_continuous(trans = scales::log_trans(),breaks=c(0,0.1,0.3,1,3,10,30,100),labels=point)+
  theme(legend.position = "none", text=element_text(size=17))


p1



p3<- ggforestplot::forestplot(
  df = fg_4c,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=5) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "none", text=element_text(size=17))


p3


p4<- ggforestplot::forestplot(
  df = fg_4d,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "none", text=element_text(size=17))

p4

p5<- ggforestplot::forestplot(
  df = fg_4e,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800")),labels=rev(c("Inverse variance weighted or Wald ratio","Weighted mode","Weighted median")))+ guides(colour=guide_legend(title="Method",reverse=TRUE))+
  theme(legend.position = "bottom", text=element_text(size=17))

p5

library(ggpubr)
#labels=c("ASAT","GFAT","Liver fat","Pancreas fat")
fg_4<-ggarrange(p1,p3,p4,p5,ncol=1,labels=c("A","B","C","D"),heights=c(19,8,5,4))
fg_4

ggsave(filename="Figure 4 wide.png", plot=last_plot(),width = 0.4*1200, height = 0.4*3*200, units = "mm",bg = "white")


# Figure 5 ----------------------------------------------------------------

dat<-data.frame(matrix(ncol=7,nrow=0))

#ASAT EC
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=3))
df[,1]<-"ASAT"
df[,2]<-"Endometrial cancer"
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#Fasting insulin
MVMR <- read.csv("results/LiverFatMR/ASAT_FI_EC_MVMR_results.csv")
ASAT_FI_EC <- MVMR[1,2]
ASAT_FI_EC_se <- MVMR[1,3]
df[1,3]<-"Fasting insulin"
df[1,6]<-ASAT_FI_EC
df[1,7]<-ASAT_FI_EC_se

#SHBG
MVMR <- read.csv("results/LiverFatMR/ASAT_SHBG_EC_MVMR_results.csv")
ASAT_SHBG_EC <- MVMR[1,2]
ASAT_SHBG_EC_se <- MVMR[1,3]
df[2,3]<-"SHBG (female)"
df[2,6]<-ASAT_SHBG_EC
df[2,7]<-ASAT_SHBG_EC_se

#Bioavailable testosterone
MVMR <- read.csv("results/LiverFatMR/ASAT_BT_EC_MVMR_results.csv")
ASAT_BT_EC <- MVMR[1,2]
ASAT_BT_EC_se <- MVMR[1,3]
df[3,3]<-"Bioavailable testosterone (female)"
df[3,6]<-ASAT_BT_EC
df[3,7]<-ASAT_BT_EC_se

dat<-rbind(dat,df)

#ASAT EEC
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Endometrioid endometrial cancer",]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=3))
df[,1]<-"ASAT"
df[,2]<-"Endometrioid endometrial cancer"
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#Fasting insulin
MVMR <- read.csv("results/LiverFatMR/ASAT_FI_EEC_MVMR_results.csv")
ASAT_FI_EEC <- MVMR[1,2]
ASAT_FI_EEC_se <- MVMR[1,3]
df[1,3]<-"Fasting insulin"
df[1,6]<-ASAT_FI_EEC
df[1,7]<-ASAT_FI_EEC_se

#SHBG
MVMR <- read.csv("results/LiverFatMR/ASAT_SHBG_EEC_MVMR_results.csv")
ASAT_SHBG_EEC <- MVMR[1,2]
ASAT_SHBG_EEC_se <- MVMR[1,3]
df[2,3]<-"SHBG (female)"
df[2,6]<-ASAT_SHBG_EEC
df[2,7]<-ASAT_SHBG_EEC_se

#Bioavailable testosterone
MVMR <- read.csv("results/LiverFatMR/ASAT_BT_EEC_MVMR_results.csv")
ASAT_BT_EEC <- MVMR[1,2]
ASAT_BT_EEC_se <- MVMR[1,3]
df[3,3]<-"Bioavailable testosterone (female)"
df[3,6]<-ASAT_BT_EEC
df[3,7]<-ASAT_BT_EEC_se

dat<-rbind(dat,df)

#ASAT-IGFBP-1-oesophageal adenocarcinoma
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Esophagus adenocarcinoma" ,]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=1))
df[,1]<-"ASAT"
df[,2]<-"Oesophageal adenocarcinoma" 
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#IGFBP-1
MVMR <- read.csv("results/LiverFatMR/ASAT_IGFBP_OAC_MVMR_results.csv")
ASAT_IGFBP_OAC <- MVMR[1,2]
ASAT_IGFBP_OAC_se <- MVMR[1,3]
df[1,3]<-"IGFBP-1"
df[1,6]<-ASAT_IGFBP_OAC
df[1,7]<-ASAT_IGFBP_OAC_se

dat<-rbind(dat,df)

#ASAT-CXCL-8-Ovarian cancer
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Ovarian cancer" ,]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=1))
df[,1]<-"ASAT"
df[,2]<-"Ovarian cancer" 
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#CXCL-8
MVMR <- read.csv("results/LiverFatMR/ASAT_CXC_OC_MVMR_results.csv")
ASAT_CXCL_OAC <- MVMR[1,2]
ASAT_CXC_OC_se <- MVMR[1,3]
df[1,3]<-"CXCL-8"
df[1,6]<-ASAT_CXCL_OAC
df[1,7]<-ASAT_CXC_OC_se

dat<-rbind(dat,df)

#ASAT-HDl cholesterol-Triple negative or basal-like breast cancer
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="ASAT",]
res<-res[res$outcome=="Triple negative or basal-like breast cancer" ,]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=1))
df[,1]<-"ASAT"
df[,2]<-"Triple negative or basal-like breast cancer" 
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#HDl cholesterol
MVMR <- read.csv("results/LiverFatMR/ASAT_HDL_TNBC_MVMR_results.csv")
ASAT_HDL_TNBC <- MVMR[1,2]
ASAT_HDL_TNBC_se <- MVMR[1,3]
df[1,3]<-"HDl cholesterol"
df[1,6]<-ASAT_HDL_TNBC
df[1,7]<-ASAT_HDL_TNBC_se

dat<-rbind(dat,df)

#GFAT
#GFAT-adiponectin-Non-endometrioid endometrial cancer
res<-fread("results/LiverFatMR/Adiposity_measures_obesity_cancers_results_new_proxies.csv")
res<-res[res$exposure=="GFAT",]
res<-res[res$outcome=="Non-endometrioid endometrial cancer" ,]
res<-res[(res$method=="Inverse variance weighted" | res$method=="Wald ratio"),]
unadjusted<-res[1,7]
unadjusted_se<-res[1,8]

df<-data.frame(matrix(ncol=7,nrow=1))
df[,1]<-"GFAT"
df[,2]<-"Non-endometrioid endometrial cancer" 
df[,4]<-unadjusted
df[,5]<-unadjusted_se


#HDl cholesterol
MVMR <- read.csv("results/LiverFatMR/GFAT_AD_NEC_MVMR_results.csv")
GFAT_AD_NEC <- MVMR[1,2]
GFAT_AD_NEC_se <- MVMR[1,3]
df[1,3]<-"Adiponectin"
df[1,6]<-GFAT_AD_NEC
df[1,7]<-GFAT_AD_NEC_se

dat<-rbind(dat,df)

colnames(dat)<-c("Exposure","Outcome","Mediator","Unadjusted","Unadjusted SE","Adjusted","Adjusted SE")

fwrite(dat,"results/LiverFatMR/Adiposity_cancers_mediation_results.csv")

plot<-data.frame(matrix(ncol=5,nrow=16))
colnames(plot)<-c("Plot","Exposure/Mediator","Beta","SE","Outcome")
plot$`Exposure/Mediator`<-c("ASAT",dat$Mediator[1:3],"ASAT",dat$Mediator[4:6],"ASAT",dat$Mediator[7],"ASAT",dat$Mediator[8],"ASAT",dat$Mediator[9],"GFAT",dat$Mediator[10])
plot$Plot<-c(rep("A",14),rep("B",2))
plot$Outcome<-c(rep("Endometrial cancer",4),rep("Endometrioid endometrial cancer",4),rep("Oesophageal adenocarcinoma",2),rep("Ovarian cancer",2),rep("Triple negative or basal-like breast cancer",2),rep("Non-endometrioid endometrial cancer",2))
plot$Beta[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Endometrial cancer"]<-dat$Unadjusted[1]
plot$Beta[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$Unadjusted[4]
plot$Beta[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Oesophageal adenocarcinoma"]<-dat$Unadjusted[7]
plot$Beta[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Ovarian cancer"]<-dat$Unadjusted[8]
plot$Beta[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Triple negative or basal-like breast cancer"]<-dat$Unadjusted[9]
plot$Beta[plot$`Exposure/Mediator`=="GFAT" & plot$Outcome=="Non-endometrioid endometrial cancer"]<-dat$Unadjusted[10]
plot$SE[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Endometrial cancer"]<-dat$`Unadjusted SE`[1:3]
plot$SE[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$`Unadjusted SE`[4:6]
plot$SE[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Oesophageal adenocarcinoma"]<-dat$`Unadjusted SE`[7]
plot$SE[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Ovarian cancer"]<-dat$`Unadjusted SE`[8]
plot$SE[plot$`Exposure/Mediator`=="ASAT" & plot$Outcome=="Triple negative or basal-like breast cancer"]<-dat$`Unadjusted SE`[9]
plot$SE[plot$`Exposure/Mediator`=="GFAT" & plot$Outcome=="Non-endometrioid endometrial cancer"]<-dat$`Unadjusted SE`[10]

plot$Beta[plot$`Exposure/Mediator`=="Fasting insulin" & plot$Outcome=="Endometrial cancer"]<-dat$Adjusted[dat$Mediator=="Fasting insulin" & dat$Outcome=="Endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="SHBG (female)" & plot$Outcome=="Endometrial cancer"]<-dat$Adjusted[dat$Mediator=="SHBG (female)" & dat$Outcome=="Endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="Bioavailable testosterone (female)" & plot$Outcome=="Endometrial cancer"]<-dat$Adjusted[dat$Mediator=="Bioavailable testosterone (female)" & dat$Outcome=="Endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="Fasting insulin" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$Adjusted[dat$Mediator=="Fasting insulin" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="SHBG (female)" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$Adjusted[dat$Mediator=="SHBG (female)" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="Bioavailable testosterone (female)" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$Adjusted[dat$Mediator=="Bioavailable testosterone (female)" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$Beta[plot$`Exposure/Mediator`=="IGFBP-1"]<-dat$Adjusted[dat$Mediator=="IGFBP-1"]
plot$Beta[plot$`Exposure/Mediator`=="CXCL-8"]<-dat$Adjusted[dat$Mediator=="CXCL-8"]
plot$Beta[plot$`Exposure/Mediator`=="HDl cholesterol"]<-dat$Adjusted[dat$Mediator=="HDl cholesterol"]
plot$`Exposure/Mediator`[plot$`Exposure/Mediator`=="HDl cholesterol"]<-"HDL cholesterol"
plot$Beta[plot$`Exposure/Mediator`=="Adiponectin"]<-dat$Adjusted[dat$Mediator=="Adiponectin"]

plot$SE[plot$`Exposure/Mediator`=="Fasting insulin" & plot$Outcome=="Endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="Fasting insulin" & dat$Outcome=="Endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="SHBG (female)" & plot$Outcome=="Endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="SHBG (female)" & dat$Outcome=="Endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="Bioavailable testosterone (female)" & plot$Outcome=="Endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="Bioavailable testosterone (female)" & dat$Outcome=="Endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="Fasting insulin" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="Fasting insulin" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="SHBG (female)" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="SHBG (female)" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="Bioavailable testosterone (female)" & plot$Outcome=="Endometrioid endometrial cancer"]<-dat$`Adjusted SE`[dat$Mediator=="Bioavailable testosterone (female)" & dat$Outcome=="Endometrioid endometrial cancer"]
plot$SE[plot$`Exposure/Mediator`=="IGFBP-1"]<-dat$`Adjusted SE`[dat$Mediator=="IGFBP-1"]
plot$SE[plot$`Exposure/Mediator`=="CXCL-8"]<-dat$`Adjusted SE`[dat$Mediator=="CXCL-8"]
plot$SE[plot$`Exposure/Mediator`=="HDL cholesterol"]<-dat$`Adjusted SE`[dat$Mediator=="HDl cholesterol"]
plot$SE[plot$`Exposure/Mediator`=="Adiponectin"]<-dat$`Adjusted SE`[dat$Mediator=="Adiponectin"]

plot$`Exposure/Mediator`<-paste("+ ",plot$`Exposure/Mediator`,sep="")
plot$`Exposure/Mediator`[plot$`Exposure/Mediator`=="+ ASAT"]<-"ASAT (Unadjusted)"
plot$`Exposure/Mediator`[plot$`Exposure/Mediator`=="+ GFAT"]<-"GFAT (Unadjusted)"

plot$LCI<-plot$Beta-1.96*plot$SE
plot$UCI<-plot$Beta+1.96*plot$SE
plot$pval<-1
plot$pval[plot$Beta>0 & plot$LCI>0]<-0.04
plot$pval[plot$Beta<0 & plot$UCI<0]<-0.04

p1<- ggforestplot::forestplot(
  df = plot,
  pvalue = pval,
  name=`Exposure/Mediator`,
  estimate = Beta,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=SE,
  colour=Plot
)+
  ggplot2::facet_wrap(vars(Outcome), scales = "free_y") +
  scale_color_manual(values=rev(c("#D11D22","#D11D22")))+theme(legend.position = "none")

p1

setwd("~/OneDrive - University of Bristol/Documents/Year 3/Liver Fat MR/Paper/")

ggsave(filename="Figure 5.png", plot=last_plot(),width = 420, height = 0.7*300, units = "mm",bg = "white")

# Supplementary figures ----
### Supplementary figure 3 ####
st_3<-fread("LiverFatMR/Repeating_no_sample_overlap.csv")
st_3$outcome<-gsub("uminal B-HER2","uminal B/HER2",st_3$outcome)
st_3$outcome<-gsub(" no UK Biobank"," (no UK Biobank)",st_3$outcome)

st_3<-st_3[!(st_3$method=="Weighted median" & st_3$nsnp<10),]
st_3<-st_3[!(st_3$method=="Weighted mode" & st_3$nsnp<10),]

st_3$analysis<-NA

st_3$analysis[st_3$exposure=="VAT" & st_3$outcome == "Thyroid cancer (FinnGen)"]<-"Adiposity-cancer"
st_3$analysis[st_3$outcome == "Liver cancer (FinnGen)"]<-"Adiposity-cancer"
st_3$analysis[st_3$outcome == "HDL cholesterol (no UK Biobank)"]<-"Adiposity-molecular"
st_3$analysis[st_3$outcome == "IGF-1 (no UK Biobank)"]<-"Adiposity-molecular"
st_3$analysis[st_3$outcome == "Triglycerides (no UK Biobank)"]<-"Adiposity-molecular"


st_3$exposure[st_3$exposure=="Percent liver fat || id:ebi-a-GCST90016673"]<-"Liver fat"
st_3$exposure[st_3$exposure=="Pancreas fat || id:ebi-a-GCST90016675"]<-"Pancreas fat"
st_3$exposure[st_3$exposure=="body mass index || id:ieu-b-40"]<-"BMI"

st_3$analysis[st_3$exposure=="Pancreas fat" & st_3$outcome == "Multiple myeloma (FinnGen)"]<-"Adiposity-cancer"

st_3$exposure<-factor(st_3$exposure,levels=c("ASAT","VAT","GFAT","Liver fat","Pancreas fat","BMI"))

fg_2_all<-fg_2_all[,1:10]

adiposity_cancer<-st_3[st_3$analysis=="Adiposity-cancer",]
adiposity_cancer<-adiposity_cancer[,1:10]
colnames(adiposity_cancer)[10]<-"BP"

st_3primarya<-fg_2_all[fg_2_all$exposure=="VAT" & fg_2_all$outcome=="Liver cancer",]
st_3primaryb<-fg_2_all[fg_2_all$exposure=="ASAT" & fg_2_all$outcome=="Liver cancer",]
st_3primaryc<-fg_2_all[fg_2_all$exposure=="Liver fat" & fg_2_all$outcome=="Liver cancer",]
st_3primaryd<-fg_2_all[fg_2_all$exposure=="VAT" & fg_2_all$outcome=="Thyroid cancer",]
st_3primarye<-fg_2_all[fg_2_all$exposure=="Pancreas fat" & fg_2_all$outcome=="Multiple myeloma",]

adiposity_cancer<-rbind(adiposity_cancer,st_3primarya,st_3primaryb,st_3primaryc,st_3primaryd,st_3primarye)

p1a<- ggforestplot::forestplot(
  df = adiposity_cancer,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=6) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#D11D22")))+ guides(colour=guide_legend(title="Method",reverse=TRUE)) +
  theme(legend.position = "none", text=element_text(size=17))+
  coord_cartesian(xlim=c(0.05,15))

p1a

subset_data <- subset(adiposity_cancer, outcome == "Liver cancer (FinnGen)" & exposure == "VAT")

p1a <- p1a +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 4, yend = 4),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#D11D22",
    data = subset_data
  ) 

p1a

subset_data <- subset(adiposity_cancer, outcome == "Thyroid cancer (FinnGen)" & exposure == "VAT")

p1a <- p1a +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 4, yend = 4),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#D11D22",
    data = subset_data
  ) 

p1a


p1a <- p1a+ theme(plot.margin = unit(c(0, 0, 0, 1.35), "cm"))



adiposity_molecular<-st_3[st_3$analysis=="Adiposity-molecular",]
adiposity_molecular<-adiposity_molecular[,1:10]
colnames(adiposity_molecular)[10]<-"BP"

fg_3<-fg_3[,1:10]


st_3primaryc<-fg_3[fg_3$exposure=="ASAT" & fg_3$outcome=="HDL cholesterol ",]
st_3primaryd<-fg_3[fg_3$exposure=="GFAT" & fg_3$outcome=="HDL cholesterol ",]
st_3primarye<-fg_3[fg_3$exposure=="ASAT" & fg_3$outcome=="IGF-1",]
st_3primaryf<-fg_3[fg_3$exposure=="GFAT" & fg_3$outcome=="IGF-1",]
st_3primaryg<-fg_3[fg_3$exposure=="ASAT" & fg_3$outcome=="Triglycerides",]
st_3primaryh<-fg_3[fg_3$exposure=="GFAT" & fg_3$outcome=="Triglycerides",]
st_3primaryi<-fg_3[fg_3$exposure=="Pancreas fat" & fg_3$outcome=="Triglycerides",]

adiposity_molecular<-rbind(adiposity_molecular,st_3primaryc,st_3primaryd,st_3primarye,st_3primaryf,st_3primaryg,st_3primaryh,st_3primaryi)

adiposity_molecular$outcome<-factor(adiposity_molecular$outcome,levels=c("SHBG ","SHBG (no UK Biobank)","HDL cholesterol ","HDL cholesterol (no UK Biobank)","IGF-1","IGF-1 (no UK Biobank)","Triglycerides","Triglycerides (no UK Biobank)"))

p1b<- ggforestplot::forestplot(
  df = adiposity_molecular,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Beta (95% CI)",
  logodds=FALSE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol = 6) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#D11D22")),labels=c("Inverse variance weighted or Wald ratio","Weighted mode","Weighted median"))+ guides(colour=guide_legend(title="Method",reverse=FALSE)) +
  theme(legend.position = "bottom", text=element_text(size=17))


p1b



p1<-ggarrange(p1a,p1b,ncol=1,labels=c("A","B"),heights=c(3.5,3))
p1
ggplot2::ggsave(filename="Supplementary figure 3.png", plot=ggplot2::last_plot(),width = 0.6*0.9*1000, height = 1.06*480, units = "mm",bg="white",dpi=1000)

### Supplementary figures 4 and 5 ####
st_4<-fread("LiverFatMR/Sex_specific_results.csv")
st_4$outcome<-gsub("\\|.*","",st_4$outcome)
st_4$exposure<-gsub("\\|.*","",st_4$exposure)
st_4$exposure<-gsub("IGF1","IGF-1",st_4$exposure)
st_4$exposure[st_4$exposure=="Body mass index "]<-"BMI"
st_4$exposure<-gsub("Female ","",st_4$exposure)
st_4$outcome<-gsub("uminal B-HER2","uminal B/HER2",st_4$outcome)
st_4<-st_4[!(st_4$method=="Weighted median" & st_4$nsnp<10),]
st_4<-st_4[!(st_4$method=="Weighted mode" & st_4$nsnp<10),]

cancerlist<-c("Endometrial cancer" ,"Endometrioid endometrial cancer","Non-endometrioid endometrial cancer","Ovarian cancer","High grade serous carcinoma ovarian cancer","Low grade serous carcinoma ovarian cancer","Invasive mucinous ovarian cancer","Endometrioid ovarian cancer","Clear cell ovarian cancer","Low malignant potential ovarian cancer","Breast cancer","Luminal A-like breast cancer" ,"Luminal B/HER2-negative-like breast cancer","Luminal B-like breast cancer" ,"HER2-enriched-like breast cancer","Triple negative or basal-like breast cancer" ,"Luminal B-HER2-negative-like breast cancer")

adiposity_cancer<-st_4[st_4$outcome %in% cancerlist,]
adiposity_cancer<-adiposity_cancer[adiposity_cancer$exposure %in% c("BMI","ASAT","VAT","GFAT","Liver fat","Pancreas fat"),]
adiposity_cancer_sexcombined<-fg_2_all[fg_2_all$outcome %in% cancerlist,]
adiposity_cancer_sexcombined<-adiposity_cancer_sexcombined[adiposity_cancer_sexcombined$exposure %in% adiposity_cancer$exposure,]

adiposity_cancer$outcome<-paste(adiposity_cancer$outcome," female-specific",sep="")

adiposity_cancer<-rbind(adiposity_cancer,adiposity_cancer_sexcombined,fill=TRUE)


adiposity_cancer$type<-"Subtype"
adiposity_cancer$type[adiposity_cancer$outcome=="Breast cancer female-specific"|adiposity_cancer$outcome=="Endometrial cancer female-specific"|adiposity_cancer$outcome=="Ovarian cancer female-specific"|adiposity_cancer$outcome=="Breast cancer"|adiposity_cancer$outcome=="Endometrial cancer"|adiposity_cancer$outcome=="Ovarian cancer"]<-"Overall"

adiposity_cancer$exposure<-factor(adiposity_cancer$exposure,levels=c("ASAT","VAT","GFAT","BMI"))
adiposity_cancer$method<-factor(adiposity_cancer$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","Wald ratio")))

adiposity_cancer<-adiposity_cancer[order(adiposity_cancer$exposure),]
adiposity_cancer<-adiposity_cancer[order(adiposity_cancer$method),]

adiposity_cancer_overall<-adiposity_cancer[type == "Overall",]
adiposity_cancer_subtype<-adiposity_cancer[type == "Subtype",]

p1a<- ggforestplot::forestplot(
  df = adiposity_cancer_overall,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE)) +
  theme(legend.position = "none", text=element_text(size=13))+
  coord_cartesian(xlim=c(0.05,15))

p1a


p1b<- ggforestplot::forestplot(
  df = adiposity_cancer_subtype,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE)) +
  theme(legend.position = "bottom", text=element_text(size=13))+
  coord_cartesian(xlim=c(0.05,15))

p1b


subset_data <- subset(adiposity_cancer_subtype, outcome == "Non-endometrioid endometrial cancer female-specific")

p1b <- p1b +
  geom_segment(
    aes(x = 20, xend = 20+0.001, y = 3, yend = 3),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#00CCCC",
    data = subset_data
  ) 

p1b

subset_data <- subset(adiposity_cancer_subtype, outcome == "Clear cell ovarian cancer female-specific")

p1b <- p1b +
  geom_segment(
    aes(x = 0.039, xend = 0.039-0.001, y = 3, yend = 3),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#00CCCC",
    data = subset_data
  ) 

p1b

subset_data <- subset(adiposity_cancer_subtype, outcome == "Low grade serous carcinoma ovarian cancer female-specific")

p1b <- p1b +
  geom_segment(
    aes(x = 0.039, xend = 0.039-0.001, y = 3, yend = 3),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#00CCCC",
    data = subset_data
  ) 

p1b

p1<-ggarrange(p1a,p1b,ncol=1,labels=c("A","B"),heights=c(1.7,5))
p1
ggplot2::ggsave(filename="Supplementary figure 4.png", plot=ggplot2::last_plot(),width = 500, height = 520, units = "mm",bg="white",dpi=1000)

molecular_cancer<-st_4[st_4$outcome %in% cancerlist,]
molecular_cancer<-molecular_cancer[!molecular_cancer$exposure %in% c("BMI","ASAT","VAT","GFAT","Liver fat","Pancreas fat"),]
molecular_cancer_sexcombined<-fg_4_all[fg_4_all$outcome %in% cancerlist,]
molecular_cancer_sexcombined<-molecular_cancer_sexcombined[molecular_cancer_sexcombined$exposure %in% molecular_cancer$exposure,]

molecular_cancer$outcome<-paste(molecular_cancer$outcome," female-specific",sep="")

molecular_cancer<-rbind(molecular_cancer,molecular_cancer_sexcombined,fill=TRUE)


molecular_cancer$type<-"Subtype"
molecular_cancer$type[molecular_cancer$outcome=="Breast cancer female-specific"|molecular_cancer$outcome=="Endometrial cancer female-specific"|molecular_cancer$outcome=="Ovarian cancer female-specific"|molecular_cancer$outcome=="Breast cancer"|molecular_cancer$outcome=="Endometrial cancer"|molecular_cancer$outcome=="Ovarian cancer"]<-"Overall"

molecular_cancer$method<-factor(molecular_cancer$method,levels=rev(c("Inverse variance weighted","Weighted mode","Weighted median","Wald ratio")))

molecular_cancer<-molecular_cancer[order(molecular_cancer$method),]

molecular_cancer_overall<-molecular_cancer[type == "Overall",]
molecular_cancer_subtype<-molecular_cancer[type == "Subtype",]

p1a<- ggforestplot::forestplot(
  df = molecular_cancer_overall,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE)) +
  theme(legend.position = "none", text=element_text(size=17))+
  coord_cartesian(xlim=c(0.05,15))

p1a

p1b<- ggforestplot::forestplot(
  df = molecular_cancer_subtype,
  name=exposure,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=se,
  colour = method
)+
  ggplot2::facet_wrap(vars(outcome),ncol=4) +
  scale_color_manual(values=rev(c("#D11D22","#4472ED","#FFB800","#00CCCC")))+ guides(colour=guide_legend(title="Method",reverse=TRUE)) +
  theme(legend.position = "bottom", text=element_text(size=17))+
  coord_cartesian(xlim=c(0.05,15))

p1b

subset_data <- subset(molecular_cancer_subtype, outcome == "Non-endometrioid endometrial cancer female-specific")

p1b <- p1b +
  geom_segment(
    aes(x = 0.039, xend = 0.039-0.001, y = 7, yend = 7),
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "#D11D22",
    data = subset_data
  ) 

p1b

p1<-ggarrange(p1a,p1b,ncol=1,labels=c("A","B"),heights=c(2.5,5))
p1
ggplot2::ggsave(filename="Supplementary figure 5.png", plot=ggplot2::last_plot(),width = 1.37*480, height = 1.35*520, units = "mm",bg="white",dpi=1000)



# Supplementary figure 6 --------------------------------------------------

st_6<-read.xlsx("~/OneDrive - University of Bristol/Documents/Year 3/Liver Fat MR/Paper/Tables/Supplementary tables.xlsx",sheet="Supplementary table 7")

plot<-data.frame(matrix(ncol=5,nrow=2))
colnames(plot)<-c("pval","Exposure/Mediator","Beta","SE","Plot")

plot$pval<-c(0.04,1)
plot$`Exposure/Mediator`<-c("ASAT (unadjusted)","+ HDL cholesterol (no UK Biobank)")
plot$Beta<-c(-0.843191848,-0.1992097)
plot$SE<-c(0.257168095,0.09570179)
plot$Plot<-"A"

p1<- ggforestplot::forestplot(
  df = plot,
  pvalue = pval,
  name=`Exposure/Mediator`,
  estimate = Beta,
  psignif = 0.05,
  xlab = "Odds ratio (95% CI)",
  logodds=TRUE,
  se=SE,
  colour=Plot
) +
  scale_color_manual(values=rev(c("#D11D22","#D11D22")))+theme(legend.position = "none")

p1


ggplot2::ggsave(filename="Supplementary figure 6.png", plot=ggplot2::last_plot(),width = 0.35*1.37*480, height = 0.35*1.35*520, units = "mm",bg="white",dpi=1000)
