library(ggplot2)
library(reshape2)


data <- read.table(text = "
        ASAT BMI GFAT VAT Liver_fat Pancreas_fat
ASAT    1.0 0.8 0.8 0.7 0.5 0.4
BMI     0.8 1.0 0.7 0.6 0.4 0.4
GFAT    0.8 0.7 1.0 0.5 0.2 0.3
VAT     0.7 0.6 0.5 1.0 0.6 0.6
Liver_fat   0.5 0.4 0.2 0.6 1.0 0.3
Pancreas_fat    0.4 0.4 0.3 0.6 0.3 1.0", header = TRUE)

rownames(data)<-c("ASAT","BMI","GFAT","VAT","Liver fat","Pancreas fat")
colnames(data)<-c("ASAT","BMI","GFAT","VAT","Liver fat","Pancreas fat")

data <- data[order(rownames(data)), order(colnames(data))]

data$row <- rownames(data)

data_melted <- melt(data, id.vars = "row")

colnames(data_melted)<-c("Column","Row","value")

data_melted$Column<-factor(data_melted$Column,levels=c("BMI","ASAT","GFAT","VAT","Liver fat","Pancreas fat"))
data_melted$Row<-factor(data_melted$Row,levels=rev(c("BMI","ASAT","GFAT","VAT","Liver fat","Pancreas fat")))

# Create the heatmap
ggplot(data_melted, aes(x = Column, y = Row, fill = value, label = round(value, 2))) +
  geom_tile() +
  geom_text(color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "red",limits=c(0,1)) +
  labs(title = " ", x = " ", y = " ",fill="RG") +
  theme_minimal()

ggsave("results/LiverFatMR/Genetic Correlations Heatmap.png",plot=last_plot(),dpi=600,width=7.5,height=6,bg = "white")

