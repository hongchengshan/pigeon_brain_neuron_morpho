library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggord)

df <- read_excel("pca.xlsx")

str(df)

head(df)

rownames(df) <- df$cellnumber

df$cellnumber <- NULL 

df

dt<-as.matrix(scale(df[,1:5]))
head(dt)

ord <-prcomp(dt[,1:5])
summary(ord)

dt<-ord$x
head(dt)
df<-data.frame(dt,df$region)
head(df)

summ<-summary(ord)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

p <- ggord(ord, df$df.region,coord_fix=F,labcol = 'black', repel = TRUE,
cols = c('#BDDEBA', '#E9BBAF', '#B7E0EA','#FFC0CB','#E98F9E'), size=2,alpha=1,
arrow=0.25, vec_ext =2,poly = FALSE,polylntyp="dashed")
p
ggsave("/data/work/AR_pca/neuron_pca.pdf", width = 6, height = 4.5)