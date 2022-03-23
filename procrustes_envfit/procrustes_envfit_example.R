library(readxl)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(dplyr)
library(tidyr)


##read in 16S feature table with features as columns and samples as rows
data = read.delim("big_soil_feature_table_CTprocrustes.txt", sep = '\t', header = TRUE)
rownames(data)=data[,1]
data = data[,-1]

##read in chemistry feature table with features as columns and samples as rows
chem = read.delim("HILIC_data_CTprocrustes.txt", sep = '\t', header = TRUE)
rownames(chem)=chem[,1]
chem = chem[,-1]

##read in factors
factors = read.delim("HILIC_factors_CTprocrustes.txt", sep = '\t', header = TRUE)
rownames(factors)=factors[,1]
factors = factors[,-1]
factors$day = as.character(factors$day)

# log10 transform both features and metabolites
log_features=log(data+1,10)
log_chem=log(chem+1,10)

##NMDS on ASV features
NMDS_Bray_data <-metaMDS(log_features, distance = "bray", k=2,
                         autotransform = FALSE, noshare = 0.1, trace = 1)
ord.data = as.data.frame(scores(NMDS_Bray_data), display="sites")
NMDS_Bray_chem <-metaMDS(log_chem, distance = "bray", k =2,
                         autotransform = FALSE, noshare = 0.1, trace = 1)
ord.chem = as.data.frame(scores(NMDS_Bray_chem), display="sites")
#plot ASV NMDS
data_plot = ggplot(ord.data, aes(x=NMDS1, y=NMDS2, color = factors$day, shape = factors$replicate)) +
  geom_point() +
  theme_classic()
data_plot
#plot metabolite NMDS
chem_plot = ggplot(ord.chem, aes(x=NMDS1, y=NMDS2, shape = factors$replicate, color = factors$day)) +
  geom_point() +
  theme_classic()
chem_plot


######
## procrustes test and plotting data
## https://john-quensen.com/wp-content/uploads/2018/10/Oksanen-Jari-vegantutor.pdf
######

#run the procrustes test= get the p-value here
pro = protest(ord.chem, ord.data)
pro
#plot the procrustes results
proc = procrustes(ord.chem, ord.data, symmetirc = TRUE)
summary(proc)
pro_df = data.frame(nmds1=proc$Yrot[,1], 
                    nmds2=proc$Yrot[,2], 
                    xnmds1=proc$X[,1],
                    xnmds2=proc$X[,2], 
                    day=c("0","1","3","5","7","10","14","20",
                          "0","1","3","5","7","10","14","20",
                          "0","1","3","5","7","10","14","20"),
                    trip=rep(c("A","B","C"), each = 8))
pro_plot = ggplot(pro_df, aes(color = day)) +
  geom_point(aes(x=nmds1, y=nmds2), shape = 1, size =5) +
  geom_point(aes(x=xnmds1, y=xnmds2), shape = 0, size =5) +
  geom_segment(aes(x=nmds1, y=nmds2, xend=xnmds1, yend=xnmds2), arrow = arrow(length = unit(0.2,"cm"))) +
  theme_classic()
pro_plot


######
## envfit to overlay significant loadings on the ordination
######
fit=envfit(ord.chem, chem)
fit.df<-as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))
fit.df$species<-rownames(fit.df)
fit.df$pvals=fit$vectors$pvals
#filter to only be significant loadings
fit.df_0.05=filter(fit.df,pvals<=0.05)

#note sometimes the vectors are too long, so i usually divide them all by the same value to scale them all to fit in the plot window
ggplot() +
  geom_point(data=ord.chem, aes(x=NMDS1, y=NMDS2, shape = factors$replicate, color = factors$day)) +
  geom_segment(data=fit.df_0.05,aes(x=0,xend=NMDS1/10,y=0,yend=NMDS2/10), arrow = arrow(length = unit(0.5, "cm")),colour="grey") + 
  geom_text(data=fit.df_0.05,aes(x=NMDS1/10,y=NMDS2/10, label=species),size=3)+
  theme_classic()


