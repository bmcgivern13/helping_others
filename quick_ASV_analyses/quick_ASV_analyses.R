library(readxl)
library(ggplot2)
library(vegan)
library(grid)
library(dplyr)
library(MASS)
library(tidyr)
library(RColorBrewer)
library(stats)
library(readxl)

########
##Read in trimmed data and create factors
########
#read in count table with features as columns and samples as rows
feature_counts = read_excel("example_table.xlsx",sheet="count_table")

#read in taxonomy tab
taxonomy= read_excel("example_table.xlsx",sheet="taxonomy")

#read in metadata
metadata=read_excel("example_table.xlsx",sheet="metadata")

########
##Convert Feature Counts to Featrue Abundance
########
feature_counts=as.data.frame(feature_counts)
row.names(feature_counts)=feature_counts[,1]
feature_counts=feature_counts[,-1]
feature_abundance = feature_counts
feature_abundance$row_sum=rowSums(feature_abundance)
feature_abundance=feature_abundance/feature_abundance$row_sum
# check that rows add up to 1
rowSums(feature_abundance[,-591])
# remove sum column
feature_abundance=feature_abundance[,-591]

#########
## remove ASVs found in two or fewer samples
#########
# find ASVs found in 3+ samples
feature_abundance_ge3samples = as.data.frame(colSums(ifelse(feature_abundance>0,1,0)))
ASV_ge3samples=filter(feature_abundance_ge3samples,feature_abundance_ge3samples[,1]>2)
ASV_ge3samples$asv=row.names(ASV_ge3samples)
feature_abundance_t = as.data.frame(t(feature_abundance))
feature_abundance_t$asv=row.names(feature_abundance_t)
feature_abundance_ge3samples=left_join(ASV_ge3samples,feature_abundance_t,by=c("asv"))
row.names(feature_abundance_ge3samples)=feature_abundance_ge3samples$asv
feature_abundance_ge3samples=feature_abundance_ge3samples[,3:41]
feature_abundance_ge3samples=as.data.frame(t(feature_abundance_ge3samples))

########
## make an NMDS. here, i am directly running on the filtered abundance table. 
## however, you may want to filter out more ASVs or log-transform the matrix.
########
set.seed(1)
# calculate Bray-Curtis distance matrix
bray_dist = metaMDSdist(feature_abundance_ge3samples, distance = "bray", 
                       autotransform = FALSE, noshare = 0.1, trace = 1)
# run NMDS using Bray-Curtis distance matrix
NMDS_Bray <-metaMDS(feature_abundance_ge3samples, distance = "bray", 
                        autotransform = FALSE, noshare = 0.1, trace = 1)
ord.scrs = as.data.frame(scores(NMDS_Bray), display="sites")
ord.scrs$sample = row.names(ord.scrs)
#join with sample metadata
ord.scrs=left_join(ord.scrs,metadata,by=c("sample"))

# define custom color pallette for the three treatments
treats=c("purple","orange","green")

#plot NMDS coloring points by treatment and labelling by day
ggplot(ord.scrs, aes(x=NMDS1, y =NMDS2, color=treatment))+
  geom_point(size = 6,alpha=0.7) +
  geom_text(aes(label = day), alpha= 0.6, color = "black") +
  scale_color_manual(values=treats)+
  theme_classic()

#plot NMDS coloring points by treatment and labelling by day,and adding ellipses by group
ggplot(ord.scrs, aes(x=NMDS1, y =NMDS2, color=treatment))+
  geom_point(size = 6,alpha=0.7) +
  geom_text(aes(label = day), alpha= 0.6, color = "black") +
  scale_color_manual(values=treats)+
  stat_ellipse(level=0.95)+
  theme_classic()

########
##plot stacked barplots
########
feature_abundance_ge3samples$sample=row.names(feature_abundance_ge3samples)
rowSums(feature_abundance_ge3samples[,-115])
tidy_features = feature_abundance_ge3samples %>%
  gather(-sample,key="ASV",value="rel_abun")

tidy_features = left_join(tidy_features,taxonomy,by="ASV")
tidy_features=left_join(tidy_features,metadata,by="sample")

#get colors to match taxonomy level, here using class
#count number of class
tidy_features %>%group_by(Class)%>%summarise(n=n())
#make this equal the number of colors you need
nb.cols=9
c_colors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(tidy_features)+
  geom_col(aes(x=sample,y=rel_abun,fill=Class))+
  scale_fill_manual(values = c_colors) +
  theme_classic()

# re-color by order
tidy_features %>%group_by(Order)%>%summarise(n=n())
nb.cols=11
o_colors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(tidy_features)+
  geom_col(aes(x=sample,y=rel_abun,fill=Order))+
  scale_fill_manual(values = o_colors) +
  theme_classic()

######
##diversity
######
shan = diversity(feature_abundance_ge3samples[,-115], index="shannon")
simp = diversity(feature_abundance_ge3samples[,-115], index="simpson")
invsimp = diversity(feature_abundance_ge3samples[,-115], index="invsimp")
#convert to presence/absence for ASV richness
ASV_rich = rowSums(ifelse(feature_abundance_ge3samples[,-115]>0,1,0))

metrics = metadata
metrics$h = shan
metrics$simp = simp
metrics$invsimp = invsimp
metrics$ASV_rich = ASV_rich
all_metrics = metrics %>%
  gather(-sample,-day,-treatment, key="index",value="value")

ggplot(all_metrics, aes(x=day, y=value,color=treatment)) +
  geom_boxplot(aes(group=interaction(day,treatment)),position = "dodge")+
  geom_point(position=position_dodge(width=0.75))+
  scale_color_manual(values=treats)+
  facet_grid(index~., scales = "free_y")+
  theme_bw()

######
##stats on the diversity metrics
######
ASV_rich_aov=aov(ASV_rich~treatment, metrics)
summary(ASV_rich_aov)
TukeyHSD(ASV_rich_aov)

h_aov=aov(h~treatment, metrics)
summary(h_aov)
TukeyHSD(h_aov)

simp_aov=aov(simp~treatment, metrics)
summary(simp_aov)
TukeyHSD(simp_aov)

invsimp_aov=aov(invsimp~treatment, metrics)
summary(invsimp_aov)
TukeyHSD(invsimp_aov)

######
##heatmap
######

library(pheatmap)
heatmap = as.matrix(feature_abundance_ge3samples[,-115])
pheatmap(heatmap,scale="column",cluster_cols = F,cluster_rows = F)
