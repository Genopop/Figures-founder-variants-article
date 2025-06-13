################################
#
# Claudia Moreau
#
# 2024
#
# MAking UMAP from PCA and plots for WGS and imputed data
#
################################


rm(list = ls());

library(ggplot2)
library(umap)
library(pals)
library("ggpubr")
library(stringr)


## umap ## WGS ## Fig suppl 3


cag_regions_file <- 'region_cartagene.30k.txt'
cag_regions_df<-read.table(cag_regions_file , header=F,sep="")
dimnames(cag_regions_df)[[2]] <- c('num','file111','Cregion')
cag_regions_df[1:10,]
length(cag_regions_df[cag_regions_df$Cregion=='Saguenay',1])

cag_birth_file <- 'cartagene_birth_country.30k.txt'
cag_birth_df<-read.table(cag_birth_file , header=T,sep="",colClasses=rep('numeric',3))
cag_birth_df$birth <- cag_birth_df$COUNTRY_BIRTH
cag_birth_df[1:10,]
cag_birth_df$birth [cag_birth_df$COUNTRY_BIRTH==99 | is.na(cag_birth_df$COUNTRY_BIRTH)] <- cag_birth_df$COUNTRY_BIRTH_OTHER[cag_birth_df$COUNTRY_BIRTH==99 | is.na(cag_birth_df$COUNTRY_BIRTH)] + 100
cag_birth_df[1:10,]
cag_birth_df$birth[!duplicated(cag_birth_df$birth)]

cag_continent_file <- 'cartagene_some_continent.txt'
cag_continent_df<-read.table(cag_continent_file , header=T,sep="\t",colClasses=c('numeric',rep('character',2),'numeric'))
cag_birth_df$continent <- cag_continent_df$continent[match(cag_birth_df$birth,cag_continent_df$country_num)]


######### UMAP of WGS #################
evecfile <- 'wgs_CaG_maf0.05_geno0.02_mind0.02_snpsonly_nomulti_allchrs_mask_hg38_prune200_5_0.1.evec'
evalfile <- 'wgs_CaG_maf0.05_geno0.02_mind0.02_snpsonly_nomulti_allchrs_mask_hg38_prune200_5_0.1.eval'

eigenvals <- read.table(evalfile , header=F)
eigenvals$variance <- eigenvals[,1] / sum(eigenvals[,1])
eigenval1 <- eigenvals[1,1] / sum(eigenvals[,1])
eigenval2 <- eigenvals[2,1] / sum(eigenvals[,1])

outfile <- (paste(evalfile,'.elbow.variance.jpg',sep=''))
jpeg(filename=outfile,width = 8.5, height = 8.5 ,res=300,units="in", pointsize = 12)
ggplot(eigenvals, aes(x=c(1:length(eigenvals[,1])),y=variance)) +
	geom_point() +
	xlim(c(0,20))
dev.off()

pc_num<-3
outfile <- (paste(evecfile,'.umap.',pc_num,'pc.txt',sep=''))
umaplayout_wg <- read.table(outfile,header=TRUE)
umaplayout_wg[1:10,]
## nbr of sag in the sag cluster
length(umaplayout_wg[umaplayout_wg$birth_country_pop=='Saguenay',9])
length(umaplayout_wg[umaplayout_wg$birth_country_pop=='Saguenay' & umaplayout_wg$clusters==1,9])

## graph
my_cols <- as.vector(glasbey())
p <- ggplot(umaplayout_wg, aes(x=Dim1,y=Dim2,color=factor(clusters))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="bottom",text = element_text(size = 20)) +
	guides(color = guide_legend(title.position = "left",override.aes = list(size=4))) +
	labs(title ='',x='Dimension 1',y='Dimension 2',color='Clusters')	

legend <- get_legend(p)

p <- ggplot(umaplayout_wg, aes(x=Dim1,y=Dim2,color=factor(clusters))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="none",text = element_text(size = 20)) +
	labs(title ='',x='Dimension 1',y='Dimension 2')	

my_cols <- as.vector(kelly())
## ok choose a more evident color for Sag
my_cols[9] <- 'red'
my_cols[7] <- 'pink'
p1 <- ggplot(umaplayout_wg, aes(x=Dim1,y=Dim2,color=(birth_country_pop))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="bottom",legend.text = element_text(size=4),text = element_text(size = 20)) +
	guides(color = guide_legend(title.position = "top",override.aes = list(size=4))) +
	labs(title ='',x='Dimension 1',y='Dimension 2',color='Birth country or recruitement place')

legend1 <- get_legend(p1)

p1 <- ggplot(umaplayout_wg, aes(x=Dim1,y=Dim2,color=(birth_country_pop))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="none",text = element_text(size = 20)) +
	labs(title ='',x='Dimension 1',y='Dimension 2')

k=3
outfile <- (paste(evecfile,'.umap.',pc_num,'pc_',k,'clusters_QcSag.png',sep=''))
png(filename=outfile,width = 11, height = 5 ,res=600,units="in")
ggarrange(p1, p,as_ggplot(legend1),as_ggplot(legend), ncol = 2, nrow = 2,labels = c("A", "B",'',''),heights = c(2, 0.7))
dev.off()


### Imputations on CaG clusters #### Fig suppl 4
setwd('imputation/PCA')

my_file <- 'CaG_imputed_allchips.INFO0.3.maf0.05_pruned50_5_0.2'
evecfile <- paste(my_file,'.eigenvec',sep='') # evec file
evalfile <- paste(my_file,'.eigenval',sep='') # eval fil

mypca <- read.table(evecfile , header=F,sep="")
mypca$sample_tmp <- sapply(str_split(mypca$V1,pattern=':'),"[[",1)
mypca$sample <- sapply(str_split(mypca$sample_tmp,pattern='_'),"[[",1)
dimnames(mypca)[[2]] <- c('fam','id',paste('PC',c(1:20),sep=''),'tmp','sample')

eigenvals <- read.table(evalfile , header=F)
eigenvals$variance <- eigenvals[,1] / sum(eigenvals[,1])
eigenval1 <- eigenvals[1,1] / sum(eigenvals[,1])
eigenval2 <- eigenvals[2,1] / sum(eigenvals[,1])

 
pc_num=5
my_umap <- (umap(mypca[,c(2:pc_num+1)]))
my_umap$layout
umaplayout <- as.data.frame(my_umap$layout)
dimnames(umaplayout)[[2]] <- c('Dim1','Dim2')

umaplayout$sample <- mypca$sample

umaplayout$pop <- cag_regions_df$Cregion[match (umaplayout$sample,cag_regions_df$file111)]
umaplayout$continent <- cag_birth_df$continent[match(umaplayout$sample,cag_birth_df$file_111)]
umaplayout$birth <- cag_birth_df$birth[match (umaplayout$sample,cag_regions_df$file111)]
umaplayout$birth_country <- cag_continent_df$country[match(umaplayout$birth,cag_continent_df$country_num)]

umaplayout$birth_country_pop<-umaplayout$birth_country
umaplayout$birth_country_pop [umaplayout$birth_country=="CANADA" & !is.na(umaplayout$birth_country)] <- umaplayout$pop [umaplayout$birth_country=="CANADA" & !is.na(umaplayout$birth_country)]

umaplayout$continent_pop<-umaplayout$continent
umaplayout$continent_pop [umaplayout$birth_country=="CANADA" & !is.na(umaplayout$birth_country_pop)] <- umaplayout$pop [umaplayout$birth_country=="CANADA" & !is.na(umaplayout$birth_country_pop)]

k=5
res.km <- kmeans(scale(umaplayout[, c(1,2)]), k, nstart = 25)
cluster5 <- res.km$cluster
umaplayout$cluster5 <- factor(res.km$cluster)
## Sag cluster = 3

umaplayout$recluster <- 0
umaplayout$recluster[umaplayout$cluster5==4 | umaplayout$cluster5==2 | umaplayout$cluster5==5] <- 2
umaplayout$recluster[umaplayout$cluster5==3] <- 1

my_cols <- as.vector(glasbey())
p <- ggplot(umaplayout, aes(x=Dim1,y=Dim2,color=factor(recluster))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="bottom",text = element_text(size = 20)) +
	guides(color = guide_legend(title.position = "left",override.aes = list(size=4))) +
	labs(title ='',x='Dimension 1',y='Dimension 2',color='Clusters')	

legend <- get_legend(p)
as_ggplot(legend)

p <- ggplot(umaplayout, aes(x=Dim1,y=Dim2,color=factor(recluster))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="none",text = element_text(size = 20)) +
	labs(title ='',x='Dimension 1',y='Dimension 2')	

my_cols <- as.vector(kelly())


my_cols[8] <- "#0067A5"
my_cols[11] <- 'red'
my_cols[10] <- "#C2B280"
my_cols[12] <- "#008856"
my_cols[7] <- 'pink'

p1 <- ggplot(umaplayout[!is.na(umaplayout$continent_pop) & umaplayout$continent_pop!='',], aes(x=Dim1,y=Dim2,color=(continent_pop))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="bottom",legend.text = element_text(size=4),text = element_text(size = 20)) +
	guides(color = guide_legend(title.position = "top",override.aes = list(size=4))) +
	labs(title ='',x='Dimension 1',y='Dimension 2',color='Birth continent or recruitement place')

legend1 <- get_legend(p1)
as_ggplot(legend1)

p1 <- ggplot(umaplayout[!is.na(umaplayout$continent_pop) & umaplayout$continent_pop!='',], aes(x=Dim1,y=Dim2,color=(continent_pop))) +
	geom_point(cex=0.5) +
	scale_color_manual(values=c(my_cols)) +
	theme(legend.position="none",text = element_text(size = 20)) +
	labs(title ='',x='Dimension 1',y='Dimension 2')

k=5
outfile <- (paste(evecfile,'.umap.',pc_num,'pc_',k,'clusters_QcSag.png',sep=''))
png(filename=outfile,width = 11, height = 5 ,res=600,units="in")
ggarrange(p1, p,as_ggplot(legend1),as_ggplot(legend), ncol = 2, nrow = 2,labels = c("A", "B",'',''),heights = c(1.5, 1))
dev.off()
