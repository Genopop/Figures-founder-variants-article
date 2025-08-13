#####################
##
## Claudia Moreau
## 
## 07/08/2025
##
## Redo fig 1 with all variants not only enriched
#####################

### Packages
rm(list = ls());
library('ggplot2')
library("ggpubr")
library("gridExtra")
library("ggpmisc")
library("stringr")
library("stringi")
library("tidyr")
library(grid)
library(cowplot)

options(scipen=999)


################################################################# Fig.1 ############################################################
######### Graph Gnomad/Imput MAF on ALL Variants in QcP
## Variants with freq of Gnomad
mat = matrix(ncol = 6, nrow = 0)
info_Gnomad=data.frame(mat)
names(info_Gnomad) <- c("CHR","SNP","A1","A2","MAF","NCHROBS")
for (i in 1:22){
  VTA_variants_Gnomad <- read.table(paste0("gnomad_clinVar_01082024_chr",i,".frq"), header =T)
  info_Gnomad <- rbind(info_Gnomad,VTA_variants_Gnomad)
}
## Get info only for CBNFE - freq
freq_NTNFE <- data.frame(SNP=info_Gnomad$SNP,Freq_NTNFE=info_Gnomad$AC_non_topmed_nfe/info_Gnomad$AN_non_topmed_nfe)

## Open freq of all chr Imput CaG QcP
freq_Imput <- read.table("mymerged_all_chr_Imput_VTA_WQ.frq", header =T)
freq_Imput <- freq_Imput[,-c(1,3,4)]
freq_Imput$SNP <- gsub(':','_', freq_Imput$SNP)

## Variants with freq of Gnomad
## Get MAF for SLSJ in Imput
freq_SLSJ_Imput <- read.table("mymerged_all_chr_Imput_VTA_SAG.afreq", header =F)
freq_SLSJ_Imput <- freq_SLSJ_Imput[,-c(1,3,4,5)]
freq_SLSJ_Imput$SNP <- gsub(':','_', freq_SLSJ_Imput$V2)
freq_SLSJ_Imput <- freq_SLSJ_Imput[,c(4,2,3)]

freq_nfe_imput <- merge(freq_NTNFE,freq_Imput,by=1,all=TRUE)
freq_nfe_imput <- merge(freq_nfe_imput,freq_SLSJ_Imput,by=1,all=TRUE)

dimnames(freq_nfe_imput)[[2]] <- c('SNP','Freq_NTNFE','MAF_WQ','n_WQ','MAF_SLSJ','n_SLSJ')
## Put all NAs to 0
freq_nfe_imput[is.na(freq_nfe_imput)] <- 0

## And if only keep imputations, non topmed at 21066 - remove freq smaller than smallest dataset
n_nfe <- max(info_Gnomad$AN_non_topmed_nfe)
freq_nfe_imput <- freq_nfe_imput[freq_nfe_imput$MAF_WQ >= 1/(n_nfe+1) | freq_nfe_imput$Freq_NTNFE >= 1/(n_nfe),]
##Remove snps not in WGS WQ or unreliable in imputations
## Unreliable in imput and absent in WGS
to_remove_imput_file <- 'ImputClaudia_1729WGSCaG_false_positive_snps.txt'
snps2remove_imput_df <- read.table(to_remove_imput_file)
## SNPs not in WGS
snps_not_in_wgs <- read.table("variants_in_imput_not_in_WGS.txt")
snps_not_in_wgs$V1 <- gsub(':','_',snps_not_in_wgs$V1)

freq_nfe_imput <- freq_nfe_imput[!is.element(freq_nfe_imput$SNP,snps2remove_imput_df$V1),]
freq_nfe_imput <- freq_nfe_imput[!is.element(freq_nfe_imput$SNP,snps_not_in_wgs$V1),]
## Also remove variants at " 5% in NFE
freq_nfe_imput <- freq_nfe_imput[freq_nfe_imput$Freq_NTNFE<=0.05,]

## Stats
length(freq_nfe_imput[freq_nfe_imput$MAF_WQ>0,1])/length(freq_nfe_imput[,1])
length(freq_nfe_imput[freq_nfe_imput$MAF_SLSJ>0,1])/length(freq_nfe_imput[,1])
length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE>0,1])/length(freq_nfe_imput[,1])
# > length(freq_nfe_imput[freq_nfe_imput$MAF_WQ>0,1])/length(freq_nfe_imput[,1])
# [1] 0.156364
# > length(freq_nfe_imput[freq_nfe_imput$MAF_SLSJ>0,1])/length(freq_nfe_imput[,1])
# [1] 0.098861
# > length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE>0,1])/length(freq_nfe_imput[,1])
# [1] 0.9288953

## Color interesting snps
max_threshold <- 0.005
min_threshold <- 0.0025
freq_nfe_imput$col_WQ <- rep(0,length(freq_nfe_imput[,1]))
freq_nfe_imput$col_WQ[freq_nfe_imput$Freq_NTNFE<min_threshold & freq_nfe_imput$MAF_WQ>max_threshold] <- 1
freq_nfe_imput$col_WQ[freq_nfe_imput$MAF_WQ<min_threshold & freq_nfe_imput$Freq_NTNFE>max_threshold] <- 2
freq_nfe_imput$col_SLSJ <- rep(0,length(freq_nfe_imput[,1]))
freq_nfe_imput$col_SLSJ[freq_nfe_imput$Freq_NTNFE<min_threshold & freq_nfe_imput$MAF_SLSJ>max_threshold] <- 1
freq_nfe_imput$col_SLSJ[freq_nfe_imput$MAF_SLSJ<min_threshold & freq_nfe_imput$Freq_NTNFE>max_threshold] <- 2

## Do figure QcP
gnomad_vs_founder_enriched_only <- ggplot(data = freq_nfe_imput,aes(x = MAF_WQ, y = Freq_NTNFE,col=factor(col_WQ))) +
  geom_point(size=0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values=c('grey','orange','purple'),labels = c("None" ,"QcP/SLSJ", "gnomAD"))+
  theme(axis.text.x = element_text( size = 10),
        axis.title=element_text(size=10),
        axis.text.y = element_text(size = 10),
        plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm")) +
  xlab("Frequency in the QcP") + ylab("Frequency in gnomAD nfe")+labs(col='More frequent in')  +xlim(0,0.02)+ylim(0,0.02) 

## Figure SLSJ
gnomad_vs_founder_SAG_enriched_only <- ggplot(data = freq_nfe_imput,aes(x = MAF_SLSJ, y = Freq_NTNFE,col=factor(col_SLSJ))) + 
  geom_point(size=0.1)+
  scale_color_manual(values=c('grey','orange','purple'),labels = c("None" ,"QcP/SLSJ", "gnomAD"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  theme(axis.text.x = element_text(size = 10),
        axis.title=element_text(size=10),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm")) +
  guides(color=guide_legend(override.aes = list(size=3)))+
  rremove("ylab") +
  xlab("Frequency in SLSJ") + ylab("Frequency in gnomAD nfe") +labs(col='More frequent in') +xlim(0,0.02)+ylim(0,0.02) 


##### --  FIGURE --
figure_1_freq_Gnomad <- ggarrange(gnomad_vs_founder_enriched_only, gnomad_vs_founder_SAG_enriched_only, ncol=2, common.legend = TRUE, legend="bottom",labels = c("A", "B"),label.y = 1,widths = c(1.2,1))## Save tiff format
tiff("Figure_1.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
figure_1_freq_Gnomad
dev.off()

## How many of these variants are more frequent in SLSJ and QcP?
length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE<min_threshold & freq_nfe_imput$MAF_WQ>max_threshold,1])
length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE<min_threshold  & freq_nfe_imput$MAF_SLSJ>max_threshold,1])
length(freq_nfe_imput[freq_nfe_imput$MAF_WQ<min_threshold  & freq_nfe_imput$Freq_NTNFE>max_threshold,1])
length(freq_nfe_imput[freq_nfe_imput$MAF_SLSJ<min_threshold  & freq_nfe_imput$Freq_NTNFE>max_threshold,1])
# > length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE<min_threshold & freq_nfe_imput$MAF_WQ>max_threshold,1])
# [1] 2
# > length(freq_nfe_imput[freq_nfe_imput$Freq_NTNFE<min_threshold  & freq_nfe_imput$MAF_SLSJ>max_threshold,1])
# [1] 40
#   > length(freq_nfe_imput[freq_nfe_imput$MAF_WQ<min_threshold  & freq_nfe_imput$Freq_NTNFE>max_threshold,1])
# [1] 5
# > length(freq_nfe_imput[freq_nfe_imput$MAF_SLSJ<min_threshold  & freq_nfe_imput$Freq_NTNFE>max_threshold,1])
# [1] 8

## total snps num:
length(freq_nfe_imput[,1])
## [1] 9043
