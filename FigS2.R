################################
#
# September 2024 by MICHEL Elisa
#
# Do Figure in article Founder variants: Supp Fig.2
################################
rm(list = ls());
### Packages
library("tidyr")
library("stringr")
library("dplyr")
library("ggplot2")

## Open freq of all chr sequencing CaG and merge them for Gnomad
mat = matrix(ncol = 6, nrow = 0)
info_Gnomad=data.frame(mat)
names(info_Gnomad) <- c("CHR","SNP","A1","A2","MAF","NCHROBS")
for (i in 1:22){
  VTA_variants_Gnomad <- read.table(paste0("/lustre03/project/6033529/saguenay_disease/results/gnomad_freq/gnomad_clinVar_01082024_chr",i,".frq"), header =T)
  info_Gnomad <- rbind(info_Gnomad,VTA_variants_Gnomad)
}
## Open freq of all chr WGS CaG QcP
freq_WQ <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Data/WGS/WQ/mymerged_all_chr_CaG_seq_VTA_WQ.frq", header =T)
freq_WQ <- freq_WQ[,-c(1,3,4)]

############################ Gnomad non_topmed_nfe = NTNFE ######################################
## Get info only for CBNFE - freq
freq_NTNFE <- as.data.frame(info_Gnomad$AC_non_topmed_nfe/info_Gnomad$AN_non_topmed_nfe)
names(freq_NTNFE) <- "Freq_NTNFE"
freq_NTNFE$SNP <- info_Gnomad$SNP
freq_NTNFE$AN_non_topmed_nfe <- info_Gnomad$AN_non_topmed_nfe
## Add the missing variants in gnomAD
m=nrow(freq_NTNFE)
for (i in 1:nrow(freq_WQ)){
  if (!freq_WQ[i,1] %in% freq_NTNFE$SNP){
    m=m+1
    freq_NTNFE[m,1] <- 1/(21066+1)
    freq_NTNFE[m,2] <- freq_WQ[i,1]
    freq_NTNFE[m,3] <- 21066
  }
}
freq_NTNFE <- freq_NTNFE[complete.cases(freq_NTNFE),]
m=nrow(freq_WQ)
for (i in 1:nrow(freq_NTNFE)){
  if (!freq_NTNFE[i,2] %in% freq_WQ$SNP){
    m=m+1
    freq_WQ[m,2] <- 1/(3704+1)
    freq_WQ[m,1] <- freq_NTNFE[i,2]
    freq_WQ[m,3] <- 3704
  }
  print(i)
}
## Merge the info
VTA_variants_NTNFE <- merge(x = freq_NTNFE, y = freq_WQ, by = c("SNP"), all = T)
colnames(VTA_variants_NTNFE)<- c("SNP","MAF_NTNFE","NCHROBs_NTNFE","MAF_WQ","NCHROBs_WQ")
## Keep with min freq in 1 of 2 data
VTA_variants_NTNFE$Status <- ifelse((VTA_variants_NTNFE$MAF_NTNFE <= 1/(3704+1) & VTA_variants_NTNFE$MAF_WQ <= 1/(3704+1)),
                                   "to_delete", "OK")
VTA_variants_NTNFE_clean <- VTA_variants_NTNFE[VTA_variants_NTNFE$Status == "OK",]
## Calculate enrichment for QcP
m <- nrow(VTA_variants_NTNFE_clean)
for (i in 1:m){
  VTA_variants_NTNFE_clean$enrich_stand_WQ[i]= (as.numeric(VTA_variants_NTNFE_clean$MAF_WQ[i]) - as.numeric(VTA_variants_NTNFE_clean$MAF_NTNFE[i]))/ as.numeric(VTA_variants_NTNFE_clean$MAF_WQ[i])
}
## Calculate enrichment for gnomAD
for (i in 1:m){
  VTA_variants_NTNFE_clean$enrich_stand_gnomAD[i]= (as.numeric(VTA_variants_NTNFE_clean$MAF_NTNFE[i]) - as.numeric(VTA_variants_NTNFE_clean$MAF_WQ[i]))/ as.numeric(VTA_variants_NTNFE_clean$MAF_NTNFE[i])
}
## Add color 
VTA_variants_NTNFE_clean$enriched_WQ <- ifelse(VTA_variants_NTNFE_clean$enrich_stand_WQ >=0.1, "enriched_WQ","not_enriched")
VTA_variants_NTNFE_clean$enriched_gnomAD <- ifelse(VTA_variants_NTNFE_clean$enrich_stand_gnomAD >=0.1, "enriched_gnomAD","not_enriched")
for (i in 1:nrow(VTA_variants_NTNFE_clean)){
  if (VTA_variants_NTNFE_clean[i,9] == "not_enriched" & VTA_variants_NTNFE_clean[i,10] == "not_enriched"){
    VTA_variants_NTNFE_clean[i,11] <- "not_enriched"
  } else if (VTA_variants_NTNFE_clean[i,9] ==  "enriched_WQ" & VTA_variants_NTNFE_clean[i,10] == "not_enriched"){
    VTA_variants_NTNFE_clean[i,11] <- "enriched_WQ" 
  } else if (VTA_variants_NTNFE_clean[i,10] ==  "enriched_gnomAD" & VTA_variants_NTNFE_clean[i,9] == "not_enriched"){
    VTA_variants_NTNFE_clean[i,11] <- "enriched_gnomAD"
  } else {
    VTA_variants_NTNFE_clean[i,11] <- "other"
  }
}
## Look at data
table(VTA_variants_NTNFE_clean$V11)
## Do the figure
figure_method_compare_gnomAD_WGS_WQ <- ggplot(data = VTA_variants_NTNFE_clean,
                  aes(x = MAF_WQ, y = MAF_NTNFE,col=V11))+ geom_point()+
  scale_shape_manual(values=c(1,5))+
  theme(axis.text.x = element_text(face="bold", size = 15),axis.title=element_text(size=18,face="bold"),
        axis.text.y = element_text(face="bold", size = 15),legend.text = element_text(size=16)) + xlim(0,0.02)+ylim(0,0.02)+
  theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
  xlab("Variant frequency in ClinVar for QcP WGS") + ylab("Variant frequency in gnomAD")+
  scale_color_manual(labels = c("Variants with RFD≥10% in gnomAD (n=389)", "Variants with RFD≥10% in QcP (n=1537)",
                                "Variants without RFD≥10% (n=111)"),
                     values = c("deeppink3", "deepskyblue","grey"))

## Save in PDF format
outfile <- paste0("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/figure_supp_3.pdf")
pdf(outfile,  width = 30, height = 10)
print(figure_method_compare_gnomAD_WGS_WQ)
dev.off()
## Save in PNG format
ggsave("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/figure_supp_2.png", width = 12, height = 8)

