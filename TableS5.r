#!/bash/env R

#######################
##
##
## Claudia
##
## Make correspondance between WGS and imputations
##
## juilet 2024
##
########################

# BiocManager::version()
# library ('BiocManager')
library('snpStats')
library(ggplot2)
library(stringr)

colfunc <- colorRampPalette(c("red", "grey",'blue'))


WGS_bedfile <- "mymerged_all_chr_CaG_seq_VTA_WQ_rfd0.1snps.bed"
WGS_bimfile <- "mymerged_all_chr_CaG_seq_VTA_WQ_rfd0.1snps.bim"
WGS_famfile <- "mymerged_all_chr_CaG_seq_VTA_WQ_rfd0.1snps.fam"

WGS_genos <- read.plink(WGS_bedfile, WGS_bimfile,WGS_famfile)
WGS_genos$genotypes

WGS_ind_sum <- row.summary(WGS_genos$genotypes)
WGS_snp_sum <- col.summary(WGS_genos$genotypes)

## Topmed imput
bedfile <- "topmed_imput_WQ_rfd0.1snps.bed"
bimfile <- "topmed_imput_WQ_rfd0.1snps.bim"
famfile <- "topmed_imput_WQ_rfd0.1snps.fam"
## Local imput on cag+schizo WGS
bedfile <- "local_imput_VTA_WQ_rfd0.1snps.bed"
bimfile <- "local_imput_VTA_WQ_rfd0.1snps.bim"
famfile <- "local_imput_VTA_WQ_rfd0.1snps.fam"

bim_df <- read.table(bimfile)
	
imput_genos <- read.plink(bedfile, bimfile,famfile)

imput_ind_sum <- row.summary(imput_genos$genotypes)
imput_snp_sum <- col.summary(imput_genos$genotypes)

compare_obj <- sm.compare(WGS_genos$genotypes, imput_genos$genotypes, row.wise = TRUE, col.wise = TRUE)
snps_df <- as.data.frame(compare_obj$col.wise)

hom_switch <- dimnames(snps_df[snps_df$Hom.switch>0,])[[1]]

snps_df$inhetWGS <- WGS_snp_sum$P.AB[match(dimnames(snps_df)[[1]],dimnames(WGS_snp_sum)[[1]])]*WGS_snp_sum$Calls[match(dimnames(snps_df)[[1]],dimnames(WGS_snp_sum)[[1]])]
snps_df$inhetimput <- imput_snp_sum$P.AB[match(dimnames(snps_df)[[1]],dimnames(imput_snp_sum)[[1]])]*imput_snp_sum$Calls[match(dimnames(snps_df)[[1]],dimnames(imput_snp_sum)[[1]])]
snps_df$hetdiff <- round(snps_df$inhetWGS-snps_df$inhetimput)
snps_df$false_pos <- round(snps_df$inhetimput-snps_df$Het.agree)
snps_df$false_neg <- round(snps_df$inhetWGS-snps_df$Het.agree)
snps_df$num_genos_agree_in_false_pos_neg_homswitch <- snps_df$Agree
snps_df$num_genos_agree_in_false_pos_neg_homswitch[snps_df$Hom.switch==0 & snps_df$Het.Hom==0] <- 0

snps_df$chr <- bim_df$V1[match(dimnames(snps_df)[[1]],bim_df$V2)]
snps_df$pos <- bim_df$V4[match(dimnames(snps_df)[[1]],bim_df$V2)]
snps_df <- snps_df[order(snps_df$chr,snps_df$pos),]

## Make table of switches
my_imput <- 'local'
num_snps_total <- 1302
num_snps_in_both <- 1047
num_inds <- 1852
num_missing_snps<-num_snps_total-num_snps_in_both
num_snps_hom_switch<-length(snps_df[snps_df$Hom.switch>0,1])
num_snps_agree<-length(snps_df[snps_df$Agree + snps_df$NA.disagree==num_inds,1])
num_genos_agree<-(num_snps_agree*num_inds) + sum(snps_df$num_genos_agree_in_false_pos_neg_homswitch)
num_genos_false_pos_hom_switch <- sum(snps_df$Hom.switch[snps_df$Hom.switch>0]) 
## without hom switch
num_snps_false_positives_only<-length(snps_df[snps_df$false_pos>0 & snps_df$false_neg==0  & snps_df$Hom.switch==0,1])
num_snps_false_negatives_only<-length(snps_df[snps_df$false_neg>0 & snps_df$false_pos==0  & snps_df$Hom.switch==0,1])
num_genos_false_positives_only<-sum(snps_df$false_pos[snps_df$false_pos>0 & snps_df$false_neg==0 & snps_df$Hom.switch==0])
num_genos_false_negatives_only<-sum(snps_df$false_neg[snps_df$false_neg>0 & snps_df$false_pos==0 & snps_df$Hom.switch==0])
num_snps_false_both<-length(snps_df[snps_df$false_pos>0 & snps_df$false_neg>0 & snps_df$Hom.switch==0,1])
num_genos_false_both_pos<-sum(snps_df$false_pos[snps_df$false_pos>0 & snps_df$false_neg>0 & snps_df$Hom.switch==0])
num_genos_false_both_neg<-sum(snps_df$false_neg[snps_df$false_pos>0 & snps_df$false_neg>0 & snps_df$Hom.switch==0])

my_stats_snps_df <- data.frame(num=c(num_missing_snps,num_snps_hom_switch,num_snps_agree,num_snps_false_positives_only,num_snps_false_negatives_only,num_snps_false_both))
dimnames(my_stats_snps_df)[[1]] <- c('num_missing_snps','num_snps_hom_switch','num_snps_agree','num_snps_false_positives_only','num_snps_false_negatives_only','num_snps_false_both')
my_stats_snps_df$prop <- my_stats_snps_df$num / num_snps_in_both
## For the overall missing
my_stats_snps_df$prop[1] <- my_stats_snps_df$num[1] / num_snps_total

outfile <- paste(my_imput,'Imput_WGSCaG_SNPs_compare.txt',sep='')
write.table(my_stats_snps_df,outfile,sep='\t',quote=FALSE)

my_stats_genos_df <- data.frame(num=c(num_genos_agree,num_genos_false_positives_only+num_genos_false_both_pos,num_genos_false_negatives_only+num_genos_false_both_neg,num_genos_false_pos_hom_switch))
dimnames(my_stats_genos_df)[[1]] <- c('num_genos_agree','false_positives','false_negatives','hom_switch')
my_stats_genos_df$prop <- my_stats_genos_df$num / (num_snps_in_both * num_inds)
outfile <- paste(my_imput,'Imput_WGSCaG_genos_compare.txt',sep='')
write.table(my_stats_genos_df,outfile,sep='\t',quote=FALSE)

false<-dimnames(snps_df[snps_df$false_pos>0 & snps_df$false_neg==0  & snps_df$Hom.switch==0,])[[1]]
fasle2<-dimnames(snps_df[snps_df$false_pos>0 & snps_df$false_neg>0 & snps_df$Hom.switch==0,])[[1]]
fasle3<-dimnames(snps_df[snps_df$Hom.switch>0,])[[1]]
outfile <- paste(my_imput,'Imput_WGSCaG_fasle_positive_snps.txt',sep='')
write.table(c(false,fasle2,fasle3),outfile,sep='\t',quote=FALSE,row.names = FALSE)

snps_df$chr_pos <- c(1:length(snps_df$pos))
snps_df$faslse_pos_prop <- snps_df$false_pos/snps_df$Het.agree
# snps_df$diff_imput_agree <- round(snps_df$inhetimput - snps_df$Het.agree)
outfile <- paste(my_imput,'Imput_WGSCaG_descriptive_table.txt',sep='')
write.table(snps_df,outfile,sep='\t',quote=FALSE,row.names = FALSE)

