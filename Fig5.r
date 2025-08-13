################################
#
# Claudia Moreau 11 juin 2025
#
# Do Figure in article Founder variants: Fig.5
################################
rm(list = ls());
### Packages
library("dplyr")
library("tidyr")
library("stringr")
library("ggpubr")
library("reshape2")
library("ggplot2")
library('Rmisc')

print("Graphs number of carriers (enriched variants)")

## Known variants
variants_known <- read.table("variants_reported_in_SLSJ.txt")
## Get enriched variants positions
variants_enriched <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Delete variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr14_94380925_T_A","chr4_102901325_AAC_A")
variants_enriched <- variants_enriched[!variants_enriched$SNP %in% variants_to_freq,]
## Remove variants with one star classification except known variants
toremove <- read.table("snps_clinvar_1star.txt", header = FALSE)
known_tokeep <- read.table("known_founder_variants_in4reviews.txt", header = FALSE)
known_tokeep$V1 <- gsub(":", "_", known_tokeep$V1)
toremove <- as.data.frame(toremove[!toremove$V1 %in% known_tokeep$V1,])
variants_enriched_2stars <- variants_enriched[!variants_enriched$SNP %in% toremove$V1,]
## Get variants known SNP
variants_known_SNP <- variants_enriched[variants_enriched$ID %in% variants_known[,1],1]
variants_known_SNP_2stars <- variants_enriched_2stars[variants_enriched_2stars$ID %in% variants_known[,1],1]

####################### known variants only and all (known and new) side by side - SLSJ ############################
## Open file with the individuals in SLSJ
## Get individuals info in SLSJ
individual_file_SLSJ <- read.table(paste0('Imput/SLSJ/CaG_Imput_VTA_SAG_enriched_in_WQ_chr_1.tfam'))
## Open file with all info
df_final <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Remove variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
df_final <- df_final[!df_final$SNP %in% variants_to_freq,]
## Keep variants founder in SLSJ, UQc and QcP
df_final$Founder_SAG <- ifelse(is.na(df_final$Ind_Imput_SAG), 
							   df_final$founder_WGS_SAG,df_final$founder_Imput_SAG)
df_final$Founder_RQ <- ifelse(is.na(df_final$Ind_Imput_RQ), 
							  df_final$founder_WGS_RQ,df_final$founder_Imput_RQ)
df_final$Founder_WQ <- ifelse(is.na(df_final$Ind_Imput_WQ), 
							  df_final$founder_WGS_WQ,df_final$founder_Imput_WQ)
## Remove variants with one star classification
df_final_2stars <- df_final[!df_final$SNP %in% toremove$V1,]
## Select founder variants
Founder_SLSJ_RQ <- df_final[df_final$Founder_SAG == "Founder" | df_final$Founder_RQ == "Founder"| df_final$Founder_WQ == "Founder", ]
Founder_SLSJ_RQ_2stars <- df_final_2stars[df_final_2stars$Founder_SAG == "Founder" | df_fina_2starsl$Founder_RQ == "Founder"| df_final_2stars$Founder_WQ == "Founder", ]

## Select carriers
carriers_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('Imput/SLSJ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
	table_ind <- as.data.frame(table_ind)
	for (icol in c(1:ncol(table_ind))){
		my_rows <- cbind(dimnames(table_ind)[[2]][icol],table_ind[!is.na(table_ind[icol]),icol])
		carriers_df <- rbind(carriers_df,my_rows)
	}
}
unique(carriers_df[,1])

## Change : in _ in variant name
carriers_df[,1] <- gsub('\\.', "_", carriers_df[,1])
## Remove NR imput
NR_var <- read.table("snps_to_remove_seuil1_allcriteria_homswitch.txt")
carriers_df <- carriers_df[!carriers_df[,1] %in% NR_var[,1],]
## Select only variants founder in SLSJ
carriers_df <- carriers_df[carriers_df[,1] %in% Founder_SLSJ_RQ$SNP,]
carriers_df_2stars <- carriers_df[carriers_df[,1] %in% Founder_SLSJ_RQ_2stars$SNP,]

## Merge all info between Imput and WGS
carriers_founder_all_SAG <- carriers_df
carriers_founder_all_SAG_2stars <- carriers_df_2stars
## Remove duplicates (homozygotes)
remove_dupes <- function(df) {
  df[!duplicated(paste(df$V1, df$V2, sep = '_')), ]
}
carriers_founder_all_SAG <- remove_dupes(carriers_df)
carriers_founder_all_SAG_2stars <- remove_dupes(carriers_df_2stars)

### Table to find the number of variants for each ind
ind_num_var_df <- as.data.frame(table(carriers_founder_all_SAG[,2]))
hist_like_df <- table(ind_num_var_df$Freq)
ind_num_var_df_2stars <- as.data.frame(table(carriers_founder_all_SAG_2stars[,2]))
hist_like_df_2stars <- table(ind_num_var_df_2stars$Freq)
## Select known SNPs only
know_slsj <- carriers_founder_all_SAG[is.element(carriers_founder_all_SAG[,1],variants_known_SNP),]
know_slsj_2stars <- carriers_founder_all_SAG_2stars[is.element(carriers_founder_all_SAG_2stars[,1],variants_known_SNP_2stars),]

ind_num_var_know_slsj_df <- as.data.frame(table(know_slsj[,2]))
hist_like_know_slsj_df <- table(ind_num_var_know_slsj_df$Freq)
ind_num_var_know_slsj_df_2stars <- as.data.frame(table(know_slsj_2stars[,2]))
hist_like_know_slsj_df_2stars <- table(ind_num_var_know_slsj_df_2stars$Freq)
################################### Founder variants - UQc ##############################################
## Open file with individuals in UQc
individual_file_RQ <- read.table(paste0('Imput/RQ/CaG_Imput_VTA_RQ_enriched_in_WQ_chr_1.tfam'))
#### Create table with all variants enriched (compared with gnomAD) - Imput
carriers_rq_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('Imput/RQ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
	table_ind <- as.data.frame(table_ind)
	for (icol in c(1:ncol(table_ind))){
		my_rows <- cbind(dimnames(table_ind)[[2]][icol],table_ind[!is.na(table_ind[icol]),icol])
		carriers_rq_df <- rbind(carriers_rq_df,my_rows)
	}
}
unique(carriers_rq_df[,1])

## Change : in _ in variant name
carriers_rq_df[,1] <- gsub('\\.', "_", carriers_rq_df[,1])
## Remove NR imput
carriers_rq_df <- carriers_rq_df[!carriers_rq_df[,1] %in% NR_var[,1],]
## Select only variants founder in SLSJ
carriers_rq_df <- carriers_rq_df[carriers_rq_df[,1] %in% Founder_SLSJ_RQ$SNP,]
## Remove variants with one star classification
carriers_rq_df_2stars <- carriers_rq_df[!carriers_rq_df$V1 %in% toremove$V1, ]

## Merge all info between Imput and WGS
carriers_enriched_all_RQ <- carriers_rq_df
carriers_enriched_all_RQ_2stars <- carriers_rq_df_2stars
## Remove duplicates (homozygotes)
carriers_enriched_all_RQ <- carriers_enriched_all_RQ[-which(duplicated(paste(carriers_enriched_all_RQ$V1,carriers_enriched_all_RQ$V2,sep='_'))),]
carriers_enriched_all_RQ_2stars <- carriers_enriched_all_RQ_2stars[-which(duplicated(paste(carriers_enriched_all_RQ_2stars$V1,carriers_enriched_all_RQ_2stars$V2,sep='_'))),]

### table to find the number of variants for each ind
ind_num_var_rq_df <- as.data.frame(table(carriers_enriched_all_RQ[,2]))
hist_like_rq_df <- table(ind_num_var_rq_df$Freq)
ind_num_var_rq_df_2stars <- as.data.frame(table(carriers_enriched_all_RQ_2stars[,2]))
hist_like_rq_df_2stars <- table(ind_num_var_rq_df_2stars$Freq)

## Select know SNPs only
know_RQ <- carriers_enriched_all_RQ[is.element(carriers_enriched_all_RQ[,1],variants_known_SNP),]
know_RQ_2stars <- carriers_enriched_all_RQ_2stars[is.element(carriers_enriched_all_RQ_2stars[,1],variants_known_SNP_2stars),]

ind_num_var_rq_known_df <- as.data.frame(table(carriers_enriched_all_RQ[is.element(carriers_enriched_all_RQ[,1],variants_known_SNP),2]))
hist_like_rq_known_df <- table(ind_num_var_rq_known_df$Freq)
ind_num_var_rq_known_df_2stars <- as.data.frame(table(carriers_enriched_all_RQ_2stars[is.element(carriers_enriched_all_RQ_2stars[,1],variants_known_SNP_2stars),2]))
hist_like_rq_known_df_2stars <- table(ind_num_var_rq_known_df_2stars$Freq)

########## Resampling #################
resampling_vec <- c(rep(0,(21472-length(ind_num_var_rq_df[,1]))),ind_num_var_rq_df[,2])
resampling_df <- data.frame(variants_num=c(0:4))
for (i in c(1:1000)){
	my_sample <- sample(resampling_vec,3589,replace = FALSE)
	my_table_sample <- as.data.frame(table(my_sample))
	resampling_df <- merge(resampling_df,my_table_sample,by=1,all=TRUE)
}
resampling_df[is.na(resampling_df)] <- 0
resampling_df <- resampling_df[,-1]
## Look at means and CI
mean_ci <- apply(resampling_df,MARGIN = 1,CI)
## 2 stars
resampling_vec_2stars <- c(rep(0,(21472-length(ind_num_var_rq_df_2stars[,1]))),ind_num_var_rq_df_2stars[,2])
resampling_df_2stars <- data.frame(variants_num=c(0:4))
for (i in c(1:1000)){
	my_sample_2stars <- sample(resampling_vec_2stars,3589,replace = FALSE)
	my_table_sample_2stars <- as.data.frame(table(my_sample_2stars))
	resampling_df_2stars <- merge(resampling_df_2stars,my_table_sample_2stars,by=1,all=TRUE)
}
resampling_df_2stars[,1:10]
resampling_df_2stars[is.na(resampling_df_2stars)] <- 0
resampling_df_2stars <- resampling_df_2stars[,-1]
## Look at means and CI
mean_ci_2stars <- apply(resampling_df_2stars,MARGIN = 1,CI)
## Do it again with only known
########## Resample #################
resampling_known_vec <- c(rep(0,(21472-length(ind_num_var_rq_known_df[,1]))),ind_num_var_rq_known_df[,2])
resampling_known_df <- data.frame(variants_num=c(0:4))
for (i in c(1:1000)){
	my_sample <- sample(resampling_known_vec,3589,replace = FALSE)
	my_table_sample <- as.data.frame(table(my_sample))
	resampling_known_df <- merge(resampling_known_df,my_table_sample,by=1,all=TRUE)
}
resampling_known_df[is.na(resampling_known_df)] <- 0
resampling_known_df <- resampling_known_df[,-1]
## Look at means and CI
mean_ci_known <- apply(resampling_known_df,MARGIN = 1,CI)
## 2 stars
resampling_known_vec_2stars <- c(rep(0,(21472-length(ind_num_var_rq_known_df_2stars[,1]))),ind_num_var_rq_known_df_2stars[,2])
resampling_known_df_2stars <- data.frame(variants_num=c(0:4))
for (i in c(1:1000)){
	my_sample_2stars <- sample(resampling_known_vec_2stars,3589,replace = FALSE)
	my_table_sample_2stars <- as.data.frame(table(my_sample_2stars))
	resampling_known_df_2stars <- merge(resampling_known_df_2stars,my_table_sample_2stars,by=1,all=TRUE)
}
resampling_known_df_2stars[,1:10]
resampling_known_df_2stars[is.na(resampling_known_df_2stars)] <- 0
resampling_known_df_2stars <- resampling_known_df_2stars[,-1]
## Look at means and CI
mean_ci_known_2stars <- apply(resampling_known_df_2stars,MARGIN = 1,CI)

###### Make figure 5 plot with all variants
### Make SLSJ graph
barplot_res_df <- data.frame(variants_num = c(0:6))
SLSJ_res <- as.numeric(hist_like_df)
SLSJ_res <- c(3589-sum(SLSJ_res),SLSJ_res)
barplot_res_df$value <- SLSJ_res
barplot_res_df$group <- 'Founder variants'

barplot_known_res_df <- data.frame(variants_num = c(0:6))
SLSJ_res <- as.numeric(hist_like_know_slsj_df)
SLSJ_res <- c(3589-sum(SLSJ_res),SLSJ_res)
barplot_known_res_df$value <- SLSJ_res
barplot_known_res_df$group <- 'Known founder variants'
barplot_res_df <- rbind(barplot_res_df,barplot_known_res_df)

graph_SAG <- ggplot(barplot_res_df, aes(x=variants_num, y=value/3589,fill=group)) +
	geom_bar(stat='identity', position = position_dodge(width = 0.9)) +
	labs(x = "Number of variants in SLSJ", y = "") + 
	scale_x_continuous( labels = c(0:6),breaks=c(0:6))+
	geom_text(aes(label=value),
			  vjust = -0.5, color="black", size=2,position = position_dodge(width = 0.9)) +
	theme(plot.title = element_text(size=10,  colour= "black" ),
		  axis.title.x = element_text(size=12,  colour = "black"),    
		  axis.title.y = element_text(size=12,  colour = "black"),    
		  axis.text.x = element_text(size=12,  colour = "black"), 
		  axis.text.y = element_text(size=12,  colour = "black"),
		  strip.text.x = element_text(size = 12,  colour = "black" ),
		  strip.text.y = element_text(size = 12,  colour = "black"),
		  plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm"),
		  legend.title = element_blank(),legend.text = element_text(size=8),
		  legend.position = "bottom")+
	scale_fill_manual(values=c("pink", "red"),labels=c("Founder variants","Known founder variants")) +
	scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1.0))


### Make UQc graph
barplot_res_df <- as.data.frame(t(mean_ci))
barplot_res_df$variants_num = c(0:4)
UQc_res <- as.numeric(hist_like_rq_df)
UQc_res <- c(21472-sum(UQc_res),UQc_res)
barplot_res_df$value <- UQc_res
barplot_res_df$group <- 'Founder variants'

barplot_known_res_df <- as.data.frame(t(mean_ci_known))
barplot_known_res_df$variants_num = c(0:4)
UQc_res <- c(as.numeric(hist_like_rq_known_df),0)
UQc_res <- c(21472-sum(UQc_res),UQc_res)
barplot_known_res_df$value <- UQc_res
barplot_known_res_df$group <- 'Known founder variants'

barplot_res_df <- rbind(barplot_res_df,barplot_known_res_df)

graph_UQc <- ggplot(barplot_res_df, aes(x=variants_num, y=mean/3589,fill=group)) +
	geom_bar(stat='identity', position = position_dodge(width = 0.9)) +
	geom_errorbar(aes(ymin = lower/3589, ymax = upper/3589),position = position_dodge(width = 0.9),width=0.5,linewidth=0.5) +
	labs(x = "Number of variants in UQc", y = "Proportion of individuals") + 
	scale_x_continuous( labels = c(0:6),breaks=c(0:6))+
	geom_text(aes(label=value),
			  vjust = -0.5, color="black", size=2,position = position_dodge(width = 0.9)) +
	theme(plot.title = element_text(size=10,  colour= "black" ),
		  axis.title.x = element_text(size=12,  colour = "black"),    
		  axis.title.y = element_text(size=12,  colour = "black"),    
		  axis.text.x = element_text(size=12,  colour = "black"), 
		  axis.text.y = element_text(size=12,  colour = "black"),
		  strip.text.x = element_text(size = 12,  colour = "black" ),
		  strip.text.y = element_text(size = 12,  colour = "black"),
		  plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm"),
		  legend.title = element_blank(),legend.text = element_text(size=8),
		  legend.position = "bottom")+
	scale_fill_manual(values=c("pink", "red"),labels=c("Founder variants","Known founder variants")) +
	scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1.0))


final_plot <- ggarrange(graph_UQc, graph_SAG,
						ncol = 2, nrow = 1,
						common.legend = TRUE,
						labels = c("A", "B"),label.y = 1,widths = c(1.2,1),
						legend = "bottom")


tiff("Figure_5.tiff",width =190, height = 132,units = 'mm', res=600 ,compression = "lzw")
final_plot
dev.off()

###### Make figure 5 subplots with variants with at least 2 stars classification 
### Make SLSJ graph
barplot_res_df_2stars <- data.frame(variants_num = c(0:6))
SLSJ_res_2stars <- rep(0, 7)
names(SLSJ_res_2stars) <- 0:6
SLSJ_res_2stars[names(hist_like_df_2stars)] <- as.numeric(hist_like_df_2stars)
SLSJ_res_2stars["0"] <- 3589 - sum(SLSJ_res_2stars[names(SLSJ_res_2stars) != "0"])
barplot_res_df_2stars$value <- SLSJ_res_2stars
barplot_res_df_2stars$group <- 'Founder variants'

barplot_known_res_df_2stars <- data.frame(variants_num = c(0:6))
SLSJ_known_res_2stars <- rep(0, 7)
names(SLSJ_known_res_2stars) <- 0:6
SLSJ_known_res_2stars[names(hist_like_know_slsj_df_2stars)] <- as.numeric(hist_like_know_slsj_df_2stars)
SLSJ_known_res_2stars["0"] <- 3589 - sum(SLSJ_known_res_2stars[names(SLSJ_known_res_2stars) != "0"])
barplot_known_res_df_2stars$value <- SLSJ_known_res_2stars
barplot_known_res_df_2stars$group <- "Known founder variants"
barplot_res_df_2stars <- rbind(barplot_res_df_2stars,barplot_known_res_df_2stars)

graph_SAG <- ggplot(barplot_res_df_2stars, aes(x=variants_num, y=value/3589,fill=group)) +
	geom_bar(stat='identity', position = position_dodge(width = 0.9), show.legend = FALSE) +
	labs(x = "", y = "", title = "> 1* ClinVar review status") + 
	scale_x_continuous( labels = c(0:6),breaks=c(0:6))+
	theme(plot.title = element_text(size=9,  colour = "black"),
		  axis.title.x = element_blank(),    
		  axis.title.y = element_blank(),    
		  axis.text.x = element_text(size=9,  colour = "black"), 
		  axis.text.y = element_text(size=9,  colour = "black"),
		  strip.text.x = element_text(size = 10,  colour = "black" ),
		  strip.text.y = element_text(size = 10,  colour = "black"),
		  legend.title = element_blank(),legend.text = element_text(size=8),
		  legend.position = "bottom")+
	scale_fill_manual(values=c("pink", "red"),labels=c("Founder variants","Known founder variants")) +
	scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1.0))
## Save in tiff format
ggsave("Figure_5_founder_SLSJ_2stars.tiff",graph_SAG,
       width =50, height = 40, units = 'mm', dpi=600,bg="white")

### Make UQc graph
barplot_res_df_2stars <- as.data.frame(t(mean_ci_2stars))
barplot_res_df_2stars$variants_num = c(0:4)
UQc_res_2stars <- as.numeric(hist_like_rq_df_2stars)
UQc_res_2stars <- c(21472-sum(UQc_res_2stars),UQc_res_2stars)
barplot_res_df_2stars$value <- UQc_res_2stars
barplot_res_df_2stars$group <- 'Founder variants'

barplot_known_res_df_2stars <- as.data.frame(t(mean_ci_known_2stars))
barplot_known_res_df_2stars$variants_num = c(0:4)
UQc_res_2stars <- c(as.numeric(hist_like_rq_known_df_2stars),0)
UQc_res_2stars <- c(21472-sum(UQc_res_2stars),UQc_res_2stars)
barplot_known_res_df_2stars$value <- UQc_res_2stars
barplot_known_res_df_2stars$group <- 'Known founder variants'

barplot_res_df_2stars <- rbind(barplot_res_df_2stars,barplot_known_res_df_2stars)

graph_UQc <- ggplot(barplot_res_df_2stars, aes(x=variants_num, y=mean/3589,fill=group)) +
	geom_bar(stat='identity', position = position_dodge(width = 0.9), show.legend = FALSE) +
	geom_errorbar(aes(ymin = lower/3589, ymax = upper/3589),position = position_dodge(width = 0.9),width=0.5,linewidth=0.05) +
	labs(x = "", y = "", title = "> 1* ClinVar review status") + 
	scale_x_continuous( labels = c(0:6),breaks=c(0:6))+
  theme(plot.title = element_text(size=9,  colour = "black"),
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.text.x = element_text(size=9,  colour = "black"), 
        axis.text.y = element_text(size=9,  colour = "black"),
		  strip.text.x = element_text(size = 10,  colour = "black" ),
		  strip.text.y = element_text(size = 10,  colour = "black"),
		  legend.title = element_blank(),legend.text = element_text(size=8),
		  legend.position = "bottom")+
	scale_fill_manual(values=c("pink", "red"),labels=c("Founder variants","Known founder variants")) +
	scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1.0))

## Save in tiff format
ggsave("Figure_5_founder_Qc_2stars.tiff",graph_UQc,
	   width =50, height = 40, units = 'mm', dpi=600,bg="white")
