################################
#
# Claudia Moreau 11 juin 2025
#
# Do Figure in article Founder variants: Fig.2
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

## Get enriched variants positions
variants_enriched <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Delete variants with freq greater than 5%
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
variants_enriched <- variants_enriched[!variants_enriched$SNP %in% variants_to_freq,]

####################### Enriched variants - SLSJ ############################
## Open file with the individuals in SLSJ
individual_file_SLSJ <- read.table(paste0('CaG_Imput_VTA_SAG_enriched_in_WQ_chr_1.tfam'))
#### Create table with all variants enriched (compared with gnomAD) - Imput

carriers_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('results/Imput/SLSJ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
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
## Select only variants enriched in SLSJ
carriers_df <- carriers_df[carriers_df[,1] %in% variants_enriched$SNP,]

## Found missing variants enriched to add with WGS
variants_missing_SAG <- variants_enriched[!variants_enriched$SNP %in% carriers_df[,1],1]
NR_var[!is.element(NR_var,variants_missing_SAG) ]
unique(carriers_df[,1])

carriers_wgs_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('WGS/SLSJ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
	table_ind <- as.data.frame(table_ind)
	for (icol in c(1:ncol(table_ind))){
		my_rows <- cbind(dimnames(table_ind)[[2]][icol],table_ind[!is.na(table_ind[icol]),icol])
		carriers_wgs_df <- rbind(carriers_wgs_df,my_rows)
	}
}

## Select only variants enriched in SLSJ
carriers_wgs_df <- carriers_wgs_df[carriers_wgs_df[,1] %in% variants_missing_SAG,]
carriers_wgs_df <- carriers_wgs_df[carriers_wgs_df[,1] %in% variants_enriched$SNP,]

unique_wgs_snps <- unique(carriers_wgs_df[,1])
unique_imput_snps <- unique(carriers_df[,1])
unique_wgs_snps[is.element(unique_wgs_snps,unique_imput_snps)]

## Merge all info between Imput and WGS
carriers_enriched_all_SAG <- rbind(carriers_wgs_df,carriers_df)
## persons who carry 2 alleles for the same disease?
carriers_enriched_all_SAG[duplicated(paste(carriers_enriched_all_SAG$V1,carriers_enriched_all_SAG$V2,sep='_')),]
carriers_enriched_all_SAG[carriers_enriched_all_SAG$V2==11133803,]
## Lets remove these...
carriers_enriched_all_SAG <- carriers_enriched_all_SAG[-which(duplicated(paste(carriers_enriched_all_SAG$V1,carriers_enriched_all_SAG$V2,sep='_'))),]

### table this to find the number of variants for each ind
ind_num_var_df <- as.data.frame(table(carriers_enriched_all_SAG[,2]))
ind_num_var_df[ind_num_var_df$Var1==11133803,]
hist_like_df <- table(ind_num_var_df$Freq)

################################### Enriched variants - UQc ##############################################
## Open file with individuals in UQc
individual_file_RQ <- read.table(paste0('CaG_Imput_VTA_RQ_enriched_in_WQ_chr_1.tfam'))
#### Create table with all variants enriched (compared with gnomAD) - Imput

carriers_rq_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('Imput/RQ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
	table_ind <- as.data.frame(table_ind)
	for (icol in c(1:ncol(table_ind))){
		print (icol)
		my_rows <- cbind(dimnames(table_ind)[[2]][icol],table_ind[!is.na(table_ind[icol]),icol])
		carriers_rq_df <- rbind(carriers_rq_df,my_rows)
	}
}
unique(carriers_rq_df[,1])

## Change : in _ in variant name
carriers_rq_df[,1] <- gsub('\\.', "_", carriers_rq_df[,1])
## Remove NR imput
carriers_rq_df <- carriers_rq_df[!carriers_rq_df[,1] %in% NR_var[,1],]
## Select only variants enriched in SLSJ
carriers_rq_df <- carriers_rq_df[carriers_rq_df[,1] %in% variants_enriched$SNP,]

## Found missing variants enriched to add with WGS
variants_missing <- variants_enriched[!variants_enriched$SNP %in% carriers_rq_df[,1],1]
NR_var[!is.element(NR_var,variants_missing) ]
unique(carriers_rq_df[,1])

## Find WGS for missing
carriers_wgs_df <- data.frame()
for (i in 1:22){
	table_ind <- read.table(paste0('results/WGS/RQ/table_of_individuals_with_VTA_enriched_variants_for_IBD_chr',i,'.txt'), header=T)
	table_ind <- as.data.frame(table_ind)
	for (icol in c(1:ncol(table_ind))){
		my_rows <- cbind(dimnames(table_ind)[[2]][icol],table_ind[!is.na(table_ind[icol]),icol])
		carriers_wgs_df <- rbind(carriers_wgs_df,my_rows)
	}
}


## Select only variants enriched in SLSJ
carriers_wgs_df <- carriers_wgs_df[carriers_wgs_df[,1] %in% variants_missing,]
carriers_wgs_df <- carriers_wgs_df[carriers_wgs_df[,1] %in% variants_enriched$SNP,]

unique_wgs_snps <- unique(carriers_wgs_df[,1])
unique_imput_snps <- unique(carriers_rq_df[,1])
unique_wgs_snps[is.element(unique_wgs_snps,unique_imput_snps)]

## Merge all info between Imput and WGS
carriers_enriched_all_RQ <- rbind(carriers_wgs_df,carriers_rq_df)
## persons who carry 2 alleles for the same disease?
carriers_enriched_all_RQ[duplicated(paste(carriers_enriched_all_RQ$V1,carriers_enriched_all_RQ$V2,sep='_')),]
carriers_enriched_all_RQ[carriers_enriched_all_RQ$V2==11104189,]
## Lets remove these...
carriers_enriched_all_RQ <- carriers_enriched_all_RQ[-which(duplicated(paste(carriers_enriched_all_RQ$V1,carriers_enriched_all_RQ$V2,sep='_'))),]

### table this to find the number of variants for each ind
ind_num_var_rq_df <- as.data.frame(table(carriers_enriched_all_RQ[,2]))
ind_num_var_rq_df[ind_num_var_rq_df$Var1==11104189,]
hist_like_rq_df <- table(ind_num_var_rq_df$Freq)


########## Resampling #################
resampling_vec <- c(rep(0,(21472-length(ind_num_var_rq_df[,1]))),ind_num_var_rq_df[,2])
resampling_df <- data.frame(variants_num=c(0:9))
for (i in c(1:1000)){
	my_sample <- sample(resampling_vec,3589,replace = FALSE)
	my_table_sample <- as.data.frame(table(my_sample))
	#print (my_table_sample)
	resampling_df <- merge(resampling_df,my_table_sample,by=1,all=TRUE)
}
resampling_df[,1:10]
resampling_df[is.na(resampling_df)] <- 0
resampling_df <- resampling_df[,-1]

## Look at means and CI
mean_ci <- apply(resampling_df,MARGIN = 1,CI)

# Table for figure
 
barplot_res_df <- as.data.frame(t(mean_ci))
barplot_res_df$variants_num <- c(0:9)
sum(barplot_res_df$mean)
barplot_res_df$group <- 'UQc'
SLSJ_res <- as.numeric(hist_like_df)
SLSJ_res <- c(3589-sum(SLSJ_res),SLSJ_res)
SLSJ_res_df <- data.frame(upper=NA,mean=SLSJ_res,lower=NA,variants_num=c(0:9),group='SLSJ')
barplot_res_df <- rbind(barplot_res_df,SLSJ_res_df)
RQ_res <- as.numeric(hist_like_rq_df)
RQ_res <- c(21472-sum(RQ_res),RQ_res)
barplot_res_df$value <- c(RQ_res,SLSJ_res)

#### Do the graph

graph_SAG <- ggplot(barplot_res_df, aes(x=variants_num, y=mean/3589,fill=group)) +
	geom_bar(stat='identity', position = position_dodge(width = 0.9)) +
	geom_errorbar(aes(ymin = lower/3589, ymax = upper/3589),position = position_dodge(width = 0.9),width=0.5,linewidth=0.05) +
	labs(x = "Number of variants", y = "Proportion of individuals") + 
	scale_x_continuous("Number of variants", labels = c(0:9),breaks=c(0:9))+
	geom_text(aes(label=value),
			   vjust = -0.5, color="black", size=2,position = position_dodge(width = 0.9)) +
	theme(plot.title = element_text(size=10,  colour= "black" ),
		  axis.title.x = element_text(size=12,  colour = "black"),    
		  axis.title.y = element_text(size=12,  colour = "black"),    
		  axis.text.x = element_text(size=12,  colour = "black"), 
		  axis.text.y = element_text(size=12,  colour = "black"),
		  strip.text.x = element_text(size = 12,  colour = "black" ),
		  strip.text.y = element_text(size = 12,  colour = "black"),
		  legend.title = element_blank(),legend.text = element_text(size=8),
		  legend.position = "bottom")+
	scale_fill_manual(values=c("blue", "green"),labels=c("Variants with RFD>=10% in SLSJ","Variants with RFD>=10% in UQc")) +
	scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.4))

## Save in tiff format
ggsave("Figure_2_enriched.tiff",graph_SAG,
	    width =132, height = 132,units = 'mm', dpi=600,bg="white")

tiff("Figure_2_enriched.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
graph_SAG
dev.off()


########## Resampling for the number of variants not found. ############

resampling_vec <- c()
for (i in c(1:1000)){
	my_sample <- sample(individual_file_RQ[,1],3589,replace = FALSE)
	num_var_present <- carriers_enriched_all_RQ[is.element(carriers_enriched_all_RQ$V2,my_sample),1]
	print (length(num_var_present))
	print (length(unique(num_var_present)))
	resampling_vec <- c(resampling_vec,1302-length(unique(num_var_present)))
}
CI(resampling_vec)

## Perform chi2

# Create a matrix of counts
counts <- matrix(c(318, 1302-318, 540, 1302-540), nrow = 2, byrow = TRUE)
colnames(counts) <- c("Absent", "Present")
rownames(counts) <- c("UQc", "SLSJ")

# Perform chi-squared test
chisq.test(counts)


