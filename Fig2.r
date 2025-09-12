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

## Get all variants present either in the QcP or in gnomAD
variants_enriched <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt", header=T)
variants_enriched <- variants_enriched[,2]

## Remove variants with frequencies greater than 5%
variants_frq_df <- read.table("mymerged_all_chr_Imput_VTA_WQ.frq", header=TRUE)
variants_to_freq <- variants_frq_df$SNP[variants_frq_df$MAF>0.05]
variants_enriched <- variants_enriched[!variants_enriched %in% variants_to_freq]

## Remove unreliable imputed
NR_var <- read.table("ImputClaudia_1729WGSCaG_false_positive_snps.txt")
NR_var$V1 <- gsub('_',':',NR_var$V1)

## SNPs not in WGS
snps_not_in_wgs <- read.table("variants_in_imput_not_in_WGS.txt")

## All snps to be included in fig 2
variants_enriched <- variants_enriched[!variants_enriched %in% NR_var$V1]
variants_included <- variants_enriched[!variants_enriched %in% snps_not_in_wgs$V1]

## Remove variants with one star classification except known variants (for subplot)
toremove <- read.table("1302_snps_clinvar_1star.txt", header = FALSE)
toremove$V1 <- gsub("_", ":", toremove$V1)
known_tokeep <- read.table("known_founder_variants_in4reviews.txt", header = FALSE)
toremove <- toremove[!toremove$V1 %in% known_tokeep$V1,]
toremove <- as.data.frame(toremove)
variants_included_2stars <- variants_included[!variants_included %in% toremove$toremove]

####################### Enriched variants - SLSJ ############################
## Open file with the individuals in SLSJ
individual_file_SLSJ <- read.table(paste0('CaG_Imput_VTA_SAG_enriched_in_WQ_chr_1.tfam'))
## Read SLSJ carrier file (SNP; carrier)
carriers_df <- read.table("carriers_mymerged_all_chr_Imput_VTA_SAG.txt", header = TRUE, sep = ";", check.names = FALSE, stringsAsFactors = FALSE)
carriers_df$variant <- sub("([A-Z])\\d+$", "\\1", carriers_df$variant)
carriers_df[,1] <- gsub('\\.', "_", carriers_df[,1])

## Select only variants included in SLSJ
carriers_df <- carriers_df[carriers_df[,1] %in% variants_included,]
carriers_df_2stars <- carriers_df[carriers_df[,1] %in% variants_included_2stars,]

carriers_enriched_all_SAG <- carriers_df
carriers_enriched_all_SAG_2stars <- carriers_df_2stars

## Remove duplicates (homozygotes)
carriers_enriched_all_SAG <- carriers_enriched_all_SAG[!duplicated(carriers_enriched_all_SAG[c("variant", "IID")]), ]
carriers_enriched_all_SAG_2stars <- carriers_enriched_all_SAG_2stars[!duplicated(carriers_enriched_all_SAG_2stars[c("variant", "IID")]), ]
## How many unique snps
length(unique(carriers_enriched_all_SAG$variant))
length(unique(carriers_enriched_all_SAG_2stars$variant))

## Find the number of variants for each individual
ind_num_var_df <- as.data.frame(table(carriers_enriched_all_SAG[,2]))
hist_like_df <- table(ind_num_var_df$Freq)
ind_num_var_df_2stars <- as.data.frame(table(carriers_enriched_all_SAG_2stars[,2]))
hist_like_df_2stars <- table(ind_num_var_df_2stars$Freq)

################################### Enriched variants - UQc ##############################################
## Open file with individuals in UQc
individual_file_RQ <- read.table(paste0('CaG_Imput_VTA_RQ_enriched_in_WQ_chr_1.tfam'))
## Read UQc carrier file (SNP; carrier)
carriers_rq_df <- read.table("carriers_mymerged_all_chr_Imput_VTA_RQ.txt", header = TRUE, sep = ";")
carriers_rq_df$variant <- sub("([A-Z])\\d+$", "\\1", carriers_rq_df$variant)
carriers_rq_df[,1] <- gsub('\\.', "_", carriers_rq_df[,1])

## Remove NR imput and not in WGS
carriers_rq_df <- carriers_rq_df[carriers_rq_df[,1] %in% variants_included,]
carriers_rq_df_2stars <- subset(carriers_rq_df, !(variant %in% toremove$toremove))

## Merge all info between Imput and WGS
carriers_enriched_all_RQ <- carriers_rq_df
carriers_enriched_all_RQ_2stars <- carriers_rq_df_2stars
## Remove duplicates (homozygotes)
carriers_enriched_all_RQ <- carriers_enriched_all_RQ[!duplicated(carriers_enriched_all_RQ[c("variant", "IID")]), ]
carriers_enriched_all_RQ_2stars <- carriers_enriched_all_RQ_2stars[!duplicated(carriers_enriched_all_RQ_2stars[c("variant", "IID")]), ]
## Number of unique variants
length(unique(carriers_enriched_all_RQ$variant))
length(unique(carriers_enriched_all_RQ_2stars$variant))

## Find the number of variants for each individual
ind_num_var_rq_df <- as.data.frame(table(carriers_enriched_all_RQ[,2]))
hist_like_rq_df <- table(ind_num_var_rq_df$Freq)
ind_num_var_rq_df_2stars <- as.data.frame(table(carriers_enriched_all_RQ_2stars[,2]))
hist_like_rq_df_2stars <- table(ind_num_var_rq_df_2stars$Freq)

########## Resampling #################
## All variants
resampling_vec <- c(rep(0,(21472-length(ind_num_var_rq_df[,1]))),ind_num_var_rq_df[,2])
resampling_df <- data.frame(variants_num=c(0:9))
for (i in c(1:1000)){
	my_sample <- sample(resampling_vec,3589,replace = FALSE)
	my_table_sample <- as.data.frame(table(my_sample))
	#print (my_table_sample)
	resampling_df <- merge(resampling_df,my_table_sample,by=1,all=TRUE)
}
resampling_df[is.na(resampling_df)] <- 0
resampling_df <- resampling_df[,-1]
## Look at means and CI
mean_ci <- apply(resampling_df,MARGIN = 1,CI)
### 2 stars
resampling_vec_2stars <- c(rep(0,(21472-length(ind_num_var_rq_df_2stars[,1]))),ind_num_var_rq_df_2stars[,2])
resampling_df_2stars <- data.frame(variants_num=c(0:9))
for (i in c(1:1000)){
  my_sample <- sample(resampling_vec_2stars,3589,replace = FALSE)
  my_table_sample <- as.data.frame(table(my_sample))
  #print (my_table_sample)
  resampling_df_2stars <- merge(resampling_df_2stars,my_table_sample,by=1,all=TRUE)
}
resampling_df_2stars[,1:10]
resampling_df_2stars[is.na(resampling_df_2stars)] <- 0
resampling_df_2stars <- resampling_df_2stars[,-1]
## Look at means and CI
mean_ci_2stars <- apply(resampling_df_2stars,MARGIN = 1,CI)

##### Prepare data with all variants for the plot
# Table for figure
barplot_res_df <- as.data.frame(t(mean_ci))
barplot_res_df$variants_num <- c(0:9)
sum(barplot_res_df$mean)
barplot_res_df$group <- 'UQc'
SLSJ_res_df <- data.frame(variants_num = 0:9)
SLSJ_res <- setNames(rep(0, 10), 0:9)
SLSJ_res[names(hist_like_df)] <- as.numeric(hist_like_df)
# Fill in count for 0 variants (total is nind from tfam)
SLSJ_res["0"] <- 3589 - sum(SLSJ_res[names(SLSJ_res) != "0"])
SLSJ_res_df$mean  <- SLSJ_res
SLSJ_res_df$upper <- NA
SLSJ_res_df$lower <- NA
SLSJ_res_df$group <- "SLSJ"
barplot_res_df <- rbind(barplot_res_df,SLSJ_res_df)
RQ_res <- setNames(rep(0, 10), 0:9)
RQ_res[names(hist_like_rq_df)] <- as.numeric(hist_like_rq_df)
# Fill in count for 0 variants (total is nind from tfam)
RQ_res["0"] <- 21472 - sum(RQ_res[names(RQ_res) != "0"])
barplot_res_df$value <- c(RQ_res, SLSJ_res)
##### Prepare data with > 1* variants for the facet plot
# Table for figure
barplot_res_df_2stars <- as.data.frame(t(mean_ci_2stars))
barplot_res_df_2stars$variants_num <- c(0:9)
sum(barplot_res_df_2stars$mean)
barplot_res_df_2stars$group <- 'UQc'
SLSJ_res_df_2stars <- data.frame(variants_num = 0:9)
SLSJ_res_2stars <- setNames(rep(0, 10), 0:9)
SLSJ_res_2stars[names(hist_like_df_2stars)] <- as.numeric(hist_like_df_2stars)
# Fill in count for 0 variants (total is nind from tfam)
SLSJ_res_2stars["0"] <- 3589 - sum(SLSJ_res_2stars[names(SLSJ_res_2stars) != "0"])
SLSJ_res_df_2stars$mean  <- SLSJ_res_2stars
SLSJ_res_df_2stars$upper <- NA
SLSJ_res_df_2stars$lower <- NA
SLSJ_res_df_2stars$group <- "SLSJ"
barplot_res_df_2stars <- rbind(barplot_res_df_2stars,SLSJ_res_df_2stars)
RQ_res_2stars <- setNames(rep(0, 10), 0:9)
RQ_res_2stars[names(hist_like_rq_df_2stars)] <- as.numeric(hist_like_rq_df_2stars)
# Fill in count for 0 variants (total is nind from tfam)
RQ_res_2stars["0"] <- 21472 - sum(RQ_res_2stars[names(RQ_res_2stars) != "0"])
barplot_res_df_2stars$value <- c(RQ_res_2stars, SLSJ_res_2stars)

#### Do the graph
## All variants
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
## 2*
graph_SAG_2stars <- ggplot(barplot_res_df_2stars, aes(x=variants_num, y=mean/3589,fill=group)) +
  geom_bar(stat='identity', position = position_dodge(width = 0.9), show.legend = FALSE) +
  geom_errorbar(aes(ymin = lower/3589, ymax = upper/3589),position = position_dodge(width = 0.9),width=0.5,linewidth=0.05) +
  labs(x = "", y = "", title = "> 1* ClinVar review status") +
  scale_x_continuous("Number of variants", labels = c(0:9),breaks=c(0:9))+
  theme(plot.title = element_text(size=9,  colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=9,  colour = "black"),
        axis.text.y = element_text(size=9,  colour = "black"),
        strip.text.x = element_text(size = 12,  colour = "black" ),
        strip.text.y = element_text(size = 12,  colour = "black"),
        legend.title = element_blank(),legend.text = element_text(size=8),
        legend.position = "bottom")+
  scale_fill_manual(values=c("blue", "green"),labels=c("Variants with RFD>=10% in SLSJ","Variants with RFD>=10% in UQc")) +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.4))
tiff("Figure_2_enriched_2stars_only.tiff",width =66, height = 66,units = 'mm', res=600 ,compression = "lzw")
graph_SAG_2stars
dev.off()

########## Resampling for the number of variants not found. ############
## for UQc
total_variants <- length(unique(variants_included))
resampling_vec <- c()
for (i in c(1:1000)){
	my_sample <- sample(individual_file_RQ[,1],3589,replace = FALSE)
	num_var_present <- carriers_enriched_all_RQ[is.element(carriers_enriched_all_RQ$V2,my_sample),1]
	print (length(num_var_present))
	print (length(unique(num_var_present)))
	resampling_vec <- c(resampling_vec,total_variants-length(unique(num_var_present)))
}
CI(resampling_vec)
## for SLSJ 
my_sample <- individual_file_SLSJ[,1]
num_var_present <- carriers_enriched_all_SAG[is.element(carriers_enriched_all_SAG$IID, my_sample), 1]
result <- total_variants - length(unique(num_var_present))
print(length(num_var_present))
print(length(unique(num_var_present)))
total_variants-length(unique(num_var_present))

## Perform chi2
# Create a matrix of counts
counts <- matrix(c(401, total_variants-401, 682, total_variants-682), nrow = 2, byrow = TRUE)
colnames(counts) <- c("Absent", "Present")
rownames(counts) <- c("UQc", "SLSJ")

# Perform chi-squared test
chisq.test(counts)

##### chi2 and resampling for individuals with at least 2 variants
resampling_vec_diff <- c()
num_indiv_more2 <- c()
num_indiv_less2 <- c()

for (i in 1:1000) {
  my_sample <- sample(individual_file_RQ[,1], 3589, replace = FALSE)
  variants_sample <- carriers_enriched_all_RQ[carriers_enriched_all_RQ$IID %in% my_sample, ]
  counts_per_indiv <- table(variants_sample$IID)
  nind_2vars <- sum(counts_per_indiv >= 2)
  nind_less <- length(my_sample) - nind_2vars
  num_indiv_more2 <- c(num_indiv_more2, nind_2vars)
  num_indiv_less2 <- c(num_indiv_less2, nind_less)
  resampling_vec_diff <- c(resampling_vec_diff, nind_2vars - nind_less)
}
mean_more <- mean(num_indiv_more2)
mean_less <- mean(num_indiv_less2)

## Same for slsj 
counts_per_indiv <- table(carriers_enriched_all_SAG$IID)
nind_2vars_slsj <- sum(counts_per_indiv >= 2)
n_lt2_slsj <- 3589 - nind_2vars_slsj

## Perform chi2
# Create a matrix of counts
counts <- matrix(c(mean_less, mean_more, n_lt2_slsj, nind_2vars_slsj), nrow = 2, byrow = TRUE)
colnames(counts) <- c("Less", "More")
rownames(counts) <- c("UQc", "SLSJ")
counts
# Perform chi-squared test
chisq.test(counts)

