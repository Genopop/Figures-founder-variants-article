################################
#
# September 2024 by ?lisa Michel
#
# Modified by Claudia Moreau 2025
#
# Do Figures in article
################################

### Packages
rm(list = ls());
library('ggplot2')
library("ggpubr")
library("gridExtra")
library("ggpmisc")
#library('ggstats')
library("stringr")
library("stringi")
library("tidyr")
library(grid)
#library("openxlsx")
#library('plotly')
library(htmlwidgets)
library(patchwork)
library(cowplot)



options(scipen=999)
## File to open with all variants enriched (1304)
df_final <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Open file of variants known
variants_known <- read.table("variants_reported_in_SLSJ.txt")
#variants_new_whose_gene_reported <- read.table("/lustre03/project/6033529/saguenay_disease/Analyses_final_version_really/Article/variants_new_with_genes_reported_in_SLSJ.txt")

################################################################# Fig.3 ############################################################
######### Founder variants in SLSJ or UQc, graphs with CR UQc vs CR SLSJ
df_final_founder_RQ_SAG <- df_final
df_final_founder_RQ_SAG$Founder_SAG <- ifelse(is.na(df_final$Ind_Imput_SAG), 
											  df_final$founder_WGS_SAG,df_final$founder_Imput_SAG)
df_final_founder_RQ_SAG$Founder_RQ <- ifelse(is.na(df_final$Ind_Imput_RQ), 
											 df_final$founder_WGS_RQ,df_final$founder_Imput_RQ)
df_final_founder_RQ_SAG$Founder_WQ <- ifelse(is.na(df_final$Ind_Imput_WQ), 
											 df_final$founder_WGS_WQ,df_final$founder_Imput_WQ)
## Select founder variants SLSJ or UQc or QcP
df_final_founder_RQ_SAG <- df_final_founder_RQ_SAG[df_final_founder_RQ_SAG$Founder_SAG == "Founder" | df_final_founder_RQ_SAG$Founder_RQ == "Founder" | df_final_founder_RQ_SAG$Founder_WQ == "Founder",]
df_final_founder_RQ_SAG <- df_final_founder_RQ_SAG[complete.cases(df_final_founder_RQ_SAG$SNP),]
## Select columns
df_final_founder_RQ_SAG_clean <- df_final_founder_RQ_SAG[,c(1:23)]
df_final_founder_RQ_SAG_clean <- df_final_founder_RQ_SAG_clean[,-c(4,5,7,8,12,13,15,16,21,22)]

## Select carrier rate for SLSJ and UQc (wether Imput or WGS if missing)
df_final_founder_RQ_SAG_clean$CR_SAG <- ifelse(is.na(df_final_founder_RQ_SAG_clean$CR_Imput_SAG), 
											   df_final_founder_RQ_SAG_clean$CR_WGS_SAG,df_final_founder_RQ_SAG_clean$CR_Imput_SAG)

df_final_founder_RQ_SAG_clean$CR_RQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$CR_Imput_RQ), 
											  df_final_founder_RQ_SAG_clean$CR_WGS_RQ,df_final_founder_RQ_SAG_clean$CR_Imput_RQ)
## Select Ind that correspond in WQ
df_final_founder_RQ_SAG_clean$Ind_WQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$Ind_Imput_WQ), 
											   df_final_founder_RQ_SAG_clean$Ind_WGS_WQ,df_final_founder_RQ_SAG_clean$Ind_Imput_WQ)
df_final_founder_RQ_SAG_clean$Source_ind_WQ <- ifelse(is.na(df_final_founder_RQ_SAG_clean$Ind_Imput_WQ), 
													  "WGS","Imputed data")

## Say if founder or not 
df_final_founder_RQ_SAG_clean$Founder_SAG <- df_final_founder_RQ_SAG$Founder_SAG
df_final_founder_RQ_SAG_clean$Founder_RQ <- df_final_founder_RQ_SAG$Founder_RQ
## Select rows for graphs
df_for_graphs_RQ_SAG <- df_final_founder_RQ_SAG_clean[,c(1,2,3,14,15,16,17,18,19)]

## Fill all NA values for peak
df_for_graphs_RQ_SAG[is.na(df_for_graphs_RQ_SAG$CR_SAG),4] <- 0
df_for_graphs_RQ_SAG[is.na(df_for_graphs_RQ_SAG$CR_RQ),5] <- 0
## Set missing as not founder
df_for_graphs_RQ_SAG$Founder_RQ = ifelse(is.na(df_for_graphs_RQ_SAG$Founder_RQ), "Not_founder",df_for_graphs_RQ_SAG$Founder_RQ)
df_for_graphs_RQ_SAG$Founder_SAG = ifelse(is.na(df_for_graphs_RQ_SAG$Founder_SAG), "Not_founder",df_for_graphs_RQ_SAG$Founder_SAG)
## Set founder both or none or only one
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Founder'] <- 'Both'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Not_founder'] <- 'Founder_RQ'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Not_founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Founder'] <- 'Founder_SAG'
df_for_graphs_RQ_SAG$Founder[df_for_graphs_RQ_SAG$Founder_RQ == 'Not_founder' & df_for_graphs_RQ_SAG$Founder_SAG == 'Not_founder'] <- 'None'

## With variants known or not
df_for_graphs_RQ_SAG$Status <-ifelse(df_for_graphs_RQ_SAG$ID %in% variants_known[,1], "Known founder variant",
									 "New founder variant")

founder_WQ_known <- ggplot(data = df_for_graphs_RQ_SAG, aes(x = CR_SAG, y = CR_RQ,colour =Status,
							shape=Source_ind_WQ,text=SNP),size=2,position=position_jitter(h=0.00,w=0.15))+
	geom_point(size=0.5) + 
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
	scale_shape_manual(values=c(19,0))+
	theme(axis.text.x = element_text( size = 6),axis.title=element_text(size=12),
		  axis.text.y = element_text( size = 6),legend.text = element_text(size=12)) + 
	scale_color_manual(values = c("red","blue"))+
	theme(legend.title = element_blank(),legend.position = "bottom",legend.box = "vertical")+ guides(color=guide_legend(override.aes = list(size=2)))+
	xlab("Carrier rate in SLSJ") + ylab("Carrier rate in UQc")+
	scale_x_continuous("Carrier rate in SLSJ",breaks=c(0,100,200,300),limits=c(0,300),labels = c(0,"1/100","1/200","1/300")) +
	scale_y_continuous("Carrier rate in UQc",breaks=c(0,2000,4000,6000,8000),limits=c(0,8000),
					   labels = c(0,"1/2,000","1/4,000","1/6,000","1/8,000"))

## Save tiff format
ggsave("Figure_3.tiff",plot = founder_WQ_known, width = 132, height = 132,dpi=600,units = "mm")
tiff("Figure_3.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
founder_WQ_known
dev.off()

################################################################# Fig.1 ############################################################
######### Graph Gnomad/(Imput if missing WGS) MAF on enriched Variants (1200) in QcP
## Variants with freq of Gnomad
NTNFE_WGS <- read.table("final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_WGS_CaG_all_chr_tresh_0.1.txt", header =T)
NTNFE_WGS$SNP <- paste0("chr",NTNFE_WGS$CHROM,"_",NTNFE_WGS$POS,"_",NTNFE_WGS$REF,"_",NTNFE_WGS$ALT)
NTNFE_WGS <- NTNFE_WGS[,c(5,18,13,15)]
NTNFE_imput <- read.table("final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_Imput_CaG_tresh_0.1.txt", header =T)
NTNFE_imput <- NTNFE_imput[,c(1,13,15)]
## Create table
mat = matrix(ncol = 0, nrow = 1304)
df_gnomad=data.frame(mat)
df_gnomad$SNP <- df_final$SNP
df_gnomad <- merge(df_gnomad,NTNFE_imput,by='SNP', all.x=T, all.y=F)
df_gnomad <- merge(df_gnomad,NTNFE_WGS,by='SNP', all.x=T, all.y=F)
df_gnomad <- df_gnomad[,-2]
names(df_gnomad) <- c("SNP","MAF_WQ_Imput","ID","MAF_GNOMAD","MAF_WQ_WGS")
## Get MAF missing
df_gnomad$MAF_WQ <- ifelse(is.na(df_gnomad$MAF_WQ_Imput),df_gnomad$MAF_WQ_WGS,df_gnomad$MAF_WQ_Imput)

## Variants founder in QcP WGS and Imput
variants_WGS_WQ_founder <- read.table('WGS/Table_founder_variants_in_WQ.txt',header=T)
variants_Imput_WQ_founder <- read.table('Imput/Table_founder_variants_in_WQ.txt',header=T)

## Get the status of the variant
df_gnomad$founder_Imput_WQ <- df_final$founder_Imput_WQ
df_gnomad$founder_WGS_WQ <- df_final$founder_WGS_WQ

## Get founder of either Imput and if missing WGS
df_gnomad$Founder_status<-  ifelse(is.na(df_gnomad$founder_Imput_WQ),df_gnomad$founder_WGS_WQ,df_gnomad$founder_Imput_WQ )
## Replace not_founder by Not founder
df_gnomad$Founder_status<-  gsub("Not_founder", "Non founder with relative difference ≥ 0.1", df_gnomad$Founder_status)
df_gnomad$Founder_status<-  ifelse(is.na(df_gnomad$Founder_status),"Non founder with relative difference ≥ 0.1",df_gnomad$Founder_status)
## Get if known or new
df_gnomad$Founder_status<-  ifelse((df_gnomad$ID %in% variants_Imput_WQ_founder$ID | (!(df_gnomad$ID %in%variants_Imput_WQ_founder) & df_gnomad$ID %in% variants_WGS_WQ_founder$ID)),
								   "New founder variant",df_gnomad$Founder_status)
df_gnomad$Founder_status<-  ifelse((df_gnomad$ID %in% variants_known[,1] & df_gnomad$Founder_status == "New founder variant"),
								   "Known founder variant",df_gnomad$Founder_status)

## Do figure
gnomad_vs_founder_enriched_only <- ggplot(data = df_gnomad,aes(x = MAF_WQ, y = MAF_GNOMAD)) +
	geom_point(size=0.1)+
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
	scale_shape_manual(values=c(1,5))+
	theme(axis.text.x = element_text( size = 12),axis.title=element_text(size=12),
		  axis.text.y = element_text(size = 12),
		  plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm")) +
	xlab("Frequency in the QcP") + ylab("Frequency in gnomAD nfe")  +xlim(0,0.03)+ylim(0,0.02) 


######### Graph Gnomad/(Imput if missing WGS) MAF on enriched Variants (1304) in SLSJ 
## Variants with freq of Gnomad
## Get MAF for SLSJ in Imput ad WGS
NTNFE_WGS_SAG <- read.table("SLSJ/mymerged_all_chr_CaG_seq_VTA_SAG_enriched.frq", header =T)
NTNFE_WGS_SAG <- NTNFE_WGS_SAG[,c(2,5)]
NTNFE_imput_SAG <- read.table("SLSJ/mymerged_all_chr_Imput_VTA_SAG_enriched.frq", header =T)
NTNFE_imput_SAG <- NTNFE_imput_SAG[,c(2,5)]
NTNFE_imput_SAG$SNP <- gsub(":", "_", NTNFE_imput_SAG$SNP)
## Get freq in gnomad
Gnomad_freq <- read.table("final_table_VTA_variant_enriched_in_WQ_compared_NTNFE_WGS_CaG_all_chr_tresh_0.1.txt", header =T)
Gnomad_freq$SNP <- paste0("chr",Gnomad_freq$CHROM,"_",Gnomad_freq$POS,"_",Gnomad_freq$REF,"_",Gnomad_freq$ALT)
Gnomad_freq <- Gnomad_freq[,c(5,13,18)]

## Create table
df_gnomad_SAG=data.frame(mat)
df_gnomad_SAG$SNP <- df_final$SNP
df_gnomad_SAG <- merge(df_gnomad_SAG,NTNFE_imput_SAG,by='SNP', all.x=T, all.y=F)
df_gnomad_SAG <- merge(df_gnomad_SAG,NTNFE_WGS_SAG,by='SNP', all.x=T, all.y=F)
df_gnomad_SAG <- merge(df_gnomad_SAG,Gnomad_freq,by='SNP', all.x=T, all.y=F)
names(df_gnomad_SAG) <- c("SNP","MAF_SAG_Imput","MAF_SAG_WGS","ID","MAF_GNOMAD")

## Select founder variants in SLSJ (Imput and WGS if missing)
## Variants founder in SLSJ WGS and Imput
variants_WGS_SAG_founder <- read.table('WGS/Table_founder_variants_in_SLSJ.txt',header=T)
variants_Imput_SAG_founder <- read.table('Imput/Table_founder_variants_in_SLSJ.txt',header=T)

## Get the status of the variant
df_gnomad_SAG$founder_Imput_SAG <- df_final$founder_Imput_SAG
df_gnomad_SAG$founder_WGS_SAG <- df_final$founder_WGS_SAG

## Get founder of either Imput and if missing WGS
df_gnomad_SAG$Founder_status<-  ifelse(is.na(df_gnomad_SAG$founder_Imput_SAG),df_gnomad_SAG$founder_WGS_SAG,df_gnomad_SAG$founder_Imput_SAG)
## Replace not_founder by Not founder
df_gnomad_SAG$Founder_status<-  gsub("Not_founder", "Non founder with relative difference ≥ 0.1", df_gnomad_SAG$Founder_status)
df_gnomad_SAG$Founder_status<-  ifelse(is.na(df_gnomad_SAG$Founder_status),"Non founder with relative difference ≥ 0.1",df_gnomad_SAG$Founder_status)
## Get if known or new
df_gnomad_SAG$Founder_status<-  ifelse((df_gnomad_SAG$ID %in% variants_Imput_SAG_founder$ID | (!(df_gnomad_SAG$ID %in%variants_Imput_SAG_founder) & df_gnomad_SAG$ID %in% variants_WGS_SAG_founder$ID)),
									   "New founder variant",df_gnomad_SAG$Founder_status)
df_gnomad_SAG$Founder_status<-  ifelse((df_gnomad_SAG$ID %in% variants_known[,1]& df_gnomad_SAG$Founder_status == "New founder variant"), 
									   "Known founder variant",df_gnomad_SAG$Founder_status)

## Get MAF (Imput or WGS if missing)
df_gnomad_SAG$MAF_SAG <- ifelse(is.na(df_gnomad_SAG$MAF_SAG_Imput),df_gnomad_SAG$MAF_SAG_WGS,df_gnomad_SAG$MAF_SAG_Imput)


## Do figure
gnomad_vs_founder_SAG_enriched_only <- ggplot(data = df_gnomad_SAG,aes(x = MAF_SAG, y = MAF_GNOMAD)) + 
	geom_point(size=0.1)+
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
	scale_shape_manual(values=c(1,5))+
	theme(axis.text.x = element_text(size = 12),axis.title=element_text(size=12),
		  axis.text.y = element_blank(),
		  plot.margin = unit(c(0.8, 0.2, 0, 0.2), "cm")) +
	theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=3)))+
	rremove("ylab") +
	# scale_x_log10() + scale_y_log10() +
	xlab("Frequency in SLSJ") + ylab("Frequency in gnomAD nfe")  +xlim(0,0.03)+ylim(0,0.02) 


##### -- FIRST FIGURE --
figure_1_freq_Gnomad <- ggarrange(gnomad_vs_founder_enriched_only, gnomad_vs_founder_SAG_enriched_only, ncol=2, common.legend = TRUE, legend="bottom",labels = c("A", "B"),label.y = 1,widths = c(1.2,1))## Save tiff format
ggsave("Figure_1.tiff",plot=figure_1_freq_Gnomad, width = 132, height = 100,dpi=600,units = "mm", bg="white")

tiff("Figure_1.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
figure_1_freq_Gnomad
dev.off()

################################################################# Fig.4 ############################################################
##### Figure Simon with CR > 1/1000 for SLSJ
## Select row I want

CR_for_graph <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11,26,27)]
## Remove variants to freq
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph <- CR_for_graph[!CR_for_graph$SNP %in% variants_to_freq,]
## Select founder in SLSJ and UQc
CR_for_graph$Founder_SAG <- ifelse(is.na(CR_for_graph$CR_Imput_SAG),
								   CR_for_graph$founder_WGS_SAG,CR_for_graph$founder_Imput_SAG)
CR_for_graph$Founder_RQ <- ifelse(is.na(CR_for_graph$CR_Imput_RQ),
								  CR_for_graph$founder_WGS_RQ,CR_for_graph$founder_Imput_RQ)
CR_for_graph$Founder_WQ <- ifelse(is.na(CR_for_graph$CR_Imput_WQ),
								  CR_for_graph$founder_WGS_WQ,CR_for_graph$founder_Imput_WQ)

CR_for_graph <- CR_for_graph[CR_for_graph$Founder_SAG == "Founder" | CR_for_graph$Founder_RQ == "Founder"| CR_for_graph$Founder_WQ == "Founder",]

## Get CR (either imput or if missing)
CR_for_graph$CR_SAG <- ifelse(is.na(CR_for_graph$CR_Imput_SAG),CR_for_graph$CR_WGS_SAG,CR_for_graph$CR_Imput_SAG)
CR_for_graph$CR_RQ <- ifelse(is.na(CR_for_graph$CR_Imput_RQ),CR_for_graph$CR_WGS_RQ,CR_for_graph$CR_Imput_RQ)

## Filter CR > 1/1000 in SLSJ
CR_for_graph <- CR_for_graph[complete.cases(CR_for_graph$SNP),]

CR_for_graph[is.na(CR_for_graph)] <- -1

# Create bins with intervals of 20 for UQc
## Claudia modified 11 juin 2025 - until 1/1000
bins_rq <- seq(0, 1000, by = 100)
# Cut the RQ data into bins
CR_for_graph$carrier_rate_intervals <- cut(CR_for_graph$CR_RQ, breaks = bins_rq, right = FALSE)
# Count the frequency of each interval for UQc
frequency_rq <- table(CR_for_graph$carrier_rate_intervals)

# Convert frequency table to dataframe for UQc
frequency_df_rq <- as.data.frame(frequency_rq)
frequency_df_rq$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_rq)))
frequency_df_rq$dataset <- "RQ"

# Create bins with intervals of 20 for SLSJ
bins_sag <- seq(0, 1000, by = 100)

# Cut the SLSJ data into bins
CR_for_graph$carrier_rate_intervals <- cut(CR_for_graph$CR_SAG, breaks = bins_sag, right = FALSE)

# Count the frequency of each interval for SLSJ
frequency_sag <- table(CR_for_graph$carrier_rate_intervals)

# Convert frequency table to dataframe for SLSJ
frequency_df_sag <- as.data.frame(frequency_sag)
frequency_df_sag$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_sag)))
frequency_df_sag$dataset <- "SAG"

# Combine the frequency dataframes for UQc and SLSJ
combined_frequency <- rbind(frequency_df_rq, frequency_df_sag)

########################## Enriched 
## Select row I want
CR_for_graph_enriched <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11)]
## Remove variants to freq
variants_to_freq <- c("chr6_26090951_C_G","chr5_126577145_T_C","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph_enriched <- CR_for_graph_enriched[!CR_for_graph_enriched$SNP %in% variants_to_freq,]

## Get CR (either imput or if missing)
CR_for_graph_enriched$CR_SAG <- ifelse(is.na(CR_for_graph_enriched$CR_Imput_SAG),CR_for_graph_enriched$CR_WGS_SAG,CR_for_graph_enriched$CR_Imput_SAG)
CR_for_graph_enriched$CR_RQ <- ifelse(is.na(CR_for_graph_enriched$CR_Imput_RQ),CR_for_graph_enriched$CR_WGS_RQ,CR_for_graph_enriched$CR_Imput_RQ)

## Filter CR > 1/1000 in SLSJ
CR_for_graph_enriched <- CR_for_graph_enriched[complete.cases(CR_for_graph_enriched$SNP),]

CR_for_graph_enriched[is.na(CR_for_graph_enriched)] <- -1

# Create bins with intervals of 20 for UQc
bins_rq_enriched <- seq(0, 1000, by = 100)
# Cut the RQ data into bins
CR_for_graph_enriched$carrier_rate_intervals <- cut(CR_for_graph_enriched$CR_RQ, breaks = bins_rq_enriched, right = FALSE)
# Count the frequency of each interval for UQc
frequency_rq_enriched <- table(CR_for_graph_enriched$carrier_rate_intervals)

# Convert frequency table to dataframe for UQc
frequency_df_rq_enriched <- as.data.frame(frequency_rq_enriched)
frequency_df_rq_enriched$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_rq_enriched)))
frequency_df_rq_enriched$dataset <- "RQ_enriched"

# Create bins with intervals of 20 for SLSJ
bins_sag_enriched <- seq(0, 1000, by = 100)

# Cut the SLSJ data into bins
CR_for_graph_enriched$carrier_rate_intervals <- cut(CR_for_graph_enriched$CR_SAG, breaks = bins_sag_enriched, right = FALSE)

# Count the frequency of each interval for SLSJ
frequency_sag_enriched <- table(CR_for_graph_enriched$carrier_rate_intervals)

# Convert frequency table to dataframe for SLSJ
frequency_df_sag_enriched <- as.data.frame(frequency_sag_enriched)
frequency_df_sag_enriched$interval <- as.numeric(gsub("\\[|,.*", "", rownames(frequency_df_sag_enriched)))
frequency_df_sag_enriched$dataset <- "SAG_enriched"

# Combine the frequency dataframes for UQc and SLSJ
combined_frequency_enriched <- rbind(frequency_df_rq_enriched, frequency_df_sag_enriched)

## Combine the 2 plots
combined_freq_all <- rbind(combined_frequency,combined_frequency_enriched)
combined_freq_all$Dots <- ifelse(stri_detect_fixed(combined_freq_all$dataset,"enriched"),"variants with RFD≥10%","Founder variants")
combined_freq_all$Region <- ifelse(stri_detect_fixed(combined_freq_all$dataset,"RQ"),"UQc","SLSJ")
## Added by Claudia for barplot
combined_freq_all$new_freq <- combined_freq_all$Freq
combined_freq_all$new_freq[combined_freq_all$Region=='UQc'] <- -combined_freq_all$new_freq[combined_freq_all$Region=='UQc']
combined_freq_all$Var1 <- gsub (pattern = ')',replacement = '[',x = combined_freq_all$Var1)

############## Histogram

plot_final_hist <- ggplot(combined_freq_all[order(-combined_freq_all$Freq),], aes(x = factor(Var1), y = new_freq, fill = interaction(Region, Dots),group=Region)) +
	geom_bar(stat = "identity", color = "black",  position = "identity") +   # Dodge pour séparer SLSJ & UQc
	scale_fill_manual(labels = c("SLSJ founder variants","UQc founder variants",
								 "SLSJ variants with RFD>=10%","UQc variants with RFD>=10%"),
					  values = c("blue","green","lightblue","darkseagreen1")) +
	theme_minimal() +
	scale_y_continuous( breaks=seq(-100,100,10),labels = abs(seq(-100,100,10))) +
	scale_x_discrete( labels = c("1/14-1/99","1/100-1/199","1/200-1/299","1/300-1/399","1/400-1/499","1/500-1/599","1/600-1/699","1/700-1/799","1/800-1/899","1/900-1/999")) +
	theme(axis.text.x = element_text(angle = 90, size = 10), 
		  axis.title=element_text(size=12),
		  axis.text.y = element_text( size = 12),
		  legend.position = "none") +
			xlab("Carrier rate") + 
			ylab("Variants' count") 

## costum legend 
# Custom legend grob
legend_grob <- grobTree(
	textGrob("Founder", x = 0.54, y = 0.95, hjust = 1),
	textGrob("SLSJ", x = 0.45, y = 0.8, hjust = 1),
	pointsGrob(x = 0.5, y = 0.8, pch = 22, gp = gpar(fill = "blue", cex = 1)),
	pointsGrob(x = 0.6, y = 0.8, pch = 22, gp = gpar(fill = "lightblue", cex = 1)),
	textGrob("RFD>=10%", x = 0.74, y =0.95, hjust = 1),
	textGrob("UQc", x = 0.45, y = 0.6, hjust = 1),
	pointsGrob(x = 0.5, y = 0.6, pch = 22, gp = gpar(fill = "green", cex = 1)),
	pointsGrob(x = 0.6, y = 0.6, pch = 22, gp = gpar(fill = "darkseagreen1", cex = 1))
	
)

# Convert grob to ggplot object
legend_plot <- ggdraw() + draw_grob(legend_grob)

plot_final <- plot_grid(plot_final_hist, legend_plot, nrow = 2,rel_heights = c(3, 1))

## Save in tiff format for plos gen
# ggsave("Figure_4_CM.tiff",plot=plot_final, width = 132, height = 132,bg="white",dpi=600, units="mm")
tiff("Figure_4.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
plot_final
dev.off()


################################### Supp Fig.1 ##############################################
df_final_CR <- df_final
## open files
References <- read.table("CR_reported_paper.txt", header=T)
References <- References[-c(35,6),]
CR_cumulated <- read.table("CR_cumulated_all_gene.txt", header=T)
df_final_known_founder <-read.xlsx("Table_known_variants.xlsx")
df_final_known_founder$Position.GRCh38 <- gsub(":","_",df_final_known_founder$Position.GRCh38)

## Add info of references
## Change format References
References$Carrier_rate_reported_SLSJ <- as.numeric(str_split_fixed(References$Carrier_rate_reported_SLSJ,"",3)[,3])
df_final_CR_ref <- merge(df_final_CR,References,by="ID",all=F)
## Select known
df_final_CR_ref <- df_final_CR_ref[df_final_CR_ref$SNP %in% df_final_known_founder$Position.GRCh38, ]
## Add info cumulated CR
df_final_CR_ref <- merge(df_final_CR_ref,CR_cumulated,by="ID",all=T)

## Select CR data
df_final_CR_ref$CR_SLSL_data <- ifelse(is.na(df_final_CR_ref$Ind_Imput_SAG), 
									   df_final_CR_ref$CR_WGS_SAG,df_final_CR_ref$CR_Imput_SAG)

## Select CR cumulated first, if missing CR data
df_final_CR_ref$CR_SLSL_final <- ifelse(is.na(df_final_CR_ref$CR_cumulated_SLSJ), 
										df_final_CR_ref$CR_SLSL_data,df_final_CR_ref$CR_cumulated_SLSJ)

df_final_CR_ref <- df_final_CR_ref[complete.cases(df_final_CR_ref$Carrier_rate_reported_SLSJ),]
## Create the equation
lm_eqn <- function(df){
	m <- lm(Carrier_rate_reported_SLSJ~ CR_SLSL_final, df);
	eq <- substitute(~~italic(R)^2~"="~r2, 
					 list(r2 = format(summary(m)$r.squared, digits = 3)))
	as.character(as.expression(eq));
}

## Remove dot not good
df_final_CR_ref <- df_final_CR_ref[df_final_CR_ref$ID != "21439",]
## Do graph
CR_graph <- ggplot(data = df_final_CR_ref,aes(x = CR_SLSL_final, y = Carrier_rate_reported_SLSJ)) + geom_point()+
	scale_shape_manual(values=c(1,5))+
	theme(axis.text.x = element_text(face="bold", size = 10,angle=45, hjust = 1),axis.title=element_text(size=10,face="bold"),
		  axis.text.y = element_text(face="bold", size = 10),legend.text = element_text(size=16)) + 
	theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
	scale_x_continuous("CR found in our analysis in SLSJ",breaks=c(0,25,50,75,100,125,150,175,200,225),limits=c(0,225),labels = c("0","1/25","1/50","1/75","1/100","1/125","1/150","1/175","1/200","1/225")) +
	scale_y_continuous("CR reported in the literature in SLSJ",breaks=c(0,25,50,75,100,125,150,175,200,225),limits=c(0,225),
					   labels = c("0","1/25","1/50","1/75","1/100","1/125","1/150","1/175","1/200","1/225"))+
	geom_text(x = 50, y = 175, label = lm_eqn(df_final_CR_ref), parse = TRUE,size=8)

## Save in PDF format
outfile <- paste0("Figure_supp_1.pdf")
cairo_pdf(outfile,  width = 20, height = 10)
print(CR_graph)
dev.off()

## Save in png format
ggsave("Figure_supp_1.pdf", width = 107, height = 107,bg="white",dpi=300, units="mm")



################################################################# Supp Fig.5 ############################################################
######### Graph Correlation Imput / WGS 
## Select data for Graph 
df_final_hist <- df_final[,c(1:23)]
## Select variants with WGS and Imput info
df_final_hist_CR_both_complete <- df_final_hist[complete.cases(df_final_hist$Ind_Imput_WQ),]
df_final_hist_CR_both_complete <- df_final_hist_CR_both_complete[complete.cases(df_final_hist_CR_both_complete$Ind_WGS_WQ),]
## Calculate the CR
df_final_hist_CR_both_complete$Prop_WGS_WQ <- df_final_hist_CR_both_complete$Ind_WGS_WQ/1852
df_final_hist_CR_both_complete$Prop_Imput_WQ <- df_final_hist_CR_both_complete$Ind_Imput_WQ/25061
## Do the equation
lm_eqn <- function(df){
	m <- lm(Prop_Imput_WQ~ Prop_WGS_WQ, df);
	eq <- substitute(~~italic(R)^2~"="~r2, 
					 list(r2 = format(summary(m)$r.squared, digits = 3)))
	as.character(as.expression(eq));
}
## Do the figure
figure_method_compare_WGS_Imput <- ggplot(data = df_final_hist_CR_both_complete,aes(x = Prop_WGS_WQ, y = Prop_Imput_WQ)) + geom_point(cex=0.5)+
	scale_shape_manual(values=c(1,5))+
	theme(axis.text.x = element_text( size = 12),axis.title=element_text(size=12),
		  axis.text.y = element_text( size = 12),legend.text = element_text(size=12)) + 
	xlim(c(0,0.02)) + ylim(c(0,0.02)) +
	theme(legend.title = element_blank())+ guides(color=guide_legend(override.aes = list(size=4)))+
	xlab("Variant frequency in WGS") + ylab("Variant frequency in imputed data")+
	geom_text(x = 0.002, y = 0.015, label = lm_eqn(df_final_hist_CR_both_complete), parse = TRUE,size=6)

## Save in tiff format
ggsave("Figure_supp_5.tiff", width = 132, height = 132,units = 'mm',dpi=600)
