################################
#
# September 2024 by ME
#
# Do Figures in article
#
# Modified by CM
################################

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



options(scipen=999)
## File to open with all variants enriched (1304)
df_final <- read.table("comparison_WGS_Imput_Gnomad_all_variants_enriched__all_info_SAG_RQ_WQ.txt",header=T)
## Open file of variants known
variants_known <- read.table("variants_reported_in_SLSJ.txt")

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
tiff("Figure_3.tiff",width =132, height = 132,units = 'mm', res=600 ,compression = "lzw")
founder_WQ_known
dev.off()



################################################################# Fig.4 ############################################################
##### Figure Simon with CR > 1/1000 for SLSJ
## Select row I want

CR_for_graph <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11,26,27)]
## Remove variants to freq
variants_to_freq <- c("chr6_26090951_C_G","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph <- CR_for_graph[!CR_for_graph$SNP %in% variants_to_freq,]
## Select founder in SLSJ and UQc
CR_for_graph$Founder_SAG[!is.na(CR_for_graph$CR_Imput_SAG)] <- CR_for_graph$founder_Imput_SAG[!is.na(CR_for_graph$CR_Imput_SAG)]
CR_for_graph$Founder_RQ[!is.na(CR_for_graph$CR_Imput_RQ)] <- CR_for_graph$founder_Imput_RQ[!is.na(CR_for_graph$CR_Imput_RQ)]
CR_for_graph$Founder_WQ[!is.na(CR_for_graph$CR_Imput_WQ)] <- CR_for_graph$founder_Imput_WQ[!is.na(CR_for_graph$CR_Imput_WQ)]

CR_for_graph <- CR_for_graph[CR_for_graph$Founder_SAG == "Founder" | CR_for_graph$Founder_RQ == "Founder"| CR_for_graph$Founder_WQ == "Founder",]
CR_for_graph <- CR_for_graph[!is.na(CR_for_graph$SNP),]


## Only take imput
CR_for_graph$CR_SAG <- CR_for_graph$CR_Imput_SAG[!is.na(CR_for_graph$CR_Imput_SAG)]
CR_for_graph$CR_RQ <- CR_for_graph$CR_Imput_RQ[!is.na(CR_for_graph$CR_Imput_RQ)]


CR_for_graph[is.na(CR_for_graph)] <- -1

# Create bins with intervals of 20 for UQc
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
sum(combined_frequency$Freq[combined_frequency$dataset=="SAG"])
########################## Enriched 
## Select row I want
CR_for_graph_enriched <- df_final[,c(1,2,6,9,14,17,22,25,18,19,10,11)]
## Remove variants to freq
variants_to_freq <- c("chr6_26090951_C_G","chr14_94380925_T_A","chr4_102901325_AAC_A")
CR_for_graph_enriched <- CR_for_graph_enriched[!CR_for_graph_enriched$SNP %in% variants_to_freq,]

## Get CR (either imput or if missing)
CR_for_graph_enriched$CR_SAG[!is.na(CR_for_graph_enriched$CR_Imput_SAG)] <- CR_for_graph_enriched$CR_Imput_SAG[!is.na(CR_for_graph_enriched$CR_Imput_SAG)]
CR_for_graph_enriched$CR_RQ[!is.na(CR_for_graph_enriched$CR_Imput_RQ)] <- CR_for_graph_enriched$CR_Imput_RQ[!is.na(CR_for_graph_enriched$CR_Imput_RQ)]
## Only take imput variants
CR_for_graph_enriched <- CR_for_graph_enriched[!is.na(CR_for_graph_enriched$CR_SAG) | !is.na(CR_for_graph_enriched$CR_RQ),]


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



################################################################# Supp Fig.4 ############################################################
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
ggsave("Figure_supp_4.tiff", width = 132, height = 132,units = 'mm',dpi=600)


## Save in tiff format
ggsave("Figure_supp_5.tiff", width = 132, height = 132,units = 'mm',dpi=600)
