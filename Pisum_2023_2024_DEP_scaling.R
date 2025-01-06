#setting up working directory
setwd("C:/AkilaY/Work/Akila/SUB968_Pisum/Pisum_2023_2024")


#### Installing and loading libraries
chooseCRANmirror(ind=1) #selecting cran mirror option
packages = c("BiocManager","tidyverse", 'DEP', 'limma','data.table','org.At.tair.db','parallel',
				'doParallel','SummarizedExperiment','ggrepel','pheatmap','ggalt','ggvenn','openxlsx','UpSetR')
				
for(lib in packages){
if(!lib %in% installed.packages()){
if(lib %in% available.packages()[,1]){
install.packages(lib)} else{ BiocManager::install(lib)}}}

suppressMessages(sapply(packages, require, character = TRUE))		#loading libraries with suppressing any messages

############################ DATA PREPARATION ##################################

#loading data
data_v1 = fread("Pisum_2023_2024_report.tsv")
dim(data_v1)

#Instruction 1. Only work with protein groups that have PG.Q.Value <0.01
data_v1_n = data_v1[data_v1$PG.Q.Value<0.01,]
dim(data_v1_n)

#Instruction 3. Filter your data to only have the following columns: Run, Protein.Group, Protein.Names, and PG.MaxLFQ.
#Here we also include Protein.IDs column as it is required for make_se_parse() function.

data_v2 = data_v1_n %>% dplyr::select(Run, Protein.Group, Protein.Names, Protein.Ids, PG.MaxLFQ)

#Instruction 4. Filter that smaller dataset to only contain unique Run & Protein.Group (where every row now represents a specific protein group for each run)

data_v2_F = data_v2[!duplicated(data_v2[,c('Run','Protein.Group')]),]
dim(data_v2_F)

#Reshape data.
# 1. Rows should be unique. for this ""Protein.Names" or "Protein.Ids" can be used.
# 2. reshape table using pivot_wider(names_from = 'Run', values_from = "PG.MaxLFQ")

data_v2_T <- data_v2_F %>% pivot_wider(names_from = 'Run', values_from = 'PG.MaxLFQ')
dim(data_v2_T)

colnames(data_v2_T)

rename_data = read.xlsx("Sample names_Pisum_2023_2024.xlsx")


#Rename columns
data <- data_v2_T %>% rename("10_Axillary_Buds"="Axillary_Buds_0h_R1",
				"20_Axillary_Buds"="Axillary_Buds_0h_R2",
				"40_Axillary_Buds"="Axillary_Buds_0h_R3",
				"20221028_Akila_sample8_Apical_Buds_0h_R3"="Apical_Buds_0h_R3",
				"20221028_Akila_sample9_Apical_Buds_2h_R4"="Apical_Buds_2h_R4",
				"20221028_Akila_sample1_2_Axillary_Buds_R2"="Axillary_Buds_2h_R2",
				"20221028_Akila_sample1_Axillary_Buds_R1"="Axillary_Buds_2h_R1",
				"20221028_Akila_sample2_Axillary_Buds_R3"="Axillary_Buds_2h_R3",
				"20221028_Akila_sample3_Apical_Buds_2h_R1"="Apical_Buds_2h_R1",
				"20221028_Akila_sample4_Apical_Buds_2h_R2"="Apical_Buds_2h_R2",
				"20221028_Akila_sample7_Apical_Buds_0h_R2"="Apical_Buds_0h_R2")
				
colnames(data)




#Make unique names using the annotation in the "Gene.Names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have a gene name.
data_unique <- make_unique(data, "Protein.Names", "Protein.Ids", delim = ";")

data_unique[,12:14] = data_unique[,12:14]/12.5		##Scaling the 5 µg Samples (5/0.4 = 12.5)

#preparing experimnetal design
experimental_design = data.frame(label = colnames(data_unique)[4:14], condition = sub("_[^_]*$", "", colnames(data_unique)[4:14]),replicate = gsub("[^_]+_", "", colnames(data_unique)[4:14]))

assay_columns <- 4:14	#getting the assay columns


# Generate a SummarizedExperiment object using an experimental design
data_se <- make_se(data_unique, assay_columns, experimental_design)

data_se



#Bar plot for protein counts in each Sample
protein_counts = data.frame(samples = colnames(data_se))
for(i in 1:11){
protein_counts$counts[i] = table(is.na(assay(data_se)[,i]))[2]
}


p <- ggplot(protein_counts,aes(x=counts,y=samples))+
		geom_bar(stat="identity",fill="darkorchid4") +
		geom_text(aes(label=counts),color="white" ,hjust=1.5) +
		ggtitle("Protein counts in each sample")+
		theme_bw()+ 
		theme(
			plot.title = element_text(size = 18, face = 'bold', hjust=0.5, vjust = 1),
			axis.text = element_text(face = 'bold', color ='black'),
			axis.title = element_text(size = 12, face = 'bold'),
			axis.title.y = element_blank(),
			axis.line = element_line(colour = 'black')) 



jpeg("DEP/Pisum_2023_2024_V2/Protein_counts_BarPlot_V2.jpg", height = 15*350, width = 25*350, res = 800)
p
dev.off()

#heatmap of common proteins among samples using raw data
common_data = na.omit(assay(data_se))
col = colorRampPalette(c("purple3", "white", "red2"))(15)

#Heatmap of all common proteins
out=pheatmap(t(common_data),col = col,scale="none",main="\nHeatmap of Common proteins among samples using raw data",
		border_color='gray60',fontsize=25, fontsize_col = 15,cluster_rows = FALSE,show_colnames = F) 


jpeg("DEP/Pisum_2023_2024_V2/heatmap_Pisum_common_proteins_V2.jpg", width = 100*350,  height = 42*350,res=1000)
out
dev.off()

							
#Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
				
# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)				

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
#The data is background corrected and normalized by variance stabilizing transformation (vsn).
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
# To explore the pattern of missing values in the data, a heatmap is plotted indicating whether values are missing (0) or not (1). 
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

set.seed(1)

# here we will compare left censored methods and mixed imputation method
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp_MinProb <- impute(data_norm, fun = "MinProb", q = 0.01)		

##Mixed imputation on proteins (rows)
# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(data_norm) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(data_norm) %in% proteins_MNAR

# Perform a mixed imputation
data_imp_mixed <- impute(
  data_norm, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "MLE", # imputation function for MAR
  mnar = "MinProb") # imputation function for MNAR
  
  

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp_MinProb,data_imp_mixed)


### Checking the data before going for differential analysis by plotting boxplot
#boxplot
colData(data_se)$label = factor(colData(data_se)$label,levels = c("Axillary_Buds_0h_R1", "Axillary_Buds_0h_R2", "Axillary_Buds_0h_R3",
									"Apical_Buds_0h_R3","Apical_Buds_2h_R4",
									"Axillary_Buds_2h_R2","Axillary_Buds_2h_R1","Axillary_Buds_2h_R3",
									"Apical_Buds_2h_R1","Apical_Buds_2h_R2","Apical_Buds_0h_R2"))
									
									
									
jpeg("DEP/Pisum_2023_2024_V2/BOXPLOT_Pisum_Comparison_V2.jpg", height = 15*350, width = 40*350,res=500)
layout(matrix(c(1,2,3,4,4,4), 2, 3, byrow = TRUE), heights=c(1,0.2))
par(mar=c(12,7,6,3),font=2,font.main=2,font.axis=2)
col=rep(c("lightskyblue","olivedrab3","limegreen","steelblue","limegreen","olivedrab3"),c(3,1,1,3,2,1))
boxplot(assay(data_se)[,order(colData(data_se)$label)],target="core",las=2,main=c("Non-Normalized Data"), outline=T, col=col,cex.main=3,cex.axis=2 )
mtext("Intensities", side=2, line=3.5,cex=2,font=2)


par(mar=c(12,1.5,6,3),font=2,font.main=2,font.axis=2)
boxplot(assay(data_norm)[,order(colData(data_se)$label)],target="core",las=2,main=c("VSN Normalized Data"), outline=T, col=col,cex.main=3,cex.axis=2)
#mtext("Expression Values", side=2, line=5.2,cex=1.5,font=2)

par(mar=c(12,1.5,6,3),font=2,font.main=2,font.axis=2)
boxplot(assay(data_imp_mixed)[,order(colData(data_se)$label)],target="core",las=2,main=c("Mixed Imputed Data"), outline=T, col=col,cex.main=3,cex.axis=2)

par(mar = c(1, 0, 0, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('bottom',legend = c("Axillary_Buds_0h","Axillary_Buds_2h","Apical_Buds_0h","Apical_Buds_2h"), col = c("lightskyblue","steelblue","olivedrab3","limegreen"), xpd = TRUE, horiz = TRUE, lwd=25,seg.len=1,cex = 2, bty = 'n')
   
dev.off()



# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# manual comparisons of samples



contrasts =  c("Axillary_Buds_2h_vs_Axillary_Buds_0h", "Apical_Buds_2h_vs_Apical_Buds_0h", "Axillary_Buds_0h_vs_Apical_Buds_0h")

###### --------------Differential enrichment analysis-------------- ######

# Differential enrichment analysis  based on linear models and empherical Bayes statistics using limma workflow


#design matrix
conditions = factor(colData(data_imp_mixed)$condition)		
design <- model.matrix(~0+conditions, data = colData(data_imp_mixed))
colnames(design) <- gsub("conditions", "", colnames(design))


fit <- lmFit(assay(data_imp_mixed), design)	#fitting model

#comparisons we want for this study

		 
cntrst <- gsub("_vs_", " - ", contrasts)		 

made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)		#making contrasts

contrast_fit <- contrasts.fit(fit, made_contrasts)		#fitting contrasts into the model

eB_fit <- eBayes(contrast_fit)		#empirical bayes statistics for differential expression

#creating a function to retrieve the statistics results for each comparisons
retrieve_fun <- function(comp, fit = eB_fit) {
        res <- topTable(fit, sort.by = "t", coef = comp, number = Inf, 
            confint = TRUE)
        res <- res[!is.na(res$t), ]
        res$comparison <- rep(comp, dim(res)[1])
        res <- rownames_to_column(res)
        return(res)
    }
limma_res <- map_df(cntrst, retrieve_fun)	#applyinf retrieve function to each comparison

#preparing the result table as in the DEP package
table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, CI.R, 
        P.Value, adj.P.Val, comparison) %>% mutate(comparison = gsub(" - ", 
        "_vs_", comparison)) %>% gather(variable, value, -c(rowname, 
        comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
        P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
        variable) %>% spread(temp, value)
		
#creating the summarizedExperiment object with differential results		
data_diff <- data_imp_mixed		
rowData(data_diff) <- merge(rowData(data_imp_mixed, use.names = FALSE), table, 
        by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
		
# Denote significant proteins based on user defined cutoffs
#here I have written a different code instead of using DEP's 'add_rejections' function
dep=data_diff		#assigning se object from differential enrichment analysis to a new variable
row_data <- rowData(dep) %>% as.data.frame()	#extracting rowData of the summarizedExperiment object	
    
cols_p <- grep("_p.adj", colnames(row_data))		#getting adjusted p-values columns
cols_diff <- grep("_diff", colnames(row_data))		#getting log fold change columns
    
p_reject <- row_data[, cols_p] <= 0.05				#adjusted pvalue threshold

diff_reject <- abs(row_data[, cols_diff]) >= log2(1.5)		#adding log fold change threshold
        
sign_df <- p_reject & diff_reject						
sign_df <- cbind(sign_df, significant = apply(sign_df, 1, function(x) any(x)))		#adding significance
colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
rowData(dep) <- merge(rowData(dep), sign_df, by = "name",sort=FALSE)		#adding significance to the rowData of the summarizedExperiment object of differentisl result
    
  
rowData(dep)$significant%>%table()	#getting number of significant proteins




###### ---------------Visualization of the results-------------- ######

######################
##### PCA PLOT #######
######################

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = nrow(dep), point_size = 4)

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = nrow(dep),indicate = "condition", point_size = 4)


#pca plot using ggplot2 package
pca <- prcomp(t(assay(data_imp_mixed)), scale = FALSE)
pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>% 
        left_join(., data.frame(colData(dep)), by = c(rowname = "ID"))
percent <- round(100 * pca$sdev^2/sum(pca$sdev^2), 1)

p = ggplot(pca_df, aes(x=PC1, y=PC2, color = condition)) + 
		geom_point(size=5) + 
		labs(title = "\nPCA plot: Pisum_2023_2024_DEP", x = paste0("PC1: ", percent[1], "% variation"), y = paste0("\nPC2: ", percent[2], "% variation")) + 
		geom_text_repel(aes(label=rowname),color='black',show.legend = FALSE,size=3, max.overlaps = Inf,force=35)+
		scale_color_manual(values=c("olivedrab3","limegreen","lightskyblue","steelblue"))+
 		theme_minimal() +
    	theme(	axis.text = element_text(size=12,face = 'bold', color ='black'),
				axis.title=element_text(size=15, face = 'bold'),
				plot.title = element_text(hjust = 0.5, size = 22, face = 'bold', vjust = 1),
				legend.text = element_text(size = 13,face = "bold"),
				title = element_text(size = 15, face = "bold"),				
				panel.border = element_rect(color='black',fill= "transparent"))

p		

jpeg("DEP/Pisum_2023_2024_V2/PCA_Pisum_DEP_V2.jpg", height = 30*350, width = 35*350,res=1200)
p
dev.off()



###### Plot a P value histogram

plot_p_hist(dep)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 3, show_row_names = TRUE)


#heatmap of aggregated intensities for each conditions
#getting the rowmeans for each condition
sig_pro_rowdata = dep[rowData(dep)$significant==TRUE,]
sig_pro_t = rowData(sig_pro_rowdata)[,c(1,grep("_p.adj", colnames(rowData(sig_pro_rowdata))))]
sig_pro_t <- melt(data.frame(sig_pro_t), id.vars = "name", variable.name = "comparison")
sig_pro_top50 = sig_pro_t%>%group_by(comparison)%>% arrange(value)%>% dplyr::slice(1:50)

sig_pro = assay(dep)[rowData(dep)$significant==TRUE,]
sig_pro_m = melt(sig_pro)	#changing to long format
sig_pro_m$condition = sub("_[^_]*$","",sig_pro_m$Var2)		#adding condition column

sig_pro_agg = aggregate(value ~ Var1+condition, data = sig_pro_m, FUN = mean)		#aggregating intensities for each condition

sig_pro_w = sig_pro_agg %>% pivot_wider(names_from = condition, values_from = value)%>% as.data.frame()
rownames(sig_pro_w) = sig_pro_w[,1]
sig_pro_w = sig_pro_w[,-1]

top50_pro = sig_pro_w[rownames(sig_pro_w)%in%sig_pro_top50$name,]
#color palette for heatmap  
col = colorRampPalette(c("royalblue3", "white", "brown2"))(15)

#Heatmap of all significant proteins in each condition
out=pheatmap(t(top50_pro),col = col,scale="none",main="\nIntensity Heatmap of Top 50 significant proteins in Each comparison",cellheight = 150,
		cellwidth = 15,border_color='gray60',fontsize=25, fontsize_col = 15,cluster_rows = FALSE) 


jpeg("DEP/Pisum_2023_2024_V2/heatmap_Pisum_DEP_top50_V2.jpg", width = 90*350,  height = 42*350,res=1000)
out
dev.off()


############################################################ heatmap of contrasts ##################################################################

filtered = dep[rowData(dep)$significant, ]
df <- rowData(filtered, use.names = FALSE) %>% data.frame() %>% 
            column_to_rownames(var = "name") %>% dplyr::select(ends_with("_diff"))
        colnames(df) <- gsub("_diff", "", colnames(df)) %>% gsub("_vs_", 
            " vs ", .)
df <- as.matrix(df)

top50_pro_logfc = df[rownames(df)%in%sig_pro_top50$name,]

col = colorRampPalette(c("darkgreen", "white", "gold"))(15)

#Heatmap of all significant proteins in each condition
out=pheatmap(t(top50_pro_logfc),col = col,scale="none",main="\nHeatmap of top 50 significant proteins in tested Contrasts (log Fold Change)",
		cellheight = 180,cellwidth = 13,border_color='gray60',fontsize=25, fontsize_col = 15) 


jpeg("DEP/Pisum_2023_2024_V2/heatmapOfcontrast_Pisum_DEP_top50_V2.jpg", width = 90*350, height = 42*350,res=1000)
out
dev.off()







### volcano plot of all comparisons using ggplot2 package

#preparing data for volcano plot
volcdat_p = rowData(dep)[c(1,grep("_p.adj", colnames(rowData(dep))))] #extracting adjusted p-values for all comparisons
volcdat_p = melt(as.data.table(volcdat_p)) 		#melting the data
volcdat_p$variable = gsub("_p.adj",'',volcdat_p$variable)		
colnames(volcdat_p)[3] = "p.adj"
volcdat_d = rowData(dep)[c(1,grep("_diff", colnames(rowData(dep))))]	#extracting lof fold change for all comparisons
volcdat_d = melt(as.data.table(volcdat_d))
volcdat_d$variable = gsub("_diff",'',volcdat_d$variable)
colnames(volcdat_d)[3] = "lfc"
volcdat_sig = rowData(dep)[c(1,grep("_significant", colnames(rowData(dep))))]			#extracting significance column for all comparisons
volcdat_sig = melt(as.data.table(volcdat_sig))
volcdat_sig$variable = gsub("_significant",'',volcdat_sig$variable)
colnames(volcdat_sig)[3] = "significant"

#merging melted p.adj, lfc and significance datasets into one
volcdat =  merge(volcdat_p,volcdat_d,sort=F)
volcdat = merge(volcdat,volcdat_sig,sort=F)

#creating a column to define up and down regulation of the compound in the comparisons
volcdat$diff <- "NO"
volcdat$diff[volcdat$lfc > log2(1.5) & volcdat$p.adj < 0.05] <- "UP"
volcdat$diff[volcdat$lfc < -log2(1.5) & volcdat$p.adj < 0.05] <- "DOWN"
volcdat$delabel <- NA
volcdat$delabel[volcdat$diff != "NO"] <- volcdat$name[volcdat$diff != "NO"] #adding name of significant compounds to show on the plot
sig_volc = volcdat[volcdat$p.adj<0.05,]
sig_volc = sig_volc[abs(sig_volc$lfc)>log2(1.5),]


labels = rbind(sig_volc[sig_volc$diff=='UP',],sig_volc[sig_volc$diff=='DOWN',])			#label data
#preparing top 10 up and down-regulated gene labels data for volcano plots

labelData<- labels %>%                                     
  arrange(p.adj) %>%
  group_by(variable,diff) %>%
  dplyr::slice(1:10)
  
#volcano plot
p = ggplot(volcdat, aes(x=lfc, y=-log10(p.adj),col=diff)) + 
		geom_hline(yintercept=-log10(0.05), col="gray", linetype = 'dashed',linewidth=1)+
		geom_vline(xintercept=c(-log2(1.5), log2(1.5)), col="gray", linetype = 'dashed',linewidth=1) +
		geom_point(size = 2) +geom_point(data=sig_volc,aes(x=lfc,y=-log10(p.adj)),size = 3)+guides(color = guide_legend(override.aes=list(size=5)))+
		facet_wrap(.~variable , scales = "free",nrow=2)+
		scale_color_manual(name = "Expression", values = c("red","gray80","limegreen"),
                     labels = c("Down-regulated","Not significant", "Up-regulated"))+
		geom_text_repel(data = labelData,aes(label=delabel),col='black',max.overlaps=Inf,fontface="bold",show.legend=FALSE,size=3,force=30) +
		labs(title="\nVolcano Plot: Biological and Statistical Significance\n", x="\nlog2 Fold change", y="-log10 P-value\n")+
		
		theme_minimal()+ 
  theme(
    panel.border = element_rect(fill= "transparent",colour='black'),
    plot.title = element_text(size = 15, face = 'bold', hjust=0.5, vjust = 1),
	strip.text = element_text(size = 12,color="black",face='bold',hjust=0.5),
    axis.text = element_text(size=8,face = 'bold', color = 'black'),
    axis.title = element_text(size = 11, face = 'bold',color = 'black'),
    axis.line = element_line(colour = 'black'),
    legend.text = element_text(size=9,face = "bold"), # Text size
	title = element_text(size = 11, face = "bold")) 
	
jpeg("DEP/Pisum_2023_2024_V2/Volcano_Pisum_DEP_V2.jpg", height = 45*350, width = 55*350, res = 1200)
p
dev.off()	



# Plot a barplot for A0A0F6NGI2_PEA and A0A9D4VLH4_PEA
plot_single(dep, proteins = c("A0A0F6NGI2_PEA", "A0A9D4VLH4_PEA"))	#selected randomly from significant proteins

# Plot a barplot for the protein A0A9D4VLH4_PEA with the data centered
plot_single(dep, proteins = "A0A9D4VLH4_PEA", type = "centered") 


###Frequency plot of significant proteins and overlap of conditions
# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)


# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

# Column names of the results table
colnames(data_results)

#significant proteins
sig_protein = data_results %>% filter(significant)

#wwriting result table
write.csv(data_results, 'DEP/Pisum_2023_2024_V2/Pisum_DEP_result_table_V2.csv',row.names=F)


####Separate results table for the comparisons you asked for

#result table for d14_vs_Col0
Apical_Buds_2h_vs_Apical_Buds_0h_result = data_results[,c(1,2,grep("Apical_Buds_2h_vs_Apical_Buds_0h",colnames(data_results)))]%>% filter(Apical_Buds_2h_vs_Apical_Buds_0h_significant)

write.csv(Apical_Buds_2h_vs_Apical_Buds_0h_result, "DEP/Pisum_2023_2024_V2/Apical_Buds_2h_vs_Apical_Buds_0h_DEP_result_table.csv",row.names=F)

#result table for max2_vs_Col0
Axillary_Buds_0h_vs_Apical_Buds_0h_result = data_results[,c(1,2,grep("Axillary_Buds_0h_vs_Apical_Buds_0h",colnames(data_results)))]%>% filter(Axillary_Buds_0h_vs_Apical_Buds_0h_significant)

write.csv(Axillary_Buds_0h_vs_Apical_Buds_0h_result, "DEP/Pisum_2023_2024_V2/Axillary_Buds_0h_vs_Apical_Buds_0h_DEP_result_table.csv",row.names=F)

#result table for smxl678_vs_Col0
Axillary_Buds_2h_vs_Axillary_Buds_0h_result = data_results[,c(1,2,grep("Axillary_Buds_2h_vs_Axillary_Buds_0h",colnames(data_results)))]%>% filter(Axillary_Buds_2h_vs_Axillary_Buds_0h_significant)

write.csv(Axillary_Buds_2h_vs_Axillary_Buds_0h_result, "DEP/Pisum_2023_2024_V2/Axillary_Buds_2h_vs_Axillary_Buds_0h_DEP_result_table.csv",row.names=F)

#venn diagram
library("ggvenn")

x = list(A = data_results$name[data_results$Apical_Buds_2h_vs_Apical_Buds_0h_significant],
		B = data_results$name[data_results$Axillary_Buds_0h_vs_Apical_Buds_0h_significant],
		C = data_results$name[data_results$Axillary_Buds_2h_vs_Axillary_Buds_0h_significant])


names(x) <- c("Apical_Buds_2h_vs_Apical_Buds_0h","Axillary_Buds_0h_vs_Apical_Buds_0h","Axillary_Buds_2h_vs_Axillary_Buds_0h")
jpeg("DEP/Pisum_2023_2024_V2/Pisum_venn_V2.jpg",height=25*350,width=35*350,res=1200)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


######################################### VENN DIAGRAM ##############################################################
#preparing data for venn diagram
labeldata=volcdat[volcdat$diff!='NO',]		#getting significant proteins
vennData=labeldata
#vennData$variable = gsub("contrast: ",'',vennData$variable)
signif_genelist <- list()
for(j in 1:length(unique(vennData$variable))){
for(k in 1:length(unique(vennData$diff))){
x = list(vennData$name[vennData$variable==unique(vennData$variable)[j]&vennData$diff==unique(vennData$diff)[k]])
name <- paste(unique(vennData$variable)[j],'@',unique(vennData$diff)[k],sep='')
signif_genelist[name] <- x
}}



#install.packages('ggVennDiagram')
library(ggVennDiagram)

venn =ggVennDiagram(signif_genelist, label_alpha = 0,set_size = 3.5) +
			labs(title="Venn Diagram",subtitle="Significant proteins in all tested contrasts")+
			scale_fill_gradient(low="cornsilk",high = "aquamarine")+
			#scale_fill_distiller(palette = "RdBu")+ 
			#scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
			guides(fill = guide_legend(title = "Count")) +
			scale_x_continuous(expand = expansion(mult = .2))+
			theme(legend.position = "bottom",
				plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
				plot.subtitle = element_text(size = 12, face = 'bold', hjust = 0.5),		
				legend.text = element_text(face = "bold"), # Text size
				title = element_text(size = 11, face = "bold")) 
  
  
jpeg("DEP/Pisum_2023_2024_V2/Pisum_DEP_vennDiagram_V2.jpg", height = 20*350, width = 30*350, res = 800)
venn
dev.off()

listInput = signif_genelist

#UpsetR plot
jpeg("DEP/Pisum_2023_2024_V2/upsetPlot_V2.jpg",height = 10*350, width = 20*350, res = 300)
upset(fromList(listInput), nsets=6,mb.ratio = c(0.55, 0.45),order.by = "freq",
	keep.order = TRUE,main.bar.color = "darkred",mainbar.y.label="Intersection Size",matrix.color = "black", 
	sets.bar.color = c("darkorchid4","coral1","darkolivegreen","darkgoldenrod","royalblue","deeppink"),
	point.size=5,line.size = 1, shade.color = "palegreen",shade.alpha = 0.5,text.scale=2)
	
	
dev.off()	


  
# getting gene names of common genes in all comparisons 
common = data.frame(Common_Genes=Reduce(intersect,x))
write.csv(common,"DEP/Pisum_2023_2024_V2/Pisum_COMMON_DEP_V2.csv")


#getting uncommon genes in each comparison
uncommon = list(Apical_2h_vs_Apical_0h = setdiff(setdiff(x$Apical_Buds_2h_vs_Apical_Buds_0h,x$Axillary_Buds_0h_vs_Apical_Buds_0h),x$Axillary_Buds_2h_vs_Axillary_Buds_0h),
Axillary_0h_vs_Apical_0h = setdiff(setdiff(x$Axillary_Buds_0h_vs_Apical_Buds_0h,x$Apical_Buds_2h_vs_Apical_Buds_0h),x$Axillary_Buds_2h_vs_Axillary_Buds_0h),
Axillary_2h_vs_Axillary_0h = setdiff(setdiff(x$Axillary_Buds_2h_vs_Axillary_Buds_0h,x$Apical_Buds_2h_vs_Apical_Buds_0h),x$Axillary_Buds_0h_vs_Apical_Buds_0h))

write.xlsx(uncommon,"DEP/Pisum_2023_2024_V2/Pisum_UNCOMMON_DEP_V2.xlsx")

