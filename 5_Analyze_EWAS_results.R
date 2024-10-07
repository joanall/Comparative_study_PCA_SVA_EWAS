#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 5: Analyze results                                #
#                                                                              #
#  <Display description>                                                       #----
#==============================================================================#

# Objective: Assess whether performing PCA and SVA on EWAS data and adding PCs 
#or SVs on the models changes the results by truly removing unknonw variation effects
# and therefore uncovering stronger/clearer biological signals. 

# STEP 1: Get Residuals to remove the effect of known covariates from DNAm 
# Residuals represent the part of DNAm not explained by the covariates. 
#   Methylation ∼ Exposure + Covariates

# STEP 2: SVA or PCA on Residuals to identify hidden patterns of variation- 

# STEP 3: Choose optimal SVs or PCs ("optional") by applying elbow method

# STEP 4: Running EWAS with and without SVs or PCs
#   Methylation ∼ Exposure + Covariates +  SVs/PCs

# STEP 5: Compare hits and statistics within approaches and with reference study.



# [INPUT]:
# - RDS with PACeanalysis object after running EWAS 
# - Dataframe with unreliable probes for EPIC and/or 450K array 

# [OUPUT]:
# - Clean EWAS results (remooval of unreliable plots)
# - Interactive number of CpGs surpassing bonferroni and suggestive threshold
# - Manhattan plot 
# - QQ plots with lambda-CI
# - Venn plot comparing obatined results with reference hits

# [CONSIDERATIONS]
#===============================================================================

#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# BISC
# - basic : 'sva_pca/results/archive/bisc_16092024_Output/maternal_smoking_smoke_basic_pace/bisc_16092024_maternal_smoking_smoke_basic_pace_allanalyses.RData'
# - pca: 'sva_pca/results/bisc_26092024_Output/maternal_smoking_6opt_pcs/bisc_26092024_maternal_smoking_opt_pcs_allanalyses.RData'
# - 22 pcs 'sva_pca/results/bisc_21092024_Output/maternal_smoking_22opt_pcs/bisc_21092024_maternal_smoking_22opt_pcs_allanalyses.RData'
# - 6 sva: 'sva_pca/results/bisc_26092024_Output/maternal_smoking_6opt_svs/bisc_26092024_maternal_smoking_opt_svs_allanalyses.RData'
# - meffil : 'sva_pca/results/meffil_ewas_21092024.rds'

# INMA 
# - basic: 'sva_pca/results/INMA_21092024_Output/Smoke_basic_ewas/INMA_21092024_Smoke_basic_ewas_allanalyses.RData'
# - 5 pcs: sva_pca/results/INMA_21092024_Output/Smoke_5opt_pcs/INMA_21092024_Smoke_5opt_pcs_allanalyses.RData
# - 16 pcs: 'sva_pca/results/INMA_21092024_Output/Smoke_16opt_pcs/INMA_21092024_Smoke_16opt_pcs_allanalyses.RData'
# - 7 sv: 'sva_pca/results/INMA_21092024_Output/Smoke_7opt_svs/INMA_21092024_Smoke_7opt_svs_allanalyses.RData'
# - meffil 

ewas_file <- 'sva_pca/results/INMA_21092024_Output/Smoke_16opt_pcs/INMA_21092024_Smoke_16opt_pcs_allanalyses.RData'
filter_EPIC_cpgs_path <- "db/final_data/filter_EPIC.RData"

array <- 'EPIC'
cohort <- 'INMA'
analysis_name <- 'opt16pcs'
current_date <- '21092024'
file_from <- 'PACE' # "PACE" or "meffil"


#==============================================================================#
#                                                                              #
#                        INIT WORKING ENVIORNEMT                               #
#                                                                              #
#==============================================================================#

# Load libraries

library(ggplot2)
library(data.table)
library(ewascatalog)
library(VennDiagram)
library(grid)
library(dplyr)

#------------------------------------------------------------------------------#
#                         Define functions                                     #----
#------------------------------------------------------------------------------#

 
do_manhattan <- function(results, array, suggestive_cutoff, filename_save){
  
  # Get annotation
  if (array == "450K") {
    library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
    annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (array == "EPIC") {
    library('IlluminaHumanMethylationEPICanno.ilm10b2.hg19')
    annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  }
  
  annotation <- as.data.frame(annotation)
  
  # Add annotation to results by merging on row names = CpGs
  results <- merge(results, annotation, by.x = 'row.names', by.y = 'row.names')
  
  # Get numeric for chromosome (e.g. from "chr6" to 6)
  results$chromosome <- as.numeric(gsub('chr', "", results$chr))
  
  # Add a column for -log10 p-values
  results$logp <- -log10(results$pval)
  
  # Prepare data to make plot
  don <- results %>%
    # Compute chromosome size
    group_by(chromosome) %>%
    summarise(chr_len = max(as.numeric(pos))) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(results, by = c("chromosome" = "chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(chromosome, pos) %>%
    mutate(BPcum = as.numeric(pos) + tot)
  
  axisdf <- don %>%
    group_by(chromosome) %>%
    summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)
  
  # Define thresholds
  bonferroni_threshold <- -log10(0.05 / nrow(results))
  suggestive_threshold <- -log10(suggestive_cutoff)
  
  # Create the Manhattan plot
  manhattan_plot <- ggplot(don, aes(x = BPcum, y = logp)) +
    geom_point(aes(color = as.factor(chromosome)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 22)) +
    scale_x_continuous(label = axisdf$chromosome, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(don$logp) + 2)) +
    geom_hline(yintercept = bonferroni_threshold, color = "#e7298a", linetype = "dashed") +
    geom_hline(yintercept = suggestive_threshold, color = "#78c679", linetype = "dashed") +
    geom_text(data = subset(don, logp > suggestive_threshold), aes(label = CpG), vjust = -1, size = 2) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x = "Chromosome", y = "-log10(p)")
  
  # Save the plot
  ggsave(filename = filename_save, plot = manhattan_plot, width = 14, height = 7, dpi = 300)
}


# Function to remove unreliable CpGs -------------------------------------------
exclude_CpGs <- function(cohort, cpgid, exclude, ethnic, artype='450K', filename = NULL, fileresume = NULL) {
  # Test if filter data exists and gets it
  if( toupper(artype) == '450K' ) {
    #try(data("filter_450K") )
    filters <- filter_450K
  }else if(toupper(artype) == 'EPIC') {
    #try(data("filter_EPIC"))
    filters <- filter_EPIC
  }else{
    stop( paste0( "Unknown array type ", toupper(artype) ) )
  }
  
  
  # In order to minimize data size we merge only the selected ethnic column
  
  
  # Get first ethnic column position in filters
  firstEthnicPosition <-  min(grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE))-1
  
  # Select ethnic columns not related with our data from filters (ONLY IF ethnic != '' or NA)
  if(ethnic!='' && !is.na(ethnic)) {
    # Get first ethnic column position in filters
    fieldstodelete <- grep(paste0("MASK_snp5_.*[^",ethnic,"]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    
  } else {
    #.Remove from code.# fieldstodelete <- grep(paste0("MASK_snp5_.*[^EUR]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
    fieldstodelete <- ""
  }
  
  
  fieldstomerge <- which(!seq(1:dim(filters)[2]) %in% fieldstodelete)
  
  # Merge cohort with CpGs filters
  cohort <- merge(cohort, filters[,fieldstomerge], by.x= cpgid, by.y = "probeID", all.x = TRUE )
  
  
  # CpGs id to remove
  if( is.null(exclude)) {
    excludeid <- NULL
  } else {
    excludeid <- cohort[eval(parse(text=getCritera(exclude, ethnic))), cpgid]
  }
  
  
  # Report descriptive exclussions to a descriptive file
  if(!is.null(fileresume)) {
    write(sprintf('\n# %s', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Remove "problematic"  CpGs : '), file = fileresume, append = TRUE)
    write(sprintf('# %s\n', strrep("-",16)), file = fileresume, append = TRUE)
    write(sprintf('# Criteria : \n\tArray type : %s \n\tEthnia : %s \n', toupper(artype), toupper(ethnic)), file = fileresume, append = TRUE)
    write(sprintf('# Mask : %s', exclude), file = fileresume, append = TRUE)
    write(sprintf('\n'), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = fileresume, append = TRUE)
    write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = fileresume, append = TRUE)
    write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = fileresume, append = TRUE)
    
    suppressWarnings(
      write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
                   filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  }
  
  # # Report CpG excluded and reason to a file
  # if(!is.null(filename)) {
  #    # write(sprintf('# Criteria : %s\n', toupper(artype)), file = filename)
  #    # write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = filename, append = TRUE)
  #    # write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = filename, append = TRUE)
  #    # write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = filename, append = TRUE)
  #    # write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = filename, append = TRUE)
  #    # Report exclusion reason
  #    suppressWarnings(
  #       write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
  #                    filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
  # }
  
  # Rmove CpGs with exclusion parameters
  cohort <- cohort[ !cohort[,cpgid] %in% excludeid,]
  
  warning( length(excludeid), ' CpGs have been excluded')
  
  return(cohort)
  
}


getCritera <- function(exclude, ethnic){
  
  criter = ''
  possible_crit <- c( 'MASK_sub25_copy', 'MASK_sub30_copy', 'MASK_sub35_copy', 'MASK_sub40_copy',
                      'MASK_mapping', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 'MASK_snp5_common', 'MASK_snp5_GMAF1p',
                      'MASK_general', 'cpg_probes', 'noncpg_probes', 'control_probes', 'Unrel_450_EPIC_blood', 'MASK_rmsk15',
                      'Sex', 'Unrel_450_EPIC_pla_restrict', 'Unrel_450_EPIC_pla', 'MASK_snp5_ethnic')
  
  if(!is.null(exclude[1]) && !is.na(exclude[1]) && exclude[1] != '') {
    # Test if all parameters are allowed
    if(length(which(! exclude %in% possible_crit))>=1)
      stop(paste0('Parameter(s) : ',paste(exclude[which(! exclude %in% possible_crit)], sep = ','),' not valid. Possible values are ',possible_crit ))
    
    # Get all exclude variables
    ##.. Works with lists -> convert list to matrix ..## exclusion.crit <- do.call(cbind, lapply( ls(patt="exclude"), get) )
    
    # Create formula with exclusion criteria
    ##.. Works with lists -> gets only parameters with value 'Exclude'..## criter <- paste(paste("cohort$",attributes(exclusion.crit)$dimnames[[1]][which(exclusion.crit=='Exclude')], sep = ""),"TRUE | ",sep=" == ", collapse = '')
    criter <- paste(paste("cohort$",exclude, sep = ""),"TRUE | ",sep=" == ", collapse = '')
    
    # If  MASK_general = 'Exclude' --> Change to :  ("MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p") == TRUE
    if( length(grep("cohort$MASK_general == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$MASK_general == TRUE", " cohort$MASK_sub30_copy == TRUE | cohort$MASK_mapping == TRUE | cohort$MASK_extBase == TRUE | cohort$MASK_typeINextBaseSwitch == TRUE | cohort$MASK_snp5_GMAF1p == TRUE " , criter, fixed = TRUE)
    
    # If Sex = 'Exclude' --> Change to : ( CpG_chrm %in% "chrX" | CpG_chrm %in% "chrY" )
    if( length(grep("cohort$Sex == TRUE",criter, fixed = TRUE))>0 )
      criter <- sub("cohort$Sex == TRUE", " cohort$CpG_chrm %in% 'chrX' | cohort$CpG_chrm %in% 'chrY' " , criter, fixed = TRUE)
    
    newCpGcond <- vector()
    if( length( grep("| cohort$cpg_probes == TRUE",criter, fixed = TRUE)) >0 )
      newCpGcond <- append(newCpGcond, c('cg'))
    if( length(grep("| cohort$noncpg_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('ch'))
    if( length(grep("| cohort$control_probes == TRUE",criter, fixed = TRUE))>0 )
      newCpGcond <- append(newCpGcond, c('rs'))  
    
    if(!is.null(vector())) {
      # Remove invalid criteria
      criter <- sub("| cohort$cpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$noncpg_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$control_probes == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("| cohort$MASK_snp5_ethnic == TRUE ", "", criter, fixed = TRUE)
      criter <- sub("cohort$cpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$noncpg_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$control_probes == TRUE |", "", criter, fixed = TRUE)
      criter <- sub("cohort$MASK_snp5_ethnic == TRUE |", "", criter, fixed = TRUE)
      
      # Add new criteria
      criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse = ''))
    }
    
    # Add ethnicity exlcusion criteria
    criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse= ''))
    if(is.na(ethnic) | ethnic == ''){
      criter <- paste("which(", substr(criter,0,nchar(criter)-3),")",sep="")
    }
  }
  
  if(!is.na(ethnic) & ethnic != ''){
    criter <- paste0(criter, " cohort$MASK_snp5_", ethnic," == TRUE | ")
    criter <- paste("which(", substr(criter, 0, nchar(criter)-2), ")", sep="")
  }
  criter
}


# Function to calculate lambda & confidence intervals
inflation <- function(pvector) {
  chisq <- qchisq(1 - pvector, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq))))
  lower_ci <-  lambda - (SE.median / qchisq(0.5, 1))
  upper_ci <-  lambda + (SE.median / qchisq(0.5, 1))
  lambda_ci <- c(lambda, lower_ci, upper_ci)
  return(lambda_ci)
}

# Function to do a qqplot norm
qq_norm_plot <- function(df, var_name,title=NULL) {
  
  vals <- df[[var_name]]
  o <- -log10(sort(vals,decreasing=F))
  e <- -log10( 1:length(o)/length(o))
  
  qqdata <- data.frame(o=o,e=e)
  
  ggplot(qqdata, aes(e,o)) +
    geom_point() +
    
    geom_abline(intercept=0,slope=1,col="#d7301f") + 
    
    # axis labels 
    xlab(expression(Expected~(-log[10](pval)))) + 
    ylab(expression(Observed~(-log[10](pval)))) +
    
    # styling
    theme_bw() +
    ggtitle(title) +
    theme(axis.title=element_text(face="bold",size=12),
          plot.title=element_text(face="bold",size=18,hjust=0))
}

# Function to add lambda and CI intervals to a qqplot
add_lambda_to_plot <- function(pvector, qqplot) {
  lambda <- inflation(pvector)
  # Add lamba + CI 
  lambda_CI <- sprintf("λ = %.2f [%.2f,%.2f]", lambda[1], lambda[2], lambda[3]) 
  print(lambda_CI)
  qqplot <- qqplot + annotate(geom = "text", x = -Inf, y = Inf,
                              hjust = -0.15, vjust = 1 + 0.15 * 3, label = lambda_CI, size = 5)
  return(qqplot)
}

#===============================================================================





#==============================================================================#
#                        LOAD & ARRANGE DATA                                   #----
#==============================================================================#

if (file_from == 'PACE') {
  load(ewas_file)
  ewas_result <- data.frame(alldataout$AdjustedwithCellType)
  colnames(ewas_result) <- gsub("^pval_.*", "pval", colnames(ewas_result))  # Keep 'pval' only
  colnames(ewas_result) <- gsub("^coef_.*", "coef", colnames(ewas_result))  # Keep 'coef' only
  colnames(ewas_result) <- gsub("^se_.*", "se", colnames(ewas_result))      # Keep 'se' only
  ewas_result <- na.omit(ewas_result)
  rm(alldataout)
}

if (file_from == 'meffil') {
  meffil <- readRDS(ewas_file)
  ewas_result <- meffil$analyses$none$table
  ewas_result$CpG <- rownames(ewas_result)
  colnames(ewas_result)[colnames(ewas_result) == 'p.value'] <- 'pval'
  colnames(ewas_result)[colnames(ewas_result) == 'coefficient'] <- 'coef'
  colnames(ewas_result)[colnames(ewas_result) == 'coefficient.se'] <- 'se'
  ewas_result <- na.omit(ewas_result)
}

##### we have p.value change to change them also with coef and se i think 


#==============================================================================#
#                                                                              #
#                            QUALITY CONTROL                                   #----
#                                                                              #
#==============================================================================#

if(array == "450K") {load(filter_450K_cpgs_path)}
if(array == "EPIC") {load(filter_EPIC_cpgs_path)}

ethnic <- 'GMAF1p'

# Prepare inputs to run EASIER functions 
exclude <-  c( 'control_probes', 'noncpg_probes', 'Sex',
               'MASK_mapping', 'MASK_sub30_copy', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 
               'MASK_snp5_common','MASK_snp5_ethnic',
               'Unrel_450_EPIC_pla_restrict')

# 'MASK_snp5_ethnic' -> remove SNPs common in European population 

cat('Number CpGs before exclusion', nrow(ewas_result)) #807201
file_path1 <- paste0("sva_pca/results/list_excluded_probes_",cohort,"_",current_date,".txt")
file_path2 <- paste0("sva_pca/results/descriptives_excluded_probes",cohort,"_",current_date,".txt")

ewas_result_filter <- exclude_CpGs(ewas_result, "CpG", 
                          exclude, 
                          ethnic = ethnic, 
                          filename = file_path1, 
                          fileresume = file_path2,
                          artype = array )

# started with #807201, 204072 CpGs have been excluded, now 603129 in BiSC
# started with 811982, 207480 CpGs exlcuded, no 604502 CpGs in INMA 


#==============================================================================#
#                                                                              #
#                    SIGNIFICANT AND SUGGESTIVE CpGs                           #----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                              BEFORE BACON                                    #----
#------------------------------------------------------------------------------#

# Filter CpG sites 
suggestive_cpgs <- subset(ewas_result_filter, pval < 1e-5) #120
ewas_result_filter$bf <- p.adjust(ewas_result_filter$pval, method = "bonferroni")
significant_cpgs <- ewas_result_filter[ewas_result_filter$bf < 0.05,] #18

cat('There are',nrow(suggestive_cpgs), 'suggestive CpGs in', analysis_name, 'model\n')
cat('There are',nrow(significant_cpgs), 'significant CpGs in', analysis_name, 'model\n')

# write.table(suggestive_cpgs, "sva_pca/results/basic_ewas_10e5.csv", row.names = FALSE, sep=";")

#------------------------------------------------------------------------------#
#                             AFTER BACON                                      #----
#------------------------------------------------------------------------------#

# Apply bacon adjustment to EWAS results
bc <- bacon::bacon(effectsizes = ewas_result_filter$coef,
                   standarderrors = ewas_result_filter$se,
                   na.exclude = TRUE) 

ewas_result_filter$pval_bacon <- bacon::pval(bc)
ewas_result_filter$bf_bacon <- p.adjust(ewas_result_filter$pval_bacon, method = "bonferroni")

# Filter CpG sites after bacon correction 
suggestive_cpgs_bacon <- subset(ewas_result_filter,  pval_bacon < 1e-5) 
significant_cpgs_bacon <- ewas_result_filter[ewas_result_filter$bf_bacon < 0.05,] 

cat('There are',nrow(suggestive_cpgs_bacon), 'suggestive CpGs after bacon correction in', analysis_name, 'model\n')
cat('There are',nrow(significant_cpgs_bacon), 'significant CpGs after bacon correction in', analysis_name, 'model\n')

saveRDS(ewas_result_filter, 
        paste0('sva_pca/results/',analysis_name, '_ewas_result_filtered_',cohort,'_',current_date, '.rds'))



#==============================================================================#
#                                                                              #
#                       QQ-lambda and MANHATTAN PLOTS                          #----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                              BEFORE BACON                                    #----
#------------------------------------------------------------------------------#

# Get QQ PLOT with lambda
qq_plot <- qq_norm_plot(ewas_result_filter, 'pval')
qq_lambda <- add_lambda_to_plot(ewas_result_filter$pval, qq_plot)
ggsave(paste0('sva_pca/results/QQplot_lambda_',analysis_name, '_ewas',cohort,'.png'), 
       plot = qq_lambda)



# GET MANHATTAN 
rownames(ewas_result_filter) <- ewas_result_filter$CpG
do_manhattan(results = ewas_result_filter, 
             array= array, 
             suggestive_cutoff = 1E-05,
             filename_save = paste0('sva_pca/results/Manhattan_',analysis_name, '_ewas',cohort,'.png'))


#------------------------------------------------------------------------------#
#                             AFTER BACON                                      #----
#------------------------------------------------------------------------------#

# Get QQ PLOT with lambda
qq_plot <- qq_norm_plot(ewas_result_filter, 'pval_bacon')
qq_lambda <- add_lambda_to_plot(ewas_result_filter$pval_bacon, qq_plot)
ggsave(paste0('sva_pca/results/bacon_QQplot_lambda_',analysis_name, '_ewas',cohort,'.png'), 
       plot = qq_lambda)

# GET MANHATTAN 
rownames(ewas_result_filter) <- ewas_result_filter$CpG
temp_ewas_result_filter <- ewas_result_filter
temp_ewas_result_filter$pval <- ewas_result_filter$pval_bacon
do_manhattan(results = temp_ewas_result_filter, 
             array = array, 
             suggestive_cutoff = 1E-05,
             filename_save = paste0('sva_pca/results/bacon_Manhattan_',analysis_name, '_ewas',cohort,'.png'))

rm(temp_ewas_result_filter)



#==============================================================================#
#                                                                              #
#                       COMPARISON WITH REFERENCE STUDY                        #----
#                                                                              #
#==============================================================================#

# Reference Study. 
# Everson, T.M., Vives-Usano, M., Seyve, E. et al. 
# Placental DNA methylation signatures of maternal smoking during pregnancy and 
# potential impacts on fetal growth. Nat Commun 12, 5095 (2021). 
# https://doi.org/10.1038/s41467-021-24558-y

ref_study <- readxl::read_excel('sva_pca/results/SuppData8_MetaRes_AfterBACON_WithAnnotations.xlsx')
ref_study <- as.data.frame(ref_study)
colnames(ref_study) <- ref_study[1:1,]

# Get CpGs for any+sust smoking during pregnancy
ref_study_cpgs <- na.omit(ref_study$CpG)[2:444] #443


#------------------------------------------------------------------------------#
#                              BEFORE BACON                                    #----
#------------------------------------------------------------------------------#
case_study_cpgs <- significant_cpgs$CpG
top_cpgs <- list(ref_study_cpgs, case_study_cpgs)

# Make venn diagram 
sink("/dev/null")

# Create the Venn diagram for two sets
venn_counts <- venn.diagram(
  x = top_cpgs,
  category.names = c("Reference Study", "Case Study"),
  filename = NULL,  # Set to NULL for further customization before saving
  output = TRUE,
  
  # Appearance customization for 2 sets
  fill = c("#D9D9D9", "#B0B0B0"),  # Grayscale colors for two groups
  alpha = 0.5,
  cex = 1.5,  # Size for the counts
  fontface = "bold",
  fontfamily = "serif",
  cat.cex = 1.5,  # Size for category labels
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25),  # Adjust positions of the category names for two sets
  cat.dist = c(0.025, 0.025),  # Bring category labels closer to the circles
  cat.fontfamily = "serif",
  lwd = 2,  # Line width of the circles
  lty = 'solid',  # Line type of the circles
  col = c("#404040", "#606060"),  # Grayscale border colors for two sets
  
  # Display counts only
  print.mode = c("raw"),  # Show raw counts in the intersections
  label.col = "black",
  label.fontface = "bold",
  label.cex = 1.2,
  label.just = "center",
  
  # Other customizations
  margin = 0.1  # Adjust margin to fit long text if necessary
)

# Save the Venn diagram
png(paste0('sva_pca/results/Venn_RefCpGs_',analysis_name, '_ewas',cohort,'.png'), 
    width = 800, height = 800)
grid.draw(venn_counts)
dev.off()
sink()

# Get CpG IDs in the intersection 
intersect_ref_case <- intersect(ref_study_cpgs, case_study_cpgs)
cat(length(intersect_ref_case),'CpGs in BOTH case and reference study:',intersect_ref_case, "\n")

# Get the CpG sites only present in the case study (not in the reference study)
only_case_study <- setdiff(case_study_cpgs, ref_study_cpgs)
cat(length(only_case_study),'CpGs in ONLY CASE study:',only_case_study, "\n")

#------------------------------------------------------------------------------#
#                             AFTER BACON                                      #----
#------------------------------------------------------------------------------#
case_study_cpgs <- na.omit(significant_cpgs_bacon$CpG) #8 or 65

top_cpgs <- list(ref_study_cpgs, case_study_cpgs)

# Make venn diagram 
sink("/dev/null")

# Create the Venn diagram for two sets
venn_counts <- venn.diagram(
  x = top_cpgs,
  category.names = c("Reference Study", "Case Study"),
  filename = NULL,  # Set to NULL for further customization before saving
  output = TRUE,
  
  # Appearance customization for 2 sets
  fill = c("#D9D9D9", "#B0B0B0"),  # Grayscale colors for two groups
  alpha = 0.5,
  cex = 1.5,  # Size for the counts
  fontface = "bold",
  fontfamily = "serif",
  cat.cex = 1.5,  # Size for category labels
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25),  # Adjust positions of the category names for two sets
  cat.dist = c(0.025, 0.025),  # Bring category labels closer to the circles
  cat.fontfamily = "serif",
  lwd = 2,  # Line width of the circles
  lty = 'solid',  # Line type of the circles
  col = c("#404040", "#606060"),  # Grayscale border colors for two sets
  
  # Display counts only
  print.mode = c("raw"),  # Show raw counts in the intersections
  label.col = "black",
  label.fontface = "bold",
  label.cex = 1.2,
  label.just = "center",
  
  # Other customizations
  margin = 0.1  # Adjust margin to fit long text if necessary
)

# Save the Venn diagram
png(paste0('sva_pca/results/bacon_Venn_RefCpGs_',analysis_name, '_ewas',cohort,'.png'), 
    width = 800, height = 800)
grid.draw(venn_counts)
dev.off()
sink()

# Get CpG IDs in the intersection 
intersect_ref_case <- intersect(ref_study_cpgs, case_study_cpgs)
cat(length(intersect_ref_case), 'CpGs in BOTH case and reference study:',intersect_ref_case, "\n")

# Get the CpG sites only present in the case study (not in the reference study)
only_case_study <- setdiff(case_study_cpgs, ref_study_cpgs)
cat(length(only_case_study), 'CpGs in ONLY CASE study:',only_case_study, "\n")









#------------------------------------------------------------------------------#
#                    Check results with EWAS catalogue                         #----
#------------------------------------------------------------------------------#

# ewascat <- data.frame()
# 
# for (i in seq_along(1:length(case_study_cpgs))) {
#   cpg <- case_study_cpgs[i]
#   print(cpg)
#   
#   res <- tryCatch({
#     ewascatalog(cpg, "cpg")
#   }, error = function(e) {
#     NULL
#   })
#   
#   # Check if the result is empty or not found
#   if (!is.null(res) && nrow(res) > 0) {
#     res$CpG <- cpg
#     ewascat <- rbind(ewascat, res)
#   } else {
#     cat("No results found for CpG:", cpg, "\n")
#   }
# }
# 
# 
# # 21 - CpG
# # 15 - tissue 
# # 3 - pmid
# # 9 - outcome
# # 10 - exposure
# # 27 - beta
# # 28 - se
# # 29 - p
# 
# ewascat  <- ewascat[, c(21,15,3,9,10,27,28,29)]
# 
# write.table(ewascat,"sva_pca/results/ewascat_basic_ewas_sig_bf_bacon.csv",
#             sep=";", row.names = FALSE, col.names = TRUE)
# 
# 
# rm(basic, basic_filter, ewascat, bc)

#===============================================================================







