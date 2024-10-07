#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 6: Compare results                                #
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



# [INPUT]: All models/files we want to compare. 
# - RDS with fitlered results after 5_Analyze_EWAS_results 

# [OUPUT]:
# - Venn plots with intersection for all input files. 
# - Table with CpGs BF < 0.05 in any of the models, table with Pval and coeff 

# [CONSIDERATIONS]
#===============================================================================

#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# INMA -> now running :) 

# BISC
basic <- 'sva_pca/results/basic_ewas_result_filtered_BISC_21092024.rds'
opt6pca <- 'sva_pca/results/opt6PCs_ewas_result_filtered_BISC_21092024.rds'
opt22pca <- 'sva_pca/results/opt22PCs_ewas_result_filtered_BISC_21092024.rds'
opt6sva <- 'sva_pca/results/opt6SV_ewas_result_filtered_BISC_21092024.rds'
opt136mef <- 'sva_pca/results/meffil_136sv_ewas_result_filtered_BISC_21092024.rds'

# INMA 
# basic <-  'sva_pca/results/basic_ewas_result_filtered_INMA_21092024.rds'
# opt5pca <- 'sva_pca/results/opt5pcs_ewas_result_filtered_INMA_21092024.rds'
# opt16pca <- 'sva_pca/results/opt16pcs_ewas_result_filtered_INMA_21092024.rds'
# opt7sva <- 'sva_pca/results/opt7svs_ewas_result_filtered_INMA_21092024.rds'
# - meffil 

list_ewas_files <- c(basic, opt6pca, opt22pca, opt6sva, opt136mef)
#list_ewas_files <- c(basic, opt5pca, opt16pca, opt7sva)

# Define names to be used in the venn plot 
list_names <- c("basic", "opt6pca", "opt22pca", "opt6sva", "opt136mef")
#list_names <- c("basic", "opt5pca", "opt16pca", "opt7sva")

cohort <- 'BISC'
current_date <- '21092024'

cutoff <- 0.05 # to get sumamry statsitcs table
#==============================================================================#
#                                                                              #
#                        INIT WORKING ENVIORNEMT                               #
#                                                                              #
#==============================================================================#

library(ggplot2)
library(data.table)
library(ewascatalog)
library(VennDiagram)
library(grid)
library(dplyr)
library(UpSetR)

#==============================================================================#
#                        LOAD & ARRANGE DATA                                   #----
#==============================================================================#
# Colnames = CpG, coef, se, pval, bonferroni_bacon, pval_bacon, bonferroni

list_ewas_results <- lapply(list_ewas_files, function(file) {readRDS(file)})

# Read the files and assign names to the resulting list
list_ewas_results <- setNames(lapply(list_ewas_files, readRDS), list_names)
list_cpgs <- lapply(list_ewas_results, function(ewas){ewas[ewas$bf < 0.05, "CpG"]})


#==============================================================================#
#                                                                              #
#                     VENN PLOT CpGs - MODELS                                  #----
#                                                                              #
#==============================================================================#
#[NOTE]: Maxmum 5 elements in list/models if not venn diagrama doesn't work 

# Make Venn plot 
sink("/dev/null")

# Generate grayscale colors based on the number of sets in list_cpgs
num_sets <- length(list_cpgs)  # Number of sets
grayscale_fill <- gray.colors(num_sets, start = 0.8, end = 0.3)  # Light to dark gray
grayscale_border <- gray.colors(num_sets, start = 0.5, end = 0.1)  # Darker grays for borders


# Create the Venn diagram
venn_counts <- venn.diagram(
  x = list_cpgs,  # Assuming this is a named list of CpGs for each model
  category.names = names(list_cpgs),
  filename = NULL,  # Set to NULL for further customization before saving
  output = TRUE,
  
  #Dynamically generated grayscale colors
  fill = grayscale_fill,  # Grayscale colors for the sets
  alpha = 0.5,
  cex = 1.5,  # Size for the counts
  fontface = "bold",
  fontfamily = "serif",
  cat.cex = 1.5,  # Size for category labels
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = seq(-25, 25, length.out = num_sets),  # Adjust category name positions dynamically
  cat.dist = rep(0.025, num_sets),  # Bring category labels closer to the circles
  cat.fontfamily = "serif",
  lwd = 2,  # Line width of the circles
  lty = 'solid',  # Line type of the circles
  col = grayscale_border,  # Grayscale border colors for the sets
  
  # Display counts only
  print.mode = c("raw"),  # Show raw counts in the intersections
  label.col = "black",
  label.fontface = "bold",
  label.cex = 1.2,
  label.just = "center",
  
  # Other customizations
  margin = 0.1  # Adjust margin to fit long text if necessary
  )

# Save the Venn diagram to a PNG file
png(paste0('sva_pca/results/Venn_RefCpGs_ewas', cohort, '.png'), 
    width = 800, height = 800)
grid.draw(venn_counts)
dev.off()

# Re-enable output to the console
sink()



#==============================================================================#
#                                                                              #
#                     UPSET CpGs - MODELS - Reference                          #----
#                                                                              #
#==============================================================================#


# Load Reference Study 
ref_study <- readxl::read_excel('sva_pca/results/SuppData8_MetaRes_AfterBACON_WithAnnotations.xlsx')
ref_study <- as.data.frame(ref_study)
colnames(ref_study) <- ref_study[1:1,]

# Get CpGs for any+sust smoking during pregnancy
ref_study_cpgs <- na.omit(ref_study$CpG)[2:444] #4

# Get list of CpGs for each model + reference
list_cpgs <- c(list_cpgs, ref_study = list(ref_study_cpgs))

# Convert list_cpgs to a data frame suitable for UpSet plot
upset_input <- fromList(list_cpgs)
png_filename <- paste0('sva_pca/results/UpSet_RefCpGs_ewas', cohort,'_', current_date, '.png')
png(png_filename)    
upset(upset_input)
dev.off()



#==============================================================================#
#                                                                              #
#                     SUMMARY TABLES CpGs - MODELS                             #----
#                                                                              #
#==============================================================================#


#------------------------------------------------------------------------------#
#                        Get summary tables CpGs  bef bacon                    #----
#------------------------------------------------------------------------------#


stats_colnames <- c('CpG', 'coef', 'pval', 'pval_bacon')
#stats_colnames <- c('CpG', 'coef', 'se', 'pval', 'bf')


# Extract top hits CpGs for each file base on suggestive cutoff
top_cpgs <- lapply(list_ewas_results, 
                        function(ewas) {subset(ewas, bf <= cutoff)$CpG})

# Get a unique list of all top CpGs
unique_top_cpgs <- unique(unlist(top_cpgs))

# Create a summary table for each EWAS result with the specified columns
summary_tables_main <- lapply(seq_along(list_ewas_results), function(i) {
  ewas <- list_ewas_results[[i]]
  summtab <- ewas[ewas$CpG %in% unique_top_cpgs, stats_colnames]
  colnames(summtab) <- c('CpG', paste0('coef_model', i), 
                         #paste0('se_model', i), 
                         paste0('pval_model', i)) 
                         #paste0('bf_model', i))
  return(summtab)
})

# Merge all summary tables by CpG
combined_summary_table_main <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), 
                                      summary_tables_main)

# Order by CpG
combined_summary_table_main <- combined_summary_table_main[order(combined_summary_table_main$CpG), ]

# Save to CSV file 
save_to <- paste0('sva_pca/results/summtab_stat_BFcpg',cohort,'_',current_date,'.csv')

write.table(combined_summary_table_main, 
            file = save_to, sep = ";", row.names = FALSE)






#==============================================================================#
#                                                                              #
#                 SUMMARY TABLES STATISTICS - MODELS                           #----
#                                                                              #
#==============================================================================#


# Function to calculate summary statistics
get_summary_stats <- function(df, colname) {
  stats <- c(
    Min = min(df[[colname]], na.rm = TRUE),
    Q1 = quantile(df[[colname]], 0.25, na.rm = TRUE),
    Median = median(df[[colname]], na.rm = TRUE),
    Mean = mean(df[[colname]], na.rm = TRUE),
    Q3 = quantile(df[[colname]], 0.75, na.rm = TRUE),
    Max = max(df[[colname]], na.rm = TRUE)
  )
  return(stats)
}


#------------------------------------------------------------------------------#
#                        For All CpGs in All models                            #----
#------------------------------------------------------------------------------#


# SLOPE stats for each df of 
slope_stats <- lapply(list_ewas_results_main, function(ewas) {get_summary_stats(ewas, 'coef')})
slope_stats <- data.frame(slope_stats)
colnames(slope_stats) <- c('slope_basic', 'slope_sva', 'slope_pca', 'slope_meffil')


# PVAL  stats for each df of 
pval_stats <- lapply(list_ewas_results_main, function(ewas) {
  ewas$log_pval <-  -log10(ewas$pval)
  get_summary_stats(ewas, 'log_pval')
})

pval_stats <- data.frame(pval_stats)
colnames(pval_stats) <- c('pval_basic', 'pval_sva', 'pval_pca', 'pval_meffil')

# DIST   stats for each df of 

dist_stats <- lapply(list_ewas_results_main, function(ewas) {
  ewas$dist <- ewas$CpG_end - ewas$CpG_beg
  get_summary_stats(ewas, 'dist')
})

dist_stats <- data.frame(dist_stats)
colnames(dist_stats) <- c('dist_basic', 'dist_sva', 'dist_pca', 'dist_meffil')


all_stats <- cbind(slope_stats, pval_stats, dist_stats)
write.table(all_stats, 'sva_pca/results/stats_all_cpgs_all_models.csv')



