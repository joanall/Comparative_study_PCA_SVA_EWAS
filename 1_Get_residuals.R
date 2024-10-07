#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 1: Get Residuals                                  #
#                                                                              #
#  <Display description>                                                       #----
#==============================================================================#

# Objective: Assess whether performing PCA and SVA on EWAS data and adding PCs 
# or SVs on the models changes the results by truly removing unknown variation 
# effects and therefore uncovering stronger/clearer biological signals. 

# STEP 1: Get Residuals to remove the effect of known covariates from DNAm 
# Residuals represent the part of DNAm not explained by the covariates. 
#   Methylation ∼ Exposure + Covariates

# STEP 2: SVA or PCA on Residuals to identify hidden patterns of variation- 

# STEP 3: Choose optimal SVs or PCs ("optional") by applying elbow method

# STEP 4: Running EWAS with and without SVs or PCs
#   Methylation ∼ Exposure + Covariates +  SVs/PCs

# STEP 5: Compare hits and statistics within approaches and with reference study.


# [INPUT]: 
# - Methylation data in DATAFRAME format (.rds), where CpGs (rownames) x samples (colnames)
# - Covariables data in DATAFRAME format (.rds), where samples x covariables
# - Covariables/Metadata dataframe needs a colum named "Basename" with the
# samples IDs, the samples IDs have to match the columnames of Methylation data


# [OUTPUT]:
# - A Large MAtrix object with the resiuals after performing the lm, where rownames 
#  are samples and colnames are CpGs


# [CONSIDERATIONS]: SVA() can't handle missing data, so we apply imputation
# to our methylation data (instead of removing that CpG). 

#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# Paths to methylation and metadata
# methylation_path <- "sva_pca/db/bVals_BISC_nosib_noout_nomask.RDS" #bisc
# metadata_path <- "sva_pca/db/sva_pca_metadata_20092024.rds" #bisc
methylation_path <- "sva_pca/db/bVals_INMA_noout.RDS" #inma
metadata_path <- "sva_pca/db/INMA_sva_pca_metadata_03102024.rds" #inma


#exposure <- 'maternal_smoking' #bisc
exposure <- 'Smoke' #inma


# covars <- c("Sex", "maternal_educational_level", "gestational_age", # bisc
#             "study_center", "covid_confinement", "child_ethnicity")
 
covars <- c("estudios3c", "ethnic_origin2c", "Cohort", "sges", "Sex")  #inma


cell_type_vars <- c('Trophoblasts', 'Hofbauer', 'Endothelial','nRBC',
                    'Syncytiotrophoblast')

cell_type_vars <- c()
# we have removed the stromal 
# Get curent date
current_date <- "03102024"
cohort_name <- 'INMA'



#==============================================================================#
#                                                                              #
#                         INIT  WORKING ENVIRONMENT                            #----
#                                                                              #
#==============================================================================#

#rm(list=ls())
#.rs.restartR()
#gc()

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))


# Load libraries 
library(dplyr)
#library(sfsmisc)R
library(readr)


#------------------------------------------------------------------------------#
#                              Define functions                                #----
#------------------------------------------------------------------------------#

# Function to impute (copied from meffil package)
impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
  if (margin == 2) x <- t(x)
  
  idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[,"row"])
    v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
    v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
    x[idx] <- v[match(idx[,"row"],na.idx)] ##
    stopifnot(all(!is.na(x)))
  }
  
  if (margin == 2) x <- t(x)
  x
}



#==============================================================================#
#                                                                              #
#                            DATA PREPARATION                                  #----
#                                                                              #
#==============================================================================#

cat('-------------------------- Prepare data ----------------------------------')

meth_data <- readRDS(methylation_path)  # Cpgs x samples
metadata <- readRDS(metadata_path) # samples x covars

#------------------------------------------------------------------------------#
#                               Imputation                                     #----
#------------------------------------------------------------------------------#

cat('Imputation of methylation values ...')
start_time <- Sys.time()
meth_data <- impute.matrix(meth_data)
write_rds(meth_data,
          paste0('sva_pca/results/imputed_methylation_',cohort_name, '_', current_date, '.rds'))
cat('Time for imputing meth data', Sys.time() - start_time, "\n")


#------------------------------------------------------------------------------#
#                               Match and Order                                #----
#------------------------------------------------------------------------------#

# Ensure that all samples in metadata are in meth_data
all(metadata$Basename %in% colnames(meth_data))  # Should return TRUE

# Ensure that all samples in meth_data are in metadata
all(colnames(meth_data) %in% metadata$Basename)  # Should return TRUE

# Subset meth_data to include only samples present in metadata
meth_data <- meth_data[,colnames(meth_data) %in% metadata$Basename]

# Match and order samples in metadata with meth_data
metadata <- metadata[match(colnames(meth_data), metadata$Basename), ]

# Check match-order again 
all(metadata$Basename %in% colnames(meth_data))  # Should return TRUE
all(colnames(meth_data) %in% metadata$Basename)  # Should return TRUE

cell_type_data <- as.matrix(metadata[,cell_type_vars])
rownames(cell_type_data) <- metadata$Basename

cat('Methylation data', dim(meth_data), '\n') #5000 CpGs x 553 samples 
cat('Metadata', dim(metadata), '\n')   # 553 samples x 15 covars 


#==============================================================================#
#                                                                              #
#                           RESIDUALIZATION                                    #----
#                                                                              #
#==============================================================================#

cat('--------------------- Start residualization -----------------------------')

# Prepare residuals matrix to store the residuals
residuals_matrix <- matrix(NA, nrow = ncol(meth_data), ncol = nrow(meth_data))  # Samples x CpGs
rownames(residuals_matrix) <- colnames(meth_data)  # Match samples
colnames(residuals_matrix) <- rownames(meth_data)  # Match CpGs

# Covariates for the model
exp_cov <- c(exposure, covars)
covariates <- metadata[, exp_cov]


start_time <- Sys.time()
# Loop through each CpG site to fit linear models and extract residuals
for (i in 1:nrow(meth_data)) {
  
  if(i/10 == 0) {print(i)}
  # Get methylation values for the CpG site (transpose, CpG i for all samples)
  y <- meth_data[i, ]
  
  # Fit the model: Methylation ~ Covariates (adjusting for covariates)
  model <- lm(y ~ ., data = cbind(covariates, y = y))
  
  # Store the residuals
  residuals_matrix[, i] <- residuals(model)
}


write_rds(residuals_matrix, 
          paste0('sva_pca/results/residuals_matrix_noCT', cohort_name, '_', current_date, '.rds'))
cat('Residualization complete. Residuals matrix dimensions:', dim(residuals_matrix), '\n')
cat('Time for obtaining the residuals', Sys.time() - start_time, "\n")

# Rownames are samples and colnames are CpGs





