#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 4: Run EWAS meffil                                #
#                                                                              #
#  <Display description>                                                       #----
#==============================================================================#

# Objective: Assess whether performing PCA and SVA on EWAS data and adding PCs 
#or SVs on the models changes the results by truly removing unknonw variation effects
# and therefore uncovering stronger/clearer biological signals. 

# Meffil perfroms sva() by default and chooses the number of optimal SVs using 
# croos validtion with 10 folds 


# [INPUT]:
# - Methylation data in a dataframe samples x cpgs 
# - Covariables data in DATAFRAME format (.rds), where samples x covariables

# Covariables/Metadata dataframe needs a colum named "Basename" with the
# samples IDs, the samples IDs have to match the columnames of Methylation data

# [OUTPUT]:


#==============================================================================#
#                           USER PARAMETERS                                    #
#==============================================================================#

# Paths to methylation and metadata
# methylation_path <- "db/bVals_BISC_nosib_noout_nomask.RDS"
# metadata_path <- "sva_pca/db/sva_pca_metadata_20092024.rds"
# exposure <- 'maternal_smoking'
# covars <- c("Sex", "maternal_educational_level", "gestational_age",
#             "study_center", "covid_confinement", "child_ethnicity")
# cell_type_vars <- c('Trophoblasts','Stromal', 'Hofbauer', 'Endothelial','nRBC', 
#                     'Syncytiotrophoblast')


methylation_path <- "sva_pca/db/bVals_INMA_noout.RDS" #inma
metadata_path <- "sva_pca/db/INMA_sva_pca_metadata_03102024.rds" #inma

exposure <- 'Smoke'
covars <- c("estudios3c", "ethnic_origin2c", "Cohort", "sges", "Sex")  #INMA
            cell_type_vars <- c('Trophoblasts', 'Hofbauer', 'Endothelial','nRBC',
                    'Syncytiotrophoblast')


cohort <- 'INMA'

current_date <- '21092024'

#==============================================================================#
#                                                                              #
#                         INIT  WORKING ENVIRONMENT                            #----
#                                                                              #
#==============================================================================#

print(paste0('YOUR WORKING DIRECTORY IS:', getwd()))

# Load libraries 
library(meffil) # ewas
library(tableone)
library(dplyr)
library(parallel)
library(readr)




#------------------------------------------------------------------------------#
#                         DATA PREPARATION                                     #----
#------------------------------------------------------------------------------#


# Load methylation and metadata
meth_data <- readRDS(methylation_path)  # CpGs x samples
metadata <- readRDS(metadata_path)      # samples x covariates

# Remove missing values from metadata
metadata <- na.omit(metadata)

# Subset methylation data to match the samples in metadata
meth_data <- meth_data[, colnames(meth_data) %in% metadata$Basename]

# Match and order metadata with the methylation data
metadata <- metadata[match(colnames(meth_data), metadata$Basename), ]

# Extract exposure variable
exposure_vector <- metadata[[exposure]]  # This should be a column in metadata

# Check if lengths match
if (ncol(meth_data) != length(exposure_vector)) {
  stop("Number of samples in methylation data and exposure variable do not match.")
}

# Create covariates data frame
covariates_data <- metadata[, covars]
covariates_data <- cbind(covariates_data, metadata[, cell_type_vars])

#==============================================================================#
#                                                                              #
#                                RUN MEFFIL EWAS                               #                  
#                                                                              #
#==============================================================================#

ewas_result <- meffil.ewas(
  beta = meth_data,                   # Methylation matrix
  variable = exposure_vector,         # Exposure variable vector
  covariates = covariates_data,       # Covariates data frame
  sva = TRUE,                         # Surrogate Variable Analysis
  rlm = TRUE                          # Robust linear model
)

# Save results
write_rds(ewas_result, 
          file = paste0('sva_pca/results/meffil_ewas_',cohort,'_',current_date,'.rds'))

# Print a summary of the EWAS result
print(summary(ewas_result))

