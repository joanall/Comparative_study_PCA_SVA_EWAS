#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 2: Perform SVA or PCA                             #
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
# - Residuals from EWAS saved in Large MAtrix format. Samples (rownames) x CpGs (colnames) 
# - Covariables data in DATAFRAME format (.rds), where samples x covariables

# Covariables/Metadata dataframe needs a colum named "Basename" with the
# samples IDs, the samples IDs have to match the columnames of Methylation data


# [OUTPUT]:
# - PCs inside Large prcomp object 
# - SVs as LAreg lsit object, returned by sva() 

# [CONSIDERATIONS]: SVA() may fail if variables included in the matrix are highly
# correlated. Uncomment the section code named "Detection of co-linearility" to 
# know which variables are causing the problem. In my case, BiSC cohort didn't
# have this problem but in INMA cohort CellTypes Proporiton were making the model
# fail. In that casewe advise you to remove the Stromal. 
# It is very important to keep Syncytiotrophoblasts and trophoblasts. 

#===============================================================================

#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# Paths to methylation and metadata
# residuals_matrix <- 'sva_pca/results/residuals_variance/residuals_matrix_21092024.rds' # bisc
# metadata_path <- "sva_pca/db/sva_pca_metadata_20092024.rds"

residuals_matrix <- 'sva_pca/results/residuals_matrix_INMA_03102024.rds' #inma
metadata_path <- "sva_pca/db/INMA_sva_pca_metadata_03102024.rds"

current_date <- "03102024"
cohort_name <- 'INMA'

# Covars to include in the model matrix for computing sva() 
covars <- c('Smoke', 'estudios3c', 'ethnic_origin2c',
  'Cohort', 'sges', 'Sex', 'Trophoblasts', 'Syncytiotrophoblast', 
  'Hofbauer', 'Endothelial', 'nRBC')
# in bisc we do add Stromal 

#===============================================================================



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
#library(sfsmisc)
library(stats)
library(readr)
library(Matrix)
library(sva)

residuals_matrix <- readRDS(residuals_matrix)


#==============================================================================#
#                                                                              #
#                    PRINCIPAL COMPONENT ANALYSIS (PCA)                        #----
#                                                                              #
#==============================================================================#

cat('---------------------- Get residuals PCs --------------------------------\n')

start_time <- Sys.time()

# Perform PCA on the residuals (samples x CpGs)
pca_residuals <- prcomp(residuals_matrix, scale = TRUE)
# sdev, rotation, center, scale,x

write_rds(pca_residuals,
          paste0('sva_pca/results/pca_residuals_',cohort_name, '_', current_date,'.rds'))

cat('Time for running PCA', Sys.time() - start_time, "\n")



#==============================================================================#
#                                                                              #
#                      SINGULAR VALUE DECOMPOSITION                            #----
#                                                                              #
#==============================================================================#


cat('---------------------- Get residuals SVs --------------------------------\n')

metadata <- readRDS(metadata_path) # samples x covars
# residual matrix: (CpGs as columns, samples as rows)

#------------------------------------------------------------------------------#
#                             Prepare data                                     #----
#------------------------------------------------------------------------------#

# Remove metadata rows with missing values
metadata <- na.omit(metadata)  

# Ensure that all samples in metadata are in meth_data
all(metadata$Basename %in% rownames(residuals_matrix))  # Should return TRUE

# Ensure that all samples in meth_data are in metadata
all(rownames(residuals_matrix) %in% metadata$Basename)  # Should return TRUE

# Subset meth_data to include only samples present in metadata
residuals_matrix <- residuals_matrix[rownames(residuals_matrix) %in% metadata$Basename,]

# Match and order samples in metadata with meth_data
metadata <- metadata[match(rownames(residuals_matrix), metadata$Basename), ]

# Check match-order again
all(metadata$Basename %in% rownames(residuals_matrix))  # Should return TRUE
all(rownames(residuals_matrix) %in% metadata$Basename)  # Should return TRUE

cat('Residuals data', dim(residuals_matrix), '\n') # samples x CpgS
cat('Metadata', dim(metadata), '\n')   # 553 samples x 15 covars


#------------------------------------------------------------------------------#
#                             Get Surrogates                                   #----
#------------------------------------------------------------------------------#

# Prepare model matrices
# mod <- model.matrix(~ maternal_smoking + Sex + maternal_educational_level +
#                       gestational_age + study_center + covid_confinement +
#                       child_ethnicity + Trophoblasts + Stromal + Hofbauer +
#                       Endothelial + nRBC + Syncytiotrophoblast,
#                       data = metadata)

metadata <- metadata[,covars]

mod <- model.matrix(~ Smoke + estudios3c + ethnic_origin2c +
                      Cohort + sges + sex  + Trophoblasts + Syncytiotrophoblast + Hofbauer +
                      Endothelial + nRBC,
                      data = metadata)


mod0 <- model.matrix(~ 1, data = metadata)

cat('Matrices defined, starting transposing methylation... \n')

# Prepare residuals matrix
sparse_residuals_matrix <- Matrix(residuals_matrix, sparse = TRUE)
t_residuals_matrix <- t(sparse_residuals_matrix)
rownames(t_residuals_matrix) <- colnames(residuals_matrix)
colnames(t_residuals_matrix) <- rownames(residuals_matrix)
t_residuals_matrix <- as.matrix(t_residuals_matrix)

cat('Tranposing matrix finished...\n')

# Run SVA 
# needs to be CpGs x samples and samples x vars
start_time <- Sys.time()

# Alternative: n.sv = num.sv (data, mod,method = "leek") anf vie n.sv as paraaet in sva()
# be default it also comutes de optimal n sv though it uses "be" instead of "leek" i think

cat('starting actual sva')
sv_obj <- sva(t_residuals_matrix, mod, mod0)  # Use imputed methylation data

write_rds(sv_obj,
          paste0('sva_pca/results/sva_woStromal_residuals',cohort_name,'_', current_date,'.rds'))

cat('Time for running SVA', Sys.time() - start_time, "\n")


# [UNCCOMENT IF SVA FAILS, CHECK MULTICOLINEARITY]

# alias(lm(as.matrix(mod) ~ ., data = metadata))
# alias_info <- alias(lm(as.matrix(mod) ~ ., data = metadata))
# print(alias_info$Complete)

# If one or more variables are multicolinear remove them. Run again the residuals
# wo including them. Run again PCA and SVA without including them. 


