#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 4:  Run EWAS                                      #
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
# - Methylation data in a dataframe samples x cpgs 
# - Covariables data with PCs/SVs in DATAFRAME format (.rds), where samples x covariables

# Covariables/Metadata dataframe needs a colum named "Basename" with the
# samples IDs, the samples IDs have to match the columnames of Methylation data

# [OUTPUT]:


#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# Paths to methylation and metadata
# methylation_path <- "db/bVals_BISC_nosib_noout_nomask.RDS"       # BISC 
# metadata_path <- "sva_pca/results/metadata_22pcs_21092024.rds"
 
methylation_path <- "sva_pca/results/imputed_methylation_INMA_03102024.rds"       # INMA
metadata_path <- "sva_pca/db/INMA_sva_pca_metadata_21092024.rds"


save_to <- '/PROJECTES/BISC_OMICS/analyses/BiSC_23/EWAS_diet_CL/sva_pca/results'

exposure <- 'Smoke'
pcs <- paste0("SV", 1:14)

# covars <- c("Sex", "maternal_educational_level", "gestational_age",   # BISC
#             "study_center", "covid_confinement", "child_ethnicity",
#             pcs)

covars <- c("Smoke", "estudios3c", "ethnic_origin2c", "Cohort", "sges", "Sex")  #INMA
           # pcs)

cell_type_vars <- c('Trophoblasts', 'Hofbauer', 'Endothelial','nRBC', 
                    'Syncytiotrophoblast')

# Get curent date
current_date <- "21092024"
analysis_name <- 'basic_ewas'
cohort <- 'INMA'

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
library(tableone)
library(dplyr)
#library(sfsmisc)
library(parallel)
library(readr)
library(PACEanalysis)

#------------------------------------------------------------------------------#
#                              Define functions                                #----
#------------------------------------------------------------------------------#



#==============================================================================#
#                                                                              #
#                            DATA PREPARATION                                  #----
#                                                                              #
#==============================================================================#

cat('-------------------------- Prepare data ----------------------------------')

meth_data <- readRDS(methylation_path) # Cpgs x samples
metadata <- readRDS(metadata_path) # samples x covars

cell_type_data <- as.matrix(metadata[,cell_type_vars])
rownames(cell_type_data) <- metadata$Basename


if(is.factor(metadata[[exposure]]) == FALSE) {metadata[[exposure]] <- as.factor(metadata[[exposure]])}
cat(colnames(metadata))
tempresults <- dataAnalysis(phenofinal=metadata,
                            betafinal=meth_data,
                            array='EPIC',
                            maxit=100,
                            Omega=cell_type_data, #matrix -> 6 Cell Type proportions (planet package) -> 565 samples x 6 CT
                            vartype= "ExposureCat",
                            robust=TRUE,
                            varofinterest=exposure, # main exposure
                            
                            Table1vars=covars,
                            StratifyTable1=FALSE,
                            StratifyTable1var=NULL,
                            
                            #which model
                            adjustmentvariables=covars, 
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=FALSE,
                            RunCellTypeInteract=FALSE,
                            
                            # susbet
                            RestrictToSubset=FALSE,
                            RestrictionVar=NULL, # change to ethnicity variable name
                            RestrictToIndicator=NULL, #write category ethnicity (i.e., 1 or "White) 
                            
                            number_cores=8,
                            runparallel=TRUE,
                            destinationfolder=save_to,
                            savelog=TRUE,
                            
                            cohort=cohort,
                            analysisdate=current_date,
                            analysisname=analysis_name)

