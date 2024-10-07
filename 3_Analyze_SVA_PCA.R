#==============================================================================#
#                                                                              #
#                    Comparative Analysis of SVA and PCA                       #
#                to deal with variance in high throughput data                 #
#                                                                              #
#                       Part 3: Analyze SVA or PCA                             #
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
# - Methylation_path <- "sva_pca/db/bVals_BISC_nosib_noout_nomask.RDS" #bisc
# - Covariables data in DATAFRAME format (.rds), where samples x covariables
# - PCs inside Large prcomp object 
# - SVs as LAreg lsit object, returned by sva() 

# Covariables/Metadata dataframe needs a colum named "Basename" with the
# samples IDs, the samples IDs have to match the columnames of Methylation data


# [OUPUT]:
# - csv with explained and cumulative variance for PCs and SVs
# . plots with cumulative variance for all PCs/SVs and for 30 first PCs/SVs
# - Interactive, different number of optimal PCs and SVs
# - csv with association of PCs and SVs with covariables
# - csv with associations between PCs and SVs

# [CONSIDERATIONS]


#==============================================================================#
#                           USER PARAMETERS                                    #----
#==============================================================================#

# Paths to methylation and metadata
# pca_res_path <- 'sva_pca/results/residuals_variance/pca_residuals_21092024.rds'    # BiSC
# sva_res_path <- 'sva_pca/results/residuals_variance/sva_residuals_21092024.rds'
# metadata_path <- "sva_pca/db/sva_pca_metadata_20092024.rds"
# methylation_path <- "sva_pca/results/imputed_methylation_21092024.rds"

pca_res_path <- 'sva_pca/results/pca_woStromal_residuals_woStromal_INMA_03102024.rds'    # INMA
sva_res_path <- 'sva_pca/results/sva_woStromal_residuals_woStromalINMA_03102024.rds'
metadata_path <- "sva_pca/db/INMA_sva_pca_metadata_21092024.rds"
methylation_path <- "sva_pca/results/imputed_methylation_INMA_03102024.rds"

current_date <- '21092024'

# covars <- c('maternal_smoking', "Sex", "maternal_educational_level", "gestational_age",    # BiSC
#             "study_center", "covid_confinement", "child_ethnicity",
#             'Trophoblasts','Stromal', 'Hofbauer', 'Endothelial','nRBC', 
#             'Syncytiotrophoblast')

# to do:  in bisc we have to add technical varibles to llok at assocaiiton with PCs and SVs 


covars <- c("Smoke", "estudios3c", "ethnic_origin2c", "Cohort", "sges", "Sex",  #INMA
            'Trophoblasts', 'Hofbauer', 'Endothelial','nRBC',
                    'Syncytiotrophoblast')

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
#library(sfsmisc)
library(stats)
library(readr)
#devtools::install_github("haotian-zhuang/findPC")
library(findPC)
#devtools::install_github("heatherjzhou/PCAForQTL")
library(ggplot2)
#install.packages("pathviewr")
library(pathviewr)
library(reshape2)

#------------------------------------------------------------------------------#
#                           Define functions                                   #----
#------------------------------------------------------------------------------#

#https://github.com/heatherjzhou/PCAForQTL/blob/master/R/22.01.04_main1.1_runElbow.R
runElbow <- function(X=NULL,
                     prcompResult=NULL){
  
  #Obtain prcompResult.
  if(is.null(prcompResult)){
    if(is.null(X)){
      stop("Please input X or prcompResult.")
    }else{
      cat("Running PCA...\n")
      prcompResult<-prcomp(X,center=TRUE,scale.=TRUE) #This may take a moment.
    }
  }else{
    if(class(prcompResult)!="prcomp"){
      stop("prcompResult must be a prcomp object returned by the function prcomp().")
    }
  }
  
  importanceTable<-summary(prcompResult)$importance
  x<-1:ncol(importanceTable) #PC indices.
  y<-importanceTable[2,] #PVEs.
  
  #Given x  and y, calculate the distance between each point and the diagonal line (the line connecting the first and last points).
  x1<-x[1] #First point.
  y1<-y[1] #First point.
  
  x2<-x[length(x)] #Last point.
  y2<-y[length(y)] #Last point.
  
  x0<-x
  y0<-y
  
  distancesDenominator<-sqrt((x2-x1)^2+(y2-y1)^2)
  distancesNumerator<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))
  distances<-distancesNumerator/distancesDenominator
  # plot(distances)
  
  numOfPCsChosen<-which.max(distances) #12.
  names(numOfPCsChosen)<-NULL
  return(numOfPCsChosen)
}



# Function to calculate R² or adjusted R² for continuous variables
calculate_r2 <- function(model) {
  summary(model)$adj.r.squared
}

# Generalized function to test associations for both PCs and SVs
test_association <- function(components, covariates, data, component_type = "PC") {
  results <- data.frame(Component = character(), Covariate = character(), P.Value = numeric(), Radj2 = numeric())
  
  for (component in components) {
    for (covar in covariates) {
      if (is.factor(data[[covar]]) || is.character(data[[covar]])) {
        # If the covariate is categorical, use ANOVA
        model <- aov(as.formula(paste(component, "~", covar)), data = data)
        p_val <- summary(model)[[1]][["Pr(>F)"]][1]
        radj2 <- calculate_r2(lm(as.formula(paste(component, "~", covar)), data = data))  # Adjusted R² for categorical
      } else {
        # If the covariate is continuous, use linear regression
        model <- lm(as.formula(paste(component, "~", covar)), data = data)
        p_val <- summary(model)$coefficients[2, 4]  # p-value for covariate
        radj2 <- calculate_r2(model)  # Adjusted R² for continuous
      }
      
      # Store the results
      results <- rbind(results, data.frame(Component = paste0(component_type, " ", component), Covariate = covar, P.Value = p_val, Radj2 = radj2))
    }
  }
  
  return(results)
}



# Function to calculate the variance explained by SVs based on posterior probabilities from SVA
# Code adpated from D.Pelegrí, who based it on: https://support.bioconductor.org/p/88553/

variance_explained_sv <- function(dat, pprob.gam, pprob.b, n.sv) {
  pprob <- pprob.gam * (1 - pprob.b)
  dats <- dat * pprob
  dats <- dats - rowMeans(dats)  # Center the data
  
  # Eigen decomposition of the scaled data
  uu <- eigen(t(dats) %*% dats)
  
  # Percent variance explained (PVE) by each SV/PC
  uu_val <- uu$values / sum(uu$values)  # Proportion of each eigenvalue over total sum
  
  sv <- data.frame(
    SV = paste0("SV", 1:length(uu_val)),           
    PVE = uu_val)
  return(sv)
}



test_association_pcs_svs <- function(pcs, svs, data) {
  results <- data.frame(PC = character(), SV = character(), P.Value = numeric(), Radj2 = numeric())
  
  for (pc in pcs) {
    for (sv in svs) {
      # Perform linear regression to test association between each PC and each SV
      model <- lm(as.formula(paste(sv, "~", pc)), data = data)
      p_val <- summary(model)$coefficients[2, 4]  # p-value for the PC-SV association
      radj2 <- calculate_r2(model)  # Adjusted R² value
      
      # Store the results
      results <- rbind(results, data.frame(PC = pc, SV = sv, P.Value = p_val, Radj2 = radj2))
    }
  }
  
  return(results)
}


#------------------------------------------------------------------------------#
#                            LOAD METADATA                                     #----
#------------------------------------------------------------------------------#

metadata <- readRDS(metadata_path)




#==============================================================================#
#                                                                              #
#                    PRINCIPAL COMPONENT ANALYSIS (PCA)                        #----
#                                                                              #
#==============================================================================#

cat('-------------------------- PCA -------------------------------------------')

pca_residuals <- readRDS(pca_res_path)


#------------------------------------------------------------------------------#
#                     Get explained and cumulative variance PCs                #----
#------------------------------------------------------------------------------#

# Get variance explaine for each PCs and cumualtive varicance 
explained_variance <- summary(pca_residuals)$importance[2, ]  # Proportion of variance explained by each PC
cumulative_variance <- cumsum(explained_variance)


# Plot cumulative var  
jpeg(paste0("sva_pca/results/explained_variance_resPCs_",cohort,"_", current_date, ".jpeg"))
plot(explained_variance , type = "b", xlab = "Principal Component", 
     ylab = "Variance Explained")
#abline(h = 0.8, col = "red", lty = 2)  # Line at 80% variance explained
dev.off()


jpeg(paste0("sva_pca/results/explained_variance_resPCs_zoom30", cohort,"_", current_date, ".jpeg"))
plot(explained_variance[0:30] , type = "b", xlab = "Principal Component", 
     ylab = "Variance Explained")
#abline(h = 0.8, col = "red", lty = 2)  # Line at 80% variance explained
dev.off()


pc <- data.frame(explained_variance, cumulative_variance)
write.table(pc, 
paste0('sva_pca/results/explained_cumualtive_variance_resPCs_',cohort,"_", current_date,'.csv'))


# Alternative (with more decimals)
#explained_variance2 <- (pca_residuals$sdev^2) / sum(pca_residuals$sdev^2)  


#------------------------------------------------------------------------------#
#                      Select optimal number of PCs                            #----
#------------------------------------------------------------------------------#


# OPTION 1: FindPC() -----------------------------------------------------------

# Methods include 'all', 'piecewise linear model', 'first derivative', 
#'second derivative', 'preceding residual', 'perpendicular line (default)', 
#''k-means clustering'

jpeg(paste0("sva_pca/results/find_optPC_plot_", cohort,"_",current_date,".jpeg"), 
     width = 800, height = 600)

opt_PCs <- findPC(sdev = pca_residuals$sdev,
                  number=30, #it's because of the number 
                  method='all', 
                  aggregate=NULL,
                  figure=TRUE)

dev.off()

# Optimal PCs are: 6

# OPTION 2: runElbow() perpendicular approach ----------------------------------

optimal_pc_runelbow <- runElbow(prcompResult = pca_residuals)
cat("Optimal number of PCs using runElbow():", optimal_pc_runelbow, "\n") #22

# OPTION 3: find_curve_elbow perpendicular approach -----------------------------

df_explained_variance <- data.frame(explained_variance)
df_explained_variance$PC <- as.numeric(gsub('PC', "", rownames(df_explained_variance)))

optimal_pc_curve <- find_curve_elbow(data_frame = df_explained_variance)
cat("Optimal number of PCs using find_curve_elbow():", optimal_pc_curve, "\n") #22


#------------------------------------------------------------------------------#
#                     Get metadata with PCs for EWAS analysis                  #----
#------------------------------------------------------------------------------#

# Add optimal number of PCs to metadata --> MANAL? 
n_pc <- 16

# Remove sample not in pca_residuals 
metadata <- metadata[match(rownames(pca_residuals$x), metadata$Basename),]
metadata_pc <- cbind(metadata, pca_residuals$x[,1:n_pc])
write_rds(metadata_pc, paste0('sva_pca/results/metadata_',cohort,'_',n_pc,'pcs_',current_date,'.rds'))
cat('Optimal PCs added to metadata\n')
rm(metadata_pc)


#------------------------------------------------------------------------------#
#               Associations PCs with covariables                              #----
#------------------------------------------------------------------------------#

pcs <- as.data.frame(pca_residuals$x[, 1:10])  # First 10 PCs
colnames(pcs) <- paste0("PC", 1:10)

# Combine PCs with covariates
metadata <- metadata[metadata$Basename %in% rownames(pcs),]
data_for_analysis <- cbind(metadata, pcs)

# Test associations between PCs and covariates
pc_association_results <- test_association(components = colnames(pcs), 
                                           covariates = covars, 
                                           data = data_for_analysis, 
                                           component_type = "PC")


write.table(pc_association_results, 
            paste0('sva_pca/results/assoc_10pca_vars', cohort, '_', current_date, '.csv'), 
            sep=";")



#==============================================================================#
#                                                                              #
#                      SINGULAR VALUE DECOMPOSITION                            #----
#                                                                              #
#==============================================================================#

cat('--------------------------   SVA  -----------------------------------------')

sva_residuals <- readRDS(sva_res_path)

# List of 4
# $ sv       : num [1:553, 1:114] -0.0781 0.0366 -0.018 0.0158 -0.0218 ...
# $ pprob.gam: num [1:807201] 1 1 1 0.996 1 ...
# $ pprob.b  : num [1:807201] 0 0 0 0 0 0 0 0 0 0 ...
# $ n.sv     : num 114


#------------------------------------------------------------------------------#
#                     Get explained and cumulative variance SVs                #----
#------------------------------------------------------------------------------#

# NOTE: Unlike PCs we don't have direclty the variance from SVs, we need to 
# compute it. !!! We need imputation meth, cause sva can't handle NAs. 

meth_data <- readRDS(methylation_path)  # Cpgs x samples


explained_variance <- variance_explained_sv(
  dat = meth_data,                       
  pprob.gam = sva_residuals$pprob.gam,   
  pprob.b = sva_residuals$pprob.b,       
  n.sv = sva_residuals$n.sv              
)

# Cumulative variance 
cumulative_variance <- cumsum(explained_variance$PVE)

sv <- data.frame(explained_variance,cumulative_variance)

# Save explained an cumulative var
write.table(sv, 
paste0('sva_pca/results/explained_cumualtive_variance_resSVs_', cohort, '_', current_date,'.csv'))


jpeg(paste0("sva_pca/results/explained_variance_resSVs_", cohort, current_date, ".jpeg"))
plot(explained_variance$PVE, type = "b", xlab = "Surrogate Variable", 
     ylab = "Variance Explained (estimate)")
dev.off()



jpeg(paste0("sva_pca/results/explained_variance_resSVs_zoom30_", cohort, current_date, ".jpeg"))
plot(explained_variance$PVE[1:30], type = "b", xlab = "Surrogate Variable", 
     ylab = "Variance Explained (estimate)")
dev.off()



#------------------------------------------------------------------------------#
#                          Select optimal number of SVs                        #----
#------------------------------------------------------------------------------#

# OPTION 1: use n.sv() parameter that sva() returns. In BiSC it was n.sv = 114

# OPTION 2: findPC() -----------------------------------------------------------

# NOTE: we need the sdev, in PCA sdev are directly linked to the eigenvalues of 
# the covariance matrix of the data.  findPC() function is designed to work with 
# standard deviations (sdev) because it expects data in the form of variances 
# (which are derived from the eigenvalues). 
# The PVE gives us normalized variance information, but to use it with findPC(), 
# we need to "reverse-engineer" it into a form the function can handle—hence,
# the need to compute the sdev from PVE.


# 1) Get eigenvalyes 
eigenvalues <- explained_variance$PVE * sum(explained_variance$PVE)

# 2) Compute the "sdev" values (analogous to PCA standard deviations)
sdev_sv <- sqrt(eigenvalues)


jpeg(paste0("sva_pca/results/find_optSC_plot_", cohort, '_', current_date,".jpeg"), 
     width = 800, height = 600)
opt_PCs <- findPC(sdev = sdev_sv,
                  number=30, #it's because of the number 
                  method='all', 
                  aggregate=NULL,
                  figure=TRUE)

dev.off()


# OPTION 4: find_curve_elbow() perpendicular approach -----------------------------





df_explained_variance <- data.frame(explained_variance$PVE)
df_explained_variance$SV <- as.numeric(gsub('SV', "", rownames(df_explained_variance)))

optimal_sv_curve <- find_curve_elbow(data_frame = df_explained_variance)
cat("Optimal number of SVs using find_curve_elbow():", optimal_sv_curve, "\n") #21


#------------------------------------------------------------------------------#
#                     Get metadata with PCs for EWAS analysis                  #----
#------------------------------------------------------------------------------#

# Add optimal number of SVs to metadata 
n_sv <- 14
svs <- data.frame(sva_residuals$sv[,1:n_sv])  
colnames(svs) <- paste0("SV", 1:n_sv)
metadata_sv <- cbind(metadata, svs)
write_rds(metadata_sv, paste0('sva_pca/results/metadata_', cohort, "_",n_sv, 'svs_',current_date,'.rds'))

#------------------------------------------------------------------------------#
#                     Associations  SCs with covariables                       #----
#------------------------------------------------------------------------------#

# Step 1: Extract the first 10 SVs from the sva_residuals
sv_matrix <- sva_residuals$sv  # Assuming this contains all SVs
svs <- as.data.frame(sv_matrix[, 1:10])  # First 10 SVs
colnames(svs) <- paste0("SV", 1:10)
data_for_analysis <- cbind(metadata, svs)

# Test associations between PCs and covariates
sv_association_results <- test_association(components = colnames(svs), 
                                           covariates = covars, 
                                           data = data_for_analysis, 
                                           component_type = "SV")

write.table(sv_association_results, 
            paste0('sva_pca/results/assoc_10sva_vars', cohort, '_', current_date, '.csv'), 
            sep=";")





#==============================================================================#
#                                                                              #
#                         ASSOCITATION SVA - PCA                               #----
#                                                                              #
#==============================================================================#

pcs <- as.data.frame(pca_residuals$x[, 1:10])  # Assuming pca_residuals$x contains the PC scores
colnames(pcs) <- paste0("PC", 1:10)

svs <- as.data.frame(sva_residuals$sv[, 1:10])  # Assuming sva_residuals$sv contains the SVs
colnames(svs) <- paste0("SV", 1:10)

# Combine PCs and SVs into one dataset
data_for_analysis <- cbind(pcs, svs)


association_results_pcs_svs <- test_association_pcs_svs(pcs = colnames(pcs), 
                                                        svs = colnames(svs), 
                                                        data = data_for_analysis)

write.table(association_results_pcs_svs, 
            paste0('sva_pca/results/assoc_pcs_svs_',cohort, '_', current_date, '.csv'),
            sep = ";", row.names = FALSE)


# Generate the heatmap for Adjusted R² values
association_results_long <- melt(association_results_pcs_svs, 
                                 id.vars = c("PC", "SV"), 
                                 measure.vars = "Radj2", 
                                 variable.name = "Metric", 
                                 value.name = "Radj2")

ggplot(association_results_long, aes(x = PC, y = SV, fill = Radj2)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "gray", name = "Adj. R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Adjusted R² Values between PCs and SVs",
       x = "Principal Components (PCs)",
       y = "Surrogate Variables (SVs)") +
  coord_fixed()


#  Generate the heatmap for P-Values
association_results_long_pvalue <- melt(association_results_pcs_svs, 
                                        id.vars = c("PC", "SV"), 
                                        measure.vars = "P.Value", 
                                        variable.name = "Metric", 
                                        value.name = "P.Value")

ggplot(association_results_long_pvalue, aes(x = PC, y = SV, fill = P.Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "red", high = "blue", na.value = "gray", name = "P-Value", trans = "log10") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of P-Values between PCs and SVs",
       x = "Principal Components (PCs)",
       y = "Surrogate Variables (SVs)") +
  coord_fixed()









