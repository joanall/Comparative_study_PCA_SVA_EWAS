#! /bin/bash
#SBATCH --job-name=basicIN
#SBATCH --partition=long
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-user=joana.llaurado@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=80gb
#SBATCH --output=/PROJECTES/BISC_OMICS/analyses/BiSC_23/EWAS_diet_CL/sva_pca/basicIN_%A_%a.txt

module purge > /dev/null 2>&1


source ~/.bashrc
source activate PACEanalysis #name of your environment 


Rscript sva_pca/script/4_Run_EWAS.R

# To run it : sbatch --array=1-6 your_script.sh