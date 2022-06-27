#######################################
#                                     #
#   Individual to Population Niches   #
#           Empirical Demos           #
#                                     #
#          Scott Yanco, PhD           #
#        scott.yanco@yale.edu         #
#                                     #
#######################################


#-----------------------#
#---- Annotate Data ----#
#-----------------------#

# Activate conda env
conda activate niche_mix

# Declare working directory
wd=~/projects/ind2pop/
cd $wd

# Run annotation script
Rscript $wd/src/workflow/annotate_data.r


#---------------------------------#
#---- Calculate Niche Metrics ----#
#---------------------------------#

# Run calculation script
Rscript $wd/src/workflow/calc_niches.r


#--------------------#
#---- Make Plots ----#
#--------------------#

# Switch conda env
conda activate plots

# Run plotting script script
Rscript $wd/src/workflow/calc_niches.r

# Run 