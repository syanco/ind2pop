#######################################
#                                     #
#   Individual to Population Niches   #
#           Empirical Demos           #
#                                     #
#          Scott Yanco, PhD           #
#        scott.yanco@yale.edu         #
#                                     #
#######################################

# DESCRIPTION:  Top-level workflow script to run empirical demonstration for the
#   Individual to Population Niches Project


# Activate conda env
# conda activate niche_mix

# Declare working directory
wd=~/projects/ind2pop/
cd $wd

#-----------------------#
#---- Annotate Data ----#
#-----------------------#

Rscript $wd/src/workflow/annotate_data.r


#---------------------------------#
#---- Calculate Niche Metrics ----#
#---------------------------------#

# Run calculation script
Rscript $wd/src/workflow/calc_niches.r


#--------------------------#
#---- Make Home Ranges ----#
#--------------------------#

Rscript $wd/src/workflow/make_akde.r


#------------------------------#
#---- Make Mean LST Raster ----#
#------------------------------#

# IMPORTANT: Do not run this from command line - run in GEE code editor
$wd/src/workflow/make_lst_mean_rast.js


#-------------------------------#
#---- Calc Vulnlerabilities ----#
#-------------------------------#

# Outputs an RDS of a df with future weights per ind
Rscript $wd/src/workflow/calc_vulnerabilities.r

#--------------------#
#---- Make Plots ----#
#--------------------#

# Switch conda env
# conda activate plots

# Run plotting script script
Rscript $wd/src/workflow/make_figs.r

# Run 