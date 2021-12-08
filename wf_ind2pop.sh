#######################################
#                                     #
#   Individual to Population Niches   #
#           Empirical Demos           #
#                                     #
#          Scott Yanco, PhD           #
#        scott.yanco@yale.edu         #
#                                     #
#######################################

#-------------------------------------#
#---- Build and Populate Mosey DB ----#
#-------------------------------------#

#TODO: write conda yml to file
conda activate mosey

# Database code derived from github.com/bencarlson/mosey_db

# Presumes project directory already exists
wd=~/projects/ind2pop/

MOSEYDB_SRC=~/projects/mosey_db

cd $wd

# Don't run this if database already exists!
# Below is commented out to prevent accidental execution

# cat $MOSEYDB_SRC/db/create_db.sql | sqlite3 data/mosey.db

csvdir=/projects/ind2pop/data/csvs

$MOSEYDB_SRC/db/load_studies.sh --csvdir $csvdir --process civ
