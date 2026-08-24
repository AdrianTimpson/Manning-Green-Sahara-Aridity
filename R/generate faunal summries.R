#---------------------------------------------------------------------
# script to generate faunal summary .csv files direct from database
# will not work without access to database!! Not for public
#---------------------------------------------------------------------
source('functions.R')
faunal <- get.animals()
make.animal.summary(faunal)