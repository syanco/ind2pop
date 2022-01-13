# Libraries
library(tidyverse)
library(rstoat)
library(lubridate)

#-- Filter out sugar maples and the write to file to avoid loading full FIA data
# trees <- read.csv("data/data/FIA_090720.csv")
# acsa <- trees %>% 
#   filter(genus == "Acer" & species == "saccharum")
# write_csv(acsa, file = "data/sugar_maple.csv")
acsa <- read.csv("data/sugar_maple.csv")


set.seed(4444)
# filter to past 2000, and indviudals with at least 2 records, and take a random sample
samp <- acsa %>% 
  filter(INVYR >= 2000) %>% 
  group_by(IDtree) %>% 
  summarise(n=n()) %>% 
  filter(n >1) %>% 
  sample_n(50) %>% 
  left_join(acsa)

####---- Generate Time series for each ind ----####

# Need to simulate a time series of occurrences between 1st and last record for
# each individual.  Will generate weekly records for each ind.

makeTS <- function(df){
  minyr <- min(df$INVYR)
  maxyr <- max(df$INVYR)
  nyr <- maxyr-minyr+1
  # TODO: there's probably a more elegant way to make the vector of years...
  yrvec <- list()
  for(i in 1:length(seq(minyr, maxyr, 1))){
    yrvec[[i]] <- rep(seq(minyr, maxyr, 1)[i], 52)
  }
  yrvec <- unlist(yrvec)
  
  # grab ist row
  r1 <- df[1,]
  
  #duplicate it nyrs*52 tims
  new <- uncount(r1, nyr*52) %>% 
    # construct dates for annotations
    mutate(week = rep(1:52, nyr),
           year = yrvec,
           timestamp = as.Date(paste(year, week, 1, sep="-"), "%Y-%U-%u"))
}

# create new timeseries for each individual
tree_ts <- samp %>% 
   split(.$IDtree) %>% 
   map(~makeTS(.x)) %>% 
    bind_rows()


write_csv(tree_ts, "data/sugar_maple_samp.csv")

  
