
#create simple not in function (negation of `%in%`)
`%notin%` <- Negate(`%in%`)

# function to calculate individual-specific time-varying niche size
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNTS <- function(ind_ID, data, interval, min_obs = 2, vars, log=F){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    group_by_(interval) %>% #group by the supplied time interval
    mutate(n=n()) %>% #calculate # of observations per interval
    filter(n >= min_obs) %>% #filter out those below the supplied minimum
    ungroup() #ungroup for next operation
  
  #calculate niches
  indMVNH <- dat_ind %>%
    group_by_(interval) %>% #group by time interval
    group_split() %>% #split into lists of DFs by interval
    lapply(., FUN = function(x){ #apply across elements of ths list
      x %>% select(all_of(vars)) %>% #select only the supplied vars
        MVNH_det(data = ., log = log) #calculate niche hypevolume
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  #create matching identifier information to link back up to the niche estimates
  indID <- dat_ind %>% 
    group_by_(interval) %>% #group by time interval
    group_split() %>% # split into lists of DFs by interval
    lapply(., FUN = function(x){#apply across elements of the list
      #summarize the identifier info (just use first record)
      x %>% summarize(WKxYR = WKxYR[1], # week X yr factor
                      ind = individual_id[1], #ind ID
                      ts = timestamp[1]) #interval starting time stamp
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  out <- cbind(indMVNH, indID) #bind the niche estimates back to the identifiers
  return(out) #return the combined DF for the target individual
}


# function to calculate individual-specific time-varying niche dissimilarity
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNDissim <- function(ind_ID, data, interval, min_obs = 2, vars){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  # .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    group_by_(interval) %>% #group by the supplied time interval
    mutate(n=n()) %>% #calculate # of observations per interval
    filter(n >= min_obs) %>% #filter out those below the supplied minimum
    ungroup() #ungroup for next operation
  
  #get individual total niche to use in calculating dissim below.
  niche0 <- dat_ind %>% 
    select(all_of(vars)) %>% 
    na.omit()
  
  #calculate niches
  indMVNH <- dat_ind %>%
    group_by_(interval) %>% #group by time interval
    group_split() %>% #split into lists of DFs by interval
    lapply(., FUN = function(x){ #apply across elements of this list
      tryCatch(x %>% select(all_of(vars)) %>% #select only the supplied vars
                 na.omit() %>% # omit NAs
                 MVNH_dissimilarity(db1 = ., db2 = niche0), #calculate niche dissim
               error=function(e) {list(Bhattacharyya_distance = NA,
                                       Mahalanobis_distance = NA,
                                       Determinant_ratio = NA)}
      ) # tryCatch
    } #f(x)
    ) %>% # lapply
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  #create matching identifier information to link back up to the niche estimates
  indID <- dat_ind %>% 
    group_by_(interval) %>% #group by time interval
    group_split() %>% # split into lists of DFs by interval
    lapply(., FUN = function(x){#apply across elements of the list
      #summarize the identifier info (just use first record)
      x %>%  
        summarize(WKxYR = WKxYR[1], # week X yr factor
                  ind = individual_id[1], #ind ID
                  ts = timestamp[1] #interval starting time stamp
        ) 
    }) %>%
    # bind_rows()
    do.call("rbind", .) %>% #bind the list elements back together
    as.data.frame() #convert to DF (not tibble)
  
  tryCatch(
    out <- cbind(indMVNH, indID) #bind the niche estimates back to the identifiers
  ) # tryCatch
  
  return(out) #return the combined DF for the target individual
}


# function to calculate individual-specific time-varying niche size
# intended to applied across a list of individuals
#
# ind_ID = unique identifier for target individual
# data = DF containing niche and timestamp info
# interval = time slice for calculating niche - 
#             should refer to a factor variable in the data
# min obs = minimum # of observations per time slice
# vars = vector of niche variables to use (declared as objects, not strings)
# log = should `MVNH_det` log transform output?
indNicheAccum <- function(ind_ID, data, interval, min_obs = 2, vars){
  
  if("individual_id" %notin% colnames(data))
  {stop("Data must contain a column called 'individual_id'.")}
  
  #convert vars vector to enquos so that select can read it
  .vars <- enquos(vars)
  
  #prep data for niche size estimation
  dat_ind <- data %>% 
    filter(individual_id==ind_ID) %>% #select target individual
    arrange({{interval}}) %>% 
    group_by({{interval}}) %>% #group by the supplied time interval
    # mutate(n=n()) %>% #calculate # of observations per interval
    # filter(n >= min_obs) %>% #filter out those below the supplied minimum
    summarise(n = n(),
              max = max({{vars}}, na.rm = T),
              min = min({{vars}}, na.rm = T),
              range = max-min) %>%
    ungroup()
  # 
  # nsd_dat <- dat_ind %>% 
  #   make_track(.x = lon, .y = lat, .t = ts, order_by_ts = F, 
  #            crs = CRS("+init=epsg:4326"), all_cols = T) %>% 
  #   mutate(netSQ=nsd(.),
  #          nsd_accum = netSQ+dplyr::lag(netSQ)) #add NSD variable to each individual DF
  # 
  
  return(dat_ind) #return the combined DF for the target individual
}


#-- Individual to Population Functions  --##

# Anticipates vector of sample means, vector of sample variances
estPopVar <- function(var, means, pop_mean, w = NULL){
  n <- length(var)
  ifelse(is.null(w), w <- 1/n, w <- w)
  
  vec <- c()
  
  for(i in 1:n){
    vec[i] <- w*(var[i] + (means[i]^2 - pop_mean^2))
  }
  
  out <- sum(vec)
  return(out)
}

# calculate an individual's total contribution to pop variance 
indContrib <- function(mu_i, var_i, mu_pop, n){
  # n <- length(mu_i)
  cont <- (1/n)*(var_i + mu_i^2 - mu_pop^2)
  return(cont)
}

# calculate the mean component of ind contribution to pop variance
#(variance contrib is just the ind variance itself)
muContrib <- function(mu_i, mu_pop, n){
  cont <- (1/n)*(mu_i^2 - mu_pop^2)
  return(cont)
}

individual_contribution <- function(x, w = NULL) {
  # x is the individual parameter data, columns are mu and sigma (note that it is the sd not the variance),
  # but the returned population sigma2 is the variance
  n = nrow(x)
  if(!is.null(w)) n = w
  mu_i = x[,"mu"]
  sigma_i = x[,"sigma"]
  mu = mean(mu_i)
  marginality_sigma2 = (mu_i^2 - mu^2)/n
  specialization_sigma2 = (sigma_i^2)/n
  sigma2 = sum(marginality_sigma2) + sum(specialization_sigma2)
  marginality_skew = (mu_i^3 - mu^3)/(n*sigma2^(3/2))
  specialization_skew = (3*mu_i*(sigma_i^2)-3*mu*sigma2)/(n*sigma2^(3/2))
  skew = sum(marginality_skew) + sum(specialization_skew)
  
  return(cbind(mu_pop = mu,marginality_sigma2 = marginality_sigma2,
               specialization_sigma2 = specialization_sigma2,
               sigma2_pop = sigma2,
               marginality_skew = marginality_skew,
               specialization_skew = specialization_skew,
               skew_pop = skew))
}


# Creates n bootstrapped samples from individual niche parameters for use in 
# calculating confidence intervals. df is the dataframe of individual niche 
# parameters and n is the number of bootstrapped samples to take. Returns a list
# of length n each element of which is a re-sampled df of length nrow(df). A
# user could call this manually, but it's intended as a sub-function called by 
# the CI function
niche_boot <- function(df, n = 1000){
  new <- list()
  for(i in 1:n){
    nr <- nrow(df)
    idx <- sample(1:nr, size = nr, replace = T)
    new[[i]] <- df[idx,]
  }  
  return(new)
}

# Generates bootstrapped CIs for population mean as estimated from individuals.  
# Calls niche_boot to generate bootsrapped samples.
mean_CIs <- function(df, n = 1000, ci_l = 0.025, ci_h = 0.975){
  new <- niche_boot(df=df, n=n)
  mean_vec <- c()
  for(i in 1:length(new)){
    mean_vec[i] <- mean(na.omit(new[[i]]$mu)) 
  }
  CIs <- quantile(mean_vec, probs = c(ci_l, ci_h)) 
  return(CIs)
}

# Generates bootstrapped CIs for population niche breadth as estimated from individuals.  
# Calls niche_boot to generate bootstrapped samples.
var_CIs <- function(df, n = 1000, ci_l = 0.025, ci_h = 0.975){
  
  # run bootstrap re-samples
  new <- niche_boot(df=df, n=n)
  
  #create empty vector for results
  var_vec <- c()
  
  # loop through bootstrap draws
  for(i in 1:length(new)){
    # get grand mean
    grand_mean <- mean(na.omit(new[[i]]$mu))
    
    #estimate pop var using ind2pop meethod (mixture dist)
    var_vec[i] <- estPopVar(var= na.omit(new[[i]]$var), 
                            means = na.omit(new[[i]]$mu),
                            pop_mean = grand_mean) 
  }
  CIs <- quantile(var_vec, probs = c(ci_l, ci_h)) 
  return(CIs)
}

# based on https://stackoverflow.com/questions/67818541/extract-confidence-interval-for-a-combination-sum-of-variance-estimates-from-a
# gets bootstrapped CIs for pooled variance from a variance components lmm
get_total_var_CIs <- function(mod) {
  varCorr_df <- as.data.frame(VarCorr(mod))
  components <- varCorr_df[['vcov']]
  names(components) <- varCorr_df[['grp']]
  
  var_total <- sum(components)
  
  c(components, "Total" = var_total)
}

# estimates CIs variance using chi-sq distribution.
# can be used to get pooled variance CI (wheer n = n_samples + n_groups [we think...])
var_CI_chi <- function(var, n, alpha){
  a <- 1-(alpha/2)
  b <- alpha/2
  lower_lim <- qchisq(b, n)
  upper_lim <- qchisq(a, n)
  return(c((n-1)*var/upper_lim,(n-1)*var/lower_lim))
  
}


individual_contribution <- function(x, ID = "individual.local.identifier", w = NULL) {
  # x is the individual parameter data, columns are mu and sigma (note that it 
  # is the sd not the variance),
  # but the returned population sigma2 is the variance
  n = nrow(x)
  if(!is.null(w)) n = w
  mu_i = x[,"mu"]
  sigma_i = x[,"sigma"]
  mu = mean(mu_i$mu)
  marginality_sigma2 = (mu_i$mu^2 - mu^2)/n
  specialization_sigma2 = (sigma_i$sigma^2)/n
  sigma2 = sum(marginality_sigma2) + sum(specialization_sigma2)
  marginality_skew = (mu_i$mu^3 - mu^3)/(n*sigma2^(3/2))
  specialization_skew = (3*mu_i$mu*(sigma_i$sigma^2)-3*mu*sigma2)/(n*sigma2^(3/2))
  skew = sum(marginality_skew) + sum(specialization_skew)
  
  return(data.frame(mu_pop = mu,
               marginality_sigma2 = marginality_sigma2,
               specialization_sigma2 = specialization_sigma2,
               sigma2_pop = sigma2,
               marginality_skew = marginality_skew,
               specialization_skew = specialization_skew,
               skew_pop = skew,
               ID = x %>% pull(ID)))
}

