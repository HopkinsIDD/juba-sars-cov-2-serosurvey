#' @title Predict IgG from the standard dose response curve
#' 
#' @param y raw value
#' @param fit dose response curve fit
#' 
#' @return IgG value
antiLL4 <- function(y, fit) {
  # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021
  
  pVec <- as.numeric(fit$parmMat)
  
  b <- pVec[1]
  c <- pVec[2]
  d <- pVec[3]
  e <- pVec[4]
  # inverse of equation 2
  x <- exp(log((d-c)/(y-c)-1)/b+log(e))
  
  return(x)
}


#' @title Run ELISA analysis using raw plate data
#' @description Pulls raw plate data and calculates IgG
#' 
#' @param save_file T/F whether or not to save a csv version of the analyzed results
#' 
#' @return a data table with IgG estimates by sample_id
run_elisa_analysis <- function(save_file) {
  
  # load packages and installs in case they aren't available
  pacman::p_load('drc', 'data.table', 'kableExtra', 'tidyverse', 'openxlsx', 
                 'EnvStats', 'ggplot2', 'ggbeeswarm', 'gridExtra')
  
  # some key paths
  rel_path <- here::here('data/raw_data/elisa_plate_metadata/')
  rel_rawdata_path <- here::here('data/raw_data/elisa_plate_data/')
  
  # load data
  all_files <- list.files(rel_rawdata_path)
  plate_files <- all_files[grep('SS_R', all_files)]
  data <- lapply(paste0(rel_rawdata_path, plate_files), fread)
  names(data) <- plate_files
  
  # load all metadata for plate templates
  all_meta <- list.files(rel_path)
  meta_files <- all_meta[grep('SS_Plate_Lookup_column_format_w_notations_', all_meta)]
  meta <- lapply(paste0(rel_path, meta_files), read.xlsx)
  names(meta) <- meta_files
  
  # pull sample/barcode lookup table
  lookup <- data.table(read.xlsx(paste0(rel_rawdata_path, 'sample_id_lookup.xlsx')))
  
  # fix some typos in lookup table
  lookup[sample_id == '00057', sample_id := '0057']
  lookup[sample_id == '00039', sample_id := '0039']
  
  # standard curve values
  igg_ng <- c(25, 12.5, 6.25, 3.125, 1.5625, 0.78125, 0.390625, 0.1953125)
  igg_ug <- igg_ng/1000
  
  # function to analyze plate data
  analyze_plate <- function(in_dat, plate_name) {
    
    # get plate ID and corresponding template
    plateID <- gsub('.*[_]([^.]+)[_].*', '\\1', plate_name)
    templateID <- unique(lookup[plate_id == plateID, template_id])
    if (length(templateID) > 1) {
      stop('There is more than one template for plate: ', plateID)
    } else {
      message('Analyzing plate: ', plateID)
    }
    
    # pull the template that corresponds to the data
    meta_plate <- meta[[paste0('SS_Plate_Lookup_column_format_w_notations_', templateID, '.xlsx')]]
    mt <- data.table(well = as.character(meta_plate[2, 3:386]),
                     well_id = as.character(meta_plate[3, 3:386]))
    
    # just keep the needed data cells
    dt <- data.table(well = as.character(in_dat[3, 3:386]), 
                     value = as.numeric(in_dat[4, 3:386]))
    
    # merge ids on and sort
    dt <- merge(mt, dt, by = 'well')
    setorderv(dt, 'well_id')
    
    # add plate id, date run, and replicate
    dt[, plate_id := plateID]
    dt[, template_id := templateID]
    dt[, date_run := gsub('.*[2020]([^.]+)[.csv].*', '\\1', plate_name)]
    dt[, replicate := gsub('_R[0-9]', '', well_id)]
    
    # average of blanks
    blank_avg <- mean(dt[grep('B_', well_id), value])
  
    # subtract the blank and get average and sd
    dt[, value_avg := mean(value - blank_avg), by = replicate]
    dt[, value_cv := cv(value - blank_avg), by = replicate]
    
    # pull out blanks and standard and keep important cols only
    keepcols <- c('plate_id', 'date_run', 'replicate', 'value_avg', 'value_cv')
    dt_std <- unique(dt[grep('S[0-9]_', well_id), keepcols, with = FALSE])
    dt_sample <- unique(dt[!grep('B_|S[0-9]_', well_id), keepcols, with = FALSE])
    
    # keep values that have less than 1 residual SD from the linear model fit
    dt_std[, true_igg := igg_ug]
    lm_fit <- lm(true_igg ~ value_avg, data = dt_std)
    outside <- which(resid(lm_fit) > sd(resid(lm_fit)))
    if (!is.na(outside[1])) {
      dt_std[outside, exclude_std := TRUE]
      dt_std[!outside, exclude_std := FALSE]
    } else {
      dt_std[, exclude_std := FALSE]
    }
    
    # flag values outside of the curve
    std_range <- range(dt_std[exclude_std == FALSE, value_avg])
    dt_sample[, above_std := ifelse(value_avg > std_range[[2]], T, F)]
    dt_sample[, below_std := ifelse(value_avg < std_range[[1]], T, F)]
    
    # fit dose response curve
    dr_fit <- drm(value_avg ~ true_igg, data = dt_std[exclude_std == FALSE], fct = LL.4())
    
    # bind sample data with standard so we can save the model fits for both
    dt_sample <- rbindlist(list(dt_sample, dt_std), use.names = TRUE, fill = TRUE)
    
    # predict using the standard curve
    preds <-antiLL4(dt_sample$value_avg, dr_fit)
    dt_sample[, prediction := preds]
    
    # set anything above or below the standard curve to NA
    dt_sample[above_std == TRUE, prediction := NA]
    dt_sample[below_std == TRUE, prediction := NA]
    
    # adjust by dilution factor
    dt_sample[, dilution := as.numeric(gsub('.*(\\d{4}).*', '\\1', replicate))]
    dt_sample[, igg := prediction * dilution]
    
    # calculate average by sample
    dt_sample[, position_id := gsub('\\_.*', '', replicate)]
    dt_sample[, position_number := as.numeric(gsub('\\_.*', '', replicate))]
    dt_sample[, igg_avg := mean(igg, na.rm = T), by = position_id]
    
    # add replicate values back in so we can report this with the CV analysis
    dt_sample <- merge(dt_sample, dt, by = keepcols, all = T)
    
    # end function
    return(dt_sample)
  }
  
  # apply function to get data from plates
  results <- mapply(analyze_plate, data, plate_name = names(data),
                    SIMPLIFY = FALSE)
  
  # bind and merge with actual sample IDs
  results <- rbindlist(results, use.names = T)
  results <- merge(lookup, results, by = c('plate_id', 'position_number', 'template_id'), 
                   all.y = TRUE)
  
  # add flag for most recent plates (assuming ascending plate_id numbers)
  results[, order_run := as.numeric(gsub('R', '', plate_id))]
  ord_tmp <- unique(results[, c('plate_id', 'order_run')])
  setorderv(ord_tmp, 'order_run', order = -1)
  recent_plates <- unique(ord_tmp)[1:8, plate_id]
  results[, most_recent := ifelse(plate_id %in% recent_plates, TRUE, FALSE)]
  
  # assume that if sample is below the standard that it is 0
  results[below_std == TRUE & is.na(igg_avg), igg_avg := 0]
  results[below_std == TRUE, igg := 0]
  
  # summarize plate calculations
  plates <- unique(results[!grep('S[0-9]', position_id), 
                           c('plate_id', 'order_run', 'position_id', 'igg_avg')])
  
  # prepare the plate data to summarize 
  lk_tmp <- copy(lookup)
  plates2 <- copy(plates)
  plates2[, order_run := NULL]
  plates2[, position_number := as.numeric(position_id)]
  lk_tmp <- merge(lk_tmp, plates2, by = c('plate_id', 'position_number'), all = T)
  lk_tmp[, position_number := NULL]
  
  # subset to just the samples
  setorderv(plates, 'order_run', order = -1)
  prelim <- plates[grep('[0-9]', position_id)]
  
  # remove the repeated samples
  prelim <- merge(prelim, lk_tmp[, c('plate_id', 'position_id', 'sample_id')],
                  by = c('plate_id', 'position_id'))
  prelim[, needed_repeat := duplicated(sample_id, fromLast = TRUE)]
  prelim <- prelim[needed_repeat == FALSE]
  prelim <- prelim[!is.na(sample_id), c('sample_id', 'igg_avg')]
  setnames(prelim, 'igg_avg', 'igg')
  
  # saving cleaned/simplified data
  if (save_file) {
    filename <- here::here('data', 'generated_data', 'serology_data', 
                           'ssd_cleaned_elisa_results.csv')
    write.csv(prelim, filename)
  }
  
  # end function
  return(prelim)
}


#' @title Get cleaned serodynamics data from the Boston cohort
#' 
#' @return a data table with cumulative and weekly cases
clean_serodynamics_data <- function() {

  serodata <- read_csv(here::here("data", "raw_data", "boston_serodynamics",
                                  "boston_serology_data_12_11_20.csv"),
                       col_types = cols(
                         sample_id=col_character()
                       )) %>%
    rename(id=sample_id,
           sample=`sample  ID`,
           IgG=cov2_igg
    ) %>%
    
    mutate(
      new_data=ifelse(cohort %in% c(3), 1, 0),
      cohort=ifelse(cohort %in% c(3), 1, cohort),
      cohort=factor(cohort,levels=c(2,1),
                    labels=c("prepandemic","case")),
      
      #week
      #4 categories
      week=factor(cut(dos,c(-1,7,14,28,80),
                      right = TRUE),
                  labels= c("\u22647 days", "8-14 days",
                            "15-28 days",">28 days")),
      #20 categories
      week2=factor(cut(dos,c(-1,seq(7,140,7)),
                       labels=1:20)),
      week2=factor(week2,levels=c(1:20)),
      
      #Severity
      Severity=ifelse(is.na(Severity),5,Severity),
      Severity=ifelse(Severity==6,1,Severity),
      Severity=factor(Severity,
                      labels=c("Not Hospitalized",
                               "Hospitalized, no ICU",
                               "Hospitalized, required ICU",
                               "Died due to COVID-19",
                               "Missing")),
      
    )
  
  # set limit of detection
  serodata <- data.table(serodata)
  serodata <- serodata[, !grep('X', names(serodata)), with = FALSE]
  min_ab <- min(serodata[new_data == 0, IgG])
  serodata[IgG < min_ab, IgG := min_ab]
  
  # select columns we want
  serodata <- serodata[, c('id', 'cohort', 'dos', 'Severity', 
                           'IgG', 'new_data', 'week')]
  setnames(serodata, names(serodata), tolower(names(serodata)))
  serodata <- serodata[cohort == 'prepandemic' | severity != 'Missing']
  serodata <- serodata[, case := ifelse(cohort == 'case', 1, 0)]
  
  return(serodata)
}


#' @title Prepare positive control validation data for Stan model
#'
#' @param run_id id for run/experiment
#' @param validation_data validation data from Boston cohort
#' @param severity_data probability weights by severity
#' @param sev_cat severity categories
#' @param threshold threshold for seropositivity
#' 
#' @return input validation values for Stan model
get_pos_val_data <- function(run_id,
                             validation_data, 
                             severity_data,
                             sev_cat,
                             threshold) {
  
  # if we're not using all the positive control data, 
  # assign severity group proportions
  if (run_id != 'allpos') {
    
    for (i in 1:length(sev_cat)) {
      validation_data[severity == sev_cat[[i]], severity_prob := severity_data[[i]]]
    }
    
    # create positive control set by severity
    pos_total <- round(nrow(validation_data[severity == 'Not Hospitalized'])/severity_data[[1]])
    all_totals <- round(pos_total*severity_data)
    keep_rows <- validation_data[severity_prob == severity_data[[which(sev_cat == 'Not Hospitalized')]], 
                             row_id]
    for (i in 2:4) {
      keep_rows <- c(keep_rows, 
                     sample(validation_data[severity == sev_cat[[i]], row_id], all_totals[[i]]))
    }
    validation_data <- validation_data[row_id %in% keep_rows]
  }
  
  # number of positive controls from validation data
  # Iyer et al, 2020 and new mild/asymptomatic cohort data
  pos_control <- nrow(validation_data)
  
  # true positive rate for cases
  control_tp <- nrow(validation_data[igg > threshold])
  
  # return
  return(c(pos_control, control_tp))
}


#' @title Find population-level proportions of age and sex
#'
#' @param pop_dt data with population counts by age/sex
#' 
#' @return data with population proportions by age/sex
get_pop_props <- function(pop_dt) {
  pop_dt %>%
    mutate(
      age_med = (age_lower + age_upper) / 2
    ) %>%
    filter(age_cat %in% age_cats) %>%
    group_by(age_cat) %>%
    summarize(
      male = sum(men),
      female = sum(women)
    ) %>%
    pivot_longer(
      cols = -age_cat,
      names_to = 'sex',
      values_to = 'pop'
    ) %>%
    mutate(
      sex = ifelse(sex == 'male', 1, 0),
      pct = pop / sum(pop)
    ) %>%
    full_join(
      expand_grid(
        sex = 0:1
      )
    ) 
}


#' @title Run Stan analysis
#' @description Runs the Stan seroprevalence model
#'
#' @param model_script Stan model script
#' @param dat list with data to run the model
#' @param analysis which analysis to run among different seropositivity thresholds
#' @param run_id which run_id to launch (generally, different proportions of mild in positive controls)
#' @param coef_eqn formula in character format expressing the probability of seropositivity on the logit scale
#' @param pos_control number of predicted positives in the control dataset
#' @param neg_control number of predicted negatives in the control dataset
#' @param control_tp number of true positives in the control dataset
#' @param control_tn number of true negatives in the control dataset
#' @param pop_age_cats population and age category matrix
#' @param n_cores number of cores to use for parallel computation
#' @param sex_ref reference category for sex
#' @param age_ref reference category for age
#' @param redo redo fit or load pre-computed posteriors if available
#' 
#' @details The analysis parameter needs to be one of:
#' - main: main analysis with 0.32 threshold
#' 
#' @return a list with parameter posteriors and results
run_analysis_stan <- function(model_script,
                              dat,
                              analysis = "main",
                              run_id,
                              coef_eqn,
                              pos_control,
                              neg_control,
                              control_tp,
                              control_fp,
                              pop_age_cats,
                              n_cores = detectCores() - 2,
                              sex_ref = 0,
                              age_ref = "[20,30)",
                              redo = F,
                              chains,
                              iter,
                              warmup,
                              control, 
                              ...) {
                                 
  ## Prepare data for analysis
  ana_suffix <- case_when(analysis == "main" ~ "",
                          T ~ as.character(NA))
  if (is.na(ana_suffix))
    stop("Unknown analysis parameter, must be one of main, ...")
  
  # Set analysis data
  ana_dat <- as_tibble(dat)
  
  for (var in c("pos", "neg")) {
    ana_dat[[var]] <- as.numeric(dat[[paste0(var, ana_suffix)]])
  }
  
  # Set model matrix
  X <- model.matrix(as.formula(paste("~", coef_eqn)), data = ana_dat)
  
  # name of the stan output file 
  # just in case results directory isn't created already
  if(!dir.exists(here::here("data", "generated_data", "model_fits"))){
    dir.create(here::here("data", "generated_data", "model_fits"))
  }
    
  stan_out_file <- here::here("data", "generated_data", "model_fits", 
                              paste0("stan_fit_", run_id,".rds"))
  
  if (!file.exists(stan_out_file) | redo) {
    
    re_model <- stan_model(model_script)
    
    # Run stan model
    stan_est <- sampling(re_model,
                         data = list(
                           N_survey = nrow(ana_dat),
                           p_vars = ncol(X),
                           X = X,
                           survey_pos = ana_dat$pos,
                           N = 1000, # times we sampled the positive control dataset
                           N_pos_control = pos_control,
                           control_tp = control_tp,
                           N_neg_control = neg_control,
                           control_fp = control_fp
                         ),
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         control = control)
    
    saveRDS(stan_est, stan_out_file)
  } else {
    cat("Model not re-run. Loading pre-computed posteriors from ", stan_out_file, "\n")
    stan_est <- readRDS(stan_out_file)
  }
  
  ## extract log-likelihoods
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  
  ## extract parameters
  beta <- rstan::extract(stan_est, pars = "beta")[[1]]
  
  pop_cat_mat <- pop_age_cats %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)
  
  ## compute estimates by age category
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  pop_cat_p <- foreach(i = 1:nrow(pop_cat_mat),
                       .combine = rbind, 
                       .inorder = F, 
                       .packages = c("tidyverse", "foreach")) %dopar% 
    { 
      foreach(j = 1:nrow(beta), 
              .combine = rbind, 
              .inorder = T) %do% 
        {
          
          # Compute probability integrating across in household random effects
          prob <- integrate(function(x) {
            plogis(qnorm(
              x, beta[j, , drop = F] %*% t(pop_cat_mat[i, , drop = F])
            ))
          }, 0, 1)[[1]]
          
          tibble(
            age_cat = pop_age_cats$age_cat[i],
            sex = pop_age_cats$sex[i],
            pop = pop_age_cats$pop[i],
            seropos = prob
          ) %>%
            mutate(sim = j)
        }
    }
  parallel::stopCluster(cl)
  
  ## overall estimate
  overall_re <- pop_cat_p %>%
    mutate(
      var = "Overall",
      val = ""
    ) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  ## find age specific probabilities in order to make relative risks
  age_re <- pop_cat_p %>%
    filter(sex == sex_ref) %>%
    mutate(var = "Age") %>%
    rename(val = age_cat) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  # sex-specific probabilities
  sex_re <- pop_cat_p %>%
    filter(age_cat == age_ref) %>%
    mutate(var = "sex") %>%
    rename(val = sex) %>%
    group_by(sim, var, val) %>%
    summarize(p = weighted.mean(seropos, pop)) %>%
    ungroup()
  
  # restuls list
  res <- list(
    beta = rstan::extract(stan_est, pars = "beta")[[1]],
    model_mtx = X,
    sens = rstan::extract(stan_est, pars = "sens")[[1]],
    spec = rstan::extract(stan_est, pars = "spec")[[1]],
    obs = nrow(ana_dat),
    pos = sum(ana_dat$pos),
    neg = sum(ana_dat$neg),
    stan_ll = stan_ll,
    pop_cat_p = pop_cat_p,
    subset_est = rbindlist(
      list(
        overall_re,
        sex_re,
        age_re
      )
    )
  )
  
  return(res)
}


#' @title Summarize relative risk of seropositivity by age and sex
#'
#' @param est seroprevalence estimats by age and sex
#' @param age_ref reference category for age
#' 
#' @return relative risk table
get_rr_by_group <- function(est, age_ref = '[20,30)') {
  est %>%
    filter(var == 'Age') %>%
    group_by(sim) %>%
    mutate(rr = ifelse(val == age_ref, NA, p / p[val == age_ref])) %>%
    ungroup() %>%
    left_join(sero_dat %>%
                group_by(age_cat) %>%
                summarize(
                  n = n(),
                  pos = sum(pos),
                  neg = sum(neg)
                ),
              by = c('val' = 'age_cat')
    ) %>%
    bind_rows(
      subset_est %>%
        filter(var == 'sex') %>%
        group_by(sim) %>%
        mutate(rr = ifelse(val == 0, NA, p / p[val == 0])) %>%
        ungroup() %>%
        mutate(val = ifelse(val == 0, 'Female', 'Male')) %>%
        left_join(sero_dat %>%
                    mutate(val = ifelse(sex == 'male', 'Male', 'Female')) %>%
                    group_by(val) %>%
                    summarize(
                      n = n(),
                      pos = sum(pos),
                      neg = sum(neg)
                    ))) %>%
    group_by(var, val, n, pos, neg) %>%
    summarize(
      `Relative risk (95% CI)` = ifelse(is.na(mean(rr)), '--',
                                        paste0(
                                          mean(rr, na.rm = T) %>%
                                            formatC(2, format = 'f'),
                                          ' (', quantile(rr, probs = .025, na.rm = T) %>%
                                            formatC(2, format = 'f'), '-',
                                          quantile(rr, probs = .975, na.rm = T) %>%
                                            formatC(2, format = 'f'), ')'
                                        )
      ),
      p = ifelse(is.na(mean(rr)), '--',
                 min(2 * c(
                   mean(rr > 1, na.rm = T),
                   mean(rr < 1, na.rm = T)
                 )) %>%
                   formatC(3, format = 'f')
      )
    ) %>%
    ungroup() %>%
    mutate(
      pos = paste0(pos, ' (', formatC(100 * pos / n, 1, format = 'f'), '%)'),
      neg = paste0(neg, ' (', formatC(100 * neg / n, 1, format = 'f'), '%)')
    ) %>%
    rename(
      `Test positive` = pos, `Test negative` = neg,
      Obs = n, Category = val
    )
}
