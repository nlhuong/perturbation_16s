# (...) - a list options of options to fdapace()
fda_pca <- function(df.long, time_column, value_column, replicate_column,
                    feat_column = NULL, parallel = FALSE, ncores = 1, 
                    num_nonzero = 10, thresh = 0, clust_min_num_replicate = 15,
                    cluster = FALSE, cmethod = "EMCluster", K = 2, 
                    fpca_optns = NULL, fclust_optns = NULL){
  require(dplyr)
  if (!time_column %in% colnames(df.long)) {
    stop("'df.long' does not have", time_column, "as time column.")
  }
  if (!value_column %in% colnames(df.long)) {
    stop("'df.long' does not have", value_column, "as value column.")
  }
  if (!replicate_column %in% colnames(df.long)) {
    stop("'df.long' does not have", replicate_column, "as replicate id column.")
  }
  if (all(!is.null(feat_column), !feat_column %in% colnames(df.long))) {
    stop("'df.long' does not have", feat_column, "as feature id column.")
  }
  if(!is.null(feat_column)){
    features <- unique(df.long[[feat_column]])
  } else {
    feat_column <- "Feature_ID"
    features <- df.long[[feat_column]] <- "Feat1"
  }
  if(parallel) {
    require(doParallel)
    ncores <- min(ncores, detectCores())
    registerDoParallel(ncores)
    res <- foreach(feat=features) %dopar% {
      fit_fdapca(df.long, time_column, value_column, 
                 replicate_column, feat_column, feat, 
                 cluster = cluster, cmethod = cmethod, K = K, 
                 num_nonzero = num_nonzero, thresh = thresh,  
                 clust_min_num_replicate = clust_min_num_replicate, 
                 fpca_optns = fpca_optns, fclust_optns = fclust_optns)
    }
  } else {
    res <- lapply(features, function(feat) {
      cat("Processing feature: ", feat, "...\n")
      fit_fdapca(df.long, time_column, value_column, 
                 replicate_column, feat_column, feat, 
                 cluster = cluster, cmethod = cmethod, K = K, 
                 num_nonzero = num_nonzero, thresh = thresh,  
                 clust_min_num_replicate = clust_min_num_replicate, 
                 fpca_optns = fpca_optns, fclust_optns = fclust_optns)
    })
  }
  names(res) <- features
  return(res)
}


fit_fdapca <- function(
  df, time_column, value_column, replicate_column,
  feat_column = NULL, feat = NULL,
  cluster = FALSE, cmethod = "EMCluster", K = 2, filter = TRUE,
  thresh = 0, num_nonzero = 10, clust_min_num_replicate = 15,
  fpca_optns = NULL, fclust_optns = NULL){
  if(!is.null(feat) & !is.null(feat_column)) {
    df <- df %>% 
      filter_at(.vars = feat_column, any_vars((.) == feat))
  }
  y_lst <- plyr::dlply(df, replicate_column, function(x) x[[value_column]])
  t_lst <- plyr::dlply(df, replicate_column, function(x) x[[time_column]])
  if (filter){
    idx <- sapply(y_lst, function(y){sum(abs(y) > thresh) >= num_nonzero})
    y_lst <- y_lst[idx]
    t_lst <- t_lst[idx]
  }  
  if(length(y_lst) == 0) return(NULL)
  feat_res <- NULL
  if(cluster & (length(y_lst) >= clust_min_num_replicate)) {
    feat_res <- try(fdapace::FClust(
      y_lst, t_lst, k = K, optnsFPCA = fpca_optns, optnsCS = fclust_optns))
  } else {
    if(is.null(fpca_optns)) {
      fpca_optns <- list()
    } 
    if (length(y_lst) < clust_min_num_replicate) {
      fpca_optns$dataType <- 'Sparse'
    }
    feat_res <- fdapace::FPCA(y_lst, t_lst, optns = fpca_optns)
    if(cluster){
      feat_res <- list("cluster" = rep(1, length(y_lst)), 
                       "fpca" = feat_res, "clusterObj" = NULL)
    }
  }
  feat_res
}

get_fpca_means <- function(object_lst, out_long_format = TRUE){
  if("cluster" %in% names(object_lst[[1]])){
    object_lst <- lapply(object_lst, function(x) x[["fpca"]])
    names(object_lst) <- names(object_lst)
  }
  mu <- plyr::ldply(object_lst, .fun = function(obj) {
    mu <- data.frame(t(obj[["mu"]]))
  }, .id = "Feature_ID") 
  colnames(mu) <- c("Feature_ID", object_lst[[1]]$workGrid)
  if(out_long_format) {
    mu <- mu %>% 
      gather(key = "time", value = "value", -Feature_ID) %>%
      mutate(
        time = as.numeric(time),
        Feature_ID = as.character(Feature_ID)) %>%
      arrange(Feature_ID, time)
  }
  return(mu)
}

get_fpca_fits <- function(object_lst, derOptns = list(p = 0)){
  fit_df <- plyr::ldply(seq_along(object_lst), function(i) {
    feat <- names(object_lst)[i]
    obj <- object_lst[[i]]
    fit <- fpca_fit(obj, derOptns) %>%
      mutate(Feature_ID = feat)
    return(fit)
  }, .id = NULL)
  return(fit_df)
}


fpca_fit <- function(obj, derOptns = list(p = 0)) {
  selectedK <- NULL
  clust_df <- NULL
  if("cluster" %in% names(obj)){
    clust_df <- data.frame(
      "Replicate_ID" = names(obj$fpca$inputData$Ly),
      "Cluster" =  as.character(obj[["cluster"]]),
      stringsAsFactors = FALSE)
    obj <- obj[["fpca"]]
  }
  if(derOptns$p > 0) {
    selectedK <- fdapace::SelectK(obj)$K
    if(!is.finite(selectedK)) selectedK <- ncol(obj$xiEst)
  }  
  fit <- fdapace:::fitted.FPCA(obj, K = selectedK, derOptns = derOptns)
  rownames(fit) <- names(obj$inputData$Ly)
  colnames(fit) <- obj$workGrid
  fit <- reshape2::melt(fit, varnames = c("Replicate_ID", "time")) %>%
    mutate(Replicate_ID = as.character(Replicate_ID))
  if(!is.null(clust_df)){
    fit <- suppressMessages(fit  %>% left_join(clust_df))
  }
  return(fit)
}


get_recovery_stats <- function(
  df, replicate_column = "Subject",
  stab_tol =1e-3, tolerance = 1.5, n_contigeous = 3){
  df$replicate <- df[[replicate_column]]
  df_lst <- lapply(unique(df$replicate), function(r){
    df_rep <- df %>% 
      filter(replicate == r) %>%
      mutate(
        sd_prior_value = sd(value[time <= -10]),
        diff = value - mean(value[time < -10]),
        stab = abs(deriv) < stab_tol) %>%
      arrange(time)
    df_rep$stab <- sapply(1:nrow(df_rep), function(i) {
      if(df_rep$time[i] <= 0) return(FALSE) 
      all(df_rep$stab[i:(min(nrow(df_rep), i+n_contigeous-1))]) 
    })
    df_rep <- df_rep %>%
      mutate(
        diff_max = diff[which.max(abs(diff))], 
        diff_t40 = diff[which.min(abs(time - 40))],
        slope_t4 = deriv[which.min(abs(time - 4))],
        stabilization_time = ifelse(any(stab), min(time[stab]), NA),
        stable = time >= stabilization_time,
        reach_stability = any(stable),
        recovered = stable & (abs(diff) < tolerance*sd_prior_value),
        recovery_time = ifelse(any(recovered), min(time[recovered]), NA),
        recovered = time >= recovery_time,
        recovery = any(recovered)
      )
    return(df_rep)
  })
  df_recov <- plyr::ldply(df_lst, function(x) x)
  return(df_recov)  
}


