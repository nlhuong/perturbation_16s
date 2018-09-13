# (...) - a list options of options to fdapace()
fda_pca <- function(df.long, time_column, value_column, replicate_column,
                    feat_column = NULL, parallel = FALSE, ncores = 1, 
                    num_nonzero = 10, thresh = 0, 
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
      fit_fdapca(df.long, time_column, value_column, replicate_column, 
                 feat_column, feat, num_nonzero = num_nonzero, thresh = thresh,  
                 cluster = cluster, cmethod = cmethod, K = K, 
                 fpca_optns = fpca_optns, fclust_optns = fclust_optns)
    }
  } else {
    res <- lapply(features, function(feat) {
      cat("Processing feature: ", feat, "...\n")
      fit_fdapca(df.long, time_column, value_column, replicate_column, 
                 feat_column, feat, num_nonzero = num_nonzero, thresh = thresh,  
                 cluster = cluster, cmethod = cmethod, K = K, 
                 fpca_optns = fpca_optns, fclust_optns = fclust_optns)
    })
  }
  names(res) <- features
  return(res)
}


fit_fdapca <- function(
  df, time_column, value_column, replicate_column,
  feat_column = NULL, feat = NULL, num_nonzero = 10, thresh = 0, 
  cluster = FALSE, cmethod = "EMCluster", K = 2, 
  fpca_optns = NULL, fclust_optns = NULL){
  if(!is.null(feat) & !is.null(feat_column)) {
    df <- df %>% 
      filter_at(.vars = feat_column, any_vars((.) == feat))
  }
  y_lst <- plyr::dlply(df, replicate_column, function(x) x[[value_column]])
  t_lst <- plyr::dlply(df, replicate_column, function(x) x[[time_column]])
  idx <- sapply(y_lst, function(y){sum(abs(y) > thresh) >= num_nonzero})
  feat_res <- NULL
  if(sum(idx) > 0) {
    y_lst <- y_lst[idx]
    t_lst <- t_lst[idx]
    if(cluster) {
      feat_res <- fdapace::FClust(
        y_lst, t_lst, k = K, optnsFPCA = fpca_optns, optnsCS = fclust_optns)
    } else {
      if(is.null(fpca_optns)) fpca_optns <- list()
      feat_res <- fdapace::FPCA(y_lst, t_lst, optns = fpca_optns)
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



