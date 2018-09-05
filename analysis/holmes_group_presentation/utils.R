plot_projection <- function(data, xname, yname, labname = "Subject", size = 3,
                            color = NULL, eigs = NULL, ...){
  if(all(!is.null(color), color %in% colnames(data))){
    plt <- ggplot(data, aes_string(x = xname, y = yname, color = color)) 
  } else {
    plt <- ggplot(data, aes_string(x = xname, y = yname)) 
  }
  if(all(!is.null(labname), labname %in% colnames(data))) {
    if(size %in% colnames(data)) {
      plt <- plt + geom_text(aes_string(label = labname, size = size), ...)
    } else {
      plt <- plt + geom_text(aes_string(label = labname), size = size, ...)
    }
  } else {
    if(size %in% colnames(data)) {
      plt <- plt + geom_point(aes_string(size = size), ...)
    } else {
      plt <- plt + geom_point(size = size, ...)
    }
  }
  if(!is.null(eigs)){
    var.explained <- round(100 * eigs/sum(eigs), 2)
    axis1 <- as.numeric(gsub("\\D", "", xname))
    axis2 <- as.numeric(gsub("\\D", "", yname))
    plt <- plt + 
      labs(
        x = sprintf("Axis%d [%s%% variance]", axis1, var.explained[axis1]),
        y = sprintf("Axis%d [%s%% variance]", axis2, var.explained[axis2])
      ) +
      coord_fixed(sqrt( eigs[axis2]/eigs[axis1])) 
  }
  return(plt + theme(text = element_text(size = 20))) 
}

plot_coefficients <- function (out.treeda, 
                               color = NULL,
                               remove.bl = TRUE, 
                               ladderize = TRUE, 
                               tree.height = 2) 
{
  predictor.names <- colnames(out.treeda$input$predictors)
  tr = out.treeda$input$tree
  if (remove.bl) {
    tr$edge.length = rep(1, length(tr$edge.length))
  }
  tree.plot = plot_tree(tr, ladderize = ladderize) + coord_flip() + 
    scale_x_reverse()
  leaf.position = get_leaf_position(out.treeda$input$tree, 
                                    ladderize = ladderize)$otu.pos
  coef = as(out.treeda$leafCoefficients$beta, "matrix")
  colnames(coef) = paste("Axis", 1:ncol(coef))
  df = data.frame(coef, leaf.position)
  df = reshape2::melt(df, id.vars = "leaf.position")
  df$predictor.names <- predictor.names 
  coef.plot = ggplot(df) +
    facet_grid(variable ~ .) + 
    ylab("Coefficient value") + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  if(!is.null(color)) {
    coef.plot <- coef.plot + 
      geom_point(aes(x = leaf.position, y = value, color = color))
  } else {
    coef.plot <- coef.plot + geom_point(aes(x = leaf.position, y = value))
  }
  p = combine_plot_and_tree(coef.plot, tree.plot, tree.height = tree.height, 
                            print = FALSE)
  grid::grid.draw(p)
  invisible(p)
  return(list("df" = df, "coef.plot" = coef.plot, "tree.plot" = tree.plot))
}

get_responses <- function(dist_df, dist_cols = c("bray", "jaccard")) {
  group_cols <- dist_df %>%
    select(-ends_with("S1"),  -dist_cols) %>%
    colnames()
  
  dist_df <- dist_df %>%
    group_by_at(group_cols) %>%
    summarise_at(
      .vars = dist_cols,
      .funs = mean) %>%
    arrange(Subject, DaysFromStart.S2) 
}


recovery_stats <- function(response_df, relday_col, 
                           perturb_day = 4, fac = 1.2,
                           thresh_num_days_below = 5,
                           dist_cols = c("jaccard", "bray")) {
  
  subject_nointv_thresh <- response_df %>%
    filter_at(.vars = relday_col, any_vars((.) < -7 )) %>%
    group_by(Subject) %>%
    summarise(jaccard = mean(jaccard), bray = mean(bray))
  
  response_df <- response_df %>% 
    as_data_frame() %>%
    left_join(subject_nointv_thresh, 
              by = "Subject", suffix = c("", "_thresh")) 
  
  recov_time <- response_df %>% 
    as_data_frame() %>%
    filter_at(.vars = relday_col, any_vars((.) > perturb_day))
  
  for(dname in dist_cols){
    dist_recov_time <- recov_time %>%
      mutate(
        below_thresh = (.)[[dname]] < fac * (.)[[paste0(dname, "_thresh")]]) %>%
      filter(below_thresh) %>%
      arrange_at(.vars = c("Subject", relday_col)) %>%
      group_by(Subject) %>%
      mutate(num_below_after = cumsum(below_thresh)) %>%
      group_by(Subject) %>%
      filter(num_below_after >= thresh_num_days_below) %>%
      summarise_at(.vars = relday_col, .funs = min)
    colnames(dist_recov_time) <- c("Subject", paste0(dname, "_recov_time"))
    
    response_df <- response_df %>%
      left_join(dist_recov_time)
  }
  
  return(response_df)
}
