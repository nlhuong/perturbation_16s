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
}



