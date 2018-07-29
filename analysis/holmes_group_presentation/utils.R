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
