#' Plot random field
#'
#' @param x \code{RandomField} object.
#' @param main title to the plot.
#' @param colors vector of colors over which a gradient is created.
#' @param name title of the legend, character string.
#' @param alpha opacity of the raster tiles. Numeric value between 0 and 1.
#' @param facet_id optional factors to divide plots into facets. Has to be a factor vector of length \code{prod(dim(x))} resp. \code{length(x)}.
#' @param p additional ggplot2 component added BEFORE this function's components.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' @examples
#' plot(genField(c(50, 50)))
#' plot(genField(c(50, 50), type = 2))
#' 
#' @return a ggplot2 object. Function is mostly called for its side effect.
#'
#' @exportS3Method
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn ggtitle geom_raster guides guide_colorbar facet_wrap geom_line theme element_rect element_line element_text 
plot.RandomField <- function(x, main = "", colors, name = "value", alpha = 1, facet_id = NULL,
                             p = NULL, ...)
{
  ### binding parameters for automatic check reasons
  y <- NA
  z <- NA
  ###
  
  n <- dim(x)
  
  if(is.null(n) || (n[1] == 1 | n[2] == 1))
  {
    if(is.null(n))
    {
      xvar <- names(x)
    } else if(n[1] == 1)
    {
      xvar <- colnames(x)
    } else if(n[2] == 1)
    {
      xvar <- rownames(x)
    }
    xvar <- as.numeric(xvar)
    
    if(length(xvar) == 0) xvar <- 1:length(x)
    data <- cbind(x = xvar, y = as.vector(x))
    
    if(missing(colors) || any(is.na(colors))) colors <- "black"
    
    if(!is.null(facet_id)) data <- cbind(data, id = facet_id)
    
    p <- ggplot(data, aes(x = x, y = y)) +
      p +
      geom_line(colour = colors) +
      ggtitle(main) + 
      theme(panel.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_line(color = "grey95"), 
            panel.grid.minor = element_line(color = "grey95"), 
            axis.line = element_line(color = "grey90"), 
            strip.background = element_rect(fill = "grey90"), 
            strip.text.y = element_text(angle = 0, size = 12),
            strip.text.x = element_text(size = 12))

    if(!is.null(facet_id)) p <- p + facet_wrap(~ id, scales = "free")
    
    return(p)
    # plot.ts(x)
    # return(invisible(NULL))
  }
  
  if(missing(colors) || any(is.na(colors)))
    colors <- c(rgb(0.95, 1, 0), rgb(1, 1, 0), rgb(1, 0.8, 0), "orange", 
                "red", "darkred", "black")
  
  xvar <- as.numeric(colnames(x))
  if(length(xvar) == 0) xvar <- 1:n[2]
  yvar <- as.numeric(rownames(x)) 
  if(length(yvar) == 0) yvar <- 1:n[1]

  data <- cbind(expand.grid(y = yvar, x = xvar), z = as.vector(x))
  
  if(!is.null(facet_id)) data <- cbind(data, id = facet_id)
  
  p <- ggplot(data, aes(x = x, y = y, fill = z)) +
    p +
    geom_raster(alpha = alpha) +
    scale_fill_gradientn(colors = colors, name = name, ...) +
    ggtitle(main) +
    guides(fill = guide_colorbar(order = 1))
  
  if(!is.null(facet_id)) p <- p + facet_wrap(~ id, scales = "free")

  return(p)
}


