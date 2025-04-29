#' Plot random field
#'
#' @param x \code{RandomField} object.
#' @param main title to the plot.
#' @param colors vector of colors over which a gradient is created.
#' @param name title of the legend, character string.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' @examples
#' plot(genField(c(50, 50)))
#' plot(genField(c(50, 50), type = 2))
#' 
#' @return a ggplot2 object. Function is mostly called for its side effect.
#'
#' @exportS3Method
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn ggtitle geom_raster guides guide_colorbar
plot.RandomField <- function(x, main = "", colors, name = "value", ...)
{
  n <- dim(x)
  
  if(is.null(n) || (n[1] == 1 | n[2] == 1))
  {
    plot.ts(x)
    return(invisible(NULL))
  }
  
  if(missing(colors)) colors <- c(rgb(0.95, 1, 0), rgb(1, 1, 0), rgb(1, 0.8, 0), 
                                  "orange", "red", "darkred", "black")
  
  ### binding parameters for automatic check reasons
  y <- NA
  z <- NA
  ###

  data <- cbind(expand.grid(y = 1:n[1], x = 1:n[2]), z = as.vector(x))
  p <- ggplot(data, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(colors = colors, name = name, ...) +
    ggtitle(main) +
    guides(fill = guide_colorbar(order = 1))

  return(p)
}


