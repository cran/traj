
#' @title Plot Median Trajectory of \code{traj} Object
#' @description Plot cluster-specific median trajectory for one or all clusters provided by a \code{traj} object.
#' @param x \code{traj} object.
#' @param clust.num Integer indicating the cluster number to plot.\code{NULL} to print all clusters. Defaults to
#'\code{NULL}.
#' @param plot.percentile Value indicating if the function should plot percentiles. Defaults to \code{TRUE}.
#' @param low.percentile Value of the lower percentile to be plotted. Must be between 0 and 1. Defaults to 0.1.
#' @param high.percentile Value of the high percentile to be plotted. Must be between 0 and 1. Defaults to 0.9.
#' @param ... Extra parameters used in the \code{plot} function.
#' 
#' @details The function plots cluster specific median trajectory calculated at each time point, in addition to 10th and 90th percentiles. By setting the \code{clust.num} argument to an integer corresponding to a cluster number, one can plot the median trajectory of that cluster only. Any other plotting arguments can be added to the function.
#'@seealso 
#'\code{\link[graphics]{plot}}
#'\code{\link[stats]{median}}
#'\code{\link[stats]{quantile}}
#'
#' @author Marie-Pierre Sylvestre, Dan Vatnik
#'
#'marie-pierre.sylvestre@umontreal.ca
#'
#' @examples 
#' \dontrun{
#' # Setup data and time
#' data = example.data$data
#' time = example.data$time
#'
#' # Run step1measures, step2factors and step3clusters with a predetermined number of clusters
#' s1 = step1measures(data,time, ID=TRUE)
#' s2 = step2factors(s1)
#' s3.4clusters = step3clusters(s2, nclust = 4)
#'
#' # Plot median trajectories
#' plotMedTraj(s3.4clusters)
#'}
#'
#' @rdname plotMedTraj
#' @export 
plotMedTraj <- function(x, clust.num = NULL, plot.percentile = TRUE, low.percentile = 0.10, high.percentile = 0.90, ...)
{

  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)


  data = data[,-1]
  time = time[,-1]

  unique.clusters = sort(unique(clust[,2]))

  num.clust = length(unique.clusters)

  # Set to default percentile
  # If plot.percentile = FALSE, they will not be plotted.
  if(is.null(low.percentile)){
    low.percentile = 0.10
  }
  if(is.null(high.percentile)){
    high.percentile = 0.90
  }


  # Create multiplot canvas
  if(is.null(clust.num))
  {
    num.plot.x = round(sqrt(num.clust))
    if(num.plot.x^2 < num.clust)
      num.plot.y = num.plot.x + 1
    else
      num.plot.y = num.plot.x

    par(mfrow=c(num.plot.y, num.plot.x), oma = c(0,0,4,0) )


  }
  else unique.clusters = clust.num


  min.y = min(data, na.rm=T)
  max.y = max(data, na.rm=T)

  min.x = min(time, na.rm=T)
  max.x = max(time, na.rm=T)

  xlab = "time"
  ylab = "y"
  args = list(...)
  params = names(args)



  # Loops through all clusters
  for(i_clust in unique.clusters)
  {

    clust.data.pos = which(clust$cluster == i_clust)

    clust.data = data[clust.data.pos,]

    clust.time = time[clust.data.pos,]

    unique.clust.time = as.numeric(names(table(clust.time)))

    comp.med = data.frame(matrix(ncol = 3, nrow = 1))

    for(i_time in unique.clust.time)
    {
      time.pos = which(clust.time == i_time)

      data.time = clust.data[time.pos]

      median.data = median(data.time)
      centiles = quantile(data.time,c(low.percentile, high.percentile))

      comp.med = rbind(comp.med, c(median.data, centiles[1], centiles[2]))
    }

    comp.med = comp.med[-1,]

    # Setup args for plotting
    args[["x"]] = unique.clust.time
    args[["y"]] = comp.med[,1]
    if(!("xlim" %in% params))
      args[["xlim"]] = c(min.x, max.x)
    if(!("ylim" %in% params))
      args[["ylim"]] = c(min.y, max.y)
    if(!("main" %in% params))
      args[["main"]] = paste("Cluster ", i_clust, sep = "")
    if(!("xlab" %in% params))
      args[["xlab"]] = xlab
    if(!("ylab" %in% params))
      args[["ylab"]] = ylab

    # Generates the cluster plot with the first sampled data.
    do.call(plot, args)

    lines(unique.clust.time, comp.med[,1])


    # Loops throuth the rest of the samples in cluster and plots them on the same plot.
    if(plot.percentile){
      for(i_centiles in 2:3)
      {
        points(unique.clust.time, comp.med[,i_centiles], col = "red")
        lines(unique.clust.time,comp.med[,i_centiles], col = "red", lty= 2)
      }
    }
  }

  if(is.null(clust.num))
  {
    text = "Median, 10% and 90% for Every Cluster"
    mtext(text, side = 3, line = 0.5, outer = TRUE, cex = 1.2)
  }

  if(is.null(clust.num)){
    par(mfrow = c(1,1))
  }

}
