
#' @title Plot Mean Trajectory
#' @description Plot cluster-specific mean trajectory for one or all clusters provided by a \code{traj} object.
#' @param x \code{traj} object.
#' @param clust.num Integer indicating the cluster number to plot.\code{NULL} to print all clusters. Defaults to
#'\code{NULL}.
#' @param ... Arguments to be passed to plot.
#' 
#' @details The function plots cluster specific mean trajectory calculated at each time point. By setting the \code{clust.num} argument to an integer corresponding to a cluster number, one can plot the mean trajectory of that cluster only. Any other plotting arguments can be added to the function.
#' 
#' @author Marie-Pierre Sylvestre, Dan Vatnik
#'
#'marie-pierre.sylvestre@umontreal.ca
#'
#' @examples 
#' \dontrun{
#' # Setup data 
#' data = example.data$data
#'
#' # Run step1measures, step2factors and step3clusters with a predetermined number of clusters
#' s1 = step1measures(data, ID=TRUE)
#' s2 = step2factors(s1)
#' s3.4clusters = step3clusters(s2, nclust = 4)
#'
#' # Plot mean trajectories
#' plotMeanTraj(s3.4clusters)
#'}
#' @rdname plotMeanTraj
#' @export 
plotMeanTraj <- function(x, clust.num = NULL,  ...)
{
  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)


  data = data[,-1]
  time = time[,-1]

  unique.clusters = sort(unique(clust[,2]))

  num.clust = length(unique.clusters)


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

    comp.mean = data.frame(matrix(ncol = 3, nrow = 1))

    for(i_time in unique.clust.time)
    {
      time.pos = which(clust.time == i_time)

      data.time = clust.data[time.pos]

      mean.data = mean(data.time)
      centiles = 0

      comp.mean = rbind(comp.mean, c(mean.data, centiles[1], centiles[2]))
    }

    comp.mean = comp.mean[-1,]

    # Setup args for plotting
    args[["x"]] = unique.clust.time
    args[["y"]] = comp.mean[,1]
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

    lines(unique.clust.time, comp.mean[,1])

  }

  if(is.null(clust.num))
  {
    text = "Mean for Every Cluster"
    mtext(text, side = 3, line = 0.5, outer = TRUE, cex = 1.2)
  }

  if(is.null(clust.num)){
    par(mfrow = c(1,1))
  }

}
