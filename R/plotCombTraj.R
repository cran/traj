
#' @title Plot Cluster-Specific Mean or Median Trajectories provided by a \code{traj} object
#' @description Plot cluster-specific mean or median trajectories.
#' @param x \code{traj} object.
#' @param stat.type Choice between "mean" or "median". The mean or the median calculated at each time point for a cluster-specific set of trajectories will be plotted.
#'Defaults to "mean."
#' @param colored Boolean in dictating if the plot should use colors. If not, the trajectory lines will be distinctively patterned.
#'Defaults to \code{FALSE} for patterns.
#' @param ... Any extra parameter used by the \code{plot} function.
#' 
#' @details The function plots the mean or the median cluster-specific trajectory, calculated at each time point. A legend is generated in the
#'top left corner of the plot. Other plotting parameter(s) can be added to the function with the use of \dots.
#'@seealso 
#'\code{\link[base]{mean}}
#'\code{\link[stats]{median}}
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
#' # Plot mean combination trajectories
#' plotCombTraj(s3.4clusters)
#'}
#'
#' @rdname plotCombTraj
#' @export 
plotCombTraj <- function(x, stat.type = "mean", colored = FALSE, ...)
{

  clust = x$clusters
  data = as.matrix(x$data)
  time = as.matrix(x$time)


  data = data[,-1]
  time = time[,-1]

  unique.clusters = sort(unique(clust[,2]))

  num.clust = length(unique.clusters)

  min.y = min(data, na.rm=T)
  max.y = max(data, na.rm=T)

  min.x = min(time, na.rm=T)
  max.x = max(time, na.rm=T)


  max.y = max.y + num.clust * 0.04 * (max.y - min.y)

  xlab = "time"
  ylab = "y"
  args = list(...)
  params = names(args)


  for(i_clust in unique.clusters)
  {

    clust.data.pos = which(clust$cluster == i_clust)

    clust.data = data[clust.data.pos,]

    clust.time = time[clust.data.pos,]

    unique.clust.time = as.numeric(names(table(clust.time)))

    comp.stat = data.frame(matrix(ncol = 3, nrow = 1))

    for(i_time in unique.clust.time)
    {
      time.pos = which(clust.time == i_time)

      data.time = clust.data[time.pos]

      if(stat.type == "mean"){
        stat.data = mean(data.time)
        centiles = 0

        comp.stat = rbind(comp.stat, c(stat.data, centiles[1], centiles[2]))
      }
      else if(stat.type == "median"){
        stat.data = median(data.time)
        centiles = 0

        comp.stat = rbind(comp.stat, c(stat.data, centiles[1], centiles[2]))
      }
    }

    comp.stat = comp.stat[-1,]

    stat.name = paste(toupper(substr(stat.type, 1, 1)), substr(stat.type, 2, nchar(stat.type)), sep = "")

    # Setup args for plotting
    args[["x"]] = unique.clust.time
    args[["y"]] = comp.stat[,1]
    if(!("xlim" %in% params))
      args[["xlim"]] = c(min.x, max.x)
    if(!("ylim" %in% params))
      args[["ylim"]] = c(min.y, max.y + num.clust * 0.05 * (max.y - min.y) + 2 * 0.05 * (max.y - min.y))
    if(!("main" %in% params))
      args[["main"]] = paste(stat.name," Trajectory of all Clusters", sep = "")
    if(!("xlab" %in% params))
      args[["xlab"]] = xlab
    if(!("ylab" %in% params))
      args[["ylab"]] = ylab
    if(colored)
      args[["col"]] = i_clust

      args[["pch"]] = i_clust
    # Generates the cluster plot with the first sampled data.
    if(i_clust == 1)
      do.call(plot, args)
    else if(colored)
      points(unique.clust.time, comp.stat[,1], col = i_clust, pch = i_clust)
    else
      points(unique.clust.time, comp.stat[,1], pch = i_clust)

    if(colored)
      lines(unique.clust.time, comp.stat[,1], col = i_clust)
    else
      lines(unique.clust.time, comp.stat[,1], lty = i_clust)
  }

  if(colored){
    legend(min.x + 0.05 * (max.x - min.x), max.y + num.clust * 0.07 * (max.y - min.y),
           unique.clusters,
           lty=rep(1, num.clust),
           lwd= rep(2.5, num.clust),
           col=unique.clusters,
           pch = unique.clusters)
  }
  else{
    legend(min.x + 0.05 * (max.x - min.x), max.y + num.clust * 0.07 * (max.y - min.y),
           unique.clusters,
           lty=unique.clusters,
           lwd= rep(2.5, num.clust),
           col=rep(1, num.clust),
           pch = unique.clusters)
  }
}
