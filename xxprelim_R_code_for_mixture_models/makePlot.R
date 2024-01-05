makePlot <- function(df, fit){
  
  mu = fit$mu
  sigma = fit$sigma
  p = nrow(mu)
  k = ncol(mu)
  
  means = as.data.frame(t(mu))
  
  ggplot(df, aes(x=y1, y=y2, col=x)) + 
    geom_point() +
    geom_point(data = means, aes(x = V1, y = V2), col = 1, shape = 17, size = 2) +  
    labs(title = paste0('Fit with K = ',k)) +
    theme_classic() 
}