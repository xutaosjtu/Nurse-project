
plot_mean_se = function(data, metabolite){
  
  mean_se = function(x){
    m = mean(x, na.rm = T)
    SD = sd(x, na.rm = T)
    count = sum(!is.na(x))
    se = SD/sqrt(count)
    return(c(m, se))
  }
  
  tmp = aggregate(data[,metabolite], by=list(data$shift, data$Kontrolle,data$date2, as.character(data$hour)), mean_se)
  
  tmp = data.frame(tmp[,1:4], tmp$x)
  colnames(tmp) = c("shift", "group", "date", "time", "m", "se")
  
  #require(ggplot2)
  
  limits = aes(ymax = m+se, ymin = m - se)
  p = ggplot(tmp, aes(x = time, y = m))
  p = p + geom_point(aes(colour = interaction(as.factor(group), as.factor(shift)))) + geom_errorbar(limits, width = 0.4) + scale_color_manual(values = c("blue", "orange", "red"))
  p = p + ylab("Relative concentration") + ggtitle(metabolite)
  #p = p + facet_grid(shift + group~date)
  return(p)
}


