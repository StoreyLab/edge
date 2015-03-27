lrtStat <- function(resNull, resFull) {
  # lrt statistic
  #
  # Args:
  #   res0: matrix- Null residuals
  #   res1: matrix- Full residuals
  #
  # Returns:
  #   stat: vector- F statistic 
  # Calculates RSS
  rss.full <- rowSums(resFull ^ 2)
  rss.null <- rowSums(resNull ^ 2)
  
  # F-statistic
  stat <- (rss.null - rss.full) / (rss.full)  
  return(stat)
}