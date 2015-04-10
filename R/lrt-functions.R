lrtStat <- function(resNull, resFull) {
  rss.full <-  rowSums(resFull ^ 2)
  rss.null <- rowSums(resNull ^ 2)

  # F-statistic
  stat <- (rss.null - rss.full) / (rss.full)
  return(stat)
}
