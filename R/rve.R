rve <- function(x, y=NULL) {
  # Performs a one-sided  Wilcox test:
  stat <- wilcox.test(x, y, alternative = "g")        # greater
  stat
}