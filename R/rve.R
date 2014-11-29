rve <- function(data) {
  # Performs a Kruskal-Wallis rank sum test.
  stat <- kruskal.test(data) 
  stat
}