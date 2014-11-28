rve <- function(table) {
  stat <- fisher.test(table)$p.value # Fisher's exact test
  stat
}