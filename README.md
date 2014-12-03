#vartools - Variant Association tools R-package

## Installation

```R CMD INSTALL vartools```

## Usage

###CMC test:
```
library(vartools)
cmc.pvalue <- cmc(table=system.file("extdata","phengen.dat",package="vartools"))
print(cmc.pvalue)
summary(cmc.pvalue)
```

###KBAC test:
```
library(vartools)
alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05
kbac.pvalue <- kbac(table=system.file("extdata","phengen.dat",package="vartools"), alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
summary(kbac.pvalue)
```

### C-alpha test:
```
library(vartools)

calpha.stat <- calpha(system.file("extdata","phengen.dat",package="vartools"))
calpha.stat
```

###RVE
```
library(vartools)

x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
rve.stat <- rve(x, y)
rve.stat
```

## Contact

Ilya Y. Zhbannikov, i.zhbannikov@mail.ru, ilyaz@uidaho.edu