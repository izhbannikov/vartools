#vartools - Variant Association tools R-package

## Installation

```R CMD INSTALL vartools```

```

## Usage

library(vartools)

#CMC test:
cmc.pvalue <- cmc(table=system.file("extdata","phengen.dat",package="vartools"))
print(cmc.pvalue)
summary(cmc.pvalue)

#KBAC test:
alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05
kbac.pvalue <- kbac(table=system.file("extdata","phengen.dat",package="vartools"), alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
summary(kbac.pvalue)
```

## Contact

Ilya Y. Zhbannikov, i.zhbannikov@mail.ru, ilyaz@uidaho.edu