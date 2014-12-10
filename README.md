#vartools - Variant Association tools R-package

## Installation

```R CMD INSTALL vartools```

## Usage

###CMC test:
```
library(vartools)
pgdata <- as.matrix(read.table(system.file("extdata","phengen.dat",package="vartools"), as.is=T, skip = 1))
cmc.pvalue <- cmc(table=pgdata)
print(cmc.pvalue)
```

###KBAC test:
```
library(vartools)
alpha <- 0.05
num.permutation <- 3000
quiet <- 1
alternative <- 1
maf.upper.bound <- 0.05

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1) 
kbac.pvalue <- kbac(table=casectrl.dat, alpha, num.permutation, quiet, maf.upper.bound, alternative)
print(kbac.pvalue)
```

### C-alpha test:
```
library(vartools)

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)    
calpha.stat <- calpha(casectrl.dat)
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

###VT
```
library(vartools)
?vt

casectrl.dat <- read.table(system.file("extdata","phengen2.dat",package="vartools"), skip = 1)
vt.stat <- vt(casectrl.dat)
vt.stat

```

###RareCover
```
library(vartools)
?rarecover

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)
rarecover.stat <- rarecover(casectrl.dat)
rarecover.stat

```

###RBT
```
library(vartools)
?rbt

casectrl.dat <- read.table(system.file("extdata","phengen.dat",package="vartools"), skip = 1)
rbt.stat <- rbt(casectrl.dat)
rbt.stat

```

## Contact

Ilya Y. Zhbannikov, i.zhbannikov@mail.ru, ilyaz@uidaho.edu