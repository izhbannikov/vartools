#vartools - Variant Association tools R-package

```vartools``` offers nine varitan association tests: 

* CMC: Combined and Multivariate Collapsing Method for Rare Variants
* KBAC: Kernel Based Adaptive Clustering Method
* RBT: Replication Based Test for Protective Variants
* VT: Variable Thresholds Test for Case Control Data Analysis
* RareCover: A "Covering Algorithm" for Rare and Low Frequency Variants
* C(alpha): C-alpha Test for Protective Variants
* WSS: Weighted Sum Statistic via Rank Test
* aSum: Data-adaptive Sum Test for Protective and Deleterious Variants
* SKAT: SNP-set (Sequence) Kernel Association Test Method

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

### Other tests

In additional, ```vartools``` offers SKAT, WSS, aSum tests.Please check the user manual to use them.

## Contact

Ilya Y. Zhbannikov, i.zhbannikov@mail.ru, ilyaz@uidaho.edu