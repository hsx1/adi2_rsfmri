# Test fdr correction

set.seed(123)
x <- rnorm(20, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))

round(p, 3)
rank=seq(1,20)
round(p.adjust(p, "none"), 3)
round(p.adjust(p, "fdr"), 3)
plot(p, p.adjust(p, "fdr"))


pvalues <- c(.002, .005, .015, .113, .222, .227, .454, .552, .663, .751)
round(p.adjust(pvalues, "fdr"), 3)
plot(pvalues, p.adjust(pvalues, "fdr"))


p=c(0.002,0.04) 
round(p.adjust(p, "fdr"), 3) #-> beides signifikant
p=c(0.003, 0.12,0.4)
round(p.adjust(p, "fdr"), 3)#-> nur 1. signifikant

p=c(0.002,0.003,0.04,0.12,0.4)#-> nur beiden kleinsten signifikant
round(p.adjust(p, "fdr"), 3)
