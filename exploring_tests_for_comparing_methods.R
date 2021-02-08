method1 <- setNames(c(20,89), c("ca", "nonca"))
method2 <- setNames(c(54,253), c("ca", "nonca"))
results <- as.matrix(rbind(method1, method2))
prop.test(results)$p.val # 0.9742808
prop.test(c(20,54), c(20+89,54+253))$p.val # 0.9742808
prop.test(results, correct=F)$p.val # 0.8587072
chisq.test(results)$p.val # 0.9742828 Pearson's chi-squared test with Yates's correction
chisq.test(results, correct=F)$p.val # 0.8587072 Pearson's chi-squared test
# prop.test takes a matrix or table with two columns - successes and failures - or a 
# vector of successes and a vector of trials.
# chisq.test performs Pearson's chi-squared test with the null hypothesis that the joint 
# distribution of the cell counts is the product of the row and column marginals.


method1 <- c(1,5,9,14,56,58,96,99,101,122,123,158,499,688,1552,1574,1699,1822,2012,2511,2600,2992)
method2 <- c(7,8,10,12,17,19,22,64,251,259,633,699,705,785,799,800,801,850,899,901,999,1110)
mean(method1)
mean(method2)
wilcox.test(method1, method2, paired=F)
wilcox.test(method1, method2, paired=T)


cases <- c(149,211,190,186,212,174,127,114,150,116)
days <- 24:33
nsw <- data.frame("New cases" = cases, 
                  "Date" = days)
library(ggplot2)
plot <- ggplot(data = nsw, 
               aes(x = Date, 
                   y = New.cases)) + 
  geom_line(size=2) + 
  ylab("New cases")
plot



unpaired <- matrix(c(41,33,44,52), 2, 2, 
                   dimnames=list(method=c("A", "B"), 
                                 result=c("+", "-"))
)
paired <- matrix(c(35,17,5,28), 2, 2, 
                 dimnames=list(methodA=c("+", "-"), 
                               methodB=c("+", "-"))
)
unpaired
paired
fisher.test(unpaired)$p.val
# 0.28; don't reject H_0 method independent of result
# valid test but not using pairing
fisher.test(paired)$p.val
# 2e-6; reject H_0 methodA independent of methodB
# valid test but not what we're interested in
mcnemar.test(unpaired)$p.val
# 0.25; don't reject H_0 (???)
# not valid test
mcnemar.test(paired)$p.val
# 0.02; reject H_0 f(result|A) = f(result|B)
# valid test, using pairing
chisq.test(unpaired)$p.val
# 0.28; don't reject H_0 method independent of result or f(result|A) = f(result|B)
chisq.test(paired)$p.val
# 8e-6; reject H_0 methodA independent of methodB


listA <- sample(1:100)
listB <- sample(1:100)
sum((listA-listB)^2)
cor.test(listA, listB, method="spearman")
listC <- listA[100:1]
sum((listA-listC)^2)
cor.test(listA, listC, method="spearman")
listD <- listA + 10
sum((listA-listD)^2)
cor.test(listA, listD, method="spearman")
pvals <- numeric(1000)
for (i in 1:1000) {
  listE <- sample(1:1000)[1:100]
  listF <- sample(1:1000)[1:100]
  pvals[i] <- cor.test(listE, listF, method="spearman")$p.val
}
hist(pvals)

library(OrderedList)
listG <- sample(1:1000)
listH <- sample(1:1000)
listI <- listG[1000:1]
GH <- compareLists(listG, listH)
GI <- compareLists(listG, listI)
getOverlap(GH)
getOverlap(GI)
cor.test(listG, listH, method="spearman", alt="less")
cor.test(listG, listI, method="spearman", alt="less")
data(OL.data)
list1 <- as.character(OL.data$map$prostate)
list2 <- c(sample(list1[1:500]), sample(list1[501:1000]))
y <- compareLists(list1, list2, two.sided=F)
z <- getOverlap(y)
z
ranks1 <- 1:length(list1)
ranks2 <- match(list1, list2)
cor.test(ranks1[1:100], ranks2[1:100], method="spearman", alt="less")

overlap <- numeric(1000)
p.overlap <- numeric(1000)
p.cor12 <- numeric(1000)
p.cor21 <- numeric(1000)
for (i in 1:1000) {
  overlap[i] <- mean(list1[1:i] %in% list2[1:i])
  perm.overlap <- numeric(1000)
  for (r in 1:1000) {
    list2null <- sample(list2)
    perm.overlap[r] <- mean(list1[1:i] %in% list2null[1:i])
  }
  p.overlap[i] <- mean(perm.overlap <= overlap[i])
  if (i > 1) {
    p.cor12[i] <- cor.test(1:i, match(list1, list2)[1:i], method="spearman", alt="less")$p.val
    p.cor21[i] <- cor.test(1:i, match(list2, list1)[1:i], method="spearman", alt="less")$p.val
  }
}
plot(overlap)
lines(p.overlap)
lines(p.cor12, col='red')
lines(p.cor21, col='blue')


table0 <- matrix(c(0,2,2,2),2,2)
table1 <- matrix(c(1,1,1,3),2,2)
table2 <- matrix(c(2,0,0,4),2,2)
fisher.test(table0, alt="less")$p # 0.4
fisher.test(table1, alt="less")$p # 0.9333333 i.e. 14/15
fisher.test(table2, alt="less")$p # 1
phyper(0,2,4,2) # 0.4
phyper(1,2,4,2) # 0.9333333 i.e. 14/15
phyper(2,2,4,2) # 1
# hypergeometric/one-sided Fisher's exact test for 1 through to whatever draws
# phyper(number of genes common to both lists, length of list, number of genes not in list, length of list)
pvals <- numeric(length(list1))
p.cor <- numeric(length(list1))
for (i in seq_len(length(list1))) {
  pvals[i] <- phyper(sum(list1[1:i] %in% list2[1:i]), i, length(list1) - i, i)
  genes <- unique(c(list1[1:i], list2[1:i]))
  p.cor[i] <- cor.test(match(genes, list1), match(genes, list2), method="spearman", alt="less")$p.val
}
plot(pvals, type='l', ylim=c(0,1))
lines(p.overlap, lty=2)
lines(p.cor, lty=3)
# p-values from hypergeometric test nearly identical to permutation-based p-values from overlap. That gives 
# me confidence that the hypergeometric test is testing what I want to test.

i<-50
phyper(sum(list1[1:i] %in% list2[1:i]), i, length(list1) - i, i)
fisher.test(matrix(c(
  sum(list1[1:i] %in% list2[1:i]), 
  i - sum(list1[1:i] %in% list2[1:i]), 
  i - sum(list1[1:i] %in% list2[1:i]), 
  length(list1) - 2*i + sum(list1[1:i] %in% list2[1:i])
),2,2), alt="less")$p.val
# same p-value for hypergeometric test and one-sided Fisher's exact test, which is also reassuring.


## Does permutation test on average squared distance between ranks of same gene in each list test the same 
## thing as correlation test?
# i.e. does permutation test give same results as 
# cor.test(match(genes, list1), match(genes, list2), method="spearman", alt="less")$p.val?
p.rankdiff <- numeric(1000)
for (i in 1:1000) {
  perm.rankdiff <- numeric(1000)
  genes <- unique(c(list1[1:i], list2[1:i]))
  for (r in 1:1000) {
    perm.rankdiff[r] <- sum((rank(match(genes, list1)) - sample(1:length(genes)))^2)
  }
  p.rankdiff[i] <- mean(perm.rankdiff >= sum((rank(match(genes, list1)) - rank(match(genes, list2)))^2))
}
plot(p.rankdiff, ylim=c(0,1), type="l")
lines(p.cor, lty=2)
# Yes.

