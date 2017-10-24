####### Multiclass Diffential Expression Analysis
set.seed(10)
G <- 5
N <- 30
M <- 1000
initmean <- 5
initvar <- 10
mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
rownames(mat) <- paste0('gene', 1:M)
colnames(mat) <- paste0('cell', 1:(N*G))
group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
names(group) <- colnames(mat)
heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


### Now, let’s simulate some differentially upregulated genes unique to each group.
set.seed(10)
upreg <- 5
upregvar <- 10
ng <- 100

diff <- lapply(1:G, function(x) {
  diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
  mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})
names(diff) <- paste0('group', 1:G)
heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


# Let’s also simulate some differentially upregulated genes affecting groups of groups. 
# This is typical of cell differentiation processes.

diff2 <- lapply(2:(G-1), function(x) {
  y <- x+G
  diff <- rownames(mat)[(((y-1)*ng)+1):(((y-1)*ng)+ng)]
  mat[diff, group %in% paste0("group", 1:x)] <<- mat[diff, group %in% paste0("group", 1:x)] + rnorm(ng, upreg, upregvar)
  return(diff)
})

heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


### 1 vs many (in group vs. not in group)
# One approach is to consider each group vs. all others. We can use a simple T-test to test whether 
# each gene is significantly upregulated in each group vs. all other groups.
pv.sig <- lapply(levels(group), function(g){
  ingroup <- names(group)[group %in% g]
  outgroup <- names(group)[!(group %in% g)]
  pv <- sapply(1:M, function(i) {
    t.test(mat[i,ingroup], mat[i,outgroup], alternative='greater')$p.value
    #t.test(mat[i,ingroup], mat[i,outgroup])$p.value
  })
  names(pv) <- rownames(mat)
  pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni
  pv.sig
})

heatmap(mat[unique(unlist(pv.sig)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)



### Indeed, we pick up many of the differentially upregulated genes we simulated that are unique
### to each group ie. in diff. However, we have a much harder time picking up the differentially
### upregulated genes marking multiple groups, ie. in diff2 as expected.


perf <- function(pv.sig) {
  predtrue <- unique(unlist(pv.sig)) ## predicted differentially expressed genes
  predfalse <- setdiff(rownames(mat), predtrue) ## predicted not differentially expressed genes
  true <- c(unlist(diff), unlist(diff2)) ## true differentially expressed genes
  false <- setdiff(rownames(mat), true) ## true not differentially expressed genes
  TP <- sum(predtrue %in% true)
  TN <- sum(predfalse %in% false)
  FP <- sum(predtrue %in% false)
  FN <- sum(predfalse %in% true)
  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  prec <- TP/(TP+FP)
  fdr <- FP/(TP+FP)
  acc <- (TP+TN)/(TP+FP+FN+TN)
  return(data.frame(sens, spec, prec, fdr, acc))
}
print(perf(pv.sig))

#####################################################
########################## ANOVA ####################
#####################################################
## Alternatively, we can use ANOVA (analysis of variance) to identify genes that are variable 
## among and between group
pv <- sapply(1:M, function(i) {
  mydataframe <- data.frame(y=mat[i,], ig=group)
  fit <- aov(y ~ ig, data=mydataframe)
  summary(fit)[[1]][["Pr(>F)"]][1]
})
names(pv) <- rownames(mat)
pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni

heatmap(mat[pv.sig,], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)
print(perf(pv.sig))




##### Testing combined groups (pairs, triplets, etc)
pv.sig.pair <- lapply(levels(group), function(g1) {
  g2 <- setdiff(levels(group), g1)
  unlist(lapply(g2, function(g) {
    ## test two groups vs. all others
    ingroup <- names(group)[group %in% c(g1, g)]
    outgroup <- names(group)[!(group %in% c(g1, g))]
    pv <- sapply(1:M, function(i) {
      t.test(mat[i,ingroup], mat[i,outgroup], alternative='greater')$p.value
    })
    names(pv) <- rownames(mat)
    pv.sig <- names(pv)[pv < 0.05/M/G] ## bonferonni
    pv.sig
  }))
})

heatmap(mat[unique(unlist(pv.sig.pair)),], Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)


print(perf(pv.sig.pair))










