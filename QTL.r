install.packages("qtl") 
library(qtl)
setwd("~/Desktop/Zurich_2017_computerlab/")
ls()
rm(list=ls()) #remove all items in your list to start fresh

data=read.cross(format=c("csv"), file="LerxCvix.csv", na.strings=c("NA"), genotypes=c("A","B")) 
#reads the data file in, tells R NA are missing data, and you've coded genotyeps as A and B
?read.cross

#tell rqtl that it is a recombinant inbred line (RIL), or a recombinant inbred selfed line
class(data)[1]="riself"
#### As it matters no of recombination in BC vs RIL

#summarize data
summary(data)

#plot all phenotypes and genetic map
plot(data)
par(mar=c(5.1,4.1,4.1,2.1))

#plot individual plots
plotPheno(data, pheno.col=3)

#make scatterplots to assess relationship between traits
plot(data$pheno$SD.FT,data$pheno$SD.RLN)

#single marker mapping#
lod1traitmr=scanone(data, pheno.col=3, model=c("normal"), method=c("mr")) #this is the first qtl model you run, using regression analysis
plot(lod1traitmr) #plots first qtl map
title("LDVer.RLN") #adds title to plot


#test for significance#
maxlod1traitmr=scanone(data, pheno.col=3, model=c("normal"), method=c("mr"), n.perm=1000) #same model as above, but using permutations to more accurately estimate the null hypothesis distribution
plot(maxlod1traitmr) #plot histogram of LOD scores from permutations
summary(maxlod1traitmr) #shows the 5% and 10% cutoff for LOD scores, ie any trait with LOD >5% cutoff would be significant at 5% level
plot(lod1traitmr)
add.threshold(lod1traitmr, perms=maxlod1traitmr, alpha=0.05) #adds a threshold line for the 5% cutoff value
summary(lod1traitmr, perms=maxlod1traitmr, alpha=0.05, pvalues=TRUE) #provides a table of QTLs that have LODs >5% value, shows their chromosome and position on the chromosome, and their significance


#Interval Mapping#
dataIM=calc.genoprob(data, step=1, map.function="kosambi") #model takes a guess at what is inbetween two known markers, and then makes a new dataset
lod1traitIM=scanone(dataIM, pheno.col=3, model=c("normal"), method=c("hk")) #run essentially the same scanone, but use this new dataset, but here we use a new method called Haley-Knott (hk)
plot(lod1traitIM)
plot(lod1traitIM, lod1traitmr) #plot both methods so far to compare. You will see there is no noticeable difference

#test for significance#
maxlod1traitIM=scanone(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=1000) #similar to above, we find out thresholds using permutation
plot(maxlod1traitIM) #plot a histogram of permutations
summary(maxlod1traitIM) #summarize your threshold results
plot(lod1traitIM) #plot qtl
add.threshold(lod1traitIM, perms=maxlod1traitIM, alpha=0.05) #now add the 5% threshold to qlt plot
summary(lod1traitIM, perms=maxlod1traitmr, alpha=0.05, pvalues=TRUE) #summarize qtls that are identified by interval mapping


#Composite Interval Mapping#
lod1traitCIM=cim(dataIM, pheno.col=3, n.marcovar=3, window=10, method=c("hk"), imp.method=c("imp"), error.prob=0.0001, map.function=c("kosambi")) #this is the sliding window method that essentially removes effects of "nearby" qtl, so as to not mask the focal qtl
plot(lod1traitCIM) #plot qtl map, and the differences look pretty big from regression or interval mapping methods

#try different numbers of background markers or window sizes#
#it's important to do this because you are trying different "window" sizes, and the number of possible other qtl to control for (n.marcovar)
lod1traitCIM6_10=cim(dataIM, pheno.col=3, n.marcovar=6, window=10, method=c("hk"), imp.method=c("imp"), error.prob=0.0001, map.function=c("kosambi"))
plot(lod1traitCIM6_10)

lod1traitCIM3_30=cim(dataIM, pheno.col=3, n.marcovar=3, window=30, method=c("hk"), imp.method=c("imp"), error.prob=0.0001, map.function=c("kosambi"))
plot(lod1traitCIM3_30)

lod1traitCIM3_80=cim(dataIM, pheno.col=3, n.marcovar=3, window=80, method=c("hk"), imp.method=c("imp"), error.prob=0.0001, map.function=c("kosambi"))
plot(lod1traitCIM3_80)

#compare results of different window and qtl numbers
plot(lod1traitCIM6_10, lod1traitCIM3_30, lod1traitCIM3_80)
legend("top", c("CIM6_10", "CIM3_30", "CIM3_80"), col=c("black", "blue","red"), lty=1)
title("testing different compositve interval options")



#test for significance#
maxlod1traitCIM=cim(dataIM, pheno.col=3, n.marcovar=3, window=20, method=c("hk"), imp.method=c("imp"), error.prob=0.0001, map.function=c("kosambi"), n.perm=1000) #run permutations using composite interval mapping to figure out error
plot(maxlod1traitCIM) #make a histogram of your permutations
summary(maxlod1traitCIM) #summarize the qtl that are found under permutation results
plot(lod1traitCIM) #make your qtl plot
add.threshold(lod1traitCIM, perms=maxlod1traitCIM, alpha=0.05) #add 5% threshold cutoff for significant qtl


#add other significance tresholds#
#so far we have only added 5% threshold, this can add further thresholds
add.threshold(lod1traitCIM, perms=maxlod1traitCIM, alpha=0.10, col="red")
add.threshold(lod1traitCIM, perms=maxlod1traitCIM, alpha=0.01, col="green")
legend("top", c("1%","5%","10%"), col=c("green", "black", "red"), lty=1) #add legend for each of these thresholds
title("different thresholds for composite interval mapping") #add plot title




##########################################################################################
##########################################################################################

#get confidence interval for QTL location#
#using LOD drop#
lodint(lod1traitIM, chr=1, drop=1.5, lodcolumn=1)

#using bayesian credible interval#

bayesint(lod1traitIM, chr=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#using nonparametric bootstrap#
bootoutputtrait=scanoneboot(dataIM, pheno.col=3, chr=1, model="normal", method="hk", n.boot=100, verbose=TRUE)
plot(bootoutputtrait)
summary(bootoutputtrait)

#plot QTL effect
CRY2<-find.marker(dataIM, 1, 7.2)
effectplot(dataIM, pheno.col=3, mname1=CRY2)

#fit multiple QTL model
summary(lod1traitIM, perms=maxlod1traitIM, alpha=0.05, pvalues=TRUE)
# take out several QTLs and make QTL object
qc <- c("1", "4", "5")
qp <- c(7.12, 34.0, 108.00)
qtl <- makeqtl(dataIM, qc, qp, what="prob")
summary(qtl)
lod.add <- fitqtl(dataIM, pheno.col=3, qtl, formula=y~Q1+Q2+Q3, method="hk")
summary(lod.add)


#complete 2-dimensional, 2 QTL scan#
lod2trait=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"))
plot(lod2trait)

#use stepwise model building#

lod2trait=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"))
perm1=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm2=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm3=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm4=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm5=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm6=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm7=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm8=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm9=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm10=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm11=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm12=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm13=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm14=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm15=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm16=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm17=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm18=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm19=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm20=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm21=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm22=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm23=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm24=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm25=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm26=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm27=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm28=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm29=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm30=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm31=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm32=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm33=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm34=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm35=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm36=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm37=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm38=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm39=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
perm40=scantwo(dataIM, pheno.col=3, model=c("normal"), method=c("hk"), n.perm=25)
maxlod2trait=c(perm1, perm2, perm3, perm4, perm5, perm6, perm7, perm8, perm9, perm10, perm11, perm12, perm13, perm14, perm15, perm16, perm17, perm18, perm19,perm20, perm21, perm22, perm23, perm24, perm25, perm26, perm27, perm28, perm29, perm30, perm31, perm32, perm33, perm34, perm35, perm36, perm37, perm38, perm39, perm40) 

save(maxlod2trait, file="permtest.R")
load("permtest.R")

summary(maxlod2trait)
pens<-calc.penalties(maxlod2trait)
traitstepwise <- stepwiseqtl(dataIM, pheno.col=3, max.qtl=6, method="hk", penalties=pens, keeplodprofile=TRUE, keeptrace=TRUE)
plotLodProfile(traitstepwise)
plotModel(traitstepwise)
