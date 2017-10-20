#library(devtools)
#install_github('sartorlab/methylSig')
library(methylSig)
require(parallel)
require(boot)


setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/methylSig/")
path <- "/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/Aphid/DSS_Analysis/methylSig"
fileList = c(paste(path, "A-01_LE_13_713.MethylKit",sep = "/" ), 
             paste(path, "A-03_LE_13_716.MethylKit",sep = "/" ),
             paste(path, "A-11_LE_4_743.MethylKit",sep = "/" ),
             paste(path, "A-12_LE_28_787.MethylKit",sep = "/" ),
             paste(path, "A-14_LE_4_748.MethylKit",sep = "/" ),
             paste(path, "A-20_LE_4_744.MethylKit",sep = "/" ),
             paste(path, "A-22_LE_28_785.MethylKit",sep = "/" ),
             paste(path, "A-23_LE_13_712.MethylKit",sep = "/" ),
             paste(path, "B-29_LE_28_795.MethylKit",sep = "/" ),
             paste(path, "B-33_LE_4_747.MethylKit",sep = "/" ),
             paste(path, "B-35_LE_13_715.MethylKit",sep = "/" ),
             paste(path, "B-36_LE_28_799.MethylKit",sep = "/" ),

             
             paste(path, "A-07_Ancestral_0_824.MethylKit",sep = "/" ),
             paste(path, "A-08_Ancestral_0_817.MethylKit",sep = "/" ),
             paste(path, "A-10_Ancestral_0_829.MethylKit",sep = "/" ),
             paste(path, "A-15_Ancestral_0_812.MethylKit",sep = "/" ),
             paste(path, "A-21_Ancestral_0_830.MethylKit",sep = "/" ),
             paste(path, "B-30_Ancestral_0_811.MethylKit",sep = "/" ),
             paste(path, "B-34_Ancestral_0_823.MethylKit",sep = "/" ))

sample.id = c("A-01_LE", "A-03_LE", "A-11_LE", "A-12_LE", "A-14_LE", "A-20_LE", "A-22_LE", "A-23_LE", "A-29_LE", "A-33_LE","A-35_LE", "A-36_LE",
              "A-07_Ancestral", "A-08_Ancestral", "A-10_Ancestral", "A-15_Ancestral", "A-21_Ancestral", "B-30_Ancestral", "B-34_Ancestral")
treatment = c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0)
#### Read Data ####
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "TAIR10",
                          treatment = treatment, context = "CpG", minCount=5, maxCount=100, quiet=TRUE, destranded=TRUE, num.cores = 4)

#### Differential methylation analysis, dispersion calculation from both groups
bothGroup_DMC <- methylSigCalc(meth, groups = c(Treatment = 1, Control = 0), dispersion = "both",
              local.disp = FALSE, winsize.disp = 50, min.per.group = c(5, 5), num.cores = 4)

#### Differential methylation analysis, dispersion calculation from one groups
oneGroup_DMC <- methylSigCalc(meth, groups = c(Treatment = 1, Control = 0), dispersion = 0,
              local.disp = FALSE, winsize.disp = 50, min.per.group = c(5, 5), num.cores = 4)

#myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=2, num.cores=4)

oneGDiffSigbothDMCs = oneGroup_DMC[oneGroup_DMC[,"qvalue"] <= 1.05 & abs(oneGroup_DMC[,"meth.diff"])>=20, ]
bothGDiffSigbothDMCs = bothGroup_DMC[bothGroup_DMC[,"qvalue"] <= 1.05 & abs(bothGroup_DMC[,"meth.diff"])>=20, ]



