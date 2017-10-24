

library(BiSeq)


fileList = c(paste(path, "A-01_LE_13_713.cov",sep = "" ), 
             paste(path, "A-03_LE_13_716.cov",sep = "" ),
             paste(path, "A-11_LE_4_743.cov",sep = "" ),
             paste(path, "A-12_LE_28_787.cov",sep = "" ),
             paste(path, "A-14_LE_4_748.cov",sep = "" ),
             paste(path, "A-20_LE_4_744.cov",sep = "" ),
             paste(path, "A-22_LE_28_785.cov",sep = "" ),
             paste(path, "A-23_LE_13_712.cov",sep = "" ),
             paste(path, "B-29_LE_28_795.cov",sep = "" ),
             paste(path, "B-33_LE_4_747.cov",sep = "" ),
             paste(path, "B-35_LE_13_715.cov",sep = "" ),
             paste(path, "B-36_LE_28_799.cov",sep = "" ),
             
             paste(path, "A-07_Ancestral_0_824.cov",sep = "" ),
             paste(path, "A-08_Ancestral_0_817.cov",sep = "" ),
             paste(path, "A-10_Ancestral_0_829.cov",sep = "" ),
             paste(path, "A-15_Ancestral_0_812.cov",sep = "" ),
             paste(path, "A-21_Ancestral_0_830.cov",sep = "" ),
             
             paste(path, "B-30_Ancestral_0_811.cov",sep = "" ),
             paste(path, "B-34_Ancestral_0_823.cov",sep = "" ))



group = c("LE", "LE","LE","LE","LE","LE","LE","LE","LE","LE","LE","LE",
          "ANS","ANS","ANS","ANS","ANS","ANS","ANS")

colData <- DataFrame(group=group, 
                     row.names=c("A-01_LE","A-03_LE","A-11_LE","A-12_LE","A-14_LE",
                                 "A-20_LE","A-22_LE","A-23_LE","A-29_LE","A-33_LE",
                                 "A-35_LE","A-36_LE","A-07_Ancestral","A-08_Ancestral",
                                 "A-10_Ancestral","A-15_Ancestral","A-21_Ancestral",
                                 "A-30_Ancestral","A-34_Ancestral"))



BSraw_ANS_LE <- readBismark(fileList, colData)



colData(BSraw_ANS_LE) <- colData
BSraw_ANS_LE.small <- BSraw_ANS_LE[1:1000]

BSraw_ANS_LE.small.unlim <- clusterSites(object = BSraw_ANS_LE.small,
                                         groups = colData(BSraw_ANS_LE.small)$group,
                                         perc.samples = 1,
                                         min.sites = 20,
                                         max.dist = 100)

























