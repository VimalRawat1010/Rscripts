#http://bis.zju.edu.cn/BNArray/#OLE9

## Install package "deal"
library(BNArray)
data(total.data)
attach(total.data)



#    Impute missing values
### The current version of BNArray (V1.0) allows using Least Local Squares algorithm to estimate the missing values as the linear 
### combination of their nearest neighbors:
ori.compact <- LLSimpute(total.data$df.all, total.data$df.ori, total.data$n.changed)
bn.data <- PrepareCompData(ori.compact)


TF_Exp_FPKM$XYZ <- factor(rep("xyz", nrow(TF_Exp_FPKM)))
nw <- network(data.matrix(TF_Exp_FPKM))

###  Construct Bayesian network
# implement Bayesian network learning by using deal package with a slight modification.

nw <- network(bn.data)
nw.prior <- jointprior(nw,20)