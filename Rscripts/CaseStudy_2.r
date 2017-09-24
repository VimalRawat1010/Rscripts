rm(list =ls()) ### Clear the workspace
.libPaths("/home/vimal/R_Library/")
library("openxlsx")
library(ggplot2)
library(reshape)
library(knitr)

###  Setting home directory
setwd("/home/vimal/Software/Github/Rscripts/")

### Reading excel file using read.xls function from library(openxlsx) is comparitively faster
All_users <- read.xlsx("Analytics_case.xlsx", sheet = 2, startRow = 1, colNames = TRUE)
Redeemers<- read.xlsx("Analytics_case.xlsx", sheet = 3, startRow = 1, colNames = TRUE)
### changing date field to date feild type
All_users$Week = as.Date(All_users$Week, origin =  "1900-01-01") -2 

#----------------------------------------------------
# Task 1:
#-----------------------------------------------------
#Adding two factor column (givenCoupon, RedeemerGroup) to data frame "All_users"
All_users$usedCoupon <- ifelse(All_users$User.ID %in% Redeemers$User.ID,"TRUE","FALSE")
All_users$RedeemerGroup <- ifelse(All_users$Group =="Control", "Control", ifelse(All_users$Group =="Test" & All_users$usedCoupon =="TRUE","Redeemer", "Non-Redeemer"))
unique_ID = data.frame(unique(All_users[c("User.ID", "Group")]))
unique_ID_groups = data.frame(unique(All_users[c("User.ID", "RedeemerGroup")]))
#Table to see the user count for each group
table(unique_ID$Group)
table(unique_ID_groups$RedeemerGroup)


#### Adding extra column to Classify data as Nuetral period (7 weeks), PreLaunch (6 weeks), Launch(2 weeks), and post 
##@# Post Launch (6 weeks)
All_users$Period <- ifelse(All_users$Week < '2014-11-24', "Neutral",ifelse(All_users$Week <'2015-01-05', "PreLaunch",
                                              ifelse(All_users$Week <='2015-01-18', "Launch", "PostLaunch")))



#### Checking for potential *********************OUTLIERS***********************
#### Using Aggregate function to take average spending by each User/weeek : Control/Test categories
#aggregate(All_users$Spend, by=list(Category=All_users$Group, Week=All_users$Week), FUN=mean)
ggplot(All_users , aes(x = All_users$Week, y = All_users$Spend, fill= Group)) +  geom_bar(stat="identity", position=position_dodge())
ggplot(All_users , aes(x = Week, y =log10(Spend+1), colour =All_users$Group , group=User.ID)) +  geom_boxplot()
ggplot(All_users , aes(x = Week, y =log10(Spend+1), colour =All_users$Group , group=User.ID)) +  geom_point()


User_Spending_perweek <-setNames(aggregate(All_users$Spend, by=list(All_users$User.ID, All_users$Group,All_users$Week), FUN=sum) ,c( "User.ID","Group","Week","Spend"))
User_Spending_perweek$Period <- ifelse(User_Spending_perweek$Week < '2014-11-24', "Neutral",ifelse(User_Spending_perweek$Week <'2015-01-05', "PreLaunch",
                                                                           ifelse(User_Spending_perweek$Week <='2015-01-18', "Launch", "PostLaunch")))


### LogNormal distribution
hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="Neutral"])+1,xlim=c(-2,8), main="Control group spending Neutral period", col="red")
hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"     & User_Spending_perweek$Period =="Neutral"])+1, xlim=c(-2,8), col="Blue",main="Test group spending Neutral period",)

hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="PreLaunch"])+1, xlim=c(-2,8), main="Control group spending PreLaunch period")
hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="PreLaunch"])+1, xlim=c(-2,8), main="Test group spending PreLaunch period")

hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="Launch"])+1, xlim=c(-2,8), main="Control group spending Launch period")
hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="Launch"])+1, xlim=c(-2,8), main="Test group spending Launch period")

hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="PostLaunch"])+1, xlim=c(-2,8), main="Control group spending PostLaunch period")
hist(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="PostLaunch"])+1, xlim=c(-2,8), main="Test group spending PostLaunch period")
############ Ploting Densities
plot(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="Neutral"])+1), xlim=c(-2,8), ylim=c(0,1.0), main="Control (Red) & Test (Blue) group spending Neutral period", col="red")
lines(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"    & User_Spending_perweek$Period =="Neutral"])+1), col="Blue")

plot(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="PreLaunch"])+1), xlim=c(-2,8), ylim=c(0,1.0), main="Control (Red) & Test (Blue) group spending PreLaunch period", col="red")
lines(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="PreLaunch"])+1), col="Blue")

plot(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="Launch"])+1), xlim=c(-2,8), ylim=c(0,1.0),main="Control (Red) & Test (Blue) group spending Launch period", col="red")
lines(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="Launch"])+1), col="Blue")

plot(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Control"  & User_Spending_perweek$Period =="PostLaunch"])+1), xlim=c(-2,8), ylim=c(0,1.0), main="Control (Red) & Test (Blue) group spending PostLaunch period", col="red")
lines(density(log10(User_Spending_perweek$Spend[User_Spending_perweek$Group=="Test"  & User_Spending_perweek$Period =="PostLaunch"])+1), col="Blue")





library(extremevalues)


ggplot(All_users , aes(x = All_users$Week, y = All_users$Spend, fill= RedeemerGroup)) +  geom_bar(stat="identity", position=position_dodge())

#Spending_perweek_peruser_2grp <-setNames(aggregate(All_users$Spend, by=list(Category=All_users$User.ID, Week=All_users$Week), FUN=sum) ,c( "RedeemGroup","Week","Spend"))




#### Using Aggregate function to take average spending by each group/weeek
#aggregate(All_users$Spend, by=list(Category=All_users$Group, Week=All_users$Week), FUN=mean)
Spending_perweek_2grp <-setNames(aggregate(All_users$Spend, by=list(Category=All_users$Group, Week=All_users$Week), FUN=sum) ,c( "RedeemGroup","Week","Spend"))
#Normalize by group size
Spending_perweek_2grp$Spend <- ifelse(Spending_perweek_2grp$RedeemGroup =="Control",Spending_perweek_2grp$Spend/1001 ,Spending_perweek_2grp$Spend/3996)

Spending_perweek_3grp <-setNames(aggregate(All_users$Spend, by=list(Category=All_users$RedeemerGroup, Week=All_users$Week), FUN=sum) ,c( "RedeemGroup","Week","Spend"))
#Normalize by group size
Spending_perweek_3grp$Spend <- ifelse(Spending_perweek_3grp$RedeemGroup =="Control",Spending_perweek_3grp$Spend/1001 ,ifelse(Spending_perweek_3grp$RedeemGroup =="Non-Redeemer",Spending_perweek_3grp$Spend/3192,Spending_perweek_3grp$Spend/804) )



Spending_perweek_2grp$Period <- ifelse(Spending_perweek_2grp$Week < '2014-11-24', "Neutral", ifelse(Spending_perweek_2grp$Week <'2015-01-05', "PreLaunch",
                                                                                                    ifelse(Spending_perweek_2grp$Week <= '2015-01-18',"Launch", "PostLaunch")))
Spending_perweek_3grp$Period <- ifelse(Spending_perweek_3grp$Week < '2014-11-24', "Neutral", ifelse(Spending_perweek_3grp$Week <'2015-01-05', "PreLaunch",
                                                                                                    ifelse(Spending_perweek_3grp$Week <= '2015-01-18',"Launch", "PostLaunch")))


ggplot(Spending_perweek_2grp , aes(x = Spending_perweek_2grp$Week, y = Spending_perweek_2grp$Spend, 
                                  colour = RedeemGroup)) + geom_line(size=4) + xlab("Weeks")+ ylab("Avg spending") + 
                                  ylim(low=0,high=200) + xlim(as.Date(c('22/11/2014', '09/02/2015'), format="%d/%m/%Y")) +
                                 ggtitle("Avg. Spending Per Week") 

ggplot(Spending_perweek_3grp , aes(x =Spending_perweek_3grp$Week, y = Spending_perweek_3grp$Spend, 
                                   colour = RedeemGroup)) + geom_line(size=4) + xlab("Weeks")+ ylab("Avg spending") + 
                                   ylim(low=0,high=200) + xlim(as.Date(c('22/11/2014', '09/02/2015'), format="%d/%m/%Y"))+ ggtitle("Avg. Spending Per Week") 


Spending_perweek_2grp$Period <- factor(Spending_perweek_2grp$Period,levels = c('Neutral','PreLaunch','Launch','PostLaunch'),ordered = TRUE)
ggplot(Spending_perweek_2grp, aes(x=Period, y=Spend, colour = RedeemGroup)) +geom_boxplot() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")

Spending_perweek_3grp$Period <- factor(Spending_perweek_3grp$Period,levels = c('Neutral','PreLaunch','Launch','PostLaunch'),ordered = TRUE)
ggplot(Spending_perweek_3grp, aes(x=Period, y=Spend, colour = RedeemGroup)) +geom_boxplot() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")
ggplot(Spending_perweek_3grp, aes(x=Period, y=Spend, colour = RedeemGroup)) +geom_violin() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")

ggplot(All_users, aes(x=Period, y=log10(Spend+1), colour = RedeemerGroup)) +geom_boxplot() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")
#ggplot(All_users, aes(x=Period, y=log10(Spend+1), colour = RedeemerGroup)) +geom_violin() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")

#library(beeswarm)
#beeswarm(All_users$Spend~All_users$RedeemerGroup, pch=16, cex=.1, pwcol=as.numeric(sqrt(All_users$Spend)))



###------------------------------------------
###  TASK 3
###-------------------------------------------

ggplot(Spending_perweek_2grp , aes(x = Spending_perweek_2grp$Week, y = Spending_perweek_2grp$Spend, 
                                   colour = RedeemGroup)) + geom_line(size=4) + xlab("Weeks")+ ylab("Avg spending") + 
                                  ylim(low=0,high=200) + xlim(as.Date(c('19/01/2015', '23/02/2015'), format="%d/%m/%Y")) +
                                  ggtitle("Avg. Spending Per Week") 


count_test_W1 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-01-19' )))
count_test_W2 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-01-26' )))
count_test_W3 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-02-02' )))
count_test_W4 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-02-09' )))
count_test_W5 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-02-16' )))
count_test_W6 = length(unique(which(All_users$Group == "Test" & All_users$Week == '2015-02-23' )))


Cont_W1 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-01-19']
Cont_W2 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-01-26']
Cont_W3 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-02-02']
Cont_W4 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-02-09']
Cont_W5 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-02-16']
Cont_W6 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Control" & Spending_perweek_2grp$Week=='2015-02-23']

Test_W1 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-01-19']
Test_W2 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-01-26']
Test_W3 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-02-02']
Test_W4 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-02-09']
Test_W5 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-02-16']
Test_W6 = Spending_perweek_2grp$Spend[Spending_perweek_2grp$RedeemGroup=="Test" & Spending_perweek_2grp$Week=='2015-02-23']


ActualSpend   =(Test_W1*count_test_W1) +(Test_W2*count_test_W2)+(Test_W3*count_test_W3)+(Test_W4*count_test_W4)+(Test_W5*count_test_W5)+(Test_W6*count_test_W6)
ExpectedSpend =(Cont_W1*count_test_W1) +(Cont_W2*count_test_W2)+(Cont_W3*count_test_W3)+(Cont_W4*count_test_W4)+(Cont_W5*count_test_W5)+(Cont_W6*count_test_W6)

perIncrease = (ActualSpend-ExpectedSpend)*100/ExpectedSpend
redeemer =length(unique(All_users$User.ID[All_users$RedeemerGroup == "Redeemer"]))


Revenue = (ActualSpend-ExpectedSpend) * 0.1 
Investment = 10 * redeemer
ROI =  (Revenue-Investment)/Investment
ROI
###### DELETEDTest for Normality and difference in spending control vs test
######DELETED Spend  = Base + Coef X Coupon
#Spending_perweek <- setNames(aggregate(All_users$Spend ~ All_users$Week + All_users$RedeemerGroup, FUN=sum),c("Week", "RedeemGroup","Spend"))
#Spending_perweek$Spend <- ifelse(Spending_perweek$RedeemGroup =="Control",Spending_perweek$Spend/1001 ,ifelse(Spending_perweek$RedeemGroup =="Non-Redeemer",Spending_perweek$Spend/3192,Spending_perweek$Spend/804) )
#Spending_perweek$Period <- ifelse(Spending_perweek$Week < '2015-11-24', "Neutral", ifelse(Spending_perweek$Week <'2015-01-05', "PreLaunch", ifelse(Spending_perweek$Week <='2015-02-18',"Launch","PostLaunch")))
#Spending_perweek_4groups <-setNames(aggregate(All_users$Spend, by=list(Category=All_users$RedeemerGroup, Week=All_users$Week), FUN=mean) ,c( "RedeemGroup","Week","Spend"))                  
#Spending_perweek_4groups$Period <- ifelse(Spending_perweek_4groups$Week < '2014-11-24', "Neutral", ifelse(Spending_perweek_4groups$Week < '2015-01-05', "PreLaunch", ifelse(Spending_perweek_4groups$Week <='2015-02-18',"Launch","PostLaunch")))
#Spending_perweek_Pre <- filter(Spending_perweek,Spending_perweek$Week < '2014-11-24')
#Spending_perweek_Post <- filter(Spending_perweek,Spending_perweek$Week >= '2014-11-24')
#### Ploting
#ggplot(Spending_perweek , aes(x = Spending_perweek$Week, y = Spending_perweek$Spend, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg spending") + ylim(low=0,high=200) + ggtitle("Avg. Spending Per Week") + scale_x_continuous(breaks = round(seq(min(Spending_perweek$Week), max(Spending_perweek$Week), by = 1),1))
#ggplot(Spending_perweek_Pre , aes(x = Spending_perweek_Pre$Week, y = Spending_perweek_Pre$Spend, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg spending")+ ylim(low=0,high=200)
#ggplot(Spending_perweek_Post, aes(x = Spending_perweek_Post$Week, y = Spending_perweek_Post$Spend, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg spending")+ ylim(low=0,high=200)
###Pre Post Launch
#ggplot(All_users , aes(x = All_users$Week, y = All_users$Spend, fill= RedeemerGroup)) +  geom_bar(stat="identity", position=position_dodge())
#ggplot(Spending_perweek , aes(x = Spending_perweek$Period, y = Spending_perweek$Spend, fill=Spending_perweek$RedeemGroup )) +  geom_bar(stat="identity", position=position_dodge())
#ggplot(Spending_perweek_4groups , aes(x = Spending_perweek_4groups$Period, y = Spending_perweek_4groups$Spend, fill=Spending_perweek_4groups$RedeemGroup)) +  geom_bar(stat="identity", position=position_dodge())
#ggplot(Spending_perweek , aes(x = Spending_perweek$Week, y = Spending_perweek$Spend, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")
#ggplot(Spending_perweek_3groups , aes(x = Spending_perweek_3groups$Week, y = Spending_perweek_3groups$Spend, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg spending") + ggtitle("Avg. Spending Per Week")
#library(beeswarm)
#beeswarm(Spending_perweek$Spend~Spending_perweek$RedeemGroup, pch=16, cex=1.5, pwcol=as.numeric(sqrt(Spending_perweek$Spend)))
#boxplot(Spending_perweek$Spend~Spending_perweek$RedeemGroup)
#ggplot(Spending_perweek, aes(x=RedeemGroup, y=Spend, colour = RedeemGroup)) +geom_boxplot() + xlab("Weeks")+ ylab("Avg Cum spending") + ggtitle("Avg. Spending Per Week")
#library(dplyr)
#Spending_perweek_cum <- dplyr::mutate(group_by(Spending_perweek,RedeemGroup), spendSum=cumsum(Spend))
##### 
#ggplot(Spending_perweek_cum , aes(x = Spending_perweek_cum$Week, y = Spending_perweek_cum$spendSum, colour = RedeemGroup)) + geom_line() + xlab("Weeks")+ ylab("Avg Cum spending") + ggtitle("Avg. Spending Per Week")
#### Aggregating customers by Ids, summing spending and putting label
#aggregate(cbind(All_users$Spend) ~ All_users$Week + All_users$User.ID + All_users$RedeemerGroup, FUN=sum)
