# prepare packages
library(metafor)
library(meta)

# import data

#library(readxl)
#data <- read_excel("YourPath/data.xlsx")

# calculate effect sizes
dataset <-escalc(measure="ZCOR", ri=r, ni=N, data=data)

# overall associations
overall <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                  tdist = TRUE, data =  dataset)

summary(overall, digits=3)

# within-study variance 
modelnovar2 <- rma.mv(yi, vi, random = list(~ 1 | EID, 
                                            ~ 1 | SID), 
                      sigma2=c(0,NA), tdist=TRUE, data=dataset)
anova(overall,modelnovar2)


# between-study variance 
modelnovar3 <- rma.mv(yi, vi, random = list(~ 1 | EID, 
                                            ~ 1 | SID), sigma2=c(NA,0), 
                      tdist=TRUE, 
                      data=dataset)
anova(overall,modelnovar3)


# Determining how the total variance is distributed over the 
# three levels of the meta-analytic model;
# Print the results in percentages on screen.
n <- length(dataset$vi)
list.inverse.variances <- 1 / (dataset$vi)
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2
list.inverse.variances.square <- 1 / (dataset$vi^2)
sum.inverse.variances.square <-
  sum(list.inverse.variances.square)
numerator <- (n - 1) * sum.inverse.variances
denominator <- squared.sum.inverse.variances - sum.inverse.variances.square

estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (overall$sigma2[1]
                                         + overall$sigma2[2] + estimated.sampling.variance)
I2_2 <- (overall$sigma2[1]) / (overall$sigma2[1]
                               + overall$sigma2[2] + estimated.sampling.variance)
I2_3 <- (overall$sigma2[2]) / (overall$sigma2[1]
                               + overall$sigma2[2] + estimated.sampling.variance)
amountvariancelevel1 <- I2_1 * 100
amountvariancelevel2 <- I2_2 * 100
amountvariancelevel3 <- I2_3 * 100
amountvariancelevel1
amountvariancelevel2
amountvariancelevel3


# moderator analysis

## measurement time (designco = concurrently, designbefore = aggression before 
##                   ToM, designafter = aggression after ToM)
designco <- rma.mv(yi, vi, mods = ~ designbefore + designafter, random = 
                     list(~ 1 | EID, ~ 1 | SID), 
                   tdist=TRUE, data=dataset)
summary(designco, digits=3)

designbefore <- rma.mv(yi, vi, mods = ~ designco + designafter, random = 
                         list(~ 1 | EID, ~ 1 | SID), 
                       tdist=TRUE, data=dataset)
summary(designbefore, digits=3)

designafter <- rma.mv(yi, vi, mods = ~ designco + designbefore, random = 
                        list(~ 1 | EID, ~ 1 | SID), 
                      tdist=TRUE, data=dataset)
summary(designafter, digits=3)

# bully or not
bully <- rma.mv(yi, vi, mods = ~ nonbully, random = 
                      list(~ 1 | EID, ~ 1 | SID), 
                    tdist=TRUE, data=dataset)
summary(bully, digits=3)

nonbully <- rma.mv(yi, vi, mods = ~ bully, random = 
                         list(~ 1 | EID, ~ 1 | SID), 
                       tdist=TRUE, data=dataset)
summary(nonbully, digits=3)

# form (aggphysical = physical aggression, aggrelational = relational aggression)
aggphysical <- rma.mv(yi, vi, mods = ~ aggrelational, random = 
                         list(~ 1 | EID, ~ 1 | SID), 
                       tdist=TRUE, data=dataset)
summary(aggphysical, digits=3)

aggrelational <- rma.mv(yi, vi, mods = ~ aggphysical, random = 
                           list(~ 1 | EID, ~ 1 | SID), 
                         tdist=TRUE, data=dataset)
summary(aggrelational, digits=3)

# function (aggreactive = reactive aggression, aggproactive = proactive aggression)
aggreactive <- rma.mv(yi, vi, mods = ~ aggproactive, random = 
                         list(~ 1 | EID, ~ 1 | SID), 
                       tdist=TRUE, data=dataset)
summary(aggreactive, digits=3)

aggproactive <- rma.mv(yi, vi, mods = ~ aggreactive, random = 
                          list(~ 1 | EID, ~ 1 | SID), 
                        tdist=TRUE, data=dataset)
summary(aggproactive, digits=3)

## aggression measurement (aggself = self report, aggteacher = teacher report,
## aggparent = parent report, aggpeer = peer report)
aggself <- rma.mv(yi, vi, mods = ~ aggparent + aggteacher + aggpeer, random = 
                    list(~ 1 | EID, ~ 1 | SID), 
                  tdist=TRUE, data=dataset)
summary(aggself, digits=3)

aggparent <- rma.mv(yi, vi, mods = ~ aggself + aggteacher + aggpeer, random = 
                      list(~ 1 | EID, ~ 1 | SID), 
                    tdist=TRUE, data=dataset)
summary(aggparent, digits=3)

aggteacher <- rma.mv(yi, vi, mods = ~ aggparent + aggself + aggpeer, random = 
                       list(~ 1 | EID, ~ 1 | SID), 
                     tdist=TRUE, data=dataset)
summary(aggteacher, digits=3)

aggpeer <- rma.mv(yi, vi, mods = ~ aggparent + aggself + aggteacher, random = 
                       list(~ 1 | EID, ~ 1 | SID), 
                     tdist=TRUE, data=dataset)
summary(aggpeer, digits=3)

## ToM type (ToMAPT = affective ToM, ToMCPT = cognitive ToM, TOMACPT = combination
##           of affective and cognitive ToM)
ToMAPT <- rma.mv(yi, vi, mods = ~ ToMCPT + ToMACPT, random = 
                   list(~ 1 | EID, ~ 1 | SID), 
                 tdist=TRUE, data=dataset)
summary(ToMAPT, digits=3)

ToMCPT <- rma.mv(yi, vi, mods = ~ ToMAPT + ToMACPT, random = 
                   list(~ 1 | EID, ~ 1 | SID), 
                 tdist=TRUE, data=dataset)
summary(ToMCPT, digits=3)

ToMACPT <- rma.mv(yi, vi, mods = ~ ToMCPT + ToMAPT, random = 
                    list(~ 1 | EID, ~ 1 | SID), 
                  tdist=TRUE, data=dataset)
summary(ToMACPT, digits=3)

## age group (agehigh = 7 years old and above, agelow = 6 years old and below,
##            age middle = 2 to 15 years old)
agelow <- rma.mv(yi, vi, mods = ~ agemiddle+agehigh, random = 
                   list(~ 1 | EID, ~ 1 | SID), 
                 tdist=TRUE, data=dataset)
summary(agelow, digits=3)

agehigh <- rma.mv(yi, vi, mods = ~ agelow + agemiddle, random = 
                    list(~ 1 | EID, ~ 1 | SID), 
                  tdist=TRUE, data=dataset)
summary(agehigh, digits=3)

agemiddle <- rma.mv(yi, vi, mods = ~ agelow+agehigh, random = 
                      list(~ 1 | EID, ~ 1 | SID), 
                    tdist=TRUE, data=dataset)
summary(agemiddle, digits=3)


## individualism of culture (c_individualism = centralized index of individualism of country)
individualism <- rma.mv(yi, vi, mods = ~ c_individualism, random = list(~ 1 | EID, ~ 1 | SID), 
                        tdist=TRUE, data=dataset)
summary(individualism, digits=3)

## moderation of publication 
published <- rma.mv(yi, vi, mods = ~ unpublish, random = 
                      list(~ 1 | EID, ~ 1 | SID), 
                    tdist=TRUE, data=dataset)
summary(published, digits=3)

notpublished <- rma.mv(yi, vi, mods = ~ publish, random = 
                         list(~ 1 | EID, ~ 1 | SID), 
                       tdist=TRUE, data=dataset)
summary(notpublished, digits=3)


## Egger's test and Trim and fill for multilevel meta-analysis
SE <-sqrt(overall$vi)
Egger <- rma.mv(yi, vi, mods=SE, random = list (~1| EID, ~1 | SID), data=dataset)
print(Egger)
Funnel <- rma.mv(yi, vi, mods=N, random= list (~1| EID, ~1 | SID), data= dataset)
print(Funnel)

h <- rma.mv(yi, vi, random= list (~1| EID, ~1 | SID), data= dataset)
pooled_d<-h$b[1]
d<-dataset$yi
R0_func<-function(d, pooled_d){ 
  d_difference1=0
  rank_positive=0
  final_rank=0
  d_difference<- d-pooled_d
  for(i in 1:length(d_difference)){
    d_difference1[i]<- if ( d_difference[i]<0) d_difference[i]*-1 
    else d_difference[i] }
  rank_positive<-rank(d_difference1)
  for(i in 1:length(d_difference)){
    final_rank[i]<- if (d_difference[i]<0) rank_positive[i]*-1
    else rank_positive[i]}
  O=length(d_difference)-(min(final_rank)*-1)
  R0_ps= O-1
  R0=if (R0_ps<0) 0 else  R0_ps
  print(R0)
}
L0_func<-function(d, pooled_d){ 
  d_difference1=0
  rank_positive=0
  d_difference<- d-pooled_d
  d_difference1=0
  only_positive=0
  final_rank=0
  N_ES<-length(d)
  for(i in 1:length(d_difference)){
    d_difference1[i]<- if ( d_difference[i]<0) d_difference[i]*-1 
    else d_difference[i] }
  rank_positive<-rank(d_difference1)
  for(i in 1:length(d_difference)){
    final_rank[i]<- if (d_difference[i]<0) rank_positive[i]*-1
    else rank_positive[i]
    only_positive[i]<- if (final_rank[i]<0) 0 
    else final_rank[i]}
  t=sum(only_positive)
  L0_ps=(4*t-N_ES*(N_ES+1))/(2*N_ES-1)
  L0= if (L0_ps<0) 0 else  L0_ps
  print(L0)
}

R0_func(d, pooled_d)
L0_func(d, pooled_d)


# funnel plot
## Aggregate data for funnel plot 
dataset.agg <- dataset[(c("SID", "N", "yi", "vi"))]
dataset.agg <- aggregate(cbind(yi, vi, N) ~ SID, data = dataset.agg, mean, na.rm = TRUE)
overall.agg <- rma(yi, vi, data =  dataset.agg)
funnel(overall.agg)

# forest plot
m1 <- metacor(data = data, r, N, studlab=paste(Study), sm = "COR",comb.fixed = FALSE)
print(summary(m1,digits=2))
jpeg("picture.jpeg",height=2700,width=800)
forest(m1,hetstat=FALSE, leftlabs=c("Study","Sample size"))
dev.off()

# heterogeneity for each effect size
datapublish <- subset(dataset, publish == 1, select = c(EID,SID,yi,vi))
qpublish <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                  tdist = TRUE, data =  datapublish)

summary(qpublish, digits=3)

dataunpublish <- subset(dataset, unpublish == 1, select = c(EID,SID,yi,vi))
qunpublish <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                   tdist = TRUE, data =  dataunpublish)

summary(qunpublish, digits=3)

dataagelow <- subset(dataset, agelow == 1, select = c(EID,SID,yi,vi))
qagelow <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                  tdist = TRUE, data =  dataagelow)

summary(qagelow, digits=3)

dataagehigh <- subset(dataset, agehigh == 1, select = c(EID,SID,yi,vi))
qagehigh <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                   tdist = TRUE, data =  dataagehigh)

summary(qagehigh, digits=3)

dataagemiddle <- subset(dataset, agemiddle == 1, select = c(EID,SID,yi,vi))
qagemiddle <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                     tdist = TRUE, data =  dataagemiddle)

summary(qagemiddle, digits=3)

datadesignco <- subset(dataset, designco == 1, select = c(EID,SID,yi,vi))
qdesignco <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                    tdist = TRUE, data =  datadesignco)

summary(qdesignco, digits=3)

datadesignbefore <- subset(dataset, designbefore == 1, select = c(EID,SID,yi,vi))
qdesignbefore <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                        tdist = TRUE, data =  datadesignbefore)

summary(qdesignbefore, digits=3)

datadesignafter <- subset(dataset, designafter == 1, select = c(EID,SID,yi,vi))
qdesignafter <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                       tdist = TRUE, data =  datadesignafter)

summary(qdesignafter, digits=3)

dataToMAPT <- subset(dataset, ToMAPT == 1, select = c(EID,SID,yi,vi))
qToMAPT <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                  tdist = TRUE, data =  dataToMAPT)

summary(qToMAPT, digits=3)

dataToMCPT <- subset(dataset, ToMCPT == 1, select = c(EID,SID,yi,vi))
qToMCPT <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                  tdist = TRUE, data =  dataToMCPT)

summary(qToMCPT, digits=3)


dataToMACPT <- subset(dataset, ToMACPT == 1, select = c(EID,SID,yi,vi))
qToMACPT <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                   tdist = TRUE, data =  dataToMACPT)

summary(qToMACPT, digits=3)

databully <- subset(dataset, bully == 1, select = c(EID,SID,yi,vi))
qbully <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                 tdist = TRUE, data =  databully)

summary(qbully, digits=3)

datanonbully <- subset(dataset, nonbully == 1, select = c(EID,SID,yi,vi))
qnonbully <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                    tdist = TRUE, data =  datanonbully)

summary(qnonbully, digits=3)

dataaggphysical <- subset(dataset, aggphysical == 1, select = c(EID,SID,yi,vi))
qaggphysical <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                       tdist = TRUE, data =  dataaggphysical)

summary(qaggphysical, digits=3)

dataaggrelational <- subset(dataset, aggrelational == 1, select = c(EID,SID,yi,vi))
qaggrelational <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                         tdist = TRUE, data =  dataaggrelational)

summary(qaggrelational, digits=3)

dataaggreactive <- subset(dataset, aggreactive == 1, select = c(EID,SID,yi,vi))
qaggreactive <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                       tdist = TRUE, data =  dataaggreactive)

summary(qaggreactive, digits=3)

dataaggproactive <- subset(dataset, aggproactive == 1, select = c(EID,SID,yi,vi))
qaggproactive <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                        tdist = TRUE, data =  dataaggproactive)

summary(qaggproactive, digits=3)

dataaggself <- subset(dataset, aggself == 1, select = c(EID,SID,yi,vi))
qaggself <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                   tdist = TRUE, data =  dataaggself)

summary(qaggself, digits=3)

dataaggparent <- subset(dataset, aggparent == 1, select = c(EID,SID,yi,vi))
qaggparent <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                     tdist = TRUE, data =  dataaggparent)

summary(qaggparent, digits=3)

dataaggteacher <- subset(dataset, aggteacher == 1, select = c(EID,SID,yi,vi))
qaggteacher <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                      tdist = TRUE, data =  dataaggteacher)

summary(qaggteacher, digits=3)

dataaggpeer <- subset(dataset, aggpeer == 1, select = c(EID,SID,yi,vi))
qaggpeer <- rma.mv(yi, vi, random = list (~1| EID, ~1 | SID), 
                   tdist = TRUE, data =  dataaggpeer)

summary(qaggpeer, digits=3)
