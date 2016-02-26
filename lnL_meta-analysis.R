rm(list=ls(all=TRUE))


setwd("Y:/Kevin/brouillons/meta-analysis vaccine")

# library(Hmisc)
# library(ggplot2)
library(nlme) #for BIC
library(lme4)
trans = function(x) exp(x) / (1 + exp(x))



tab = read.csv("meta-analyse vaccine.csv", h=T)
tab

# if we exclude the 0.75% efficacy study (3rd line)
tab= tab[-3,]


tab$efficacy=  tab$responder/ tab$sample
tab$conf = binconf(tab$responder, tab$sample)
  


tab$lnL = lgamma(tab$sample+1) - lgamma(tab$responder+1) - lgamma(tab$sample-tab$responder+1) + 
  tab$responder*log(tab$efficacy) + 
  ifelse(tab$efficacy == 1, 0, (tab$sample-tab$responder)*log(1-tab$efficacy))

sum(tab$lnL)

tot_sample= sum(tab$sample)
tot_responder = sum(tab$responder)
tot_efficacy=  tot_responder/ tot_sample
binconf(tot_responder, tot_sample)

tab$one_param_lnL = lgamma(tab$sample+1) - lgamma(tab$responder+1) - lgamma(tab$sample-tab$responder+1) + 
  tab$responder*log(tot_efficacy) + 
  ifelse(tab$efficacy == 1, 0, (tab$sample-tab$responder)*log(1-tot_efficacy))





-2*(sum(tab$one_param_lnL) -sum(tab$lnL) )

1 - pchisq(-2*(sum(tab$one_param_lnL) -sum(tab$lnL) ) , df= 15)


binconf(tot_responder, tot_sample)


tab$conf = binconf(tab$responder, tab$sample)


####




###################################
### Plotting a forrest plot

tab$conf[,1]
d = data.frame(x=paste(tab$Study, tab$Vaccine, sep=", "), y=tab$conf[,1], ylo = tab$conf[,2], yhi = tab$conf[,3])
d
                 
credplot.gg <- function(d){
  # d is a data frame with 4 columns
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+
    coord_flip() + geom_hline(aes(x=0.97), lty=2)+ xlab('Variable')
  return(p)
}

credplot.gg(d)

#reorder 
d$x <- factor(d$x, levels=d$x[order(d$y)])

plot_theme<-theme(axis.text = element_text(size=11),
                  axis.title.x = element_text(vjust=-0.5, size=20),                      
                  axis.title.y = element_text(vjust=0.9, size=20)   )
p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) + geom_pointrange(size=1.01)+
  coord_flip() + geom_hline(aes(yintercept=0.97), lty=2, lwd=1.01)  + geom_hline(aes(yintercept=1), lty=1, lwd=1, col="red") +
  scale_y_continuous(limits=c(0.72, 1.01)) 
p + plot_theme + ylab("Vaccine efficacy") + xlab("Study") 

png("forest_plot_meta.png", width=10,height=7,res=120,units="in")
p + plot_theme + ylab("Vaccine efficacy") + xlab("Study") 
dev.off()




###################################
### meta-regression

tab = read.csv("meta-analyse vaccine.csv", h=T)

# if we exclude the 0.75% efficacy study (3rd line)
tab= tab[-3,]
tab

tab$Study=paste(tab$Study, tab$Vaccine, sep=", ")
tab=tab[,-3]

sum(tab$sample)
sum(tab$responder)
sum(tab$responder)/sum(tab$sample)


# creating individual table
resp=Study=wider_study=type=vaccine_type=NULL
datind=dat.tmp=data.frame()

for (i in 1: nrow(tab)){
  study = rep(tab$Study[i], tab$sample[i])
  wider_study = rep(tab$wider_study[i], tab$sample[i])
  type = rep(tab$type[i], tab$sample[i])
  vaccine_type = rep(tab$vaccine_type[i], tab$sample[i])
  resp = c(rep(1, tab$responder[i]), rep(0, (tab$sample[i]- tab$responder[i])))
  
  dat.tmp = data.frame(study,wider_study,type,vaccine_type,resp)
  datind = rbind(datind,dat.tmp)
}



datind$resp= as.factor(datind$resp)
datind$wider_study= as.factor(datind$wider_study)
summary(datind)
tab


?glmer

fglm =  glm(resp ~  study, data=datind, family="binomial" )
summary(fglm)





#######################################
##################################
# excluding 100% studies
##################################
#######################################
dim(tab)
tab2=tab[-c(3,13,15),]

# creating individual table
resp=Study=wider_study=type=vaccine_type=NULL
datind2=dat.tmp=data.frame()

for (i in 1: nrow(tab2)){
  study = rep(tab2$Study[i], tab2$sample[i])
  wider_study = rep(tab2$wider_study[i], tab2$sample[i])
  type = rep(tab2$type[i], tab2$sample[i])
  vaccine_type = rep(tab2$vaccine_type[i], tab2$sample[i])
  resp = c(rep(1, tab2$responder[i]), rep(0, (tab2$sample[i]- tab2$responder[i])))
  
  dat.tmp = data.frame(study,wider_study,type,vaccine_type,resp)
  datind2 = rbind(datind2,dat.tmp)
}

datind2$resp= as.factor(datind2$resp)
datind2$wider_study= as.factor(datind2$wider_study)
summary(datind2)
tab


glm1bis = glmer(resp ~  1+ (1|study), data=datind2, family="binomial" )
summary(glm1bis)
vavar = sqrt( 0.3262^2 + 0.909)
trans(3.5584)
trans(3.5584 - 1.96*(vavar))
trans(3.5584 + 1.96*(vavar))


#######################################
##################################
# trying to add 1 neg in the 100% studies
##################################
#######################################
#dim(tab)
#tab3=tab[-3,]

# creating individual table
resp=Study=wider_study=type=vaccine_type=NULL
datind3=dat.tmp=data.frame()

for (i in 1: nrow(tab3)){
  
  if(i ==12 | i==14) {
    study = rep(tab3$Study[i], tab3$sample[i]+1)
    wider_study = rep(tab3$wider_study[i], tab3$sample[i]+1)
    type = rep(tab3$type[i], tab3$sample[i]+1)
    vaccine_type = rep(tab3$vaccine_type[i], tab3$sample[i]+1)
    resp = c(rep(1, tab3$responder[i]),0)
  }
  else {  
    resp = c(rep(1, tab3$responder[i]), rep(0, (tab3$sample[i]- tab3$responder[i]))) 
    study = rep(tab3$Study[i], tab3$sample[i])
    wider_study = rep(tab3$wider_study[i], tab3$sample[i])
   type = rep(tab3$type[i], tab3$sample[i])
   vaccine_type = rep(tab3$vaccine_type[i], tab3$sample[i])}
  
  dat.tmp = data.frame(study,wider_study,type,vaccine_type,resp)
  datind3 = rbind(datind3,dat.tmp)
}

datind3$resp= as.factor(datind3$resp)
datind3$wider_study= as.factor(datind3$wider_study)
summary(datind3)
table(datind3$study)
tab


glm1ter = glmer(resp ~  1+ (1|study), data=datind3, family="binomial" )
summary(glm1ter)
vavar = sqrt( 0.3262^2 + 1.192)
trans(3.8271)
trans(3.8271 - 1.96*(vavar))
trans(3.8271 + 1.96*(vavar))

ranef(glm1ter)
var(ranef(glm1ter)$study$"(Intercept)")



###########################################
# trying with covariates
glm2 = glmer(resp ~ type +  (1|study), data=datind, family="binomial" )
summary(glm2)
# on retombe sur 97%, OK
# comment interpreter le reste??


glm2bis = glmer(resp ~ vaccine_type +  (1|study), data=datind, family="binomial" )
summary(glm2bis)
meta_prop = exp(3.4742) / (1 + exp(3.4742))
# on retombe sur 97%, OK
# comment interpreter le reste??

glm3 = glmer(resp ~ type +  (1|wider_study) + (1|study:wider_study), data=datind, family="binomial" )
summary(glm3)
trans(3.5206)
trans(3.5206-1.96*0.5346)
trans(3.5206+1.96*0.5346)


glm4 = glmer(resp ~ 1 +  (1|wider_study) + (1|study:wider_study), data=datind, family="binomial" )
summary(glm4)
trans(3.9368)
trans(3.9368-1.96*0.4529)
trans(3.9368+1.96*0.4529)


##############################################
#### plotiing on a logit scale


glm1 = glmer(resp ~  1+ (1|study), data=datind, family="binomial" )
summary(glm1)
# I sum the variance of fixed effect and random effect
0.3782^2 + 1.613
# and take the sqrt to have the sd
sqrt(0.3782^2 + 1.613)

trans(3.9923)
trans(3.9923 - 1.96*(1.325))
trans(3.9923 + 1.96*(1.325))


tab$conf = binconf(tab$responder, tab$sample)
tab$conf_logit = log(tab$conf/(1-tab$conf))
tab$conf_logit = ifelse(tab$conf_logit== Inf | tab$conf_logit>10, 10, tab$conf_logit)


d = data.frame(x=paste(tab$Study, tab$Vaccine, sep=", "), y=tab$conf_logit[,1], ylo = tab$conf_logit[,2], yhi = tab$conf_logit[,3])
d[,1]=as.character(d[,1])

# when adding fake negative results in 1 study
#add = c("Meta-estimate", 3.8271, 1.5937,  6.06047 )

# without post-hoc adjustement
sqrt(0.3782^2 + 1.613)
3.9923 - 1.96*(1.325)
3.9923 + 1.96*(1.325)
add = c("Meta-estimate", 3.9923, 1.3953,  6.5893 )

trans(3.9923)
trans(3.9923 - 1.96*(1.325))
trans(3.9923 + 1.96*(1.325))

d=rbind(d,add)
d

d$y = as.numeric(d$y)
d$ylo = as.numeric(d$ylo)
d$yhi = as.numeric(d$yhi)
d$x = as.factor(d$x)

print(levels(d$x))
d$x <- factor(d$x, levels=d$x[order(d$y)])
print(levels(d$x))
d$x <- factor(d$x, levels(d$x)[c(9,1:8,10:16)] )



plot_theme<-theme(axis.text = element_text(size=11),
                  axis.title.x = element_text(vjust=-0.5, size=20),                      
                  axis.title.y = element_text(vjust=0.9, size=20)   )

p = ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) + geom_pointrange(size=0.8)+
  coord_flip() +
  scale_y_continuous(limits=c(0,10))
p



png("forest_plot_meta_logit_scale.png", width=10,height=7,res=120,units="in")
p + plot_theme + ylab("Vaccine efficacy") + xlab("Study") 
dev.off()



########
## Vaccine prior distribution

### results of the meta-analysis gave a estimated efficacy distributed on a
### Normal distribution with parameters mean=3.9923 and sd=1.325

y= rnorm(10000, mean=3.9923, sd=1.325)
hist(y)


x <- seq(-2, 10, length=100000)
y <- dnorm(x, mean=3.9923, sd=1.325)
plot(x, y)

y_logit = exp(y)/(1 + exp(y))
hist(y_logit)

x_logit = exp(x)/(1 + exp(x))
hist(x_logit)

png("prior_distribution_efficacy.png", width=7,height=7,res=120,units="in")
plot(x_logit, y_logit, type="l", lty=1, lwd=2, xlab="Density", ylab="", col="blue", xlim=c(0.5,1), main="Prior distribution" )
dev.off()

library("pracma")
?trapz
trapz(x_logit, y_logit)
trapz(x,y)


### try with logis function
logi_trans = function (x) 1/(1+exp(-x))
y= rnorm(100000, mean=3.9923, sd=1.325)
sd(y)

logis_y = logi_trans(y)
mean(logis_y)
sd(logis_y)
quantile(logis_y, probs=c(0.025, 0.5, 0.975))

y2 = rlogis(1000000, location=mean(logis_y), scale =sqrt(3)*sd(logis_y)/pi  )
logitt = function()
hist(y2)
mean(y2)
sd(y2)
quantile(y2, probs=c(0.025, 0.5, 0.975))



scale =sqrt(3)*sd(logis_y)/pi
m = median(logis_y)

pr = function(x) 1/(1+exp(-(x-m)/scale))
pr(0.5)
pr(1)
pr(0.98)


x= 0.98



x=-1000:1000
fr = function(x) 1/(1+exp(-x))
plot(x, fr(x))

y2 = rlogis(1000000, location=mean(logis_y), scale =sqrt(3)*sd(logis_y)/pi  )
hist(y2)


x = seq(0,1, by=1e-5)
y = dlogis(x, location=mean(logis_y) , scale =1 )
hist(y)
plot(x, y, type="l", lty=1, lwd=2, xlab="Density", ylab="", col="blue")
# 
# 
# 
# ##################################
# ## On ne se sert plus de ce qu il y a ci dessous
# 
# 
# 
# ### values of the beta distrib
# # mean
# trans(3.9923)
# trans(3.9923 - 1.96*(1.325))
# trans(3.9923 + 1.96*(1.325))
# 
# range=trans(3.9923 + 1.96*(1.325)) - trans(3.9923 - 1.96*(1.325))
# range / 1.96
# (range / 1.96)^2
# 
# mu= trans(3.9923)
# var = ( ( trans(3.9923 + 1.96*(1.325)) - trans(3.9923 - 1.96*(1.325)) ) / 1.96)^2
# 
# estBetaParams <- function(mu, var) {
#   alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
#   beta <- alpha * (1 / mu - 1)
#   return(params = list(alpha = alpha, beta = beta))
# }
# alpha = estBetaParams(mu, var)$alpha
# beta = estBetaParams(mu, var)$beta
# alpha
# beta
# 
# z=rbeta(10000,alpha, beta)
# mean(z)
# quantile(z, probs=c(0.025, 0.5, 0.975))
# 
# x <- seq(0.90, 0.999, length=10000)
# #beta =  (1/0.99 - 1)*alpha
# y <- dbeta(x, alpha, beta)
# plot(x, y, type="l", lty=1, lwd=2, xlab="Density", ylab="", col="blue", xlim=c(0.94,1), main="Beta distribution" )
# 
# 
# 
# alpha=605
# beta=9.3
# z=rbeta(10000,alpha, beta)
# mean(z)
# var(z)
# x <- seq(0.90, 0.999, length=10000)
# #beta =  (1/0.99 - 1)*alpha
# y <- dbeta(x, alpha, beta)
# plot(x, y, type="l", lty=1, lwd=2, xlab="Density", ylab="", col="blue", xlim=c(0.94,1), main="Beta distribution" )
# 
# 
# 
# #### with a logistic distrib
# m = 3.9923
# v = 1.325
# 
# logi_trans = function (x) 1/(1+exp(-x))
# x_val = logi_trans(rnorm(10000, mean=m, sd=v))
# sd(x_val)
# mean(x_val)
# quantile(x_val, probs=c(0.025, 0.5, 0.975))
# 
# 
