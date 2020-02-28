# date: 2020/2/28
# Example code for implementing the method proposed in Zhang and Zhang: A stable and more efficient doubly robust estimator.
# simulated data, simulation setting 1. 



################################################# functions###########################

# see manuscript for descritions of the data and correct and wrong working models. 
# this function implements the proposed estimators with different specifications of working models. 
# dataframe: data set name (Yobs: outcome; r: non-missing indicator; z1-z4, x1-x4: covariates)
# outcomeX=T or F: indicates whether the outcome model is correctly specified or not.
# propensity score model: indicates whether the propensity score model is correct or not.
# lambda: bandwidth

############## estimate ##################
estimate<-function(dataframe, outcomeX,propensityX,lambda){
n<-length(dataframe$r)

if (outcomeX==T) omodel<-lm(Yobs~z1+z2+z3+z4,data=dataframe)
else omodel<-lm(Yobs~x1+x2+x3+x4,data=dataframe)

predicty<-predict(omodel,dataframe)

if (propensityX==T)  pimodel<-glm(r~z1+z2+z3+z4, family=binomial,data=dataframe)
else pimodel<-glm(r~x1+x2+x3+x4, family=binomial,data=dataframe)

predictpi<-predict(pimodel,dataframe)
hatpi<-exp(predictpi)/(1+exp(predictpi))
spi<-(hatpi-mean(hatpi))/sd(hatpi)

dryx<-c(rep(0,n))
loc<-(dataframe$r==1)
for(i in 1:n){
kij<-exp(-(spi[i]-spi[loc])^2/(lambda^2))
dryx[i]=sum(kij*(dataframe$Y[loc]-predicty[loc]))/sum(kij)
}

drsmthwt=sum(dryx+predicty)/n
return(drsmthwt)
}


################# bootstrap for variance estimation ################
bootstrap<-function(dataframe,outcomeX,propensityX,lambda){
n<-length(dataframe$r)
bootest<-matrix(0,100,1)
id<-seq(1,n)
for(i in 1:100){
bootid<-sample(id,replace=T)
bootdata<-dataframe[bootid,]
names(bootdata)<-names(dataframe)
bootest[i]<-estimate(dataframe=bootdata,outcomeX=outcomeX,propensityX=propensityX,lambda)
}

se<-sqrt(apply(bootest, 2, var))
return(se)
}





##########################################################################
### generate data###

# simu: number of Monte Carlo replicates
# n: sample size

simu<-10  
n<-1000
simudata<-NULL
for (j in 1:simu) {

z1<-rnorm(n)
z2<-rnorm(n)
z3<-rnorm(n) 
z4<-rnorm(n)
  
y<-210+27.4*z1+13.7*z2+13.7*z3+13.7*z4+rnorm(n)

pi<-exp(-z1+0.5*z2-0.25*z3-0.1*z4)/(1+exp(-z1+0.5*z2-0.25*z3-0.1*z4))
r<-NULL
for (i in 1:n) r<-c(r,rbinom(1,1,pi[i]))

x1<-exp(z1/2);
x2<-z2/(1+exp(z1))+10;
x3<-(z1*z3/25+0.6); x3=x3*x3*x3;
x4<-(z2+z4+20)*(z2+z4+20);
  
Yobs<-as.numeric(r==1)*y
Yobs[r==0]<-NA

simudata<- rbind(simudata, cbind(Yobs, r, z1,z2,z3,z4, pi, x1,x2,x3,x4, j))
}

simudata<-data.frame(simudata)
names(simudata)<-c("Yobs", "r","z1","z2","z3","z4","pi","x1","x2","x3","x4","j")


estimatorTT<-matrix(0,simu,1)
seTT<-matrix(0,simu,1)

estimatorTW<-matrix(0,simu,1)
seTW<-matrix(0,simu,1)

estimatorWT<-matrix(0,simu,1)
seWT<-matrix(0,simu,1)

estimatorWW<-matrix(0,simu,1)
seWW<-matrix(0,simu,1)


for(j in 1:simu){
print(j)
simudataj<-simudata[simudata$j==j,]

n<-length(simudataj$r)
lambda=n^(-1/3)
estimatorTT[j,]<-estimate(dataframe=simudataj,outcomeX=T,propensityX=T,lambda=lambda)
seTT[j,]<-bootstrap(dataframe=simudataj,outcomeX=T,propensityX=T,lambda=lambda)

estimatorTW[j,]<-estimate(dataframe=simudataj,outcomeX=T,propensityX=F,lambda=lambda)
seTW[j,]<-bootstrap(dataframe=simudataj,outcomeX=T,propensityX=F,lambda=lambda)


estimatorWT[j,]<-estimate(dataframe=simudataj,outcomeX=F,propensityX=T,lambda=lambda)
seWT[j,]<-bootstrap(dataframe=simudataj,outcomeX=F,propensityX=T,lambda=lambda)

estimatorWW[j,]<-estimate(dataframe=simudataj,outcomeX=F,propensityX=F,lambda=lambda)
seWW[j,]<-bootstrap(dataframe=simudataj,outcomeX=F,propensityX=F,lambda=lambda)
}



#### summarize results
mean(estimatorTT)-210
sd(estimatorTT)
mean(seTT)

mean(estimatorTW)-210
sd(estimatorTW)
mean(seTW)

mean(estimatorWT)-210
sd(estimatorWT)
mean(seWT)

mean(estimatorWW)-210
sd(estimatorWW)
mean(seWW)

