# date: 2020/07/14
# Example code for implementing the method proposed in Zhang and Zhang: A stable and more efficient doubly robust estimator.
 



################################################# functions ###########################

# see manuscript and code below for descritions of the data 
# dataframe: data set name (Yobs: outcome; r: non-missing indicator; x1-x4: observed covariates)
# lambda: bandwidth


############## estimate ##################
estimate<-function(dataframe,lambda){
n<-length(dataframe$r)

    ### build  outcome and propensity score models; users need to modify this part accordingly ###
 omodel<-lm(Yobs~x1+x2+x3+x4,data=dataframe)

 pimodel<-glm(r~x1+x2+x3+x4, family=binomial,data=dataframe)
    ### end of buliding outcome and propensity score models ###



  ### omodel is the model for outcome, and pimodel is the model for propensity, dataframe$Y is outcome ###
predicty<-predict(omodel,dataframe)
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
bootstrap<-function(dataframe,lambda){
n<-length(dataframe$r)
bootest<-matrix(0,100,1)
id<-seq(1,n)
for(i in 1:100){
bootid<-sample(id,replace=T)
bootdata<-dataframe[bootid,]
names(bootdata)<-names(dataframe)
bootest[i]<-estimate(dataframe=bootdata,lambda)
}

se<-sqrt(apply(bootest, 2, var))
return(se)
}





##########################################################################
### generate data###

# n: sample size

n<-1000
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
simudata<-  cbind(Yobs, r, z1,z2,z3,z4, pi, x1,x2,x3,x4)
simudata<-data.frame(simudata)
names(simudata)<-c("Yobs", "r","z1","z2","z3","z4","pi","x1","x2","x3","x4")

################## end of generating data#############3







################ call the functions to estimate mu and use bootstrap for se############3
n<-length(simudataj$r)
lambda=n^(-1/3)
estimator<-estimate(dataframe=simudata,lambda=lambda)
se<-bootstrap(dataframe=simudata, lambda=lambda)







