setwd("C:/Users/Ivona Voroneckaja/Desktop/MoE/MoE")
data<-read.csv("data_sin_stand.csv",header = FALSE)
y<-data[,1]
X<-data[,2:3]
x<-data[,3]
z_before<-read.csv("z_before.csv",header=FALSE)
z_before<-z_before[,1]
plot(x,y,col=z_before,pch=20, main="Starting point")
abline(lm(y[z_before==1]~x[z_before==1]))
abline(lm(y[z_before==2]~x[z_before==2]),col=2)

###After 10 runs of MCMC
z_afterMCMC<-read.csv("z_afterMCMC.csv",header=FALSE)
z_afterMCMC<-z_afterMCMC[,1]

plot(x,y,col=z_afterMCMC,pch=20, main="After 10 MCMC Runs")

beta1<-c(-1.8061,-0.6532)
beta2<-c(1.0757,-0.4515)
lines(x,beta1[1]+beta1[2]*x)
lines(x,beta2[1]+beta2[2]*x,col=2)

X_sub<-X[z_afterMCMC==2,]
y_sub<-y[z_afterMCMC==2,]

###After Split

z_afterSplit<-read.csv("z_afterSplit.csv",header=FALSE)
z_afterSplit<-z_afterSplit[,1]
abline(v=1.66667,col="grey")
gamma<-c(15.0697,-8.87489)

lines(x,(exp(gamma[1]+gamma[2]*x)/(1+exp(gamma[1]+gamma[2]*x)))*3.2-1.5,lty=2,col="grey")

plot(x,y,col=z_afterSplit,pch=20, main="After Split")
abline(v=1.66667,col="grey")
lines(x,(exp(gamma[1]+gamma[2]*x)/(1+exp(gamma[1]+gamma[2]*x)))*3.2-1.5,lty=2,col="grey")
beta1<-c(-1.8061,-0.6532)
beta2<-c(2.8086,-0.9263)
beta3<-c(0.4112,0.5580)
lines(x,beta1[1]+beta1[2]*x)
lines(x,beta2[1]+beta2[2]*x,col=2)
lines(x,beta3[1]+beta3[2]*x,col=3)

###After another 100 runs of MCMC

z_afterMCMC2<-read.csv("z_afterMCMC2.csv",header=FALSE)
z_afterMCMC2<-z_afterMCMC2[,1]

plot(x,y,col=z_afterMCMC2,pch=20, main="After Another 100 MCMC Runs")
beta1<-c(-2.3084,-0.7966)
beta2<-c(2.4667,-0.8286)
beta3<-c(0.1605,0.9285)
lines(x,beta1[1]+beta1[2]*x)
lines(x,beta2[1]+beta2[2]*x,col=2)
lines(x,beta3[1]+beta3[2]*x,col=3)

###Record of acceptance of splits
beta1<-read.csv("beta1.csv",header=FALSE)
beta1<-t(beta1[,-(ncol(beta1))])
beta2<-read.csv("beta2.csv",header=FALSE)
beta2<-t(beta2[,-(ncol(beta2))])
beta3<-read.csv("beta3.csv",header=FALSE)
beta3<-t(beta3[,-(ncol(beta3))])
gamma1<-read.csv("gamma1.csv",header=FALSE)
gamma1<-t(gamma1[,-(ncol(gamma1))])
gamma2<-read.csv("gamma2.csv",header=FALSE)
gamma2<-t(gamma2[,-(ncol(gamma2))])
z<-read.csv("z.csv",header=FALSE)
z<-t(z[,-(ncol(z))])

beta1MCMC<-read.csv("beta1MCMC.csv",header=FALSE)
beta1MCMC<-t(beta1MCMC[,-(ncol(beta1MCMC))])
beta2MCMC<-read.csv("beta2MCMC.csv",header=FALSE)
beta2MCMC<-t(beta2MCMC[,-(ncol(beta2MCMC))])
beta3MCMC<-read.csv("beta3MCMC.csv",header=FALSE)
beta3MCMC<-t(beta3MCMC[,-(ncol(beta3MCMC))])
gamma1MCMC<-read.csv("gamma1MCMC.csv",header=FALSE)
gamma1MCMC<-t(gamma1MCMC[,-(ncol(gamma1MCMC))])
gamma2MCMC<-read.csv("gamma2MCMC.csv",header=FALSE)
gamma2MCMC<-t(gamma2MCMC[,-(ncol(gamma2MCMC))])
zMCMC<-read.csv("zMCMC.csv",header=FALSE)
zMCMC<-t(zMCMC[,-(ncol(zMCMC))])

acc<-read.csv("split_record.csv",header=FALSE)
acc<-acc[,1]
mean(acc)

set.seed(38)
choice<-sample(nrow(beta1),1)
for(i in 1:nrow(beta1)){
  if(i%%20==0){
  par(mfrow=c(1,2))
  choice<-i
  plot(x,y,col=z[choice,],pch=20, main="An Example of Accepted Jump")
  lines(x,beta1[choice,1]+beta1[choice,2]*x)
  lines(x,beta2[choice,1]+beta2[choice,2]*x,col=2)
  lines(x,beta3[choice,1]+beta3[choice,2]*x,col=3)
  #lines(x,(exp(gamma1[choice,1]+gamma1[choice,2]*x)/(1+exp(gamma1[choice,1]+gamma1[choice,2]*x)))*3.2-1.5,lty=2,col="grey")
  #lines(x,(exp(gamma2[choice,1]+gamma2[choice,2]*x)/(1+exp(gamma2[choice,1]+gamma2[choice,2]*x)))*3.2-1.5,lty=2,col="grey")


  plot(x,y,col=zMCMC[choice,],pch=20, main="An Example of Accepted Jump Post MCMC")
  lines(x,beta1MCMC[choice,1]+beta1MCMC[choice,2]*x)
  lines(x,beta2MCMC[choice,1]+beta2MCMC[choice,2]*x,col=2)
  lines(x,beta3MCMC[choice,1]+beta3MCMC[choice,2]*x,col=3)
  #lines(x,(exp(gamma1MCMC[choice,1]+gamma1MCMC[choice,2]*x)/(1+exp(gamma1MCMC[choice,1]+gamma1MCMC[choice,2]*x)))*3.2-1.5,lty=2,col="grey")
  #lines(x,(exp(gamma2MCMC[choice,1]+gamma2MCMC[choice,2]*x)/(1+exp(gamma2MCMC[choice,1]+gamma2MCMC[choice,2]*x)))*3.2-1.5,lty=2,col="grey")
  }
}

###Record of acceptance of merges

mer<-read.csv("merge_record.csv",header=FALSE)
mer<-mer[,1]
mean(mer)

#####TESTING BOTH DIRECTIONS#######

data<-read.csv("data_sin_stand.csv",header = FALSE)
y<-data[,1]
X<-data[,2:3]
x<-data[,3]
z_before<-read.csv("z_before.csv",header=FALSE)
z_before<-z_before[,1]
plot(x,y,col=z_before,pch=20, main="Starting point")
abline(lm(y[z_before==1]~x[z_before==1]))
abline(lm(y[z_before==2]~x[z_before==2]),col=2)

###After 10 runs of MCMC
z_afterMCMC<-read.csv("z_afterMCMC_both.csv",header=FALSE)
z_afterMCMC<-z_afterMCMC[,1]

plot(x,y,col=z_afterMCMC,pch=20, main="After 10 MCMC Runs")

###Propose a split

z_afterSplit<-read.csv("z_afterSplit_both.csv",header=FALSE)
z_afterSplit<-z_afterSplit[,1]
plot(x,y,col=z_afterSplit,pch=20, main="After Split")

###MCMC post split

z_afterMCMC<-read.csv("z_afterMCMC_both2.csv",header=FALSE)
z_afterMCMC<-z_afterMCMC[,1]

plot(x,y,col=z_afterMCMC,pch=20, main="After 100 MCMC Runs post split")

###AFTER RJ MCMC##
z_before<-read.csv("z_before.csv",header=FALSE) #case 1
z_before<-read.csv("z_for_merge.csv",header=FALSE) #case 2
z_before<-read.csv("z_4exp_merge.csv",header=FALSE) #case 3
z_before<-z_before[,1]
plot(x,y,col=z_before,pch=20, main="Starting point")

z_afterMCMC<-read.csv("z_MCMCRJ.csv",header=FALSE)
z_afterMCMC<-z_afterMCMC[,1]
plot(x,y,col=z_afterMCMC,pch=20, main="After RJ MCMC")

