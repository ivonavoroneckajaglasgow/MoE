setwd("C:/Users/Ivona Voroneckaja/Desktop/MoE/MoE")
#x<-seq(from=-5,to=5,length.out = 100)
#y<-sin(x)
#y<-y+rnorm(n = 100,mean = 0,sd = 0.25)

data<-read.csv("data_sin.csv")
y<-data[,1]
X<-data[,2:3]
x<-data[,3]

m1<-lm(y[1:33]~x[1:33])
m2<-lm(y[34:66]~x[34:66])
m3<-lm(y[67:100]~x[67:100])

coef(m1)
coef(m2)
coef(m3)

z_before<-read.csv("z_before.csv",header=FALSE)
z_before<-z_before[,1]

#z_before<-c(rep(1,times=33),rep(2,times=33),rep(3,times=34))
# z_before[10:15]<-2
# z_before[45:50]<-3
# z_before[80:85]<-1
# write.csv(file="z_before.csv",z_before,row.names = FALSE)
# write.table(z_before,file="z_before.csv",sep=",",col.names=FALSE,row.names = FALSE)

par(mfrow=c(1,1))
#plot(y~x,pch=col=c(rep(1,times=33),rep(2,times=33),rep(3,times=34)),main="Estimates by R")
plot(y~x,pch=20,col=z_before,main="Estimates by R")
abline(m1)
abline(m2,col=2)
abline(m3,col=3)


###Results after MCMC###
coefs<-read.csv("coefs.csv",header = FALSE)
logsigma_sqs<-as.numeric(read.csv("sigmas.csv",header = FALSE))
betas <-as.numeric(unlist(coefs[,1:3]))
gammas<-as.numeric(unlist(coefs[,4:5]))
sigmas<-exp(logsigma_sqs)
z_after<-read.csv("z_afterMCMC.csv",header=FALSE)
z_after<-z_after[,1]
P<-read.csv("P.csv",header=FALSE)
x_new<-read.csv("x_new.csv",header=FALSE)

par(mfrow=c(1,1))
plot(y~x,col=z_after,pch=20,main="After MCMC")
lines(x, betas[1]+betas[2]*x)
lines(x, betas[1]+betas[2]*x+2*sqrt(sigmas[1]),lty=2)
lines(x, betas[1]+betas[2]*x-2*sqrt(sigmas[1]),lty=2)
lines(x, betas[3]+betas[4]*x,col=2)
lines(x, betas[3]+betas[4]*x+2*sqrt(sigmas[2]),lty=2,col=2)
lines(x, betas[3]+betas[4]*x-2*sqrt(sigmas[2]),lty=2,col=2)
lines(x, betas[5]+betas[6]*x,col=3)
lines(x, betas[5]+betas[6]*x+2*sqrt(sigmas[3]),lty=2,col=3)
lines(x, betas[5]+betas[6]*x-2*sqrt(sigmas[3]),lty=2,col=3)

par(mfrow=c(1,1))
plot(x_new[,2],P[nrow(P),],main="Predictions",ylab="Predicted value",xlab="x new")

par(mfrow=c(1,2))
plot(x,exp(gammas[1]+gammas[2]*x)/(1+exp(gammas[1]+gammas[2]*x)),type="l",main="Gate 1",ylab="pi")
plot(x,exp(gammas[3]+gammas[4]*x)/(1+exp(gammas[3]+gammas[4]*x)),type="l",main="Gate 2",ylab="pi")

###After Swap###
z_swap<-read.csv("z_swap.csv",header=FALSE)
z_swap<-z_swap[,1]
coefs_swap <-read.csv("coefs_swap.csv",header = FALSE)
betas_swap <-as.numeric(unlist(coefs[,1:3]))
gammas_swap<-as.numeric(unlist(coefs[,4:5]))
par(mfrow=c(1,1))
plot(x,y,pch=20,col=z_after)
plot(x,y,pch=20,col=z_swap)
lines(x, betas_swap[1]+betas_swap[2]*x,col=1)
lines(x, betas_swap[3]+betas_swap[4]*x,col=2)
lines(x, betas_swap[5]+betas_swap[6]*x,col=3)

plot(x,exp(gammas[1]+gammas[2]*x)/(1+exp(gammas[1]+gammas[2]*x)),type="l",main="Gate 2",ylab="pi")
plot(x,exp(gammas[3]+gammas[4]*x)/(1+exp(gammas[3]+gammas[4]*x)),type="l",main="Gate 2",ylab="pi")



