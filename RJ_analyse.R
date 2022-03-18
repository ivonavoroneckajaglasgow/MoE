setwd("C:/Users/Ivona Voroneckaja/Desktop/MoE/MoE")

###Choose which case and subcase
scenarios<- c("sin","cube")
case     <- scenarios[2]
subcase  <- 2
last_iter<- 1

###Load data
if(case=="sin")  data<-read.csv("data_sin.csv",header = FALSE)
if(case=="cube") data<-read.csv("data_cube.csv",header = FALSE)
y<-data[,1]
X<-data[,2:3]
x<-data[,3]
plot(x,y,pch=20)

###Load initial allocations
if(case=="sin"){
  if(subcase==1) z_before<-read.csv("z_sin1.csv",header=FALSE)
  if(subcase==2) z_before<-read.csv("z_sin2.csv",header=FALSE)
  if(subcase==3) z_before<-read.csv("z_sin3.csv",header=FALSE)
}
if(case=="cube"){
  if(subcase==1) z_before<-read.csv("z_cube1.csv",header=FALSE)
  if(subcase==2) z_before<-read.csv("z_cube2.csv",header=FALSE)
  if(subcase==3) z_before<-read.csv("z_cube3.csv",header=FALSE)
}
z_before<-z_before[,1]
plot(x,y,col=z_before,pch=20, main="Starting point")
n_exp_before<-length(unique(z_before))
for(i in 1: n_exp_before) abline(lm(y[z_before==i]~x[z_before==i]),col=i)

###Load predictions and acceptance record
accept<- read.csv("accept_RJ.csv",header=FALSE)
predictions <-read.csv("predictions.csv", header=FALSE)
no_expt<-read.csv("no_expt.csv",header=FALSE)
###If interested in the last iteration
if(last_iter==1){
  z_last<-read.csv("z_last_iter.csv",header=FALSE)
  z_last<-z_last[,1]
  par(mfrow=c(1,1))
  plot(x,y,col=z_last,pch=20, main="Last Iteration")
  params<-read.csv("all_params.csv",header=FALSE)
  exp_unique<-unique(z_last)
  for(i in exp_unique) lines(x,params[1,i]+params[2,i]*x,col=i)
  lines(x,predictions[,ncol(predictions)],col="grey")
}

###Analyse acceptance
accept_split<- accept[accept[,1]==1,]
accept_merge<- accept[accept[,1]==2,]
mean(accept_split[,2])
mean(accept_merge[,2])
mean(accept[,2])
barplot(table(accept[,1])/nrow(accept))

###Plot predictions
plot(x,y,pch=20,main="Predictions")
mycol <- rgb(0, 0, 255, max = 255, alpha = 50, names = "blue50")
for(i in 1:ncol(predictions)) if(i%%100==0) lines(x,predictions[,i],pch=2,col=mycol)

###Analyse autocorrelations
autocors<- c() 
for(i in 1:nrow(predictions)){
  autocors<-c(autocors, cor(c(unlist(predictions[i,2:ncol(predictions)])),c(unlist(predictions[i,1:(ncol(predictions)-1)]))))
}
summary(abs(autocors))
rbPal <- colorRampPalette(c('blue','red'))
Col <- rbPal(10)[as.numeric(cut(autocors,breaks = 10))]
plot(x,y,pch = 20,col=Col,main="Autocorrelation Heat Map")

###Analyse # of experts distribution
no_expt_all  <-no_expt[,1]
no_expt_empty<-no_expt[,2]
m<-rbind(c(0,1,1,0),c(2,2,3,3))
layout(m)
barplot(table(no_expt_all)/length(no_expt_all),main="All experts")
barplot(table(no_expt_empty)/length(no_expt_empty),main="Empty experts")
barplot(table(no_expt_all-no_expt_empty)/length(no_expt_all-no_expt_empty),main="Full Experts")
mean(no_expt_all)
median(no_expt_all)
mean(no_expt_all-no_expt_empty)
median(no_expt_all-no_expt_empty)
###Analyse convergence of predictions
par(mfrow=c(2,2))
rand_point<-sample(1:100,4)
for(i in 1:length(rand_point)) plot(c(unlist(predictions[rand_point[i],])),type="l",ylab=paste0("Observation ",rand_point[i]))

     