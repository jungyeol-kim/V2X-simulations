# [Vehicular Messaging Simulation using Clustered Epidemiological Differential Equations]
## Required packages
list.of.packages <- c("deSolve")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)


## Setting parameters values
total.clusters<- 120      #the total number of clusters
sim.time<- 100     #end time of the time sequence
step.size<- 1         #step, increment of time
### e.g., if step.size=1 and sim.time=100, the result will be generated every 1 time unit, from 0 to 100 units.


## Mobility and Communication Network
setwd("/Users/jy/V2Vproject") # set working directory
mobility_network<-read.csv("mobility-network.csv", header=T, as.is=T)
mobility_network$lambda_from_to<-mobility_network$routing_prob*mobility_network$lambda

communication_network<-read.csv("communication-network.csv", header=T, as.is=T)
communication_network<-communication_network[order(communication_network$cluster_i,communication_network$cluster_j),]
rownames(communication_network) <- NULL # reset row names


## Defining a neighborhood of a cluster
mob.edge.out<-vector(mode='list',length=total.clusters) # outgoing edges from node i in the mobility network
mob.edge.out.rate<-vector(mode='list',length=total.clusters) # mobility parameter of outgoing edges from node i in the mobility network
mob.edge.in<-vector(mode='list',length=total.clusters) # incoming edges to node i in the mobility network
mob.edge.in.rate<-vector(mode='list',length=total.clusters) # mobility parameter of incoming edges from node i in the mobility network
comm.edge<-vector(mode='list',length=total.clusters) # undirected edges from or to node i in the communication network
comm.edge.rate<-vector(mode='list',length=total.clusters) # communication parameter of undirected edges from or to node i in the communication network
for (i in 1:total.clusters) {
  temp1<-mobility_network[mobility_network$from_clust==i,]
  temp2<-communication_network[communication_network$cluster_i==i,]
  mob.edge.out[[i]]<-temp1$to_clust
  mob.edge.out.rate[[i]]<-temp1$lambda_from_to
  comm.edge[[i]]<-temp2$cluster_j
  comm.edge.rate[[i]]<-temp2$beta_ij
  for(j in 1:length(mob.edge.out[[i]])){
    mob.edge.in[[mob.edge.out[[i]][j]]]<-c(mob.edge.in[[mob.edge.out[[i]][j]]],i)
    mob.edge.in.rate[[mob.edge.out[[i]][j]]]<-c(mob.edge.in.rate[[mob.edge.out[[i]][j]]],mob.edge.out.rate[[i]][j])
  }
}


## Generation of Differential Equations
mob.out.text<-vector(mode='character',length=total.clusters*2)
for(i in 1:total.clusters){
  mob.out.text.temp.1<-c();  mob.out.text.temp.2<-c();
  mob.out.text.temp.3<-c();  mob.out.text.temp.4<-c()
  for(j in 1:length(mob.edge.out[[i]])){
    mob.out.text.temp.1<-paste("- ",mob.edge.out.rate[[i]][j],"*y[",i,"]", sep="")
    mob.out.text.temp.2<-paste(mob.out.text.temp.2,mob.out.text.temp.1,sep=" ")
    mob.out.text.temp.3<-paste("- ",mob.edge.out.rate[[i]][j],"*y[",total.clusters+i,"]", sep="")
    mob.out.text.temp.4<-paste(mob.out.text.temp.4,mob.out.text.temp.3,sep=" ")
  }
  mob.out.text[i]<-mob.out.text.temp.2
  mob.out.text[total.clusters+i]<-mob.out.text.temp.4
}

comm.text<-vector(mode='character',length=total.clusters*2)
for(i in 1:total.clusters){
  comm.text.temp.1<-c();  comm.text.temp.2<-c();
  comm.text.temp.3<-c();  comm.text.temp.4<-c()
  for(j in 1:length(comm.edge[[i]])){
    comm.text.temp.1<-paste(comm.edge.rate[[i]][j],"*y[",comm.edge[[i]][j],"]*y[",total.clusters+i,"]", sep="")
    comm.text.temp.2<-paste(comm.text.temp.2," + ",comm.text.temp.1,sep="")
    comm.text.temp.3<-paste(comm.edge.rate[[i]][j],"*y[",comm.edge[[i]][j],"]*y[",total.clusters+i,"]", sep="")
    comm.text.temp.4<-paste(comm.text.temp.4," - ",comm.text.temp.3,sep="")
  }
  comm.text[i]<-comm.text.temp.2
  comm.text[total.clusters+i]<-comm.text.temp.4
}

mob.in.text<-vector(mode='character',length=total.clusters*2)
for(i in 1:total.clusters){
  mob.in.text.temp.1<-c();  mob.in.text.temp.2<-c();
  mob.in.text.temp.3<-c();  mob.in.text.temp.4<-c()
  for(j in 1:length(mob.edge.in[[i]])){
    mob.in.text.temp.1<-paste(mob.edge.in.rate[[i]][j],"*y[",mob.edge.in[[i]][j],"]", sep="")
    mob.in.text.temp.2<-paste(mob.in.text.temp.2," + ",mob.in.text.temp.1,sep=" ")
    mob.in.text.temp.3<-paste(mob.edge.in.rate[[i]][j],"*y[",total.clusters+mob.edge.in[[i]][j],"]", sep="")
    mob.in.text.temp.4<-paste(mob.in.text.temp.4," + ",mob.in.text.temp.3,sep=" ")
  }
  mob.in.text[i]<-mob.in.text.temp.2
  mob.in.text[total.clusters+i]<-mob.in.text.temp.4
}

dy<-vector(mode='character',length=total.clusters*2)
for (i in 1:(total.clusters*2)) {
   dy[i]<-paste("dy",i," <- ",mob.out.text[i],mob.in.text[i],comm.text[i],sep="")
} 

dy_name<-c()
for (i in 1:(total.clusters*2)) {
  if(i==1){dy_name<-paste(dy_name,"list(c(dy1",sep="")}
  else if(i==(total.clusters*2)){dy_name<-paste(dy_name,",dy",total.clusters*2,"))}",sep="")}
  else{dy_name<-paste(dy_name,",",paste("dy",i,sep=""),sep="")}  
}
set.diff.eqn<-c("f <- function(t, y, parms) {",dy,dy_name)
write(set.diff.eqn, file = "set_diff_eqn.R")
source("set_diff_eqn.R")


## Solving Differential Equations
initial_condition<-read.csv("initial-condition.csv", header=T, as.is=T) # import preset initial condition
yini<-c(initial_condition$I_ini,initial_condition$S_ini) # initial condition: 2J-dimensional vector

times <- seq(from = 0, to = sim.time, by = step.size) # output wanted at these time intervals
out <- ode(times = times, y=yini, func = f, parms = NULL) # numerically solve the set of differential equations 
solution<-out[,-1]
rownames(solution)<-times

frac.inf.clust<-solution[,1:total.clusters] # fraction of informed vehicles per cluster
frac.non.inf.clust<-solution[,(1+total.clusters):(2*total.clusters)]; colnames(frac.non.inf.clust)<-1:total.clusters # fraction of non-informed vehicles per cluster

write.table(frac.inf.clust, file = "fraction_of_informed_vehicles_per_cluster.csv",row.names=TRUE,col.names=TRUE, sep=",") # export a matrix to a file.
write.table(frac.non.inf.clust, file = "fraction_of_non_informed_vehicles_per_cluster.csv",row.names=TRUE,col.names=TRUE, sep=",") # export a matrix to a file.

write.table(frac.inf.clust, file = "fraction_of_informed_vehicles_per_cluster.csv",row.names=TRUE,col.names=TRUE, sep=",") # export a matrix to a file.
write.table(frac.non.inf.clust, file = "fraction_of_non_informed_vehicles_per_cluster.csv",row.names=TRUE,col.names=TRUE, sep=",") # export a matrix to a file.


## Generating Figures
frac.inf.veh<-rowSums(solution[,1:total.clusters]) #fraction of overall vehicles that are informed over time.
plot(times,frac.inf.veh, xlab="Time", ylab="Fraction of informed vehicles")

cluster.specific<-10 # determine the specific cluster of interest.
frac.inf.veh.clust<-frac.inf.clust[,cluster.specific] # fraction of vehicles that are informed over time in the particular cluster.
frac.non.inf.veh.clust<-frac.non.inf.clust[,cluster.specific] # fraction of vehicles that are informed over time in the particular cluster.

plot(times, frac.inf.veh.clust, xlab="Time", col="black", ylab=paste("Fraction of (non)informed vehicles in cluster ",cluster.specific,sep = ""))
par(new=T)
plot(times, frac.non.inf.veh.clust, xlab='', ylab='', col="red", axes=F)
par(new=F)
legend(0, 0.0025, legend=c("Infomred","Non-infomred"), pch = c(1, 1), col=c("black","red"))