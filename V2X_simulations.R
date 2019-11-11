
# 0. R packages --------------------------------------------------------------
# install Packages
install.packages("igraph")
install.packages("qgraph")
install.packages("plyr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("zoo")
install.packages("matrixStats")
install.packages("MASS")
install.packages("RColorBrewer")
install.packages("oce")
install.packages("scales")
install.packages("pdist")
install.packages("reshape2")
install.packages("reshape")
install.packages("GGally")
install.packages("akima")
install.packages("rgl")
install.packages("data.table")
install.packages("sqldf")
install.packages("matlib")
install.packages("deSolve")
install.packages("tidyverse")
install.packages("plotly")

# load the packages for use
library(GGally)
library(igraph)
library(qgraph)
library(plyr)
library(dplyr)
library(ggplot2) 
library(zoo)
library(MASS)
library(RColorBrewer)
# Source the colorRampAlpha file
source ('colorRampPaletteAlpha.R')
library(oce)
library(scales)
library(pdist)
library(reshape)
library(reshape2)
library(network)
library(akima)
library(rgl)
library(data.table)
library(sqldf)
library(matlib)
library(deSolve)
library(readr)
library(tidyverse)
library(plotly)


# 1. Transportation network --------------------------------------------------------------
### 1-1. Setting a size of grid topology --------------------------------------------------------------
system_size<-2000 # size of whole area (system_size x system_size)
block<-400  # size of block (block x block)
sim_time<-150  # simulation time

### 1-2. Creating directed mobility network and corresponding adjacency metrix --------------------------------------------------------------
temp1<-seq(as.integer(block/2),system_size,as.integer(block))
temp2<-seq(0,system_size,as.integer(block))
x_h<-rep(temp1,times=length(temp2))
y_h<-rep(temp2,each=length(temp1))
x_v<-rep(temp2,each=length(temp1))
y_v<-rep(temp1,times=length(temp2))
intersec1<-data.frame("x"=c(x_h,x_v), "y"=c(y_h,y_v), "from_x"=c(x_h-as.integer(block/2),x_v), "from_y"=c(y_h,y_v-as.integer(block/2)), "to_x"=c(x_h+as.integer(block/2),x_v), "to_y"=c(y_h,y_v+as.integer(block/2)))
intersec2<-data.frame("x"=c(x_h,x_v), "y"=c(y_h,y_v), "from_x"=c(x_h+as.integer(block/2),x_v), "from_y"=c(y_h,y_v+as.integer(block/2)), "to_x"=c(x_h-as.integer(block/2),x_v), "to_y"=c(y_h,y_v-as.integer(block/2)))
intersec<-rbind(intersec1,intersec2)
intersec$clust<-seq(1,nrow(intersec),1)

### boundaries of the CBD (Central Business District)
downtown.xl<-400
downtown.xr<-1200
downtown.yb<-400
downtown.yt<-1200

intersec$DO<-ifelse(intersec$x>=downtown.xl & intersec$x<=downtown.xr & intersec$y>=downtown.yb & intersec$y<=downtown.yt,"D","O")
num_clust<-max(intersec$clust)

e<-c()
for(i in 1:num_clust){
  temp_x<-intersec[intersec$clust==i,]$to_x
  temp_y<-intersec[intersec$clust==i,]$to_y
  clust_temp<-intersec[intersec$from_x==temp_x & intersec$from_y==temp_y,]$clust
  e_temp<-data.frame("from"=rep(i,length(clust_temp)), "to"=as.integer(clust_temp), stringsAsFactors = FALSE)
  e<-rbind(e,e_temp)
}
e<-e[abs(e$from-e$to)!=nrow(intersec1),]
g <- graph.data.frame(e, directed=TRUE)

e1<-e
colnames(e1)<-c("from_clust","to_clust")
e2<-sqldf("SELECT e1.*, intersec.to_x as from_to_x, intersec.to_y as from_to_y, intersec.DO as from_DO
          FROM e1
          LEFT OUTER JOIN intersec ON e1.from_clust = intersec.clust")
e3<-sqldf("SELECT e2.*, intersec.to_x as to_to_x, intersec.to_y as to_to_y, intersec.DO as to_DO
          FROM e2
          LEFT OUTER JOIN intersec ON e2.to_clust = intersec.clust")
e4<-e3



temp1<-sqldf("select a.from_clust, a.to_clust, b.x as from_x, b.y as from_y
      from e4 a
      left outer join intersec b 
      on a.from_clust=b.clust")

temp1<-sqldf("select a.*, b.x as to_x, b.y as to_y
      from temp1 a
      left outer join intersec b 
      on a.to_clust=b.clust")
temp2<-intersec[intersec$DO=="D",c(1,2,7)]
temp3<-merge(temp1, temp2, by = NULL)
temp3$dist_from<-sqrt((temp3$from_x-temp3$x)^2+(temp3$from_y-temp3$y)^2)
temp3$dist_to<-sqrt((temp3$to_x-temp3$x)^2+(temp3$to_y-temp3$y)^2)

temp4<-sqldf("select a.from_clust, a.to_clust, a.from_x, a.from_y, a.to_x, a.to_y, a.x, a.y, a.clust, min(a.dist_from) as min_dist_from
      from temp3 a
      group by a.from_clust, a.to_clust")
temp5<-sqldf("select a.from_clust, a.to_clust, a.from_x, a.from_y, a.to_x, a.to_y, a.x, a.y, a.clust, min(a.dist_to) as min_dist_to
      from temp3 a
      group by a.from_clust, a.to_clust")
temp6<-sqldf("select a.from_clust, a.to_clust, a.from_x, a.from_y, a.to_x, a.to_y, a.x, a.y, a.min_dist_from, b.min_dist_to
      from temp4 a
      inner join temp5 b
      on a.from_clust=b.from_clust and a.to_clust=b.to_clust")
sel<-temp6$from_clust %in% temp2$clust
temp6$from_DO[sel]<-"D"
temp6$from_DO[!sel]<-"O"
sel<-temp6$to_clust %in% temp2$clust
temp6$to_DO[sel]<-"D"
temp6$to_DO[!sel]<-"O"

temp6$direction<-ifelse(temp6$from_DO=="D" & temp6$to_DO=="D" & temp6$min_dist_from==temp6$min_dist_to,"D",
                        ifelse(temp6$from_DO=="O" & temp6$to_DO=="O" & temp6$min_dist_from==temp6$min_dist_to,"O",
                               ifelse(temp6$min_dist_from>temp6$min_dist_to,"D","O")))

dir.final<-temp6


### 1-3. Temporal variation of traffic density and routing: setting routing probability pertaining to vehicle movement --------------------------------------------------------------

gamma_val=1 # set gamma
# gamma > 1 indicates a movement pattern in which peripheral dwellers move into the city center, e.g., during morning commute time; 
# 0 < gamma < 1 indicates a movement pattern in which workers leave the CBD and go back to the periphery, e.g., during evening commute time

neighbor<-vector(mode='list',length=num_clust) # the neighborhood of cluster j is defined as the set of clusters connected from j ∈ V through a directed edge, denoted N_G(j).
neighbor_dir<-vector(mode='list',length=num_clust)
neighbor_dir.prob<-vector(mode='list',length=num_clust)
neighbor_which<-vector(mode='list',length=num_clust)
neighbor_which_dir.prob<-vector(mode='list',length=num_clust)
neighbor_num<-vector(mode='integer',length=num_clust) # the number of neighborhood of cluster j:=|N_G(j)|.
neighbor_inf<-vector(mode='list',length=num_clust)
neighbor_inf_num<-vector(mode='integer',length=num_clust)
for (i in 1:num_clust) {
  temp1<-dir.final[dir.final$from_clust==i,][c("to_clust","direction")]
  neighbor[[i]]<-temp1$to_clust
  neighbor_dir[[i]]<-temp1$direction
  temp2<-ifelse(neighbor_dir[[i]]=="D",gamma_val,1)
  neighbor_dir.prob[[i]]<-temp2/sum(temp2)
  for(j in 1:length(neighbor[[i]])){
    neighbor_which[[neighbor[[i]][j]]]<-c(neighbor_which[[neighbor[[i]][j]]],i)
    neighbor_which_dir.prob[[neighbor[[i]][j]]]<-c(neighbor_which_dir.prob[[neighbor[[i]][j]]],neighbor_dir.prob[[i]][j])
  }
  neighbor_num[i]<-length(neighbor[[i]])
}
for (i in 1:as.integer(num_clust/2)) {
  neighbor_inf[[i]]<-(i+as.integer(num_clust/2))
  neighbor_inf_num[i]<-length(neighbor[[i]])
}
for (i in ((as.integer(num_clust/2)+1):num_clust)) {
  neighbor_inf[[i]]<-(i-as.integer(num_clust/2))
  neighbor_inf_num[i]<-length(neighbor[[i]])
}

# 2. Simulation of Markov chain --------------------------------------------------------------
### 2-1. Creating function: Simulations for information propagation in V2V-enabled transportation network --------------------------------------------------------------
params <- list(beta=3/total_nodes, lambda=0.1, clust_size=num_clust) # parameters - beta:=commnunication rates; lambda:=mobility rates
num_of_ini_inf_veh<-10 # the number of initially informaed vehicles
ini_inf_clust<-1 # location where the initially informaed vehicles are located

Sini_temp<-rep(1,num_clust)*100 # the number of vehicles per cluster at initial time
total_nodes<-sum(Sini_temp) # total number of vehicles
trace_ini<-data.frame("time"=rep(0,num_clust), "clust"=seq(1,num_clust,1), "S"=Sini_temp, "I"=rep(0,num_clust), stringsAsFactors = FALSE) 

epidemic<-function(iteration){  
  ### initial condition
  trace_ini<-data.frame("time"=rep(0,num_clust), "clust"=seq(1,num_clust,1), "S"=Sini_temp, "I"=rep(0,num_clust), stringsAsFactors = FALSE)
  X<-c();Y<-c();tau<-c();tau_cum<-c();temp_type<-c();
  tau_temp_mob.S<-c();tau_temp_mob.I<-c();tau_temp_epi<-c();mob_in_clust<-c()
  next_clust_order<-c(); mob_next_clust<-c(); neigh_out<-c(); neigh_in<-c(); neigh<-c();
  
  X<-trace_ini$S; Y<-trace_ini$I
  X[ini_inf_clust]=X[ini_inf_clust]-num_of_ini_inf_veh; Y[ini_inf_clust]=Y[ini_inf_clust]+num_of_ini_inf_veh; # X and Y denote the number of non-informed and informed vehicles per cluster respectively
  
  output<-array(dim=c(20000000,4))
  time_step<-0; time_cum<-0; time_step_min<-0;
  for(i in 1:nrow(trace_ini)){
    output[i,]<-c(clust=i, time=0, X=X[i], Y=Y[i])
    tau_temp_mob.S<-if(X[i]==0) Inf else rexp(n=1,rate=params$lambda*(X[i])) # exponential delay for mobility of non-inforemd vehicle
    tau_temp_mob.I<-if(Y[i]==0) Inf else rexp(n=1,rate=params$lambda*(Y[i])) # exponential delay for mobility of inforemd vehicle
    tau_temp_epi<-ifelse(X[i]*Y[i]==0,Inf,rexp(n=1,rate=params$beta*(X[i]*Y[i]))) # exponential delay for successful communication within cluster
    tau_temp_epi_inter<-ifelse(X[i]*Y[neighbor_inf[[i]]]==0,Inf,rexp(n=1,rate=params$beta*(X[i]*Y[neighbor_inf[[i]]]))) # exponential delay for successful communication across cluster
    temp_type[i]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter)) # type for minimum exponential delay 
    tau[i] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter)) # delay for minimum exponential delay 
    tau_cum[i]<-tau[i] # cumulative time
  }
  
  ### information propagation simulation
  step<-nrow(trace_ini)
  while(time_step_min<=sim_time){
    time_step_min<-min(tau_cum) # minimum cumulative time
    mob_out_clust<-which(min(tau_cum)==tau_cum) # origin cluster for mobility
    mob_out_type<-temp_type[mob_out_clust] # type for minimum cumulative time
    neigh_out<-neighbor[[mob_out_clust]] # adjacent clusters of the origin cluster in mobile network
    mob_in_clust<-neigh_out[sample(neighbor_num[mob_out_clust], 1, prob=neighbor_dir.prob[[mob_out_clust]])] # destination cluster for mobility

    if(mob_out_type==1){ # mobility of a non-informed vehicle
      X[mob_out_clust] <- X[mob_out_clust]-1; Y[mob_out_clust] <- Y[mob_out_clust];
      X[mob_in_clust] <- X[mob_in_clust]+1; Y[mob_in_clust] <- Y[mob_in_clust];
      step<-step+1
      output[step,]<-c(clust=mob_out_clust, time=time_step_min, X=X[mob_out_clust], Y=Y[mob_out_clust])
      step<-step+1
      output[step,]<-c(clust=mob_in_clust, time=time_step_min, X=X[mob_in_clust], Y=Y[mob_in_clust])

      inf_out_clust<-neighbor_inf[[mob_out_clust]]
      inf_in_clust<-neighbor_inf[[mob_in_clust]]
      
      tau_temp_mob.S<-if(X[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[mob_out_clust]))
      tau_temp_mob.I<-if(Y[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[mob_out_clust]))
      tau_temp_epi<-ifelse(X[mob_out_clust]*Y[mob_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[mob_out_clust])))
      tau_temp_epi_inter<-ifelse(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]])))
      temp_type[mob_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[mob_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[mob_out_clust]<-time_step_min+tau[mob_out_clust]
      
      tau_temp_mob.S<-if(X[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[inf_out_clust]))
      tau_temp_mob.I<-if(Y[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[inf_out_clust]))
      tau_temp_epi<-ifelse(X[inf_out_clust]*Y[inf_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[inf_out_clust])))
      tau_temp_epi_inter<-ifelse(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]])))
      temp_type[inf_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[inf_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[inf_out_clust]<-time_step_min+tau[inf_out_clust]
      
      tau_temp_mob.S<-if(X[mob_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[mob_in_clust]))
      tau_temp_mob.I<-if(Y[mob_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[mob_in_clust]))
      tau_temp_epi<-ifelse(X[mob_in_clust]*Y[mob_in_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_in_clust]*Y[mob_in_clust])))
      tau_temp_epi_inter<-ifelse(X[mob_in_clust]*Y[neighbor_inf[[mob_in_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_in_clust]*Y[neighbor_inf[[mob_in_clust]]])))
      temp_type[mob_in_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[mob_in_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[mob_in_clust]<-time_step_min+tau[mob_in_clust]
      
      tau_temp_mob.S<-if(X[inf_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[inf_in_clust]))
      tau_temp_mob.I<-if(Y[inf_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[inf_in_clust]))
      tau_temp_epi<-ifelse(X[inf_in_clust]*Y[inf_in_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_in_clust]*Y[inf_in_clust])))
      tau_temp_epi_inter<-ifelse(X[inf_in_clust]*Y[neighbor_inf[[inf_in_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_in_clust]*Y[neighbor_inf[[inf_in_clust]]])))
      temp_type[inf_in_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[inf_in_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[inf_in_clust]<-time_step_min+tau[inf_in_clust]
    }
    else if(mob_out_type==2){ # mobility of an informed vehicle
      X[mob_out_clust] <- X[mob_out_clust]; Y[mob_out_clust] <- Y[mob_out_clust]-1;
      X[mob_in_clust] <- X[mob_in_clust]; Y[mob_in_clust] <- Y[mob_in_clust]+1;
      step<-step+1
      output[step,]<-c(clust=mob_out_clust, time=time_step_min, X=X[mob_out_clust], Y=Y[mob_out_clust])
      step<-step+1
      output[step,]<-c(clust=mob_in_clust, time=time_step_min, X=X[mob_in_clust], Y=Y[mob_in_clust])

      inf_out_clust<-neighbor_inf[[mob_out_clust]]
      inf_in_clust<-neighbor_inf[[mob_in_clust]]
      
      tau_temp_mob.S<-if(X[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[mob_out_clust]))
      tau_temp_mob.I<-if(Y[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[mob_out_clust]))
      tau_temp_epi<-ifelse(X[mob_out_clust]*Y[mob_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[mob_out_clust])))
      tau_temp_epi_inter<-ifelse(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]])))
      temp_type[mob_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[mob_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[mob_out_clust]<-time_step_min+tau[mob_out_clust]
      
      tau_temp_mob.S<-if(X[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[inf_out_clust]))
      tau_temp_mob.I<-if(Y[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[inf_out_clust]))
      tau_temp_epi<-ifelse(X[inf_out_clust]*Y[inf_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[inf_out_clust])))
      tau_temp_epi_inter<-ifelse(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]])))
      temp_type[inf_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[inf_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[inf_out_clust]<-time_step_min+tau[inf_out_clust]
      
      tau_temp_mob.S<-if(X[mob_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[mob_in_clust]))
      tau_temp_mob.I<-if(Y[mob_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[mob_in_clust]))
      tau_temp_epi<-ifelse(X[mob_in_clust]*Y[mob_in_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_in_clust]*Y[mob_in_clust])))
      tau_temp_epi_inter<-ifelse(X[mob_in_clust]*Y[neighbor_inf[[mob_in_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_in_clust]*Y[neighbor_inf[[mob_in_clust]]])))
      temp_type[mob_in_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[mob_in_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[mob_in_clust]<-time_step_min+tau[mob_in_clust]
      
      tau_temp_mob.S<-if(X[inf_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[inf_in_clust]))
      tau_temp_mob.I<-if(Y[inf_in_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[inf_in_clust]))
      tau_temp_epi<-ifelse(X[inf_in_clust]*Y[inf_in_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_in_clust]*Y[inf_in_clust])))
      tau_temp_epi_inter<-ifelse(X[inf_in_clust]*Y[neighbor_inf[[inf_in_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_in_clust]*Y[neighbor_inf[[inf_in_clust]]])))
      temp_type[inf_in_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[inf_in_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[inf_in_clust]<-time_step_min+tau[inf_in_clust]
    }
    else{ # successful communication between a non-informed vehicle and an informed vehicle
      X[mob_out_clust] <- X[mob_out_clust]-1; Y[mob_out_clust] <- Y[mob_out_clust]+1;
      step<-step+1
      output[step,]<-c(clust=mob_out_clust, time=time_step_min, X=X[mob_out_clust], Y=Y[mob_out_clust])

      inf_out_clust<-neighbor_inf[[mob_out_clust]]
      
      tau_temp_mob.S<-if(X[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[mob_out_clust]))
      tau_temp_mob.I<-if(Y[mob_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[mob_out_clust]))
      tau_temp_epi<-ifelse(X[mob_out_clust]*Y[mob_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[mob_out_clust])))
      tau_temp_epi_inter<-ifelse(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[mob_out_clust]*Y[neighbor_inf[[mob_out_clust]]])))
      temp_type[mob_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[mob_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[mob_out_clust]<-time_step_min+tau[mob_out_clust]
      
      tau_temp_mob.S<-if(X[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(X[inf_out_clust]))
      tau_temp_mob.I<-if(Y[inf_out_clust]==0) Inf else rexp(n=1,rate=params$lambda*(Y[inf_out_clust]))
      tau_temp_epi<-ifelse(X[inf_out_clust]*Y[inf_out_clust]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[inf_out_clust])))
      tau_temp_epi_inter<-ifelse(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]]==0,Inf,rexp(n=1,rate=params$beta*(X[inf_out_clust]*Y[neighbor_inf[[inf_out_clust]]])))
      temp_type[inf_out_clust]<-which.min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau[inf_out_clust] <- min(c(tau_temp_mob.S,tau_temp_mob.I,tau_temp_epi,tau_temp_epi_inter))
      tau_cum[inf_out_clust]<-time_step_min+tau[inf_out_clust]
    }
  }
  
  ### data cleansing
  # summary: number of non-informed and informed vehicles over time (every 0.1 sec)
  # history: number of non-informed and informed vehicles per cluster over time (every 0.1 sec)

  output<-output[complete.cases(output), ]
  output<-output[order(output[,1],output[,2]),]
  
  floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
  ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
  
  output[,2]<-ceiling_dec(output[,2],1)
  output_final<-output[which(diff(output[,2])!=0),]
  output_final<-data.frame("clust"=output_final[,1],"time"=output_final[,2],"S"=output_final[,3],"I"=output_final[,4],stringsAsFactors = FALSE)
  
  output_final$time<-as.integer(output_final$time*10)

  seq_time_temp<-as.integer(seq(0,sim_time,0.1)*10)
  seq_time<-data.frame("clust"=rep(1:nrow(trace_ini),each=length(seq_time_temp)), "time"=rep(seq_time_temp,times=nrow(trace_ini)), stringsAsFactors = FALSE)
  final_temp <- sqldf("SELECT seq_time.clust as clust, seq_time.time as time, output_final.S as S, output_final.I as I 
                      FROM seq_time
                      LEFT OUTER JOIN output_final USING (clust,time)")
  
  final_temp<-final_temp %>% do(na.locf(.))
  
  final <- sqldf("SELECT time, sum(S) as S, sum(I) as I 
                 FROM final_temp
                 GROUP BY time")
  
  final_temp$time<-final_temp$time/10
  final$time<-final$time/10
  
  list(summary=final,history=final_temp)
}


### 2-2. Running multiple simulations to calculate ensemble average --------------------------------------------------------------
num_iteration<-100 # number of simulation runs
summary<-vector(mode='list',length=num_iteration)
history<-vector(mode='list',length=num_iteration)
for(pp in 1:num_iteration){
  set.seed(pp) # set seed
  temp<-epidemic(pp)
  summary[[pp]] <- temp$summary # summary[[pp]]: number of non-informed and informed vehicles over time (every 0.1 sec) for pp-th simulation run
  history[[pp]] <- temp$history # history[[pp]]: number of non-informed and informed vehicles per cluster over time (every 0.1 sec) for pp-th simulation run
}


# 3. Model solution; set of ordinary differential equations
# 3. Model solution; solution of the set of ordinary differential equations --------------------------------------------------------------
### 3-1. Automatic generation of the set of differential equations --------------------------------------------------------------
dy<-vector(mode='character',length=num_clust)
for (i in 1:num_clust) {
  neigh<-neighbor_which[[i]]
  neigh.prob<-neighbor_which_dir.prob[[i]]
  if(length(neigh)==1){
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + (",neigh.prob[1],")*lambda*y[",neigh[1],"] + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + (",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
  else if(length(neigh)==2){
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + ((",neigh.prob[1],")*lambda*y[",neigh[1],"] + (",neigh.prob[2],")*lambda*y[",neigh[2],"]) + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] + (",neigh.prob[2],")*lambda*y[",(neigh[2])+num_clust,"]) - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
  else if(length(neigh)==3){
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1]),"] + (",neigh.prob[2],")*lambda*y[",(neigh[2]),"] + (",neigh.prob[3],")*lambda*y[",(neigh[3]),"]) + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] + (",neigh.prob[2],")*lambda*y[",(neigh[2])+num_clust,"] + (",neigh.prob[3],")*lambda*y[",(neigh[3])+num_clust,"]) - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
  else if(length(neigh)==4){
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1]),"] + (",neigh.prob[2],")*lambda*y[",(neigh[2]),"] + (",neigh.prob[3],")*lambda*y[",(neigh[3]),"] + (",neigh.prob[4],")*lambda*y[",(neigh[4]),"]) + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] + (",neigh.prob[2],")*lambda*y[",(neigh[2])+num_clust,"] + (",neigh.prob[3],")*lambda*y[",(neigh[3])+num_clust,"] + (",neigh.prob[4],")*lambda*y[",(neigh[4])+num_clust,"]) - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
  else if(length(neigh)==5){
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1]),"] + (",neigh.prob[2],")*lambda*y[",(neigh[2]),"] + (",neigh.prob[3],")*lambda*y[",(neigh[3]),"] + (",neigh.prob[4],")*lambda*y[",(neigh[4]),"] + (",neigh.prob[5],")*lambda*y[",(neigh[5]),"]) + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] + (",neigh.prob[2],")*lambda*y[",(neigh[2])+num_clust,"] + (",neigh.prob[3],")*lambda*y[",(neigh[3])+num_clust,"] + (",neigh.prob[4],")*lambda*y[",(neigh[4])+num_clust,"] + (",neigh.prob[5],")*lambda*y[",(neigh[5])+num_clust,"]) - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
  else{
    dy[i]<-paste("dy",i," <- ","-(lambda*y[",i,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1]),"] + (",neigh.prob[2],")*lambda*y[",(neigh[2]),"] + (",neigh.prob[3],")*lambda*y[",(neigh[3]),"] + (",neigh.prob[4],")*lambda*y[",(neigh[4]),"] + (",neigh.prob[5],")*lambda*y[",(neigh[5]),"] + (",neigh.prob[6],")*lambda*y[",(neigh[6]),"]) + beta*y[",i,"]*y[",i+num_clust,"]","+ beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
    dy[i+num_clust]<-paste("dy",i+num_clust," <- ","-(lambda*y[",i+num_clust,"]) + ((",neigh.prob[1],")*lambda*y[",(neigh[1])+num_clust,"] + (",neigh.prob[2],")*lambda*y[",(neigh[2])+num_clust,"] + (",neigh.prob[3],")*lambda*y[",(neigh[3])+num_clust,"] + (",neigh.prob[4],")*lambda*y[",(neigh[4])+num_clust,"] + (",neigh.prob[5],")*lambda*y[",(neigh[5])+num_clust,"] + (",neigh.prob[6],")*lambda*y[",(neigh[6])+num_clust,"]) - beta*y[",i,"]*y[",i+num_clust,"]","- beta*y[",neighbor_inf[[i]],"]*y[",i+num_clust,"]", sep = "")
  }
}  # automatically generate the set of differential equations, for given gamma, beta, lambda, and arbitrary size of grid topology.
cat(noquote(dy),sep='\n');

dy_name<-c()
for (i in 1:(num_clust*2)) {
  dy_name<-c(dy_name,paste("dy",i,sep=""))
}
cat(noquote(dy_name),sep=',');
### 3-2. Creating function: Set of ordinary differential equations --------------------------------------------------------------
f <- function(t, y, parms) {
  dy1 <- -(lambda*sum(neighbor_dir.prob[[1]])*y[1]) + (neighbor_which_dir.prob[[1]][1])*lambda*y[91] + beta*y[1]*y[121]+ beta*y[61]*y[121]
  dy2 <- -(lambda*sum(neighbor_dir.prob[[2]])*y[2]) + ((neighbor_which_dir.prob[[2]][1])*lambda*y[1] + (neighbor_which_dir.prob[[2]][2])*lambda*y[96]) + beta*y[2]*y[122]+ beta*y[62]*y[122]
  dy3 <- -(lambda*sum(neighbor_dir.prob[[3]])*y[3]) + ((neighbor_which_dir.prob[[3]][1])*lambda*y[2] + (neighbor_which_dir.prob[[3]][2])*lambda*y[101]) + beta*y[3]*y[123]+ beta*y[63]*y[123]
  dy4 <- -(lambda*sum(neighbor_dir.prob[[4]])*y[4]) + ((neighbor_which_dir.prob[[4]][1])*lambda*y[3] + (neighbor_which_dir.prob[[4]][2])*lambda*y[106]) + beta*y[4]*y[124]+ beta*y[64]*y[124]
  dy5 <- -(lambda*sum(neighbor_dir.prob[[5]])*y[5]) + ((neighbor_which_dir.prob[[5]][1])*lambda*y[4] + (neighbor_which_dir.prob[[5]][2])*lambda*y[111]) + beta*y[5]*y[125]+ beta*y[65]*y[125]
  dy6 <- -(lambda*sum(neighbor_dir.prob[[6]])*y[6]) + ((neighbor_which_dir.prob[[6]][1])*lambda*y[31] + (neighbor_which_dir.prob[[6]][2])*lambda*y[92]) + beta*y[6]*y[126]+ beta*y[66]*y[126]
  dy7 <- -(lambda*sum(neighbor_dir.prob[[7]])*y[7]) + ((neighbor_which_dir.prob[[7]][1])*lambda*y[6] + (neighbor_which_dir.prob[[7]][2])*lambda*y[36] + (neighbor_which_dir.prob[[7]][3])*lambda*y[97]) + beta*y[7]*y[127]+ beta*y[67]*y[127]
  dy8 <- -(lambda*sum(neighbor_dir.prob[[8]])*y[8]) + ((neighbor_which_dir.prob[[8]][1])*lambda*y[7] + (neighbor_which_dir.prob[[8]][2])*lambda*y[41] + (neighbor_which_dir.prob[[8]][3])*lambda*y[102]) + beta*y[8]*y[128]+ beta*y[68]*y[128]
  dy9 <- -(lambda*sum(neighbor_dir.prob[[9]])*y[9]) + ((neighbor_which_dir.prob[[9]][1])*lambda*y[8] + (neighbor_which_dir.prob[[9]][2])*lambda*y[46] + (neighbor_which_dir.prob[[9]][3])*lambda*y[107]) + beta*y[9]*y[129]+ beta*y[69]*y[129]
  dy10 <- -(lambda*sum(neighbor_dir.prob[[10]])*y[10]) + ((neighbor_which_dir.prob[[10]][1])*lambda*y[9] + (neighbor_which_dir.prob[[10]][2])*lambda*y[51] + (neighbor_which_dir.prob[[10]][3])*lambda*y[112]) + beta*y[10]*y[130]+ beta*y[70]*y[130]
  dy11 <- -(lambda*sum(neighbor_dir.prob[[11]])*y[11]) + ((neighbor_which_dir.prob[[11]][1])*lambda*y[32] + (neighbor_which_dir.prob[[11]][2])*lambda*y[93]) + beta*y[11]*y[131]+ beta*y[71]*y[131]
  dy12 <- -(lambda*sum(neighbor_dir.prob[[12]])*y[12]) + ((neighbor_which_dir.prob[[12]][1])*lambda*y[11] + (neighbor_which_dir.prob[[12]][2])*lambda*y[37] + (neighbor_which_dir.prob[[12]][3])*lambda*y[98]) + beta*y[12]*y[132]+ beta*y[72]*y[132]
  dy13 <- -(lambda*sum(neighbor_dir.prob[[13]])*y[13]) + ((neighbor_which_dir.prob[[13]][1])*lambda*y[12] + (neighbor_which_dir.prob[[13]][2])*lambda*y[42] + (neighbor_which_dir.prob[[13]][3])*lambda*y[103]) + beta*y[13]*y[133]+ beta*y[73]*y[133]
  dy14 <- -(lambda*sum(neighbor_dir.prob[[14]])*y[14]) + ((neighbor_which_dir.prob[[14]][1])*lambda*y[13] + (neighbor_which_dir.prob[[14]][2])*lambda*y[47] + (neighbor_which_dir.prob[[14]][3])*lambda*y[108]) + beta*y[14]*y[134]+ beta*y[74]*y[134]
  dy15 <- -(lambda*sum(neighbor_dir.prob[[15]])*y[15]) + ((neighbor_which_dir.prob[[15]][1])*lambda*y[14] + (neighbor_which_dir.prob[[15]][2])*lambda*y[52] + (neighbor_which_dir.prob[[15]][3])*lambda*y[113]) + beta*y[15]*y[135]+ beta*y[75]*y[135]
  dy16 <- -(lambda*sum(neighbor_dir.prob[[16]])*y[16]) + ((neighbor_which_dir.prob[[16]][1])*lambda*y[33] + (neighbor_which_dir.prob[[16]][2])*lambda*y[94]) + beta*y[16]*y[136]+ beta*y[76]*y[136]
  dy17 <- -(lambda*sum(neighbor_dir.prob[[17]])*y[17]) + ((neighbor_which_dir.prob[[17]][1])*lambda*y[16] + (neighbor_which_dir.prob[[17]][2])*lambda*y[38] + (neighbor_which_dir.prob[[17]][3])*lambda*y[99]) + beta*y[17]*y[137]+ beta*y[77]*y[137]
  dy18 <- -(lambda*sum(neighbor_dir.prob[[18]])*y[18]) + ((neighbor_which_dir.prob[[18]][1])*lambda*y[17] + (neighbor_which_dir.prob[[18]][2])*lambda*y[43] + (neighbor_which_dir.prob[[18]][3])*lambda*y[104]) + beta*y[18]*y[138]+ beta*y[78]*y[138]
  dy19 <- -(lambda*sum(neighbor_dir.prob[[19]])*y[19]) + ((neighbor_which_dir.prob[[19]][1])*lambda*y[18] + (neighbor_which_dir.prob[[19]][2])*lambda*y[48] + (neighbor_which_dir.prob[[19]][3])*lambda*y[109]) + beta*y[19]*y[139]+ beta*y[79]*y[139]
  dy20 <- -(lambda*sum(neighbor_dir.prob[[20]])*y[20]) + ((neighbor_which_dir.prob[[20]][1])*lambda*y[19] + (neighbor_which_dir.prob[[20]][2])*lambda*y[53] + (neighbor_which_dir.prob[[20]][3])*lambda*y[114]) + beta*y[20]*y[140]+ beta*y[80]*y[140]
  dy21 <- -(lambda*sum(neighbor_dir.prob[[21]])*y[21]) + ((neighbor_which_dir.prob[[21]][1])*lambda*y[34] + (neighbor_which_dir.prob[[21]][2])*lambda*y[95]) + beta*y[21]*y[141]+ beta*y[81]*y[141]
  dy22 <- -(lambda*sum(neighbor_dir.prob[[22]])*y[22]) + ((neighbor_which_dir.prob[[22]][1])*lambda*y[21] + (neighbor_which_dir.prob[[22]][2])*lambda*y[39] + (neighbor_which_dir.prob[[22]][3])*lambda*y[100]) + beta*y[22]*y[142]+ beta*y[82]*y[142]
  dy23 <- -(lambda*sum(neighbor_dir.prob[[23]])*y[23]) + ((neighbor_which_dir.prob[[23]][1])*lambda*y[22] + (neighbor_which_dir.prob[[23]][2])*lambda*y[44] + (neighbor_which_dir.prob[[23]][3])*lambda*y[105]) + beta*y[23]*y[143]+ beta*y[83]*y[143]
  dy24 <- -(lambda*sum(neighbor_dir.prob[[24]])*y[24]) + ((neighbor_which_dir.prob[[24]][1])*lambda*y[23] + (neighbor_which_dir.prob[[24]][2])*lambda*y[49] + (neighbor_which_dir.prob[[24]][3])*lambda*y[110]) + beta*y[24]*y[144]+ beta*y[84]*y[144]
  dy25 <- -(lambda*sum(neighbor_dir.prob[[25]])*y[25]) + ((neighbor_which_dir.prob[[25]][1])*lambda*y[24] + (neighbor_which_dir.prob[[25]][2])*lambda*y[54] + (neighbor_which_dir.prob[[25]][3])*lambda*y[115]) + beta*y[25]*y[145]+ beta*y[85]*y[145]
  dy26 <- -(lambda*sum(neighbor_dir.prob[[26]])*y[26]) + (neighbor_which_dir.prob[[26]][1])*lambda*y[35] + beta*y[26]*y[146]+ beta*y[86]*y[146]
  dy27 <- -(lambda*sum(neighbor_dir.prob[[27]])*y[27]) + ((neighbor_which_dir.prob[[27]][1])*lambda*y[26] + (neighbor_which_dir.prob[[27]][2])*lambda*y[40]) + beta*y[27]*y[147]+ beta*y[87]*y[147]
  dy28 <- -(lambda*sum(neighbor_dir.prob[[28]])*y[28]) + ((neighbor_which_dir.prob[[28]][1])*lambda*y[27] + (neighbor_which_dir.prob[[28]][2])*lambda*y[45]) + beta*y[28]*y[148]+ beta*y[88]*y[148]
  dy29 <- -(lambda*sum(neighbor_dir.prob[[29]])*y[29]) + ((neighbor_which_dir.prob[[29]][1])*lambda*y[28] + (neighbor_which_dir.prob[[29]][2])*lambda*y[50]) + beta*y[29]*y[149]+ beta*y[89]*y[149]
  dy30 <- -(lambda*sum(neighbor_dir.prob[[30]])*y[30]) + ((neighbor_which_dir.prob[[30]][1])*lambda*y[29] + (neighbor_which_dir.prob[[30]][2])*lambda*y[55]) + beta*y[30]*y[150]+ beta*y[90]*y[150]
  dy31 <- -(lambda*sum(neighbor_dir.prob[[31]])*y[31]) + (neighbor_which_dir.prob[[31]][1])*lambda*y[61] + beta*y[31]*y[151]+ beta*y[91]*y[151]
  dy32 <- -(lambda*sum(neighbor_dir.prob[[32]])*y[32]) + ((neighbor_which_dir.prob[[32]][1])*lambda*y[31] + (neighbor_which_dir.prob[[32]][2])*lambda*y[66]) + beta*y[32]*y[152]+ beta*y[92]*y[152]
  dy33 <- -(lambda*sum(neighbor_dir.prob[[33]])*y[33]) + ((neighbor_which_dir.prob[[33]][1])*lambda*y[32] + (neighbor_which_dir.prob[[33]][2])*lambda*y[71]) + beta*y[33]*y[153]+ beta*y[93]*y[153]
  dy34 <- -(lambda*sum(neighbor_dir.prob[[34]])*y[34]) + ((neighbor_which_dir.prob[[34]][1])*lambda*y[33] + (neighbor_which_dir.prob[[34]][2])*lambda*y[76]) + beta*y[34]*y[154]+ beta*y[94]*y[154]
  dy35 <- -(lambda*sum(neighbor_dir.prob[[35]])*y[35]) + ((neighbor_which_dir.prob[[35]][1])*lambda*y[34] + (neighbor_which_dir.prob[[35]][2])*lambda*y[81]) + beta*y[35]*y[155]+ beta*y[95]*y[155]
  dy36 <- -(lambda*sum(neighbor_dir.prob[[36]])*y[36]) + ((neighbor_which_dir.prob[[36]][1])*lambda*y[1] + (neighbor_which_dir.prob[[36]][2])*lambda*y[62]) + beta*y[36]*y[156]+ beta*y[96]*y[156]
  dy37 <- -(lambda*sum(neighbor_dir.prob[[37]])*y[37]) + ((neighbor_which_dir.prob[[37]][1])*lambda*y[6] + (neighbor_which_dir.prob[[37]][2])*lambda*y[36] + (neighbor_which_dir.prob[[37]][3])*lambda*y[67]) + beta*y[37]*y[157]+ beta*y[97]*y[157]
  dy38 <- -(lambda*sum(neighbor_dir.prob[[38]])*y[38]) + ((neighbor_which_dir.prob[[38]][1])*lambda*y[11] + (neighbor_which_dir.prob[[38]][2])*lambda*y[37] + (neighbor_which_dir.prob[[38]][3])*lambda*y[72]) + beta*y[38]*y[158]+ beta*y[98]*y[158]
  dy39 <- -(lambda*sum(neighbor_dir.prob[[39]])*y[39]) + ((neighbor_which_dir.prob[[39]][1])*lambda*y[16] + (neighbor_which_dir.prob[[39]][2])*lambda*y[38] + (neighbor_which_dir.prob[[39]][3])*lambda*y[77]) + beta*y[39]*y[159]+ beta*y[99]*y[159]
  dy40 <- -(lambda*sum(neighbor_dir.prob[[40]])*y[40]) + ((neighbor_which_dir.prob[[40]][1])*lambda*y[21] + (neighbor_which_dir.prob[[40]][2])*lambda*y[39] + (neighbor_which_dir.prob[[40]][3])*lambda*y[82]) + beta*y[40]*y[160]+ beta*y[100]*y[160]
  dy41 <- -(lambda*sum(neighbor_dir.prob[[41]])*y[41]) + ((neighbor_which_dir.prob[[41]][1])*lambda*y[2] + (neighbor_which_dir.prob[[41]][2])*lambda*y[63]) + beta*y[41]*y[161]+ beta*y[101]*y[161]
  dy42 <- -(lambda*sum(neighbor_dir.prob[[42]])*y[42]) + ((neighbor_which_dir.prob[[42]][1])*lambda*y[7] + (neighbor_which_dir.prob[[42]][2])*lambda*y[41] + (neighbor_which_dir.prob[[42]][3])*lambda*y[68]) + beta*y[42]*y[162]+ beta*y[102]*y[162]
  dy43 <- -(lambda*sum(neighbor_dir.prob[[43]])*y[43]) + ((neighbor_which_dir.prob[[43]][1])*lambda*y[12] + (neighbor_which_dir.prob[[43]][2])*lambda*y[42] + (neighbor_which_dir.prob[[43]][3])*lambda*y[73]) + beta*y[43]*y[163]+ beta*y[103]*y[163]
  dy44 <- -(lambda*sum(neighbor_dir.prob[[44]])*y[44]) + ((neighbor_which_dir.prob[[44]][1])*lambda*y[17] + (neighbor_which_dir.prob[[44]][2])*lambda*y[43] + (neighbor_which_dir.prob[[44]][3])*lambda*y[78]) + beta*y[44]*y[164]+ beta*y[104]*y[164]
  dy45 <- -(lambda*sum(neighbor_dir.prob[[45]])*y[45]) + ((neighbor_which_dir.prob[[45]][1])*lambda*y[22] + (neighbor_which_dir.prob[[45]][2])*lambda*y[44] + (neighbor_which_dir.prob[[45]][3])*lambda*y[83]) + beta*y[45]*y[165]+ beta*y[105]*y[165]
  dy46 <- -(lambda*sum(neighbor_dir.prob[[46]])*y[46]) + ((neighbor_which_dir.prob[[46]][1])*lambda*y[3] + (neighbor_which_dir.prob[[46]][2])*lambda*y[64]) + beta*y[46]*y[166]+ beta*y[106]*y[166]
  dy47 <- -(lambda*sum(neighbor_dir.prob[[47]])*y[47]) + ((neighbor_which_dir.prob[[47]][1])*lambda*y[8] + (neighbor_which_dir.prob[[47]][2])*lambda*y[46] + (neighbor_which_dir.prob[[47]][3])*lambda*y[69]) + beta*y[47]*y[167]+ beta*y[107]*y[167]
  dy48 <- -(lambda*sum(neighbor_dir.prob[[48]])*y[48]) + ((neighbor_which_dir.prob[[48]][1])*lambda*y[13] + (neighbor_which_dir.prob[[48]][2])*lambda*y[47] + (neighbor_which_dir.prob[[48]][3])*lambda*y[74]) + beta*y[48]*y[168]+ beta*y[108]*y[168]
  dy49 <- -(lambda*sum(neighbor_dir.prob[[49]])*y[49]) + ((neighbor_which_dir.prob[[49]][1])*lambda*y[18] + (neighbor_which_dir.prob[[49]][2])*lambda*y[48] + (neighbor_which_dir.prob[[49]][3])*lambda*y[79]) + beta*y[49]*y[169]+ beta*y[109]*y[169]
  dy50 <- -(lambda*sum(neighbor_dir.prob[[50]])*y[50]) + ((neighbor_which_dir.prob[[50]][1])*lambda*y[23] + (neighbor_which_dir.prob[[50]][2])*lambda*y[49] + (neighbor_which_dir.prob[[50]][3])*lambda*y[84]) + beta*y[50]*y[170]+ beta*y[110]*y[170]
  dy51 <- -(lambda*sum(neighbor_dir.prob[[51]])*y[51]) + ((neighbor_which_dir.prob[[51]][1])*lambda*y[4] + (neighbor_which_dir.prob[[51]][2])*lambda*y[65]) + beta*y[51]*y[171]+ beta*y[111]*y[171]
  dy52 <- -(lambda*sum(neighbor_dir.prob[[52]])*y[52]) + ((neighbor_which_dir.prob[[52]][1])*lambda*y[9] + (neighbor_which_dir.prob[[52]][2])*lambda*y[51] + (neighbor_which_dir.prob[[52]][3])*lambda*y[70]) + beta*y[52]*y[172]+ beta*y[112]*y[172]
  dy53 <- -(lambda*sum(neighbor_dir.prob[[53]])*y[53]) + ((neighbor_which_dir.prob[[53]][1])*lambda*y[14] + (neighbor_which_dir.prob[[53]][2])*lambda*y[52] + (neighbor_which_dir.prob[[53]][3])*lambda*y[75]) + beta*y[53]*y[173]+ beta*y[113]*y[173]
  dy54 <- -(lambda*sum(neighbor_dir.prob[[54]])*y[54]) + ((neighbor_which_dir.prob[[54]][1])*lambda*y[19] + (neighbor_which_dir.prob[[54]][2])*lambda*y[53] + (neighbor_which_dir.prob[[54]][3])*lambda*y[80]) + beta*y[54]*y[174]+ beta*y[114]*y[174]
  dy55 <- -(lambda*sum(neighbor_dir.prob[[55]])*y[55]) + ((neighbor_which_dir.prob[[55]][1])*lambda*y[24] + (neighbor_which_dir.prob[[55]][2])*lambda*y[54] + (neighbor_which_dir.prob[[55]][3])*lambda*y[85]) + beta*y[55]*y[175]+ beta*y[115]*y[175]
  dy56 <- -(lambda*sum(neighbor_dir.prob[[56]])*y[56]) + (neighbor_which_dir.prob[[56]][1])*lambda*y[5] + beta*y[56]*y[176]+ beta*y[116]*y[176]
  dy57 <- -(lambda*sum(neighbor_dir.prob[[57]])*y[57]) + ((neighbor_which_dir.prob[[57]][1])*lambda*y[10] + (neighbor_which_dir.prob[[57]][2])*lambda*y[56]) + beta*y[57]*y[177]+ beta*y[117]*y[177]
  dy58 <- -(lambda*sum(neighbor_dir.prob[[58]])*y[58]) + ((neighbor_which_dir.prob[[58]][1])*lambda*y[15] + (neighbor_which_dir.prob[[58]][2])*lambda*y[57]) + beta*y[58]*y[178]+ beta*y[118]*y[178]
  dy59 <- -(lambda*sum(neighbor_dir.prob[[59]])*y[59]) + ((neighbor_which_dir.prob[[59]][1])*lambda*y[20] + (neighbor_which_dir.prob[[59]][2])*lambda*y[58]) + beta*y[59]*y[179]+ beta*y[119]*y[179]
  dy60 <- -(lambda*sum(neighbor_dir.prob[[60]])*y[60]) + ((neighbor_which_dir.prob[[60]][1])*lambda*y[25] + (neighbor_which_dir.prob[[60]][2])*lambda*y[59]) + beta*y[60]*y[180]+ beta*y[120]*y[180]
  dy61 <- -(lambda*sum(neighbor_dir.prob[[61]])*y[61]) + ((neighbor_which_dir.prob[[61]][1])*lambda*y[62] + (neighbor_which_dir.prob[[61]][2])*lambda*y[96]) + beta*y[61]*y[181]+ beta*y[1]*y[181]
  dy62 <- -(lambda*sum(neighbor_dir.prob[[62]])*y[62]) + ((neighbor_which_dir.prob[[62]][1])*lambda*y[63] + (neighbor_which_dir.prob[[62]][2])*lambda*y[101]) + beta*y[62]*y[182]+ beta*y[2]*y[182]
  dy63 <- -(lambda*sum(neighbor_dir.prob[[63]])*y[63]) + ((neighbor_which_dir.prob[[63]][1])*lambda*y[64] + (neighbor_which_dir.prob[[63]][2])*lambda*y[106]) + beta*y[63]*y[183]+ beta*y[3]*y[183]
  dy64 <- -(lambda*sum(neighbor_dir.prob[[64]])*y[64]) + ((neighbor_which_dir.prob[[64]][1])*lambda*y[65] + (neighbor_which_dir.prob[[64]][2])*lambda*y[111]) + beta*y[64]*y[184]+ beta*y[4]*y[184]
  dy65 <- -(lambda*sum(neighbor_dir.prob[[65]])*y[65]) + (neighbor_which_dir.prob[[65]][1])*lambda*y[116] + beta*y[65]*y[185]+ beta*y[5]*y[185]
  dy66 <- -(lambda*sum(neighbor_dir.prob[[66]])*y[66]) + ((neighbor_which_dir.prob[[66]][1])*lambda*y[36] + (neighbor_which_dir.prob[[66]][2])*lambda*y[67] + (neighbor_which_dir.prob[[66]][3])*lambda*y[97]) + beta*y[66]*y[186]+ beta*y[6]*y[186]
  dy67 <- -(lambda*sum(neighbor_dir.prob[[67]])*y[67]) + ((neighbor_which_dir.prob[[67]][1])*lambda*y[41] + (neighbor_which_dir.prob[[67]][2])*lambda*y[68] + (neighbor_which_dir.prob[[67]][3])*lambda*y[102]) + beta*y[67]*y[187]+ beta*y[7]*y[187]
  dy68 <- -(lambda*sum(neighbor_dir.prob[[68]])*y[68]) + ((neighbor_which_dir.prob[[68]][1])*lambda*y[46] + (neighbor_which_dir.prob[[68]][2])*lambda*y[69] + (neighbor_which_dir.prob[[68]][3])*lambda*y[107]) + beta*y[68]*y[188]+ beta*y[8]*y[188]
  dy69 <- -(lambda*sum(neighbor_dir.prob[[69]])*y[69]) + ((neighbor_which_dir.prob[[69]][1])*lambda*y[51] + (neighbor_which_dir.prob[[69]][2])*lambda*y[70] + (neighbor_which_dir.prob[[69]][3])*lambda*y[112]) + beta*y[69]*y[189]+ beta*y[9]*y[189]
  dy70 <- -(lambda*sum(neighbor_dir.prob[[70]])*y[70]) + ((neighbor_which_dir.prob[[70]][1])*lambda*y[56] + (neighbor_which_dir.prob[[70]][2])*lambda*y[117]) + beta*y[70]*y[190]+ beta*y[10]*y[190]
  dy71 <- -(lambda*sum(neighbor_dir.prob[[71]])*y[71]) + ((neighbor_which_dir.prob[[71]][1])*lambda*y[37] + (neighbor_which_dir.prob[[71]][2])*lambda*y[72] + (neighbor_which_dir.prob[[71]][3])*lambda*y[98]) + beta*y[71]*y[191]+ beta*y[11]*y[191]
  dy72 <- -(lambda*sum(neighbor_dir.prob[[72]])*y[72]) + ((neighbor_which_dir.prob[[72]][1])*lambda*y[42] + (neighbor_which_dir.prob[[72]][2])*lambda*y[73] + (neighbor_which_dir.prob[[72]][3])*lambda*y[103]) + beta*y[72]*y[192]+ beta*y[12]*y[192]
  dy73 <- -(lambda*sum(neighbor_dir.prob[[73]])*y[73]) + ((neighbor_which_dir.prob[[73]][1])*lambda*y[47] + (neighbor_which_dir.prob[[73]][2])*lambda*y[74] + (neighbor_which_dir.prob[[73]][3])*lambda*y[108]) + beta*y[73]*y[193]+ beta*y[13]*y[193]
  dy74 <- -(lambda*sum(neighbor_dir.prob[[74]])*y[74]) + ((neighbor_which_dir.prob[[74]][1])*lambda*y[52] + (neighbor_which_dir.prob[[74]][2])*lambda*y[75] + (neighbor_which_dir.prob[[74]][3])*lambda*y[113]) + beta*y[74]*y[194]+ beta*y[14]*y[194]
  dy75 <- -(lambda*sum(neighbor_dir.prob[[75]])*y[75]) + ((neighbor_which_dir.prob[[75]][1])*lambda*y[57] + (neighbor_which_dir.prob[[75]][2])*lambda*y[118]) + beta*y[75]*y[195]+ beta*y[15]*y[195]
  dy76 <- -(lambda*sum(neighbor_dir.prob[[76]])*y[76]) + ((neighbor_which_dir.prob[[76]][1])*lambda*y[38] + (neighbor_which_dir.prob[[76]][2])*lambda*y[77] + (neighbor_which_dir.prob[[76]][3])*lambda*y[99]) + beta*y[76]*y[196]+ beta*y[16]*y[196]
  dy77 <- -(lambda*sum(neighbor_dir.prob[[77]])*y[77]) + ((neighbor_which_dir.prob[[77]][1])*lambda*y[43] + (neighbor_which_dir.prob[[77]][2])*lambda*y[78] + (neighbor_which_dir.prob[[77]][3])*lambda*y[104]) + beta*y[77]*y[197]+ beta*y[17]*y[197]
  dy78 <- -(lambda*sum(neighbor_dir.prob[[78]])*y[78]) + ((neighbor_which_dir.prob[[78]][1])*lambda*y[48] + (neighbor_which_dir.prob[[78]][2])*lambda*y[79] + (neighbor_which_dir.prob[[78]][3])*lambda*y[109]) + beta*y[78]*y[198]+ beta*y[18]*y[198]
  dy79 <- -(lambda*sum(neighbor_dir.prob[[79]])*y[79]) + ((neighbor_which_dir.prob[[79]][1])*lambda*y[53] + (neighbor_which_dir.prob[[79]][2])*lambda*y[80] + (neighbor_which_dir.prob[[79]][3])*lambda*y[114]) + beta*y[79]*y[199]+ beta*y[19]*y[199]
  dy80 <- -(lambda*sum(neighbor_dir.prob[[80]])*y[80]) + ((neighbor_which_dir.prob[[80]][1])*lambda*y[58] + (neighbor_which_dir.prob[[80]][2])*lambda*y[119]) + beta*y[80]*y[200]+ beta*y[20]*y[200]
  dy81 <- -(lambda*sum(neighbor_dir.prob[[81]])*y[81]) + ((neighbor_which_dir.prob[[81]][1])*lambda*y[39] + (neighbor_which_dir.prob[[81]][2])*lambda*y[82] + (neighbor_which_dir.prob[[81]][3])*lambda*y[100]) + beta*y[81]*y[201]+ beta*y[21]*y[201]
  dy82 <- -(lambda*sum(neighbor_dir.prob[[82]])*y[82]) + ((neighbor_which_dir.prob[[82]][1])*lambda*y[44] + (neighbor_which_dir.prob[[82]][2])*lambda*y[83] + (neighbor_which_dir.prob[[82]][3])*lambda*y[105]) + beta*y[82]*y[202]+ beta*y[22]*y[202]
  dy83 <- -(lambda*sum(neighbor_dir.prob[[83]])*y[83]) + ((neighbor_which_dir.prob[[83]][1])*lambda*y[49] + (neighbor_which_dir.prob[[83]][2])*lambda*y[84] + (neighbor_which_dir.prob[[83]][3])*lambda*y[110]) + beta*y[83]*y[203]+ beta*y[23]*y[203]
  dy84 <- -(lambda*sum(neighbor_dir.prob[[84]])*y[84]) + ((neighbor_which_dir.prob[[84]][1])*lambda*y[54] + (neighbor_which_dir.prob[[84]][2])*lambda*y[85] + (neighbor_which_dir.prob[[84]][3])*lambda*y[115]) + beta*y[84]*y[204]+ beta*y[24]*y[204]
  dy85 <- -(lambda*sum(neighbor_dir.prob[[85]])*y[85]) + ((neighbor_which_dir.prob[[85]][1])*lambda*y[59] + (neighbor_which_dir.prob[[85]][2])*lambda*y[120]) + beta*y[85]*y[205]+ beta*y[25]*y[205]
  dy86 <- -(lambda*sum(neighbor_dir.prob[[86]])*y[86]) + ((neighbor_which_dir.prob[[86]][1])*lambda*y[40] + (neighbor_which_dir.prob[[86]][2])*lambda*y[87]) + beta*y[86]*y[206]+ beta*y[26]*y[206]
  dy87 <- -(lambda*sum(neighbor_dir.prob[[87]])*y[87]) + ((neighbor_which_dir.prob[[87]][1])*lambda*y[45] + (neighbor_which_dir.prob[[87]][2])*lambda*y[88]) + beta*y[87]*y[207]+ beta*y[27]*y[207]
  dy88 <- -(lambda*sum(neighbor_dir.prob[[88]])*y[88]) + ((neighbor_which_dir.prob[[88]][1])*lambda*y[50] + (neighbor_which_dir.prob[[88]][2])*lambda*y[89]) + beta*y[88]*y[208]+ beta*y[28]*y[208]
  dy89 <- -(lambda*sum(neighbor_dir.prob[[89]])*y[89]) + ((neighbor_which_dir.prob[[89]][1])*lambda*y[55] + (neighbor_which_dir.prob[[89]][2])*lambda*y[90]) + beta*y[89]*y[209]+ beta*y[29]*y[209]
  dy90 <- -(lambda*sum(neighbor_dir.prob[[90]])*y[90]) + (neighbor_which_dir.prob[[90]][1])*lambda*y[60] + beta*y[90]*y[210]+ beta*y[30]*y[210]
  dy91 <- -(lambda*sum(neighbor_dir.prob[[91]])*y[91]) + ((neighbor_which_dir.prob[[91]][1])*lambda*y[66] + (neighbor_which_dir.prob[[91]][2])*lambda*y[92]) + beta*y[91]*y[211]+ beta*y[31]*y[211]
  dy92 <- -(lambda*sum(neighbor_dir.prob[[92]])*y[92]) + ((neighbor_which_dir.prob[[92]][1])*lambda*y[71] + (neighbor_which_dir.prob[[92]][2])*lambda*y[93]) + beta*y[92]*y[212]+ beta*y[32]*y[212]
  dy93 <- -(lambda*sum(neighbor_dir.prob[[93]])*y[93]) + ((neighbor_which_dir.prob[[93]][1])*lambda*y[76] + (neighbor_which_dir.prob[[93]][2])*lambda*y[94]) + beta*y[93]*y[213]+ beta*y[33]*y[213]
  dy94 <- -(lambda*sum(neighbor_dir.prob[[94]])*y[94]) + ((neighbor_which_dir.prob[[94]][1])*lambda*y[81] + (neighbor_which_dir.prob[[94]][2])*lambda*y[95]) + beta*y[94]*y[214]+ beta*y[34]*y[214]
  dy95 <- -(lambda*sum(neighbor_dir.prob[[95]])*y[95]) + (neighbor_which_dir.prob[[95]][1])*lambda*y[86] + beta*y[95]*y[215]+ beta*y[35]*y[215]
  dy96 <- -(lambda*sum(neighbor_dir.prob[[96]])*y[96]) + ((neighbor_which_dir.prob[[96]][1])*lambda*y[6] + (neighbor_which_dir.prob[[96]][2])*lambda*y[67] + (neighbor_which_dir.prob[[96]][3])*lambda*y[97]) + beta*y[96]*y[216]+ beta*y[36]*y[216]
  dy97 <- -(lambda*sum(neighbor_dir.prob[[97]])*y[97]) + ((neighbor_which_dir.prob[[97]][1])*lambda*y[11] + (neighbor_which_dir.prob[[97]][2])*lambda*y[72] + (neighbor_which_dir.prob[[97]][3])*lambda*y[98]) + beta*y[97]*y[217]+ beta*y[37]*y[217]
  dy98 <- -(lambda*sum(neighbor_dir.prob[[98]])*y[98]) + ((neighbor_which_dir.prob[[98]][1])*lambda*y[16] + (neighbor_which_dir.prob[[98]][2])*lambda*y[77] + (neighbor_which_dir.prob[[98]][3])*lambda*y[99]) + beta*y[98]*y[218]+ beta*y[38]*y[218]
  dy99 <- -(lambda*sum(neighbor_dir.prob[[99]])*y[99]) + ((neighbor_which_dir.prob[[99]][1])*lambda*y[21] + (neighbor_which_dir.prob[[99]][2])*lambda*y[82] + (neighbor_which_dir.prob[[99]][3])*lambda*y[100]) + beta*y[99]*y[219]+ beta*y[39]*y[219]
  dy100 <- -(lambda*sum(neighbor_dir.prob[[100]])*y[100]) + ((neighbor_which_dir.prob[[100]][1])*lambda*y[26] + (neighbor_which_dir.prob[[100]][2])*lambda*y[87]) + beta*y[100]*y[220]+ beta*y[40]*y[220]
  dy101 <- -(lambda*sum(neighbor_dir.prob[[101]])*y[101]) + ((neighbor_which_dir.prob[[101]][1])*lambda*y[7] + (neighbor_which_dir.prob[[101]][2])*lambda*y[68] + (neighbor_which_dir.prob[[101]][3])*lambda*y[102]) + beta*y[101]*y[221]+ beta*y[41]*y[221]
  dy102 <- -(lambda*sum(neighbor_dir.prob[[102]])*y[102]) + ((neighbor_which_dir.prob[[102]][1])*lambda*y[12] + (neighbor_which_dir.prob[[102]][2])*lambda*y[73] + (neighbor_which_dir.prob[[102]][3])*lambda*y[103]) + beta*y[102]*y[222]+ beta*y[42]*y[222]
  dy103 <- -(lambda*sum(neighbor_dir.prob[[103]])*y[103]) + ((neighbor_which_dir.prob[[103]][1])*lambda*y[17] + (neighbor_which_dir.prob[[103]][2])*lambda*y[78] + (neighbor_which_dir.prob[[103]][3])*lambda*y[104]) + beta*y[103]*y[223]+ beta*y[43]*y[223]
  dy104 <- -(lambda*sum(neighbor_dir.prob[[104]])*y[104]) + ((neighbor_which_dir.prob[[104]][1])*lambda*y[22] + (neighbor_which_dir.prob[[104]][2])*lambda*y[83] + (neighbor_which_dir.prob[[104]][3])*lambda*y[105]) + beta*y[104]*y[224]+ beta*y[44]*y[224]
  dy105 <- -(lambda*sum(neighbor_dir.prob[[105]])*y[105]) + ((neighbor_which_dir.prob[[105]][1])*lambda*y[27] + (neighbor_which_dir.prob[[105]][2])*lambda*y[88]) + beta*y[105]*y[225]+ beta*y[45]*y[225]
  dy106 <- -(lambda*sum(neighbor_dir.prob[[106]])*y[106]) + ((neighbor_which_dir.prob[[106]][1])*lambda*y[8] + (neighbor_which_dir.prob[[106]][2])*lambda*y[69] + (neighbor_which_dir.prob[[106]][3])*lambda*y[107]) + beta*y[106]*y[226]+ beta*y[46]*y[226]
  dy107 <- -(lambda*sum(neighbor_dir.prob[[107]])*y[107]) + ((neighbor_which_dir.prob[[107]][1])*lambda*y[13] + (neighbor_which_dir.prob[[107]][2])*lambda*y[74] + (neighbor_which_dir.prob[[107]][3])*lambda*y[108]) + beta*y[107]*y[227]+ beta*y[47]*y[227]
  dy108 <- -(lambda*sum(neighbor_dir.prob[[108]])*y[108]) + ((neighbor_which_dir.prob[[108]][1])*lambda*y[18] + (neighbor_which_dir.prob[[108]][2])*lambda*y[79] + (neighbor_which_dir.prob[[108]][3])*lambda*y[109]) + beta*y[108]*y[228]+ beta*y[48]*y[228]
  dy109 <- -(lambda*sum(neighbor_dir.prob[[109]])*y[109]) + ((neighbor_which_dir.prob[[109]][1])*lambda*y[23] + (neighbor_which_dir.prob[[109]][2])*lambda*y[84] + (neighbor_which_dir.prob[[109]][3])*lambda*y[110]) + beta*y[109]*y[229]+ beta*y[49]*y[229]
  dy110 <- -(lambda*sum(neighbor_dir.prob[[110]])*y[110]) + ((neighbor_which_dir.prob[[110]][1])*lambda*y[28] + (neighbor_which_dir.prob[[110]][2])*lambda*y[89]) + beta*y[110]*y[230]+ beta*y[50]*y[230]
  dy111 <- -(lambda*sum(neighbor_dir.prob[[111]])*y[111]) + ((neighbor_which_dir.prob[[111]][1])*lambda*y[9] + (neighbor_which_dir.prob[[111]][2])*lambda*y[70] + (neighbor_which_dir.prob[[111]][3])*lambda*y[112]) + beta*y[111]*y[231]+ beta*y[51]*y[231]
  dy112 <- -(lambda*sum(neighbor_dir.prob[[112]])*y[112]) + ((neighbor_which_dir.prob[[112]][1])*lambda*y[14] + (neighbor_which_dir.prob[[112]][2])*lambda*y[75] + (neighbor_which_dir.prob[[112]][3])*lambda*y[113]) + beta*y[112]*y[232]+ beta*y[52]*y[232]
  dy113 <- -(lambda*sum(neighbor_dir.prob[[113]])*y[113]) + ((neighbor_which_dir.prob[[113]][1])*lambda*y[19] + (neighbor_which_dir.prob[[113]][2])*lambda*y[80] + (neighbor_which_dir.prob[[113]][3])*lambda*y[114]) + beta*y[113]*y[233]+ beta*y[53]*y[233]
  dy114 <- -(lambda*sum(neighbor_dir.prob[[114]])*y[114]) + ((neighbor_which_dir.prob[[114]][1])*lambda*y[24] + (neighbor_which_dir.prob[[114]][2])*lambda*y[85] + (neighbor_which_dir.prob[[114]][3])*lambda*y[115]) + beta*y[114]*y[234]+ beta*y[54]*y[234]
  dy115 <- -(lambda*sum(neighbor_dir.prob[[115]])*y[115]) + ((neighbor_which_dir.prob[[115]][1])*lambda*y[29] + (neighbor_which_dir.prob[[115]][2])*lambda*y[90]) + beta*y[115]*y[235]+ beta*y[55]*y[235]
  dy116 <- -(lambda*sum(neighbor_dir.prob[[116]])*y[116]) + ((neighbor_which_dir.prob[[116]][1])*lambda*y[10] + (neighbor_which_dir.prob[[116]][2])*lambda*y[117]) + beta*y[116]*y[236]+ beta*y[56]*y[236]
  dy117 <- -(lambda*sum(neighbor_dir.prob[[117]])*y[117]) + ((neighbor_which_dir.prob[[117]][1])*lambda*y[15] + (neighbor_which_dir.prob[[117]][2])*lambda*y[118]) + beta*y[117]*y[237]+ beta*y[57]*y[237]
  dy118 <- -(lambda*sum(neighbor_dir.prob[[118]])*y[118]) + ((neighbor_which_dir.prob[[118]][1])*lambda*y[20] + (neighbor_which_dir.prob[[118]][2])*lambda*y[119]) + beta*y[118]*y[238]+ beta*y[58]*y[238]
  dy119 <- -(lambda*sum(neighbor_dir.prob[[119]])*y[119]) + ((neighbor_which_dir.prob[[119]][1])*lambda*y[25] + (neighbor_which_dir.prob[[119]][2])*lambda*y[120]) + beta*y[119]*y[239]+ beta*y[59]*y[239]
  dy120 <- -(lambda*sum(neighbor_dir.prob[[120]])*y[120]) + (neighbor_which_dir.prob[[120]][1])*lambda*y[30] + beta*y[120]*y[240]+ beta*y[60]*y[240]
  dy121 <- -(lambda*sum(neighbor_dir.prob[[1]])*y[121]) + (neighbor_which_dir.prob[[1]][1])*lambda*y[211] - beta*y[1]*y[121]- beta*y[61]*y[121]
  dy122 <- -(lambda*sum(neighbor_dir.prob[[2]])*y[122]) + ((neighbor_which_dir.prob[[2]][1])*lambda*y[121] + (neighbor_which_dir.prob[[2]][2])*lambda*y[216]) - beta*y[2]*y[122]- beta*y[62]*y[122]
  dy123 <- -(lambda*sum(neighbor_dir.prob[[3]])*y[123]) + ((neighbor_which_dir.prob[[3]][1])*lambda*y[122] + (neighbor_which_dir.prob[[3]][2])*lambda*y[221]) - beta*y[3]*y[123]- beta*y[63]*y[123]
  dy124 <- -(lambda*sum(neighbor_dir.prob[[4]])*y[124]) + ((neighbor_which_dir.prob[[4]][1])*lambda*y[123] + (neighbor_which_dir.prob[[4]][2])*lambda*y[226]) - beta*y[4]*y[124]- beta*y[64]*y[124]
  dy125 <- -(lambda*sum(neighbor_dir.prob[[5]])*y[125]) + ((neighbor_which_dir.prob[[5]][1])*lambda*y[124] + (neighbor_which_dir.prob[[5]][2])*lambda*y[231]) - beta*y[5]*y[125]- beta*y[65]*y[125]
  dy126 <- -(lambda*sum(neighbor_dir.prob[[6]])*y[126]) + ((neighbor_which_dir.prob[[6]][1])*lambda*y[151] + (neighbor_which_dir.prob[[6]][2])*lambda*y[212]) - beta*y[6]*y[126]- beta*y[66]*y[126]
  dy127 <- -(lambda*sum(neighbor_dir.prob[[7]])*y[127]) + ((neighbor_which_dir.prob[[7]][1])*lambda*y[126] + (neighbor_which_dir.prob[[7]][2])*lambda*y[156] + (neighbor_which_dir.prob[[7]][3])*lambda*y[217]) - beta*y[7]*y[127]- beta*y[67]*y[127]
  dy128 <- -(lambda*sum(neighbor_dir.prob[[8]])*y[128]) + ((neighbor_which_dir.prob[[8]][1])*lambda*y[127] + (neighbor_which_dir.prob[[8]][2])*lambda*y[161] + (neighbor_which_dir.prob[[8]][3])*lambda*y[222]) - beta*y[8]*y[128]- beta*y[68]*y[128]
  dy129 <- -(lambda*sum(neighbor_dir.prob[[9]])*y[129]) + ((neighbor_which_dir.prob[[9]][1])*lambda*y[128] + (neighbor_which_dir.prob[[9]][2])*lambda*y[166] + (neighbor_which_dir.prob[[9]][3])*lambda*y[227]) - beta*y[9]*y[129]- beta*y[69]*y[129]
  dy130 <- -(lambda*sum(neighbor_dir.prob[[10]])*y[130]) + ((neighbor_which_dir.prob[[10]][1])*lambda*y[129] + (neighbor_which_dir.prob[[10]][2])*lambda*y[171] + (neighbor_which_dir.prob[[10]][3])*lambda*y[232]) - beta*y[10]*y[130]- beta*y[70]*y[130]
  dy131 <- -(lambda*sum(neighbor_dir.prob[[11]])*y[131]) + ((neighbor_which_dir.prob[[11]][1])*lambda*y[152] + (neighbor_which_dir.prob[[11]][2])*lambda*y[213]) - beta*y[11]*y[131]- beta*y[71]*y[131]
  dy132 <- -(lambda*sum(neighbor_dir.prob[[12]])*y[132]) + ((neighbor_which_dir.prob[[12]][1])*lambda*y[131] + (neighbor_which_dir.prob[[12]][2])*lambda*y[157] + (neighbor_which_dir.prob[[12]][3])*lambda*y[218]) - beta*y[12]*y[132]- beta*y[72]*y[132]
  dy133 <- -(lambda*sum(neighbor_dir.prob[[13]])*y[133]) + ((neighbor_which_dir.prob[[13]][1])*lambda*y[132] + (neighbor_which_dir.prob[[13]][2])*lambda*y[162] + (neighbor_which_dir.prob[[13]][3])*lambda*y[223]) - beta*y[13]*y[133]- beta*y[73]*y[133]
  dy134 <- -(lambda*sum(neighbor_dir.prob[[14]])*y[134]) + ((neighbor_which_dir.prob[[14]][1])*lambda*y[133] + (neighbor_which_dir.prob[[14]][2])*lambda*y[167] + (neighbor_which_dir.prob[[14]][3])*lambda*y[228]) - beta*y[14]*y[134]- beta*y[74]*y[134]
  dy135 <- -(lambda*sum(neighbor_dir.prob[[15]])*y[135]) + ((neighbor_which_dir.prob[[15]][1])*lambda*y[134] + (neighbor_which_dir.prob[[15]][2])*lambda*y[172] + (neighbor_which_dir.prob[[15]][3])*lambda*y[233]) - beta*y[15]*y[135]- beta*y[75]*y[135]
  dy136 <- -(lambda*sum(neighbor_dir.prob[[16]])*y[136]) + ((neighbor_which_dir.prob[[16]][1])*lambda*y[153] + (neighbor_which_dir.prob[[16]][2])*lambda*y[214]) - beta*y[16]*y[136]- beta*y[76]*y[136]
  dy137 <- -(lambda*sum(neighbor_dir.prob[[17]])*y[137]) + ((neighbor_which_dir.prob[[17]][1])*lambda*y[136] + (neighbor_which_dir.prob[[17]][2])*lambda*y[158] + (neighbor_which_dir.prob[[17]][3])*lambda*y[219]) - beta*y[17]*y[137]- beta*y[77]*y[137]
  dy138 <- -(lambda*sum(neighbor_dir.prob[[18]])*y[138]) + ((neighbor_which_dir.prob[[18]][1])*lambda*y[137] + (neighbor_which_dir.prob[[18]][2])*lambda*y[163] + (neighbor_which_dir.prob[[18]][3])*lambda*y[224]) - beta*y[18]*y[138]- beta*y[78]*y[138]
  dy139 <- -(lambda*sum(neighbor_dir.prob[[19]])*y[139]) + ((neighbor_which_dir.prob[[19]][1])*lambda*y[138] + (neighbor_which_dir.prob[[19]][2])*lambda*y[168] + (neighbor_which_dir.prob[[19]][3])*lambda*y[229]) - beta*y[19]*y[139]- beta*y[79]*y[139]
  dy140 <- -(lambda*sum(neighbor_dir.prob[[20]])*y[140]) + ((neighbor_which_dir.prob[[20]][1])*lambda*y[139] + (neighbor_which_dir.prob[[20]][2])*lambda*y[173] + (neighbor_which_dir.prob[[20]][3])*lambda*y[234]) - beta*y[20]*y[140]- beta*y[80]*y[140]
  dy141 <- -(lambda*sum(neighbor_dir.prob[[21]])*y[141]) + ((neighbor_which_dir.prob[[21]][1])*lambda*y[154] + (neighbor_which_dir.prob[[21]][2])*lambda*y[215]) - beta*y[21]*y[141]- beta*y[81]*y[141]
  dy142 <- -(lambda*sum(neighbor_dir.prob[[22]])*y[142]) + ((neighbor_which_dir.prob[[22]][1])*lambda*y[141] + (neighbor_which_dir.prob[[22]][2])*lambda*y[159] + (neighbor_which_dir.prob[[22]][3])*lambda*y[220]) - beta*y[22]*y[142]- beta*y[82]*y[142]
  dy143 <- -(lambda*sum(neighbor_dir.prob[[23]])*y[143]) + ((neighbor_which_dir.prob[[23]][1])*lambda*y[142] + (neighbor_which_dir.prob[[23]][2])*lambda*y[164] + (neighbor_which_dir.prob[[23]][3])*lambda*y[225]) - beta*y[23]*y[143]- beta*y[83]*y[143]
  dy144 <- -(lambda*sum(neighbor_dir.prob[[24]])*y[144]) + ((neighbor_which_dir.prob[[24]][1])*lambda*y[143] + (neighbor_which_dir.prob[[24]][2])*lambda*y[169] + (neighbor_which_dir.prob[[24]][3])*lambda*y[230]) - beta*y[24]*y[144]- beta*y[84]*y[144]
  dy145 <- -(lambda*sum(neighbor_dir.prob[[25]])*y[145]) + ((neighbor_which_dir.prob[[25]][1])*lambda*y[144] + (neighbor_which_dir.prob[[25]][2])*lambda*y[174] + (neighbor_which_dir.prob[[25]][3])*lambda*y[235]) - beta*y[25]*y[145]- beta*y[85]*y[145]
  dy146 <- -(lambda*sum(neighbor_dir.prob[[26]])*y[146]) + (neighbor_which_dir.prob[[26]][1])*lambda*y[155] - beta*y[26]*y[146]- beta*y[86]*y[146]
  dy147 <- -(lambda*sum(neighbor_dir.prob[[27]])*y[147]) + ((neighbor_which_dir.prob[[27]][1])*lambda*y[146] + (neighbor_which_dir.prob[[27]][2])*lambda*y[160]) - beta*y[27]*y[147]- beta*y[87]*y[147]
  dy148 <- -(lambda*sum(neighbor_dir.prob[[28]])*y[148]) + ((neighbor_which_dir.prob[[28]][1])*lambda*y[147] + (neighbor_which_dir.prob[[28]][2])*lambda*y[165]) - beta*y[28]*y[148]- beta*y[88]*y[148]
  dy149 <- -(lambda*sum(neighbor_dir.prob[[29]])*y[149]) + ((neighbor_which_dir.prob[[29]][1])*lambda*y[148] + (neighbor_which_dir.prob[[29]][2])*lambda*y[170]) - beta*y[29]*y[149]- beta*y[89]*y[149]
  dy150 <- -(lambda*sum(neighbor_dir.prob[[30]])*y[150]) + ((neighbor_which_dir.prob[[30]][1])*lambda*y[149] + (neighbor_which_dir.prob[[30]][2])*lambda*y[175]) - beta*y[30]*y[150]- beta*y[90]*y[150]
  dy151 <- -(lambda*sum(neighbor_dir.prob[[31]])*y[151]) + (neighbor_which_dir.prob[[31]][1])*lambda*y[181] - beta*y[31]*y[151]- beta*y[91]*y[151]
  dy152 <- -(lambda*sum(neighbor_dir.prob[[32]])*y[152]) + ((neighbor_which_dir.prob[[32]][1])*lambda*y[151] + (neighbor_which_dir.prob[[32]][2])*lambda*y[186]) - beta*y[32]*y[152]- beta*y[92]*y[152]
  dy153 <- -(lambda*sum(neighbor_dir.prob[[33]])*y[153]) + ((neighbor_which_dir.prob[[33]][1])*lambda*y[152] + (neighbor_which_dir.prob[[33]][2])*lambda*y[191]) - beta*y[33]*y[153]- beta*y[93]*y[153]
  dy154 <- -(lambda*sum(neighbor_dir.prob[[34]])*y[154]) + ((neighbor_which_dir.prob[[34]][1])*lambda*y[153] + (neighbor_which_dir.prob[[34]][2])*lambda*y[196]) - beta*y[34]*y[154]- beta*y[94]*y[154]
  dy155 <- -(lambda*sum(neighbor_dir.prob[[35]])*y[155]) + ((neighbor_which_dir.prob[[35]][1])*lambda*y[154] + (neighbor_which_dir.prob[[35]][2])*lambda*y[201]) - beta*y[35]*y[155]- beta*y[95]*y[155]
  dy156 <- -(lambda*sum(neighbor_dir.prob[[36]])*y[156]) + ((neighbor_which_dir.prob[[36]][1])*lambda*y[121] + (neighbor_which_dir.prob[[36]][2])*lambda*y[182]) - beta*y[36]*y[156]- beta*y[96]*y[156]
  dy157 <- -(lambda*sum(neighbor_dir.prob[[37]])*y[157]) + ((neighbor_which_dir.prob[[37]][1])*lambda*y[126] + (neighbor_which_dir.prob[[37]][2])*lambda*y[156] + (neighbor_which_dir.prob[[37]][3])*lambda*y[187]) - beta*y[37]*y[157]- beta*y[97]*y[157]
  dy158 <- -(lambda*sum(neighbor_dir.prob[[38]])*y[158]) + ((neighbor_which_dir.prob[[38]][1])*lambda*y[131] + (neighbor_which_dir.prob[[38]][2])*lambda*y[157] + (neighbor_which_dir.prob[[38]][3])*lambda*y[192]) - beta*y[38]*y[158]- beta*y[98]*y[158]
  dy159 <- -(lambda*sum(neighbor_dir.prob[[39]])*y[159]) + ((neighbor_which_dir.prob[[39]][1])*lambda*y[136] + (neighbor_which_dir.prob[[39]][2])*lambda*y[158] + (neighbor_which_dir.prob[[39]][3])*lambda*y[197]) - beta*y[39]*y[159]- beta*y[99]*y[159]
  dy160 <- -(lambda*sum(neighbor_dir.prob[[40]])*y[160]) + ((neighbor_which_dir.prob[[40]][1])*lambda*y[141] + (neighbor_which_dir.prob[[40]][2])*lambda*y[159] + (neighbor_which_dir.prob[[40]][3])*lambda*y[202]) - beta*y[40]*y[160]- beta*y[100]*y[160]
  dy161 <- -(lambda*sum(neighbor_dir.prob[[41]])*y[161]) + ((neighbor_which_dir.prob[[41]][1])*lambda*y[122] + (neighbor_which_dir.prob[[41]][2])*lambda*y[183]) - beta*y[41]*y[161]- beta*y[101]*y[161]
  dy162 <- -(lambda*sum(neighbor_dir.prob[[42]])*y[162]) + ((neighbor_which_dir.prob[[42]][1])*lambda*y[127] + (neighbor_which_dir.prob[[42]][2])*lambda*y[161] + (neighbor_which_dir.prob[[42]][3])*lambda*y[188]) - beta*y[42]*y[162]- beta*y[102]*y[162]
  dy163 <- -(lambda*sum(neighbor_dir.prob[[43]])*y[163]) + ((neighbor_which_dir.prob[[43]][1])*lambda*y[132] + (neighbor_which_dir.prob[[43]][2])*lambda*y[162] + (neighbor_which_dir.prob[[43]][3])*lambda*y[193]) - beta*y[43]*y[163]- beta*y[103]*y[163]
  dy164 <- -(lambda*sum(neighbor_dir.prob[[44]])*y[164]) + ((neighbor_which_dir.prob[[44]][1])*lambda*y[137] + (neighbor_which_dir.prob[[44]][2])*lambda*y[163] + (neighbor_which_dir.prob[[44]][3])*lambda*y[198]) - beta*y[44]*y[164]- beta*y[104]*y[164]
  dy165 <- -(lambda*sum(neighbor_dir.prob[[45]])*y[165]) + ((neighbor_which_dir.prob[[45]][1])*lambda*y[142] + (neighbor_which_dir.prob[[45]][2])*lambda*y[164] + (neighbor_which_dir.prob[[45]][3])*lambda*y[203]) - beta*y[45]*y[165]- beta*y[105]*y[165]
  dy166 <- -(lambda*sum(neighbor_dir.prob[[46]])*y[166]) + ((neighbor_which_dir.prob[[46]][1])*lambda*y[123] + (neighbor_which_dir.prob[[46]][2])*lambda*y[184]) - beta*y[46]*y[166]- beta*y[106]*y[166]
  dy167 <- -(lambda*sum(neighbor_dir.prob[[47]])*y[167]) + ((neighbor_which_dir.prob[[47]][1])*lambda*y[128] + (neighbor_which_dir.prob[[47]][2])*lambda*y[166] + (neighbor_which_dir.prob[[47]][3])*lambda*y[189]) - beta*y[47]*y[167]- beta*y[107]*y[167]
  dy168 <- -(lambda*sum(neighbor_dir.prob[[48]])*y[168]) + ((neighbor_which_dir.prob[[48]][1])*lambda*y[133] + (neighbor_which_dir.prob[[48]][2])*lambda*y[167] + (neighbor_which_dir.prob[[48]][3])*lambda*y[194]) - beta*y[48]*y[168]- beta*y[108]*y[168]
  dy169 <- -(lambda*sum(neighbor_dir.prob[[49]])*y[169]) + ((neighbor_which_dir.prob[[49]][1])*lambda*y[138] + (neighbor_which_dir.prob[[49]][2])*lambda*y[168] + (neighbor_which_dir.prob[[49]][3])*lambda*y[199]) - beta*y[49]*y[169]- beta*y[109]*y[169]
  dy170 <- -(lambda*sum(neighbor_dir.prob[[50]])*y[170]) + ((neighbor_which_dir.prob[[50]][1])*lambda*y[143] + (neighbor_which_dir.prob[[50]][2])*lambda*y[169] + (neighbor_which_dir.prob[[50]][3])*lambda*y[204]) - beta*y[50]*y[170]- beta*y[110]*y[170]
  dy171 <- -(lambda*sum(neighbor_dir.prob[[51]])*y[171]) + ((neighbor_which_dir.prob[[51]][1])*lambda*y[124] + (neighbor_which_dir.prob[[51]][2])*lambda*y[185]) - beta*y[51]*y[171]- beta*y[111]*y[171]
  dy172 <- -(lambda*sum(neighbor_dir.prob[[52]])*y[172]) + ((neighbor_which_dir.prob[[52]][1])*lambda*y[129] + (neighbor_which_dir.prob[[52]][2])*lambda*y[171] + (neighbor_which_dir.prob[[52]][3])*lambda*y[190]) - beta*y[52]*y[172]- beta*y[112]*y[172]
  dy173 <- -(lambda*sum(neighbor_dir.prob[[53]])*y[173]) + ((neighbor_which_dir.prob[[53]][1])*lambda*y[134] + (neighbor_which_dir.prob[[53]][2])*lambda*y[172] + (neighbor_which_dir.prob[[53]][3])*lambda*y[195]) - beta*y[53]*y[173]- beta*y[113]*y[173]
  dy174 <- -(lambda*sum(neighbor_dir.prob[[54]])*y[174]) + ((neighbor_which_dir.prob[[54]][1])*lambda*y[139] + (neighbor_which_dir.prob[[54]][2])*lambda*y[173] + (neighbor_which_dir.prob[[54]][3])*lambda*y[200]) - beta*y[54]*y[174]- beta*y[114]*y[174]
  dy175 <- -(lambda*sum(neighbor_dir.prob[[55]])*y[175]) + ((neighbor_which_dir.prob[[55]][1])*lambda*y[144] + (neighbor_which_dir.prob[[55]][2])*lambda*y[174] + (neighbor_which_dir.prob[[55]][3])*lambda*y[205]) - beta*y[55]*y[175]- beta*y[115]*y[175]
  dy176 <- -(lambda*sum(neighbor_dir.prob[[56]])*y[176]) + (neighbor_which_dir.prob[[56]][1])*lambda*y[125] - beta*y[56]*y[176]- beta*y[116]*y[176]
  dy177 <- -(lambda*sum(neighbor_dir.prob[[57]])*y[177]) + ((neighbor_which_dir.prob[[57]][1])*lambda*y[130] + (neighbor_which_dir.prob[[57]][2])*lambda*y[176]) - beta*y[57]*y[177]- beta*y[117]*y[177]
  dy178 <- -(lambda*sum(neighbor_dir.prob[[58]])*y[178]) + ((neighbor_which_dir.prob[[58]][1])*lambda*y[135] + (neighbor_which_dir.prob[[58]][2])*lambda*y[177]) - beta*y[58]*y[178]- beta*y[118]*y[178]
  dy179 <- -(lambda*sum(neighbor_dir.prob[[59]])*y[179]) + ((neighbor_which_dir.prob[[59]][1])*lambda*y[140] + (neighbor_which_dir.prob[[59]][2])*lambda*y[178]) - beta*y[59]*y[179]- beta*y[119]*y[179]
  dy180 <- -(lambda*sum(neighbor_dir.prob[[60]])*y[180]) + ((neighbor_which_dir.prob[[60]][1])*lambda*y[145] + (neighbor_which_dir.prob[[60]][2])*lambda*y[179]) - beta*y[60]*y[180]- beta*y[120]*y[180]
  dy181 <- -(lambda*sum(neighbor_dir.prob[[61]])*y[181]) + ((neighbor_which_dir.prob[[61]][1])*lambda*y[182] + (neighbor_which_dir.prob[[61]][2])*lambda*y[216]) - beta*y[61]*y[181]- beta*y[1]*y[181]
  dy182 <- -(lambda*sum(neighbor_dir.prob[[62]])*y[182]) + ((neighbor_which_dir.prob[[62]][1])*lambda*y[183] + (neighbor_which_dir.prob[[62]][2])*lambda*y[221]) - beta*y[62]*y[182]- beta*y[2]*y[182]
  dy183 <- -(lambda*sum(neighbor_dir.prob[[63]])*y[183]) + ((neighbor_which_dir.prob[[63]][1])*lambda*y[184] + (neighbor_which_dir.prob[[63]][2])*lambda*y[226]) - beta*y[63]*y[183]- beta*y[3]*y[183]
  dy184 <- -(lambda*sum(neighbor_dir.prob[[64]])*y[184]) + ((neighbor_which_dir.prob[[64]][1])*lambda*y[185] + (neighbor_which_dir.prob[[64]][2])*lambda*y[231]) - beta*y[64]*y[184]- beta*y[4]*y[184]
  dy185 <- -(lambda*sum(neighbor_dir.prob[[65]])*y[185]) + (neighbor_which_dir.prob[[65]][1])*lambda*y[236] - beta*y[65]*y[185]- beta*y[5]*y[185]
  dy186 <- -(lambda*sum(neighbor_dir.prob[[66]])*y[186]) + ((neighbor_which_dir.prob[[66]][1])*lambda*y[156] + (neighbor_which_dir.prob[[66]][2])*lambda*y[187] + (neighbor_which_dir.prob[[66]][3])*lambda*y[217]) - beta*y[66]*y[186]- beta*y[6]*y[186]
  dy187 <- -(lambda*sum(neighbor_dir.prob[[67]])*y[187]) + ((neighbor_which_dir.prob[[67]][1])*lambda*y[161] + (neighbor_which_dir.prob[[67]][2])*lambda*y[188] + (neighbor_which_dir.prob[[67]][3])*lambda*y[222]) - beta*y[67]*y[187]- beta*y[7]*y[187]
  dy188 <- -(lambda*sum(neighbor_dir.prob[[68]])*y[188]) + ((neighbor_which_dir.prob[[68]][1])*lambda*y[166] + (neighbor_which_dir.prob[[68]][2])*lambda*y[189] + (neighbor_which_dir.prob[[68]][3])*lambda*y[227]) - beta*y[68]*y[188]- beta*y[8]*y[188]
  dy189 <- -(lambda*sum(neighbor_dir.prob[[69]])*y[189]) + ((neighbor_which_dir.prob[[69]][1])*lambda*y[171] + (neighbor_which_dir.prob[[69]][2])*lambda*y[190] + (neighbor_which_dir.prob[[69]][3])*lambda*y[232]) - beta*y[69]*y[189]- beta*y[9]*y[189]
  dy190 <- -(lambda*sum(neighbor_dir.prob[[70]])*y[190]) + ((neighbor_which_dir.prob[[70]][1])*lambda*y[176] + (neighbor_which_dir.prob[[70]][2])*lambda*y[237]) - beta*y[70]*y[190]- beta*y[10]*y[190]
  dy191 <- -(lambda*sum(neighbor_dir.prob[[71]])*y[191]) + ((neighbor_which_dir.prob[[71]][1])*lambda*y[157] + (neighbor_which_dir.prob[[71]][2])*lambda*y[192] + (neighbor_which_dir.prob[[71]][3])*lambda*y[218]) - beta*y[71]*y[191]- beta*y[11]*y[191]
  dy192 <- -(lambda*sum(neighbor_dir.prob[[72]])*y[192]) + ((neighbor_which_dir.prob[[72]][1])*lambda*y[162] + (neighbor_which_dir.prob[[72]][2])*lambda*y[193] + (neighbor_which_dir.prob[[72]][3])*lambda*y[223]) - beta*y[72]*y[192]- beta*y[12]*y[192]
  dy193 <- -(lambda*sum(neighbor_dir.prob[[73]])*y[193]) + ((neighbor_which_dir.prob[[73]][1])*lambda*y[167] + (neighbor_which_dir.prob[[73]][2])*lambda*y[194] + (neighbor_which_dir.prob[[73]][3])*lambda*y[228]) - beta*y[73]*y[193]- beta*y[13]*y[193]
  dy194 <- -(lambda*sum(neighbor_dir.prob[[74]])*y[194]) + ((neighbor_which_dir.prob[[74]][1])*lambda*y[172] + (neighbor_which_dir.prob[[74]][2])*lambda*y[195] + (neighbor_which_dir.prob[[74]][3])*lambda*y[233]) - beta*y[74]*y[194]- beta*y[14]*y[194]
  dy195 <- -(lambda*sum(neighbor_dir.prob[[75]])*y[195]) + ((neighbor_which_dir.prob[[75]][1])*lambda*y[177] + (neighbor_which_dir.prob[[75]][2])*lambda*y[238]) - beta*y[75]*y[195]- beta*y[15]*y[195]
  dy196 <- -(lambda*sum(neighbor_dir.prob[[76]])*y[196]) + ((neighbor_which_dir.prob[[76]][1])*lambda*y[158] + (neighbor_which_dir.prob[[76]][2])*lambda*y[197] + (neighbor_which_dir.prob[[76]][3])*lambda*y[219]) - beta*y[76]*y[196]- beta*y[16]*y[196]
  dy197 <- -(lambda*sum(neighbor_dir.prob[[77]])*y[197]) + ((neighbor_which_dir.prob[[77]][1])*lambda*y[163] + (neighbor_which_dir.prob[[77]][2])*lambda*y[198] + (neighbor_which_dir.prob[[77]][3])*lambda*y[224]) - beta*y[77]*y[197]- beta*y[17]*y[197]
  dy198 <- -(lambda*sum(neighbor_dir.prob[[78]])*y[198]) + ((neighbor_which_dir.prob[[78]][1])*lambda*y[168] + (neighbor_which_dir.prob[[78]][2])*lambda*y[199] + (neighbor_which_dir.prob[[78]][3])*lambda*y[229]) - beta*y[78]*y[198]- beta*y[18]*y[198]
  dy199 <- -(lambda*sum(neighbor_dir.prob[[79]])*y[199]) + ((neighbor_which_dir.prob[[79]][1])*lambda*y[173] + (neighbor_which_dir.prob[[79]][2])*lambda*y[200] + (neighbor_which_dir.prob[[79]][3])*lambda*y[234]) - beta*y[79]*y[199]- beta*y[19]*y[199]
  dy200 <- -(lambda*sum(neighbor_dir.prob[[80]])*y[200]) + ((neighbor_which_dir.prob[[80]][1])*lambda*y[178] + (neighbor_which_dir.prob[[80]][2])*lambda*y[239]) - beta*y[80]*y[200]- beta*y[20]*y[200]
  dy201 <- -(lambda*sum(neighbor_dir.prob[[81]])*y[201]) + ((neighbor_which_dir.prob[[81]][1])*lambda*y[159] + (neighbor_which_dir.prob[[81]][2])*lambda*y[202] + (neighbor_which_dir.prob[[81]][3])*lambda*y[220]) - beta*y[81]*y[201]- beta*y[21]*y[201]
  dy202 <- -(lambda*sum(neighbor_dir.prob[[82]])*y[202]) + ((neighbor_which_dir.prob[[82]][1])*lambda*y[164] + (neighbor_which_dir.prob[[82]][2])*lambda*y[203] + (neighbor_which_dir.prob[[82]][3])*lambda*y[225]) - beta*y[82]*y[202]- beta*y[22]*y[202]
  dy203 <- -(lambda*sum(neighbor_dir.prob[[83]])*y[203]) + ((neighbor_which_dir.prob[[83]][1])*lambda*y[169] + (neighbor_which_dir.prob[[83]][2])*lambda*y[204] + (neighbor_which_dir.prob[[83]][3])*lambda*y[230]) - beta*y[83]*y[203]- beta*y[23]*y[203]
  dy204 <- -(lambda*sum(neighbor_dir.prob[[84]])*y[204]) + ((neighbor_which_dir.prob[[84]][1])*lambda*y[174] + (neighbor_which_dir.prob[[84]][2])*lambda*y[205] + (neighbor_which_dir.prob[[84]][3])*lambda*y[235]) - beta*y[84]*y[204]- beta*y[24]*y[204]
  dy205 <- -(lambda*sum(neighbor_dir.prob[[85]])*y[205]) + ((neighbor_which_dir.prob[[85]][1])*lambda*y[179] + (neighbor_which_dir.prob[[85]][2])*lambda*y[240]) - beta*y[85]*y[205]- beta*y[25]*y[205]
  dy206 <- -(lambda*sum(neighbor_dir.prob[[86]])*y[206]) + ((neighbor_which_dir.prob[[86]][1])*lambda*y[160] + (neighbor_which_dir.prob[[86]][2])*lambda*y[207]) - beta*y[86]*y[206]- beta*y[26]*y[206]
  dy207 <- -(lambda*sum(neighbor_dir.prob[[87]])*y[207]) + ((neighbor_which_dir.prob[[87]][1])*lambda*y[165] + (neighbor_which_dir.prob[[87]][2])*lambda*y[208]) - beta*y[87]*y[207]- beta*y[27]*y[207]
  dy208 <- -(lambda*sum(neighbor_dir.prob[[88]])*y[208]) + ((neighbor_which_dir.prob[[88]][1])*lambda*y[170] + (neighbor_which_dir.prob[[88]][2])*lambda*y[209]) - beta*y[88]*y[208]- beta*y[28]*y[208]
  dy209 <- -(lambda*sum(neighbor_dir.prob[[89]])*y[209]) + ((neighbor_which_dir.prob[[89]][1])*lambda*y[175] + (neighbor_which_dir.prob[[89]][2])*lambda*y[210]) - beta*y[89]*y[209]- beta*y[29]*y[209]
  dy210 <- -(lambda*sum(neighbor_dir.prob[[90]])*y[210]) + (neighbor_which_dir.prob[[90]][1])*lambda*y[180] - beta*y[90]*y[210]- beta*y[30]*y[210]
  dy211 <- -(lambda*sum(neighbor_dir.prob[[91]])*y[211]) + ((neighbor_which_dir.prob[[91]][1])*lambda*y[186] + (neighbor_which_dir.prob[[91]][2])*lambda*y[212]) - beta*y[91]*y[211]- beta*y[31]*y[211]
  dy212 <- -(lambda*sum(neighbor_dir.prob[[92]])*y[212]) + ((neighbor_which_dir.prob[[92]][1])*lambda*y[191] + (neighbor_which_dir.prob[[92]][2])*lambda*y[213]) - beta*y[92]*y[212]- beta*y[32]*y[212]
  dy213 <- -(lambda*sum(neighbor_dir.prob[[93]])*y[213]) + ((neighbor_which_dir.prob[[93]][1])*lambda*y[196] + (neighbor_which_dir.prob[[93]][2])*lambda*y[214]) - beta*y[93]*y[213]- beta*y[33]*y[213]
  dy214 <- -(lambda*sum(neighbor_dir.prob[[94]])*y[214]) + ((neighbor_which_dir.prob[[94]][1])*lambda*y[201] + (neighbor_which_dir.prob[[94]][2])*lambda*y[215]) - beta*y[94]*y[214]- beta*y[34]*y[214]
  dy215 <- -(lambda*sum(neighbor_dir.prob[[95]])*y[215]) + (neighbor_which_dir.prob[[95]][1])*lambda*y[206] - beta*y[95]*y[215]- beta*y[35]*y[215]
  dy216 <- -(lambda*sum(neighbor_dir.prob[[96]])*y[216]) + ((neighbor_which_dir.prob[[96]][1])*lambda*y[126] + (neighbor_which_dir.prob[[96]][2])*lambda*y[187] + (neighbor_which_dir.prob[[96]][3])*lambda*y[217]) - beta*y[96]*y[216]- beta*y[36]*y[216]
  dy217 <- -(lambda*sum(neighbor_dir.prob[[97]])*y[217]) + ((neighbor_which_dir.prob[[97]][1])*lambda*y[131] + (neighbor_which_dir.prob[[97]][2])*lambda*y[192] + (neighbor_which_dir.prob[[97]][3])*lambda*y[218]) - beta*y[97]*y[217]- beta*y[37]*y[217]
  dy218 <- -(lambda*sum(neighbor_dir.prob[[98]])*y[218]) + ((neighbor_which_dir.prob[[98]][1])*lambda*y[136] + (neighbor_which_dir.prob[[98]][2])*lambda*y[197] + (neighbor_which_dir.prob[[98]][3])*lambda*y[219]) - beta*y[98]*y[218]- beta*y[38]*y[218]
  dy219 <- -(lambda*sum(neighbor_dir.prob[[99]])*y[219]) + ((neighbor_which_dir.prob[[99]][1])*lambda*y[141] + (neighbor_which_dir.prob[[99]][2])*lambda*y[202] + (neighbor_which_dir.prob[[99]][3])*lambda*y[220]) - beta*y[99]*y[219]- beta*y[39]*y[219]
  dy220 <- -(lambda*sum(neighbor_dir.prob[[100]])*y[220]) + ((neighbor_which_dir.prob[[100]][1])*lambda*y[146] + (neighbor_which_dir.prob[[100]][2])*lambda*y[207]) - beta*y[100]*y[220]- beta*y[40]*y[220]
  dy221 <- -(lambda*sum(neighbor_dir.prob[[101]])*y[221]) + ((neighbor_which_dir.prob[[101]][1])*lambda*y[127] + (neighbor_which_dir.prob[[101]][2])*lambda*y[188] + (neighbor_which_dir.prob[[101]][3])*lambda*y[222]) - beta*y[101]*y[221]- beta*y[41]*y[221]
  dy222 <- -(lambda*sum(neighbor_dir.prob[[102]])*y[222]) + ((neighbor_which_dir.prob[[102]][1])*lambda*y[132] + (neighbor_which_dir.prob[[102]][2])*lambda*y[193] + (neighbor_which_dir.prob[[102]][3])*lambda*y[223]) - beta*y[102]*y[222]- beta*y[42]*y[222]
  dy223 <- -(lambda*sum(neighbor_dir.prob[[103]])*y[223]) + ((neighbor_which_dir.prob[[103]][1])*lambda*y[137] + (neighbor_which_dir.prob[[103]][2])*lambda*y[198] + (neighbor_which_dir.prob[[103]][3])*lambda*y[224]) - beta*y[103]*y[223]- beta*y[43]*y[223]
  dy224 <- -(lambda*sum(neighbor_dir.prob[[104]])*y[224]) + ((neighbor_which_dir.prob[[104]][1])*lambda*y[142] + (neighbor_which_dir.prob[[104]][2])*lambda*y[203] + (neighbor_which_dir.prob[[104]][3])*lambda*y[225]) - beta*y[104]*y[224]- beta*y[44]*y[224]
  dy225 <- -(lambda*sum(neighbor_dir.prob[[105]])*y[225]) + ((neighbor_which_dir.prob[[105]][1])*lambda*y[147] + (neighbor_which_dir.prob[[105]][2])*lambda*y[208]) - beta*y[105]*y[225]- beta*y[45]*y[225]
  dy226 <- -(lambda*sum(neighbor_dir.prob[[106]])*y[226]) + ((neighbor_which_dir.prob[[106]][1])*lambda*y[128] + (neighbor_which_dir.prob[[106]][2])*lambda*y[189] + (neighbor_which_dir.prob[[106]][3])*lambda*y[227]) - beta*y[106]*y[226]- beta*y[46]*y[226]
  dy227 <- -(lambda*sum(neighbor_dir.prob[[107]])*y[227]) + ((neighbor_which_dir.prob[[107]][1])*lambda*y[133] + (neighbor_which_dir.prob[[107]][2])*lambda*y[194] + (neighbor_which_dir.prob[[107]][3])*lambda*y[228]) - beta*y[107]*y[227]- beta*y[47]*y[227]
  dy228 <- -(lambda*sum(neighbor_dir.prob[[108]])*y[228]) + ((neighbor_which_dir.prob[[108]][1])*lambda*y[138] + (neighbor_which_dir.prob[[108]][2])*lambda*y[199] + (neighbor_which_dir.prob[[108]][3])*lambda*y[229]) - beta*y[108]*y[228]- beta*y[48]*y[228]
  dy229 <- -(lambda*sum(neighbor_dir.prob[[109]])*y[229]) + ((neighbor_which_dir.prob[[109]][1])*lambda*y[143] + (neighbor_which_dir.prob[[109]][2])*lambda*y[204] + (neighbor_which_dir.prob[[109]][3])*lambda*y[230]) - beta*y[109]*y[229]- beta*y[49]*y[229]
  dy230 <- -(lambda*sum(neighbor_dir.prob[[110]])*y[230]) + ((neighbor_which_dir.prob[[110]][1])*lambda*y[148] + (neighbor_which_dir.prob[[110]][2])*lambda*y[209]) - beta*y[110]*y[230]- beta*y[50]*y[230]
  dy231 <- -(lambda*sum(neighbor_dir.prob[[111]])*y[231]) + ((neighbor_which_dir.prob[[111]][1])*lambda*y[129] + (neighbor_which_dir.prob[[111]][2])*lambda*y[190] + (neighbor_which_dir.prob[[111]][3])*lambda*y[232]) - beta*y[111]*y[231]- beta*y[51]*y[231]
  dy232 <- -(lambda*sum(neighbor_dir.prob[[112]])*y[232]) + ((neighbor_which_dir.prob[[112]][1])*lambda*y[134] + (neighbor_which_dir.prob[[112]][2])*lambda*y[195] + (neighbor_which_dir.prob[[112]][3])*lambda*y[233]) - beta*y[112]*y[232]- beta*y[52]*y[232]
  dy233 <- -(lambda*sum(neighbor_dir.prob[[113]])*y[233]) + ((neighbor_which_dir.prob[[113]][1])*lambda*y[139] + (neighbor_which_dir.prob[[113]][2])*lambda*y[200] + (neighbor_which_dir.prob[[113]][3])*lambda*y[234]) - beta*y[113]*y[233]- beta*y[53]*y[233]
  dy234 <- -(lambda*sum(neighbor_dir.prob[[114]])*y[234]) + ((neighbor_which_dir.prob[[114]][1])*lambda*y[144] + (neighbor_which_dir.prob[[114]][2])*lambda*y[205] + (neighbor_which_dir.prob[[114]][3])*lambda*y[235]) - beta*y[114]*y[234]- beta*y[54]*y[234]
  dy235 <- -(lambda*sum(neighbor_dir.prob[[115]])*y[235]) + ((neighbor_which_dir.prob[[115]][1])*lambda*y[149] + (neighbor_which_dir.prob[[115]][2])*lambda*y[210]) - beta*y[115]*y[235]- beta*y[55]*y[235]
  dy236 <- -(lambda*sum(neighbor_dir.prob[[116]])*y[236]) + ((neighbor_which_dir.prob[[116]][1])*lambda*y[130] + (neighbor_which_dir.prob[[116]][2])*lambda*y[237]) - beta*y[116]*y[236]- beta*y[56]*y[236]
  dy237 <- -(lambda*sum(neighbor_dir.prob[[117]])*y[237]) + ((neighbor_which_dir.prob[[117]][1])*lambda*y[135] + (neighbor_which_dir.prob[[117]][2])*lambda*y[238]) - beta*y[117]*y[237]- beta*y[57]*y[237]
  dy238 <- -(lambda*sum(neighbor_dir.prob[[118]])*y[238]) + ((neighbor_which_dir.prob[[118]][1])*lambda*y[140] + (neighbor_which_dir.prob[[118]][2])*lambda*y[239]) - beta*y[118]*y[238]- beta*y[58]*y[238]
  dy239 <- -(lambda*sum(neighbor_dir.prob[[119]])*y[239]) + ((neighbor_which_dir.prob[[119]][1])*lambda*y[145] + (neighbor_which_dir.prob[[119]][2])*lambda*y[240]) - beta*y[119]*y[239]- beta*y[59]*y[239]
  dy240 <- -(lambda*sum(neighbor_dir.prob[[120]])*y[240]) + (neighbor_which_dir.prob[[120]][1])*lambda*y[150] - beta*y[120]*y[240]- beta*y[60]*y[240]
  list(c(dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17,dy18,dy19,dy20,dy21,dy22,dy23,dy24,dy25,dy26,dy27,dy28,dy29,dy30,dy31,dy32,dy33,dy34,dy35,dy36,dy37,dy38,dy39,dy40,dy41,dy42,dy43,dy44,dy45,dy46,dy47,dy48,dy49,dy50,dy51,dy52,dy53,dy54,dy55,dy56,dy57,dy58,dy59,dy60,dy61,dy62,dy63,dy64,dy65,dy66,dy67,dy68,dy69,dy70,dy71,dy72,dy73,dy74,dy75,dy76,dy77,dy78,dy79,dy80,dy81,dy82,dy83,dy84,dy85,dy86,dy87,dy88,dy89,dy90,dy91,dy92,dy93,dy94,dy95,dy96,dy97,dy98,dy99,dy100,dy101,dy102,dy103,dy104,dy105,dy106,dy107,dy108,dy109,dy110,dy111,dy112,dy113,dy114,dy115,dy116,dy117,dy118,dy119,dy120,dy121,dy122,dy123,dy124,dy125,dy126,dy127,dy128,dy129,dy130,dy131,dy132,dy133,dy134,dy135,dy136,dy137,dy138,dy139,dy140,dy141,dy142,dy143,dy144,dy145,dy146,dy147,dy148,dy149,dy150,dy151,dy152,dy153,dy154,dy155,dy156,dy157,dy158,dy159,dy160,dy161,dy162,dy163,dy164,dy165,dy166,dy167,dy168,dy169,dy170,dy171,dy172,dy173,dy174,dy175,dy176,dy177,dy178,dy179,dy180,dy181,dy182,dy183,dy184,dy185,dy186,dy187,dy188,dy189,dy190,dy191,dy192,dy193,dy194,dy195,dy196,dy197,dy198,dy199,dy200,dy201,dy202,dy203,dy204,dy205,dy206,dy207,dy208,dy209,dy210,dy211,dy212,dy213,dy214,dy215,dy216,dy217,dy218,dy219,dy220,dy221,dy222,dy223,dy224,dy225,dy226,dy227,dy228,dy229,dy230,dy231,dy232,dy233,dy234,dy235,dy236,dy237,dy238,dy239,dy240))
}
### 3-3. Getting numerical solution of the set of the differential equations --------------------------------------------------------------
ini_inf_clust<-1
num_of_ini_inf_veh<-10
Iini<-trace_ini$I; Sini<-trace_ini$S;
Iini[ini_inf_clust]<-Iini[ini_inf_clust]+num_of_ini_inf_veh; Sini[ini_inf_clust]<-Sini[ini_inf_clust]-num_of_ini_inf_veh
yini<-c(Iini/total_nodes,Sini/total_nodes)
  
times <- seq(from = 0, to = sim_time, by = .1)
out <- ode(times = times, y=yini, func = f, parms = NULL) # solve the set of differential equations numerically
analysis<-rowSums(out[,2:(1+num_clust)])
plot(times,analysis)

# 4. Geographical representation of traffic density and information propagation for both simulation results and model solutions --------------------------------------------------------------
history<-get(load("history_2Dgrid_2way_downtown_400_1200_OD5_dep_arr_a1_b70_clust120_mobrate0_1_infrate_3_nodes12000_avg200_pz_10veh_in_clust1.RData")) # import "history" dataset
clust_avg<-data.frame("clust"=history[[1]]$clust, "time"=history[[1]]$time, "S"=Reduce('+',history)$S/length(history), "I"=Reduce('+',history)$I/length(history))

clust_avg_1<-clust_avg[(1:(nrow(clust_avg)/2)),]
clust_avg_2<-clust_avg[((nrow(clust_avg)/2+1):nrow(clust_avg)),]
clust_avg_final<-clust_avg_1
clust_avg_final$S<-clust_avg_1$S+clust_avg_2$S
clust_avg_final$I<-clust_avg_1$I+clust_avg_2$I

out_1<-out[,2:(1+num_clust)]
out_2<-out[,(2+num_clust):(1+2*num_clust)]

out_I<-out_1[,1:(ncol(out_1)/2)]+out_1[,(ncol(out_1)/2+1):(ncol(out_1))]
out_S<-out_2[,1:(ncol(out_2)/2)]+out_2[,(ncol(out_2)/2+1):(ncol(out_2))]
out_I<-cbind(out[,1],out_I)
out_S<-cbind(out[,1],out_S)

intersec_final<-intersec[1:(num_clust/2),]
intersec_final$from_clust<-paste(intersec_final$clust,"f",sep="")
intersec_final$to_clust<-paste(intersec_final$clust,"t",sep="")

for(t in 0:sim_time){
  time_temp<-t
  out.file.name <- paste("/Users/JY/Dropbox/OD5_avg_relden_a1_b40_ana_t_", t, ".png", sep="")
  png(out.file.name, width=396, height=396)
  
  ################################################################################################
  simulation_ratio<-(clust_avg_final[clust_avg_final$time==time_temp,]$I/(clust_avg_final[clust_avg_final$time==time_temp,]$I+clust_avg_final[clust_avg_final$time==time_temp,]$S))
  analysis_ratio<-(as.vector(out_I[out_I[,1]==time_temp,2:(ncol(out_I))]))/(as.vector(out_I[out_I[,1]==time_temp,2:(ncol(out_I))])+as.vector(out_S[out_S[,1]==time_temp,2:(ncol(out_S))]))
  
  simulation_I<-(clust_avg_final[clust_avg_final$time==time_temp,]$I)
  analysis_I<-(as.vector(out_I[out_I[,1]==time_temp,2:(ncol(out_I))]*total_nodes))
  num_of_nodes_clust_sim<-(clust_avg_final[clust_avg_final$time==time_temp,]$I+clust_avg_final[clust_avg_final$time==time_temp,]$S)
  num_of_nodes_clust_ana<-(as.vector(out_I[out_I[,1]==time_temp,2:(ncol(out_I))]*total_nodes))+(as.vector(out_S[out_S[,1]==time_temp,2:(ncol(out_S))]*total_nodes))
  
  ft<-data.frame("from"=intersec_final$from_clust,"to"=intersec_final$to_clust, "simulation_ratio"=simulation_ratio,"analysis_ratio"=analysis_ratio, "simulation_I"=simulation_I,"analysis_I"=analysis_I, "num_of_nodes_sim"=num_of_nodes_clust_sim, "num_of_nodes_ana"=num_of_nodes_clust_ana, stringsAsFactors=FALSE)
  ft<-ft[order(ft$num_of_nodes_ana),]
  loc1 <- data.frame("clust" = intersec_final$from_clust, 
                     "x_val" = intersec_final$from_x, 
                     "y_val" = intersec_final$from_y) 
  loc2 <- data.frame("clust" = intersec_final$to_clust, 
                     "x_val" = intersec_final$to_x, 
                     "y_val" = intersec_final$to_y) 
  loc <- rbind(loc1,loc2)
  
  g_snap <- graph.data.frame(ft[,1:2], directed = TRUE, vertices = loc)
  lo <- as.matrix(loc[,2:3])
  
#  E(g_snap)$simulation_ratio <- ft$simulation_ratio
#  E(g_snap)$analysis_ratio <- ft$analysis_ratio
  E(g_snap)$simulation_I <- ft$simulation_I # the number of informed vehicles per road segment from Markov simulation results
  E(g_snap)$analysis_I <- ft$analysis_I # the number of informed vehicles per road segment from the set of differential equations
  E(g_snap)$num_of_nodes_sim <- ft$num_of_nodes_sim # traffic density per road segment from Markov simulation results
  E(g_snap)$num_of_nodes_ana <- ft$num_of_nodes_ana # traffic density per road segment from the set of differential equations
  
  fine <- 10
  colfunc <- colorRampPalette(c("gray80","black","blue","red"))
  
  palette = colfunc(fine)
  
  color_temp1<-ifelse(ceiling(E(g_snap)$simulation_I/100)==0,1,ceiling(E(g_snap)$simulation_I/100))
  color_temp2<-ifelse(ceiling(E(g_snap)$analysis_I/100)==0,1,ceiling(E(g_snap)$analysis_I/100))
  color_temp3<-ifelse(ceiling(E(g_snap)$num_of_nodes_sim/100)==0,1,ceiling(E(g_snap)$num_of_nodes_sim/100))
  color_temp4<-ifelse(ceiling(E(g_snap)$num_of_nodes_ana/100)==0,1,ceiling(E(g_snap)$num_of_nodes_ana/100))
  
  E(g_snap)$color1 = palette[color_temp1]
  E(g_snap)$color2 = palette[color_temp2]
  E(g_snap)$color3 = palette[color_temp3]
  E(g_snap)$color4 = palette[color_temp4]
  
  par(mar=c(0,0,0,0))
  # plot the number of informed vehicles per road segment from Markov simulation results: the thickness of the road segment is linearly proportional to the number of vehicles.
  plot.igraph(g_snap, layout = lo, edge.arrow.size=0, edge.width = E(g_snap)$simulation_I/10, edge.color="red", vertex.color=NA, vertex.frame.color = NA, vertex.size = 0., vertex.label= NA, edge.label=ifelse(round(E(g_snap)$simulation_I)>0,round(E(g_snap)$simulation_I),NA), edge.label.color="black", edge.label.cex=1.3, asp=0)
  # plot the number of informed vehicles per road segment from the set of differential equations: the thickness of the road segment is linearly proportional to the number of vehicles.
  plot.igraph(g_snap, layout = lo, edge.arrow.size=0, edge.width = E(g_snap)$analysis_I/10, edge.color="red", vertex.color=NA, vertex.frame.color = NA, vertex.size = 0., vertex.label= NA, edge.label=ifelse(round(E(g_snap)$analysis_I)>0,round(E(g_snap)$analysis_I),NA), edge.label.color="black", edge.label.cex=1.3, asp=0)
  # plot traffic density per road segment from Markov simulation results: the thickness of the road segment is linearly proportional to the number of vehicles.
  plot.igraph(g_snap, layout = lo, edge.arrow.size=0, edge.width = E(g_snap)$num_of_nodes_sim/10, edge.color="azure4", vertex.color=NA, vertex.frame.color = NA, vertex.size = 0., vertex.label= NA, edge.label=ifelse(round(E(g_snap)$num_of_nodes_sim)>0,round(E(g_snap)$num_of_nodes_sim),NA), edge.label.color="black", edge.label.cex=1.3, asp=0)
  # plot traffic density per road segment from the set of differential equations: the thickness of the road segment is linearly proportional to the number of vehicles.
  plot.igraph(g_snap, layout = lo, edge.arrow.size=0, edge.width = E(g_snap)$num_of_nodes_ana/10, edge.color="azure4", vertex.color=NA, vertex.frame.color = NA, vertex.size = 0., vertex.label= NA, edge.label=ifelse(round(E(g_snap)$num_of_nodes_ana)>0,round(E(g_snap)$num_of_nodes_ana),NA), edge.label.color="black", edge.label.cex=1.3, asp=0)

  dev.off()
}



