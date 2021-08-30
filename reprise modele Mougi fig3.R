## remake of Mougi's model

#### parameters definition ####
# li rate of nutrient loss from the system
# mi mortality rate
# I nutrient input rate
# am decomposition rate
# acg uptake rate of the producers by the consumers
# acb uptake rate of microbes by predators
# em conversion efficiency of detritus into microbes
# ecg conversion efficiency of producers into consumers
# ecb conversion efficiency of microbes into predators
# deltai fraction of nutrients released
# r nutrient uptake rate of the producers
# rm nutrient uptake rate of decomposers
# q ratio C:N ratio of detritus on C:N ration of decomposers

# em = acg = acb = e
# deltai = delta
# ln = ld = l


#require(tidyverse)
#require(rootSolve)
require(pracma)


#### biomass equilibrium ####
# C-limited case
biomass_equilibrium_C <- function(param){
  with(as.list(c(param)),{
    # plants biomass equilibrium
    Peq <- (lcg + mcg)/(e * acg)  
    
    # decomposers biomass equilibrium
    Meq <- (lcb + mcb)/(e * acb)
    
    # detritus biomass equilibrium
    b <- (1 - e)*(Peq*(-mp - lp) + Meq*(-mm - lm)) # part of Neq
    c <- am*Meq*(delta*(1 - e)*(e*q + 1) - e*(q - 1)) # part of Neq
    d <- r*Peq*(1 - delta*(1 - e)) + ln # part of Neq
    
    f <- ld + am*Meq*(1-(1 - delta)*(1 - e)*(e*q + 1)) - mcb*e*am*q/acb # part of Deq
    g <- (1 - delta)*(1 - e)*Peq*r + mcg*r/acg # part of Deq
    h <- mp*Peq + mm*Meq + mcg*(-mp -lp)/acg + mcb*(-mm - lm)/acb # part of Deq
    
    Deq <- ((delta*b + I)*g/d + h + (1 - delta)*b)/(f - c*g/d)
    
    # nutrients biomass equilibrium
    Neq <- (Deq*c + delta*b + I)/d
    
    # consumers biomass equilibrium
    Cgeq <- (r*Neq - lp - mp)/acg
  
    # predators biomass equilibrium
    Cbeq <- (q*e*am*Deq - lm - mm)/acb
    
    return(list(N=Neq,P=Peq,Cg=Cgeq,D=Deq,M=Meq,Cb=Cbeq))
  })
}

# N-limited case
biomass_equilibrium_N <- function(param){
  with(as.list(c(param)),{
    # plants biomass equilibrium
    Peq <- (lcg + mcg)/(e * acg)  
    
    # decomposers biomass equilibrium
    Meq <- (lcb + mcb)/(e * acb)
    
    # detritus biomass equilibrium
    b <- (1 - e)*(Peq*(-mp - lp) + Meq*(-mm - lm)) # part of Neq
    k <- delta*(1 - e)*am*Meq*(e + 1) # part of Neq
    n <- ln + rm*Meq + r*Peq - delta*(1 - e)*(r*Peq + rm*Meq) # part of Neq
    
    m <- mp*Peq + mm*Meq # part of Deq
    u <- mcb*(-mm - lm)/acb + mcg*(-mp - lp)/acg # part of Deq
    t <- (1 - delta)*(1 - e)*(Peq*r + rm*Meq) + mcb*rm/acb + mcg*r/acg # part of Deq
    d <- ld + am*Meq*(1 - (1 - delta)*(1 - e)*(1 + e)) - mcb*e*am/acb # part of Deq
    
    Deq <- (t*(I + delta*b)/n + (1 - delta)*b + m + u)/(d - k*t/n)
    
    # nutrients biomass equilibrium
    Neq <- (Deq*k + I + delta*b)/n
    
    # consumers biomass equilibrium
    Cgeq <- (r*Neq - lp - mp)/acg
    
    # predators biomass equilibrium
    Cbeq <- (e*am*Deq + rm*Neq - lm - mm)/acb
    
    return(list(N=Neq,P=Peq,Cg=Cgeq,D=Deq,M=Meq,Cb=Cbeq))
  })
}


#### jacobian matrix ####
# C-limited case
jacobian_C <- function(B,param){
  with(as.list(c(B,param)),{
    
    ## N partial derivatives
    d1N <- -r*P -ln  # for the first equation
    d2N <- r*P       # second equation
    d3N <- 0         # 3 equation
    d4N <- 0         # 4 equation
    d5N <- 0         # 5 equation
    d6N <- 0         # 6 equation
    dN <- c(d1N,d2N,d3N,d4N,d5N,d6N)
    
    ## P partial derivatives
    d1P <- delta*(1 - e)*acg*Cg - r*N
    d2P <- r*N - mp - lp -acg*Cg
    d3P <- e*acg*Cg
    d4P <- (1 - delta)*(1 - e)*acg*Cg + mp
    d5P <- 0
    d6P <- 0
    dP <- c(d1P,d2P,d3P,d4P,d5P,d6P)
    
    ## Cg partial derivatives
    d1Cg <- delta*(1 - e)*acg*P
    d2Cg <- -acg*P
    d3Cg <- e*acg*P - mcg - lcg
    d4Cg <- (1 - delta)*(1 - e)*acg*P + mcg
    d5Cg <- 0
    d6Cg <- 0
    dCg <- c(d1Cg,d2Cg,d3Cg,d4Cg,d5Cg,d6Cg)
    
    ## D partial derivatives
    d1D <- delta*(1 - e)*am*M - e*am*M*(q - 1)
    d2D <- 0
    d3D <- 0
    d4D <- (1 - delta)*(1 - e)*am*M - ld - am*M
    d5D <- e*am*M*q
    d6D <- 0
    dD <- c(d1D,d2D,d3D,d4D,d5D,d6D)
    
    ## M partial derivatives
    d1M <- delta*(1 - e)*am*D - e*am*D*(q - 1) + delta*(1 - e)*acb*Cb
    d2M <- 0
    d3M <- 0
    d4M <- (1 - delta)*(1 - e)*acb*Cb + (1 - delta)*(1 - e)*am*D + mm - am*D
    d5M <- e*am*D*q - mm - lm - acb*Cb
    d6M <- e*acb*Cb
    dM <- c(d1M,d2M,d3M,d4M,d5M,d6M)
    
    ## Cb partial derivatives
    d1Cb <- delta*(1 - e)*acb*M
    d2Cb <- 0
    d3Cb <- 0
    d4Cb <- (1 - delta)*(1 - e)*acb*M + mcb
    d5Cb <- -acb*M
    d6Cb <- e*acb*M - mcb - lcb
    dCb <- c(d1Cb,d2Cb,d3Cb,d4Cb,d5Cb,d6Cb)
    
    ## jacobian
    Jac <- matrix(c(dN,dP,dCg,dD,dM,dCb), ncol = 6, byrow=FALSE)
    
    return(Jac)
  })
}

# N-limited case
jacobian_N <- function(B,param){
  with(as.list(c(B,param)),{
    
    ## N partial derivatives
    d1N <- -r*P - rm*M -ln
    d2N <- r*P
    d3N <- 0
    d4N <- 0
    d5N <- rm*M
    d6N <- 0
    dN <- c(d1N,d2N,d3N,d4N,d5N,d6N)
    
    ## P partial derivatives
    d1P <- delta*(1 - e)*acg*Cg - r*N
    d2P <- r*N - mp - lp -acg*Cg
    d3P <- e*acg*Cg
    d4P <- (1 - delta)*(1 - e)*acg*Cg + mp
    d5P <- 0
    d6P <- 0
    dP <- c(d1P,d2P,d3P,d4P,d5P,d6P)
    
    ## Cg partial derivatives
    d1Cg <- delta*(1 - e)*acg*P
    d2Cg <- -acg*P
    d3Cg <- e*acg*P - mcg - lcg
    d4Cg <- (1 - delta)*(1 - e)*acg*P + mcg
    d5Cg <- 0
    d6Cg <- 0
    dCg <- c(d1Cg,d2Cg,d3Cg,d4Cg,d5Cg,d6Cg)
    
    ## D partial derivatives
    d1D <- delta*(1 - e)*am*M
    d2D <- 0
    d3D <- 0
    d4D <- (1 - delta)*(1 - e)*am*M - ld - am*M
    d5D <- e*am*M
    d6D <- 0
    dD <- c(d1D,d2D,d3D,d4D,d5D,d6D)
    
    ## M partial derivatives
    d1M <- delta*(1 - e)*(acb*Cb + am*D) - rm*N
    d2M <- 0
    d3M <- 0
    d4M <- (1 - delta)*(1 - e)*(acb*Cb + am*D) + mm - am*D
    d5M <- e*am*D + rm*N - mm - lm - acb*Cb
    d6M <- e*acb*Cb
    dM <- c(d1M,d2M,d3M,d4M,d5M,d6M)
    
    ## Cb partial derivatives
    d1Cb <- delta*(1 - e)*acb*M
    d2Cb <- 0
    d3Cb <- 0
    d4Cb <- (1 - delta)*(1 - e)*acb*M + mcb
    d5Cb <- -acb*M
    d6Cb <- e*acb*M - mcb - lcb
    dCb <- c(d1Cb,d2Cb,d3Cb,d4Cb,d5Cb,d6Cb)
    
    ## jacobian
    Jac <- matrix(c(dN,dP,dCg,dD,dM,dCb), ncol = 6, byrow=FALSE)
    
    return(Jac)
  })
} 




#### function ####
model <- function(param){
  if(param$limitation=="C-limited"){
    biomass<-biomass_equilibrium_C(param) # biomass at equilibrium
    Jac <- jacobian_C(unlist(biomass),param) # jacobian matrix
  }
  if(param$limitation=="N-limited"){
    biomass<-biomass_equilibrium_N(param) # biomass at equilibrium
    Jac <- jacobian_N(unlist(biomass),param) # jacobian matrix
  }
  
  # matrix Ve disturbance variances
  Ve <- diag(x=c(zetaI=0.0100,zetaN=0.0025,zetaP=0.0025,zetaCg=0.0025,zetaD=0.0025,zetaM=0.0025,zetaCb=0.0025))

  # identity matrix
  In <- diag(6)

  #### distrubance ###
  M <- diag(x=c(biomass$N**param$z,biomass$P**param$z,biomass$Cg**param$z,biomass$D**param$z,biomass$M**param$z,biomass$Cb**param$z)) # avec C_limited[] qui represente les valeurs a lequilibre
  M <- cbind(c(param$I**param$z,zeros(5,1)),M)

  #### Lyapunov equation ####
  MVeM <- M%*%Ve%*%t(M)
  dim(MVeM) <- c(6*6,1)
  C <- -solve(Jac%x%In + In%x%Jac) %*% MVeM
  dim(C) <- c(6,6)

  #### invariability ####
  somme_var <- 0
  somme_Cov <- 0
  for (k in 1:6){
    somme_Cov <- somme_Cov + sum(C[k,]) - C[k,k] # sum of covariances
    somme_var <- somme_var + C[k,k]  # sum of variances
  }
  invariability <- sum(unlist(biomass))/(somme_var + 2*somme_Cov) # invariability total
  
  inv_species <- list(N=biomass$N/sqrt(C[1,1]),
                      P=biomass$P/sqrt(C[2,2]),
                      Cg=biomass$Cg/sqrt(C[3,3]),
                      D=biomass$D/sqrt(C[4,4]),
                      M=biomass$M/sqrt(C[5,5]),
                      Cb=biomass$Cb/sqrt(C[6,6]))
  
  variance <- list(N=C[1,1],
                   P=C[2,2],
                   Cg=C[3,3],
                   D=C[4,4],
                   M=C[5,5],
                   Cb=C[6,6])

  #### Resilience ####
  eigenvalues <- eigen(Jac, only.values = TRUE)
  Resilience <- abs(max(Re(eigenvalues$values))) # cest la valeur absolue de la valeur propre ayant la partie reelle la plus grande

  results <- list(Biomass=biomass,
                  Invariability=invariability,
                  Resilience=Resilience,
                  Inv_species=inv_species,
                  Variance=variance)
  return(results)
}


####################
##### figure 3 #####
####################

dim=6 # dimension of the system (number of compartments)
parametres=list(I=2,
                r=2,
                rm=1,
                mp=0.1,
                mcg=0.1,
                mm=0.1,
                mcb=0.1,
                e=0.25,
                acg=1,
                acb=1,
                delta=0.5,
                lp=0.1,
                lcg=0.1,
                lm=0.1,
                lcb=0.1,
                q=1.2)
am <- seq(1,50,0.1)
ln=c(0.001,0.07,0.1,0.2,0.5,1.0,2.0) ; ld=ln
limitation=c("C-limited","N-limited")
z <- c(0,0.5,1)

# data frame creation with all parameters
param_data<-expand.grid(parametres=list(parametres),
                         am=am,
                         ln=ln,
                         limitation=limitation,
                         z=z)
param_data$simulationID<-1:dim(param_data)[1]

# results storage
resilience<-param_data
resilience$resilience=0

biomass<-param_data
biomass<-cbind(biomass,matrix(0,dim(biomass)[1],dim))
names(biomass)[(dim(param_data)[2]+1):dim(biomass)[2]]<-c("N","P","Cg","D","M","Cb")

invariability <- param_data
invariability$invariability=0

inv_species <- param_data
inv_species <-cbind(inv_species,matrix(0,dim(inv_species)[1],dim))
names(inv_species)[(dim(param_data)[2]+1):dim(inv_species)[2]]<-c("N","P","Cg","D","M","Cb")

variance <- param_data
variance <-cbind(variance,matrix(0,dim(variance)[1],dim))
names(variance)[(dim(param_data)[2]+1):dim(variance)[2]]<-c("N","P","Cg","D","M","Cb")

# results
for (i in 1:dim(param_data)[1]){
  param<-param_data$parametres[[i]]
  param$am=param_data$am[i]
  param$ln=param_data$ln[i]
  param$ld=param_data$ln[i]
  param$limitation=param_data$limitation[i]
  param$z=param_data$z[i]
  results<-model(param)
  resilience$resilience[i]=results$Resilience
  invariability$invariability[i]=results$Invariability
  biomass[i,(dim(param_data)[2]+1):dim(biomass)[2]]=unlist(results$Biomass)
  inv_species[i,(dim(param_data)[2]+1):dim(inv_species)[2]]=unlist(results$Inv_species)
  variance[i,(dim(param_data)[2]+1):dim(variance)[2]]=unlist(results$Variance)
}


# results recording
setwd("~/M2/Stage/tableaux données/fig3")
resilienc <- resilience[,2:7]
write.csv2(resilienc, file="resiliencefig3.csv",row.names = FALSE)
invariabilit <- invariability[,2:7]
write.csv2(invariabilit, file="invariabilityfig3.csv",row.names = FALSE)
biomas <- biomass[,2:12]
write.csv2(biomas, file="biomassfig3.csv",row.names = FALSE)
inv_specie <- inv_species[,2:12]
write.csv2(inv_specie, file="inv_speciesfig3.csv",row.names = FALSE)
varianc <- variance[,2:12]
write.csv2(varianc, file="variance.csv",row.names = FALSE)

require('ggplot2')
require('reshape2')
library(scales)
library(cowplot)


#################
#### figures ####

#### resilience ####
# z=0
ggplot(data = resilience[resilience$z == 0,]) +
  geom_line(aes(am,resilience,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

# z=0.5
ggplot(data = resilience[resilience$z == 0.5,]) +
  geom_line(aes(am,resilience,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

# z=1
graph_res <- ggplot(data = resilience[resilience$z == 1,]) +
  geom_line(aes(am,resilience,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")
ggsave("fig3_resilience.pdf",graph_res,width = 8,height = 5)


#### biomass ####
biomass<-melt(biomass,
              id.vars = names(param_data),
              variable.name = "species",
              value.name = "biomass")

# z=0
ggplot(data = biomass[biomass$z == 0,]) +
  geom_line(aes(am,biomass,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")

# z=0.5
ggplot(data = biomass[biomass$z == 0.5,]) +
  geom_line(aes(am,biomass,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")

# z=1
graph_biomass <- ggplot(data = biomass[biomass$z == 1,]) +
  geom_line(aes(am,biomass,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")+
  scale_y_continuous(limits = c(-0.4,0.1))

ggsave("fig3_biomass.pdf",graph_biomass,width = 8,height = 5)


#### invariability ####
inv_species<-melt(inv_species,
                  id.vars = names(param_data),
                  variable.name = "species",
                  value.name = "invariability")
III <- invariability[invariability$z==0.5 & invariability$limitation=="C-limited",]
max(III$invariability[III$ln==1])
min(III$invariability[III$ln==1])
which.max(III$invariability[III$ln==2])
III[14,]
min(III$invariability[III$ln==2])
# z=0
pi1 <- ggplot(data = invariability[invariability$z == 0,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  scale_y_continuous(limits = c(-10,300))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")

pisp1 <- ggplot(data = inv_species[inv_species$z == 0,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  scale_y_continuous(limits = c(0,70))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")
ggsave("fig3_inv_sp1_cut.pdf",pisp1,width = 9,height = 6)

# z=0.5
pi2 <- ggplot(data = invariability[invariability$z == 0.5,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  scale_y_continuous(limits = c(0,200))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")

pisp2 <- ggplot(data = inv_species[inv_species$z == 0.5,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")

# z=1
pi3 <- ggplot(data = invariability[invariability$z == 1,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_wrap(~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  #scale_y_continuous(limits = c(0,225))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")

pisp3 <- ggplot(data = inv_species[inv_species$z == 1,]) +
  geom_line(aes(am,invariability,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")+
  theme(legend.position = "none")

legend <- get_legend(pi3)
# graphs grouping
graph_inv <- ggdraw(xlim=c(0,2),ylim=c(0,2))+
  draw_plot(pi1,0,1,1,1)+
  draw_plot(pi2,1,1,1,1)+
  draw_plot(pi3,0,0,1,1)+
  draw_plot(legend, 1,0,1,1)+
  draw_plot_label(c("A","B","C"),c(0,1,0),c(2,2,1),size=10)
setwd("~/M2/Stage/graphes/fig3")
ggsave("fig3_invariability.pdf",graph_inv,width = 9,height = 6)

graph_inv_sp <- ggdraw(xlim=c(0,4),ylim=c(0,1))+
  draw_plot(pisp1,0,0,1,1)+
  draw_plot(pisp2,1,0,1,1)+
  draw_plot(pisp3,2,0,1,1)+
  draw_plot(legend,3,0,1,1)+
  draw_plot_label(c("A","B","C"),c(0,1,2),c(1,1,1),size=12)
ggsave("fig3_inv_sp_cut.pdf",graph_inv_sp,width = 9,height = 6)
graph_inv_sp1 <- ggdraw(xlim=c(0,2),ylim=c(0,1))+
  draw_plot(pisp2,0,0,1,1)+
  draw_plot(pisp3,1,0,1,1)+
  draw_plot_label(c("B","C"),c(0,1),c(1,1),size=10)
ggsave("fig3_inv_sp2.pdf",graph_inv_sp1,width = 9,height = 6)

#### variance ####
variance<-melt(variance,
                  id.vars = names(param_data),
                  variable.name = "species",
                  value.name = "variance")

# z=0
pv1 <- ggplot(data = variance[variance$z == 0,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  scale_y_continuous(limits = c(0,0.5))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")+
  theme(legend.position = "none")
ggsave("fig3_variance1_cut.pdf",pv1,width = 9,height = 6)

# z=0.5
pv2 <- ggplot(data = variance[variance$z == 0.5,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")+
  theme(legend.position = "none")

# z=1
pv3 <- ggplot(data = variance[variance$z == 1,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")+
  theme(legend.position = "none")

graph_var <- ggdraw(xlim=c(0,4),ylim=c(0,1))+
  draw_plot(pv1,0,0,1,1)+
  draw_plot(pv2,1,0,1,1)+
  draw_plot(pv3,2,0,1,1)+
  draw_plot(legend,3,0,1,1)+
  draw_plot_label(c("A","B","C"),c(0,1,2),c(1,1,1),size=10)
ggsave("fig3_variance_cut.pdf",graph_var,width = 9,height = 6)

#### variance + biomass ####
# z=0
pvb1 <- ggplot(data = variance[variance$z == 0,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  geom_line(data=biomass[biomass$z == 0,], aes(am,biomass,colour=as.factor(ln)),size=1, linetype=4) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  scale_y_continuous(limits = c(0,15))+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance + biomass")+
  theme(legend.position = "none")
ggsave("fig3_variance_biomass1_cut.pdf",pvb1,width = 9,height = 6)

# z=0.5
pvb2 <- ggplot(data = variance[variance$z == 0.5,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  geom_line(data=biomass[biomass$z == 0.5,], aes(am,biomass,colour=as.factor(ln)),size=1, linetype=4) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance + biomass")+
  theme(legend.position = "none")

# z=1
pvb3 <- ggplot(data = variance[variance$z == 1,]) +
  geom_line(aes(am,variance,colour=as.factor(ln)),size=1) +
  geom_line(data=biomass[biomass$z == 1,], aes(am,biomass,colour=as.factor(ln)),size=1, linetype=4) +
  facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Openness")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance + biomass")+
  theme(legend.position = "none")

graph_var_b <- ggdraw(xlim=c(0,4),ylim=c(0,1))+
  draw_plot(pvb1,0,0,1,1)+
  draw_plot(pvb2,1,0,1,1)+
  draw_plot(pvb3,2,0,1,1)+
  draw_plot(legend,3,0,1,1)+
  draw_plot_label(c("A","B","C"),c(0,1,2),c(1,1,1),size=10)
ggsave("fig3_variance_biomass_cut.pdf",graph_var_b,width = 9,height = 6)


#### comparaison Haegeman ####
bio <- biomass[biomass$ln == 0.5,]

bio<-melt(bio,
              id.vars = names(param_data),
              variable.name = "species",
              value.name = "biomass")

b1 <- ggplot(data = bio[bio$z == 1 & bio$limitation == "C-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")

##
var <- variance[variance$ln == 0.5,]

var<-melt(var,
               id.vars = names(param_data),
               variable.name = "species",
               value.name = "variance")

v1 <- ggplot(data = var[var$z == 1 & var$limitation == "C-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

##
res <- resilience[resilience$ln == 0.5,]

r1 <- ggplot(data = res[res$z == 1 & res$limitation == "C-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

##
inv <- invariability[invariability$ln==0.5,]

i1 <- ggplot(data = inv[inv$z == 1 & inv$limitation == "C-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

invsp <- inv_species[inv_species$ln ==0.5,]

invsp<-melt(invsp,
                  id.vars = names(param_data),
                  variable.name = "species",
                  value.name = "invariability")

is1 <- ggplot(data = invsp[invsp$z == 1 & invsp$limitation == "C-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graph <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b1,0,1,1,1)+
  draw_plot(v1,1,1,1,1)+
  draw_plot(r1,0,0,1,1)+
  draw_plot(i1,1,0,1,1)+
  draw_plot(is1,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z1.pdf",graph,width = 9,height = 6)

##
b1n <- ggplot(data = bio[bio$z == 1 & bio$limitation == "N-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")


v1n <- ggplot(data = var[var$z == 1 & var$limitation == "N-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

r1n <- ggplot(data = res[res$z == 1 & res$limitation == "N-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

i1n <- ggplot(data = inv[inv$z == 1 & inv$limitation == "N-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")


is1n <- ggplot(data = invsp[invsp$z == 1 & invsp$limitation == "N-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graphn <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b1n,0,1,1,1)+
  draw_plot(v1n,1,1,1,1)+
  draw_plot(r1n,0,0,1,1)+
  draw_plot(i1n,1,0,1,1)+
  draw_plot(is1n,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z1N.pdf",graphn,width = 9,height = 6)

##
##
b0n <- ggplot(data = bio[bio$z == 0 & bio$limitation == "N-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")


v0n <- ggplot(data = var[var$z == 0 & var$limitation == "N-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

r0n <- ggplot(data = res[res$z == 0 & res$limitation == "N-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

i0n <- ggplot(data = inv[inv$z == 0 & inv$limitation == "N-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")


is0n <- ggplot(data = invsp[invsp$z == 0 & invsp$limitation == "N-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graph0n <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b0n,0,1,1,1)+
  draw_plot(v0n,1,1,1,1)+
  draw_plot(r0n,0,0,1,1)+
  draw_plot(i0n,1,0,1,1)+
  draw_plot(is0n,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z0N.pdf",graph0n,width = 9,height = 6)

b0c <- ggplot(data = bio[bio$z == 0 & bio$limitation == "C-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")


v0c <- ggplot(data = var[var$z == 0 & var$limitation == "C-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

r0c <- ggplot(data = res[res$z == 0 & res$limitation == "C-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

i0c <- ggplot(data = inv[inv$z == 0 & inv$limitation == "C-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")


is0c <- ggplot(data = invsp[invsp$z == 0 & invsp$limitation == "C-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graph0c <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b0c,0,1,1,1)+
  draw_plot(v0c,1,1,1,1)+
  draw_plot(r0c,0,0,1,1)+
  draw_plot(i0c,1,0,1,1)+
  draw_plot(is0c,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z0C.pdf",graph0c,width = 9,height = 6)

##
##
b05n <- ggplot(data = bio[bio$z == 0.5 & bio$limitation == "N-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")


v05n <- ggplot(data = var[var$z == 0.5 & var$limitation == "N-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

r05n <- ggplot(data = res[res$z == 0.5 & res$limitation == "N-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

i05n <- ggplot(data = inv[inv$z == 0.5 & inv$limitation == "N-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")


is05n <- ggplot(data = invsp[invsp$z == 0.5 & invsp$limitation == "N-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graph05n <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b05n,0,1,1,1)+
  draw_plot(v05n,1,1,1,1)+
  draw_plot(r05n,0,0,1,1)+
  draw_plot(i05n,1,0,1,1)+
  draw_plot(is05n,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z05N.pdf",graph05n,width = 9,height = 6)

b05c <- ggplot(data = bio[bio$z == 0.5 & bio$limitation == "C-limited",]) +
  geom_line(aes(am,biomass,colour=species),size=1) +
  #facet_grid(species~limitation)+
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Biomass")


v05c <- ggplot(data = var[var$z == 0.5 & var$limitation == "C-limited",]) +
  geom_line(aes(am,variance,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Variance")

r05c <- ggplot(data = res[res$z == 0.5 & res$limitation == "C-limited",]) +
  geom_line(aes(am,resilience),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Resilience")

i05c <- ggplot(data = inv[inv$z == 0.5 & inv$limitation == "C-limited",]) +
  geom_line(aes(am,invariability),colour="deeppink",size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")


is05c <- ggplot(data = invsp[invsp$z == 0.5 & invsp$limitation == "C-limited",]) +
  geom_line(aes(am,invariability,colour=species),size=1) +
  scale_colour_manual(values = rainbow(7),
                      name = "Species")+
  theme_minimal()+
  xlab("Decomposition rate")+
  ylab("Invariability")

graph05c <- ggdraw(xlim=c(0,3),ylim=c(0,2))+
  draw_plot(b05c,0,1,1,1)+
  draw_plot(v05c,1,1,1,1)+
  draw_plot(r05c,0,0,1,1)+
  draw_plot(i05c,1,0,1,1)+
  draw_plot(is05c,2,1,1,1)+
  draw_plot_label(c("(a)","(b)","(c)","(d)","(e)"),c(0,1,2,0,1),c(2,2,2,1,1),size=10)
ggsave("fig3_bio_var_res_inv_z05C.pdf",graph05c,width = 9,height = 6)

