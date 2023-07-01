https://github.com/MeriSoto/SOTO_2023_MODELLING_DRIVERS_OF_DISCARDS_THROUGH_BAYESIAN_SPATIO_TEMPORAL_MODELS_IN-TRAWLING_FI-

################################################################################
############################# SETTING THE SCENE ################################
################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the R script used for the data analysis presented in the manuscript

# Modelling drivers of trawl fisheries discards through Bayesian spatio-temporal models 

# Soto, M. , Fernández-Peralta, L., Rey, J., Czerwisnki, I., García-Cancela, R., 
# Llope, M., Cabrera-Busto, J.,  Liébana, M. and Pennino, M. G.

# FISHERIES RESEARCH (submitted 17th April 2023)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOAD FILES AND PACKAGES ~~~~~~~~~~~~~~~~~~~~~~~~

# CLEAR SCREEN AND SET WORKING DIRECTORY 

ls()
rm(list=ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################################################
################################# DATA ANALYSIS ################################
################################################################################

# LOAD THE DATA .RData 

load("SOTO_DATA.RData")
df <- SOTO_DATA

# LOAD PACKAGES ----

library(sp)
library(INLA)
library(rgeos)
library(raster)
library(fields)
library(ggplot2)
library(ggspatial)
library("rnaturalearth")


# STUDY AREA OF MAURITANIA  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mau <- getData('GADM',country="MRT",level=0)

ext<-extent(-19.01852, -15.98148, 15.98154, 21.34012) 
Mau <- crop(Mau, ext) 

xym<- as.matrix(data.frame(x =c(-19.01852, -15.98148, -15.98148,-19.01852),
                           y = c(15.98154,15.98154,21.34012,21.34012)))

p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

Mau_rec<-crop(Mau, sps)

proj4string(sps)<-proj4string(Mau_rec)

p2 = Polygon(xym)
ps2 = Polygons(list(p2),1)
sps2 = SpatialPolygons(list(ps2))

coast <- gDifference(sps2, Mau_rec)

Loc <- cbind(df$lon,df$lat)
D <- dist(Loc)

bound=bound<-inla.sp2segment(coast)


# MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mesh <- inla.mesh.2d(loc=as.matrix(Loc),boundary = bound,
                       offset = c(0.5, 0.5),
                       cutoff = 0.05, max.edge = c(0.17, 0.4))

# STACK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

datacolumns <- c("response1","response2","response3","lon","lat")
data <- df[,datacolumns]

# SPDE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spde  <- inla.spde2.pcmatern(mesh = mesh, alpha = 1.5,
                             prior.range = c(2, 0.6),
                             prior.sigma = c(20, 0.01))

# A MATRIX  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A.est <- inla.spde.make.A(mesh, loc=cbind(df$lon, df$lat)) 

# STACK FOR TD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stk.est1 <- inla.stack(data=list(y=data$response1),
                     A=list(A.est, 1),
                     effects=list(spatial=1:spde$n.spde,
                                  data.frame(beta0=1, df)),
                     tag="est",compress = TRUE, remove.unused = TRUE)

# STACK FOR DPUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stk.est2 <- inla.stack(data=list(y=data$response2),
                     A=list(A.est, 1),
                     effects=list(spatial=1:spde$n.spde,
                                  data.frame(beta0=1, df)),
                     tag="est",compress = TRUE, remove.unused = TRUE)
# STACK FOR TDR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stk.est3 <- inla.stack(data=list(y=data$response3),
                     A=list(A.est, 1),
                     effects=list(spatial=1:spde$n.spde,
                                  data.frame(beta0=1, df)),
                     tag="est",compress = TRUE, remove.unused = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMULA FOR INITIAL SPATIAL MODELS WITH ALL COVARIATES 
# CHOOSEN FROM RANDOM FOREST SELECTION
# .std = variable standardized
# skill =skipper experience

f1 <- y~-1 + beta0 + 
  offset(trawling.std) + 
  f(spatial, model=spde) +
  f(depth.std,model="rw1") + 
  f(year.std,model = "rw1") + 
  mixed.std + 
  hake.std + 
  seasonday.std + 
  strategy + 
  haul_order.std + 
  skill + 
  SST.std + 
  vessel_size.std + 
  lunar_phase.std + 
  f(vessel, model ='iid') 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INLA INITIAL MODELS ALL COVARIATES

TDMODEL <- inla(f1, 
              data=inla.stack.data(stk.est1), family="gamma",
              control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
              control.predictor=list(A=inla.stack.A(stk.est1), compute=TRUE, 
                                     quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
              num.threads = 3,
              verbose=T)


DPUEMODEL <- inla(f1, 
               data=inla.stack.data(stk.est2), family="gamma",
               control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
               control.predictor=list(A=inla.stack.A(stk.est2), compute=TRUE, 
                                      quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
               num.threads = 3,
               verbose=T)

TDRMODEL <- inla(f1, 
                 data=inla.stack.data(stk.est3), family="gamma",
                 control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
                 control.predictor=list(A=inla.stack.A(stk.est3), compute=TRUE, 
                                        quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
                 num.threads = 3,
                 verbose=TRUE)




