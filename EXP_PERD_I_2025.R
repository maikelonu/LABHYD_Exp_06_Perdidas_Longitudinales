# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: PERDIDAS LONGITUDINALES + METODOS NUMERICOS

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////
# INFO:
# Bloomfield, V.A. Using R for Numerical Analysis in Science
# and Engineering. Taylor & Francis. New York
# //////////////////////////////////////////////////////////////////////////////////
# Custom_Functions_R (simples + anidadas)
# For-Loop-R (vectores y listas)
# Resolucion_ecuaciones_No_lineales_Unidimensional
# Depuraci?n y control de outliers
# Optimizaci?n_Inversa_Unidimensional(Colebrook + Swamme)
# An?lisis_grafico_ggplot2"
# //////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# Working directory is selected
setwd("~/Downloads/LABHYD_Exp_06_Perdidas_Longitudinales-master")

# CRAN libraries are loaded
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Custom functions 
# ////////////////////////////////////////////////////////

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Custom function cole.white (Colebrook White) is created
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cole.white <- function(f, ee, DD, RRe) {
  ((1/(sqrt(f))) + 2*(log10 ((ee / (3.7 * DD)) + (2.51/(RRe*sqrt(f))))))
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Custom function f.teo is created
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
f.teo <- function(ee, DD, RRe) {
  
  # list type container object is created
  L.Container <- list()
  root <- vector() # root
  f_root <- vector() # value of the function evaluated at that point
  iter <- vector() # number of iterations used 
  init_it <- vector() # NA
  estim_prec <- vector() # approximate estimated precision for root
  
  # list container is filled up using function uniroot{stats}
  for (i in 1:length(RRe)) {
    L.Container[[i]] <- uniroot(cole.white,c(0.001,1), ee=ee, DD=DD, RRe=RRe[i])
    root[i] <- unlist(L.Container[[i]][1])
    f_root[i] <- unlist(L.Container[[i]][2])
    iter[i] <- unlist(L.Container[[i]][3])
    init_it[i] <- unlist(L.Container[[i]][4])
    estim_prec[i] <- unlist(L.Container[[i]][5])
  }
  
  # list object is returned as vector
  return(data.frame(root,f_root,iter,init_it,estim_prec))
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Custom function f.opti is created
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
f.opti <- function(par, DD, RRe, fexp) { # par f(e)
  
  # list type container object is created
  L.Container <- list()
  root <- vector() # root
  
  # list container is filled up using function uniroot{stats}
  for (i in 1:length(RRe)) {
    L.Container[[i]] <- uniroot(cole.white,c(0.001,1), e=par, DD=DD, RRe=RRe[i])
    root[i] <- unlist(L.Container[[i]][1])
  }
  
  # Weighted sum of squares residuals (WSSR) Objective Function (OF) is returned
  return(sum((fexp - root))^2)
  
}
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ////////////////////////////////////////////////////////
# BLOCK: Data input
# ////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
df.obs <- read.table("Observaciones_Test_2025.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.obs, plotit = TRUE)

# names {base} function is requested
names(df.obs)

# ////////////////////////////////////////////////////////
# BLOCK: Problem solution
# ////////////////////////////////////////////////////////

# A theoretical Reynolds repetition is generated
Re <- c(1:10 %o% 10 ^ (3 : 7))

# -------------------------------------------------------
# Parameter values are defined
# -------------------------------------------------------

# Absolute rugosity [m]
e <- 0.0015/1000 # smooth pipe material

# Internal diameters [m]
D1 <- 20.13/1000 # Nominal 18 mm copper
D2 <- 18.20/1000 # Nominal 12 mm pvc
D3 <- 23.53/1000 # Nominal 18 mm pvc

# Names on df.obs$pipe variable are requested
summary(df.obs$pipe)

# Theoretical f values are generated for each pipe based on Swamee & Jain equation
cup18.SJ <- (0.25 / ((log10 ((e / (3.7 * D1)) + (5.74 / (Re ^ 0.9)))) ^ 2))
pvc12.SJ <- (0.25 / ((log10 ((e / (3.7 * D2)) + (5.74 / (Re ^ 0.9)))) ^ 2))
pvc18.SJ <- (0.25 / ((log10 ((e / (3.7 * D3)) + (5.74 / (Re ^ 0.9)))) ^ 2))

# Theoretical f values are generated for each pipe based on Colebrook White equation
df.cup18.CW <- f.teo(ee=e, DD=D1, RRe=Re)
df.pvc12.CW <- f.teo(ee=e, DD=D2, RRe=Re)
df.pvc18.CW <- f.teo(ee=e, DD=D3, RRe=Re)

# Vector values are extracted from data.frames
cup18.CW <- df.cup18.CW[,"root"]
pvc12.CW <- df.pvc12.CW[,"root"]
pvc18.CW <- df.pvc18.CW[,"root"]

# Individual data.frames are crested for each pipe
df.cup18 <- data.frame(cup18.SJ, cup18.CW, rep(c("copper18"), 50*2))
df.pvc12 <- data.frame(pvc12.SJ, pvc12.CW, rep(c("pvc12"), 50*2))                 
df.pvc18 <- data.frame(pvc18.SJ, pvc18.CW, rep(c("pvc18"), 50*2))                      

# data.frames names are normalized
names(df.cup18) <- (c("SJ","CW", "pipe"))
names(df.pvc12) <- (c("SJ","CW", "pipe"))
names(df.pvc18) <- (c("SJ","CW", "pipe"))

# df.factor data.frames is created using function rbind (vertical)
df.factor <- rbind(df.cup18,
                   df.pvc12,
                   df.pvc18)                     

# A ggplot2 object is created
fg01 <- ggplot() +
  geom_point(aes(x = CW,y = SJ,shape = pipe,colour = pipe),data=df.factor,size = 4.0) +
  geom_abline(data=df.cup18,linewidth = 0.95,linetype = 2) +
  facet_wrap(facets = ~pipe, nrow = 2, ncol = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(min.n = 5.0),limits = c(0.01,0.07)) +
  scale_y_continuous(breaks = scales::pretty_breaks(min.n = 5.0),limits = c(0.01,0.07)) +
  ggtitle("Comparacion de factores de friccion. Swamee-Jain vs Colebrook White)") +
  xlab("Factor f Colebrook White (-) ") +
  ylab("Factor f Swamee-Jain (-)") +
  theme_bw(base_size = 16.0) 

# A ggplot object is requested
fg01

# Correlation tests are requested                
cor.test(cup18.SJ, cup18.CW)
cor.test(pvc12.SJ, pvc12.CW)
cor.test(pvc18.SJ, pvc18.CW)

# a data.frame is created
df.base <- data.frame(cup18.CW,
                      pvc12.CW,
                      pvc18.CW)

# a data.frame is created and Reynolds theoretical values are added
   df.melted <- melt(df.base)
         Re2 <- rep(Re,3)
df.melted$Re <- Re2

# A ggplot object is created
fg02 <- ggplot() +
  geom_line(aes(x = Re,y = value,linetype = variable),data=df.melted,size = 0.85) +
  geom_point(aes(x = Re,y = f,shape = feeding,colour = pipe),data=df.obs,size = 5.0) +
  scale_x_continuous(trans="log10", limits = c(1e3,1e8),
                     breaks = c(c(1e2,1e3,5e3,1e4,2e4,3e4,5e4,1e5,5e5,1e6,1e7,1e8,1e9,1e10))) +
  scale_y_continuous(trans="log10", 
                     breaks = c(c(0.008, 0.010,0.012,0.014,0.016,0.018,0.020, 0.025,0.030,0.035,
                                  0.040,0.045,0.050,0.055,0.060,0.07, 0.1))) +
  geom_text(aes(x = Re,y = f,label = q_l_s),size = 3.5,data=df.obs,parse = FALSE) +
  ggtitle("Diagrama de Moody Teorico (Colebrook White)") +
  xlab("Numero de Reynolds (-)") +
  ylab("Factor f (-)") +
  theme_bw(base_size = 16.0) 

# A ggplot object is requested
fg02

# ////////////////////////////////////////////////////////
# BLOCK: Optimization
# ////////////////////////////////////////////////////////

# Subset data.frames are created
df.subset.cup18 <- subset(df.obs, pipe == "cup18" & feeding == "Pump_FluxAnal_ManDig")
df.subset.pvc12 <- subset(df.obs, pipe == "pvc12")
df.subset.pvc18 <- subset(df.obs, pipe == "pvc18")
df.subset.cup18.02 <- subset(df.obs, pipe == "cup18" & feeding == "Pump_Rot_ManDig")

vector.c.cup18 <- df.subset.cup18[,"f"]
vector.c.pvc12 <- df.subset.pvc12[,"f"]
vector.c.pvc18 <- df.subset.pvc18[,"f"]
vector.c.cup18.02 <- df.subset.cup18.02[,"f"]

vector.f.cup18 <- df.subset.cup18[,"Re"]
vector.f.pvc12 <- df.subset.pvc12[,"Re"]
vector.f.pvc18 <- df.subset.pvc18[,"Re"]
vector.f.cup18.02 <- df.subset.cup18.02[,"Re"]

# Custom function f.opti is manually evaluated
f.opti(par=0.0015/1000, DD=D1, RRe=vector.f.cup18, fexp=vector.c.cup18)
f.opti(par=0.0015/1000, DD=D2, RRe=vector.f.pvc12, fexp=vector.c.pvc12)
f.opti(par=0.0015/1000, DD=D3, RRe=vector.f.pvc18, fexp=vector.c.pvc18)
f.opti(par=0.0015/1000, DD=D1, RRe=vector.f.cup18.02, fexp=vector.c.cup18.02)

# Simple boxplots are manually requested to deal with outlier using ggplot2

# A ggplot object is created
fg03 <- ggplot() +
  geom_boxplot(aes(y = f,x = pipe,colour = group),data=df.obs,outlier.shape = 17,outlier.size = 3, size = 0.95) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("Boxplot de dispersion. Factor f experimental") +
  xlab("tuberia") +
  ylab("Factor f experimental (-)") +
  theme_grey(base_size = 16.0) 

# A ggplot object is requested
fg03

# Outliers are manually deleted (unselected) !!!!!!!!!!!!!!!!!!!
which(vector.c.cup18 > 0.024)
vector.c.cup18 <- vector.c.cup18[-c(1,25,26)]
vector.f.cup18 <- vector.f.cup18[-c(1,25,26)]

# Outliers are manually deleted (unselected)) !!!!!!!!!!!!!!!!!!!
which(vector.c.pvc12 > 0.040)
vector.c.pvc12 <- vector.c.pvc12[-c(11)]
vector.f.pvc12 <- vector.f.pvc12[-c(11)]

# Outliers are manually deleted (unselected)) !!!!!!!!!!!!!!!!!!!
which(vector.c.pvc18 > 0.030)
vector.c.pvc18 <- vector.c.pvc18[-c(11)]
vector.f.pvc18 <- vector.f.pvc18[-c(11)]

# Outliers are manually deleted (unselected) !!!!!!!!!!!!!!!!!!!
which(vector.c.cup18.02 > 0.024)
vector.c.cup18.02 <- vector.c.cup18.02[-c(1,25,26)]
vector.f.cup18.02 <- vector.f.cup18.02[-c(1,25,26)]

# **********************************************************
# Package optim{stats} is requested for optimization
# **********************************************************

# CUPPER 18 -------------------------------

# List container is reset
result.exp.D1 <- NULL

# optim {stats} is requested
result.exp.D1 <- optim(par = 0.010/1000,
                       f.opti,
                       lower = 0.0015/1000,
                       upper = 0.50/1000,
                       DD=D1,
                       RRe=vector.f.cup18,
                       fexp=vector.c.cup18,
                       method ="Brent",
                       control = list(maxit = 500))

# Return list is requested
result.exp.D1

# Parameter value-ONLY is requested
Exp_e1 <- result.exp.D1$par

# PVC 12 -------------------------------

# List container is reset
result.exp.D2 <- NULL

# optim {stats} is requested
result.exp.D2 <- optim(par = 0.010/1000,
                       f.opti,
                       lower = 0.0015/1000,
                       upper = 0.50/1000,
                       DD=D2,
                       RRe=vector.f.pvc12,
                       fexp=vector.c.pvc12,
                       method ="Brent",
                       control = list(maxit = 500))

# Return list is requested
result.exp.D2

# Parameter value-ONLY is requested
Exp_e2 <- result.exp.D2$par

# PVC 18 -------------------------------

# List container is reset
result.exp.D3 <- NULL

# optim {stats} is requested
result.exp.D3 <- optim(par = 0.010/1000,
                       f.opti,
                       lower = 0.0015/1000,
                       upper = 0.50/1000,
                       DD=D3,
                       RRe=vector.f.pvc18,
                       fexp=vector.c.pvc18,
                       method ="Brent",
                       control = list(maxit = 500))

# Return list is requested
result.exp.D3

# Parameter value-ONLY is requested
Exp_e3 <- result.exp.D3$par

# CUPPER 18 (Esquema02)-------------------------------

# List container is reset
result.exp.D1.02 <- NULL

# optim {stats} is requested
result.exp.D1.02 <- optim(par = 0.010/1000,
                       f.opti,
                       lower = 0.0015/1000,
                       upper = 0.50/1000,
                       DD=D1,
                       RRe=vector.f.cup18.02,
                       fexp=vector.c.cup18.02,
                       method ="Brent",
                       control = list(maxit = 500))

# Return list is requested
result.exp.D1.02

# Parameter value-ONLY is requested
Exp_e4 <- result.exp.D1.02$par

# ----------------------------------------------------------
# A new adjusted ggplot2 object must be created

# Theoretical f values are generated for each pipe based on Colebrook White equation
df.cup18.CW2 <- f.teo(ee=Exp_e1, DD=D1, RRe=Re)
df.pvc12.CW2 <- f.teo(ee=Exp_e2, DD=D2, RRe=Re)
df.pvc18.CW2 <- f.teo(ee=Exp_e3, DD=D3, RRe=Re)
df.cup18.CW2.02 <- f.teo(ee=Exp_e4, DD=D1, RRe=Re)

# Vector values are extracted from data.frames
cup18.CW2 <- df.cup18.CW2[,"root"]
pvc12.CW2 <- df.pvc12.CW2[,"root"]
pvc18.CW2 <- df.pvc18.CW2[,"root"]
cup18.CW2.02 <- df.cup18.CW2.02[,"root"]

# a data.frame is created
df.base2 <- data.frame(cup18.CW2,
                       pvc12.CW2,
                       pvc18.CW2,
                       cup18.CW2.02)

# a data.frame is created and Reynolds theoretical values ared added
df.melted2 <- melt(df.base2)
Re3 <- rep(Re,4)
df.melted2$Re <- Re3

# A ggplot object is created
fg04 <- ggplot() +
  geom_line(aes(x = Re,y = value,linetype = variable),data=df.melted2,size = 0.85) +
  geom_point(aes(x = Re,y = f,shape = feeding,colour = pipe),data=df.obs,size = 5.0) +
  scale_x_continuous(trans="log10", limits = c(1e3,1e8),
                     breaks = c(c(1e2,1e3,5e3,1e4,2e4,3e4,5e4,1e5,5e5,1e6,1e7,1e8,1e9,1e10))) +
  scale_y_continuous(trans="log10", 
                     breaks = c(c(0.008, 0.010,0.012,0.014,0.016,0.018,0.020, 0.025,0.030,0.035,
                                  0.040,0.045,0.050,0.055,0.060,0.07, 0.1))) +
  geom_text(aes(x = Re,y = f,label = q_l_s),size = 3.5,data=df.obs,parse = FALSE) +
  ggtitle("Diagrama de Moody Experimental (Colebrook White)") +
  xlab("Numero de Reynolds (-)") +
  ylab("Factor f (-)") +
  theme_bw(base_size = 16.0) 

# A ggplot object is requested
fg04

# round_df function is applied to relevant data.frames
df.obs <- round_df(df=df.obs, digits=5)

# Objects to export:
# fg01, fg02, fg03, fg04
# Exp_e1, Exp_e2, Exp_e3
# OUTLIERS:
# which(vector.c.cup18 > xxx)???
# which(vector.c.pvc12 > xxx)???
# which(vector.c.pvc18 > xxx)???
write.csv(df.obs, file = "df.obs.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
