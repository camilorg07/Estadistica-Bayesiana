#Autores:


#---------- Juan Andrés Camacho Zaráte (jucamachoz@unal.edu.co)
#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

#----------------------------- Primer Modelo -----------------------------

library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(metRology)

#1. Tratamiento de datos -----------------------------------------------------

sb2022_2=read_delim("datos/Saber 11 2022-2.TXT", delim = ";")

#1. Nacionalidad colombiana
#2. Residencia en Colombia
#3. Estado en publicar
#4. No tomar colegios de San Andrés
#5. Eliminar NAS

datos=sb2022_2%>%filter(ESTU_NACIONALIDAD=="COLOMBIA" & ESTU_PAIS_RESIDE=="COLOMBIA" & ESTU_ESTADOINVESTIGACION=="PUBLICAR" & COLE_COD_DEPTO_UBICACION!=88 & !is.na(COLE_COD_DEPTO_UBICACION) & !is.na(COLE_COD_MCPIO_UBICACION) & !is.na(PUNT_GLOBAL))

X=datos%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(nj=n())
nj=X$nj #numero de estudiantes por dpto
X=cbind(COLE_COD_DEPTO_UBICACION=X$COLE_COD_DEPTO_UBICACION, dpto=c(1:nrow(X)))

X=as.data.frame(X)

y=datos%>%left_join(X, by=c("COLE_COD_DEPTO_UBICACION"="COLE_COD_DEPTO_UBICACION"))%>%
  select(COLE_COD_DEPTO_UBICACION,dpto,PUNT_GLOBAL)%>%
  arrange(COLE_COD_DEPTO_UBICACION)


J <- nrow(X) #numero de dptos

#2. Funciones para samplear los parametros --------------------------------

#2.1 THETA ----------------------------------------------------------------

sample_theta <- function(y, mu, tau2, sij, theta, J, dpto)
{
  dt <- data.table(y = y, sij = sij, dpto = dpto)
  dt[, c("c", "d") := .(y / sij, 1 / sij)]
  
  # Agrupar por dpto y calcular sumas
  res <- dt[, .(
    sum_c = sum(c),
    sum_d = sum(d)
  ), by = dpto]
  
  # Calcular media y varianza para todos los grupos a la vez
  res[, m := (mu/tau2 + sum_c) / (1/tau2 + sum_d)]
  res[, v2 := 1 / (1/tau2 + sum_d)]
  
  # Muestrear todos los theta_j simultáneamente
  theta <- rnorm(J, res$m, sqrt(res$v2))
  return(theta)}

#2.2 Variables auxiliares --------------------------------------------------

sample_sij <- function(y, nj, v, sig2,theta, dpto, sij, n)
{
  b=((y-rep(theta,nj))^2)/2
  
  sij=1/rgamma(n=n, shape = (v+1)/2, rate = b+v*sig2/2)
  return(sij)
}

#2.3 MU ---------------------------------------------------------------------

sample_mu <- function(mu0, tau2, g20, J, theta, mu)
{
  m=(J*mean(theta)/tau2 + mu0/g20)/(J/tau2 + 1/g20)
  v2=1/(J/tau2 + 1/g20)
  
  mu=rnorm(n=1, mean = m, sd = sqrt(v2))
  return(mu)
}

#2.4 TAU^2 --------------------------------------------------------------------

sample_tau2 <- function(eta0, J, t20, theta, mu, tau2)
{
  a=(eta0+J)/2
  b=(sum((theta-mu)^2)+eta0*t20)/2
  
  tau2=1/rgamma(n=1, shape = a, rate = b)
  return(tau2)
}

#2.5 SIGMA^2 ------------------------------------------------------------------

sample_sig2 <- function(n, v, nu0, sij, s20, sig2)
{
  a=(n*v+nu0)/2
  b=(v*sum(1/sij)+nu0*s20)/2
  
  sig2=rgamma(n=1, shape = a, rate = b)
  return(sig2)
}

#3. Muestreador de Gibbs -------------------------------------------------------

#3.1 Hiperparametros -----------------------------------------------------------

v <- 3 #grados de libertad
mu0  <- 250 #mu_0
g20  <- 50^2 #gamma_0^2
eta0 <- 1 #eta_0 
t20  <- 50^2 #tau_0^2
nu0  <- 1  #nu_0
s20  <- 50^2 #sigma_0^2

#3.2 Muestreador ---------------------------------------------------------------

mcM1 <- function (y, nj, J, v, mu0, g20, eta0, t20, nu0, s20, n_sams, n_burn, n_skip, verbose = TRUE) 
  {
  
  n=nrow(y)
  a<-y%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(xbar=mean(PUNT_GLOBAL))%>%arrange(COLE_COD_DEPTO_UBICACION)
  b<-y%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(s=var(PUNT_GLOBAL))%>%arrange(COLE_COD_DEPTO_UBICACION)
  dpto=y$dpto
  y=y$PUNT_GLOBAL
  
  # número de iteraciones
  
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.1*B)
  
  # valores iniciales
  
  sig2  <- var(y)
  #sij <- 1/rgamma(n=n, shape =v/2, rate= v*sig2/2)
  mu <- mean(a$xbar)
  theta <- a$xbar
  tau2 <- var(a$xbar)
  
  Progress <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (i in 1:B) {
    
    # actualizar parámetros
    
    sij <- sample_sij(y, nj, v, sig2,theta, dpto, sij, n)
    sig2   <- sample_sig2 (n, v, nu0, sij, s20, sig2)
    theta  <- sample_theta(y, mu, tau2, sij, theta, J, dpto)
    mu  <- sample_mu(mu0, tau2, g20, J, theta, mu)
    tau2     <- sample_tau2(eta0, J, t20, theta, mu, tau2)
    
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      LL  <- sum(dt.scaled(y, v, mean = rep(theta,nj), sd = sqrt(sig2), log = TRUE))
      fwrite(as.data.table(t(theta)), file = "THETA_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig2)), file = "SIG2_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(tau2)), file = "TAU2_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(mu)), file = "MU_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(LL)), file = "LL_iterations-1.txt", append = TRUE, col.names = FALSE)
    }
    # progreso
    setTxtProgressBar(Progress, i)
  }
}

#4. Simulacion y almacenamiento ------------------------------------------------
  
n_sams=10000
n_burn=10000
n_skip=10

set.seed(123)

tictoc::tic()
mcM1(y, nj, J, v, mu0, g20, eta0, t20, nu0, s20, n_sams, n_burn, n_skip, verbose = TRUE) 
tictoc::toc()


MU<-fread("MU_iterations-1.txt")
LL<-fread("LL_iterations-1.txt")
TAU2<-fread("TAU2_iterations-1.txt")
SIG2<-fread("SIG2_iterations-1.txt")
THETA<-fread("THETA_iterations-1.txt")

M1 <- list(THETA = THETA, SIG2 = SIG2, TAU2 = TAU2, MU = MU, LL = LL)

save(M1, file =  "M1.RData")

load("M1.RData")
