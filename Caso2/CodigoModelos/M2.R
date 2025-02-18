#Autores:


#---------- Juan Andrés Camacho Zaráte (jucamachoz@unal.edu.co)
#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

#----------------------------- Segundo Modelo --------------------------------

library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(metRology)

#1. Tratamiento de datos --------------------------------------------

sb2022_2=read_delim("datos/Saber 11 2022-2.TXT", delim = ";")

#1. Nacionalidad colombiana
#2. Residencia en Colombia
#3. Estado en publicar
#4. No tomar colegios de San Andrés
#5. Eliminar NAS

datos=sb2022_2%>%filter(ESTU_NACIONALIDAD=="COLOMBIA" & ESTU_PAIS_RESIDE=="COLOMBIA" & ESTU_ESTADOINVESTIGACION=="PUBLICAR" & COLE_COD_DEPTO_UBICACION!=88 & !is.na(COLE_COD_DEPTO_UBICACION) & !is.na(COLE_COD_MCPIO_UBICACION) & !is.na(PUNT_GLOBAL))

rm(sb2022_2)

X=datos%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(nj=n())
X=cbind(COLE_COD_DEPTO_UBICACION=X$COLE_COD_DEPTO_UBICACION, dpto=c(1:nrow(X)))
X=as.data.frame(X)

y=datos%>%left_join(X, by=c("COLE_COD_DEPTO_UBICACION"="COLE_COD_DEPTO_UBICACION"))%>%
  select(COLE_COD_DEPTO_UBICACION,dpto,PUNT_GLOBAL)%>%
  arrange(COLE_COD_DEPTO_UBICACION)

J<-y%>%summarise(n=n_distinct(COLE_COD_DEPTO_UBICACION))%>%select(n)%>%as.numeric() #numero de dptos
nj=(y%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(nj=n()))$nj #estudiantes por dpto

#2. Funciones para samplear los parametros --------------------------------------------

#2.1 THETA ------------------------------------------------------------

sample_theta <- function(y, mu, tau2, sij, theta, J, dpto) {
  
  # Crear un data.table en lugar de data.frame
  dt <- data.table(c = y / sij, d = 1 / sij, dpto = dpto)
  
  # Calcular sum_c y sum_d de forma eficiente con data.table
  sums <- dt[, .(sum_c = sum(c), sum_d = sum(d)), by = dpto]
  
  # Extraer valores para cálculos posteriores
  sum_c <- sums$sum_c
  sum_d <- sums$sum_d
  
  # Calcular los parámetros de la distribución normal
  m <- (mu / tau2 + sum_c) / (1 / tau2 + sum_d)
  v2 <- 1 / (1 / tau2 + sum_d)
  
  # Generar las muestras de theta
  theta <- rnorm(J, mean = m, sd = sqrt(v2))
  
  return(theta)
}
#2.2 MU ---------------------------------------------------------------

sample_mu <- function(mu0, tau2, g20, J, theta, mu)
{
  m=(J*mean(theta)/tau2 + mu0/g20)/(J/tau2 + 1/g20)
  v2=1/(J/tau2 + 1/g20)
  
  mu=rnorm(n=1, mean = m, sd = sqrt(v2))
  return(mu)
}

#2.3 TAU^2 ------------------------------------------------------------

sample_tau2 <- function(eta0, J, t20, theta, mu, tau2)
  {
  a=(eta0+J)/2
  b=(sum((theta-mu)^2)+eta0*t20)/2
  
  tau2=1/rgamma(n=1, shape = a, rate = b)
  return(tau2)
}

#2.4 Variables auxiliares --------------------------------------------

sample_sij <- function(y, nj, v, sig2j,theta, dpto, sij, n)
{
  b=((y-rep(theta,nj))^2+v*rep(sig2j,nj))/2
  
  sij=1/rgamma(n=n, shape = (v+1)/2, rate = b)
  return(sij)
}

#2.5 SIGMA_j^2 ---------------------------------------------------------

sample_sig2j <- function(v, nj, nu, sig2, dpto, sij, sig2j,J) {
  
  # Crear un data.table en lugar de data.frame
  dt <- data.table(c = 1 / sij, dpto = dpto)
  
  # Calcular la suma por grupo (más eficiente que dplyr)
  res <- dt[, .(sum_c = sum(c)), by = dpto]
  
  # Calcular los parámetros de la distribución gamma
  a <- (nj * v + nu) / 2
  b <- (v * res$sum_c + nu * sig2) / 2
  
  # Generar las muestras de sig2j
  sig2j <- rgamma(n = J, shape = a, rate = b)
  
  return(sig2j)
}
#2.6 SIGMA^2 -----------------------------------------------------------

sample_sig2 <- function(a0, b0, nu, J, sig2j, sig2)
{
  a=(a0+nu*J)/2
  b=(b0+nu*sum(sig2j))/2
  
  sig2=rgamma(n=1, shape = a, rate = b)
  return(sig2)
}

#2.7 NU ----------------------------------------------------------------


sample_nu <- function(nus0, sig2, J, sig2j, lam0, nu) {
  lpnu <- 0.5 * nus0 * J * log(0.5 * nus0 * sig2) - 
    J * lgamma(0.5 * nus0) + 
    0.5 * nus0 * sum(log(sig2j)) - 
    nus0 * (lam0 + 0.5 * sig2 * sum(sig2j))
  prob <- exp(lpnu - max(lpnu))
  nu <- sample(nus0, size = 1, prob = prob)
  return(nu)
}



#3. Muestreador de Gibbs -----------------------------------------------

#3.1 Hiperparametros ---------------------------------------------------

v <- 3 #grados de libertad
mu0 <- 250 #mu_0
g20 <- 50^2 #gamma_0^2
eta0 <- 1 #eta_0
t20 <- 50^2 #tau_0^2
lam0 <- 1 #lambda_0
a0 <- 1 #alpha_0
b0 <- 1/(50^2) #beta_0
nus0 <- 1:50 #rango de valores para nu

#3.2 Muestreador ----------------------------------------------------------

mcM2 <- function (y, nj, J, v, mu0, g20, eta0, t20, lam0, a0, b0, nus0, n_sams, n_burn, n_skip, verbose = TRUE) 
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
  
  sig2  <- rgamma(n = 1, shape = a0/2, rate = b0/2)
  nu <- 1
  sig2j <- b$s
  sij <- 1/rgamma(n=n, shape =v/2, rate= v*rep(sig2j,nj)/2)
  mu <- mean(a$xbar)
  theta <- a$xbar
  tau2 <- var(a$xbar)
  
  for (i in 1:B) {
    
    # actualizar parámetros
    sij <- sample_sij(y, nj, v, sig2j,theta, dpto, sij, n)
    sig2j <- sample_sig2j(v, nj, nu, sig2, dpto, sij,sig2j, J)
    sig2   <- sample_sig2 (a0, b0, nu, J, sig2j, sig2)
    theta  <- sample_theta(y, mu, tau2, sij, theta, J, dpto)
    mu  <- sample_mu(mu0, tau2, g20, J, theta, mu)
    tau2<- sample_tau2(eta0, J, t20, theta, mu, tau2)
    nu<-sample_nu(nus0, sig2, J, sig2j, lam0, nu)
    
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      LL  <- sum(dt.scaled(y, v, mean = rep(theta,nj), sd = sqrt(rep(sig2j,nj)), log = TRUE))
      fwrite(as.data.table(t(theta)), file = "THETA_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig2j)), file = "SIG2J_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig2)), file = "SIG2_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(tau2)), file = "TAU2_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(mu)), file = "MU_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(nu)), file = "NU_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(LL)), file = "LL_iterations-2.txt", append = TRUE, col.names = FALSE)
    }
    
    # progreso
    
    if (verbose && i %% ncat == 0)
      cat(sprintf("%.1f%% completado\n", 100*i/B))
  }
}

#4. Simulacion y almacenamiento ---------------------------------------------------

n_sams=10000
n_burn=10000
n_skip=10

set.seed(123)

mcM2(y, nj, J, v, mu0, g20, eta0, t20, lam0, a0, b0, nus0, n_sams, n_burn, n_skip, verbose = TRUE) 

#Cargar las muestras

THETA <- fread("THETA_iterations-2.txt")
SIG2 <- fread("SIG2_iterations-2.txt")
SIG2J <- fread("SIG2J_iterations-2.txt")
TAU2 <- fread("TAU2_iterations-2.txt")
MU <- fread("MU_iterations-2.txt")
NU <- fread("NU_iterations-2.txt")
LL <- fread("LL_iterations-2.txt")

#Compilar las muestras en una lista

M2 <- list(THETA = THETA, SIG2 = SIG2, SIG2J = SIG2J, TAU2 = TAU2, MU = MU, NU = NU, LL = LL)

save(M2, file = "M2.RData")
