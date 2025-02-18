#---------------------------- Tercer modelo--------------------------------------
#Autores:


#---------- Juan Andrés Camacho Zaráte (jucamachoz@unal.edu.co)
#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

#--------------------------- Tercer Modelo -------------------------------------

library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(metRology)
library(knitr)

#1. Tratamiento de datos ----------------------------------------------
sb2022_2=read_delim("Saber 11 2022-2.TXT", delim = ";")

#1. Nacionalidad colombiana
#2. Residencia en Colombia
#3. Estado en publicar
#4. No tomar colegios de San AndrC)s
#5. Eliminar NAS

datos=sb2022_2%>%filter(ESTU_NACIONALIDAD=="COLOMBIA" & ESTU_PAIS_RESIDE=="COLOMBIA" & ESTU_ESTADOINVESTIGACION=="PUBLICAR" & COLE_COD_DEPTO_UBICACION!=88 & !is.na(COLE_COD_DEPTO_UBICACION) & !is.na(COLE_COD_MCPIO_UBICACION) & !is.na(PUNT_GLOBAL))

rm(sb2022_2)

# Crear identificador único para departamento y municipio
datos <- datos %>%
  mutate(COLE_COD_DEPTO_MCPIO = paste(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION, sep = "-"))

# Agrupar correctamente por departamento y municipio
X_deptos <- datos %>%
  group_by(COLE_COD_DEPTO_UBICACION) %>%
  summarise(nj_depto = n_distinct(COLE_COD_MCPIO_UBICACION)) %>%
  mutate(dpto = 1:n())

X=as.data.frame(X_deptos)

y <- datos %>%
  left_join(X_deptos, by = "COLE_COD_DEPTO_UBICACION") %>%
  select(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION, dpto, PUNT_GLOBAL) %>%
  arrange(dpto, COLE_COD_MCPIO_UBICACION) %>%
  mutate(mun_id = group_indices(., COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION))

# Número de municipios (K) y estudiantes por municipio (njk)
muni <- y %>%
  group_by(dpto, mun_id) %>%
  summarise(njk = n(), .groups = 'drop')

K <- nrow(muni)
nj <- muni %>% group_by(dpto) %>% summarise(nj = n()) %>% pull(nj)
njk <- muni$njk
J <- length(unique(y$dpto))  # Número de departamentos
#2. Funciones para samplear los parametros ---------------------------

#2.1 Medias por municipio ----------------------------------------------

library(data.table)

sample_epsjk <- function(y, mun, sijk, theta, nj, sig2, K, epsjk) {
  
  # Crear un data.table con las variables necesarias
  dt <- data.table(mun = mun, c = y / sijk, d = 1 / sijk)
  
  # Agrupar por municipio y calcular las sumas
  sum_values <- dt[, .(sum_c = sum(c), sum_d = sum(d)), by = mun]
  
  # Calcular m y V2
  m <- (sum_values$sum_c + rep(theta, nj) / sig2) / (sum_values$sum_d + 1 / sig2)
  V2 <- 1 / (sum_values$sum_d + 1 / sig2)
  
  # Generar las medias por municipio
  epsjk <- rnorm(n = K, mean = m, sd = sqrt(V2))
  
  return(epsjk)
}


#2.2 Variables auxiliares -------------------------------------------

sample_sijk <- function(v, k2, epsjk, njk, y, n, sijk)
{
  c = (y - rep(epsjk, njk))^2
  b = (c + v*k2)/2
  
  #Generar las latentes
  
  sijk <- 1/rgamma(n = n, shape = (v+1)/2, rate = b)
  
  return(sijk)
}

#2.3 THETA ----------------------------------------------------------

sample_theta <- function(epsjk, sig2, tau2, mu, J, dpto, nj, theta) {
  
  # Crear un data.table con epsjk y el identificador del departamento
  dt <- data.table(dpto = rep(1:J, nj), epsjk = epsjk)
  
  # Agrupar por departamento y calcular la suma de epsjk
  res <- dt[, .(c = sum(epsjk)), by = dpto]
  
  # Calcular m y V2
  m <- (res$c / sig2 + mu / tau2) / (nj / sig2 + 1 / tau2)
  V2 <- 1 / (nj / sig2 + 1 / tau2)
  
  # Generar las medias por departamento
  theta <- rnorm(n = J, mean = m, sd = sqrt(V2))
  
  return(theta)
}

#2.4 SIGMA^2 -------------------------------------------------------

#2.1 sample_sig2 Corregido
sample_sig2 <- function(K, nu0, epsjk, theta, s20, nj, sig2) {
  c <- sum((epsjk - rep(theta, nj))^2)
  a <- (K + nu0)/2  # Usar K en lugar de sum(nk)
  b <- (c + nu0 * s20)/2
  sig2 <- 1/rgamma(1, shape = a, rate = b)
  return(sig2)
}

#2.5 MU ------------------------------------------------------------

sample_mu <- function(theta, mu0, g20, tau2, J, mu)
{
  m <- (J*mean(theta)/tau2 + mu0/g20)/(J/tau2 + 1/g20)
  V2 <- 1/(J/tau2 + 1/g20)
  
  #Generar mu
  
  mu <- rnorm(n=1, mean = m, sd = sqrt(V2))
  return(mu)
}

#2.6 TAU^2 ---------------------------------------------------------

sample_tau2 <- function(eta0, J, theta, mu, t20, tau2)
{
  a <- (eta0 + J)/2
  b <- (sum((theta-mu)^2) + eta0*t20)/2
  
  #Generar tau^2
  
  tau2 <- 1/rgamma(n = 1, shape = a, rate = b)
  return(tau2)
}

#2.7 K^2 -----------------------------------------------------------

sample_k2 <- function(v, nk, eps0, sijk, k20)
{
  a = (v*sum(nk) + eps0)/2
  b = (v*sum(1/sijk) + eps0*k20)/2
  
  #Generar k^2
  
  k2 <- rgamma(n= 1, shape = a, rate = b)
  return(k2)
}

#3. Muestreador de Gibbs -------------------------------------------

#3.1 Hiperparametros------------------------------------------------

v <- 3 #grados de libertad
eps0 <- 1 #epsilon_0
k20 <- 50^2 #kappa_0^2
mu0 <- 250 #mu_0
g20 <- 50^2 #gamma_0^2
eta0 <- 1 #eta_0
t20 <- 50^2 #tau_0^2
nu0 <- 1 #nu_0
s20 <- 50^2 #sigma_0^2

#3.2 Muestreador ----------------------------------------------------

mcM3 <- function (y, J, K, nj, nk, njk, v, eps0, k20, mu0, g20, eta0, t20, nu0, s20, n_sams, n_burn, n_skip, verbose = TRUE) 
{
  n=nrow(y)
  dpto=y$dpto
  mun=y$COLE_COD_MCPIO_UBICACION
  y=y$PUNT_GLOBAL
  
  # nC:mero de iteraciones
  
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.1*B)
  
  # valores iniciales
  
  sig2  <- 1/rgamma(n = 1, shape = nu0/2, rate = nu0*s20/2)
  k2 <- rgamma(n = 1, shape = eps0/2, rate = eps0*k20/2)
  sijk <- 1/rgamma(n=n, shape =v/2, rate= v*k2/2)
  mu <- rnorm(n = 1, mean = mu0, sd= sqrt(g20))
  theta <- rnorm(n = J, mean = mu, sd = sqrt(t20))
  tau2 <- 1/rgamma(n=1, shape = eta0/2, rate = eta0*t20/2)
  epsjk <- rnorm(n = K, mean = rep(theta,nj), sd = sqrt(sig2))
  
  # almacenamiento
  
  for (i in 1:B) {
    
    # actualizar parC!metros
    
    sijk <- sample_sijk(v, k2, epsjk, njk, y, n, sijk)
    k2 <- sample_k2(v, nk, eps0, sijk, k20)
    sig2 <- sample_sig2(K, nu0, epsjk, theta, s20, nj, sig2)
    
    #sig2 <- sample_sig2(nk, nu0, epsjk, theta, s20, nj, sig2)
    theta <- sample_theta(epsjk, sig2, tau2, mu, J, dpto, nj, theta)
    mu <- sample_mu(theta, mu0, g20, tau2, J, mu)
    tau2<- sample_tau2(eta0, J, theta, mu, t20, tau2)
    epsjk <-sample_epsjk(y, mun, sijk, theta, nj, sig2, K, epsjk)
    
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      LL  <- sum(dt.scaled(y, v, mean = rep(epsjk,njk), sd = sqrt(k2), log = TRUE))
      write.table(as.data.frame(t(epsjk)), file = "EPSJK_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(theta)), file = "THETA_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(sig2)), file = "SIG2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(tau2)), file = "TAU2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(LL)), file = "LL_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(k2)), file = "K2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(mu)), file = "MU_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
    
    # progreso
    
    if (verbose && i %% ncat == 0)
      cat(sprintf("%.1f%% completado\n", 100*i/B))
  }
}

#4. Simulacion y almacenamiento -----------------------------------------


n_sams=10000
n_burn=10000
n_skip=10

set.seed(123)

mcM3(y, J, K, nj, nk, njk, v, eps0, k20, mu0, g20, eta0, t20, nu0, s20, n_sams, n_burn, n_skip, verbose = TRUE)

#Carga de muestras


EPSJK <- fread("EPSJK_iterations.txt")
THETA <- fread("THETA_iterations.txt")
MU <- fread("MU_iterations.txt")
SIG2 <- fread("SIG2_iterations.txt")
TAU2 <- fread("TAU2_iterations.txt")
K2 <- fread("K2_iterations.txt")
LL <- fread("LL_iterations.txt")

#Compilar todas las muestras en una lista

M3 <- list(EPSJK = EPSJK, THETA = THETA, MU = MU, SIG2 = SIG2, TAU2 = TAU2, K2 = K2, LL = LL)

save(M3, file = "M3.RData")

#load("M3.RData")

