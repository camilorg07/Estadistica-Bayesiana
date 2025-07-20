#Autor:

#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

#-----------------------Segundo Modelo----------------------------------

setwd("C:/Users/Lenovo/OneDrive - Universidad Nacional de Colombia/Documentos/Caso Bayesiana 2025 - I")

#0. Librerías ----

library(dplyr)
library(data.table)
library(tidyr)
library(Rcpp)
library(RcppArmadillo)

#1. Carga de datos y organización por Departamento y Municipio ----

#Base de datos
sb11<-read.csv("Datos.txt", sep=";")

#Covariables
Dep<-read.csv("Departamentales.txt", sep = ";")
Mun<-read.csv("Municipales.txt", sep = ";")

sb11 <- left_join(sb11, Mun, by=c("cole_cod_mcpio_ubicacion"="codmpio"))

sb11 <- left_join(sb11, Dep, by=c("cole_cod_depto_ubicacion"="Cod"))

sb11 <- sb11 %>% drop_na(edu_madre, comp, internet, libros, estrato, etnia,
                         PIB, porc_riesgo, porc_rul, doc_alum, nbi, RISK_VICTIM_2022)

#PIB en unidades de millones
sb11 <- sb11 %>% mutate(PIB = PIB/1e6)

sb11 <- sb11 %>% mutate(
                        edu_madre = as.factor(edu_madre),
                        comp = as.factor(comp),
                        internet = as.factor(internet),
                        libros = as.factor(libros),
                        estrato = as.factor(estrato),
                        etnia = as.factor(etnia)
)

##1.1 Medidas por departamento (k) ----

dptos <- sb11 %>% group_by(cole_cod_depto_ubicacion) %>% 
  summarise(yk_barra = mean(punt_global), #Media por dptos
            s2_k = var(punt_global), #Var por dptos
            nk = n_distinct(cole_cod_mcpio_ubicacion),  #Nro municipios
            nijk = n()) #Nro de estudiantes)

##1.2  Medidas por municipio (j) ----

mpios <- sb11 %>% group_by(cole_cod_mcpio_ubicacion, cole_cod_depto_ubicacion) %>% 
  summarise(yjk_barra = mean(punt_global), #Media por mpios
            s2_jk = var(punt_global), #Var por mpios
            njk = n(), #Nro de estudiantes
            .groups = "drop")

## 1.3 Matriz de diseño ----

fit <-  lm(punt_global~edu_madre+comp+internet+libros+estrato+etnia+
             PIB+porc_riesgo+porc_rul+doc_alum+nbi+RISK_VICTIM_2022 ,data = sb11)

X <- model.matrix(fit)

e <- length(c(2:14))
l <- length(c(15:17))
d <- length(c(18:20))

# Modelo 2

#2. Funciones para muestrear parámetros ----

##2.1 \beta y \kappa_{jk}^2 ----

sourceCpp("samplesfunctions.cpp")

##2.2 \sigma_{\beta} ----

sample_sig_beta <- function(mu_beta, nu_beta, g_beta, theta, sig_beta){
  
  beta <- theta[1]
  
  a = ( nu_beta + 1 ) / 2
  
  b = ( nu_beta*g_beta + sum( (beta-mu_beta)^2 ) ) / 2
  
  sig_beta <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_beta)
}

## 2.3 \sigma_{E} ----

sample_sig_E <- function(mu_E, nu_E, g_E, theta, e, sig_E){
  
  beta_E <- theta[2:14]
  
  a = ( nu_E + e ) / 2
  
  b = ( nu_E*g_E + t(beta_E-mu_E)%*%(beta_E-mu_E) ) / 2
  
  sig_E <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_E)
}

## 2.4 \sigma_{M} ----

sample_sig_M <- function(mu_M, nu_M, g_M, theta, l, sig_M){
  
  beta_M <- theta[15:17]
  
  a = ( nu_M + l ) / 2
  
  b = ( nu_M*g_M + t(beta_M-mu_M)%*%(beta_M-mu_M) ) / 2
  
  sig_M <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_M)
}

## 2.5 \sigma_{D} ----

sample_sig_D <- function(mu_D, nu_D, g_D, theta, d, sig_D){
  
  beta_D <- theta[18:20]
  
  a = ( nu_D + d ) / 2
  
  b = ( nu_D*g_D + t(beta_D-mu_D)%*%(beta_D-mu_D) ) / 2
  
  sig_D <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_D)
}

## 2.6 \kappa_{k}^2 ----

sample_k2k <- function(alpha_k, m, nk, nu_k, beta_k, k2jk, deptos, k2k){
  
  dt <- data.table(k2jk=k2jk, deptos=deptos)
  res <- dt[,.(pj =sum(1/k2jk)), by=deptos]
  setorder(res, deptos)
  
  a = ( alpha_k + nk*nu_k ) / 2
  
  b = ( beta_k + nu_k*res$pj ) / 2
  
  k2k <- rgamma(m, shape = a, rate = b)
  
  return(k2k)
}

## 2.7 \beta_{k} ----

sample_beta_k <- function(a_bk, m, alpha_k, b_bk, k2k, beta_k){
  
  a = ( 2*a_bk + m*alpha_k ) / 2
  
  b = ( 2*b_bk + sum(k2k) ) / 2
  
  beta_k = rgamma(1, shape = a, rate = b)
  
  return(beta_k)
}

## 2.8 \alpha_{k}

### 2.8.1 log(p(\alpha_k | resto)) ----

post <- function(alpha_k, m, beta_k, k2k){
  (m*alpha_k/2)*log(beta_k/2) -m*lgamma(alpha_k/2) + (alpha_k/2)*sum(log(k2k)) + (a_ak-1)*log(alpha_k) - b_ak*alpha_k
}

### 2.8.2 Tasa de aceptación adaptativa ----

tunning <- function(delta, n_tun, mix_rate,i) {
  target_rate <- 0.35
  tolerance <- 0.05
  max_delta <- 10.0
  min_delta <- 1e-10
  
  if (i %% n_tun == 0) {
    diff <- mix_rate - target_rate
    if (abs(diff) > tolerance) {
      delta <- delta * exp(0.1 * diff)
      delta <- max(min(delta, max_delta), min_delta)  # clamp
      n_tun <- 100
    } else {
      n_tun <- n_tun + 100
    }
  }
  
  return(list(delta = delta, n_tun = n_tun,ac = ac))
}

### 2.8.3 Metropolis ----

sample_alpha_k <- function(detla, m, i, n_tun, ac, beta_k, k2k, n_burn, alpha_k){
  
  #1. Proponer nuevo valor
  gamma_c <- log(alpha_k)
  gamma_p <- rnorm(1, mean = gamma_c, sd = sqrt(delta))
  alpha_p <- exp(gamma_p)
  
  #2. Calculo tasa de acpetación
  
  r = exp(post(alpha_p, m, beta_k, k2k)- post(alpha_k, m, beta_k, k2k)
          + log(alpha_p)- log(alpha_k))
  
  #3. Actualización del parámetro
  
  if (runif(1) < r) {
    alpha_k <- alpha_p
    ac <- ac + 1
  }
  
  #4. Adaptación solo durante el periodo de calentamiento
  
  if (i <= n_burn) {
    mix_rate <- ac/i  
    ajuste <- tunning(delta, n_tun, mix_rate,i)
    delta <- ajuste$delta
    n_tun <- ajuste$n_tun
  }
  
  return(list(ac = ac, delta = delta, n_tun = n_tun, alpha_k = alpha_k))
}

# 3. Cadena de Montecarlo ----

## 3.1 Hiperparametros ----

mu_beta<-mu_M <- mu_E <- mu_D <- 0
nu_beta <- nu_M <- nu_E <- nu_D <- nu_k <- 1
g_beta <- g_M <- g_E <- g_D <- 10^2
a_ak <- 1
b_ak <- 1
a_bk <-1
b_bk <-1/50^2
delta <- 0.1
n_tun <- 100

n_sams <- 10000
n_burn <- 10000
n_skip <- 10
ac <- 0

## 3.2 Gibs sampler ----

M2MC <- function(dptos, mpios, sb11,
                 X, Z, W, e, l, d, fit,
                 mu_beta, mu_M, mu_E, mu_D, nu_beta, nu_M, nu_E, nu_D, nu_k, a_ak, b_ak,a_bk, b_bk,
                 delta, n_tun, ac,
                 n_sams, n_burn, n_skip, verbose=T){
  
  y <- sb11$punt_global
  y <- as.numeric(y)
  X <- as.matrix(X)
  tX <- t(X)
  deptos<-mpios$cole_cod_depto_ubicacion
  

  #Estadísticos suficientes
  n <- nrow(sb11)
  
  m <- nrow(dptos)
  nk <- dptos$nk
  yk_barra <- dptos$yk_barra
  s2_k <- dptos$s2_k
  
  nj <- nrow(mpios)
  njk <- mpios$njk
  yjk_barra <- mpios$yjk_barra
  s2_jk <- mpios$s2_jk
  mpios_id <- rep(c(0:(nj-1)), njk)
  deptos_levels <- unique(dptos$cole_cod_depto_ubicacion)
  deptos_id <- match(mpios$cole_cod_depto_ubicacion, deptos_levels) - 1
  beta_0 <- c(mu_beta, rep(mu_E,e), rep(mu_M,l), rep(mu_D,d))

  
  #Valores iniciales
  
  theta <- fit$coefficients
  k2jk <- s2_jk
  alpha_k <- rgamma(1, shape = a_ak, b_ak)
  beta_k <- rgamma(1, shape = a_bk, b_bk)
  k2_k <- rgamma(m, shape = alpha_k/2, rate = beta_k/2)
  sig_beta <- 1/rgamma(1, shape=nu_beta/2, rate=nu_beta*g_beta/2)
  sig_E <- 1/rgamma(1, shape=nu_E/2, rate=nu_E*g_E/2)
  sig_M <- 1/rgamma(1, shape=nu_M/2, rate=nu_M*g_M/2)
  sig_D <- 1/rgamma(1, shape=nu_D/2, rate=nu_D*g_D/2)
  

  #Cadena y progreso
  
  B <- n_burn + n_sams*n_skip
  
  Progress <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (i in 1:B) {
    
    # actualizar parámetros
    theta <- sample_theta_cpp(k2jk, njk, sig_beta, sig_E, sig_M, sig_D, beta_0, X, tX, y, e, l, d)
    theta <- matrix(theta, nrow = 1)
    k2jk <- sample_k2jk_cpp(y, nu_k, deptos_id, njk, nj, theta, X, mpios_id, k2_k)
    sig_beta <- sample_sig_beta(mu_beta, nu_beta, g_beta, theta, sig_beta)
    sig_E <- sample_sig_E(mu_E, nu_E, g_E, theta, e, sig_E)
    sig_M <- sample_sig_M(mu_M, nu_M, g_M, theta, l, sig_M)
    sig_D <- sample_sig_D(mu_D, nu_D, g_D, theta, d, sig_D)
    k2k <- sample_k2k(alpha_k, m, nk, nu_k, beta_k, k2jk, deptos, k2k)
    beta_k <- sample_beta_k(a_bk, m, alpha_k, b_bk, k2k, beta_k)
    a <- sample_alpha_k(delta, m, i, n_tun, ac, beta_k, k2k, n_burn, alpha_k)
    
    ac <- a$ac
    delta <- a$delta
    n_tun <- a$n_tun
    alpha_k <- a$alpha_k
    sjk <- X%*%t(theta)
    
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      LL  <- sum(dnorm(y, mean = sjk, sd = sqrt(rep(k2jk, njk)), log=T))
      fwrite(as.data.table(t(LL)), file = "LL_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(theta), file = "THETA_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig_beta)), file = "SIG_B_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig_E)), file = "SIG_E_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig_M)), file = "SIG_M_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig_D)), file = "SIG_D_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(k2jk)), file = "KAPPA_JK_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(k2k)), file = "KAPPA_K_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(alpha_k)), file = "ALPHA_iterations-2.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(beta_k)), file = "BETA_K_iterations-2.txt", append = TRUE, col.names = FALSE)
    }
    # progreso
    setTxtProgressBar(Progress, i)
  }
}

#4. Simulacion ----

## 4.1 Ejecución del algoritmo ----

set.seed(7)

tictoc::tic()

M2MC(dptos, mpios, sb11,
     X, Z, W, e, l, d, fit,
     mu_beta, mu_M, mu_E, mu_D, nu_beta, nu_M, nu_E, nu_D, nu_k, a_ak, b_ak,a_bk, b_bk,
     delta, n_tun, ac,
     n_sams, n_burn, n_skip, verbose=T)

tictoc::toc()

# 5 Almacenamiento ----

THETA <- fread("THETA_iterations-2.txt")
KAPPA_JK<-fread("KAPPA_JK_iterations-2.txt")
KAPPA_K<-fread("KAPPA_K_iterations-2.txt")
SIG_B<-fread("SIG_B_iterations-2.txt")
SIG_E<-fread("SIG_E_iterations-2.txt")
SIG_M<-fread("SIG_M_iterations-2.txt")
SIG_D<-fread("SIG_D_iterations-2.txt")
ALPHA<-fread("ALPHA_iterations-2.txt")
BETA_K<-fread("BETA_K_iterations-2.txt")
LL <- fread("LL_iterations-2.txt")

M2 <- list(THETA = THETA,
           KAPPA_JK = KAPPA_JK,
           SIG_B = SIG_B,
           SIG_E = SIG_E,
           SIG_M = SIG_M,
           SIG_D = SIG_D,
           KAPPA_K = KAPPA_K,
           BETA_K = BETA_K,
           ALPHA = ALPHA,
           LL = LL)

save(M2, file = "M2.RData")