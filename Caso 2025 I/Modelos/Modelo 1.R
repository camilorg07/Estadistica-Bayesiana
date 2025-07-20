#Autor:

#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

#-----------------------Primer Modelo----------------------------------

setwd("C:/Users/Lenovo/OneDrive - Universidad Nacional de Colombia/Documentos/Caso Bayesiana 2025 - I")

# 0. Librerías -----

library(dplyr)
library(tidyr)
library(data.table)

#1. Carga de datos y organizaciC3n por Departamento y Municipio ----

sb11<-read.csv("Datos.txt", sep=";")



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


#2. Funciones para muestrear los parC!metros ----

##2.1 \zeta_{j,k} ----

sample_sjk <- function(njk, nj, k2jk, sig_k, yjk_barra, theta, nk, sjk){
  
  v = 1/( (njk/k2jk) + (1/rep(sig_k,nk)) ) 
            
  mean = ( njk*yjk_barra/k2jk + rep(theta/sig_k,nk) ) * v

  sjk = rnorm(nj, mean =  mean, sd = sqrt(v))
  
  return(sjk)
}

##2.2 \kappa_{j,k}^2 ----

sample_k2jk <- function(nu_k, nk, njk, nj, sjk,k2_k, s2_jk, yjk_barra, k2jk){
  
  a <- ( nu_k + njk ) / 2

  b <- ( nu_k*rep(k2_k,nk) + (njk-1)*s2_jk + njk*(yjk_barra-sjk)^2 ) / 2
  
  k2jk <- 1/rgamma(nj, shape = a, rate = b )
  
  return(k2jk)
}

##2.3 \theta_k ----

sample_theta <- function(nk, m, sig_k, tau2, sjk, mu, deptos, theta){
  
  dt <- data.table(sjk=sjk, deptos=deptos)
  res <- dt[,.(sj =sum(sjk)), by=deptos]
  setorder(res, deptos)
  
  v = 1/( nk/sig_k + 1/tau2 )
  
  mean = ( res$sj/sig_k + mu/tau2 ) * v
  
  theta = rnorm(m, mean = mean, sd = sqrt(v))
  
  return(theta)
}

##2.4 \sigma_{k}^2 ----

sample_sig_k <- function(nu_s, m, nk, sig2, sjk, theta, deptos, sig_k){
  
  dt <- data.table(sjk=sjk, theta=rep(theta,nk), deptos=deptos)
  res <- dt[,.(s2j =sum((sjk-theta)^2)), by=deptos]
  setorder(res, deptos)
  
  a = ( nu_s + nk ) / 2
  
  b = ( nu_s*sig2 + res$s2j ) / 2
  
  sig_k <- 1/rgamma(m, shape = a, rate = b)
  
  return(sig_k)
}

##2.5 \k^2_{k} ----

sample_k2k <- function(alpha_k, m, nk, nu_k, beta_k, k2jk, deptos, k2k){
  
  dt <- data.table(k2jk=k2jk, deptos=deptos)
  res <- dt[,.(pj =sum(1/k2jk)), by=deptos]
  setorder(res, deptos)

  a = ( alpha_k + nk*nu_k ) / 2
  
  b = ( beta_k + nu_k*res$pj ) / 2
  
  k2k <- rgamma(m, shape = a, rate = b)
  
  return(k2k)
}

##2.6 \mu ----

sample_mu <- function(m, tau2, sig_mu, theta, mu_mu, mu){
  
  v = 1 / ( m/tau2 + 1/sig_mu )
  
  mean = ( m*mean(theta)/tau2 + mu_mu/sig_mu ) * v
  
  mu <- rnorm(1, mean = mean, sd = sqrt(v))
  
  return(mu)
}

##2.7 \tau^2 ----

sample_tau2 <- function(nu_tau, m, sig_tau, theta, mu, tau2){
  
  a = ( nu_tau + m ) / 2
  
  b = ( nu_tau*sig_tau + sum( (theta-mu)^2 ) ) / 2
  
  tau2 <- 1/rgamma(1, shape = a, rate = b)
  
  return(tau2)
}

##2.8 \sigma^2 ----

sample_sig2 <- function(a_sig, m, nu_s, b_sig, sig_k, sig2){

  a = ( a_sig + m*nu_s ) / 2
  
  b = ( b_sig + nu_s*sum(1/sig_k) ) / 2
  
  sig2 <- rgamma(1, shape = a, rate = b)
  
  return(sig2)
}

##2.9 \beta_{k} ----

sample_beta_k <- function(a_bk, m, alpha_k, b_bk, k2k, beta_k){

  a = ( 2*a_bk + m*alpha_k ) / 2
  
  b = ( 2*b_bk + sum(k2k) ) / 2
  
  beta_k = rgamma(1, shape = a, rate = b)
  
  return(beta_k)
}

##2.10 \alpha_{k} ----

## 2.11 log(p(\alpha_k | resto)) ----

post <- function(alpha_k, m, beta_k, k2k){
  (m*alpha_k/2)*log(beta_k/2) -m*lgamma(alpha_k/2) + (alpha_k/2)*sum(log(k2k)) + (a_ak-1)*log(alpha_k) - b_ak*alpha_k
}

## 2.12 Tasa de aceptaciC3n adaptativa ----

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

## 2.13 Metropolis ----

sample_alpha_k <- function(delta, m, i, n_tun, ac, beta_k, k2k, n_burn, alpha_k){
  
  #1. Proponer nuevo valor
  gamma_c <- log(alpha_k)
  gamma_p <- rnorm(1, mean = gamma_c, sd = sqrt(delta))
  alpha_p <- exp(gamma_p)
  
  #2. Calculo tasa de acpetaciC3n
  
  r = exp(post(alpha_p, m, beta_k, k2k)- post(alpha_k, m, beta_k, k2k)
          + log(alpha_p)- log(alpha_k))
  
  #3. ActualizaciC3n del parC!metro
  
  if (runif(1) < r) {
    alpha_k <- alpha_p
    ac <- ac + 1
  }
  
  #4. AdaptaciC3n solo durante el periodo de calentamiento
  
  if (i <= n_burn) {
    mix_rate <- ac/i  
    ajuste <- tunning(delta, n_tun, mix_rate,i)
    delta <- ajuste$delta
    n_tun <- ajuste$n_tun
  }
  
  return(list(ac = ac, delta = delta, n_tun = n_tun, alpha_k = alpha_k))
}


#3. Cadena de Montecarlo ----

## 3.1 HiperparC!metros ----

nu_k <- 1
nu_s <- 1
mu_mu <- 250
sig_mu <- 50^2
nu_tau <- 1
sig_tau <- 50^2
a_sig <- 1
b_sig <- 1/50^2
a_ak <- 1
b_ak <- 1
a_bk <-1
b_bk <-1/50^2
delta <- 0.1
n_tun <- 100

## 3.3 Iteraciones y cadena

n_sams <- 10000
n_burn <- 10000
n_skip <- 10
ac <- 0

M1MC <- function(dptos, mpios, sb11,
                 nu_k,nu_s,mu_mu,sig_mu,nu_tau,sig_tau,a_sig,b_sig,a_ak,b_ak,a_bk,b_bk,
                 delta, n_tun, ac,
                 n_sams, n_burn, n_skip, verbose=T){
  
  y <- sb11$punt_global
  deptos<-mpios$cole_cod_depto_ubicacion
  
  #EstadC-sticos suficientes
  
  m <- nrow(dptos)
  nk <- dptos$nk
  yk_barra <- dptos$yk_barra
  s2_k <- dptos$s2_k
  
  nj <- nrow(mpios)
  njk <- mpios$njk
  yjk_barra <- mpios$yjk_barra
  s2_jk <- mpios$s2_jk
  
  #Valores iniciales
  
  sjk <- yjk_barra
  k2jk <- s2_jk
  theta <- yk_barra
  sig_k <- (mpios %>% group_by(cole_cod_depto_ubicacion) %>%
                      summarise(s_k=var(yjk_barra)) %>%
                      mutate(s_k=replace_na(s_k, mean(s_k, na.rm=T))) )$s_k
  mu <- mean(theta)
  tau2 <- var(theta)
  sig2 <- rgamma(1, shape = a_sig/2, rate = b_sig/2)
  alpha_k <- rgamma(1, shape = a_ak, b_ak)
  beta_k <- rgamma(1, shape = a_bk, b_bk)
  k2_k <- rgamma(m, shape = alpha_k/2, rate = beta_k/2)
  
  #Cadena y progreso

  B <- n_burn + n_sams*n_skip
  
  Progress <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (i in 1:B) {
    
    # actualizar parC!metros
    
    sjk <- sample_sjk(njk, nj, k2jk, sig_k, yjk_barra, theta, nk, sjk)
    k2jk <- sample_k2jk(nu_k, nk, njk, nj, sjk,k2_k, s2_jk, yjk_barra, k2jk)
    theta <- sample_theta(nk, m, sig_k, tau2, sjk, mu, deptos, theta)
    sig_k <- sample_sig_k(nu_s, m, nk, sig2, sjk, theta, deptos, sig_k)
    k2k <- sample_k2k(alpha_k, m, nk, nu_k, beta_k, k2jk, deptos, k2k)
    mu <- sample_mu(m, tau2, sig_mu, theta, mu_mu, mu)
    tau2 <- sample_tau2(nu_tau, m, sig_tau, theta, mu, tau2)
    sig2 <- sample_sig2(a_sig, m, nu_s, b_sig, sig_k, sig2)
    beta_k <- sample_beta_k(a_bk, m, alpha_k, b_bk, k2k, beta_k)
    a <- sample_alpha_k(delta, m, i, n_tun, ac, beta_k, k2k, n_burn, alpha_k)
    
    ac <- a$ac
    delta <- a$delta
    n_tun <- a$n_tun
    alpha_k <- a$alpha_k
    
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      LL  <- sum(dnorm(y, mean = rep(sjk,njk), sd = sqrt(rep(k2jk, njk)), log=T))
      fwrite(as.data.table(t(LL)), file = "LL_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sjk)), file = "ZETA_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(k2jk)), file = "KAPPA_JK_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(theta)), file = "THETA_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig_k)), file = "SIGMA_K_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(k2k)), file = "KAPPA_K_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(mu)), file = "MU_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(tau2)), file = "TAU_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(sig2)), file = "SIGMA_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(alpha_k)), file = "ALPHA_iterations-1.txt", append = TRUE, col.names = FALSE)
      fwrite(as.data.table(t(beta_k)), file = "BETA_iterations-1.txt", append = TRUE, col.names = FALSE)
    }
    # progreso
    setTxtProgressBar(Progress, i)
  }
}

# 4. Simulación -----

## 4.1 Ejecución del algoritmo ----

set.seed(7)

tictoc::tic()

M1MC(dptos, mpios, sb11,
     nu_k,nu_s,mu_mu,sig_mu,nu_tau,sig_tau,a_sig,b_sig,a_ak,b_ak,a_bk,b_bk,
     delta, n_tun, ac,
     n_sams, n_burn, n_skip, verbose=T)

tictoc::toc()

# 5. Almacenamiento ----

ZETA<-fread("ZETA_iterations-1.txt")
KAPPA_JK<-fread("KAPPA_JK_iterations-1.txt")
THETA<-fread("THETA_iterations-1.txt")
SIGMA_K<-fread("SIGMA_K_iterations-1.txt")
KAPPA_K<-fread("KAPPA_K_iterations-1.txt")
MU<-fread("MU_iterations-1.txt")
TAU<-fread("TAU_iterations-1.txt")
SIGMA<-fread("SIGMA_iterations-1.txt")
ALPHA<-fread("ALPHA_iterations-1.txt")
BETA<-fread("BETA_iterations-1.txt")
LL <- fread("LL_iterations-1.txt")


M1 <- list(
  ZETA=ZETA, KAPPA_JK=KAPPA_JK, THETA=THETA, SIGMA_K = SIGMA_K,
  KAPPA_K=KAPPA_K, MU = MU, TAU = TAU, SIGMA = SIGMA, ALPHA = ALPHA,
  BETA = BETA, LL = LL
)

save(M1, file="M1.RData")