library(Matrix)
library(dplyr)
library(data.table)
library(tidyr)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
setwd("C:/Users/Lenovo/OneDrive - Universidad Nacional de Colombia/Documentos/Caso Bayesiana 2025 - I")

sourceCpp("samplesfunctions.cpp")



sb11<-read.csv("Datos.txt", sep=";")

sb11 <- sb11 %>% filter(cole_cod_mcpio_ubicacion  !=15401 &
                          cole_cod_mcpio_ubicacion!=94883) %>% 
  arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)

Dep<-read.csv("Departamentales.txt", sep = ";")
Mun<-read.csv("Municipales.txt", sep = ";")

sb11 <- left_join(sb11, Mun, by=c("cole_cod_mcpio_ubicacion"="codmpio"))

sb11 <- left_join(sb11, Dep, by=c("cole_cod_depto_ubicacion"="Cod"))

sb11 <- sb11 %>% drop_na(edu_madre, comp, internet, libros, estrato, etnia,
                         PIB, porc_riesgo, porc_rul, doc_alum, nbi, RISK_VICTIM_2022)

sb11 <- sb11 %>% mutate(PIB = PIB/1e6)

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
            .groups = "drop") %>% 
  mutate(s2_jk=replace_na(s2_jk, mean(s2_jk, na.rm=T)))




fit <-  lm(punt_global~edu_madre+comp+internet+libros+estrato+etnia+
             PIB+porc_riesgo+porc_rul+doc_alum+nbi+RISK_VICTIM_2022 ,data = sb11)

X <- model.matrix(fit)

c <- X[,c(2:14)]
Z <- X[,c(15:17)]
W <- X[,c(18:20)]

e <- ncol(c)
l <- ncol(Z)
d <- ncol(W)

rm(c)
rm(Z)
rm(W)
# Modelo 2

sample_theta<- function(k2jk, njk, sig_beta, sig_E, sig_M, sig_D, beta_0, X, tX, y, theta){
  
  sigma_0 <- diag(c(sig_beta, rep(sig_E,e), rep(sig_M,l), rep(sig_D, d)))
  sigma_inv <- chol2inv(chol(sigma_0))
  sigma <- 1/rep(k2jk, njk)
  
  s1 <- sigma_inv%*%beta_0
  v = chol2inv(chol(sigma_inv + tX %*% (X*sigma) ))
  m = v%*% (s1+tX%*%(y*sigma))
  theta <- rmvnorm(1, mean = m, sigma = v)
  return(theta)
}

sample_k2jk <- function(y,nu_k, nk, njk, nj, theta, X, mpios,k2_k, k2jk){
  
  residuos <- (y-X%*%t(theta))^2
  
  dt <- data.table(sj=residuos, mpios=mpios)
  
  res <- dt[,.(rjk = mean(sj)), by=mpios]
  
  a <- ( nu_k + njk ) / 2
  
  b <- ( nu_k*rep(k2_k,nk) +  njk*res$rjk ) / 2
  
  k2jk <- 1/rgamma(nj, shape = a, rate = b )
  
  return(k2jk)
}

sample_sig_beta <- function(mu_beta, nu_beta, g_beta, theta, sig_beta){
  
  beta <- theta[1]
  
  a = ( nu_beta + 1 ) / 2
  
  b = ( nu_beta*g_beta + sum( (beta-mu_beta)^2 ) ) / 2
  
  sig_beta <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_beta)
}

sample_sig_E <- function(mu_E, nu_E, g_E, theta, e, sig_E){
  
  beta_E <- theta[2:14]
  
  a = ( nu_E + e ) / 2
  
  b = ( nu_E*g_E + t(beta_E-mu_E)%*%(beta_E-mu_E) ) / 2
  
  sig_E <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_E)
}

sample_sig_M <- function(mu_M, nu_M, g_M, theta, l, sig_M){
  
  beta_M <- theta[15:17]
  
  a = ( nu_M + l ) / 2
  
  b = ( nu_M*g_M + t(beta_M-mu_M)%*%(beta_M-mu_M) ) / 2
  
  sig_M <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_M)
}

sample_sig_D <- function(mu_D, nu_D, g_D, theta, d, sig_D){
  
  beta_D <- theta[18:20]
  
  a = ( nu_D + d ) / 2
  
  b = ( nu_D*g_D + t(beta_D-mu_D)%*%(beta_D-mu_D) ) / 2
  
  sig_D <- 1/rgamma(1, shape = a, rate = b)
  
  return(sig_D)
}


sample_k2k <- function(alpha_k, m, nk, nu_k, beta_k, k2jk, deptos, k2k){
  
  dt <- data.table(k2jk=k2jk, deptos=deptos)
  res <- dt[,.(pj =sum(1/k2jk)), by=deptos]
  setorder(res, deptos)
  
  a = ( alpha_k + nk*nu_k ) / 2
  
  b = ( beta_k + nu_k*res$pj ) / 2
  
  k2k <- rgamma(m, shape = a, rate = b)
  
  return(k2k)
}


sample_beta_k <- function(a_bk, m, alpha_k, b_bk, k2k, beta_k){
  
  a = ( 2*a_bk + m*alpha_k ) / 2
  
  b = ( 2*b_bk + sum(k2k) ) / 2
  
  beta_k = rgamma(1, shape = a, rate = b)
  
  return(beta_k)
}


#---------------------------------------------------------------------------

## 2.11 log(p(\alpha_k | resto)) ----

post <- function(alpha_k, m, beta_k, k2k){
  (m*alpha_k/2)*log(beta_k/2) -m*lgamma(alpha_k/2) + (alpha_k/2)*sum(log(k2k)) + (a_ak-1)*log(alpha_k) - b_ak*alpha_k
}

## 2.12 Tasa de aceptación adaptativa ----

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


#Hiperparametros

 
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

M1MC <- function(dptos, mpios, sb11,
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

Rprof("perfil.out")

set.seed(7)

tictoc::tic()

M1MC(dptos, mpios, sb11,
     X, Z, W, e, l, d, fit,
     mu_beta, mu_M, mu_E, mu_D, nu_beta, nu_M, nu_E, nu_D, nu_k, a_ak, b_ak,a_bk, b_bk,
     delta, n_tun, ac,
     n_sams, n_burn, n_skip, verbose=T)

tictoc::toc()

Rprof(NULL)

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

colMeans(sqrt(KAPPA_JK))


Matrix::chol2inv(chol())
plot(LL$V1)
matplot(ALPHA, type = "l")
matplot(BETA_K, type = "l")


plot(colMeans(as.matrix(KAPPA_JK)))

coda::effectiveSize(THETA)
summary(coda::effectiveSize(KAPPA_JK))
summary(coda::effectiveSize(KAPPA_K))
coda::effectiveSize(SIG_B)
coda::effectiveSize(SIG_E)
coda::effectiveSize(SIG_M)
coda::effectiveSize(SIG_D)
coda::effectiveSize(ALPHA)
coda::effectiveSize(BETA_K)


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


mean(sqrt(SIG_B$V1))
mean(sqrt(SIG_E$V1))
mean(sqrt(SIG_M$V1))
mean(SIG_D$V1)
colMeans(THETA)
