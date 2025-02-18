#Autores:


#---------- Juan Andrés Camacho Zaráte (jucamachoz@unal.edu.co)
#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)

# ---------------------------- Cuarto Modelo ------------------------------

library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(metRology)

#1. Tratamiento de datos -------------------------------------------------

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
  select(COLE_COD_DEPTO_UBICACION,COLE_COD_MCPIO_UBICACION,dpto,PUNT_GLOBAL)%>%
  arrange(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION)

muni<-y%>%group_by(COLE_COD_MCPIO_UBICACION)%>%summarise(njk=n(), dpto=sum(as.numeric(dpto))/njk)

X=cbind(X,muni%>%group_by(dpto)%>%summarise(nj=n_distinct(COLE_COD_MCPIO_UBICACION))
        %>%select(nj))


J <- y%>%summarise(nk=n_distinct(dpto))%>%select(nk)%>%as.numeric() #numero de dptos
nj <- X$nj #Numero de municipios por dpto
njk <- muni$njk #Estudiantes munij por municipio
K <- muni%>%summarise(k=n_distinct(COLE_COD_MCPIO_UBICACION))%>%select(k)%>%as.numeric() #numero de municipios
nk <- (y%>%group_by(COLE_COD_DEPTO_UBICACION)%>%summarise(nk=n())%>%select(nk))$nk #estudiantes por dpto

#2. Funciones para samplear los parametros-------------------------

#2.1 Variables auxiliares ---------------------------------------------

sample_sijk <- function(v, k2, epsjk, njk, y, n, sijk)
{
  c = (y - rep(epsjk, njk))^2
  b = (c + v*k2)/2
  
  #Generar las latentes
  
  sijk <- 1/rgamma(n = n, shape = (v+1)/2, rate = b)
  
  return(sijk)
}

#2.2 Medias por municipios -----------------------------------------

sample_epsjk <- function(y, mun, sijk, theta, nj, sig2k, K, epsjk)
{
  dt <- data.frame(c = y/sijk, d = 1/sijk, mun = mun)
  
  sum_c <- (dt%>%group_by(mun)%>%summarise(sum_c=sum(c))%>%select(sum_c))$sum_c
  sum_d <- (dt%>%group_by(mun)%>%summarise(sum_d=sum(d))%>%select(sum_d))$sum_d
  
  m = (sum_c + rep(theta, nj)/rep(sig2k, nj))/(sum_d + 1/rep(sig2k,nj))
  V2 = 1/(sum_d + 1/rep(sig2k,nj))
  
  #Generar las medias por municipio
  
  epsjk <- rnorm(n = K, mean = m, sd = sqrt(V2))
  
  return(epsjk)
}

#2.3 MU ---------------------------------------------------------------

sample_mu <- function(theta, mu0, g20, tau2, J, mu)
{
  m <- (J*mean(theta)/tau2 + mu0/g20)/(J/tau2 + 1/g20)
  V2 <- 1/(J/tau2 + 1/g20)
  
  #Generar mu
  
  mu <- rnorm(n=1, mean = m, sd = sqrt(V2))
  return(mu)
}

#2.4 TAU^2 ---------------------------------------------------------

sample_tau2 <- function(eta0, J, theta, mu, t20, tau2)
{
  a <- (eta0 + J)/2
  b <- (sum((theta-mu)^2) + eta0*t20)/2
  
  #Generar tau^2
  
  tau2 <- 1/rgamma(n = 1, shape = a, rate = b)
  return(tau2)
}

#2.5 SIG^2 ---------------------------------------------------------

sample_sig2 <- function(J, a0, b0, nu, sig2k, sig2)
{
  a <- (a0 + J*nu)/2
  b <- (b0 + nu * sum(1/sig2k))/2   #Generar el sigma
  
  sig2 <- rgamma(n = 1, shape = a, rate = b)
  return(sig2)
}

#2.6 K^2 -----------------------------------------------------------

sample_k2 <- function(v, nk, eps0, sijk, k20, k2)
{
  a = (v*sum(nk) + eps0)/2
  b = (v*sum(1/sijk) + eps0*k20)/2
  
  #Generar k^2
  
  k2 <- rgamma(n= 1, shape = a, rate = b)
  return(k2)
}

#2.7 SIG_K^2 --------------------------------------------------------

sample_sig2k <- function(nk, nu, epsjk, theta, nj, J, sig2, sig2k)
{
  dt <- data.frame(c = (epsjk - rep(theta, nj))^2, dpto = rep(1:J, nj))
  
  res <- dt %>% group_by(dpto) %>% summarise( d = sum(c))
  
  a <- (nj + nu)/2
  b <- (res$d + nu * sig2)/2   #Generar el sigma
  
  sig2k <- 1/rgamma(n = J, shape = a, rate = b)
  return(sig2k)
}

#2.8 THETA ----------------------------------------------------------

sample_theta <- function(epsjk, sig2k, tau2, mu, J, dpto, nj, theta)
{
  dt <- data.frame( epsjk = epsjk/rep(sig2k,nj), dpto = rep(1:J, nj) )
  
  res <- dt%>%group_by(dpto)%>%summarise(c=sum(epsjk))
  
  m = (res$c + mu/tau2)/(nj/sig2k + 1/tau2)
  V2 = 1/(nj/sig2k + 1/tau2)
  
  #Generar las medias por departamento
  
  theta <- rnorm(n = J, mean = m, sd = sqrt(V2) )
  return(theta)
}

#2.9 NU ---------------------------------------------------------------

sample_nu <- function(nus0, sig2, J, sig2k, lam0, nu)
{
  lpnu=0.5*nus0*J*log(0.5*nus0*sig2)-J*lgamma(0.5*nus0)-(0.5*nus0)*sum(log(sig2k))-nus0*(lam0+0.5*sig2*sum(1/sig2k))
  
  prob <- exp(lpnu - max(lpnu))  # Restar max(lpnu) para estabilidad numérica
  
  nu=sample(x=nus0, size=1,  prob = prob)
  return(nu)
}

#3. Muestreador de Gibbs ----------------------------------------------

#3.1 Hiperparametros -----------------------

v <- 3 #Grados de libertad
eps0 <- 1 #epsilon_0
k20 <- 50^2 #kappa_0^2
mu0 <- 250 #mu_0
g20 <- 50^2 #gamma_0^2
eta0 <- 1 #eta_0
t20 <- 50^2 #tau_0^2
lam0 <- 1 #lambda_0
a0 <- 1 #alpha_0
b0 <- 1/50^2 #beta_0
postsijk <- 0 #Para guardar las medias posteriores de las auxliares
nus0 <- 1:50 #rando de valores para nu

#3.2 Muestreador---------------------------

mcM4 <- function (y, J, K, nj, nk, njk, v, eps0, k20, mu0, g20, eta0, t20, lam0, a0, b0, postsijk, nus0, n_sams, n_burn, n_skip, verbose = TRUE) 
{
  n=nrow(y)
  dpto=y$dpto
  mun=y$COLE_COD_MCPIO_UBICACION
  y=y$PUNT_GLOBAL
  
  # número de iteraciones
  
  B <- n_burn + n_sams*n_skip
  ncat <- floor(0.1*B)
  
  # valores iniciales
  
  sig2  <- rgamma(n = 1, shape = a0/2, rate = b0/2)
  k2 <- rgamma(n = 1, shape = eps0/2, rate = eps0*k20/2)
  sijk <- 1/rgamma(n=n, shape =v/2, rate= v*k2/2)
  mu <- rnorm(n = 1, mean = mu0, sd= sqrt(g20))
  theta <- rnorm(n = J, mean = mu, sd = sqrt(t20))
  tau2 <- 1/rgamma(n=1, shape = eta0/2, rate = eta0*t20/2)
  nu <- 1
  sig2k <- 1/rgamma(n = J, shape = nu/2, rate = nu*sig2/2)
  epsjk <- rnorm(n = K, mean = rep(theta,nj), sd = sqrt(rep(sig2k,nj)))
  
  # almacenamiento
  
  for (i in 1:B) {
    
    # actualizar parámetros
    
    sijk <- sample_sijk(v, k2, epsjk, njk, y, n, sijk)
    k2 <- sample_k2(v, nk, eps0, sijk, k20, k2)
    sig2 <- sample_sig2(J, a0, b0, nu, sig2k, sig2)
    theta <- sample_theta(epsjk, sig2k, tau2, mu, J, dpto, nj, theta)
    mu <- sample_mu(theta, mu0, g20, tau2, J, mu)
    tau2<- sample_tau2(eta0, J, theta, mu, t20, tau2)
    epsjk <-sample_epsjk(y, mun, sijk, theta, nj, sig2k, K, epsjk)
    sig2k <- sample_sig2k(nk, nu, epsjk, theta, nj, J, sig2, sig2k)
    nu <- sample_nu(nus0, sig2, J, sig2k, lam0, nu)
  
    # almacenar y log-verosimilitud
    
    if (i > n_burn && (i - n_burn) %% n_skip == 0)  {
      postsijk <- postsijk + sijk/n_sams
      LL  <- sum(dt.scaled(y, v, mean = rep(epsjk,njk), sd = sqrt(k2), log = TRUE))
      write.table(as.data.frame(t(epsjk)), file = "EPSJK_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(theta)), file = "THETA_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(sig2)), file = "SIG2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(tau2)), file = "TAU2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(LL)), file = "LL_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(k2)), file = "K2_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(mu)), file = "MU_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(sig2k)), file = "SIG2K_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(as.data.frame(t(nu)), file = "NU_iterations.txt", append = TRUE, col.names = FALSE, row.names = FALSE)
      
    }
    
    # progreso
    
    if (verbose && i %% ncat == 0)
      cat(sprintf("%.1f%% completado\n", 100*i/B))
  }
  # salida
  
  return(POSTSIJK = postsijk)
}

n_sams=10000
n_burn=10000
n_skip=10

set.seed(123)

M4 <- mcM4(y, J, K, nj, nk, njk, v, eps0, k20, mu0, g20, eta0, t20, lam0, a0, b0, postsijk, nus0, n_sams, n_burn, n_skip, verbose = TRUE) 

fwrite(as.data.table(M4), file = "POSTSIJK_iterations.txt", append = TRUE, col.names = FALSE) #guardar en disco

#Carga de muestras

EPSJK <- fread("EPSJK_iterations.txt")
THETA <- fread("THETA_iterations.txt")
MU <- fread("MU_iterations.txt")
SIG2 <- fread("SIG2_iterations.txt")
TAU2 <- fread("TAU2_iterations.txt")
K2 <- fread("K2_iterations.txt")
NU <- fread("NU_iterations.txt")
SIG2K <- fread("SIG2K_iterations.txt")
POSTSIJK <- fread("POSTSIJK_iterations.txt")
LL <- fread("LL_iterations.txt")

#Compilar todas las muestras en una lista

M4 <- list(EPSJK = EPSJK, THETA = THETA, MU = MU, SIG2 = SIG2, 
           TAU2 = TAU2, K2 = K2, NU = NU, SIG2K = SIG2K, 
           POSTSIJK = POSTSIJK, LL = LL)

save(M4, file = "M4.RData")

load("M4.RData")
