setwd("C:/Users/Lenovo/OneDrive - Universidad Nacional de Colombia/Documentos/Caso Bayesiana 2025 - I")

library(dplyr)
library(sf)
library(ggplot2)
library(readxl)
library(stringr)
library(ggpubr)
library(knitr)
library(kableExtra)
library(mclust)

#0. Carga de bases ----

## 0.1 Datos saber 11 ----

sb11<-read.csv("Datos.txt", sep=";")

## 0.2 Pobreza monetaria ---- 

pobreza_m=read_excel("anexo_pobreza_monetaria_18_departamento.xls",sheet = 2, range = "A16:P40") 
colnames(pobreza_m)[1]="Departamento"

pobreza_m <- pobreza_m %>% select(Departamento,`2018`)

## 0.3 Shape por municipios ----

sh_mpios <- sf::st_read("MGN_MPIO_POLITICO.shp", quiet = T)
sh_mpios = sh_mpios[sh_mpios$DPTO_CCDGO != '88',]
colnames(sh_mpios)[6] <- "MUNI" 

## 0.4 Shape por departamentos ----

sh_deptos <- st_read("MGN_DPTO_POLITICO.shp", quiet = T)
colnames(sh_deptos)[1]="DPTO"
sh_deptos = sh_deptos[sh_deptos$DPTO != '88',]

## 0.5 Shape con info global ----

mundoshp <- st_read("admin00.shp",quiet=TRUE)
mundocol <- mundoshp %>% 
  filter(CNTRY_NAME %in% c("Peru","Brazil","Venezuela","Ecuador","Panama"))
box=st_bbox(sh_deptos)

## 0.6 Cobertura neta secundaria ----

CNS <- read.csv("MEN_ESTADISTICAS_EN_EDUCACION_EN_PREESCOLAR__B_SICA_Y_MEDIA_POR_MUNICIPIO_20250623.csv")
CNS <- CNS %>% filter(AÑO==2022) %>% select(AÑO, MUNI = CÓDIGO_MUNICIPIO, CNS = COBERTURA_NETA_SECUNDARIA)
CNS <- CNS %>% mutate(MUNI = str_pad(MUNI, width = 5, pad = "0"))

## 0.7 Modelos ----

load("M1.RData")
load("M2.RData")


## 0.8 Medidas por departamento (k) ----

dptos <- sb11 %>% group_by(cole_cod_depto_ubicacion) %>% 
  summarise(yk_barra = mean(punt_global), #Media por dptos
            s2_k = var(punt_global), #Var por dptos
            nk = n_distinct(cole_cod_mcpio_ubicacion),  #Nro municipios
            nijk = n()) #Nro de estudiantes)

## 0.9  Medidas por municipio (j) ----

mpios <- sb11 %>% group_by(cole_cod_mcpio_ubicacion, cole_cod_depto_ubicacion) %>% 
  summarise(yjk_barra = mean(punt_global), #Media por mpios
            s2_jk = var(punt_global), #Var por mpios
            njk = n(), #Nro de estudiantes
            .groups = "drop")

## 0.10 Datos de prueba ----

prueba<-read.csv("Prueba.txt", sep=";")

## 0.11 Covariables ----

Dep<-read.csv("Departamentales.txt", sep = ";")
Mun<-read.csv("Municipales.txt", sep = ";")


#2. Mapa de Colombia y Pobreza monetaria ----

departamentos <- data.frame(
  Codigo = c(5, 8, 11, 13, 15, 17, 18, 19, 20, 23, 25, 27, 41, 44, 47, 50, 52, 54, 63, 66, 68, 70, 73, 76, 81, 85, 86, 88, 91, 94, 95, 97, 99),
  Nombre = c("Antioquia", "Atlántico", "Bogotá D.C.", "Bolívar", "Boyacá", "Caldas", "Caquetá", "Cauca", "Cesar", "Córdoba", 
             "Cundinamarca", "Chocó", "Huila", "La Guajira", "Magdalena", "Meta", "Nariño", "Norte de Santander", 
             "Quindío", "Risaralda", "Santander", "Sucre", "Tolima", "Valle del Cauca", "Arauca", "Casanare", "Putumayo", 
             "San Andrés", "Amazonas", "Guainía", "Guaviare", "Vaupés", "Vichada"))

pobreza_m <- inner_join(departamentos, pobreza_m, by =c("Nombre"="Departamento"))
pobreza_m <- pobreza_m %>% mutate(DPTO = str_pad(Codigo, width = 2, pad = "0"))

dptos <- dptos %>% mutate(DPTO = str_pad(cole_cod_depto_ubicacion, width = 2, pad = "0"))

## 2.1 Mapa con la media del puntaje global ----

m1 <- inner_join(x = sh_deptos, y = dptos, by = c("DPTO")) %>% 
      select(DPTO, yk_barra, geometry) %>%
      ggplot() +
      geom_sf(data=mundocol, col="white")+
      geom_sf(aes(fill = yk_barra), size = 0.125, color = "#1E1E1E", linetype="solid") +
      coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
      scale_fill_gradient(
        high = "#00FF00",   
        low = "#A50021"
      )  +
      labs(fill="Media del \n Puntaje \n Global")+
      theme(
        plot.title = element_text(hjust = 0.5,),
        panel.background=element_rect(fill="lightblue")
      )

## 2.2 Mapa de pobreza monetaria ----

m2 <- left_join(x = sh_deptos, y = pobreza_m, by = c("DPTO"))%>%
      select(DPTO ,`2018`, geometry)%>%
      ggplot() +
      geom_sf(data=mundocol, col="white")+
      geom_sf(aes(fill = `2018`), size = 0.125, color = "#1E1E1E", linetype="solid") +
      coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
      scale_fill_gradient(
        low = "#00FF00",   
        high = "#A50021",   
        na.value = "#EEE9E9"
      )  +
      labs(fill="IPM \n en \n 2018")+
      theme(
        plot.title = element_text(hjust = 0.5,),
        panel.background=element_rect(fill="lightblue")
      )

## 2.3 Unión de los mapas

m3 <- ggarrange(m1,m2, ncol = 2, nrow = 1, widths = c(1,1), heights = c(1,1))

m3 <- annotate_figure(m3, top = text_grob("Colombia", face = "bold", size = 14),
                left = text_grob("Latitud", rot = 90),
                bottom = text_grob("Longitud"))
m3

# 3. Mapa de Colombia y cobertura neta ----

## 3.1 Mapa de los municipios ----

mpios <- mpios %>% mutate(MUNI = str_pad(cole_cod_mcpio_ubicacion, width = 5, pad = "0"))

m1 <- left_join(x = sh_mpios, y = mpios, by = c("MUNI")) %>% 
      select(MUNI, yjk_barra, geometry) %>%
      ggplot() +
      geom_sf(data=mundocol, col="white")+
      geom_sf(aes(fill = yjk_barra), size = 0.0125, color = "gray30", linetype="solid") +
      coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
      scale_fill_gradient(
        high = "#85CF00",  
        low = "#FFB2B6",   
        na.value = "#EEE9E9"
      ) +
      labs(fill="Media del \n Puntaje \n Global")+
      theme(
        plot.title = element_text(hjust = 0.5,),
        panel.background=element_rect(fill="lightblue")
      )

## 3.2 Mapa cobertura neta ----

m2 <- left_join(x = sh_mpios, y = CNS, by = c("MUNI"))%>%
      select(MUNI ,CNS, geometry)%>%
      ggplot() +
      geom_sf(data=mundocol, col="white")+
      geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
      coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
      scale_fill_gradient(
        high = "#85CF00",  
        low = "#FFB2B6",   
        na.value = "#EEE9E9"
      )  +
      labs(fill="CNS \n en 2022")+
      theme(
        plot.title = element_text(hjust = 0.5,),
        panel.background=element_rect(fill="lightblue")
      )

## 3.3 Unión de los mapas ----

m4 <- ggarrange(m1,m2, ncol = 2, nrow = 1, widths = c(1,1), heights = c(1,1))

m4 <- annotate_figure(m4, top = text_grob("Colombia", face = "bold", size = 14),
                      left = text_grob("Latitud", rot = 90),
                      bottom = text_grob("Longitud"))
m4
# 4. Análisis descriptivo ----

col <- RColorBrewer::brewer.pal(9, "Set1")[1:9]

boxplot(sb11_2$punt_global~sb11_2$edu_madre,
        col = NA, outline = FALSE, ylim = range(sb11_2$punt_global),
        ylab = "Puntaje Global", main = expression(italic("Puntaje por educación de la madre")),
        names=c("Sin Educación Superior", "Educación superior"),
        border = c(col[1], col[2]),
        boxwex = 0.4,
        xlab="Nivel educativo de la madre")
points(x = jitter(as.numeric(sb11_2$edu_madre) +1, factor = 6),
       y = sb11_2$punt_global,
       pch = 16, col = adjustcolor("blue", alpha.f = 0.3))


boxplot(sb11_2$punt_global~sb11_2$estrato,
        col = NA, outline = FALSE, ylim = range(sb11_2$punt_global),
        ylab = "Puntaje Global", main = expression(italic("Puntaje por estrato")),
        border = col[1:7],
        names = c(1,2,3,4,5,6,"Sin estrato"),
        boxwex = 0.4,
        xlab="Estrato")

boxplot(sb11_2$punt_global~sb11_2$libros,
        col = NA, outline = FALSE, ylim = range(sb11_2$punt_global),
        ylab = "Puntaje Global", main = expression(italic("")),
        border = col[1:4],
        names = c("0 a 10","11 a 25","26 a 100","+100"),
        boxwex = 0.4,
        xlab="Número de libros en casa")

boxplot(sb11_2$punt_global~sb11_2$etnia,
        col = NA, outline = FALSE, ylim = range(sb11_2$punt_global),
        ylab = "Puntaje Global", main = expression(italic("")),
        border = col[1:2],
        names = c("No","Sí"),
        boxwex = 0.4,
        xlab="Pertenencia étnica")

boxplot(sb11_2$punt_global~sb11_2$internet,
        col = NA, outline = FALSE, ylim = range(sb11_2$punt_global),
        ylab = "Puntaje Global", main = expression(italic("")),
        border = col[1:2],
        names = c("No","Sí"),
        boxwex = 0.4,
        xlab="Acceso a internet en casa")

sb11_2$

sb11_2$edu_madre

points(x = jitter(rep(1, m + 4), factor = 6), y = CVMC1,
       pch = 16, col = adjustcolor(4, 0.3))


# 7. Verosimilitudes ----

png("log_verosimilitud.png", width = 1800, height = 600, res = 150) 

par(mfrow=c(1,2))

min 

lim_y <- range(c(M1$LL$V1, M2$LL$V1))

# Gráfico para el Modelo 1
plot(x = 1:10000, y = M1$LL$V1, type = "p", pch = ".", col = "#7497CD",
     cex.axis = 0.8, cex = 1.5, main = "Modelo 1", xlab = "Iteración",
     ylab = "Log-verosimilitud", ylim = lim_y)
abline(h = mean(M1$LL$V1), col = "darkblue", lty = 1, cex = 2)

# Gráfico para el Modelo 2
plot(x = 1:10000, y = M2$LL$V1, type = "p", pch = ".", col = "#FF8288",
     cex.axis = 0.8, cex = 1.5, main = "Modelo 2", xlab = "Iteración",
     ylab = "Log-verosimilitud", ylim = lim_y)
abline(h = mean(M2$LL$V1), col = "#4B112D", lty = 1, cex = 2)


dev.off()

# 8. Tamaños efectivos de muestra ----

## 8.1 Primer modelo ----

efsize1 <- rbind( summary(coda::effectiveSize(M1$ZETA)),
                  summary(coda::effectiveSize(M1$KAPPA_JK)),
                  summary(coda::effectiveSize(M1$THETA)),
                  summary(coda::effectiveSize(M1$SIGMA_K)),
                  summary(coda::effectiveSize(M1$KAPPA_K)),
                  summary(coda::effectiveSize(M1$MU)),
                  summary(coda::effectiveSize(M1$TAU)),
                  summary(coda::effectiveSize(M1$SIGMA)),
                  summary(coda::effectiveSize(M1$ALPHA)),
                  summary(coda::effectiveSize(M1$BETA)))

p1 <- c("$\\zeta_{j,k}$", "$\\kappa_{j,k}^2$",
        "$\\theta_{k}$", "$\\sigma_{k}^2$","$\\kappa_{k}^2$",
        "$\\mu$", "$\\tau^2$", "$\\sigma^2$",
        "$\\alpha_{\\kappa}$", "$\\beta_{\\kappa}$")

kable(cbind(p1, round(efsize1,0)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
      kable_styling(position = "center", latex_options = "HOLD_position")

## 8.2 Segundo modelo ----

efsize2 <- rbind(summary(coda::effectiveSize(M2$THETA$V1)),
                summary(coda::effectiveSize(M2$THETA[,c(2:14)])),
                summary(coda::effectiveSize(M2$THETA[,c(15:17)])),
                summary(coda::effectiveSize(M2$THETA[,c(18:20)])),
                summary(coda::effectiveSize(M2$KAPPA_JK)),
                summary(coda::effectiveSize(M2$KAPPA_K)),
                coda::effectiveSize(M2$SIG_B),
                coda::effectiveSize(M2$SIG_E),
                coda::effectiveSize(M2$SIG_M),
                coda::effectiveSize(M2$SIG_D),
                coda::effectiveSize(M2$ALPHA),
                coda::effectiveSize(M2$BETA_K))

p2 <- c("$\\beta$",
        "$\\beta^E$","$\\beta^M$","$\\beta^D$",
        "$\\kappa_{j,k}^2$",
        "$\\kappa_{k}^2$",
        "$\\sigma_{\\beta}^2$", "$\\sigma_{E}^2$", "$\\sigma_{M}^2$", "$\\sigma_{D}^2$",
        "$\\alpha_{\\kappa}$","$\\beta_{\\kappa}$")

kable(cbind(p2, round(efsize2,0)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                                                 "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
      kable_styling(position = "center", latex_options = "HOLD_position")



# 9. Error estándar de Montecarlo ----

## 9.1 Primer modelo ----

semc1 <- rbind( summary(apply(M1$ZETA, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$ZETA))),
                summary(apply(M1$KAPPA_JK, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$KAPPA_JK))),
                summary(apply(M1$THETA, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$THETA))),
                summary(apply(M1$SIGMA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$SIGMA_K))),
                summary(apply(M1$KAPPA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$KAPPA_K))),
                summary(sd(M1$MU$V1)/sqrt(coda::effectiveSize(M1$MU))),
                summary(sd(M1$TAU$V1)/sqrt(coda::effectiveSize(M1$TAU))),
                summary(sd(M1$SIGMA$V1)/sqrt(coda::effectiveSize(M1$SIGMA))),
                summary(sd(M1$ALPHA$V1)/sqrt(coda::effectiveSize(M1$ALPHA))),
                summary(sd(M1$BETA$V1)/sqrt(coda::effectiveSize(M1$BETA))))

kable(cbind(p1, round(semc1,3)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                                                 "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

## 9.2 Segundo Modelo ----

semc2 <- rbind(summary(sd(M2$THETA$V1)/sqrt(coda::effectiveSize(M2$THETA$V1))),
               summary(apply(M2$THETA[,c(2:14)], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[,c(2:14)]))),
               summary(apply(M2$THETA[,c(15:17)], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[,c(15:17)]))),
               summary(apply(M2$THETA[,c(18:20)], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[,c(18:20)]))),
               summary(apply(M2$KAPPA_JK, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$KAPPA_JK))),
               summary(apply(M2$KAPPA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$KAPPA_K))),
               sd(M2$SIG_B$V1)/sqrt(coda::effectiveSize(M2$SIG_B)),
               sd(M2$SIG_E$V1)/sqrt(coda::effectiveSize(M2$SIG_E)),
               sd(M2$SIG_M$V1)/sqrt(coda::effectiveSize(M2$SIG_M)),
               sd(M2$SIG_D$V1)/sqrt(coda::effectiveSize(M2$SIG_D)),
               sd(M2$ALPHA$V1)/sqrt(coda::effectiveSize(M2$ALPHA)),
               sd(M2$BETA_K$V1)/sqrt(coda::effectiveSize(M2$BETA_K)))

kable(cbind(p2, round(semc2,3)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                                                 "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")


sd(M2$SIG_B$V1)/sqrt(coda::effectiveSize(M2$SIG_B))

# 10. Criterios de informacion ----

## 10.1 DIC ----
### Primer modelo ----
njk <- mpios$njk
y <- sb11$punt_global
LL1        <- M1$LL$V1
zeta_hat  <- colMeans(M1$ZETA)
kappa_hat <- colMeans(M1$KAPPA_JK)
lpy_m1     <- sum(dnorm(x = y, mean = rep(zeta_hat, njk), sd = sqrt(rep(kappa_hat,njk)), log = TRUE))
pDIC_m1    <-  2 * (lpy_m1 - mean(LL1))
dic_m1     <- -2 * lpy_m1 + 2 * pDIC_m1

### Segundo modelo ----
library(tidyr)

sb11_2 <- left_join(sb11, Mun, by=c("cole_cod_mcpio_ubicacion"="codmpio"))

sb11_2 <- left_join(sb11_2, Dep, by=c("cole_cod_depto_ubicacion"="Cod"))

sb11_2 <- sb11_2 %>% drop_na(edu_madre, comp, internet, libros, estrato, etnia,
                         PIB, porc_riesgo, porc_rul, doc_alum, nbi, RISK_VICTIM_2022)

sb11_2 <- sb11_2 %>% mutate(PIB = PIB/1e6)

sb11_2 <- sb11_2 %>% mutate(
  edu_madre = as.factor(edu_madre),
  comp = as.factor(comp),
  internet = as.factor(internet),
  libros = as.factor(libros),
  estrato = as.factor(estrato),
  etnia = as.factor(etnia)
)

fit <-  lm(punt_global~edu_madre+comp+internet+libros+estrato+etnia+
             PIB+porc_riesgo+porc_rul+doc_alum+nbi+RISK_VICTIM_2022 ,data = sb11_2)

X <- model.matrix(fit)

mpios_2 <- sb11_2 %>% group_by(cole_cod_mcpio_ubicacion, cole_cod_depto_ubicacion) %>% 
  summarise(yjk_barra = mean(punt_global), #Media por mpios
            s2_jk = var(punt_global), #Var por mpios
            njk = n(), #Nro de estudiantes
            .groups = "drop")

njk <- mpios_2$njk

y <- sb11_2$punt_global
LL2        <- M2$LL$V1
zeta_hat  <- colMeans(M2$THETA)
zeta_hat <- as.matrix(zeta_hat)
zeta_hat <- X%*%(zeta_hat)
kappa_hat <- colMeans(M2$KAPPA_JK)
lpy_m2     <- sum(dnorm(x = y, mean = zeta_hat, sd = sqrt(rep(kappa_hat,njk)), log = TRUE))
pDIC_m2    <-  2 * (lpy_m2 - mean(LL2))
dic_m2     <- -2 * lpy_m2 + 2 * pDIC_m2

## 10.2 WAIC ----
### Modelo 1 ----
mpios <- mpios %>% mutate(id = row_number())

g <- sb11 %>%
  left_join(mpios, by = c("cole_cod_mcpio_ubicacion", "cole_cod_depto_ubicacion")) %>%
  pull(id)
y <- sb11$punt_global
n <- length(y)
lppd_m1  <- 0
pWAIC_m1 <- 0
tictoc::tic()
for (i in 1:n) {
  # lppd
  tmp1    <- dnorm(x = y[i], mean = M1$ZETA[, g[i]], sd = sqrt(M1$KAPPA_JK[,g[i]]))
  lppd_m1 <- lppd_m1 + log(mean(tmp1))
  
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = M1$ZETA[, g[i]], sd = sqrt(M1$KAPPA_JK[,g[i]]), log = T)
  pWAIC_m1 <- pWAIC_m1 + 2 * (log(mean(tmp1)) - mean(tmp2))
}
tictoc::toc()
waic_m1 <- -2 * lppd_m1 + 2 * pWAIC_m1

### Modelo 2 ----
mpios_2 <- mpios_2 %>% mutate(id = row_number())

g <- sb11_2 %>%
  left_join(mpios_2, by = c("cole_cod_mcpio_ubicacion", "cole_cod_depto_ubicacion")) %>%
  pull(id)
y <- sb11_2$punt_global
n <- length(y)
lppd_m2  <- 0
pWAIC_m2 <- 0
B<-10000

library(parallel)
library(pbapply)  

# Crear clúster con el número de núcleos disponibles
cores <- detectCores() - 1  # Usa todos los núcleos menos 1 para no sobrecargar

# Iniciar clúster
cl <- makeCluster(cores)

# Exportar funciones y objetos al clúster
clusterExport(cl, varlist = c("X", "y", "g", "M2", "B"))

WAIC <- function(i){
  # lppd
  mu_samples <- rowSums(M2$THETA * matrix(rep(X[i, ], B), ncol = ncol(X), byrow = TRUE))
  tmp1    <- dnorm(x = y[i], mean = mu_samples, sd = sqrt(M2$KAPPA_JK[,g[i]]))
  lppd_i <- log(mean(tmp1))
  
  # pWAIC
  tmp2 <- dnorm(x = y[i], mean = mu_samples, sd = sqrt(M2$KAPPA_JK[,g[i]]), log = T)
  pWAIC_i <-  2 * (log(mean(tmp1)) - mean(tmp2))
  
  return(c(lppd_i = lppd_i, pWAIC_i = pWAIC_i))
}

pboptions(type = "txt")  # Puedes usar "timer" o "win" si estás en Windows
resultados <- pblapply(1:n, WAIC, cl = cl)

# Detener clúster
stopCluster(cl)

# Procesar resultados
resultados_mat <- do.call(rbind, resultados)
lppd_m2 <- sum(resultados_mat[, "lppd_i"])
pWAIC_m2 <- sum(resultados_mat[, "pWAIC_i"])
WAIC_m2 <- -2 * (lppd_m2 - pWAIC_m2)

  ## 10.3 Tabla ----

res <- rbind(c(dic_m1, waic_m1), c(dic_m2, waic_m2))
rownames(res) <- c("M1", "M2")

kable(res, col.names = c("DIC", "WAIC"), row.names = T,format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

#11. Ranking Bayesiano ----

est <- cbind(apply(M1$THETA, MARGIN=2, FUN=function(x) quantile(x, probs = 0.025)),
             colMeans(M1$THETA),
             apply(M1$THETA, MARGIN=2, FUN=function(x) quantile(x, probs = 0.975)))

est <- as.data.frame(est)
est <- cbind(est, departamentos$Nombre[which(departamentos$Codigo!=88)])
colnames(est) <- c("q1", "mean", "q2", "DPTO")

est <- est %>% mutate(DPTO = case_when(
                       DPTO == "Valle del Cauca" ~ "V. del Cauca",
                       DPTO == "Norte de Santander" ~ "N. Santander",
                       .default = DPTO
))

est <- est %>% arrange(mean)

colores = rep('#838B83', 32)
colores[which(est[,"q2"] < 250)] = '#FF0000' 
colores[which(est[,"q1"] > 250) ] = '#00FF00'

png("ranking_puntaje.png", width = 600, height = 1000, res = 150)

par(mfrow = c(1, 1), mar = c(4, 5, 1.5, 1), mgp = c(2.5, 0.75, 0))

plot(x = c(180, 280), y = c(1, 32), type = "n",
     xlab = "Puntaje", ylab = "",
     main = expression(italic("Ranking Bayesiano")),
     yaxt = "n")

abline(h = 1:32, col = "lightgray", lwd = 1)
abline(v = 250, col = "gray", lwd = 3)
axis(side = 2, at = 1:32, labels = est$DPTO,
     las = 2, cex.axis = 0.7)
for (l in 1:32) {
  segments(x0 = est[l,'q2'], y0 = l, x1 = est[l,'q1'], y1 = l, col = colores[l], lwd=2)
  lines(x = est[l,'mean'], y = l, type = "p", pch = 16, cex = 0.8, col = colores[l])
}

#dev.off()

#12. Agrupamiento ----

indicadora_cluster <- function(cluster1, cluster2) {
  ifelse(cluster1 == cluster2, 1, 0)} #Para hallar la probabilidad de que dos deptos pertenezcan a un mismo grupo

clust=function(theta, sigma, n.clust){
  
  # Inicializar una lista para almacenar los clusters
  Omega <- vector("list", nrow(theta))
  Progress <- txtProgressBar(min = 1, max = nrow(theta), style = 3)
  # Iterar por cada fila
  for (i in 1:nrow(theta)) {
    fila_data <- matrix(cbind(theta[i, ], sigma[i,]), ncol = 2)
    resultado_kmeans <- kmeans(fila_data, centers = n.clust)
    Omega[[i]] <- resultado_kmeans$cluster
    setTxtProgressBar(Progress, i)
  }
  return(Omega)  # Retornar la lista de clusters
}

set.seed(7)
tictoc::tic()
cl <- clust(M1$THETA, M1$SIGMA_K, 5)
tictoc::toc()

cor<-matrix(0, ncol=ncol(M1$THETA), nrow=ncol(M1$SIGMA_K))

tictoc::tic()

Progress = txtProgressBar(min = 2, max = nrow(M1$THETA)+1, style = 3)

for (i in 1:10000)
{
  cor <- cor + (1/10000)*(outer(cl[[i]],cl[[i]],indicadora_cluster))
  setTxtProgressBar(Progress,i)
}

tictoc::toc()

departamentos <- departamentos %>% mutate(Nombre = case_when(
  Nombre == "Valle del Cauca" ~ "V. del Cauca",
  Nombre == "Norte de Santander" ~ "N. Santander",
  .default = Nombre
))

colnames(cor)<-departamentos$Nombre[which(departamentos$Codigo!=88)]
rownames(cor)<-colnames(cor)

orden_ranking <- list(rev(est$DPTO))

cor <- cor[unlist(orden_ranking),unlist(orden_ranking)]

cor <- Mclust(cor, G=5)

grupos_finales <- cor$classification

df <- data.frame(
  depto = names(grupos_finales),
  grupo = grupos_finales
)

png("matriz_incidencia_dpto.png", width = 1000, height = 1100, res = 150)

par(mfrow = c(1, 1), mar = c(4, 5, 3.2, 1), mgp = c(2.5, 0.75, 0))

corrplot::corrplot(corr = cor$data, 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color",
                   tl.pos = "lt", 
                   tl.cex = 0.8,
                   tl.col = 'black')
title(main = 'Matriz de incidencia por departamentos', line = 1.8, cex.main = 1.6)

#dev.off()

df <- inner_join(departamentos, df,  by = c("Nombre"="depto"))

df <- df %>% mutate(DPTO=str_pad(Codigo, 2, pad = "0"))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

inner_join(x = sh_deptos, y = df, by = c("DPTO")) %>% 
  select(DPTO, grupo) %>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = as.factor(grupo)), size = 0.125, color = "gray30", linetype="solid") +
  scale_fill_manual(values = c("#CD950C", "#C590B3", "#E9967A", "#FFEBCD", "#66CDAA")) +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  labs(x="Longitud",y="Latitud",title="Agrupamiento por Departamentos",fill="Cluster")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )


# 13. Prediccion pobreza monetaria ----

THETA <- M1$THETA
colnames(THETA)<-dptos$DPTO

predictivas <- THETA[,pobreza_m$DPTO, with = F]

resto <- THETA[,setdiff(colnames(THETA), pobreza_m$DPTO), with = F]

Predecir <- function(predictivas, resto,pobreza_m, predicciones){
  Progress = txtProgressBar(min = 1, max = nrow(resto), style = 3)
  predicciones <- matrix(NA, ncol = ncol(resto), nrow = nrow(resto))
  pm <- pobreza_m$`2018`
  for (i in 1:nrow(resto)){
  faltante <- as.data.frame(t(resto[i,])) 
  colnames(faltante) <- "Media"
  
  X <- cbind(t(predictivas[i,]),pm)
  colnames(X)<-c("Media", "PM")
  X <- as.data.frame(X)
  
  fit <- lm(PM ~ Media, data = X)
  
  predicciones[i,] <- predict(fit, faltante)
  setTxtProgressBar(Progress,i)
  }
  colnames(predicciones) <- colnames(resto)
  predicciones <- as.data.frame(predicciones)
  return(predicciones)
}

  
predicciones<-Predecir(predictivas, resto, pobreza_m, predicciones)  

r_predic <- cbind(apply(predicciones, MARGIN = 2, FUN = function(x)
  quantile(x, 0.025)), colMeans(predicciones),
  apply(predicciones, MARGIN = 2, FUN = function(x)
    quantile(x, 0.975))
)

r_predic <- as.data.frame(r_predic)

r_predic <- r_predic %>% arrange(desc(V2))

r_predic <- cbind(r_predic, DPTO=rownames(r_predic))

r_predic <- r_predic %>% mutate(DPTO=as.numeric(DPTO))

r_predic <- inner_join(r_predic, departamentos, by=c("DPTO"="Codigo"))

r_predic <- r_predic %>% mutate(DPTO=as.factor(DPTO))

## 13.1 Tabla de predicciones ----

kable(r_predic[,c(5,1:3)], col.names = c("Departamento","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

## 13.2 Mapa ----

left_join(x = sh_deptos, y = r_predic, by = c("DPTO"))%>%
  select(DPTO ,V2, geometry)%>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = V2), size = 0.125, color = "#1E1E1E", linetype="solid") +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    low = "#EEA2AD",   
    high = "#A50021",  
    na.value = "#EEE9E9"
  )  +
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Estimación \n del IPM")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

# 14. Medias por municipio ----

## 14.1 Ranking Bayesiano ----

est <- cbind(apply(M1$ZETA, MARGIN=2, FUN=function(x) quantile(x, probs = 0.025)),
             colMeans(M1$ZETA),
             apply(M1$ZETA, MARGIN=2, FUN=function(x) quantile(x, probs = 0.975)))

est <- as.data.frame(est)

est <- cbind(est, MUNI=mpios$cole_cod_mcpio_ubicacion)
est <- est %>% mutate(MUNI=str_pad(MUNI, width = 5, pad="0"))

colnames(est) <- c("q1", "mean", "q2", "MUNI")

est <- est %>% arrange(desc(mean))

sh_mpios_unico <- sh_mpios %>% distinct(MUNI, .keep_all = TRUE)
rank100 <- head(est,60)
rank100 <- inner_join(rank100, sh_mpios_unico, by=c("MUNI")) %>% 
          select(everything(rank100), MPIO_CNMBR)

colnames(rank100)[5] <- "NOMBRE"

rank100 <- rank100 %>% mutate(NOMBRE=case_when(
  NOMBRE == "SAN ANTONIO DEL TEQUENDAMA" ~ "S.A. DEL TEQUENDAMA",
  NOMBRE == "SANTA ROSA DE VITERBO" ~ "S.R. DE VITERBO",
  .default = NOMBRE
))

rank100 <- rank100 %>% arrange(mean)

colores = rep('#838B83', nrow(rank100))
colores[which(rank100[,"q2"] < 250)] = '#FF0000' 
colores[which(rank100[,"q1"] > 250) ] = '#00FF00'

#png("ranking_puntaje2.png", width = 600, height = 1000, res = 150)

par(mfrow = c(1, 1), mar = c(4, 5, 1.5, 1), mgp = c(2.5, 0.75, 0))

plot(x = c(240, 310), y = c(1, nrow(rank100)), type = "n",
     xlab = "Puntaje", ylab = "",
     main = expression(italic("Ranking Bayesiano (60 primeros)")),
     yaxt = "n")

abline(h = 1:nrow(rank100), col = "lightgray", lwd = 1)
abline(v = 250, col = "gray", lwd = 3)
axis(side = 2, at = 1:nrow(rank100), labels = rank100$NOMBRE,
     las = 2, cex.axis = 0.4)
for (l in 1:nrow(rank100)) {
  segments(x0 = rank100[l,'q2'], y0 = l, x1 = rank100[l,'q1'], y1 = l, col = colores[l], lwd=2)
  lines(x = rank100[l,'mean'], y = l, type = "p", pch = 16, cex = 0.5, col = colores[l])
}

#dev.off( )

## 14.2 Matriz ----

clust1=function(zeta, n.clust){
  
  # Inicializar una lista para almacenar los clusters
  Omega <- vector("list", nrow(zeta))
  Progress <- txtProgressBar(min = 1, max = nrow(zeta), style = 3)
  # Iterar por cada fila
  for (i in 1:nrow(zeta)) {
    fila_data <- matrix(zeta[i, ])
    resultado_kmeans <- kmeans(fila_data, centers = n.clust)
    Omega[[i]] <- resultado_kmeans$cluster
    setTxtProgressBar(Progress, i)
  }
  return(Omega)  # Retornar la lista de clusters
}

ZETA <- M1$ZETA

set.seed(7)
tictoc::tic()
cl <- clust1(ZETA, 8)
tictoc::toc()

cor<-matrix(0, ncol=ncol(ZETA), nrow=ncol(ZETA))

tictoc::tic()

Progress = txtProgressBar(min = 2, max = nrow(M1$THETA)+1, style = 3)

for (i in 1:10000)
{
  cor <- cor + (1/10000)*(outer(cl[[i]],cl[[i]],indicadora_cluster))
  setTxtProgressBar(Progress,i)
}

tictoc::toc()

colnames(cor) <- str_pad(mpios$cole_cod_mcpio_ubicacion,
                         width = 5, pad = "0")

rownames(cor) <- colnames(cor)

nuevo_orden <- list(est$MUNI)

cor <- cor[unlist(nuevo_orden),unlist(nuevo_orden)]

png("matriz_incidencia_muni.png", width = 1000, height = 1100, res = 150)

corrplot::corrplot(corr = cor, 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color",
                   tl.pos = "n", 
                   tl.cex = 0.8,
                   tl.col = 'black')
title(main = 'Matriz de incidencia por municipios', line = 1.8, cex.main = 1.6)

dev.off()

## 14.3 K-means ----

est <- est %>% mutate(cluster=kmeans(mean,8)$cluster)

  left_join(x = sh_mpios, y = est, by = c("MUNI"="MUNI")) %>% 
    select(MUNI, cluster, geometry) %>%
    ggplot() +
    geom_sf(data=mundocol, col="white")+
    geom_sf(aes(fill = as.factor(cluster)), size = 0.0125, color = "gray30", linetype="solid") +
    coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
    scale_fill_manual(values=c("#F5F5DC", "#AB82FF", "#008B8B", "#BDB76B", "#E9967A","#228B22","#D02090","#FFFF0B")
    ) +  
    labs(x="Longitud",y="Latitud",title="Agrupamiento por municipios",fill="Clusters")+
    theme(
      plot.title = element_text(hjust = 0.5,),
      panel.background=element_rect(fill="lightblue")
    )

# 15. Predicción por municipio ----
  
ZETA <- M1$ZETA
colnames(ZETA)<-str_pad(mpios$cole_cod_mcpio_ubicacion,
                        width = 5, pad = "0")
  
predictivas <- CNS$MUNI[CNS$MUNI %in% colnames(ZETA)]

predictivas <- ZETA[,predictivas]
  
resto <- ZETA[, setdiff(colnames(ZETA), CNS$MUNI), with = F]

cobertura <- CNS%>%filter(MUNI %in% colnames(predictivas))

Predecir1 <- function(predictivas, resto,CNS, predicciones){
  Progress = txtProgressBar(min = 1, max = nrow(resto), style = 3)
  predicciones <- matrix(NA, ncol = ncol(resto), nrow = nrow(resto))
  CNS <- cobertura$CNS
  for (i in 1:nrow(resto)){
    faltante <- as.data.frame(t(resto[i,])) 
    colnames(faltante) <- "Media"
    
    X <- cbind(t(predictivas[i,]),CNS)
    colnames(X)<-c("Media", "CNS")
    X <- as.data.frame(X)
    
    fit <- lm(CNS ~ Media, data = X)
    
    predicciones[i,] <- predict(fit, faltante)
    setTxtProgressBar(Progress,i)
  }
  colnames(predicciones) <- colnames(resto)
  predicciones <- as.data.frame(predicciones)
  return(predicciones)
}

predicciones<-Predecir1(predictivas, resto, CNS, predicciones)  

r_predic <- cbind(apply(predicciones, MARGIN = 2, FUN = function(x)
  quantile(x, 0.025)), colMeans(predicciones),
  apply(predicciones, MARGIN = 2, FUN = function(x)
    quantile(x, 0.975))
)

r_predic <- as.data.frame(r_predic)

r_predic <- r_predic %>% arrange(desc(V2))

r_predic <- cbind(r_predic, MUNI=rownames(r_predic))

## 15.1 Tabla resumen ----

kable(cbind(c("Belén de Bajirá","Mapiripana"),r_predic[,c(1:3)]), col.names = c("Municipio","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

## 15.2 Mapa ----

m1<-left_join(x = sh_mpios, y = CNS , by = c("MUNI"="MUNI"))%>%
  filter(DPTO_CCDGO=="27")%>%
  select(MUNI ,CNS, geometry)%>%
  ggplot() +
  geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
  scale_fill_gradient(
    high = "#3B8C4D",   # Color claro para valores bajos
    low = "#D4C95F",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  ) +
  labs(title="Chocó",fill="Cobertura \n Neta \n Secundaria \n en 2022")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="gray90")
  )

m2<-left_join(x = sh_mpios, y = CNS , by = c("MUNI"="MUNI"))%>%
  filter(DPTO_CCDGO=="94")%>%
  select(MUNI ,CNS, geometry)%>%
  ggplot() +
  geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
  scale_fill_gradient(
    high = "#3B8C4D",   # Color claro para valores bajos
    low = "#D4C95F",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  ) +
  labs(title="Guainía",fill="Cobertura \n Neta \n Secundaria \n en 2022")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="gray90")
  )

m3 <- ggarrange(m1,m2, ncol = 2, nrow = 1, widths = c(1,1), heights = c(1,1))

m3 <- annotate_figure(m3, left = text_grob("Latitud", rot = 90),
                      bottom = text_grob("Longitud"))

m3

# 16. Análisis de los coeficientes ----

## 16.1 Intercepto ----

beta_01 <- c(quantile(M2$THETA$V1, 0.025),mean(M2$THETA$V1),quantile(M2$THETA$V1, 0.975)) 

kable(matrix(c("$\\beta_0$", round(beta_01,3)), nrow = 1), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")


## 16.2 Covariables individuales ----

beta_E1 <- cbind(apply(M2$THETA[,c(2:14)] ,MARGIN=2, FUN = function(x) quantile(x, 0.025)),
             colMeans(M2$THETA[,c(2:14)]),
             apply(M2$THETA[,c(2:14)] ,MARGIN=2, FUN = function(x) quantile(x, 0.975))) 

### 16.2.1 Var dicotomicas ----

kable(cbind(paste0("$\\beta_", c(1:3,12),"$"), round(beta_E1[c(1:3,13),],3)), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")


png("dicotomicas.png", width = 1200, height = 800, res = 150)

par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

colo <- c("blue","black")[as.numeric((beta_E1[c(1:3,13),1] < 0) & (0 < beta_E1[c(1:3,13),3]))+1]

plot(NA, NA, 
     xlab = "Var. Dicotómicas", 
     ylab = "Efecto tratamiento", 
     main = "Modelo 2", 
     xlim = c(0, 5), ylim=range(beta_E1[c(1:3,13),]),
     cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:4, 
     labels = c("Madre con \n educación superior", 
                "Tiene \n computador", 
                "Tiene \n internet", 
                "Tiene \n etnia"), 
     cex.axis = 0.65)
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:4, col = "gray95", lwd = 1, lty = 2)
a <- c(1:3,13)
for (j in 1:4) {
  segments(x0 = j, y0 = beta_E1[a[j],1], x1 = j, y1 = beta_E1[a[j],3], lwd = 3, col = colo[j])
  lines(x = j, y = beta_E1[a[j],2], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

dev.off()

### 16.2.2 Num libros ----

sb11$libros <- as.factor(sb11$libros)

kable(cbind(paste0("$\\beta_", c(4:6),"$"), round(beta_E1[c(4:6),],3)
            ), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

png("libros.png", width = 1200, height = 800, res = 150)

par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

colo <- c("blue","black")[as.numeric((beta_E1[c(4:6),1] < 0) & (0 < beta_E1[c(4:6),3]))+1]

plot(NA, NA, 
     xlab = "Número de libros en casa", 
     ylab = "Efecto tratamiento", 
     main = "Modelo 2", 
     xlim = c(0, 4),ylim=c(-3,22),
     cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:3, 
     labels = levels(sb11$libros)[-1], 
     cex.axis = 0.65)
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:3, col = "gray95", lwd = 1, lty = 2)
a <- c(4:6)
for (j in 1:3) {
  segments(x0 = j, y0 = beta_E1[a[j],1], x1 = j, y1 = beta_E1[a[j],3], lwd = 3, col = colo[j])
  lines(x = j, y = beta_E1[a[j],2], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

dev.off()

### 16.2.3 Estrato ----

kable(cbind(paste0("$\\beta_", c(7:12),"$"), round(beta_E1[c(7:12),],3)), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

sb11$estrato <- as.factor(sb11$estrato)

png("estrato.png", width = 1200, height = 800, res = 150)

par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

colo <- c("blue","black")[as.numeric((beta_E1[c(7:12),1] < 0) & (0 < beta_E1[c(7:12),3]))+1]

plot(NA, NA, 
     xlab = "Estrato", 
     ylab = "Efecto tratamiento", 
     main = "Modelo 2", 
     xlim = c(0, 7),ylim=range(beta_E1[c(7:12),]),
     cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:6, 
     labels = levels(sb11$estrato)[-1], 
     cex.axis = 0.65)
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:6, col = "gray95", lwd = 1, lty = 2)
a <- c(7:12)
for (j in 1:6) {
  segments(x0 = j, y0 = beta_E1[a[j],1], x1 = j, y1 = beta_E1[a[j],3], lwd = 3, col = colo[j])
  lines(x = j, y = beta_E1[a[j],2], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

dev.off()

## 16.3 Covariables municipales ----

beta_M1 <- cbind(apply(M2$THETA[,c(18:20)] ,MARGIN=2, FUN = function(x) quantile(x, 0.025)),
                 colMeans(M2$THETA[,c(18:20)]),
                 apply(M2$THETA[,c(18:20)] ,MARGIN=2, FUN = function(x) quantile(x, 0.975))) 

kable(cbind(paste0("$\\beta_", c(17:19),"$"), round(beta_M1,3)), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

png("municipales.png", width = 1200, height = 800, res = 150)

par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

colo <- c("blue","black")[as.numeric((beta_M1[,1] < 0) & (0 < beta_M1[,3]))+1]

plot(NA, NA, 
     xlab = "Covariables municipales", 
     ylab = "Efecto tratamiento", 
     main = "Modelo 2", 
     xlim = c(0, 4),ylim=range(beta_M1),
     cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:3, 
     labels = c("Docentes por \n estudiante",
                "NBI", 
                "Índice de riesgo"), 
     cex.axis = 0.65)
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:3, col = "gray95", lwd = 1, lty = 2)
for (j in 1:3) {
  segments(x0 = j, y0 = beta_M1[j,1], x1 = j, y1 = beta_M1[j,3], lwd = 3, col = colo[j])
  lines(x = j, y = beta_M1[j,2], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

dev.off()

## 16.4 Covariables departamentales ----

beta_D1 <- cbind(apply(M2$THETA[,c(15:17)] ,MARGIN=2, FUN = function(x) quantile(x, 0.025)),
                 colMeans(M2$THETA[,c(15:17)]),
                 apply(M2$THETA[,c(15:17)] ,MARGIN=2, FUN = function(x) quantile(x, 0.975))) 

kable(cbind(paste0("$\\beta_", c(14:16),"$"), round(beta_D1,3)), 
      col.names = c("","Q1 (2.5\\%)", "Estimación", "Q3 (97.5\\%)"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = F) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

png("departamentales.png", width = 1200, height = 800, res = 150)

par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

colo <- c("blue","black")[as.numeric((beta_D1[,1] < 0) & (0 < beta_D1[,3]))+1]

plot(NA, NA, 
     xlab = "Covariables departamentales", 
     ylab = "Efecto tratamiento", 
     main = "Modelo 2", 
     xlim = c(0, 4), ylim=range(beta_D1),
     cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = 1:3, 
     labels = c("PIB", 
                "Porcentaje de mun \n en riesgo",
                "Porcentaje  de \n población rural"), 
     cex.axis = 0.65)
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:3, col = "gray95", lwd = 1, lty = 2)
for (j in 1:3) {
  segments(x0 = j, y0 = beta_D1[j,1], x1 = j, y1 = beta_D1[j,3], lwd = 3, col = colo[j])
  lines(x = j, y = beta_D1[j,2], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

dev.off()

# 17. Validación de predicción ----

## 17.1 MSE, MAE r R^2 ----

### Modelo 1 ----

g <- prueba %>%
  left_join(mpios, by = c("cole_cod_mcpio_ubicacion", "cole_cod_depto_ubicacion")) %>%
  pull(id)
y <- prueba$punt_global

ZETA <- colMeans(M1$ZETA)  
n <- length(y)
y_hat <- ZETA[g]

MSE_1 <- mean((y-y_hat)^2)

MAE_1 <- mean(abs(y-y_hat)) 

R2_1 <- 1 - n*MSE_1/((n-1)*var(y)) 

### Modelo 2 ----

prueba_2 <- left_join(prueba, Mun, by=c("cole_cod_mcpio_ubicacion"="codmpio"))

prueba_2 <- left_join(prueba_2, Dep, by=c("cole_cod_depto_ubicacion"="Cod"))

prueba_2 <- prueba_2 %>% drop_na(edu_madre, comp, internet, libros, estrato, etnia,
                             PIB, porc_riesgo, porc_rul, doc_alum, nbi, RISK_VICTIM_2022)

prueba_2 <- prueba_2 %>% mutate(PIB = PIB/1e6)

prueba_2 <- prueba_2 %>% mutate(
  edu_madre = as.factor(edu_madre),
  comp = as.factor(comp),
  internet = as.factor(internet),
  libros = as.factor(libros),
  estrato = as.factor(estrato),
  etnia = as.factor(etnia)
)

fit <-  lm(punt_global~edu_madre+comp+internet+libros+estrato+etnia+
             PIB+porc_riesgo+porc_rul+doc_alum+nbi+RISK_VICTIM_2022 ,data = prueba_2)

X <- model.matrix(fit)

mpios_2 <- prueba_2 %>% group_by(cole_cod_mcpio_ubicacion, cole_cod_depto_ubicacion) %>% 
  summarise(yjk_barra = mean(punt_global), #Media por mpios
            s2_jk = var(punt_global), #Var por mpios
            njk = n(), #Nro de estudiantes
            .groups = "drop")

g <- prueba_2 %>%
  left_join(mpios_2, by = c("cole_cod_mcpio_ubicacion", "cole_cod_depto_ubicacion")) %>%
  pull(id)

y <- prueba_2$punt_global
n <- length(y)
THETA <- M2$THETA

y_hat2 <- rowMeans(X %*% t(THETA))

MSE_2 <- mean((y-y_hat2)^2)

MAE_2 <- mean(abs(y-y_hat2)) 

R2_2 <- 1 - n*MSE_2/((n-1)*var(y)) 

### Tabla ----

pred_m1 <- c(MSE_1, MAE_1, R2_1)
pred_m2 <- c(MSE_2, MAE_2, R2_2)

pred <- rbind(pred_m1,pred_m2)

rownames(pred) <- c("M1", "M2")

kable(pred, 
      col.names = c("MSE","MAE", "$R^2$"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 3, row.names = T) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")



# Anexo ----

cvmc1 <- rbind( summary((apply(M1$ZETA, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$ZETA)))*100/colMeans(M1$ZETA)),
                summary((apply(M1$KAPPA_JK, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$KAPPA_JK)))*100/colMeans(M1$KAPPA_JK)),
                summary((apply(M1$THETA, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$THETA)))*100/colMeans(M1$THETA)),
                summary((apply(M1$SIGMA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$SIGMA_K)))*100/colMeans(M1$SIGMA_K)),
                summary((apply(M1$KAPPA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M1$KAPPA_K)))*100/colMeans(M1$KAPPA_K)),
                summary((sd(M1$MU$V1)/sqrt(coda::effectiveSize(M1$MU)))*100/mean(M1$MU$V1)),
                summary((sd(M1$TAU$V1)/sqrt(coda::effectiveSize(M1$TAU)))*100/mean(M1$TAU$V1)),
                summary((sd(M1$SIGMA$V1)/sqrt(coda::effectiveSize(M1$SIGMA)))*100/mean(M1$SIGMA$V1)),
                summary((sd(M1$ALPHA$V1)/sqrt(coda::effectiveSize(M1$ALPHA)))*100/mean(M1$ALPHA$V1)),
                summary((sd(M1$BETA$V1)/sqrt(coda::effectiveSize(M1$BETA))))*100/mean(M1$BETA$V1))

cvmc2 <- rbind(
  summary((sd(M2$THETA$V1)/sqrt(coda::effectiveSize(M2$THETA)))*100/abs(mean(M2$THETA$V1))),
  summary((apply(M2$THETA[, 2:14], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[, 2:14])))*100/abs(colMeans(M2$THETA[, 2:14]))),
  summary((apply(M2$THETA[, 15:17], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[, 15:17])))*100/abs(colMeans(M2$THETA[, 15:17]))),
  summary((apply(M2$THETA[, 18:20], MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$THETA[, 18:20])))*100/abs(colMeans(M2$THETA[, 18:20]))),
  summary((apply(M2$KAPPA_JK, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$KAPPA_JK)))*100/colMeans(M2$KAPPA_JK)),
  summary((apply(M2$KAPPA_K, MARGIN = 2, FUN = sd)/sqrt(coda::effectiveSize(M2$KAPPA_K)))*100/colMeans(M2$KAPPA_K)),
  summary((sd(M2$SIG_B$V1)/sqrt(coda::effectiveSize(M2$SIG_B)))*100/mean(M2$SIG_B$V1)),
  summary((sd(M2$SIG_E$V1)/sqrt(coda::effectiveSize(M2$SIG_E)))*100/mean(M2$SIG_E$V1)),
  summary((sd(M2$SIG_M$V1)/sqrt(coda::effectiveSize(M2$SIG_M)))*100/mean(M2$SIG_M$V1)),
  summary((sd(M2$SIG_D$V1)/sqrt(coda::effectiveSize(M2$SIG_D)))*100/mean(M2$SIG_D$V1)),
  summary((sd(M2$ALPHA$V1)/sqrt(coda::effectiveSize(M2$ALPHA)))*100/mean(M2$ALPHA$V1)),
  summary((sd(M2$BETA_K$V1)/sqrt(coda::effectiveSize(M2$BETA_K)))*100/mean(M2$BETA_K$V1))
)


kable(cbind(p1, round(cvmc1,3)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                                                 "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")

kable(cbind(p2, round(cvmc2,3)), col.names = c("Paramétro","Mínimo", "Cuartíl 1", "Mediana",
                                               "Media", "Cuartíl 3", "Máximo"), 
      format = "latex", 
      booktabs = TRUE, 
      escape = FALSE, digits = 0) %>%
  kable_styling(position = "center", latex_options = "HOLD_position")
