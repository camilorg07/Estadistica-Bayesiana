#Autores:


#---------- Juan Andrés Camacho Zaráte (jucamachoz@unal.edu.co)
#---------- Camilo Alejandro Raba Gómez (craba@unal.edu.co)


#-----------------------Solución puntos caso 2----------------------------------

#0. Librerias y carga de datos -------------------------------------------------

Sys.setlocale("LC_ALL", "es_ES.UTF-8") #poner R en español para evitar probleas con caracteres

#0.1 Librerias -----------------------------------------------------------------

library(gridExtra)
library(ggpubr)
library(ggplot2)
library(sf)
library(latex2exp)
library(igraph)
library(magick)
library(readr)
library(stringr)
library(readxl)
library(dplyr)
library(knitr)
library(coda)
library(sqldf)
library(data.table)
library(progress)

#0.2 Modelos--------------------------------------------------------------------

load("M1.RData")
load("M2.RData")
load("M3.RData")
load("M4.RData")
M2<-M21
rm(M21)

#0.3 Datos ---------------------------------------------------------------------

pobreza_m=read_excel("datos/pobreza monetaria.xls",sheet = 2, range = "A16:P41") 
colnames(pobreza_m)[1]="Departamento"


est_edu=read_delim("datos/estadísticas educación.csv", delim = ",")

sb2022_2=read_delim("datos/Saber 11 2022-2.TXT", delim = ";")

#1. Nacionalidad colombiana
#2. Residencia en Colombia
#3. Estado en publicar
#4. No tomar colegios de San Andrés
#5. Eliminar NAS

datos=sb2022_2%>%filter(ESTU_NACIONALIDAD=="COLOMBIA" & ESTU_PAIS_RESIDE=="COLOMBIA" & ESTU_ESTADOINVESTIGACION=="PUBLICAR" & COLE_COD_DEPTO_UBICACION!=88 & !is.na(COLE_COD_DEPTO_UBICACION) & !is.na(COLE_COD_MCPIO_UBICACION) & !is.na(PUNT_GLOBAL))

y=datos%>%select(COLE_COD_DEPTO_UBICACION,
                 COLE_COD_MCPIO_UBICACION,
                 PUNT_GLOBAL)%>%
  arrange(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION)

shp <- sf::st_read("MGN_DPTO_POLITICO.shp", quiet = T)
colnames(shp)[1]="DPTO"
shp = shp[shp$DPTO != '88',]

mundoshp <- st_read("admin00.shp",quiet=TRUE)
mundocol <- mundoshp %>% 
  filter(CNTRY_NAME %in% c("Peru","Brazil","Venezuela","Ecuador","Panama"))
box=st_bbox(shp)

#-------------------------------------------------------------------------------

# Crear identificador C:nico para departamento y municipio
datos1 <- datos %>%
  mutate(COLE_COD_DEPTO_MCPIO = paste(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION, sep = "-"))

# Agrupar correctamente por departamento y municipio
X_deptos <- datos1 %>%
  group_by(COLE_COD_DEPTO_UBICACION) %>%
  summarise(nj_depto = n_distinct(COLE_COD_MCPIO_UBICACION)) %>%
  mutate(dpto = 1:n())

X=as.data.frame(X_deptos)

y <- datos1 %>%
  left_join(X_deptos, by = "COLE_COD_DEPTO_UBICACION") %>%
  select(COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION, dpto, PUNT_GLOBAL) %>%
  arrange(dpto, COLE_COD_MCPIO_UBICACION) %>%
  mutate(mun_id = group_indices(., COLE_COD_DEPTO_UBICACION, COLE_COD_MCPIO_UBICACION))

# NC:mero de municipios (K) y estudiantes por municipio (njk)
muni <- y %>%
  group_by(dpto, mun_id) %>%
  summarise(njk = n(), .groups = 'drop')

K <- nrow(muni)
nj <- muni %>% group_by(dpto) %>% summarise(nj = n()) %>% pull(nj)
njk <- muni$njk
J <- length(unique(y$dpto))  # NC:mero de departamentos

# Se define un 'diccionario' con el código de departamentos DANE para posterior presentación:
Departamentos = cbind('CodDepartamento' = c('05', '08', '11', '13', '15', '17', '18', '19', '20',
                                            '23', '25', '27', '41', '44', '47', '50', '52', '54', 
                                            '63', '66', '68', '70', '73', '76', '81', '85','86',
                                            '91', '94', '95', '97', '99'),
                      'departamento' = c("ANTIOQUIA", "ATLANTICO", "BOGOTA", "BOLIVAR",
                                         "BOYACA", "CALDAS", "CAQUETA", "CAUCA", "CESAR",
                                         "CORDOBA", "CUNDINAMARCA", "CHOCO", "HUILA", 
                                         "LA GUAJIRA", "MAGDALENA", "META", "NARIÑO", 
                                         "NORTE SANTANDER", "QUINDIO", "RISARALDA", "SANTANDER", 
                                         "SUCRE", "TOLIMA", "VALLE", "ARAUCA", 
                                         "CASANARE", "PUTUMAYO", "AMAZONAS",
                                         "GUAINIA", "GUAVIARE", "VAUPES", "VICHADA"))
Departamentos = as.data.frame(Departamentos)

# LEER EL ARCHIVO SEPARADO POR ;

data = datos1

# QUITANDO LOS DATOS DESCRITOS EN EL CASO                                       # Únicamente los estudiantes:
data = subset(data, ESTU_NACIONALIDAD == 'COLOMBIA')                            # - Con nacionalidad colombiana
data = subset(data,ESTU_PAIS_RESIDE == 'COLOMBIA')                              # - Que reside en Colombia
data = subset(data, ESTU_ESTADOINVESTIGACION == 'PUBLICAR')                     # - Resultados para publicación
data = subset(data, COLE_DEPTO_UBICACION != 'SAN ANDRES')                       # - Colegios que no están en San Andrés
# colnames((head(FilteredData[, grepl("PUNT", names(FilteredData))])))

# Seleccionando las columnas de intéres de la base de datos1:
data = data[, c('COLE_DEPTO_UBICACION', 'COLE_COD_DEPTO_UBICACION', 
                'COLE_MCPIO_UBICACION','COLE_COD_MCPIO_UBICACION',
                'PUNT_GLOBAL')]
colnames(data) = c('departamento','CodDepartamento','municipio','CodMunicipio','puntaje')

# Quitando los datos1 sin datos1 en algunas de las columnas seleccionadas
data = data[complete.cases(data),]                                              # - Sin datos1 faltantes en ninguna de las columnas seleccionadas.

# Corrección código 5 a 05 y 8 a 08
data$CodDepartamento = as.character(data$CodDepartamento)
data$CodDepartamento[data$CodDepartamento == '5'] = '05'
data$CodDepartamento[data$CodDepartamento == '8'] = '08'

# CALCULANDO LAS ESTADÍSTICAS SUFICIENTES POR DEPARTAMENTOS Y POR MUNICIPIOS:
# install.packages('dplyr')
dataDepartamentos = data %>%
  group_by(CodDepartamento) %>%
  summarise(CodDepartamento = unique(CodDepartamento),
            departamento = unique(departamento),
            n = n(),
            mean = mean(puntaje),
            var = var(puntaje),
            nMunicipios = n_distinct(CodMunicipio),)

dataMunicipios = data %>%
  group_by(CodMunicipio) %>%
  summarise(CodMunicipio = unique(CodMunicipio),
            municipio = unique(municipio),
            n = n(),
            mean = mean(puntaje),
            var = var(puntaje),
            min = min(puntaje),
            max = max(puntaje),
            IQR = IQR(puntaje),
            mediana = median(puntaje))

v=3

#1. Punto 1 --------------------------------------------------------------------


#1.1 Mapa por media del puntaje global (Departamentos) -------------------------

sample_m <- y%>%group_by(COLE_COD_DEPTO_UBICACION)%>%
  summarise(xbar=mean(PUNT_GLOBAL))
colnames(sample_m)[1] <- "DPTO"

m1<-inner_join(x = shp, y = sample_m, by = c("DPTO"="DPTO")) %>% 
  select(DPTO, xbar, geometry) %>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = xbar), size = 0.125, color = "gray30", linetype="solid") +
  #scale_fill_manual(values = c("#F38D30", "#C590B3", "#FFA07A", "#FFFACD", "#95B958")) +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    high = "#87B560",   # Color claro para valores bajos
    low = "#D44D44",   # Color oscuro para valores altos
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  #geom_sf_text(aes(label=Nombre),col="black",
  #            fontface="bold",size=3,fun.geometry=function(x) sf::st_centroid(x))+
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Media del \n Puntaje \n Global")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

#1.2 Mapa por indice de pobreza monetaria --------------------------------------

departamentos <- data.frame(
  Codigo = c(5, 8, 11, 13, 15, 17, 18, 19, 20, 23, 25, 27, 41, 44, 47, 50, 52, 54, 63, 66, 68, 70, 73, 76, 81, 85, 86, 88, 91, 94, 95, 97, 99),
  Nombre = c("Antioquia", "Atlántico", "Bogotá D.C.", "Bolívar", "Boyacá", "Caldas", "Caquetá", "Cauca", "Cesar", "Córdoba", 
             "Cundinamarca", "Chocó", "Huila", "La Guajira", "Magdalena", "Meta", "Nariño", "Norte de Santander", 
             "Quindío", "Risaralda", "Santander", "Sucre", "Tolima", "Valle del Cauca", "Arauca", "Casanare", "Putumayo", 
             "San Andrés", "Amazonas", "Guainía", "Guaviare", "Vaupés", "Vichada")
)

#Convertir a factor 

departamentos$Codigo=factor(departamentos$Codigo)


pobreza_m <- pobreza_m%>%select(Departamento, '2018')

pobreza_m <- inner_join(departamentos, pobreza_m, by = c("Nombre"="Departamento"))%>%
  select(Codigo, pm = '2018') 

pobreza_m$Codigo <- as.character(pobreza_m$Codigo)

pobreza_m$Codigo[pobreza_m$Codigo=="5"] = "05"
pobreza_m$Codigo[pobreza_m$Codigo=="8"] = "08"

m2<-left_join(x = shp, y = pobreza_m, by = c("DPTO"="Codigo"))%>%
  select(DPTO ,pm, geometry)%>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = pm), size = 0.125, color = "gray30", linetype="solid") +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    low = "#87B560",   # Color claro para valores bajos
    high = "#D44D44",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  labs(x="Longitud",y="Latitud",title="Colombia",fill="IPM \n en 2018")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

#1.3 Panel de los mapas --------------------------------------------------------

m3<-ggarrange(m1, m2, ncol = 2)
m3
#ggsave("mapas_punto1.png", plot = m3, width = 7, height = 4, dpi = 300)

#2. Punto 2 --------------------------------------------------------------------

#2.1 Mapa por media del puntaje global (municipios) ----------------------------

sample_m <- y%>%group_by(COLE_COD_MCPIO_UBICACION)%>%
  summarise(xbar=mean(PUNT_GLOBAL))
colnames(sample_m)[1] <- "MUNI"

shp1 <- sf::st_read("MGN_MPIO_POLITICO.shp", quiet = T)
shp1 = shp1[shp1$DPTO_CCDGO != '88',]
colnames(shp1)[6] <- "MUNI" 


m1<-left_join(x = shp1, y = sample_m, by = c("MUNI"="MUNI")) %>% 
  select(MUNI, xbar, geometry) %>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = xbar), size = 0.0125, color = "gray30", linetype="solid") +
  #scale_fill_manual(values = c("#F38D30", "#C590B3", "#FFA07A", "#FFFACD", "#95B958")) +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    high = "#AB3F24",   # Color claro para valores bajos
    low = "#F4D166",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  #geom_sf_text(aes(label=Nombre),col="black",
  #            fontface="bold",size=3,fun.geometry=function(x) sf::st_centroid(x))+
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Media del \n Puntaje \n Global")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

#2.2 Mapa por cobertura neta secundaria ----------------------------------------

est_edu <- est_edu%>%filter(AÑO == "2022")
est_edu <- est_edu[est_edu$CÓDIGO_DEPARTAMENTO!= '88',]
est_edu <- est_edu%>%select(MUNI = CÓDIGO_MUNICIPIO,
                            CNS = COBERTURA_NETA_SECUNDARIA)

m2<-left_join(x = shp1, y = est_edu, by = c("MUNI"="MUNI"))%>%
  select(MUNI ,CNS, geometry)%>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    high = "#AB3F24",   # Color claro para valores bajos
    low = "#F4D166",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  labs(x="Longitud",y="Latitud",title="Colombia",fill="CNS \n en 2022")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

#2.3 Panel de los dos mapas ----------------------------------------------------

m3<-ggarrange(m1, m2, ncol = 2)
m3
#ggsave("mapas_punto2.png", plot = m3, width = 7, height = 4, dpi = 300)

#3. Punto 5 --------------------------------------------------------------------

#3.1 DAG del primer modelo -----------------------------------------------------

#png("DAG1.png", width = 1000, height = 800, res = 150)

Modelo1 <- graph(edges = c(2, 1, 
                           3, 4,
                           4, 1,
                           5, 2,
                           6, 2,
                           7, 3,
                           8, 3,
                           13, 4,
                           9, 5,
                           10, 5,
                           11, 6,
                           12, 6
), directed = TRUE)

V(Modelo1)$label = c(TeX("$y_{i,j}$"), 
                     TeX("$\\theta_{j}$"), TeX("$\\sigma^2$"),TeX("$\\xi_{i,j}$"),
                     TeX('$\\mu$'), TeX('$\\tau^2$') ,
                     TeX('$\\nu_0$'), TeX('$\\sigma^2_0$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'),
                     TeX('$\\eta_0$$'), TeX('$\\tau^2_0$'),
                     TeX('$v$'))

V(Modelo1)$label.color = "black"; V(Modelo1)$color = "white"

V(Modelo1)[1:6]$shape = 'circle'; V(Modelo1)[7:13]$shape = "square"

V(Modelo1)$size = 25

# Atributos de las flechas del DAG
E(Modelo1)$arrow.size = 0.5;E(Modelo1)$color = 'black'

distributions <- c(
  TeX("$N$"),  # y_ij
  TeX("$N$"),            # theta_j
  TeX("$GI$"),      # sigma^2
  TeX("$G$"),             # xi_ij
  TeX("$N$"),      # mu
  TeX("$GI$"),      # tau^2
  "", "", "", "", "", "", ""  # Hiperparámetros (sin distribución mostrada)
)


custom_layout <- matrix(c(
  0, -2,  # y_ij (abajo)
  -1.7, -1,  # theta_j
  1.7, -1,  # sigma^2
  0, -1,  # xi_ij
  -2.52, 0,  # mu
  -0.65, 0,  # tau^2
  1.75, 1,  # nu_0
  2.7, 1,  # sigma^2_0
  -3, 1,  # mu_0
  -2.05, 1,  # gamma^2_0
  -1.1, 1,  # eta_0
  -.2, 1,  # tau^2_0
  0.8, 1   # v
), ncol = 2, byrow = TRUE)


plot(Modelo1, 
     layout = custom_layout,
     vertex.label = V(Modelo1)$label, 
     vertex.shape = V(Modelo1)$shape, 
     vertex.color = V(Modelo1)$color, 
     vertex.size = V(Modelo1)$size,
     edge.arrow.size = E(Modelo1)$arrow.size,
     edge.color = E(Modelo1)$color,
     main = "Modelo 1: DAG")
text(c(0.15, -.45,  0.15,  .75, -.75, -.12),  # Ajuste en X
     c(-1.125, -.4, -.4, -.4,  .2,  .2),  # Ajuste en Y
     labels = distributions[1:6], 
     col = "black", cex = 0.7, pos = 4)

#dev.off()


#3.2 DAG del segundo modelo-----------------------------------------------------

#png("DAG2.png", width = 1000, height = 800, res = 150)

Modelo2 =  graph(edges = c(2, 1,                     # Primer nivel 
                           3, 1, 
                           4, 3,                        # Segundo nivel
                           5, 2,
                           6, 2,
                           7, 4,
                           8, 4,
                           9, 5,
                           10, 5,
                           11, 6,
                           12, 6,
                           13, 3,
                           14, 7,
                           15, 8,
                           16, 8
), directed = TRUE)

# Atributos de los nodos:
V(Modelo2)$label = c(TeX("$y_{i,j}$"), 
                     TeX("$\\theta_{j}$"),TeX("$\\xi_{i,j}$") ,TeX("$\\sigma^2_j$"),
                     TeX('$\\mu$'), TeX('$\\tau^2$'),
                     TeX('$\\nu$'), TeX('$\\sigma^2$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'),
                     TeX('$\\eta_0$$'), TeX('$\\tau^2_0$'), 
                     TeX("$\\v$"),
                     TeX("$\\lambda_0$"),
                     TeX('$\\alpha_0$'), TeX('$\\beta_0$') )

V(Modelo2)$label.color = "black"
V(Modelo2)$color = "white"

# Formas de los nodos: Variables aleatorias → círculo, Hiperparámetros → cuadrado
V(Modelo2)[1:8]$shape = 'circle'  # Variables aleatorias
V(Modelo2)[9:16]$shape = 'square'  # Hiperparámetros

# Tamaño de nodos
V(Modelo2)$size = 25

# Configuración de flechas
E(Modelo2)$arrow.size = 0.5
E(Modelo2)$color = 'black'

# Definir posiciones manualmente
custom_layout <- matrix(c(
  0, -2,   # y_ij (abajo)
  -2, -1,  # θ_j (medio)
  0, -1,  # ξ_ij (medio)
  2, -1,  # σ_j^2 (medio)
  -3.5, 0,   # μ (encima de θ_j)
  -.75, 0,   # τ^2 (encima de θ_j)
  1.5, 0,   # ν (encima de σ_j^2)
  3, 0,   # σ^2 (encima de σ_j^2)
  -4, 1,   # μ_0 (hiperparámetro arriba)
  -3, 1,   # γ^2_0 (hiperparámetro arriba)
  -.7, 1,   # η_0 (hiperparámetro arriba)
  -1.8, 1,   # τ^2_0 (hiperparámetro arriba)
  .5, 1,   # v (hiperparámetro arriba)
  1.8, 1,   # λ_0 (hiperparámetro arriba)
  3, 1,   # α_0 (hiperparámetro arriba)
  4, 1    # β_0 (hiperparámetro arriba)
), ncol = 2, byrow = TRUE)
distributions <- c(
  TeX("$N$"),  # y_ij
  TeX("$N$"),  # θ_j
  TeX("$GI$"), # ξ_ij
  TeX("$G$"),  # σ_j^2
  TeX("$N$"),  # μ
  TeX("$GI$"), # τ^2
  TeX("Exp"),  # ν
  TeX("$G$")   # σ^2
)
# Dibujar el DAG con posiciones personalizadas
plot(Modelo2, 
     layout = custom_layout,
     vertex.label = V(Modelo2)$label, 
     vertex.shape = V(Modelo2)$shape, 
     vertex.color = V(Modelo2)$color, 
     vertex.size = V(Modelo2)$size,
     edge.arrow.size = E(Modelo2)$arrow.size,
     edge.color = E(Modelo2)$color,
     main = "Modelo 2: DAG") 
text(c(0.1, -.4,  0.1,  .6, -.8, -0.15, 0.4, 0.8),  # Ajuste en X
     c(-1, -.4, -.4, -.4,  .2,  .2, .2 ,.2),  # Ajuste en Y
     labels = distributions, 
     col = "black", cex = 0.7, pos = 4)

#dev.off()

#3.3 DAG del tercer modelo -----------------------------------------------------

#png("DAG3.png", width = 1000, height = 800, res = 150)

Modelo3 = graph(edges = c(2,1,
                          3,1,
                          4,3,
                          5,2,
                          6,2,
                          7,5,
                          8,5,
                          9, 7,
                          10, 7,
                          11, 8,
                          12, 8,
                          13, 6,
                          14, 6,
                          15, 3,
                          16, 4,
                          17, 4), directed = TRUE)

V(Modelo3)$label = c(TeX("$y_{_{i, j, k}}$"), 
                     TeX("$\\zeta_{j, k}$"), TeX("$\\xi_{_{i, j, k}}$") ,TeX("$\\kappa^2$"),
                     TeX('$\\theta_k$'), TeX('$\\sigma^2$'),
                     TeX('$\\mu$$'), TeX('$\\tau^2$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'),
                     TeX('$\\eta_0$'), TeX('$\\tau^2_0$'),
                     TeX('$\\nu_0$'), TeX('$\\sigma^2_0$'),
                     TeX('$\\v$'),
                     TeX('$\\epsilon_0$'), TeX('$\\kappa^2_0$')
)

distributions <- c(
  TeX("$N$"),  # y_{i,j,k}
  TeX("$N$"),  # ζ_{j,k}
  TeX("$GI$"), # ξ_{i,j,k}
  TeX("$G$"),  # κ^2
  TeX("$N$"),  # θ_k
  TeX("$GI$"),  # σ^2
  TeX("$N$"),  # μ
  TeX("$GI$")  # τ^2
)

V(Modelo3)$label.color = "black"; V(Modelo3)$color = "white"

V(Modelo3)[1:8]$shape = 'circle'; V(Modelo3)[9:17]$shape = "square"

V(Modelo3)$size = 25


# Atributos de las flechas del DAG
E(Modelo3)$arrow.size = 0.5;E(Modelo3)$color = 'black'

custom_layout <- matrix(c(
  0, -2,   # y_{i,j,k} (abajo)
  -2, -1,  # ζ_{j,k} (medio)
  0, -1,  # ξ_{i,j,k} (medio)
  2, -1,  # κ^2 (medio)
  -3, 0,   # θ_k (encima de ζ_{j,k})
  -.75, 0,   # σ^2 (encima de ζ_{j,k})
  -4, 1,   # μ (encima de θ_k)
  -2, 1,   # τ^2 (encima de θ_k)
  -4.5, 2,   # μ_0 (hiperparámetro arriba)
  -3.5, 2,   # γ^2_0 (hiperparámetro arriba)
  -1.5, 2,   # η_0 (hiperparámetro arriba)
  -2.5, 2,   # Τ^2_0 (hiperparámetro arriba)
  -.5, 2,   # ν_0 (hiperparámetro arriba)
  .5, 2,   # σ^2_0 (hiperparámetro arriba)
  1.5, 2,   # v (hiperparámetro arriba)
  2.5, 2,   # ε_0 (hiperparámetro arriba)
  3.5, 2    # κ^2_0 (hiperparámetro arriba)
), ncol = 2, byrow = TRUE)

plot(Modelo3, 
     layout = custom_layout,
     vertex.label = V(Modelo3)$label, 
     vertex.shape = V(Modelo3)$shape, 
     vertex.color = V(Modelo3)$color, 
     vertex.size = V(Modelo3)$size,
     edge.arrow.size = E(Modelo3)$arrow.size,
     edge.color = E(Modelo3)$color,
     main = "Modelo 3: DAG") 
text(c(0.2, -.28,  0.2,  .68, -.53, 0.03, -.78, -0.28),
     c(-1.1, -.6, -.6, -.6,  -.1,  -.1, .4, .4),
     labels = distributions, 
     col = "black", cex = 0.8, pos = 4) 
#dev.off()

#3.4 DAG del cuarto modelo -----------------------------------------------------

#png("DAG4.png", width = 1000, height = 800, res = 150)

Modelo4 <- graph(edges = c(2,1,
                           3,1,
                           4,3,
                           5,2,
                           6,2,
                           7,5,
                           8,5,
                           9, 6,
                           10, 6,
                           11, 7,
                           12, 7,
                           13, 8,
                           14, 8,
                           15, 9,
                           16, 10,
                           17, 10,
                           18, 3,
                           19, 4,
                           20, 4), directed = TRUE)

V(Modelo4)$label = c(TeX("$y_{_{i, j, k}}$"), 
                     TeX("$\\zeta_{j, k}$"), TeX("$\\xi_{_{i, j, k}}$") ,TeX("$\\kappa^2$"),
                     TeX('$\\theta_k$'), TeX('$\\sigma^2_k$'),
                     TeX('$\\mu$$'), TeX('$\\tau^2$'),
                     TeX('$\\nu$'), TeX('$\\sigma^2$'),
                     TeX('$\\mu_0$'), TeX('$\\gamma^2_0$'),
                     TeX('$\\eta_0$'), TeX('$\\tau^2_0$'),
                     TeX('$\\lambda_0$'),
                     TeX('$\\alpha_0$'),TeX('$\\beta_0$'),
                     TeX('$\\v$'),
                     TeX('$\\epsilon_0$'), TeX('$\\kappa^2_0$')
)

V(Modelo4)$label.color = "black"; V(Modelo4)$color = "white"

V(Modelo4)[1:10]$shape = 'circle'; V(Modelo4)[11:20]$shape = "square"

V(Modelo4)$size = 25
# Atributos de las flechas del DAG
E(Modelo4)$arrow.size = 0.5;E(Modelo4)$color = 'black'

custom_layout <- matrix(c(
  0, -2,   # y_{i,j,k} (abajo)
  -3, -1,  # ζ_{j,k} (medio)
  0, -1,   # ξ_{i,j,k} (medio)
  3, -1,   # κ^2 (medio)
  -4.5, 0,   # θ_k (encima de ζ_{j,k})
  -1.5, 0,    # σ_k^2 (encima de ζ_{j,k})
  -6, 1,   # μ (encima de θ_k)
  -3.6, 1,   # τ^2 (encima de θ_k)
  -1.5, 1,   # ν (encima de σ_k^2)
  1, 1,  # σ^2 (encima de σ_k^2)
  -6.75, 2, # μ_0 (hiperparámetro arriba)
  -5.25, 2, # γ^2_0 (hiperparámetro arriba)
  -3.75, 2, # η_0 (hiperparámetro arriba)
  -2.25, 2, # τ^2_0 (hiperparámetro arriba)
  -0.75, 2, # λ_0 (hiperparámetro arriba)
  0.75, 2,  # α_0 (hiperparámetro arriba)
  2.25, 2,  # β_0 (hiperparámetro arriba)
  3.75, 2,  # v (hiperparámetro arriba)
  5.25, 2,  # ε_0 (hiperparámetro arriba)
  6.75, 2   # κ^2_0 (hiperparámetro arriba)
), ncol = 2, byrow = TRUE)

distributions <- c(
  TeX("$N$"),  # y_{i,j,k}
  TeX("$N$"),  # ζ_{j,k}
  TeX("$GI$"), # ξ_{i,j,k}
  TeX("$G$"),  # κ^2
  TeX("$N$"),  # θ_k
  TeX("$GI$"),  # σ_k^2
  TeX("$N$"),  # μ
  TeX("$GI$"),  # τ^2
  TeX("$Exp$"), # nu
  TeX("$G$") # σ^2
)


plot(Modelo4, 
     layout = custom_layout,
     vertex.label = V(Modelo4)$label, 
     vertex.shape = V(Modelo4)$shape, 
     vertex.color = V(Modelo4)$color, 
     vertex.size = V(Modelo4)$size,
     edge.arrow.size = E(Modelo4)$arrow.size,
     edge.color = E(Modelo4)$color,
     main = "Modelo 4: DAG") 
text(c(0.05, -.37,  0.05,  .5, -.605, -0.175, -.83, -0.48, -.15, .2),
     c(-1.1, -.6, -.6, -.6,  -.1,  -.1, .4, .4, .4,.4),
     labels = distributions, 
     col = "black", cex = 0.8, pos = 4) 

#dev.off()

#3.4 Gráfico de los cuatro DAG -------------------------------------------------

#img1 <- image_read("DAG1.png")
#img2 <- image_read("DAG2.png")
#img3 <- image_read("DAG3.png")
#img4 <- image_read("DAG4.png")

# Crear un collage 2x2
final_plot <- image_append(c(
  image_append(c(img1, img2), stack = FALSE),
  image_append(c(img3, img4), stack = FALSE)
), stack = TRUE)

final_plot

#image_write(final_plot, path = "DAGS.png", format = "png")

#4. Punto 6 --------------------------------------------------------------------

#png("log_verosimilitud.png", width = 1000, height = 800, res = 150) 
par(mfrow = c(2,2))

#4.1 Log-verosimilitud M1-------------------------------------------------------

plot(x = 1:10000, y = M1$LL$V1, type = "p", pch = ".", cex.axis = 0.8, main = "Modelo 1", xlab = "Iteración", ylab = "Log-verosimilitud", col="#458B00")

abline(h=mean(M1$LL$V1), col="#008B00")

#4.2 Log-verosimilitud M2-------------------------------------------------------

plot(x = 1:10000, y = M2$LL$V1, type = "p", pch = ".", cex.axis = 0.8, main = "Modelo 2", xlab = "Iteración", ylab = "Log-verosimilitud", col="#E66820")

abline(h=mean(M2$LL$V1), col="#A94322")

#4.3 Log-verosimilitud M3-------------------------------------------------------

plot(x = 1:10000, y = M3$LL$V1, type = "p", pch = ".", cex.axis = 0.8, main = "Modelo 3", xlab = "Iteración", ylab = "Log-verosimilitud", col="#5A8BB7")

abline(h=mean(M3$LL$V1), col="#2E5B88")

#4.4 Log-verosimilitud M4-------------------------------------------------------

plot(x = 1:10000, y = M4$LL$V1, type = "p", pch = ".", cex.axis = 0.8, main = "Modelo 4", xlab = "Iteración", ylab = "Log-verosimilitud", col="#E84746")

abline(h=mean(M4$LL$V1), col="#B5153A")

#dev.off()

#5. Punto 7 --------------------------------------------------------------------

# 5.1 Primer modelo ------------------------------------------------------------

THETA <- M1$THETA
MU <- M1$MU
SIG2 <- M1$SIG2
TAU2 <- M1$TAU2

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)


knitr::kable(t((c("$$\\theta_j$$",round(c(min(eff_theta), quantile(eff_theta, .25),
                                          quantile(eff_theta,.5), mean(eff_theta),
                                          quantile(eff_theta, .75,), max(eff_theta)),0)))), 
             col.names = c("Estadístico","Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"), digits = 0)


knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$"), round(c(eff_sig2, eff_mu, eff_tau2),0)), row.names = FALSE)

# 5.2 Segundo modelo -----------------------------------------------------------

THETA <- M2$THETA
MU <- M2$MU
SIG2 <- M2$SIG2
TAU2 <- M2$TAU2
SIG2J <- M2$SIG2J
NU <- M2$NU

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_sig2j <- coda::effectiveSize(SIG2J)
eff_nu <- coda::effectiveSize(NU)

knitr::kable(t(round(c(min(eff_theta), quantile(eff_theta, .25),
                       quantile(eff_theta,.5), mean(eff_theta),
                       quantile(eff_theta, .75,), max(eff_theta)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(eff_sig2j), quantile(eff_sig2j, .25),
                       quantile(eff_sig2j,.5), mean(eff_sig2j),
                       quantile(eff_sig2j, .75,), max(eff_sig2j)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\sigma_j^2$$", digits = 0)

knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\nu$$"), round(c(eff_sig2, eff_mu, eff_tau2, eff_nu),0)), row.names = FALSE)

# 5.3 Tercer modelo ------------------------------------------------------------

THETA <- M3$THETA
MU <- M3$MU
SIG2 <- M3$SIG2
TAU2 <- M3$TAU2
EPSJK <- M3$EPSJK
K2 <- M3$K2

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_epsjk <- coda::effectiveSize(EPSJK)
eff_k2 <- coda::effectiveSize(K2)

knitr::kable(t(round(c(min(eff_theta), quantile(eff_theta, .25),
                       quantile(eff_theta,.5), mean(eff_theta),
                       quantile(eff_theta, .75,), max(eff_theta)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(eff_epsjk), quantile(eff_epsjk, .25),
                       quantile(eff_epsjk,.5), mean(eff_epsjk),
                       quantile(eff_epsjk, .75,), max(eff_epsjk)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\xi_{j,k}$$")

knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\kappa^2$$"), round(c(eff_sig2, eff_mu, eff_tau2, eff_k2),0)), row.names = FALSE)

# 5.4 Cuarto modelo ------------------------------------------------------------

THETA <- M4$THETA
MU <- M4$MU
SIG2 <- M4$SIG2
TAU2 <- M4$TAU2
EPSJK <- M4$EPSJK
K2 <- M4$K2
SIG2K <- M4$SIG2K
NU <- M4$NU

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_epsjk <- coda::effectiveSize(EPSJK)
eff_k2 <- coda::effectiveSize(K2)
eff_sig2k <- coda::effectiveSize(SIG2K)
eff_nu <- coda::effectiveSize(NU)

knitr::kable(t(round(c(min(eff_theta), quantile(eff_theta, .25),
                       quantile(eff_theta,.5), mean(eff_theta),
                       quantile(eff_theta, .75,), max(eff_theta)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(eff_epsjk), quantile(eff_epsjk, .25),
                       quantile(eff_epsjk,.5), mean(eff_epsjk),
                       quantile(eff_epsjk, .75,), max(eff_epsjk)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\xi_{j,k}$$")

knitr::kable(t(round(c(min(eff_sig2k), quantile(eff_sig2k, .25),
                       quantile(eff_sig2k,.5), mean(eff_sig2k),
                       quantile(eff_sig2k, .75,), max(eff_sig2k)),0)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\sigma_k^2$$")

knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\kappa^2$$", "$$\\nu$$"), round(c(eff_sig2, eff_mu, eff_tau2, eff_k2, eff_nu),0)), row.names = FALSE)

# 6. Punto 8 -------------------------------------------------------------------

# 6.1 Primer modelo ------------------------------------------------------------

THETA <- M1$THETA
MU <- M1$MU
SIG2 <- M1$SIG2
TAU2 <- M1$TAU2

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)

EMC1 <- apply(X = THETA, MARGIN = 2, FUN = sd)/sqrt(eff_theta)

EMC2 <- apply(X = MU, MARGIN = 2, FUN = sd)/sqrt(eff_mu)

EMC3 <- apply(X = SIG2, MARGIN = 2, FUN = sd)/sqrt(eff_sig2)

EMC4 <- apply(X = TAU2, MARGIN = 2, FUN = sd)/sqrt(eff_tau2)

knitr::kable(t(round(c(min(EMC1), quantile(EMC1, .25),
                       quantile(EMC1,.5), mean(EMC1),
                       quantile(EMC1, .75,), max(EMC1)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")


knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$"), round(c(EMC2, EMC3, EMC4),3)), row.names = FALSE)

# 6.2 Segundo modelo -----------------------------------------------------------

THETA <- M2$THETA
MU <- M2$MU
SIG2 <- M2$SIG2
TAU2 <- M2$TAU2
SIG2J <- M2$SIG2J
NU <- M2$NU

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_sig2j <- coda::effectiveSize(SIG2J)
eff_nu <- coda::effectiveSize(NU)

EMC1 <- apply(X = THETA, MARGIN = 2, FUN = sd)/sqrt(eff_theta)

EMC2 <- apply(X = MU, MARGIN = 2, FUN = sd)/sqrt(eff_mu)

EMC3 <- apply(X = SIG2, MARGIN = 2, FUN = sd)/sqrt(eff_sig2)

EMC4 <- apply(X = TAU2, MARGIN = 2, FUN = sd)/sqrt(eff_tau2)

EMC5 <- apply(X = SIG2J, MARGIN = 2, FUN = sd)/sqrt(eff_sig2j)

EMC6 <- apply(X = NU, MARGIN = 2, FUN = sd)/sqrt(eff_nu)


knitr::kable(t(round(c(min(EMC1), quantile(EMC1, .25),
                       quantile(EMC1,.5), mean(EMC1),
                       quantile(EMC1, .75,), max(EMC1)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(EMC5), quantile(EMC5, .25),
                       quantile(EMC5,.5), mean(EMC5),
                       quantile(EMC5, .75,), max(EMC5)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\sigma_j^2$$")

knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\nu$$"), round(c(EMC3, EMC2, EMC4, EMC6),3)), row.names = FALSE)


# 6.3 Tercer modelo ------------------------------------------------------------

THETA <- M3$THETA
MU <- M3$MU
SIG2 <- M3$SIG2
TAU2 <- M3$TAU2
EPSJK <- M3$EPSJK
K2 <- M3$K2

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_epsjk <- coda::effectiveSize(EPSJK)
eff_k2 <- coda::effectiveSize(K2)

EMC1 <- apply(X = THETA, MARGIN = 2, FUN = sd)/sqrt(eff_theta)

EMC2 <- apply(X = MU, MARGIN = 2, FUN = sd)/sqrt(eff_mu)

EMC3 <- apply(X = SIG2, MARGIN = 2, FUN = sd)/sqrt(eff_sig2)

EMC4 <- apply(X = TAU2, MARGIN = 2, FUN = sd)/sqrt(eff_tau2)

EMC5 <- apply(X = EPSJK, MARGIN = 2, FUN = sd)/sqrt(eff_epsjk)

EMC6 <- apply(X = K2, MARGIN = 2, FUN = sd)/sqrt(eff_k2)



knitr::kable(t(round(c(min(EMC1), quantile(EMC1, .25),
                       quantile(EMC1,.5), mean(EMC1),
                       quantile(EMC1, .75,), max(EMC1)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(EMC5), quantile(EMC5, .25),
                       quantile(EMC5,.5), mean(EMC5),
                       quantile(EMC5, .75,), max(EMC5)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\xi_{j,k}$$")

knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\kappa^2$$"), round(c(EMC3, EMC2, EMC4, EMC6),3)), row.names = FALSE)

# 6.4 Cuarto modelo ------------------------------------------------------------

THETA <- M4$THETA
MU <- M4$MU
SIG2 <- M4$SIG2
TAU2 <- M4$TAU2
EPSJK <- M4$EPSJK
K2 <- M4$K2
SIG2K <- M4$SIG2K
NU <- M4$NU

eff_theta <- coda::effectiveSize(THETA)
eff_mu <- coda::effectiveSize(MU)
eff_sig2 <- coda::effectiveSize(SIG2)
eff_tau2 <- coda::effectiveSize(TAU2)
eff_epsjk <- coda::effectiveSize(EPSJK)
eff_k2 <- coda::effectiveSize(K2)
eff_sig2k <- coda::effectiveSize(SIG2K)
eff_nu <- coda::effectiveSize(NU)

EMC1 <- apply(X = THETA, MARGIN = 2, FUN = sd)/sqrt(eff_theta)

EMC2 <- apply(X = MU, MARGIN = 2, FUN = sd)/sqrt(eff_mu)

EMC3 <- apply(X = SIG2, MARGIN = 2, FUN = sd)/sqrt(eff_sig2)

EMC4 <- apply(X = TAU2, MARGIN = 2, FUN = sd)/sqrt(eff_tau2)

EMC5 <- apply(X = EPSJK, MARGIN = 2, FUN = sd)/sqrt(eff_epsjk)

EMC6 <- apply(X = K2, MARGIN = 2, FUN = sd)/sqrt(eff_k2)

EMC7<- apply(X = SIG2K, MARGIN = 2, FUN = sd)/sqrt(eff_sig2k)

EMC8 <- apply(X = NU, MARGIN = 2, FUN = sd)/sqrt(eff_nu)


knitr::kable(t(round(c(min(EMC1), quantile(EMC1, .25),
                       quantile(EMC1,.5), mean(EMC1),
                       quantile(EMC1, .75,), max(EMC1)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\theta_j$$")

knitr::kable(t(round(c(min(EMC5), quantile(EMC5, .25),
                       quantile(EMC5,.5), mean(EMC5),
                       quantile(EMC5, .75,), max(EMC5)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\xi_{j,k}$$")

knitr::kable(t(round(c(min(EMC7), quantile(EMC7, .25),
                       quantile(EMC7,.5), mean(EMC7),
                       quantile(EMC7, .75,), max(EMC7)),3)), 
             col.names = c("Mínimo", "Cuartíl 1", "Mediana",
                           "Media", "Cuartíl 3", "Máximo"),
             caption ="$$\\sigma_k^2$$")



knitr::kable(cbind(c("$$\\sigma^2$$","$$\\mu$$", "$$\\tau^2$$",
                     "$$\\kappa^2$$","$$\\nu$$"),  round(c(EMC3, EMC2, EMC4, EMC6, EMC8),3)), row.names = FALSE)

# 7. Punto 9 -------------------------------------------------------------------

# 7.1 DIC (Deviance Information Criterion) -------------------------------------

# 7.1.1 Modelo 1 ---------------------------------------------------------------
Log2Est = sum(dt.scaled(y$PUNT_GLOBAL, v, mean = rep(colMeans(M1$THETA[,1:32]), dataDepartamentos$n), 
                        sd = sqrt(mean(M1$SIG2$V1)), log = TRUE))
pDIC2 = 2 * (Log2Est - mean(M1$LL$V1))
DIC2 = -2 * Log2Est + 2 * pDIC2

# 7.1.2 Modelo 2 ---------------------------------------------------------------
Log3Est = sum(dt.scaled(y$PUNT_GLOBAL, v, mean = rep(colMeans(M2$THETA[,1:32]), dataDepartamentos$n), 
                        sd = rep(sqrt(colMeans(M2$SIG2J[,1:32])), dataDepartamentos$n), log = TRUE))
pDIC3 = 2 * (Log3Est - mean(M2$LL$V1))
DIC3 = -2 * Log3Est + 2 * pDIC3

# 7.1.3 Modelo 3 ---------------------------------------------------------------
Log4Est = sum(dt.scaled(y$PUNT_GLOBAL, v, mean = rep(colMeans(M3$EPSJK[,1:1112]), dataMunicipios$n), 
                        sd = sqrt(mean(M3$K2$V1)), log = TRUE))
pDIC4 = 2 * (Log4Est - mean(M3$LL$V1))
DIC4 = -2 * Log4Est + 2 * pDIC4

# 7.1.4 Modelo 4 ---------------------------------------------------------------
Log5Est = sum(dt.scaled(y$PUNT_GLOBAL, v, mean = rep(colMeans(M4$EPSJK[,1:1112]), dataMunicipios$n), 
                        sd = sqrt(mean(M4$K2$V1)), log = TRUE))
pDIC5 = 2 * (Log5Est - mean(M4$LL$V1))
DIC5 = -2 * Log5Est + 2 * pDIC5

DIC = rbind('Modelo 1' = c(pDIC2,DIC2),
            'Modelo 2' = c(pDIC3,DIC3),'Modelo 3' = c(pDIC4,DIC4),
            'Modelo 4' = c(pDIC5,DIC5))

colnames(DIC) = c('pDIC','DIC')

# 7.2 WAIC (Watanabe-Akaike Information Criterion) -----------------------------

# 7.2.1 Modelo 1 ---------------------------------------------------------------

colnames(M1$THETA) <- c('05', '08', '11', '13', '15', '17', '18', '19', '20',
                        '23', '25', '27', '41', '44', '47', '50', '52', '54', 
                        '63', '66', '68', '70', '73', '76', '81', '85', '86',
                        '91', '94', '95', '97', '99')



# Inicialización de variables
lppd2 = 0
pWAIC2 = 0

# Barra de progreso
Progress = txtProgressBar(min = 1, max = nrow(y), style = 3)
tictoc::tic()

for (i in 1:nrow(y)) {
  dep = y$COLE_COD_DEPTO_UBICACION[i]  # Código del departamento
  
  # Extraer valores de media y varianza para el departamento correspondiente
  mean_dep = M1$THETA[[dep]]  # Media: columna correspondiente de THETA
  sd_dep = sqrt(M1$SIG2$V1) # Desviación estándar: SIG2 (misma para todos los Departamentos)
  
  # Evaluar la densidad t de Student escalada (ajustado desde dnorm)
  dt1 = dt.scaled(y$PUNT_GLOBAL[i], v, mean = mean_dep, sd = sd_dep, log = TRUE)
  
  # Cálculo de la log-verosimilitud ponderada
  lppd2 = lppd2 + log(mean(exp(dt1)))
  
  # Cálculo del penalizador de complejidad pWAIC2
  pWAIC2 = pWAIC2 + 2 * (log(mean(exp(dt1))) - mean(dt1))
  
  # Actualizar barra de progreso
  setTxtProgressBar(Progress, i)
}

# Cierre de barra de progreso y tiempo de ejecución
close(Progress)
tictoc::toc()

# Cálculo final de WAIC para el modelo 2
WAIC2 = -2 * lppd2 + 2 * pWAIC2

# Resultado final
WAIC2

# 7.2.2 Modelo 2 ---------------------------------------------------------------

# Inicialización de variables
lppd3 = 0
pWAIC3 = 0

colnames(M2$THETA) <- c('05', '08', '11', '13', '15', '17', '18', '19', '20',
                        '23', '25', '27', '41', '44', '47', '50', '52', '54', 
                        '63', '66', '68', '70', '73', '76', '81', '85', '86',
                        '91', '94', '95', '97', '99')

colnames(M2$SIG2J) <- c('05', '08', '11', '13', '15', '17', '18', '19', '20',
                        '23', '25', '27', '41', '44', '47', '50', '52', '54', 
                        '63', '66', '68', '70', '73', '76', '81', '85', '86',
                        '91', '94', '95', '97', '99')
# Barra de progreso
Progress = txtProgressBar(min = 1, max = nrow(y), style = 3)
tictoc::tic()

for (i in 1:nrow(y)) {
  dep = y$COLE_COD_DEPTO_UBICACION[i]  # Código del departamento
  
  # Extraer valores de media y varianza para el departamento correspondiente
  mean_dep = M2$THETA[[dep]]  # Media: columna correspondiente de THETA
  sd_dep = sqrt(M2$SIG2J[[dep]]) # Desviación estándar: SIG2 (misma para todos los Departamentos)
  
  dt1 = dt.scaled(y$PUNT_GLOBAL[i], v, mean = mean_dep, sd = sd_dep, log=TRUE)
  # Cálculo de la log-verosimilitud ponderada
  lppd3 = lppd3 + log(mean(exp(dt1)))
  # Cálculo del penalizador de complejidad pWAIC2
  pWAIC3 = pWAIC3 + 2 * (log(mean(exp(dt1))) - mean(dt1))
  # Actualizar barra de progreso
  setTxtProgressBar(Progress, i)
}

# Cierre de barra de progreso y tiempo de ejecución
close(Progress)
tictoc::toc()

# Cálculo final de WAIC para el modelo 2
WAIC3 = -2 * lppd3 + 2 * pWAIC3

# Resultado final
WAIC3

# 7.2.3 Modelo 3 ---------------------------------------------------------------

# Extraer los códigos de municipio de la base 'muni'
codigos_municipios <- unique(y$COLE_COD_MCPIO_UBICACION)  # Suponiendo que la columna se llama 'CodigoMunicipio'

# Verificar que la cantidad de códigos coincide con el número de columnas en 'otra_base'
if (length(codigos_municipios) == ncol(M3$EPSJK)) {
  colnames(M3$EPSJK) <- as.character(codigos_municipios)  # Asignar nombres de columnas
} else {
  warning("El número de códigos de municipios no coincide con el número de columnas en 'otra_base'")
}

# Inicialización de variables
lppd4 = 0
pWAIC4 = 0

# Barra de progreso
Progress = txtProgressBar(min = 1, max = nrow(y), style = 3)
tictoc::tic()

for (i in 1:nrow(y)) {
  dep = y$COLE_COD_MCPIO_UBICACION[i]  # Código del departamento
  
  # Extraer valores de media y varianza para el departamento correspondiente
  mean_dep = M3$EPSJK[[dep]]  # Media: columna correspondiente de THETA
  sd_dep = sqrt(M3$K2$V1) # Desviación estándar: SIG2 (misma para todos los Departamentos)
  
  # Evaluar la densidad t de Student escalada (ajustado desde dnorm)
  dt1 = dt.scaled(y$PUNT_GLOBAL[i], v, mean = mean_dep, sd = sd_dep, log = TRUE)
  
  # Cálculo de la log-verosimilitud ponderada
  lppd4 = lppd4 + log(mean(exp(dt1)))
  
  # Cálculo del penalizador de complejidad pWAIC2
  pWAIC4 = pWAIC4 + 2 * (log(mean(exp(dt1))) - mean(dt1))
  
  # Actualizar barra de progreso
  setTxtProgressBar(Progress, i)
}

# Cierre de barra de progreso y tiempo de ejecución
close(Progress)
tictoc::toc()

# Cálculo final de WAIC para el modelo 2
WAIC4 = -2 * lppd4 + 2 * pWAIC4

# Resultado final
WAIC4

# 7.2.4 Modelo 4 ---------------------------------------------------------------

# Extraer los códigos de municipio de la base 'muni'
codigos_municipios <- unique(y$COLE_COD_MCPIO_UBICACION)  # Suponiendo que la columna se llama 'CodigoMunicipio'

# Verificar que la cantidad de códigos coincide con el número de columnas en 'otra_base'
if (length(codigos_municipios) == ncol(M3$EPSJK)) {
  colnames(M4$EPSJK) <- as.character(codigos_municipios)  # Asignar nombres de columnas
} else {
  warning("El número de códigos de municipios no coincide con el número de columnas en 'otra_base'")
}

# Inicialización de variables
lppd5 = 0
pWAIC5 = 0

# Barra de progreso
Progress = txtProgressBar(min = 1, max = nrow(y), style = 3)
tictoc::tic()

for (i in 1:nrow(y)) {
  dep = y$COLE_COD_MCPIO_UBICACION[i]  # Código del departamento
  
  # Extraer valores de media y varianza para el departamento correspondiente
  mean_dep = M4$EPSJK[[dep]]  # Media: columna correspondiente de THETA
  sd_dep = sqrt(M4$K2$V1) # Desviación estándar: SIG2 (misma para todos los Departamentos)
  
  # Evaluar la densidad t de Student escalada (ajustado desde dnorm)
  dt1 = dt.scaled(y$PUNT_GLOBAL[i], v, mean = mean_dep, sd = sd_dep, log = TRUE)
  
  # Cálculo de la log-verosimilitud ponderada
  lppd5 = lppd5 + log(mean(exp(dt1)))
  
  # Cálculo del penalizador de complejidad pWAIC2
  pWAIC5 = pWAIC5 + 2 * (log(mean(exp(dt1))) - mean(dt1))
  
  # Actualizar barra de progreso
  setTxtProgressBar(Progress, i)
}

# Cierre de barra de progreso y tiempo de ejecución
close(Progress)
tictoc::toc()

# Cálculo final de WAIC para el modelo 2
WAIC5 = -2 * lppd5 + 2 * pWAIC5

# Resultado final
WAIC5

#8. Punto 10 -------------------------------------------------------------------

Inferencia=function(a){
  estpunt=mean(a)
  cv=sqrt(var(a))/abs(mean(a))
  li=quantile(a,.025)
  ls=quantile(a,.975)
  return(t(c(estpunt,cv,li,ls)))
} #Estimacion, cv e intervalos

resultados <- data.frame(
  Modelo = c("Modelo 1", "Modelo 2", "Modelo 3", "Modelo 4"),
  do.call(rbind, lapply(list(M1$MU$V1, M2$MU$V1, M3$MU$V1, M4$MU$V1), Inferencia))
)

resultados$V2<-resultados$V2*100

# Mostrar la tabla en formato Kable
knitr::kable(resultados, digits = 3, 
             col.names = c("", "Media Posterior", "CV (%)", "L.I. (2.5%)", "L.S. (97.5%)"),
             caption = "Inferencia Bayesiana de $$\\mu$$ para los 4 modelos")

#9. Punto 11 -------------------------------------------------------------------

#9.1 Organización del ranking --------------------------------------------------

#Convertir a factor 
departamentos$Codigo=factor(departamentos$Codigo)

departamentos[departamentos == 'Norte de Santander'] = 'N. Santander'
departamentos[departamentos == 'Valle del Cauca'] = 'V. del Cauca'
nombres <- departamentos[-28,]
m <- 32
theta_hat <- colMeans(M4$THETA)
q1 <- apply(M4$THETA, 2, quantile, probs = 0.025)
q2 <- apply(M4$THETA, 2, quantile, probs = .975)

media_post <- cbind(nombres,q1, theta_hat, q2)
media_post=media_post%>%arrange(theta_hat)

#Colores de acuerdo  a la significancia:
colores = rep('#838B83', m)
colores[which(media_post[,"q2"] < 250)] = '#B22222' 
colores[which(media_post[,"q1"] > 250) ] = '#008B00'

# 9.2 Gráfico del ranking ------------------------------------------------------

#png("ranking_puntaje.png", width = 600, height = 800, res = 150)

plot(x = c(180,280), y = c(1,m), type = "n", xlab = "Puntaje", ylab = "", main = expression(italic("Ranking Bayesiano")), yaxt = "n")
axis(side = 2, at = 1:m, labels = media_post$Nombre, las = 2, cex.axis = 0.5)
abline(h = 1:m, col = "lightgray", lwd = 1)
abline(v = 250,col = "gray", lwd = 3)
for (l in 1:m) {
  segments(x0 = media_post[l,'q2'], y0 = l, x1 = media_post[l,'q1'], y1 = l, col = colores[l])
  lines(x = media_post[l,'theta_hat'], y = l, type = "p", pch = 16, cex = 0.6, col = colores[l])
}


#dev.off()

# 10. Punto 12 -----------------------------------------------------------------

# 10.1 Calculo de la matriz ---------------------------------------------------- 

indicadora_cluster <- function(cluster1, cluster2) {
  ifelse(cluster1 == cluster2, 1, 0)} #Para hallar la probabilidad de que dos deptos pertenezcan a un mismo grupo

clust=function(theta, n.clust){
  
  # Inicializar una lista para almacenar los clusters
  Omega <- vector("list", nrow(theta))
  Progress <- txtProgressBar(min = 1, max = nrow(theta), style = 3)
  # Iterar por cada fila
  for (i in 1:nrow(theta)) {
    fila_data <- matrix(theta[i, ], ncol = 1)
    resultado_kmeans <- kmeans(fila_data, centers = n.clust)
    Omega[[i]] <- resultado_kmeans$cluster
    setTxtProgressBar(Progress, i)
  }
  return(Omega)  # Retornar la lista de clusters
}


set.seed(123)
tictoc::tic()
cl <- clust(M4$THETA, 5)
tictoc::toc()  

indicadora_todas <- lapply(cl, function(clusters) {
  outer(clusters, clusters, indicadora_cluster)
})

#Calcular la matriz de similitudes

M_v <- (1 / 10000) * Reduce("+", indicadora_todas)
M_v=as.matrix(M_v)
colnames(M_v) <- nombres$Nombre
rownames(M_v) <- nombres$Nombre

media_post=media_post%>%arrange(desc(theta_hat))
nuevo_orden <- list(media_post$Nombre)


M_v <- M_v[unlist(nuevo_orden), unlist(nuevo_orden)]

# 10.2 K-means -----------------------------------------------------------------

media_post$Codigo <- as.character(media_post$Codigo)
media_post$Codigo[media_post$Codigo=="5"]="05"
media_post$Codigo[media_post$Codigo=="8"]="08"

media_post<-media_post%>%mutate(cl= kmeans(theta_hat,5)$cluster)

# 10.3 Gráfico matriz ----------------------------------------------------------

#png("matriz_incidencia_dpto.png", width = 1000, height = 1100, res = 150)
corrplot::corrplot(corr = M_v, 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color", 
                   tl.pos = "lt",
                   tl.cex = 0.8,
                   tl.col = 'black')
title(main = 'Matriz de incidencia por departamentos', line = 1.8, cex.main = 1.6)

#dev.off()

# 10.4 Gráfico mapa k-means ----------------------------------------------------

m1<-inner_join(x = shp, y = media_post, by = c("DPTO"="Codigo")) %>% 
  select(DPTO, theta_hat, geometry, Nombre, cl) %>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = as.factor(cl)), size = 0.125, color = "gray30", linetype="solid") +
  scale_fill_manual(values = c("#F38D30", "#C590B3", "#FFA07A", "#FFFACD", "#95B958")) +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  # geom_sf_text(aes(label=Nombre),col="black",
  #fontface="bold",size=3,fun.geometry=function(x) sf::st_centroid(x))+
  labs(x="Longitud",y="Latitud",title="Agrupamiento por media posterior",fill="Cluster")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

m1

#ggsave("agrupamiento_posterior_dpto.png", 
#       plot = m1, 
#       width = 10, height = 8, dpi = 300)

# 11. Punto 13 -----------------------------------------------------------------

# 11.1 Predicción --------------------------------------------------------------

THETA <- M4$THETA

colnames(THETA)<-nombres$Codigo
colnames(THETA)[1]="05"
colnames(THETA)[2]="08"

predictivas <- THETA[,pobreza_m$Codigo]

resto <- THETA[, setdiff(colnames(THETA), pobreza_m$Codigo)]


predicciones <- matrix(NA, ncol = ncol(resto), nrow = nrow(resto))

Progress = txtProgressBar(min = 2, max = nrow(THETA)+1, style = 3)
for (i in 1:nrow(THETA))
{
  newdata <- as.data.frame(t(resto[i,]))
  colnames(newdata) <- "Media"
  
  X <- cbind(t(predictivas[i,]),pobreza_m$pm)
  colnames(X)<-c("Media", "PM")
  X <- as.data.frame(X)
  
  fit <- lm(PM ~ Media, data = X)
  
  predicciones[i,] <- predict(fit, newdata)
  setTxtProgressBar(Progress,i)
} 

colnames(predicciones) <- colnames(resto)

predicciones <- as.data.frame(predicciones)

tabla<-matrix(NA, nrow = ncol(predicciones), ncol = 4)
colnames(tabla)<- c("Departamento", "Media Posterior","LI",  "LS")

# 11.2 Inferencia -------------------------------------------------------------- 

tabla[,1]<- departamentos$Nombre[departamentos$Codigo %in% colnames(predicciones)]
tabla[,3]<- t(apply(predicciones,2, quantile, prob=.025))
tabla[,4]<- t(apply(predicciones,2, quantile, prob=.975))
tabla[,2]<- t(apply(predicciones,2, mean))

tabla <- as.data.frame(tabla)
tabla$LI<-as.numeric(tabla$LI)
tabla$LS<-as.numeric(tabla$LS)
tabla$`Media Posterior`<-as.numeric(tabla$`Media Posterior`)

tabla[,c(2:4)] <- round(tabla[,c(2:4)],3) 
tabla<-tabla%>%arrange(desc(`Media Posterior`))
knitr::kable(tabla)

# 11.3 Mapa --------------------------------------------------------------------

tabla<- inner_join(tabla, departamentos, by =c("Departamento"="Nombre"))%>%select(Codigo, everything(tabla))

m1<-left_join(x = shp, y = tabla, by = c("DPTO"="Codigo"))%>%
  select(DPTO ,`Media Posterior`, geometry)%>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = `Media Posterior`), size = 0.125, color = "gray30", linetype="solid") +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    low = "#FC7D69",   # Color claro para valores bajos
    high = "#9C0824",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Media \n Posterior de \n Departamentos \n sin IPM")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

m1

#ggsave("mapa_punto13.png", 
#       plot = m1, 
#       width = 10, height = 8, dpi = 300)

# 12. Punto 14 -----------------------------------------------------------------

# 12.1 Calculo de la matriz ----------------------------------------------------

theta_hat <- colMeans(M4$EPSJK)
muni<-y%>%distinct(COLE_COD_MCPIO_UBICACION)

media_post <- cbind(muni, theta_hat)
media_post=media_post%>%arrange(desc(theta_hat))

set.seed(123)
tictoc::tic()
cl <- clust(M4$EPSJK, 8)
tictoc::toc()  

a<-matrix(0, ncol=ncol(M4$EPSJK), nrow=ncol(M4$EPSJK))

output_file <- "M_incidencia.csv"

Progress = txtProgressBar(min = 2, max = nrow(M4$EPSJK)+1, style = 3)

tictoc::tic()

for (i in 1:10000)
{
  a <- a + (1/10000)*(outer(cl[[i]],cl[[i]],indicadora_cluster))
  setTxtProgressBar(Progress,i)
}

tictoc::toc()

colnames(a)<-muni$COLE_COD_MCPIO_UBICACION
rownames(a)<-colnames(a)

nuevo_orden <- list(media_post$COLE_COD_MCPIO_UBICACION)

M_v <- a[unlist(nuevo_orden), unlist(nuevo_orden)]

# 12.2 Gráfico matriz ----------------------------------------------------------

#M_v<-read_delim("M_incidencia.csv", col_names = TRUE)

#M_v <- as.matrix(M_v)

#png("matriz_incidencia_muni.png", width = 1000, height = 1100, res = 150)

corrplot::corrplot(corr = M_v, 
                   is.corr = FALSE,
                   addgrid.col = NA, 
                   method = "color",
                   tl.pos = "n", 
                   tl.cex = 0.8,
                   tl.col = 'black')
title(main = 'Matriz de incidencia por municipios', line = 1.8, cex.main = 1.6)

#dev.off()

#img1 <- image_read("matriz_incidencia_muni.png")

#img1

# 12.3 Gráfico Mapa k-means-----------------------------------------------------

media_post<-media_post%>%mutate(cluster=kmeans(theta_hat,8)$cluster)
colnames(media_post)[1]<- "MUNI"

m1<-left_join(x = shp1, y = media_post, by = c("MUNI"="MUNI")) %>% 
  select(MUNI, cluster, geometry) %>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = as.factor(cluster)), size = 0.0125, color = "gray30", linetype="solid") +
  #scale_fill_manual(values = c("#F38D30", "#C590B3", "#FFA07A", "#FFFACD", "#95B958")) +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_manual(values=c("#F38D30", "#C590B3", "#FFA07A", "#FFFACD", "#95B958","#00CD66", "#EEEE00","#00868B")
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  #geom_sf_text(aes(label=Nombre),col="black",
  #            fontface="bold",size=3,fun.geometry=function(x) sf::st_centroid(x))+
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Clusters")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )

m1

#ggsave("agrupamiento_posterior_muni.png", 
#       plot = m1, 
#       width = 10, height = 8, dpi = 300)

# 13. Punto 15 -----------------------------------------------------------------

# 13.1 Prediccion --------------------------------------------------------------

EPSJK <- M4$EPSJK

colnames(EPSJK)<-muni$COLE_COD_MCPIO_UBICACION

predictivas <- est_edu$MUNI[est_edu$MUNI %in% colnames(EPSJK)]

resto <- EPSJK[, setdiff(colnames(EPSJK), est_edu$MUNI)]

predictivas <- EPSJK[,predictivas]

cobertura <- est_edu%>%filter(MUNI %in% colnames(predictivas))

predicciones <- matrix(NA, ncol = ncol(resto), nrow = nrow(resto))

Progress = txtProgressBar(min = 2, max = nrow(EPSJK)+1, style = 3)
for (i in 1:nrow(EPSJK))
{
  newdata <- as.data.frame(t(resto[i,]))
  colnames(newdata) <- "Media"
  
  X <- cbind(t(predictivas[i,]),cobertura$CNS)
  colnames(X)<-c("Media", "CNS")
  X <- as.data.frame(X)
  
  fit <- lm(CNS ~ Media, data = X)
  
  predicciones[i,] <- predict(fit, newdata)
  setTxtProgressBar(Progress,i)
} 

colnames(predicciones) <- colnames(resto)

predicciones <- as.data.frame(predicciones)

# 13.2 Inferencia --------------------------------------------------------------

tabla<-matrix(NA, nrow = ncol(predicciones), ncol = 4)
colnames(tabla)<- c("Municipio", "LI", "Media Posterior", "LS")

tabla[,1]<- c("Belén de Bajirá","Mapiripana")
tabla[,2]<- t(apply(predicciones,2, quantile, prob=.025))
tabla[,4]<- t(apply(predicciones,2, quantile, prob=.975))
tabla[,3]<- t(apply(predicciones,2, mean))

tabla <- as.data.frame(tabla)
tabla$LI<-as.numeric(tabla$LI)
tabla$LS<-as.numeric(tabla$LS)
tabla$`Media Posterior`<-as.numeric(tabla$`Media Posterior`)

tabla[,c(2:4)] <- round(tabla[,c(2:4)],3) 
tabla<-tabla%>%arrange(desc(`Media Posterior`))
knitr::kable(tabla)

# 13.3 Mapa --------------------------------------------------------------------

tabla<-cbind(tabla, Codigo=colnames(resto))

# Mapa de Chocó 

m1<-left_join(x = shp1, y = est_edu, by = c("MUNI"="MUNI"))%>%
  filter(DPTO_CCDGO=="27")%>%
  select(MUNI ,CNS, geometry)%>%
  ggplot() +
  geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
  scale_fill_gradient(
    high = "#3B8C4D",   # Color claro para valores bajos
    low = "#D4C95F",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  ) +
  labs(x="Longitud",y="Latitud",title="Chocó",fill="Cobertura \n Neta \n Secundaria \n en 2022")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="gray90"
    ))

# Mapa de Guainía 

m2<- left_join(x = shp1, y = est_edu, by = c("MUNI"="MUNI"))%>%
  filter(DPTO_CCDGO=="94")%>%
  select(MUNI ,CNS, geometry)%>%
  ggplot() +
  geom_sf(aes(fill = CNS), size = 0.0125, color = "gray30", linetype="solid") +
  scale_fill_gradient(
    high = "#3B8C4D",   # Color claro para valores bajos
    low = "#D4C95F",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  ) +
  labs(x="Longitud",y="Latitud",title="Guainía",fill="Cobertura \n Neta \n Secundaria \n en 2022")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="gray90")
  )

# Mapa combinado 

m3<-ggarrange(m1,m2)

m3

#ggsave("mapas_punto15.png", plot = m3, width = 7, height = 4, dpi = 300)

# 14. Punto 16 -----------------------------------------------------------------

# 14.1 Calculo outliers --------------------------------------------------------

# Paso 1: Determinar el percentil 95 de POSTSIJK
umbral_atipico <- quantile(M4$POSTSIJK$V1, 0.95)

# Paso 2: Clasificar observaciones como atípicas
M4$POSTSIJK[, atipico := ifelse(V1 > umbral_atipico, 1, 0)]

# Paso 3: Agregar códigos de departamento desde la tabla `y`
M4$POSTSIJK[, departamento := y$COLE_COD_DEPTO_UBICACION]

# Paso 4: Calcular la proporción de observaciones atípicas por departamento
proporcion_atipicos <- M4$POSTSIJK[, .(prop_outliers = mean(atipico)*100), by = departamento]

# Paso 5: Obtener el top 5 Departamentos con mayor proporción de outliers
top5_Departamentos <- as.data.frame(proporcion_atipicos[order(-prop_outliers)][1:5])

Departamentos$CodDepartamento <- as.numeric(Departamentos$CodDepartamento)
top5_Departamentos$departamento <- as.numeric(top5_Departamentos$departamento)

top5_Departamentos1 <- sqldf("
    SELECT t.departamento as cod_dpto, t.prop_outliers, d.Departamento
    FROM top5_Departamentos as t
    INNER JOIN Departamentos as d ON t.departamento = d.CodDepartamento
")
# Convertir cod_dpto a carácter y asegurarse de que tenga dos dígitos
top5_Departamentos1$cod_dpto <- sprintf("%02d", as.numeric(top5_Departamentos1$cod_dpto))


top5_Departamentos1$cod_dpto <- as.character(top5_Departamentos1$cod_dpto)


# 14.2 Mapa --------------------------------------------------------------------

m1<-left_join(x = shp, y = top5_departamentos1, by = c("DPTO"="cod_dpto"))%>%
  select(DPTO ,prop_outliers, geometry)%>%
  ggplot() +
  geom_sf(data=mundocol, col="white")+
  geom_sf(aes(fill = prop_outliers), size = 0.125, color = "gray30", linetype="solid") +
  coord_sf(xlim=c(box$xmin,box$xmax),ylim=c(box$ymin,box$ymax),expand=FALSE) +
  scale_fill_gradient(
    low = "#FC7D69",   # Color claro para valores bajos
    high = "#BD1100",   # Color oscuro para valores altos
    na.value = "#EEE9E9"
  )  +
  annotate("text", x=c(-74.5,-68,-78,-69,-78.5), y=c(-2.5,0,-1,9,9), colour="blue",
           label=c("Perú","Brasil","Ecuador","Venezuela","Panamá")) +
  labs(x="Longitud",y="Latitud",title="Colombia",fill="Proporcion \n de outliers (%)")+
  theme(
    plot.title = element_text(hjust = 0.5,),
    panel.background=element_rect(fill="lightblue")
  )


ggsave("Top_5Prop.png", plot = m1,width = 10, height = 8, dpi = 300)

# 15. Punto 17 -----------------------------------------------------------------

# 15.1 Calculo -----------------------------------------------------------------

pppValues = matrix(0, nrow = nrow(dataMunicipios), ncol = 1)
colnames(pppValues) = c('Media')

desv = sqrt(M4$K2)

Progress = txtProgressBar(min = 0, max = 10000, style = 3)
tictoc::tic()
for (i in 1:10000){
  desvi = desv[i]
  for (j in 1:nrow(dataMunicipios)){
    muestra = rnorm(n = dataMunicipios$n[j],
                    mean = M4$EPSJK[[i, dataMunicipios$CodMunicipio[j]]],
                    sd = desvi$V1)
    pppValues[j,"Media"] =  pppValues[j,"Media"] + (mean(muestra) < dataMunicipios$mean[j])
    setTxtProgressBar(Progress,i)
  }
}
close(Progress)
tictoc::toc()
beepr::beep()

pppValues1 = pppValues/10000

# 15.2 Visualización -----------------------------------------------------------

xd=ggplot(pppValues1, aes(x = "Media", y = Media)) +  # "Media" para que tenga etiqueta
  geom_boxplot(fill = "#C1CDC1", color = "black", alpha = 0.6, width = 0.4) +  # Color caja y contorno
  geom_jitter(width = 0.2, alpha = 0.3, color = "#838B83", size = 1.2) +  # Puntos dispersos en color verde azulado
  theme_minimal() +
  labs(
    title = "Bondad de ajuste",
    x = "Valores PPP para el modelo 4 por Departamentos",
    y = "ppp"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11)
  )

#ggsave("xd.png", 
#       plot = xd, 
#       width = 14, height = 8, dpi = 300)
