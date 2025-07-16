setwd("C:/Users/Lenovo/OneDrive - Universidad Nacional de Colombia/Documentos/Caso Bayesiana 2025 - I")

library(readxl)
library(dplyr)

sb11_22 <- read.csv("Examen_Saber_11_2022_2_Entrenamiento.txt", sep = ";")

educ <- read_excel("PANEL_DE_EDUCACION(2022).xlsx")

conf_vio <- read_excel("PANEL_CONFLICTO_Y_VIOLENCIA(2022).xlsx")

c_gen <- read_excel("PANEL_CARACTERISTICAS_GENERALES(2022).xlsx")

PIB <- read_excel("DANE - PIB.xlsx", sheet = 4, range = "A9:U43")

LAFT <- read_excel("LAFT.xlsx")

MOE <- data.frame(
  Depto = c(
    "Antioquia", "Cauca", "Chocó", "Norte de Santander", "Nariño", "Arauca", "Bolívar", "Putumayo", "Córdoba", "Valle del Cauca", 
    "Caquetá", "Meta", "Santander", "Cesar", "Atlántico", "Guaviare", "Sucre", "Tolima", "La Guajira", "Magdalena", "Casanare", 
    "Risaralda", "Huila", "Guainía", "Vichada", "Bogotá D.C.", "Caldas", "Quindío", "Boyacá", "Cundinamarca", "Amazonas", 
    "San Andrés", "Vaupés"
  ),
  Riesgo_Extremo = c(
    18, 14, 12, 11, 7, 6, 5, 5, 4, 3,
    3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0
  ),
  Riesgo_Alto = c(
    19, 12, 13, 7, 10, 1, 12, 4, 5, 8,
    6, 11, 3, 7, 2, 2, 13, 5, 4, 3, 3,
    3, 2, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 0
  ),
  Riesgo_Medio = c(
    23, 6, 2, 1, 13, 0, 4, 1, 7, 9,
    6, 5, 1, 8, 1, 0, 5, 4, 5, 3, 1,
    1, 5, 4, 1, 0, 0, 0, 2, 2, 0,
    0, 0
  ),
  Total_en_riesgo = c(
    60, 32, 27, 19, 30, 7, 21, 10, 16, 20,
    15, 18, 6, 16, 4, 3, 18, 9, 9, 6, 4,
    4, 7, 5, 2, 1, 1, 1, 2, 2, 0,
    0, 0
  ),
  Porcentaje_municipios_en_riesgo = c(
    48.00, 76.19, 90.00, 47.50, 46.88, 100.00, 45.65, 76.92, 53.33, 47.62,
    93.75, 62.07, 6.90, 64.00, 17.39, 75.00, 69.23, 19.15, 60.00, 20.00, 21.05,
    28.57, 18.92, 55.56, 50.00, 100.00, 3.70, 8.33, 1.63, 1.72, 0.00,
    0.00, 0.00
  )
)

departamentos <- tibble::tibble(
  Nombre = c("Amazonas", "Antioquia", "Arauca", "Atlántico", "Bogotá D.C.", "Bolívar",
             "Boyacá", "Caldas", "Caquetá", "Casanare", "Cauca", "Cesar", "Chocó",
             "Córdoba", "Cundinamarca", "Guainía", "Guaviare", "Huila", "La Guajira",
             "Magdalena", "Meta", "Nariño", "Norte de Santander", "Putumayo", "Quindío",
             "Risaralda", "San Andrés", "Santander", "Sucre", "Tolima", "Valle del Cauca",
             "Vaupés", "Vichada"),
  Código = c("91", "5", "81", "8", "11", "13", "15", "17", "18", "85", "19", "20", "27",
             "23", "25", "94", "95", "41", "44", "47", "50", "52", "54", "86", "63", "66",
             "88", "68", "70", "73", "76", "97", "99")
)


# Filtrado de datos


sb11<- sb11_22 %>% filter(estu_nacionalidad == "COLOMBIA" &
                   estu_pais_reside == "COLOMBIA" &
                   cole_cod_depto_ubicacion != 88 &
                   if_all(c(fami_educacionmadre,
                            fami_tienecomputador,
                            fami_tieneinternet,
                            fami_numlibros,
                            fami_estratovivienda,
                            estu_tieneetnia,
                            punt_global,
                            cole_cod_depto_ubicacion,
                            cole_cod_mcpio_ubicacion),
                          ~ !is.na(.))) %>% 
              select(cole_cod_depto_ubicacion,
                     cole_cod_mcpio_ubicacion,
                     punt_global,
                     fami_educacionmadre,
                     fami_tienecomputador,
                     fami_tieneinternet,
                     fami_numlibros,
                     fami_estratovivienda,
                     estu_tieneetnia)

#Variables a Nivel Municipal ---------------------------------------------------


LAFT <- LAFT %>% arrange(CODIGO)

LAFT$codmpio <- (c_gen %>% filter(ano==2018) %>%  
                arrange(depto,municipio) %>% distinct(codmpio))

LAFT$codmpio <- as.integer(unlist(LAFT$codmpio))

c_gen <- c_gen %>% mutate(codmpio = as.double(codmpio))
LAFT <- LAFT %>% mutate(codmpio = as.double(unlist(codmpio)))

municipales <- inner_join(LAFT, (c_gen %>% filter(ano==2018)), by=c("codmpio")) %>% 
               inner_join(educ %>% filter(ano==2021) %>% mutate(doc_alum = docen_total/alumn_total), by=("codmpio")) %>%
               select(codmpio, doc_alum, nbi, RISK_VICTIM_2022)

#Variables a nivel departamental -----------------------------------------------

MOE <-inner_join(MOE, departamentos, by=c("Depto"="Nombre")) %>% dplyr::select(`Código`,Porcentaje_municipios_en_riesgo)

proprul<-c_gen %>% filter(ano==2018) %>% group_by(coddepto) %>% summarise(porc_rul = 100*sum(pobl_rur)/sum(pobl_tot)) 


PIB <- PIB %>% mutate(`Código Departamento (DIVIPOLA)`= case_when(
                      `Código Departamento (DIVIPOLA)`=="05" ~"5",
                      `Código Departamento (DIVIPOLA)`=="08" ~"8",
                      .default = `Código Departamento (DIVIPOLA)`
))

proprul$coddepto <- as.character(proprul$coddepto)

departamentales <-inner_join(MOE, PIB, by=c("Código"="Código Departamento (DIVIPOLA)")) %>%
  inner_join(proprul, by=c("Código"="coddepto")) %>% 
  dplyr::select(Cod=`Código`, PIB=`2022`,porc_riesgo=Porcentaje_municipios_en_riesgo, porc_rul)


#Reescritura de datos ---------------------------------------------------------

sb11<- sb11 %>% mutate(edu_madre = case_when(
  fami_educacionmadre %in% c("Educación profesional completa",
                             "Postgrado") ~ 1,
  fami_educacionmadre == "" ~NA,
  .default = 0),
  edu_madre = as.factor(edu_madre)
)

sb11 <- sb11 %>% mutate(comp = case_when(
  fami_tienecomputador == "Si" ~ 1,
  fami_tienecomputador == "No" ~ 0,
  .default = NA),
  comp = as.factor(comp)
)

sb11 <- sb11 %>% mutate(internet = case_when(
  fami_tieneinternet == "Si" ~ 1,
  fami_tieneinternet == "No" ~ 0,
  .default = NA),
  internet = as.factor(internet)
)

sb11 <- sb11 %>% mutate(libros = case_when(
  fami_numlibros =="" ~ NA,
  .default = fami_numlibros
),
libros = as.factor(libros)
)

sb11 <- sb11 %>% mutate(estrato = case_when(
  fami_estratovivienda =="" ~ NA,
  .default = fami_estratovivienda
),
estrato = as.factor(estrato)

)

sb11 <- sb11 %>% mutate(etnia = case_when(
  estu_tieneetnia == "Si" ~ 1,
  estu_tieneetnia == "No" ~ 0,
  .default = NA),
  etnia = as.factor(etnia)
)

sb11<-sb11 %>% filter(
  if_all(c(edu_madre,
           comp,
           internet,
           libros,
           estrato,
           etnia,
           punt_global,
           cole_cod_depto_ubicacion,
           cole_cod_mcpio_ubicacion),
         ~ !is.na(.))
)


sb11 <- sb11[,c(1:3,10:15)]

sb11 <- sb11 %>% arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)


library(data.table)

fwrite(sb11, "Datos.txt", sep = ";")
fwrite(departamentales, "Departamentales.txt", sep = ";")
fwrite(municipales, "Municipales.txt", sep=";")



# Para la base de prueba ----

prueba <- read.csv("Examen_Saber_11_2022_2_Prueba.txt", sep = ";")

prueba<- prueba %>% filter(estu_nacionalidad == "COLOMBIA" &
                            estu_pais_reside == "COLOMBIA" &
                            cole_cod_depto_ubicacion != 88 &
                            if_all(c(fami_educacionmadre,
                                     fami_tienecomputador,
                                     fami_tieneinternet,
                                     fami_numlibros,
                                     fami_estratovivienda,
                                     estu_tieneetnia,
                                     punt_global,
                                     cole_cod_depto_ubicacion,
                                     cole_cod_mcpio_ubicacion),
                                   ~ !is.na(.))) %>% 
  select(cole_cod_depto_ubicacion,
         cole_cod_mcpio_ubicacion,
         punt_global,
         fami_educacionmadre,
         fami_tienecomputador,
         fami_tieneinternet,
         fami_numlibros,
         fami_estratovivienda,
         estu_tieneetnia)

prueba<- prueba %>% mutate(edu_madre = case_when(
  fami_educacionmadre %in% c("Educación profesional completa",
                             "Técnica o tecnológica completa",
                             "Postgrado") ~ 1,
  fami_educacionmadre == "" ~NA,
  .default = 0),
  edu_madre = as.factor(edu_madre)
)

prueba <- prueba %>% mutate(comp = case_when(
  fami_tienecomputador == "Si" ~ 1,
  fami_tienecomputador == "No" ~ 0,
  .default = NA),
  comp = as.factor(comp)
)

prueba <- prueba %>% mutate(internet = case_when(
  fami_tieneinternet == "Si" ~ 1,
  fami_tieneinternet == "No" ~ 0,
  .default = NA),
  internet = as.factor(internet)
)

prueba <- prueba %>% mutate(libros = case_when(
  fami_numlibros =="" ~ NA,
  .default = fami_numlibros
),
libros = as.factor(libros)
)

prueba <- prueba %>% mutate(estrato = case_when(
  fami_estratovivienda =="" ~ NA,
  .default = fami_estratovivienda
),
estrato = as.factor(estrato)

)

prueba <- prueba %>% mutate(etnia = case_when(
  estu_tieneetnia == "Si" ~ 1,
  estu_tieneetnia == "No" ~ 0,
  .default = NA),
  etnia = as.factor(etnia)
)

prueba<-prueba %>% filter(
  if_all(c(edu_madre,
           comp,
           internet,
           libros,
           estrato,
           etnia,
           punt_global,
           cole_cod_depto_ubicacion,
           cole_cod_mcpio_ubicacion),
         ~ !is.na(.))
)

prueba <- prueba[,c(1:3,10:15)]

sb11 <- sb11 %>% arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)

fwrite(prueba, "Prueba.txt", sep=";")
