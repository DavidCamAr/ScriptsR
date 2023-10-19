###################################################################################
# David Campos Arellano                                                           #
# 27/07/23                                                                        #
# Obejtivo: Representar la proporcion de inserciones de TEs en los 23 cromosomas  #
# por medio de un percent stacked bar                                             #
###################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("tidyverse", #la funcion p_load viene del paquete "pacman"
       "ggplot2",
       "ggsci",
       "ggrepel",
       "gridExtra")

##directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

# Contamos cuantos TEs hay en cada cromosoma, usando la funcion aggregate(), con el 
# parametro FUN = length. El conteo se hara solo con las variables de interes

##Contamos los TEs por cromosoma 
c1 <- aggregate(TEs_data.rc$TE, by = list(TEs_data.rc$Cromosoma, 
                                          TEs_data.rc$TE,
                                          TEs_data.rc$Condicion, 
                                          TEs_data.rc$Tejido, 
                                          TEs_data.rc$RC), 
                FUN = length)

## Etiquetamos columnas
colnames(c1) <- c("Cromosoma","TE",
                  "Condicion","Tejido",
                  "RC",
                  "Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(c1)

# Cromosoma  TE Condicion Tejido      RC Count_TE
# 1         1 ALU     Tumor  Liver 1549108        4
# 2        10 ALU     Tumor  Liver 1549108        4
# 3        11 ALU     Tumor  Liver 1549108        3
# 4        12 ALU     Tumor  Liver 1549108        8
# 5        13 ALU     Tumor  Liver 1549108        5
# 6        14 ALU     Tumor  Liver 1549108        1

# Ordenamos por cromosoma ya que no se dejan ordenar nombrando como factor a la columna, 
# primero lo tomamos como caracter
orden_chrs <- c(as.character(1:22), "X", "Y")
c1$Cromosoma <- factor(c1$Cromosoma, levels = orden_chrs)
c1 <- c1[order(c1$Cromosoma), ]

# Calculamos las RPM en funcion del total de TEs en cada cromosoma contados  
c1 <- c1 %>% 
  mutate(RPM = (Count_TE / RC) * 10^6)

# Calculamos el promedio de RPM 
c2 <- aggregate(RPM ~ Cromosoma + TE + Condicion + Tejido, data = c1, FUN = mean)

##Para anular observacion sin valores

#Sacamos combinaciones entre observaciones
tejidos <- unique(c2$Tejido)
condicion <- unique(c2$Condicion)
cromosomas <- unique(c2$Cromosoma)
tipos_TE <- unique(c2$TE)
combinaciones <- expand.grid(Cromosoma = cromosomas, Tejido = tejidos,
                             Condicion = condicion, TE = tipos_TE)

#Juntamos combinaciones al d.f de trabajo
c2 <- merge(combinaciones, c2, by = c("Cromosoma","Tejido","Condicion","TE"),
                         all.x = TRUE, all.y = FALSE)
##Agregamos valores =NA a los valores faltantes
c2$RPM[is.na(c2$RPM)] <- NA

# Ajustamos Levels
c2$Tejido <- recode_factor(c2$Tejido,
                                        "Placenta" = "Trofoblastos",
                                        "Breast" = "Mama",
                                        "Colon" = "Colon",
                                        "Liver" = "Hígado")
c2$TE <- recode_factor(c2$TE,
                                    "ALU" = "Alu",
                                    "LINE1" = "L1",
                                    "SVA" = "SVA")


##Filtro para quedarnos con cáncer y trofoblastos
c3 <- subset(c2, !(Tejido %in% c("Mama", "Colon", "Hígado") & Condicion == "Sano"))

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(c3)

# 'data.frame':	360 obs. of  5 variables:
#   $ Cromosoma: Factor w/ 24 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Tejido   : Factor w/ 4 levels "Trofoblastos",..: 2 2 2 3 3 3 4 4 4 1 ...
# $ Condicion: Factor w/ 2 levels "Sano","Tumor": 2 2 2 2 2 2 2 2 2 1 ...
# $ TE       : Factor w/ 3 levels "Alu","L1","SVA": 1 2 3 1 2 3 1 2 3 1 ...
# $ RPM      : num  0.628 0.182 0.226 1.234 0.287 ...

##graficamos

# Calculamos los valores de RPM en porcentaje
c3 <- c3 %>%
  group_by(Tejido, Cromosoma) %>%
  mutate(Porcentaje = RPM / sum(RPM, na.rm = TRUE))


bar1 <- c3 %>%
  ggplot(aes(x = Cromosoma,
                          y = RPM,
                          fill = TE))+
  geom_col(color = "white", 
           position = "fill",
           width = 1)+
  geom_text(show.legend = NULL,
             data=subset(c3, Tejido == "Hígado" & Cromosoma %in% c("X","20")),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                             accuracy = 1)),
            position = position_fill(vjust = 0.5),  # Ajusta la posición de las etiquetas
            size = 4,
            color = "white",
            fontface = "bold") +
  geom_text(show.legend = NULL,
             data=subset(c3, Tejido == "Mama" & Cromosoma == "15"),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                      accuracy = 1)),
            position = position_fill(vjust = 0.5),
            size = 4,
            color = "white",
            fontface = "bold")+
  geom_text(show.legend = NULL,
             data=subset(c3, Tejido == "Trofoblastos" & Cromosoma == "21"),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                      accuracy = 1)),
            position = position_fill(vjust = 0.5),
            size = 4,
            color = "white",
            fontface = "bold")+
  facet_grid(Tejido ~ .)+
  labs(x = "Cromosomas",
       y = "Proporción")+
  theme_classic()+
  theme(text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 20, vjust = 1.5),
        axis.title.x = element_text(size = 20, hjust = 0.5),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_blank(),
         legend.key.size = unit(1, "cm"),
         legend.position = "bottom",
  legend.text = element_text(size = 20))+
  guides(fill = guide_legend(title = NULL))+
  scale_y_continuous(labels = scales::percent_format(scale = 100, suffix = "%",
                                                     accuracy = 1),
                     breaks = c(0.25, 0.5, 0.75, 1))+
  scale_fill_d3()
bar1 
  
##Guardamos en el directorio de trabajo establecido
ggsave(filename = "stackedbar_chrs_tumor.jpg",
       plot = bar1,
       width = 34, height = 19, units = "cm",
       dpi = 300)
