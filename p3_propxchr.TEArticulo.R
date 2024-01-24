###################################################################################
# David Campos Arellano                                                           #
# 27/07/23                                                                        #
# Obejtivo: Representar la proporcion de inserciones de TEs en los 23 cromosomas  #
# por medio de un percent stacked bar                                             #
###################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("tidyverse", #la funcion p_load viene del paquete "pacman"
       "ggplot2",
       "tidyr",
       "ggsci",
       "ggrepel",
       "gridExtra")

##directorio de trabajo
setwd("~/Documents/Scripts_articulo/ScriptsR/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

# Contamos cuantos TEs hay en cada cromosoma, usando la funcion aggregate(), con el 
# parametro FUN = length. El conteo se hara solo con las variables de interes

c0 <- TEs_data.rc
c0$Tejido_Condicion <- paste(c0$Tejido, c0$Condicion, sep = "_")
c0$Tejido_Condicion <- factor(c0$Tejido_Condicion, levels = unique(c0$Tejido_Condicion))

##Contamos los TEs por cromosoma 
c1 <- aggregate(c0$TE, by = list(c0$Cromosoma, 
                                          c0$TE,
                                          c0$Tejido_Condicion,
                                          c0$RC), 
                FUN = length)

## Etiquetamos columnas
colnames(c1) <- c("Cromosoma","TE",
                  "Tejido_Condicion",
                  "RC",
                  "Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(c1)

# Cromosoma  TE Tejido_Condicion      RC Count_TE
# 1         1 ALU      Liver_Tumor 1549108        4
# 2        10 ALU      Liver_Tumor 1549108        4
# 3        11 ALU      Liver_Tumor 1549108        3
# 4        12 ALU      Liver_Tumor 1549108        8
# 5        13 ALU      Liver_Tumor 1549108        5
# 6        14 ALU      Liver_Tumor 1549108        1

# Ordenamos por cromosoma ya que no se dejan ordenar nombrando como factor a la columna, 
# primero lo tomamos como caracter
orden_chrs <- c(as.character(1:22), "X", "Y")
c1$Cromosoma <- factor(c1$Cromosoma, levels = orden_chrs)
c1 <- c1[order(c1$Cromosoma), ]

# Calculamos las RPM en funcion del total de TEs en cada cromosoma contados  
c2 <- c1 %>% 
  mutate(RPM = (Count_TE / RC) * 10^6)

# Calculamos el promedio de RPM 
# c2 <- aggregate(RPM ~ Cromosoma + TE + Condicion + Tejido, data = c1, FUN = mean)

##Para anular observacion sin valores

#Sacamos combinaciones entre observaciones
tejido_condicion <- unique(c2$Tejido_Condicion)
cromosomas <- unique(c2$Cromosoma)
tipos_TE <- unique(c2$TE)
combinaciones <- expand.grid(Cromosoma = cromosomas,
                             Tejido_Condicion = tejido_condicion, TE = tipos_TE)

#Juntamos combinaciones al d.f de trabajo
c2 <- merge(combinaciones, c2, by = c("Cromosoma","Tejido_Condicion","TE"),
                         all.x = TRUE, all.y = FALSE)
##Agregamos valores =NA a los valores faltantes
c2$RPM[is.na(c2$RPM)] <- 0

# Ajustamos Levels
c2$Tejido_Condicion <- recode_factor(c2$Tejido_Condicion,
                                        "Trofoblastos_Sano" = "Trophoblast",
                                        "Breast_Tumor" = "Breast Tumor",
                                        "Breast_Sano" = "Healthy Breast",
                                        "Colon_Tumor" = "Colon Tumor",
                                        "Colon_Sano" = "Healthy Colon",
                                        "Liver_Tumor" = "Liver Tumor",
                                        "Liver_Sano" = "Healthy Liver")
c2$TE <- recode_factor(c2$TE,
                                    "ALU" = "Alu",
                                    "LINE1" = "L1",
                                    "SVA" = "SVA")


##Filtro para quedarnos con cáncer y trofoblastos
# c3 <- subset(c2, !(Tejido %in% c("Mama", "Colon", "Hígado") & Condicion == "Sano"))

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(c2)

# 'data.frame':	360 obs. of  5 variables:
#   $ Cromosoma: Factor w/ 24 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Tejido   : Factor w/ 4 levels "Trofoblastos",..: 2 2 2 3 3 3 4 4 4 1 ...
# $ Condicion: Factor w/ 2 levels "Sano","Tumor": 2 2 2 2 2 2 2 2 2 1 ...
# $ TE       : Factor w/ 3 levels "Alu","L1","SVA": 1 2 3 1 2 3 1 2 3 1 ...
# $ RPM      : num  0.628 0.182 0.226 1.234 0.287 ...

##graficamos

# Calculamos los valores de RPM en porcentaje
# c3 <- c2 %>%
#   group_by(Tejido, Cromosoma) %>%
#   mutate(Porcentaje = RPM / sum(RPM, na.rm = TRUE))

# c3$Cromosoma = as.factor(c3$Cromosoma)
# c3$Condicion = as.factor(c3$Condicion)

# bar1 <- ggplot(c2) 
# pd <- position_dodge(.9)  

# c2_group <- c2 %>%
#   group_by(Cromosoma, Tejido, Condicion) %>% summarize(avg_RPM = mean(RPM))
# c2$Tejido_Condicion <- paste(c2$Tejido, c2$Condicion, sep = "_")
# c2$Tejido_Condicion <- factor(c2$Tejido_Condicion, levels = unique(c2$Tejido_Condicion))  
# c2_group <- c2_group %>% 
  #   group_by(Tejido) %>% 
  #   mutate(Porcentaje = avg_RPM / sum(avg_RPM, na.rm = TRUE) * 100)
  # pivot_wider(names_from = c("Tejido","Condicion"), values_from = "Porcentaje")

colores <- c("Trophoblast" = "#0DF305", 
             "Breast Tumor" = "#008CCF","Healthy Breast" = "#D52B30",
             "Colon Tumor" = "#F78472","Healthy Colon" = "#A200FA",
             "Liver Tumor" = "#BDB631","Healthy Liver" = "#0B2B9C")
  
  c3 <- c2 %>% group_by(Cromosoma, Tejido_Condicion, TE) %>% 
  summarize(avg_RPM = mean(RPM)) %>% 
  mutate(Porcentaje = avg_RPM / sum(avg_RPM, na.rm = TRUE) * 100)
bar1 <- ggplot(c3, aes(x = Cromosoma, y = Porcentaje,
             group = Tejido_Condicion,
             color = Tejido_Condicion),
             size = 5) +
  geom_line() +
  geom_point(alpha = 0.4, size = 4) +
  facet_grid(TE ~.) +
  labs(x = "Chromosomes", y = "Proportion") +
  theme_bw()+
  theme(text = element_text(family = "Arial"),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 24, angle = 45, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 30, vjust = 1.5),
        axis.title.x = element_text(size = 30, vjust = 0.1, hjust = 0.5),
        strip.text = element_text(size = 28, color = "black"),
        strip.background = element_blank(),
        legend.key.size = unit(3, "cm"),
        legend.text = element_text(size = 28))+
  guides(color = guide_legend(title = NULL))+
  scale_y_continuous(trans = "log2", 
                     breaks = c(2, 5, 12, 35, 100),
                     labels = scales::percent_format(scale = 1)) +
  scale_color_manual(values = colores)
bar1

bar1 <- bar1 +
  theme(legend.position="bottom",   # Posición de la leyenda
        legend.box="horizontal",    # Orientación horizontal
        legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))   # Ajuste de márgenes

##Guardamos en el directorio de trabajo establecido
ggsave(filename = "stackedbar_chrs_Articulo.jpg",
       plot = bar1,
       width = 45, height = 35, units = "cm",
       dpi = 300)
