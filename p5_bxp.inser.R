#####################################################################################
# David Campos Arellano                                                             #
# 17/05/23                                                                          #
# Objetivo: Realizar un boxplot de la distribucion de las inserciones de TEs        #
# en los archivos analizados para los diferentes tejidos                            #
#####################################################################################

#Paquetes requeridos
library(pacman) #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "ggsci","tidyverse",
       "ggthemes")

#directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

##Filtro para quedarnos solo con cáncer y trofoblastos
e0 <- subset(TEs_data.rc,
             !(Tejido %in% c("Breast", "Colon", "Liver") & Condicion == "Sano"))
##Sumatoria de TEs por tejido
e1 <- aggregate(e0$TE, by = list(e0$TE,
                                 e0$Tejido,
                                 e0$RC), FUN = length)

colnames(e1) <- c("TE", "Tejido", "RC", "Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(e1)

#      TE Tejido      RC Count_TE
# 1   ALU  Liver 1549108       74
# 2 LINE1  Liver 1549108       11
# 3   SVA  Liver 1549108        1
# 4   ALU  Liver 2056862      105
# 5 LINE1  Liver 2056862        8
# 6   ALU  Liver 2274425      157


# Calculamos las RPM en funcion del total de TEs en cada tejido
e1 <- e1 %>%
  mutate(RPM = (Count_TE/RC) * 10^6)

#Levels
e1$Tejido <- recode_factor(e1$Tejido,
                           "Trofoblastos" = "Trofoblastos",
                           "Breast" = "Cáncer de mama",
                           "Colon" = "Cáncer de colon",
                           "Liver" = "Cáncer de hígado")
e1$TE <- recode_factor(e1$TE,
                           "ALU" = "Alu",
                           "LINE1" = "L1",
                           "SVA" = "SVA")

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(e1)

# 'data.frame':	165 obs. of  5 variables:
#   $ TE      : Factor w/ 3 levels "Alu","L1","SVA": 1 2 3 1 2 1 2 3 1 2 ...
# $ Tejido  : Factor w/ 4 levels "Trofoblastos",..: 4 4 4 4 4 4 4 4 4 4 ...
# $ RC      : num  1549108 1549108 1549108 2056862 2056862 ...
# $ Count_TE: int  74 11 1 105 8 157 14 1 227 14 ...
# $ RPM     : num  47.769 7.101 0.646 51.049 3.889 ...

# Grafico
bxp0 <- ggplot(e1, 
               aes(x=TE, 
                   y=RPM, 
                   color= TE)) +
  geom_point(size = -1)+ ##usamos estas lineas para cambiar la forma
  geom_boxplot(show.legend = FALSE) + ##de las leyendas que por default nos da boxplot
  guides(color = guide_legend(keywidth = unit(1, "cm"),
                              keyheight = unit(1, "cm"),
                              override.aes = list(size = 8,
                                                  shape = 20)))+
  facet_wrap(~Tejido) +
  geom_jitter(
    position = position_jitter(width = 0.1),
    size=1, 
    alpha=0.5, 
    color = "black")+
  labs(
       y ="Abundancia (RPM)")+
  theme_linedraw()+
  theme(
    text = element_text(family = "Arial"),
    strip.background = element_blank(),
        strip.text.x = element_text(size = 20, color = "black"),
    panel.border = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 20, hjust = 0.5),
        axis.text.y = element_text(size = 15, color = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(size =20))+ 
  scale_y_continuous(trans = 'log2', 
                           n.breaks = 7, 
                           labels = function(RPM) sprintf("%.1f",   ##Para redondear los numeros a solo un decimal
                                                         round(RPM, 1))) + 
  scale_colour_few()
bxp0

#Guardamos
ggsave(filename = "boxplotTEs_Tejido.jpg",
       plot = bxp0,
       width = 30, height = 18, units = "cm",
       dpi = 300)
    
