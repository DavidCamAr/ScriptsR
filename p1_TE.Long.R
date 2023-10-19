##############################################################################
# David Campos Arellano                                                      #
# 07/07/23                                                                   #
# Objetivo: Graficar la longitud de los TEs y de las terminaciones TSD       #
# correspondientes a cada familia, en los diferentes tejidos analizados,     #
# por medio de un poligono de frecuencias                                    #
##############################################################################

#Paquetes requeridos
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2","tidyverse", #la funcion p_load viene del paquete "pacman"
       "ggsci",
       "forcats",
       "ggpubr")

##Directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 



 ######## LONGITUD TE #################

## Seleccionamos unicamente las variables con las que nos interesa trabajar, para
## reducir errores 
a1 <- TEs_data.rc %>% 
  select(TE,Archivo,Tejido,Condicion,Longitud_TE)

##Hacemos un subset, eliminando datos sanos de aquellos tejidos con datos 
##pareados "Tumor/Sano", quedandonos con un data frame donde solo tendremos 
##los datos del trofoblasto y los tejidos con Tumor.

a1 <- subset(a1, !(Tejido %in% c("Liver", "Colon", "Breast") & Condicion == "Sano"))


##Acomodamos el orden en el que queremos visualizar nuestras variables o "Levels"
##con recode_factor
a1$Tejido <- recode_factor(a1$Tejido,
                           "Trofoblastos" = "Trofoblastos",
                           "Breast" = "Mama",
                           "Colon" = "Colon",
                           "Liver" = "Hígado")
a1$TE <- recode_factor(a1$TE,
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")

# Para graficar, tenemos que tener un d.f con la siguiente estructura y encabezados:
str(a1)

# 'data.frame':	16419 obs. of  5 variables:
#   $ TE         : Factor w/ 3 levels "Alu","L1","SVA": 1 1 1 1 1 2 2 1 1 1 ...
# $ Archivo    : chr  "CRR101242" "CRR101242" "CRR101242" "CRR101242" ...
# $ Tejido     : Factor w/ 4 levels "Trofoblastos",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Condicion  : chr  "Sano" "Sano" "Sano" "Sano" ...
# $ Longitud_TE: num  265 277 151 264 259 ...


######## LONGITUD TSD #################

## Seleccionamos unicamente las variables con las que nos interesa trabajar, para
## reducir errores 
a2 <- TEs_data.rc %>% 
  select(TE,Archivo,Tejido,Condicion,Longitud_TSD)

##Hacemos un subset, eliminando datos sanos de los aquellos tejidos con datos 
##pareados "Tumor/Sano", quedandonos con un data frame donde solo tendremos 
##los datos del trofoblasto y los tejidos con Tumor.
a2 <- subset(a2, !(Tejido %in% c("Liver", "Colon", "Breast") & Condicion == "Sano"))

##Ordenamos los "Levels"
a2$Tejido <- recode_factor(a2$Tejido,
                           "Trofoblastos" = "Trofoblastos",
                           "Breast" = "Mama",
                           "Colon" = "Colon",
                           "Liver" = "Hígado")
a2$TE <- recode_factor(a2$TE,
                           "ALU" = "Alu",
                           "LINE1" = "L1",
                           "SVA" = "SVA")

# Para graficar, tenemos que tener un d.f con la siguiente estructura y encabezados:
str(a2)

# 'data.frame':	16419 obs. of  5 variables:
#   $ TE          : Factor w/ 3 levels "Alu","L1","SVA": 1 1 1 1 1 2 2 1 1 1 ...
# $ Archivo     : chr  "CRR101242" "CRR101242" "CRR101242" "CRR101242" ...
# $ Tejido      : Factor w/ 4 levels "Trofoblastos",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Condicion   : chr  "Sano" "Sano" "Sano" "Sano" ...
# $ Longitud_TSD: num  16 15 18 14 10 14 17 14 18 5 ...

##Graficamos la longitud de TEs con ggplot()
freqpol_LTE <- ggplot(a1, aes(x = Longitud_TE, color = TE)) + 
  geom_freqpoly(binwidth=0.25)+  ##primera geometria
  geom_point(stat="bin", aes(y=..count..), binwidth = 0.25, size = 1) +  ##segunda geometria
  xlab("Longitud del TE (pb)")+ 
  ylab("Frecuencia")+
  theme_linedraw()+
  theme(panel.border = element_blank(), ##con la funcion theme podemos ajustar toda la 
        strip.background = element_blank(),  ##visualizacion de ejes, titulos o legendas
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.y = element_text(size = 12, vjust = 2),    
        axis.title.x = element_text(size = 12, vjust = 0.1),
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom",
        legend.text = element_text(size =12))+
  labs(color = NULL)+
  facet_grid(Tejido~.)+
  scale_x_continuous(trans = 'log2',   ##Modificamos la escala del eje x
                     breaks =          ##a log2 (se modifica el eje, no los valores)
                       c(80,130,200,300,400,600,
                                         1000,2000,3000,4000,7000))+
  scale_y_continuous(trans = 'log2',   ##Modificamos la escala del eje y 
                     n.breaks = 6)+    ##a log2 (se modifica el eje, no los valores)
  scale_color_d3()
freqpol_LTE

##Graficamos longitud de TSDs

freqpol_LTSD <- ggplot(a2, aes(x = Longitud_TSD, color = TE)) +
  geom_freqpoly(binwidth=0.25)+  ##primera geometria
  geom_point(stat="bin", aes(y=..count..), binwidth = 0.25, size =1) + ##segunda geometria
  xlab("Longitud TSD (pb)")+ 
  ylab("Frecuencia")+
  theme_linedraw()+
  theme(panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 9, color = "black", angle = 270),  
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 12, vjust = 2),    
        axis.title.x = element_text(size = 12, vjust = 0.1),
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom",
        legend.text = element_text(size =12))+
  labs(color = NULL)+
  facet_grid(Tejido~.)+
  scale_x_continuous(trans = 'log2',    ##Modificamos la escala del eje x
                     n.breaks = 8)+     ##a log2 (se modifica el eje, no los valores)
  scale_y_continuous(trans = 'log2',   ##Modificamos la escala del eje y
                     breaks = c(1,4,16,64))+  ##a log2 (se modifica el eje, no los valores)
  scale_color_d3()
freqpol_LTSD 

################################################################################################3
##Unimos ambos graficos

combined_plot <- ggarrange(
  freqpol_LTE,
  freqpol_LTSD, 
  nrow = 1,
  common.legend = TRUE, # Combina las leyendas en una sola
  labels = c("A", "B"),
  font.label = list(size = 10, family = "Arial"),
  label.x = -0.01, #posicion de las etiquetas en el eje x
  align = c("hv"), #alinea la altura de los graficos
  legend = "bottom")
combined_plot 

##Guardamos en el directorio de trabajo establecido
ggsave(filename = "LongitudTEs.jpg",
       plot = combined_plot,
       width = 23, height = 12, units = "cm",
       dpi = 300)
