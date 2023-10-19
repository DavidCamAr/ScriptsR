#########################################################################################
# David Campos Arellano                                                                 #
# 19/07/23                                                                              #
# Objetivo: Observar la orientacion de las inserciones de TEs que caen exones e intrones#
# para todos los tejidos en un lollipop plot                                            #
#########################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "tidyverse",
       "ggsci")

#directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

##Filtro para quedarnos solo con cáncer y trofoblastos
g0 <- subset(TEs_data.rc,
             !(Tejido %in% c("Breast", "Colon", "Liver") & Condicion == "Sano"))

# Arreglamos texto
g0$Info_gen <- gsub(":.*","",g0$Info_gen)

##Extraemos unicamente las inserciones en exones e intrones
g0 <- g0 %>% filter(Info_gen %in% c("intron","exon"))
##Eliminamos valores nulos
g0 <-g0 %>% filter(!Orientacion == ".")

##Sumatoria de inserciones de TEs en cada orientacion
g1 <- aggregate(g0$TE, by = list(g0$TE,
                                          g0$Tejido,
                                          g0$Orientacion,
                                          g0$RC,
                                          g0$Info_gen),
                                          FUN = length)
colnames(g1) <- c("TE","Tejido","Orientacion","RC","Ubicacion","Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(g1)

#    TE       Tejido Orientacion      RC Ubicacion Count_TE
# 1 ALU        Liver           - 1549108      exon        1
# 2 ALU        Liver           - 3283222      exon        1
# 3 ALU Trofoblastos           - 3487007      exon        1
# 4 ALU        Liver           - 3708838      exon        1
# 5 ALU        Liver           + 3708838      exon        1
# 6 ALU Trofoblastos           - 3849528      exon        1

# Calculamos las RPM en funcion del total de TEs en cada orientacion
g1 <- g1 %>% 
  mutate(RPM = (Count_TE/RC) * 10^6)

##Promedio RPM
g2 <- aggregate(RPM ~ TE + Tejido + Orientacion + Ubicacion, data = g1, FUN = mean)

# Ajustamos Levels
g2$Orientacion <- recode_factor(g2$Orientacion,
                                "-" = "Antisentido",
                                "+" = "Sentido")

g2$Ubicacion <- recode_factor(g2$Ubicacion,
                                "intron" = "Inserciones en intrones",
                                "exon" = "Inserciones en exones")

g2$Tejido <- recode_factor(g2$Tejido,
                           "Liver" = "Hígado",
                           "Colon" = "Colon",
                           "Breast" = "Mama",
                           "Trofoblastos" = "Trofoblastos")
g2$TE <- recode_factor(g2$TE,
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(g2)

# 'data.frame':	34 obs. of  5 variables:
#   $ TE         : Factor w/ 3 levels "Alu","L1","SVA": 1 1 2 1 1 1 1 2 1 2 ...
# $ Tejido     : Factor w/ 4 levels "Hígado","Colon",..: 3 2 2 1 4 3 2 2 1 1 ...
# $ Orientacion: Factor w/ 2 levels "Antisentido",..: 1 1 1 1 1 2 2 2 2 2 ...
# $ Ubicacion  : Factor w/ 2 levels "Inserciones en intrones",..: 2 2 2 2 2 2 2 2 2 2 ...
# $ RPM        : num  0.0768 0.1278 0.0402 0.1937 0.1524 ...


# Grafica
lollipop <- ggplot(g2, aes(x=Tejido, y=RPM, color = Orientacion)) +
  geom_linerange(aes(x = Tejido, ymin = 0, ymax = RPM, color = Orientacion), ##Se usa linearange para separar las lineas
                 position = position_dodge(0.5), show.legend = FALSE)+       ##de la lollipop, a diferencia de  ##segment, donde junta ambas lineas en una misma
  geom_point(aes (x= Tejido, y=RPM, color = Orientacion), size=4,            ##Ademas de que los RPM se deben pivotear
             position = position_dodge(0.5))+                                ## para cada una de las orientaciones
  coord_flip() +
  theme_linedraw() +
  facet_grid(TE ~ Ubicacion, scales = "free_x") +
  xlab("") +
  ylab("Abundancia (RPM)")+
  theme(
    text = element_text(family = "Arial"),
    legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.y = element_text(size = 20, color = "black", angle = 0),
        strip.text.x = element_text(size = 20, color = "black"),
        strip.background.y = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 20, vjust = 0.1),
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#2E8B57", "#BF3EFF"))+
  scale_y_continuous(trans = 'log2',
                     n.breaks = 6)
lollipop

#Guardamos
ggsave(filename = "lollipop_orientacion.jpg",
       plot = lollipop,
       width = 30, height = 18, units = "cm",
       dpi = 300)
