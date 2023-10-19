#################################################################################
# David Eduardo Campos Arellano                                                 #
# 10/07/23                                                                      #
# Objetivo: Representar la proporcion y abundancia de las inserciones de los    #
# 3 TEs encontrados en los diferentes tejidos en graficas de barras             #
#################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "tidyverse",
       "ggsci",
       "ggpubr",
       "ggrepel")

#directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 


##Contamos cuantos TEs hay por tejido, por archivo, por condicion y por tipo de TE
d1 <- aggregate(TEs_data.rc$TE, by = list(TEs_data.rc$TE,
                                                        TEs_data.rc$Condicion,
                                                        TEs_data.rc$Tejido,
                                                        TEs_data.rc$Archivo,
                                                        TEs_data.rc$RC),
                FUN = length)

# Etiquetamos columnas
colnames(d1) <- c("TE","Condicion","Tejido","Archivo","RC","Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(d1)

#      TE Condicion Tejido   Archivo      RC Count_TE
# 1   ALU     Tumor  Liver ERR232277 1549108       74
# 2 LINE1     Tumor  Liver ERR232277 1549108       11
# 3   SVA     Tumor  Liver ERR232277 1549108        1
# 4   ALU     Tumor  Liver ERR232275 2056862      105
# 5 LINE1     Tumor  Liver ERR232275 2056862        8
# 6   ALU      Sano  Liver ERR232270 2250969      129

##Filtro para quedarnos solo con cáncer y trofoblastos
d1 <- subset(d1,
             !(Tejido %in% c("Breast", "Colon", "Liver") & Condicion == "Sano"))

# Calculamos las RPM en funcion del total de TEs en cada tejido
d1 <- d1 %>%
  mutate(RPM = (Count_TE / RC) * 10^6)

##Promedio de RPM
d2 <- aggregate(d1$RPM, by = list(d1$TE, d1$Tejido, d1$Condicion), FUN = mean)
# Etiquetamos columnas
colnames(d2) <- c("TE","Tejido","Condicion","RPM.prom")

##Calculo de SD de RPM
d3 <- aggregate(d1$RPM, by = list(d1$TE, d1$Tejido, d1$Condicion), FUN = sd)
colnames(d3) <- c("TE","Tejido","Condicion","sd")

##Unimos promedio RPM con SD
d4 <- merge(d2, d3, by = c("TE", "Tejido", "Condicion"))

##Levels
d4$Tejido <- recode_factor(d4$Tejido,
                           "Trofoblastos" = "Trofoblastos",
                           "Breast" = "Mama",
                           "Colon" = "Colon",
                           "Liver" = "Hígado")
d4$TE <- recode_factor(d4$TE,
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(d4)

# 'data.frame':	12 obs. of  5 variables:
#   $ TE       : Factor w/ 3 levels "Alu","L1","SVA": 1 1 1 1 2 2 2 2 3 3 ...
# $ Tejido   : Factor w/ 4 levels "Trofoblastos",..: 2 3 4 1 2 3 4 1 2 3 ...
# $ Condicion: chr  "Tumor" "Tumor" "Tumor" "Sano" ...
# $ RPM.prom : num  7.09 18.46 31.06 3.09 2.53 ...
# $ sd       : num  6.01 4.44 22.9 1.49 3.09 ...

# Grafico
bar0 <- ggplot(d4, aes(x = Tejido, y = RPM.prom, fill = TE))+
  geom_bar(stat = "identity",
                     position = position_dodge()) +
  geom_errorbar(aes(ymin = RPM.prom - sd, ymax = RPM.prom + sd), width = 0.2,
                position = position_dodge(0.9), color = "black")+
  labs(
  y ="Abundancia (RPM)")+
  theme_linedraw()+
  theme(
    axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5, vjust = 1.5),
        legend.title = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size =20),
    panel.border = element_blank())+
  scale_y_continuous(n.breaks = 10)+
  scale_fill_d3()
bar0
# 

############# PROPORCION #################

##Calculamos el porcentaje para las etiquetas en cada proporcion
d4 <- d4 %>%
  group_by(Tejido) %>%
  mutate(Porcentaje = RPM.prom / sum(RPM.prom))

bar <- ggplot(d4, aes(x = Tejido, y = RPM.prom, fill =TE))+
  geom_bar(stat = "identity",
           position = "fill")+
  geom_label(show.legend = FALSE,  ##Para eliminar las etiquetas "a" de la legenda de colores
    aes(label = scales::percent(Porcentaje)),
    position = position_fill(vjust = 0.5),  # Ajusta la posición de las etiquetas
    size = 5,  
    color = "white",
    fontface = "bold"  # Estilo de fuente
  ) +
  labs(
    y ="Proporción")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size =20))+
  scale_y_continuous(labels = scales::percent_format(scale = 100, suffix = "%",
                                                     accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1))+
  scale_fill_d3()
bar

##Unimos plots
combined_plot <- ggarrange(
  bar0, bar,
  ncol = 2,
  common.legend = TRUE, # Combina las leyendas en una sola
  widths = c(1, 1),  ##modifica el ancho de los gráficos
  hjust = -0.2,
  labels = c("A", "B"),
  font.label = list(size = 20, family = "Arial"),
  align = "h", #alinea la altura de los graficos
  legend = "bottom") 
combined_plot 

##Guardamos en el directorio de trabajo establecido
ggsave(filename = "prop.tejido.TEs.jpg",
       plot = combined_plot,
       width = 30, height = 18, units = "cm",
       dpi = 300)