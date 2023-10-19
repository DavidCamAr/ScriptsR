####################################################################################
# David Campos Arellano                                                            #
# 17/06/23                                                                         #
# Objetivo: Analizar cuantas inserciones caen en la referencia                     #
# por medio de stacked percentage bars, para observar la proporcion de estas        # 
####################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "tidyverse",
       "ggsci",
       "scales")

#directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

##Arreglamos texto
TEs_data.rc$Copia_referencia <- gsub('[0-9]+','', TEs_data.rc$Copia_referencia)   ## en esta linea recortamos todos los numeros de la columna
TEs_data.rc$Copia_referencia <- gsub('_\\.$','', TEs_data.rc$Copia_referencia) ## las lineas \\ escapan al . que literal esta en el texto (para que no sea tomado
TEs_data.rc$Copia_referencia <- gsub("SVA_","", TEs_data.rc$Copia_referencia)                 ##como cualquier caracter), el $ indica que el patron que se busca                                 
TEs_data.rc$Copia_referencia <- gsub("LINE_","", TEs_data.rc$Copia_referencia)                ##esta al final de la cadena de texto
TEs_data.rc$Copia_referencia <- gsub("Alu_","", TEs_data.rc$Copia_referencia)
TEs_data.rc$Copia_referencia <- gsub("Left_is_copy","Fall_in_copy", TEs_data.rc$Copia_referencia)
TEs_data.rc$Copia_referencia <- gsub("Right_is_copy","Fall_in_copy", TEs_data.rc$Copia_referencia)


# Contamos cuantos TEs hay con referencia, usando la funcion aggregate(), con el 
# parametro FUN = length. El conteo se hara solo con las variables de interes

##Conteo de TEs con referencia
b1 <- aggregate(TEs_data.rc$TE, by = list(TEs_data.rc$TE,
                                                        TEs_data.rc$Condicion,
                                                        TEs_data.rc$Tejido,
                                                        TEs_data.rc$RC,
                                                        TEs_data.rc$Copia_referencia),
                                                        FUN = length)
# Renombramos columnas
colnames(b1) <- c("TE","Condicion","Tejido","RC",
                  "Copia_referencia","Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(b1)

#      TE Condicion Tejido      RC Copia_referencia Count_TE
# 1   ALU     Tumor  Liver 1549108     Fall_in_copy        9
# 2 LINE1     Tumor  Liver 1549108     Fall_in_copy        2
# 3   ALU     Tumor  Liver 2056862     Fall_in_copy       10
# 4 LINE1     Tumor  Liver 2056862     Fall_in_copy        2
# 5   ALU      Sano  Liver 2250969     Fall_in_copy       15
# 6 LINE1      Sano  Liver 2250969     Fall_in_copy        3

# Calculamos las RPM en funcion del total de TEs con Copia_referencia contados  
b1 <- b1 %>%
  mutate(RPM = (Count_TE / RC) * 10^6)

# Calculamos el promedio de RPM con aggregate() y el parametro FUN = mean
b2 <- aggregate(RPM ~ Copia_referencia + TE + Condicion + Tejido, data = b1, FUN = mean)

# Ajustamos los Levels
b2$TE <- recode_factor(b2$TE, 
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")
b2$Tejido <- recode_factor(b2$Tejido,
                                        "Trofoblastos" = "Trofoblastos",
                                        "Breast" = "Mama",
                                        "Colon" = "Colon",
                                        "Liver" = "Hígado")
b2$Copia_referencia <- recode_factor(b2$Copia_referencia,
                           "Fall_in_copy" = "Con Ref.",
                           "not_in_copy" = "Sin Ref.")


##Hacemos un subset, eliminando datos sanos de aquellos tejidos con datos 
##pareados "Tumor/Sano", quedandonos con un data frame donde solo tendremos 
##los datos del trofoblasto y los tejidos con Tumor.

b3 <- subset(b2, !(Tejido %in% c("Hígado", "Colon", "Mama") & Condicion == "Sano"))


#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(b3)

# 'data.frame':	21 obs. of  5 variables:
#   $ Copia_referencia: Factor w/ 2 levels "Con Ref.","Sin Ref.": 1 2 1 2 2 1 2 1 2 1 ...
# $ TE              : Factor w/ 3 levels "Alu","L1","SVA": 1 1 2 2 3 1 1 2 2 3 ...
# $ Condicion       : chr  "Tumor" "Tumor" "Tumor" "Tumor" ...
# $ Tejido          : Factor w/ 4 levels "Trofoblastos",..: 2 2 2 2 2 3 3 3 3 3 ...
# $ RPM             : num  0.899 6.188 0.599 2.037 0.19 ...

##Grafica

##Modificacion de las etiquetas de la leyenda de colores
etiquetas <- c("Con \nreferencia",
               "Sin \nreferencia")

st.br1 <- ggplot(b3, aes(x = TE,
                         y = RPM,
                         fill = Copia_referencia,
                         shape = Copia_referencia)) +
  geom_bar(color = "black", 
           position = "fill", stat = "identity",
           width = 0.8)+
  facet_grid(~Tejido) + 
       ylab("Proporción") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    text = element_text(family = "Arial"),
        strip.text.x = element_text(size = 23, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 20, vjust = 1.9),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 19),
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 19))+
  scale_y_continuous(labels = scales::percent_format(scale = 100, suffix = "%", #Para poner el eje Y en porcentaje
                                                     accuracy = 1),         ##agregando el sufijo %
                    breaks = c(0, 0.25, 0.5, 0.75, 1))+
  scale_fill_manual(labels = etiquetas, values = c("#FFAEB9","#BCD2EE"))
st.br1  

##Guardamos en el directorio de trabajo establecido
ggsave(filename = "bar_cp.ref.rm.jpg",
       plot = st.br1,
       width = 34, height = 18, units = "cm",
       dpi = 300) 

