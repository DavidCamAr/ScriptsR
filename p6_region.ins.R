#####################################################################################
# David Campos Arellano                                                             #
# 17/07/23                                                                          #
# Objetivo: Visualizar el sitio de insercion de TEs (intron, exon, etc) para todos  #
# los tejidos                                                                       #
#####################################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "tidyverse",
       "ggsci",
       "ggthemes")

#directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE########### 

# Arreglamos texto

TEs_data.rc$Info_gen <- gsub(":.*","",TEs_data.rc$Info_gen)

##Sumatoria de inserciones de TEs en cada ubicacion genomica
f1 <- aggregate(TEs_data.rc$TE, by = list(TEs_data.rc$TE,
                                       TEs_data.rc$Tejido,
                                       TEs_data.rc$Condicion,
                                       TEs_data.rc$Info_gen,
                                       TEs_data.rc$RC),
                FUN = length)

colnames(f1) <- c("TE",
                  "Tejido", "Condicion",
                  "Region_insercion", 
                  "RC","Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(f1)

#      TE Tejido Region_insercion      RC Count_TE
# 1   ALU  Liver      down_stream 1549108        2
# 2   ALU  Liver             exon 1549108        1
# 3   ALU  Liver           intron 1549108       52
# 4 LINE1  Liver           intron 1549108        7
# 5   SVA  Liver           intron 1549108        1
# 6   ALU  Liver  not_gene_region 1549108       19

# Calculamos las RPM en funcion del total de TEs en cada region
f1 <- f1 %>%
  mutate(RPM = (Count_TE / RC) * 10^6)

##Promedio de RPM
# f2 <- aggregate(f1$RPM, by = list(f1$TE, f1$Region_insercion, f1$Tejido), FUN = mean)
f2 <- aggregate(RPM ~ Region_insercion + TE + Condicion + Tejido, data = f1, FUN = mean)

#Sacamos combinaciones entre observaciones
tejidos <- unique(f2$Tejido)
condicion <- unique(f2$Condicion)
region <- unique(f2$Region_insercion)
tipos_TE <- unique(f2$TE)
combinaciones <- expand.grid(Region_insercion = region, Tejido = tejidos,
                             Condicion = condicion, TE = tipos_TE)

#Juntamos combinaciones al d.f de trabajo
f2 <- merge(combinaciones, f2, by = c("Region_insercion","Tejido","Condicion","TE"),
            all.x = TRUE, all.y = FALSE)
##Agregamos valores =NA a los valores faltantes
f2$RPM[is.na(f2$RPM)] <- NA

##Levels
f2$Tejido <- recode_factor(f2$Tejido,
                           "Trofoblastos" = "Trofoblasto",
                           "Breast" = "Mama",
                           "Colon" = "Colon",
                           "Liver" = "Hígado")
f2$TE <- recode_factor(f2$TE,
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")
f2$Region_insercion <- recode_factor(f2$Region_insercion,
                               "exon" = "Exón",
                               "intron" = "Intrón",
                               "not_gene_region" = "Región intergénica",
                               "down_stream" = "Downstream",
                               "up_stream" = "Upstream",
                               "UTR5" = "UTR5",
                               "UTR3" = "UTR3")

##Filtro para quedarnos solo con cáncer y trofoblastos
f3 <- subset(f2, !(Tejido %in% c("Mama", "Colon", "Hígado") & Condicion == "Sano"))

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(f3)

# 'data.frame':	105 obs. of  5 variables:
# $ Region_insercion: Factor w/ 7 levels "Exón","Intrón",..: 4 4 4 4 4 4 4 4 4 4 ...
# $ Tejido          : Factor w/ 4 levels "Trofoblasto",..: 2 2 2 3 3 3 4 4 4 1 ...
# $ Condicion       : Factor w/ 2 levels "Sano","Tumor": 2 2 2 2 2 2 2 2 2 1 ...
# $ TE              : Factor w/ 3 levels "Alu","L1","SVA": 1 2 3 1 2 3 1 2 3 1 ...
# $ RPM             : num  0.2682 0.0806 0.1454 0.3999 0.0636 ...

# Calculamos los valores de RPM en porcentaje
f3 <- f3 %>%
  group_by(Tejido, TE) %>%
  mutate(Porcentaje = RPM / sum(RPM, na.rm = TRUE))

# Grafica
bar.reg <- f3 %>% 
  ggplot(aes(x = TE, y = RPM, fill = Region_insercion))+
  geom_col(position = "fill")+
  geom_label(show.legend = NULL,
            data=subset(f3, Tejido == "Trofoblasto" 
                        & Region_insercion == "Exón" 
                        & TE == "Alu"),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                      accuracy = 1)),
            position = position_fill(vjust = 0.98),  # Ajusta la posición de las etiquetas
            size = 5,
            color = "black",
            fontface = "bold") +
  geom_label(show.legend = NULL,
            data=subset(f3, Tejido %in% c("Mama","Colon","Hígado") 
                        & Region_insercion == "Exón" 
                        & TE %in% c("Alu","L1")),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                      accuracy = 1)),
            position = position_fill(vjust = 0.99),  # Ajusta la posición de las etiquetas
            size = 5,
            color = "black",
            fontface = "bold") +
  geom_label(show.legend = NULL,
            data=subset(f3, Region_insercion %in% c("Intrón","Región intergénica","Downstream") 
                        & Tejido %in% c("Mama","Hígado","Trofoblasto") & TE == "SVA"),
            aes(label=scales::percent(Porcentaje, scale = 100, suffix = "%",
                                      accuracy = 1)),
            position = position_fill(vjust = 0.5),  # Ajusta la posición de las etiquetas
            size = 5,
            color = "black",
            fontface = "bold") +
  facet_grid(~Tejido)+
       ylab("Proporción")+
  theme_bw()+
  theme(
    text = element_text(family = "Arial"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =17))+
  scale_fill_few()+
  scale_y_continuous(labels = scales::percent_format(scale = 100, suffix = "%",
                                                     accuracy = 1),
                     breaks = c(0.1, 0.2, 0.3, 0.4,
                                0.5, 0.6, 0.7, 0.8, 0.9, 1))
bar.reg

##guardamos
ggsave(filename = "Region.TEs.jpg",
       plot = bar.reg,
       width = 30, height = 18, units = "cm",
       dpi = 300)

