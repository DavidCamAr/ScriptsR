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

##Filtro para quedarnos solo con cáncer y trofoblastos
f0 <- subset(TEs_data.rc,
             !(Tejido %in% c("Breast", "Colon", "Liver") & Condicion == "Sano"))

# Arreglamos texto

f0$Info_gen <- gsub(":.*","",f0$Info_gen)

##Sumatoria de inserciones de TEs en cada ubicacion genomica
f1 <- aggregate(f0$TE, by = list(f0$TE,
                                       f0$Tejido,
                                       f0$Info_gen,
                                       f0$RC),
                FUN = length)

colnames(f1) <- c("TE",
                  "Tejido",
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
f2 <- aggregate(f1$RPM, by = list(f1$TE, f1$Region_insercion, f1$Tejido), FUN = mean)

# Etiquetamos columnas
colnames(f2) <- c("TE","Region_insercion","Tejido","RPM.prom")

##Levels
f2$Tejido <- recode_factor(f2$Tejido,
                           "Trofoblastos" = "Trofoblastos",
                           "Breast" = "Mama",
                           "Colon" = "Colon",
                           "Liver" = "Hígado")
f2$TE <- recode_factor(f2$TE,
                       "ALU" = "Alu",
                       "LINE1" = "L1",
                       "SVA" = "SVA")
f2$Region_ins <- recode_factor(f2$Region_ins,
                               "exon" = "Exón",
                               "intron" = "Intrón",
                               "not_gene_region" = "Región intergénica",
                               "down_stream" = "Downstream",
                               "up_stream" = "Upstream",
                               "UTR5" = "UTR5",
                               "UTR3" = "UTR3")
#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(f2)

# 'data.frame':	54 obs. of  5 variables:
#   $ TE              : Factor w/ 3 levels "Alu","L1","SVA": 1 2 3 1 1 2 3 1 2 3 ...
# $ Region_insercion: chr  "down_stream" "down_stream" "down_stream" "exon" ...
# $ Tejido          : Factor w/ 4 levels "Trofoblastos",..: 2 2 2 2 2 2 2 2 2 2 ...
# $ RPM.prom        : num  0.2682 0.0806 0.1454 0.1305 3.8922 ...
# $ Region_ins      : Factor w/ 7 levels "Exón","Intrón",..: 4 4 4 1 2 2 2 3 3 3 ...


# Grafica
bar.reg <- ggplot(f2, aes(x = TE, y = RPM.prom, fill = Region_ins))+
             geom_bar(stat = "identity",
                      position = "fill")+
  facet_grid( ~ Tejido)+
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

