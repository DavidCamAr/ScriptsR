####################################################################
# David Campos Arellano                                            #
# 07/08/23                                                         #
# Objetivo: Representar las lecturas y la calidad de cada archivo  #
# WGS analizado                                                    # 
####################################################################
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "ggrepel",
       "tidyverse",
       "ggsci",
       "ggpubr")
       # "lessR")

#directorio de trabajo
setwd("~/Documents/Datos_xTea")

##Pie plot para la proporcion de archivos para cada tejido

#### Data

Archivos <- c(16,17,17,15)
Tejidos <- c("Trofoblasto","Cáncer \nde mama","Cáncer \nde colon","Cáncer \nde hígado")

##Creamos data frame con Archivos y Tejidos
m0 <- data.frame(Tejidos, Archivos)

#Ajustamos levels
m0$Tejidos <- factor(m0$Tejidos, levels = c(
                            "Trofoblasto",
                            "Cáncer \nde mama",
                            "Cáncer \nde colon",
                           "Cáncer \nde hígado"))

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(m0)

# 'data.frame':	4 obs. of  2 variables:
#   $ Tejidos : Factor w/ 4 levels "Trofoblasto",..: 1 2 3 4
#   $ Archivos: num  16 17 17 15

##Calculamos un valor total para etiquetar despues
valor_total <- sum(m0$Archivos)

##Graficamos

### Con ggplot

# hsize <- 2
# mp <- ggplot(m0, aes(x = hsize , y= Archivos, fill=Tejidos))+
#   geom_col(color = "black")+
#   geom_text(aes(label = Archivos),
#             size = 12,
#             position = position_stack(vjust = 0.5))+
#   geom_text(aes(x = hsize / 7.8, y = 0, label = valor_total),
#             color = "black",
#             size = 20) +  # Agregar el texto con el valor total
#   coord_polar(theta = "y") +
#   scale_fill_d3(palette = "category20c") +
#   xlim(c(0.2, hsize + 2))+
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.title = element_blank(),
#         legend.position = "none",
#         legend.key.size = unit(1, "cm"),
#         legend.text = element_text(size = 24))+
#   scale_fill_tron()
# 
# donut_labels <- mp +
#   geom_label_repel(data = m00,
#                    aes(y = pos, label = Tejidos, fill = Tejidos),
#                    size = 7, position = position_nudge_repel(x = 1.8, y = 0))
# donut_labels
##Con PieChart... interesante
# PieChart(Tejido, data = j0,
#          values = "input",
#          main = NULL)


##Con ggdonutchart
hsize <- 1
donut <- ggdonutchart(m0, "Archivos", label = NULL,
  fill = "Tejidos", color = "black", size = 1)+
  scale_fill_tron()+
  theme_void()+
  geom_text(aes(label = Archivos),
            size = 17,
            position = position_stack(vjust = 0.5))+
  geom_text(aes(x = hsize / 5, y = 0, label = valor_total),  # Agregar el texto con el valor total
            color = "black",
            size = 32)  +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none")+
  xlim(c(0.2, hsize + 3)) ##IMPORTANTE delimita las dimensiones del plot para poder usar nudge_x

# establecer posiciones para las etiquetas por fuera
m00 <- m0 %>% 
  mutate(csum = rev(cumsum(rev(Archivos))), 
         pos = (Archivos)/2 + lead(csum, 1),
         pos = if_else(is.na(pos), (Archivos)/2, pos))

# Añadir etiquetas por fuera del donut
donut_labels <- donut +
  geom_label_repel(data = m00,
                   aes(y = pos, label = Tejidos, fill = Tejidos),
                   size = 12, position = position_nudge_repel(x = 1.7, y=0))
donut_labels
#Guardamos
ggsave(filename = "donut.archivos.jpg",
       plot = donut_labels,
       width = 10, height = 10,
       dpi = 300)


##Promedio de reads por tejido
##Leemos datos de conteo de reads despues del control de calidad
rc_af.trim <- read_tsv("~/Documents/Datos_xTea/reads-af.trim_04-08-23.tsv", col_skip())

##Renombrar columnas
colnames(rc_af.trim) <- c("pre_rc","Archivo")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(rc_af.trim)

# # A tibble: 6 × 2
# pre_rc Archivo  
# <dbl> <chr>    
#   1 14256979 CRR101242
# 2 13978968 CRR101243
# 3 13948029 CRR101244
# 4 13942635 CRR101245
# 5 15398110 CRR101250
# 6 15082853 CRR101251

##Dividimos entre 4 para conocer los Read Counts de los archivos fastq
# ya que los archivos fastq contiene 4 lineas de texto y solo
# una de estas pertece a las reads
rc_af.trim <- rc_af.trim %>%
  mutate(RC.af = pre_rc/4)

##Filtramos solo con cancer y trofoblastos
## Reemplazar en cada caso, incluso este paso debe ser ignorado si no hay nada que 
## filtrar
tejidos_sanos <- c("SRR2155313","SRR2155316","SRR13215448","SRR13215449","SRR13215450","SRR13215451",
                   "SRR13215452","SRR13215453","SRR13215454","SRR13215351","SRR13215352","SRR13215353",
                   "SRR13215354","SRR13215355","SRR13215356","SRR13215357","SRR13215359","ERR232278",
                   "ERR232276","ERR232274","ERR232272","ERR232270","ERR1916595","ERR1916597",
                   "ERR1916599","ERR1916601","ERR1916603","ERR1916605","ERR1916607","ERR1916609",
                   "ERR1916611","ERR1916613","SRR22797220","CRR033294","CRR033295","SRR12017692",
                   "SRR12017694","SRR12017696","SRR12017698","SRR12017700","SRR22797221",
                   "SRR22797223","SRR7707726","SRR12017702","SRR12017706")

j0 <- subset(rc_af.trim, !(Archivo %in% tejidos_sanos))

#Asignamos etiquetas de tejidos

tejidos <- c("Trofoblastos","Mama","Colon","Hígado")
tamaño_muestras <- c(16,17,17,15)
 
##Repetimos para cada observacion
repet_tejidos <- rep(tejidos, times = tamaño_muestras)             

##Y agregamos columna como vector
j0$Tejido <- repet_tejidos

##Calculamos el promedio de RC por tejido
j1 <- aggregate(RC.af ~ Tejido, data = j0, FUN = mean)

# Ajustamos Levels
j1$Tejido <- factor(j1$Tejido, levels = c(
  "Hígado",
  "Colon",
  "Mama",
  "Trofoblastos"))
                           
#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(j1)

# 'data.frame':	4 obs. of  2 variables:
#   $ Tejido: Factor w/ 4 levels "Hígado","Colon",..: 2 1 3 4
#   $ RC.af : num  22521903 11545857 12487005 13484496

# Grafica
jm <- ggplot(j1, aes(x = Tejido, y = RC.af))+
  geom_segment(aes(x = Tejido,
                   xend = Tejido,
                   y = 0, yend = RC.af),
               color = "grey", lwd =2)+
  geom_point(size= 6, pch= 21, bg= "#458B00", color = "black")+
  xlab("")+
  ylab("Reads (x10⁶)")+
  coord_flip()+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 17, vjust = 0.1),
        axis.text = element_text(size = 18, color= "black"))+
        #axis.ticks.y = element_blank())+
  scale_y_continuous(n.breaks = 10,
                     labels = function(x) paste0(format(x / 1e6, nsmall = 0.5), ""))
jm

##Grafico de reads antes y despues de pasar los archivos fastq por el control de calidad
##Leemos los datos del conteo de reads antes del control de calidad
rc_bf.trim <- read_tsv("~/Documents/Datos_xTea/reads-bf.trim.tsv", col_skip())

##renombrar columnas
colnames(rc_bf.trim) <- c("pre_rc","Archivo")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(rc_bf.trim)

# # A tibble: 6 × 2
# pre_rc Archivo    
# <dbl> <chr>      
#   1 27886814 SRR12017692
# 2 27602307 SRR12017693
# 3 27751221 SRR12017694
# 4 27824958 SRR12017695
# 5 27542971 SRR12017696
# 6 21348916 SRR12017697

##Dividimos entre 4 para conocer los Read Counts de los archivos fastq
rc_bf.trim <- rc_bf.trim %>%
  mutate(RC.bf = pre_rc/4)

##eliminamos columnas del pre_rc
rc_bf.trim <- rc_bf.trim[ ,-1]
rc_af.trim <- rc_af.trim[ ,-1]

##Unimos por Archivo tanto los datos antes como despues del control de calidad
RC_tt <- merge(rc_bf.trim, rc_af.trim, by = "Archivo")

##Filtramos tejido sano 
RC_tt.fil <- subset(RC_tt, !(Archivo %in% tejidos_sanos))

##Etiquetamos archivos con tejidos
## Reemplazar en cada caso con el nombre de los archivos correspondientes
Trofoblastos <- c("SRR6941195","SRR6941196","SRR6941197","SRR6941198",
                  "SRR6954088","SRR6954089","Embrio1","Embrio2",
                  "CRR101242","CRR101243","CRR101244","CRR101245",
                  "CRR101250","CRR101251","CRR101252","CRR101253")
Mama <- c("SRR22797222","SRR22797224","SRR22797225","SRR7707725","SRR7707726",
          "SRR7707727","SRR7707728","SRR7707729","SRR12017693","SRR12017695",
          "SRR12017697","SRR12017699","SRR12017701","SRR12017703","SRR12017704",
          "SRR12017705","SRR12017707")
Colon <- c("SRR2155320","SRR2155329","SRR13215402","SRR13215403","SRR13215404",
           "SRR13215406","SRR13215407","SRR13215408","SRR13215409","SRR13215410",
           "SRR13215411","SRR13215412","SRR13215413","SRR13215414","SRR13215415",
           "SRR13215417","SRR13215418")
Hígado <- c("ERR232279","ERR232277","ERR232275","ERR232273","ERR232271",
            "ERR1916596","ERR1916598","ERR1916600","ERR1916602","ERR1916604",
            "ERR1916606","ERR1916608","ERR1916610","ERR1916612","ERR1916614")

etiquetas_tejidos <- data.frame(
  Archivo = c(Trofoblastos, Mama, Colon, Hígado),
  Tejido = c(rep("Trofoblastos", length(Trofoblastos)),
             rep("Mama", length(Mama)),
             rep("Colon", length(Colon)),
             rep("Hígado", length(Hígado))))

##Unimos etiquetas
RC_tt.fil <- merge(RC_tt.fil, etiquetas_tejidos, by = "Archivo")

# Ajustamos Levels
RC_tt.fil$Tejido <- factor(RC_tt.fil$Tejido, levels = c(
  "Trofoblastos",
  "Mama",
  "Colon",
  "Hígado"))

#Para graficar tenemos que tener un d.f con la siguiente estructura y encabezados
str(RC_tt.fil)

# 'data.frame':	64 obs. of  4 variables:
#   $ Archivo: chr  "CRR101242" "CRR101243" "CRR101244" "CRR101245" ...
# $ RC.bf  : num  7278062 7276681 7291236 7285401 7820617 ...
# $ RC.af  : num  3564245 3494742 3487007 3485659 3849528 ...
# $ Tejido : Factor w/ 4 levels "Trofoblastos",..: 1 1 1 1 1 1 1 1 1 1 ...

# Grafica
RC.lollipop <- ggplot(RC_tt.fil)+
  geom_segment( aes(x=Archivo, xend=Archivo, y=0, yend=RC.bf), color= "grey")+
  geom_point(aes (x=Archivo, y=RC.bf, color="Antes de \nTrimmomatic"), size=4) +
  geom_segment( aes(x=Archivo, xend=Archivo, y=0, yend=RC.af), color="grey") +
  geom_point( aes(x=Archivo, y=RC.af, color="Después de \nTrimmomatic"), size=2) +
  theme_linedraw()+
  facet_grid(~Tejido, scales = "free_x")+
  xlab("") +
  ylab("Reads (x10⁶)")+
  theme(
    legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 17),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 18, color = "black", vjust = 1),
        axis.text.x = element_text(size = 12, color = "black", angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 17, vjust = 1.5),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank())+
  scale_color_manual(values = c("Antes de \nTrimmomatic" = "#27408B", "Después de \nTrimmomatic" = "#458B00"))+
  scale_y_continuous(labels = function(x) paste0(format(x / 1e6, nsmall = 0.5), ""))

RC.lollipop

##Podemos unir graficos ahora en un mismo panel 

reads_plot <- ggarrange(RC.lollipop, jm,
          labels = c("A", "B"),
          font.label = list(size = 20, family = "Arial"),
          hjust = -0.1,
          vjust = 1.2,
          ncol = 1,
          nrow = 2)  # Ajusta el tamaño de fuente de los labels "A" y "B"

reads_plot

#Guardamos
ggsave(filename = "reads_tejido.jpg",
       plot = reads_plot,
       width = 30, height = 18, units = "cm",
       dpi = 300)

