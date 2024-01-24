######################################################
##    Prueba estadistica   ###########################
######################################################

#Paquetes requeridos
library(pacman) #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("ggplot2", #la funcion p_load viene del paquete "pacman"
       "ggsci","tidyverse",
       "ggthemes")

##directorio de trabajo
setwd("~/Documents/Scripts_articulo/ScriptsR/")

####  CARGAMOS AL OBJETO "TEs_data.rc" EN NUESTRO AMBIENTE   ########### 

# Unimos Tejido y Condicion

c0 <- TEs_data.rc
c0$Tejido_Condicion <- paste(c0$Tejido, c0$Condicion, sep = "_")
c0$Tejido_Condicion <- factor(c0$Tejido_Condicion, levels = unique(c0$Tejido_Condicion))

## Sumatoria de TEs en cada tejido y condicion diferentes
## Agregamos la columna RC o "conteo de lecturas"
c1 <- aggregate(c0$TE, by = list(c0$TE,
                                 c0$Tejido_Condicion,
                                 c0$RC), FUN = length)

colnames(c1) <- c("TE", "Tejido_Condicion", "RC", "Count_TE")

##Debemos obtener como resultado un data frame con la siguiente
##estructura y encabezados:
head(c1)

#      TE Tejido      RC Count_TE
# 1   ALU  Liver 1549108       74
# 2 LINE1  Liver 1549108       11
# 3   SVA  Liver 1549108        1
# 4   ALU  Liver 2056862      105
# 5 LINE1  Liver 2056862        8
# 6   ALU  Liver 2274425      157


# Calculamos las RPM en funcion del total de TEs en cada tejido
c1 <- c1 %>%
  mutate(RPM = (Count_TE/RC) * 10^6)

## shapiro.test
shapiro.test(c1$RPM)
##Resultado:
# Shapiro-Wilk normality test
# 
# data:  c1$RPM
# W = 0.65681, p-value < 2.2e-16

# el valor de p es menor a 0.05 usado comunmente, por lo que se rechaza
# la hipotesis nula que dice que los datos siguen una distribucion normal

p_load("fitdistrplus")
p_load("logspline")
p_load("dunn.test")

# buscamos que tipo de distribucion es mas parecida 
# a la de nuestros datos 

descdist(c1$RPM, discrete = FALSE)

fit.norm <- fitdist(c1$RPM, "norm")
plot(fit.norm)

fit.expo <- fitdist(c1$RPM, "exp")
plot(fit.expo)

fit.gamma <- fitdist(c1$RPM, "gamma")
plot(fit.gamma)

fit.log <- fitdist(c1$RPM, "lnorm")
plot(fit.log)


# parece que la distribucion de los datos que mejor estima cada plot 
# es la exponencial o la gamma 

# Prueba de kruskal-wallis para datos de grupos independientes sin
# distribucion normal

resultado_kruskal <- kruskal.test(RPM ~ Tejido_Condicion, data = c1)

print(resultado_kruskal)

# segun la prueba de kruskal-wallis la hipotesis nula es que no existen
# diferencias significativas entre los grupos a un valor p de 0.05,
# por lo que se rechaza la hipotesis nula con el resultado de la prueba,
# donde se obtuvo un valor de p menor a 0.05, concluyendo que
# si existen diferencias significativas entre el valor de RPM entre 
# al menos dos tejidos y condiciones diferentes

# dunn.test para identificar que posibles pares de grupos 
# difieren entre si 
## El metodo bonferroni elimina el error de tipo I o de resultados al azar
## diviendo el valor de significancia 0.05 entre el total de comparaciones
## realizadas
resultados.comparaciones <- dunn.test(c1$RPM, c1$Tejido_Condicion, method = "bonferroni")

## Se observa que los grupos que difieren entre si con significancia son 
## Troboflastos----------Colon sano
## Trofoblastos----------Higado sano

p_load("ggsignif")
p_load("ggpubr")
p <- ggplot(c1, aes(x = Tejido_Condicion, y = RPM)) +
  geom_boxplot(fill = as.factor(Tejido_Condicion)) +
  labs(title = "Boxplot de RPM por Tejido Condición") +
  xlab("Tejido Condición") +
  ylab("RPM") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p
# p <- ggboxplot(c1, x = "Tejido_Condicion", y = "RPM", fill = "lightblue") +
#   stat_compare_means(comparisons = list(c("Trofoblastos_Sano", "Breast_Tumor"),
#                                         c("Trofoblastos_Sano", "Breast_Sano"),
#                                         c("Trofoblastos_Sano", "Liver_Tumor"),
#                                         c("Trofoblastos_Sano", "Liver_Sano"),
#                                         c("Trofoblastos_Sano", "Colon_Tumor"),
#                                         c("Trofoblastos_Sano", "Colon_Sano")),
#                      method.args = list("two.sided","bonferroni")) +
#   labs(title = "Boxplot de RPM por Tejido Condición") +
#   xlab("Tejido Condición") +
#   ylab("RPM")

p
