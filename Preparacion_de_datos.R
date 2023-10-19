###################################################################################
# David Eduardo Campos Arellano                                                   #
# 05/06/23                                                                        #
#                                                                                 #
# Objetivo: Este script juntara todos los archivos arrojados por xTea en uno solo,#
# se le data formato y se dejara listo para poder usar al graficar                #
###################################################################################

#Paquetes requeridos
library("pacman") #el paquete "pacman" ayudara a gestionar otros paquetes
p_load("readr", "dplyr") #la funcion p_load viene del paquete "pacman"

# Dentro del directorio Documents creamos la carpeta Datos_xTea, de ser
# necesario, crear el directorio completo ~/Documents/Datos_xTea/

# Establecemos el directorio de trabajo
setwd("~/Documents/Datos_xTea/")

##Leemos los datos externos con terminacion .tsv
##Se debe editar el codigo en cada caso con la ruta hacia el archivo externo

##Como ejemplo se usaran los datos de Colon

## bst = Breast
## lvr = Liver
## tfb = Trofoblasto
## cln = Colon

#Ejemplos:
# bst.tsv <- read_tsv("/home/david/Documents/files/Breast/xTea_bst.tsv", 
#                          col_skip())
# bst.tsv <- read_tsv("~/Downloads/analisisTEsfinal/xTea_bst2.tsv",
#                     col_skip())
# lvr.tsv <- read_tsv("/home/david/Documents/files/Liver/xTea_lvr2.tsv",
#                     col_skip())
# tfb.tsv <- read_tsv("/home/david/Documents/files/Placenta/xTea_pcta.tsv",
#                      col_skip())
# tfb.tsv <- read_tsv("~/Downloads/analisisTEsfinal/xTea_tfbs.tsv",
#                      col_skip())
cln.tsv <- read_tsv("~/Documents/Datos_xTea/files/colon/xTea_cln2.tsv",
                    col_skip())

##En cada caso tenemos que obtener un data frame con la siguiente estructura y encabezados:
#(Los datos van a variar para cada set)
head(cln.tsv)

# # A tibble: 6 × 15
# X1          X2 X3      X4    X5    X6    X7    X8    X9    X10   X11   X12   X13   X14  
# <chr>    <dbl> <chr>   <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#   1 chr7  95149737 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    1/1  
# 2 chr12 20755736 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    1/1  
# 3 chr12 66617790 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    1/1  
# 4 chr16   239175 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    0/1  
# 5 chr17 27268029 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    1/1  
# 6 chrX  88773290 SVTYPE… SVLE… END=… TSD=… TSDL… TD_S… STRA… INS_… REF_… GENE… GT    1/1  
# # ℹ 1 more variable: X15 <chr>

##separamos en sano y tumor para etiquetar
########################################################################################
#  SE DEBEN REMPLAZAR LOS NOMBRES DE LOS ARCHIVOS Y DE LOS OBJETOS EN CADA CASO PARA   # 
#                             CADA TEJIDO Y CADA CONDICION SEA EL CASO  
#   
########################################################################################

##los siguientes valores contienen a los archivos con Tumor correspondientes a colon
valores_filtro <- c("SRR13215448", "SRR13215449", "SRR13215450", "SRR13215451",
                    "SRR13215452", "SRR13215453", "SRR13215454", "SRR13215351",
                    "SRR13215352", "SRR13215353", "SRR13215354", "SRR13215355",
                    "SRR13215356", "SRR13215357", "SRR13215359")

##Separamos Tumor de Sanos para etiquetar en cada caso por Condicion y volver a unir
##posteriormente
cln.tsv_sano <- cln.tsv %>% filter(cln.tsv$X15 %in% valores_filtro )
cln.tsv_sano$Condicion <- "Sano"

##El data frame para tumor se llamara cln.tsv_tumor y para que se genere se debe de
##filtrar con los archivos correspondientes siguiendo la misma estrucura de las lineas
##de codigo de arriba, editanto los valores de los archivos y la Condicion como "Tumor"

#unimos sanos y tumor para etiquetar con el tejido correspondiente
cln_df <- rbind(cln.tsv_sano,
                cln.tsv_tumor)
cln_df$Tejido <- "Colon"

#Para trofoblasto todo es sano, no hay que separar para etiquetar
# tfb_df <- tfb.tsv
# tfb_df$Condicion <- "Sano"
# tfb_df$Tejido <- "Placenta"

####################################################################################
##Los siguientes pasos se deben repetir para cada tejido, editando el nombre       #
##de los objetos correspondientes en cada caso                                     #                                                             #
####################################################################################

##renombramos columnas
colnames(cln_df) <- c("Cromosoma","posicion","TE","Longitud_TE",
                        "Termino","TSD","Longitud_TSD","Transduction_source",
                        "Orientacion","Inversion_5","Copia_referencia",
                        "Info_gen","GT","Genotipo","Archivo","Condicion","Tejido")
                        #"Posicion","File_TE")   ##Esta linea se usa despues de limpiar y pegar columnas

##Hasta este punto debemos tener un data frame con la siguiente estructura y encabezados:
#(Los datos van a variar para cada set)
head(cln_df)

# # A tibble: 6 × 17
# Cromosoma posicion TE         Longitud_TE Termino TSD   Longitud_TSD Transduction_source
# <chr>        <dbl> <chr>      <chr>       <chr>   <chr> <chr>        <chr>              
#   1 chr7      95149737 SVTYPE=IN… SVLEN=265   END=95… TSD=… TSDLEN=16    TD_SRC=not_transdu…
# 2 chr12     20755736 SVTYPE=IN… SVLEN=277   END=20… TSD=… TSDLEN=15    TD_SRC=not_transdu…
# 3 chr12     66617790 SVTYPE=IN… SVLEN=151   END=66… TSD=… TSDLEN=18    TD_SRC=not_transdu…
# 4 chr16       239175 SVTYPE=IN… SVLEN=186   END=23… TSD=… TSDLEN=-1    TD_SRC=not_transdu…
# 5 chr17     27268029 SVTYPE=IN… SVLEN=264   END=27… TSD=… TSDLEN=14    TD_SRC=not_transdu…
# 6 chrX      88773290 SVTYPE=IN… SVLEN=259   END=88… TSD=… TSDLEN=10    TD_SRC=not_transdu…
# # ℹ 9 more variables: Orientacion <chr>, Inversion_5 <chr>, Copia_referencia <chr>,
# #   Info_gen <chr>, GT <chr>, Genotipo <chr>, Archivo <chr>, Condicion <chr>,
# #   Tejido <chr>

##Unimos posicion con Cromosoma
cln_df$position <- paste(cln_df$Cromosoma, cln_df$posicion, sep = "_")
##Unimos filename con TE
cln_df$file_TE <- paste(cln_df$Archivo, cln_df$TE, sep = "_")

##Limpiamos
cln_df$Cromosoma <- gsub("chr", "", cln_df$Cromosoma)
cln_df$TE <- gsub("SVTYPE=INS:ME:", "", cln_df$TE)
cln_df$Longitud_TE <- gsub("SVLEN=", "", cln_df$Longitud_TE)
cln_df$Termino <- gsub("END=", "", cln_df$Termino)
cln_df$TSD <- gsub("TSD=", "", cln_df$TSD)
cln_df$Longitud_TSD <- gsub("TSDLEN=", "", cln_df$Longitud_TSD)
cln_df$Transduction_source <- gsub("TD_SRC=", "", cln_df$Transduction_source)
cln_df$Orientacion <- gsub("STRAND=", "", cln_df$Orientacion)
cln_df$Inversion_5 <- gsub("INS_INV=", "", cln_df$Inversion_5)
cln_df$Copia_referencia <- gsub("REF_REP=", "", cln_df$Copia_referencia)
cln_df$Info_gen <- gsub("GENE_INFO=", "", cln_df$Info_gen)
cln_df$file_TE <- gsub("SVTYPE=INS:ME:", "", cln_df$file_TE)

#Eliminamos columnas irrelevantes para este caso
# ("posicion","Trasduction_source","Inversion_5" y "GT")
cln_df <- cln_df[ ,-c(2,8,10,13)]

##renombramos columnas 
colnames(cln_df) <- c("Cromosoma","TE","Longitud_TE",
                      "Termino","TSD","Longitud_TSD","Orientacion",
                      "Copia_referencia", "Info_gen","Genotipo",
                      "Archivo","Condicion","Tejido",
                      "Posicion","File_TE")

##Reordenamos columnas
cln_df <- cln_df[ , c(1,14,15,2,3,4,5,6,7,8,9,10,12,13,11)]

##Hasta este punto debemos tener un data frame con la siguiente estructura y encabezados:
#(Los datos van a variar para cada set)
head(cln_df)

# # A tibble: 6 × 15
# Cromosoma Posicion      File_TE TE    Longitud_TE Termino TSD   Longitud_TSD Orientacion
# <chr>     <chr>         <chr>   <chr> <chr>       <chr>   <chr> <chr>        <chr>      
#   1 1         chr1_8386567  SRR132… ALU   270         8386575 +GAA… 8            -          
#   2 1         chr1_27878440 SRR132… ALU   47          278784… NULL  -1           .          
# 3 1         chr1_64881262 SRR132… ALU   272         648812… +AAA… 18           +          
#   4 1         chr1_65590243 SRR132… ALU   457         655902… +AAA… 15           +          
#   5 1         chr1_74682360 SRR132… ALU   450         746823… +TGA… 11           -          
#   6 1         chr1_78643143 SRR132… ALU   87          786431… +GAA… 17           -          
#   # ℹ 6 more variables: Copia_referencia <chr>, Info_gen <chr>, Genotipo <chr>,
#   #   Condicion <chr>, Tejido <chr>, Archivo <chr>

##El siguiente paso unira los data frame individuales para cada tejido en uno solo
TEs_data <- rbind(
  #tfb_df, 
  #bst_df,
  #lvr_df,
  cln_df)


##Una vez que se crea TEs_data, las siguientes adiciones usaran el siguiente codigo
##para unir todo a la tabla madre "TEs_data"

TEs_data <- rbind(
  #bst_df,
  #tfb_df,
  TEs_data)

##Preparamos formatos                               
TEs_data$Longitud_TE <- as.numeric(TEs_data$Longitud_TE)
TEs_data$Longitud_TSD <- as.numeric(TEs_data$Longitud_TSD)

#Filtramos conforme a la literatura
TEs_data <- TEs_data %>% filter(Longitud_TSD > 2)    
TEs_data <- TEs_data %>% filter(Longitud_TE > 80)   


save(TEs_data, file = "TEs_data_18_10.RData")

############################################
##        PARA PODER CALCULAR RPM          #  
############################################

##lectura de datos
##Se leera un archivo .tsv con el numero de lineas por archivo fastq

#rc <- read_tsv("/home/david/Documents/rc_20.06.tsv", col_skip())
rc <- read_tsv("~/Documents/Datos_xTea/preparacion_tesis/reads-af.trim_04-08-23.tsv",
               col_skip())

##renombrar columnas
colnames(rc) <- c("pre_rc","Archivo")

##Dividimos entre 4, ya que de las 4 lineas del archivo fastq, solo nos importa 1
##que es la que corresponde a la secuencia biologica
rc <- rc %>%
  mutate(RC = pre_rc/4)


##unimos read counts a TEs_data por la columna de Archivo
TEs_data.rc <- merge(TEs_data, rc, by = "Archivo",
                     all.x =TRUE)

##Por ultimo guardamos 
# save(TEs_data.rc, file = "TEs_data.rc_08_08.RData")

##La estructura final junto con los encabezados deberia ser la siguiente:
glimpse(TEs_data.rc)

# Rows: 29,778
# Columns: 17
# $ Archivo          <chr> "CRR033294", "CRR033294", "CRR033294", "CRR033294", "CRR033294"…
# $ Cromosoma        <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"…
# $ Posicion         <chr> "chr1_29100313", "chr1_36750194", "chr1_43743239", "chr1_629962…
# $ File_TE          <chr> "CRR033294_ALU", "CRR033294_ALU", "CRR033294_ALU", "CRR033294_A…
# $ TE               <chr> "ALU", "ALU", "ALU", "ALU", "ALU", "ALU", "ALU", "ALU", "ALU", …
# $ Longitud_TE      <dbl> 297, 309, 247, 270, 379, 283, 275, 443, 450, 368, 236, 301, 255…
# $ Termino          <chr> "29100329", "36750209", "43743253", "62996294", "69625389", "79…
# $ TSD              <chr> "+AAAAAAATTAGCTGGG", "+AAAAAATTAAATTTA", "+AAAACAGAGCCCCC", "+C…
# $ Longitud_TSD     <dbl> 16, 15, 14, 15, 16, 7, 11, 18, 11, 17, 20, 17, 13, 14, 14, 15, …
# $ Orientacion      <chr> "+", "+", "+", "-", "-", "-", "+", "+", "+", "-", "-", "+", "+"…
# $ Copia_referencia <chr> "Fall_in_Alu_copy_12.1", "not_in_Alu_copy", "not_in_Alu_copy", …
# $ Info_gen         <chr> "intron:ENSG00000159023.21:EPB41", "not_gene_region", "intron:E…
# $ Genotipo         <chr> "1/1", "1/1", "0/1", "1/1", "0/1", "1/1", "0/1", "1/1", "1/1", …
# $ Condicion        <chr> "Sano", "Sano", "Sano", "Sano", "Sano", "Sano", "Sano", "Sano",…
# $ Tejido           <chr> "Breast", "Breast", "Breast", "Breast", "Breast", "Breast", "Br…
# $ pre_rc           <dbl> 69085228, 69085228, 69085228, 69085228, 69085228, 69085228, 690…
# $ RC               <dbl> 17271307, 17271307, 17271307, 17271307, 17271307, 17271307, 172…

##TEs_data.rc se actualizara cada vez que se agreguen nuevos datos, de ser asi
##lo comveniente es colocar la fecha en el nombre del archivo donde se guardara
##la nueva version de TEs_data.rc

##TEs_data.rc contiene las inserciones de TEs detectadas junto con todas las
## caracteristicas genotipicas de cada insercion y el valor del conteo de reads
## (RC) correspondiente a cada archivo, esta tabla nos servira para graficar
## en los scripts "p1...", p2...", "p3...", "p4...", "p5...", "p6..." y "p7..."

##Para cada caso, se hara un subset con las caracteristicas especificas para el
##objetivo del script usado