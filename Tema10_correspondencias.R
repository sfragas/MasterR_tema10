file.edit(".Rprofile")
# -------------------------------------------------------- #
## Tema 10.- Análisis de correspondencias (CA).  
# -------------------------------------------------------- #

### Introducción ------------------------------------------- 
  # Método de ordenación (unimodal) no supervisada o no restringida (CA),
  # o restringido (CCA).
  # Sensible a outliers (individuos con poca ocurrencia -especies raras-)
  # 

### Objetivos. --------------------------
  # Técnica de ordenación y reducción de dimensión.
  # Representar gráficamente la relación (o correspondencia) entre
  # muestras (sitios) y/o individuos (especies) mediante proximidad. 
  # Además, podemos restringir esta representación según ciertas variables explicativas
  # con lo cual la representación maximiza la relación con estas variables.

### Consideraciones ----------------------
  # Asume una relación unimodal con el gradiente (útil cuando los gradientes son largos). 
  # Si el gradiente es corto se recomiendan otras técnicas lineales (e.g. PCA o RDA). 
  # Si la relación no es lineal ni unimodal, sino que es monotónica, se recomienda el uso del método NMDS. 

### Clasificación ------------------------
  # CA simple y múltiple según analicemos 2 o más variables, respectivamente
  # CA canónico o restringido (CCA), con la restricción de que las dimensiones sean combinaciones lineales de un conjunto de variables explicativas.
  # CCA parcial (pCCA), podemos eliminar el efecto de un grupo de variables, antes de realizar el CCA.  
  # CCA sin tendencia (DCA) permite eliminar efectos de herradura creados por la presencia de ceros en la base de datos o no linealidades entre el 1º y 2º eje.  

  # Matriz Y respuesta (VD, nominal, distirbución unimodal): filas=muestras, columnas=individuos
  # Matriz X explicativas (VI, cuantitativas, rel. lineal con v.respuesta): filas=muestras, columnas=variables (para el CCA, pCA)
  # Matriz Z de control (cuantitativas, efectos de variables que deseamos controlar o eliminar): filas=muestras, columnas=variables (para el pCA). 
  
  # Para la matriz X (CA)
  # 1. Correspondencia entre niveles de filas o columnas.
  # 2. Correspondencia entre filas y columnas.
  # Relación de la matriz X con la matriz Y
  # 3. Correspondencia entre individuos y variables explicativas. (CCA)
  # 4. Ídem, controlando el efecto de la matriz Z (pCA).  


# -------------------------------------------------------- #
# Datos dune.   ------------------------------------------
# -------------------------------------------------------- #

#  Los datos corresponden a una investigación sobre los efectos del
#  manejo del suelo sobre la vegetación de dunas. 
#  Tenemos datos de cobertura de 30 especies en 20 sitios. 
#  El marco ambiental correspondiente se encuentra en los datos "dune.env", 
#     - A1, espesor del horizonte A1 del suelo; 
#     - Moisture, un factor ordinal con los niveles: 1<2<4<5; 
#     - Management, BF (agricultura biológica), HF (agricultura Hobby), NM (Gestión de Conservación de la Naturaleza) y SF (cultivo estándar); 
#     - Use, uso de la tierra Hayfield < Haypastu < Pasture; 
#     - Manure, uso de abonos 0<1<2<3<4.

# Objetivos ---------------
  # Qué sitios y para qué especies se obtiene una capacidad óptima de separar 
  # los sitios según su vegetación. 
  # Interpretar la tipología espacial y evaluar el potencial de las especies 
  # vegetales de las dunas como descriptores biológicos.  



#### Análisis de correspondencia (CA). -------------
# -------------------------------------------------------- #
  # Activamos la librería y los datos.
  library(ade4)
  library(vegan)
  library(gclus)
  library(ape)
  # library(scatterplot3d)
  
  data(dune)
  data(dune.env)
  
  str(dune[,1:5])
  str(dune.env)
  
# Idoneidad de los datos.  ---------------------------
  # - Debemos tener un conjunto de variables interdependientes (e.g. de especies), 
  # sin distinguir entre variables dependientes e independientes.
  # - Las variables deben medirse en escalas similares 
  # - No deben existir datos perdidos.
  # 
  # Observamos la frecuencia de ocurrencia de cada especie.
  apply(dune,2,summary)
  apply(dune,2,function(x) any(x=="NA"))#No hay ningun NA
  apply(dune,2,function(x) sum(x==0,na.rm=TRUE))
  # Vemos que no hay valores perdidos (NA), aunque sí hay muchos ceros.
  # Eliminamos aquellas especies que tienen baja frecuencia 
  var<-colnames(dune) #variables a analizar, todas las especies en este caso
  min.fo<-5 #umbral de frecuencia que queremos seleccionar
  z<-NULL;for(i in var){
    y1<-subset(dune,select=eval(parse(text=i))) 
    zi<-as.matrix(apply(y1,2,function(x,na.rm) sum(!x==0,na.rm=TRUE))<min.fo)
    z<-cbind(z,zi) 
  }
  
  y2<-dune[,z[1,]==FALSE]
  dim(dune)
  dim(y2)
  #se eliminaron 11 especies
  
  chisq.test(dune) # rechazamos la independencia, podemos continuar con el análisis

# Análisis de correspondencia para los datos de abundancia de especies de vegetación de dunas
  library(vegan3d)
  (dune.ca <- cca(dune))
  # Inercia total: es una medida de variabilidad total explicada.
  # Como es un CA clásico, no está restringido (variabilidad total de los datos)
 

  summary(dune.ca) 
  # la heterogeneidad de los datos (inercia total) es 2.115.
  # el primer eje representa el 25.3% de la variación total en la composición de especies
  # puntuaciones para especies y sitios

# Seleccionar el número de ejes
  evplot <- function(ev){
    # Broken stick model (MacArthur 1957)
    n <- length(ev)
    bsm <- data.frame(j=seq(1:n), p=0)
    bsm$p[1] <- 1/n
    for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
    bsm$p <- 100*bsm$p/n
    # Plot eigenvalues and % of variation for each axis
    op <- par(mfrow=c(2,1))
    barplot(ev, main="Eigenvalues", col="bisque", las=2)
    abline(h=mean(ev), col="red")
    legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
    barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
            main="% variation", col=c("bisque",2), las=2)
    legend("topright", c("% eigenvalue", "Broken stick model"), 
           pch=15, col=c("bisque",2), bty="n")
    par(op)
  }
  (ev2 <- dune.ca$CA$eig)
  #source("~/Dropbox/NEwR functions/evplot.R")
  evplot(ev2)

# Biplot: diagrama basado en la asociación entre objetos y variables (ambos como puntos).
  # proximidad (o distancia) = mayor asociación.
  # Peor representación en el origen (perfil promedio)
  # Escalado.  
    # 1. scaling=1 (perfiles de sitios o casos), es útil para interpretar las proximidades y encontrar gradientes o grupos de sitios. 
    # 2. scaling=2 (perfiles de especies o variables), permite identificar grupos de especies. 
    # 3. scaling=3, corresponde a un escalado simétrico.
  plot(dune.ca, main="default plot") #(por defecto scaling=2)
  #si quisiéramos un gráfico 3d
  par(mfrow=c(1,1))
  ordiplot3d(dune.ca,type="h", main="3D plot")
  # ¿Cuáles son las categorías más importantes? Las categorías que se localizan alrededor del centro del mapa (punto (0,0)) són aquellas que presentan un comportamiento "medio'', que no se diferencian tanto como las categorías que están en las partes más extremas del gráfico. Las categorías que se sitúan cerca del borde del gráfico son las más diferentes o las que muestran una variación más acusada respecto a las otras.
  # ¿Cómo se observa la relación entre las categorías de una variable? Viene dado por el ángulo que queda entre ambas categorías. Un ángulo de 0º representa una correlación positiva del 100%, uno de 180º muestra una correlación negativa del 100%, y uno de 90º (o 270º) significa que no hay correlación entre ellas.

  # Aquellas especies que están cercanas al orígen (punto (0,0)) puede significar que estas especies están en su óptimo en el medio del gradiente ecológico representado por los ejes, o que están presenten en cualquier sitio del gradiente.

# Gráficos un poco más elaborados  ------------
  
  # 1) graficar especies y sitios con símbolo proporcional a la abundancia de especies o 
  # el total del sitio, respectivamente.
  #las especies raras se observan con círculos pequeños
  par(mfrow=c(1,2))
  (sp.total<-colSums(dune))
  ordiplot(dune.ca, main="rare species",choices=c(1,2), type="none") # gráfico para diagramas de ordenación, mirar ordiplot(dune.ca, main="rare species",choices=c(1,2))
  orditorp(dune.ca, choices=c(1,2),display="species", col="blue") # agrega etiquetas que no se solapan
  points(dune.ca,  display="species", col="red", cex=100*sp.total/sum(sp.total))
  
  (site.total<-colSums(t(dune)))
  ordiplot(dune.ca, main="outlier samples", choices=c(1,2), type="none")
  orditorp(dune.ca, display="sites", col="blue")
  points(dune.ca, display="sites", col="red", cex=100*site.total/sum(site.total))

  # 2) identificar especies con un ajuste pobre
  #los círculos mayores indican un mejor ajuste
  par(mfrow=c(1,1))
  (dune.ca.gof<-goodness(dune.ca,choices=c(1,2), display="species",model="CA"))
  p<-ordiplot(dune.ca,  main="species - goodness",choices=c(1,2),display="species",type="n")
  points(dune.ca, display="species", pch=19, col="blue", cex=5*dune.ca.gof)
  text(p, "species",col="red", cex=0.5)

  # 3) identificar sitios con un ajuste pobre
  #los círculos mayores indican un peor ajuste
  (dune.ca.gof<-goodness(dune.ca, choices=c(1,2), display="sites",model="CA"))
  p<-ordiplot(dune.ca, main="sites - goodness", choices=c(1,2),display="sites",type="n")
  points(dune.ca, display="sites", pch=19, col="blue", cex=5*dune.ca.gof)
  text(p, "sites",col="red", cex=0.5)
  identify(p,"sites")

  # 4) La función "envfit" que ajusta los vectores o factores ambientales a la ordenación. 
  (ef <- envfit(dune.ca, dune.env, permutations = 999))
  # R2 para cada factor (numérico y categórico)
  # Centroides para las variables ambientales
  plot(dune.ca, display = "sites")
  plot(ef,cex=.8)

  # elipses por nivel de manejo
  p<-plot(dune.ca, display = "sites", type = "p")
  with(dune.env, ordiellipse(dune.ca, Management, kind = "se", conf = 0.95)) # elipses de dispersión mediante error estándar (se) de las puntuaciones promedio (ponderadas). La correlación (ponderada) define la dirección del eje principal de la elipse. 
  with(dune.env, ordispider(dune.ca, Management, col = "blue", label= TRUE)) #diagrama de araña donde cada punto está conectado con el centroide (ponderado) del grupo. 
  with(dune.env, ordihull(dune.ca, Management, col="blue", lty=2)) # une los ítems dentro de cada grupo
 

## Análisis de correspondencia sin tendencia o segmentado (DCA).  -------------
# -------------------------------------------------------- #

  # Podemos tener efectos en herradura por la no-linealidad en los gradientes cuando:
  # - Gradientes simples pero prolongados, aparecen como curvas o arcos en la ordenación. La solución es segmentar los ejes.   
  # - Las unidades muestrales se agrupan en los extremos del gradiente, para evitar esto habría que reescalar los ejes e igualar la varianza.   
  # - Las especies raras influyen mucho en el resultado, es necesario quitarles peso.
  
  # El DCA soluciona este problema mediante: segmentar, reescalar, y disminuir la influencia (peso) de especies raras

  par(mfrow=c(1,2))
  vare.cca <- cca(dune)
  plot(vare.cca, display="sites", main="CA")
  
  vare.dca <- decorana(dune)
  plot(vare.dca, display="sites", main="DCA")
  
  vare.dca 
  # Hemos dicho que el CA asume una relación unimodal con el gradiente (útil cuando los gradientes son largos). 
  # Lepš & Šmilauer 2001 indican que si el largo del eje es <3 mejor utilizar PCA, si es >4 utilicemos CA, y si tma un valor intermedio podemos utilizar cualquiera de los dos métodos.
  # pero hay que tener cuidado en el tipo de variables que tenemos
  # summary(vare.dca)
  
  # En este caso no se observan grandes diferencias entre ambas aproximaciones. 


## Análisis de correspondencia canónico (CCA). ---------- 
# -------------------------------------------------------- #

# 1) Para un conjunto determinado de variables (selección manual).---------------
  attach(dune.env)
  (ord <- cca(dune ~ A1 + Management, data=dune.env) )
  # Partición de la varibilidad (inercia):
  # Constrained: variabilidad de la matriz de vegetación que puede ser explicada por los ejes del CCA
  # Unconstrained: variabilidad de los residuales.
  # Aquí solo el 36,86% de la variabilidad total es capturada por el CCA, no es muy útil
  
# Importancia de los componentes:
  summary(ord) # si queremos datos de las coordenadas
  # el primer eje explica el 15.1% de la variabilidad restringida, 
  # el segundo el 11.2%...

# Contribuciones Absolutas. 
# - de una categoría a la inercia de un eje: representa el porcentaje de la Inercia Total en ese eje que se debe a cada nivel de la fila o columna. 
# - de una variable a la inercia de un eje: indica las variables que más contribuyen a la construcción de los ejes y ayudan a dar una interpretación a los ejes principales.


# gráficos ---------------------
# Triplot: ejes canónicos=combinación lineal de las v. explicativas
  # comparación CCA-CA
  par(mfrow=c(1,2))
  plot(ord, main="CCA") #método directo
  # Cómo se observa la relación entre variables?
  # Trazamos una línea desde una categoría de la variable fila (A en el gráfico) a través del origen.

  plot(vare.cca, main="envfit") #método indirecto
  plot(envfit(vare.cca, dune.env[,c("A1","Management")]))
  
  # CCA
  # display: sp=especies, wa=sitios (también puede ser "lc" pero es menos recomendado)
  # bp=flechas biplot o cn=centroides de los factores de restricción en lugar de flechas.
  # ver ?plot.cca o vegandocs("decision")
  par(mfrow=c(2,2))
  plot(ord, dis=c("wa","bp")) #por defecto se muestra las puntuaciones sp, wa y cn 
  par(mfrow=c(1,1))
  plot(ord, dis=c("wa","lc"))
  ordispider(ord) 
  plot(procrustes(vare.cca, ord))
  plot(vare.cca)

# Pruebas de permutación.---------------
  anova.cca(ord) #global
  anova.cca(ord, by="term", permu=200) #para cada término
  anova.cca(ord, by="mar") #para los efectos marginales
  anova.cca(ord, by="axis", perm=500) #para cada eje
  permutest(ord, first=T, permutations=1000) #determina si el primer eje del CCA es significativo, respecto a lo esperado por azar
  permutest(ord, permutations=1000)  # determina si hay relación total entre especies y ambiente

# Bondad de ajuste y selección del modelo. ---------------
  RsquareAdj(ord) #=ord; todavía no está disponible la opción R2 ajustado.
  goodness(ord,display = "species") #bondad de ajuste para especies
  goodness(ord,display = "sites")#bondad de ajuste para sitios
  par(mfrow=c(1,2))
  #identificar especies con un ajuste pobre
  #los círculos mayores indican un mejor ajuste
  (ord.gof<-goodness(ord, choices=c(1,2), display="species"))
  p<-ordiplot(ord,main="species - goodness", choices=c(1,2),display="species",type="n")
  points(ord, display="species",  pch=19, col="blue", cex=5*ord.gof)
  text(p, "species",col="red", cex=0.5)
  # identify(p,"species")
  
  #identificar sitios con un ajuste pobre
  #los círculos mayores indican un peor ajuste
  (ord.gof<-goodness(ord, choices=c(1,2), display="sites"))
  p<-ordiplot(ord, main="sites - goodness",choices=c(1,2),display="sites",type="n")
  points(ord, display="sites",  pch=19, col="blue", cex=5*ord.gof)
  text(p, "sites",col="red", cex=0.5)
  # identify(p,"sites")
  
  inertcomp(ord, prop = TRUE) # componentes de inertia, en proporción al total
  
  vif.cca(ord) # vif , todos son <3
  spenvcor(ord) # correlación especies-ambiente, sensitive to extreme scores (like correlations are), and very sensitive to overfitting or using too many constraints


# 2) Selección automática de variables explicativas.---------------
  m0 <- cca(dune ~  1, dune.env)
  m1 <- cca(dune ~ ., dune.env)
  (m <- ordistep(m0, scope = formula(m1),trace=0)) #basado en AIC, p-valores por permutación MC
  # variables seleccionadas Moisture + Management
  m$anova
  
  #VIF (>20, >10, >3)
  vif.cca(m1)
  vif.cca(m)

# comparación de los dos procedimientos
 RsquareAdj(ord)
 RsquareAdj(m)

## Análisis de correspondencia parcial (pCA).---------------
  
  # Los métodos de ordenación restringidos (CCA) pueden tener términos que se vuelvan no significativos al condicionarlos a otras variables.  
  # Aquí vamos a controlar la influencia de la humedad antes de analizar el efecto del A1 y del manejo.   
  
  (ord3 <- cca(dune ~ A1 + Management + Condition(Moisture), 
              data=cbind(dune,dune.env)))
  # summary(ord3)
  par(mfrow=c(1,1))
  plot(ord3)
  # flecha = dirección del gradiente.   
  # longitud = fuerza del gradiente.  

  # test de permutacion 
  anova.cca(ord3) #global
  anova.cca(ord3, by="term", permu=500) #para cada término
  # A1 deja de ser significativo al controlar la variable humedad.  

# Partición de la varianza (se basa en un rda y usa R2)
  (mod <- varpart(dune, ~ Management, ~ A1 + Moisture, data = dune.env))
  par(mfrow=c(1,2))
  plot(mod)
  showvarparts(2)
