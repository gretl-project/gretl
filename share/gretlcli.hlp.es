#
add
@Tests
Uso:            add listavar
Ejemplos:       add 5 7 9	          add xx yy zz

Se añaden las variables de la lista 'listavar' al modelo anterior y se
estima el nuevo modelo. Si se añade más de una variable, se presenta el
estadístico F para las variables añadidas (sólo para el método MCO) junto
con su valor-p. Un valor-p inferior a 0.05 indica que los coeficientes son
conjuntamente significativos al nivel del 5 por ciento.

#
addto
@Tests
Uso:            addto model_ID listavar
Ejemplo:        addto 2 5 7 9

Funciona como la instrucción "add", pero aquí hay que especificar un modelo
anterior (lo cual se hace utilizando su número de ID, que se presenta al
principio de los resultados del modelo) para tomarlo como base al añadir las
variables. En el ejemplo de arriba se añaden las variables con números 5, 7
y 9 al modelo 2.

#
adf
@Tests
Uso:	        adf orden nombrevar
Ejemplo:        adf 2 x1

Calcula los estadísticos para dos contrastes de Dickey-Fuller. En ambos casos
la hipótesis nula es que la variable en cuestión presenta una raíz unitaria.

El primero es un contraste t basado en el modelo
(1 - L)x(t) = m + g * x(t-1) + e(t).
La hipótesis nula es que g = 0.

En el segundo contraste (aumentado) se estima una regresión no restringida
(teniendo como regresores una constante, una tendencia temporal, el primer
retardo de la variable y "orden" retardos de la primera diferencia) y una
regresión restringida (quitando la tendencia temporal y el primer retardo).
El estadístico de contraste es F calculado como

[(SCRr - SCRnr)/2]/[SCRnr/(T - k)]

donde T es el tamaño muestral y k el número de parámetros del modelo no
restringido. Es necesario tener en cuenta que los valores críticos para
estos estadísticos no son los usuales.

#
ar
@Estimation
Uso:            ar retardos ; vardep varindep      o
                ar retardos ; -o vardep varindep
Ejemplo:        ar 1 3 4 ; y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento iterativo
generalizado de Cochrane-Orcutt (ver Ramanathan, sección 9.5). Las iteraciones
se terminan cuando las sucesivas sumas de cuadrados residuales no varían más del
0.005 por ciento o cuando han transcurrido 20 iteraciones. 'retardos' es una
lista de retardos de los residuos, que acaba con un punto y coma. En el
ejemplo de arriba el término de error se especifica como

u(t) = rho1 u(t-1) + rho3 u(t-3) + rho4 u(t-4) + et

'vardep' es la variable dependiente y 'varindep' es la lista de variables
independientes que van separadas con espacios. Utilice el número 0 para un
término constante. Si se usa la opción -o se presentará la matriz de varianzas
y covarianzas de los coeficientes de regresión. Los residuos de la regresión
transformada se guardan bajo el nombre 'uhat' y pueden recuperarse
usando la instrucción 'genr'.

#
arch
@Tests
Uso:            arch retardo vardep varindep
Ejemplos:       arch 4 1 0 2 4 6 7      o  arch 4 y 0 x1 x2 x3 x4

Esta instrucción contrasta la posibilidad de un ARCH en el modelo, del orden
especificado en "retardo" (que debe de ser entero). Si el estadístico de
contraste LM tiene un valor-p inferior a 0.10 también se realiza la estimación
ARCH. Si la varianza predicha de cualquier observación en la regresión auxiliar
no es positiva, se utiliza en su lugar el correspondiente 'uhat' (residuo)
al cuadrado. Luego, se estima por mínimos cuadrados ponderados el modelo 
original.


#
chow
@Tests
Uso:           chow obs
Ejemplos:      chow 25
               chow 1988.1

Primero debe ejecutarse una regresión por MCO. Crea una variable ficticia
que es igual a 1 desde el punto de corte especificado en "obs" hasta el final
de la muestra y 0 en el resto. También crea los términos de interacción entre
esta variable ficticia y las variables independientes originales. Se ejecuta
una regresión aumentada incluyendo estos términos y se calcula un estadístico F,
tomando la regresión aumentada como 'no restringida' y la original como
'restringida'. Este estadístico es adecuado para contrastar la hipótesis nula de
que no hay cambio estructural en el punto de ruptura indicado.


#
coint
@Tests
Uso:	        coint orden vardep varindep
Ejemplos:       coint 2 y x
                coint 4 y x1 x2

Realiza los contrastes de Dickey-Fuller para cada una de las variables listadas.
Para cada variable, considera la hipótesis nula de que la variable tiene una 
raíz
unitaria y utiliza para el contraste el orden de retardos dado.
Se estima la regresión cointegrante y se realiza un contraste ADF sobre
los residuos de esta regresión. También se proporciona el estadístico de 
Durbin-Watson para la regresión cointegrante.
	Hay que señalar que para ninguno de estos estadísticos de contraste
se pueden aplicar las tablas estadísticas usuales.

#
corc
@Estimation
Uso:          corc vardep varindep       o    corc -o vardep varindep
Ejemplos:       corc 1 0 2 4 6 7                corc -o 1 0 2 4 6 7
                corc y 0 x1 x2 x3               corc -o y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento iterativo de
Cochrane-Orcutt (ver Ramanathan, Sección 9.4) siendo 'vardep' la variable
dependiente y siendo 'varindep' una lista de variables independientes separadas
por espacios y acabando con un ;. Utilice el número 0 para un término constante. 




El proceso iterativo se detiene cuando valores sucesivos de rho no difieren en 
más
de 0.001 o cuando han transcurrido 20 iteraciones. Si se utiliza la opción -o, 
se muestra la matriz de covarianzas de los coeficientes de regresión. La 
regresión transformada final se calcula para el rango de observación
primobs+1 ultobs que esté actualmente en efecto. Los residuos de esta regresión
transformada se guardan con el nombre 'uhat'.


#
corr
@Statistics
Uso:         corr
             corr listavar

'corr' muestra los coeficientes de correlación para todos los pares de variables
que hay en el conjunto de datos (los valores perdidos, que se denotan por -999,
no se tienen en cuenta). 'corr listavar' muestra los coeficientes de correlación
para las variables listadas.


#
corrgm
@Statistics
Uso:          corrgm nombrevar o numerovar
              corrgm nombrevar o numerovar maxretardo

Muestra los valores de la función de autocorrelación para la variable
especificada (ver Ramanathan, sección 11.7). Es, por tanto, corr[u(t), u(t-s)],
donde u(t) es la observación t-ésima de la variable u y s es el número
de retardos. También se muestran las correlaciones parciales: estas tienen
descontado el efecto de los retardos que intervienen. La instrucción también
representa el correlograma y calcula el estadístico Q de Box-Pierce para
contrastar la hipótesis nula de que la serie es 'ruido blanco'. Este se
distribuye asintóticamente como una chi-cuadrado con un número de grados
de libertad igual al número de retardos utilizados.

Si se suministra un número entero 'maxretardo' la largura del correlograma
se limita a, como máximo, ese número de retardos, en caso contrario la largura
se determina automáticamente.


#
criteria
@Utilities
Uso:          criteria scr T k        p.ej. criteria 23.45 45 8

Dados scr (suma de cuadrados de los residuos), el número de observaciones (T) 
y el número de coeficientes (k), calcula los estadísticos de selección de
modelos (ver Ramanathan, Sección 4.4). T, k y scr pueden ser valores numéricos
o los nombres de variables definidas previamente.

#
critical
@Utilities
Uso:            critical t gl           p.ej. critical t 20
                critical X gl           p.ej. critical X 5
                critical F gln gld      p.ej. critical F 3 37
                critical d n            p.ej. critical d 24

Si el primer parámetro es t, X o F, muestra los valores críticos para la 
distribución t de student, chi-cuadrado o F respectivamente, para los niveles
de significación más comunes, utilizando los grados de libertad especificados.
Si el primer parámetro es d, muestra los valores superior e inferior del 
estadístico de Durbin-Watson al nivel de significación del 5 por ciento para el 
valor de n (número de observaciones) dado y para el rango de 1 a 5 variables 
explicativas.


#
cusum
@Tests
Uso:          cusum

Debe ejecutarse después de una estimación MCO. Desarrolla el contraste 
CUSUM de estabilidad de los parámetros. Se obtiene una serie de errores de
predicción (escalados) un paso hacia adelante, ejecutando para ello una serie
de regresiones: la primera regresión se realiza con las primeras k
observaciones (donde k es el número de parámetros en el modelo original)
y se usa para generar una predicción de la variable dependiente para la
observación k+1; la segunda utiliza las primeras k+1 observaciones y genera
una predicción para la observación k+2 y así sucesivamente. Se muestra y 
representa gráficamente la suma acumulada de los errores de predicción 
escalados. Se rechaza, al nivel de significación del 5 por ciento, la 
hipótesis nula de estabilidad de los parámetros si la suma acumulada 
sale fuera de las bandas del 95 por ciento de confianza.

También se proporciona el estadístico t de Harvey-Collier para contrastar 
la hipótesis nula de estabilidad de los parámetros. Ver Capítulo 7 del libro 
de Greene 'Econometric Analysis' para más detalles.

#
delete
@Dataset
Uso:          delete

Elimina la última (la de número más alto) variable del conjunto de datos actual.
Utilícela con precaución: no se pide confirmación. Puede ser útil para eliminar 
variables ficticias temporales. Sólo se puede eliminar la última variable.

#
diff
@Transformations
Uso:          diff listavar

Se toma la primera diferencia de cada variable en 'listavar' y el resultado se 
guarda en una nueva variable con el prefijo "d_". Así, "diff x y" crea las 
nuevas variables d_x=x(t)-x(t-1) y d_y = y(t) - y(t-1).

#
endloop
@Programming

Termina un bucle de simulación. Ver "loop".

#
eqnprint
@Printing
Uso:            eqnprint
                eqnprint -o

Debe ejecutarse después de una estimación por MCO. Imprime el modelo estimado,
en forma de ecuación LaTeX, a un fichero cuyo nombre tiene la estructura
"equation_N.tex", donde N es el número de modelos estimados hasta el momento
en la sesión actual. Este puede incorporarse en un documento LaTeX. Ver
también la instrucción 'tabprint'.

Si se utiliza la opción -o, el fichero LaTeX es un documento completo
(con preámbulo LaTeX) listo para procesar; en caso contrario debe incluirse
en un documento (que ya tenga el preámbulo).


#
fcast
@Prediction
Uso:            fcast primobs ultobs nuevonombrevar
                fcast nuevonombrevar
Ejemplo:        fcast ajustados

Los valores ajustados en la última regresión ejecutada se guardan bajo 
'nuevonombrevar'. Estos valores pueden mostrarse e representarse 
gráficamente. Las variables del lado derecho de la ecuación son las del modelo
original. No hay posibilidad de cambiar a otras variables. Si se especifican 
las observaciones inicial ('primobs') y final ('ultobs') la predicción se
restringe al rango especificado. Si se ha especificado un término de error 
autorregresivo (en 'hilu','corc' y 'ar'), la predicción es condicional un paso 
adelante e incorpora el proceso del error.


#
fcasterr
@Prediction
Uso:            fcasterr primobs ultobs
                fcasterr primobs ultobs -o

Después de estimar un modelo MCO que incluya una constante y al menos una 
variable independiente (estas restricciones pueden relajarse en algún punto),
se puede usar esta instrucción para mostrar los valores ajustados sobre el 
rango de observación especificado, junto con las desviaciones típicas estimadas
de estas predicciones y los intervalos de 95 por ciento de confianza. Si se 
utiliza la opción -o también se mostrarán los resultados en gráficos gnuplot.


#
fit
@Prediction
Uso:		fit

Esta orden (que debe seguir a una instrucción de estimación) es una atajo para
la instrucción 'fcast'. Genera valores ajustados para la muestra actual, 
basados en la última regresión y los guarda en una serie denominada "autofit".
En modelos de series temporales también muestra un gráfico gnuplot de los 
valores
actual y estimado de la variable dependiente contra el tiempo.

#
freq
@Statistics
Uso:          freq nombrevar (o numerovar)

Muestra la distribución de frecuencias de 'nombrevar' o 'numerovar';
también se proporcionan los resultados de una contraste chi-cuadrado de 
normalidad. El estadístico para este último es:

  tamaño_muestral * [asimetría^2/6 + (CURTOSIS - 3.0)^2/24.0]

Bajo la hipótesis nula de normalidad se distribuye como una chi-cuadrado
con 2 grados de libertad.

En modo interactivo se genera un gráfico gnuplot de la distribución.


#
genr
@Dataset
Uso:          genr nuevonombrevar = formula

Crea nuevas variables, normalmente por medio de transformaciones de
variables ya existentes. Ver también 'diff', 'logs', 'lags', 'ldiff',
'multiply' y 'square' como atajos.

Las operaciones aritméticas permitidas son, en orden de precedencia:
^(exponenciación); *, / y % (módulo o resto); + y -.

Los operadores booleanos (de nuevo en orden de precedencia) son:
! (NO lógico [NOT]), & (Y lógico [AND]), | (O lógico [OR]), >, <, =
y los símbolos compuestos != (no igual), >= (mayor o igual que) y
<= (menor o igual que).  Los operadores booleanos pueden usarse al
definir variables ficticias: por ejemplo (x > 10)
produce 1 si x(t) > 10, 0 en los demás casos.

Las funciones permitidas pertenecen a estos grupos:

-Funciones matemáticas standard: abs, cos, exp, int (parte entera), ln
(logaritmo natural: log es un sinónimo), sin (seno), sqrt (raíz cuadrada).

-Funciones estadísticas: mean (media aritmética), median (mediana),
var (varianza), sd(desviación típica o standard), sum, cov (covarianza),
corr (coeficiente de correlación), min (mínimo) y max (máximo).

-Funciones de series temporales: lag (retardo), lead (adelanto),
diff (primera diferencia), ldiff (log-diferencia, o primera diferencia del
logaritmo natural).

-Misceláneas: cum (acumulación), sort (ordenar), uniform, normal,
misszero (reemplazar los códigos de 'observación perdida'  por ceros),
zeromiss (operación inversa de misszero), pvalue (valor de probabilidad
para un estadístico dado contra una distribución especificada)
y mpow (elevar una serie a un exponente entero utilizando aritmética de 
precisión múltiple).

Todas las funciones anteriores, con las excepciones de cov, corr, uniform, 
normal, pvalue y mpow toman como único argumento o el nombre de una variable
(hay que tener en cuenta que en 'genr' no podemos referirnos a una variable por
su número de ID) o una expresión compuesta que se evalúa en una variable (p.ej.
ln((x1+x2)/2)). 'cov' y 'corr' necesitan dos argumentos y producen
respectivamente la covarianza y el coeficiente de correlación entre dos
variables especificadas. uniform() y normal(), que no tienen argumentos,
producen series pseudo-aleatorias obtenidas a partir de la distribución uniforme
(0-100) y la distribución normal standard respectivamente (ver también la
instrucción seed). La función pvalue() toma los mismos argumentos que la orden
pvalue (ver más abajo), pero en este contexto deben situarse comas entre sus
argumentos. La función mpow toma como argumentos el nombre de una serie de datos 
y un número entero positivo, que es exponente al cual se desea elevar la serie.

Además de los operadores y funciones mencionados hay algunos usos especiales
de 'genr':

* genr time crea una variable de tendencia temporal (1,2,3,...) denominada time.
* genr index crea una variable índice (1,2,3,...) denominada index.
* genr dummy crea variables ficticias hasta la periodicidad de los datos.
  P.ej. en el caso de datos trimestrales (periodicidad 4), el programa crea
  dummy_1 = 1 para el primer trimestre y 0 en los otros trimestres, dummy_2 = 1
  para el segundo trimestre y 0 en los otros trimestres, y así sucesivamente.
* genr paneldum crea un conjunto de variables ficticias especiales para el uso
  con un panel de datos (ver el manual de gretl para más detalle)
* Pueden recuperarse usando 'genr' varias variables internas que se definen
  mientras se ejecuta una regresión. Esto se hace de la siguiente forma:

  $ess         suma de cuadrados de los residuos
  $rsq         R-cuadrado no corregido
  $T           número de observaciones utilizadas en el modelo
  $df          grados de libertad
  $trsq        TR^2 (tamaño muestral por el R-cuadrado)
  $sigma       desviación típica de los residuos
  $lnl         log-verosimilitud (en modelos logit y probit)
  coeff(var)   coeficiente estimado para var
  stderr(var)  desviación típica estimada para var
  rho(i)       coeficiente de autorregresión de orden i-ésimo de los residuos
  vcv(xi,xj)   covarianza entre los coeficientes de las variables xi y xj

La variable interna $nobs contiene el número de observaciones en el rango
muestral actual, que puede ser o puede no ser igual al $T del último modelo.

La variable interna $pd contiene la periodicidad o frecuencia de los datos
(por ejemplo, 4 para datos trimestrales, 12 para mensuales).

La variable interna t sirve para referirse a las observaciones, comenzando
en 1. Así uno puede hacer "genr dum15 = (t=15)" para generar una variable
ficticia con valor 1 para la observación 15 y 0 en el resto.

Ejemplos de instrucciones 'genr':

  genr y = x1^3          [x1 al cubo]
  genr y=ln((x1+x2)/x3)  [argumento compuesto para una función ln]
  genr z=x>y             [asigna z(t) a 1 si x(t) > y(t), en otro caso a 0]
  genr y=x(-2)           [x retardada 2 periodos]
  genr y=x(2)            [x adelantada 2 periodos]
  genr y = mean(x)       [media aritmética]
  genr y = diff(x)       [y(t) = x(t) - x(t-1)]
  genr y = ldiff(x)      [y = ln(x(t)) - ln(x(t-1))]
                          ldiff(x) es la tasa de crecimiento instantánea de x.
  genr y = sort(x)       [ordena x en orden creciente y lo guarda en y]
  genr y = - sort(-x)    [ordena x en orden decreciente]
  genr y = int(x)        [trunca x y guarda su valor entero como y]
  genr y = abs(x)        [guarda los valores absolutos de x]
  genr y = sum(x)        [suma los valores de x excluyendo los valores perdidos 
-999]
  genr y = cum(x)        [acumula x: y(t) es la suma de x hasta t]
  genr aa = $ess         [aa = suma de cuadrados de los residuos de la última 
regresión]
  genr x = coeff(sqft)   [guarda en x el coeficiente de la variable sqft 
obtenido 
                          en el último modelo]
  genr rho4 = rho(4)     [guarda en rho4 el coeficiente autorregresivo de cuarto 



orden
                          obtenido del último modelo (supone un modelo ar)]
  genr cv=vcv(x1, x2)    [covarianza entre los coeficientes de x1 y x2 en el 
último modelo]
  genr x=uniform()/100   [variable pseudo-aleatoria uniforme, de rango 0 a 1]
  genr x=3*normal()      [variable pseudo-aleatoria normal, de media 0 y desv. 
típica 3]
  genr x=pvalue(t,20,1.4)[valor p para 1.4, bajo la distribución t con 20 grados 



de libertad]

Sugerencias sobre variables ficticias:

* Supongamos que x se codifica con los valores 1, 2, o 3 y Vd desea tres 
variables ficticias d1 = 1 si x = 1, 0 en otro caso, d2 = 1 si x = 2, y así 
sucesivamente. Para crear estas, utilice las instrucciones 
genr d1 = (x=1), genr d2 = (x=2), y genr d3 = (x=3).

* Para obtener la serie z = máx(x,y) haga genr d=x>y y genr z=(x*d)+(y*(1-d))

#
gnuplot
@Graphs
Uso:            gnuplot yvar1 xvar [ opción ]
                gnuplot yvar1 yvar2 xvar [ opción ]
		gnuplot yvar xvar ficticia -z

En los dos primeros casos las variables yvars se representan contra xvar.
Si se proporciona la opción -o el gráfico utilizará líneas; si se da la
opción -m el gráfico utiliza impulsos (líneas verticales); en los demás
casos se usan puntos.

En el tercer caso, yvar se representa contra xvar mostrando los puntos en
diferentes colores dependiendo de si el valor de 'ficticia' es 1 o 0.

Para crear un gráfico de serie temporal, pedir "gnuplot yvars time".
Si no existe la variable "time", se generará automáticamente.  Se crearán
variables ficticias especiales para representar datos trimestrales y mensuales.

En modo interactivo, el resultado se pasa a gnuplot para que lo muestre
en pantalla.  En modo 'batch' se graba un fichero de nombre gpttmp<n>.plt,
donde <n> es un número entre 1 y 99. Más tarde, pueden generarse los gráficos
usando la instrucción de consola "gnuplot gpttmp<n>.plt".

#
graph
@Graphs
Uso:            graph var1 var2
                graph -o var1 var2
                graph var1 var2 var3

En los dos primeros ejemplos, la variable var1 (que puede ser un nombre o un
número) se representa (eje y) contra var2 (eje x). La opción -o hará el gráfico
con 40 filas y 60 columnas, sin ella el gráfico será de 20 por 60 (salida de
pantalla). En el tercer ejemplo, las dos, var1 y var2 se representarán (sobre
el eje y) contra var3. Esto es especialmente útil para representar los valores
observados y predichos contra el tiempo.

#
hausman
@Tests
Uso:          hausman

Este contraste sólo está disponible después de haber estimado un modelo 
utilizando la orden "pooled" (ver también las instrucciones "panel" y 
"setobs"). Contrasta el modelo combinado simple contra las alternativas 
principales, los modelos de efectos fijos y de efectos aleatorios.

En el modelo de efectos fijos se añade una variable ficticia para todas las
unidades de sección cruzada menos una, permitiendo al término constante de la
regresión que varíe a través de las unidades. En el modelo de efectos 
aleatorios se descompone la varianza residual en dos partes, una parte 
específica de la unidad de sección cruzada y la otra específica de la 
observación particular. (Este estimador sólo puede calcularse si el número 
de unidades de sección cruzada en el conjunto de datos es mayor que el número
de parámetros a estimar. Se presenta el estadístico LM de Breusch-Pagan para 
contrastar la hipótesis nula (de que el estimador MCO combinado es adecuado)
contra la alternativa de efectos aleatorios.

El modelo de MCO combinados puede ser rechazado contra ambas alternativas, 
efectos fijos y  efectos aleatorios. Si el error específico de grupo 
-o de unidad- está incorrelacionado con las variables independientes, el 
estimador de efectos aleatorios es más eficiente que el estimador de efectos 
fijos; en caso contrario, el estimador de efectos aleatorios es inconsistente
y será preferible utilizar el estimador de efectos fijos. La hipótesis nula 
para el contraste de Hausman es que el error específico de grupo no está 
correlacionado con las variables explicativas (y por tanto, que es preferible 
el modelo de efectos aleatorios). Un valor p bajo para este contraste es 
una indicación en contra del modelo de efectos aleatorios y a favor del de 
efectos fijos.


#
hccm
@Estimation
Uso:           hccm vardep varindep
            o  hccm -o vardep varindep

(Heteroskedasticity Consistent Covariance Matrix)
Presenta las estimaciones de MCO con las desviaciones típicas de los
coeficientes obtenidas por medio de una estimación de la matriz de varianzas
y covarianzas consistente ante heterocedasticidad. Utiliza para ello el
método "jacknife" de MacKinnon-White.


#
help
@Utilities
help              proporciona una lista de instrucciones gretl
help nombreinst   describe la instrucción 'nombreinst' (p.ej. help smpl)

#
hilu
@Estimation
Uso:            hilu vardep varindep       o    hilu -o vardep varindep
Ejemplos:       hilu 1 0 2 4 6 7                hilu -o 1 0 2 4 6 7
                hilu y 0 x1 x2 x3               hilu -o y 0 x1 x2 x3

Calcula las estimaciones de un modelo utilizando el procedimiento de
búsqueda de Hildreth-Lu (se hace el ajuste fino usando el método de
Cochrane-Orcutt) siendo 'vardep' la variable dependiente y 'varindep'
una lista de variables independientes separadas por espacios. Utilice el
número 0 para incluir un término constante. Se representa gráficamente la
suma de cuadrados de los residuos del modelo transformado contra los
valores de rho desde -0.99 hasta 0.99. Si se usa la opción -o se mostrará
también la matriz de varianzas y covarianzas de los coeficientes. Finalmente,
la última regresión transformada se estima para el rango de observación
primobs+1 ultobs que esté en efecto. Los residuos de esta regresión
transformada se guardan bajo el nombre 'uhat'.


#
hsk
@Estimation
Uso:            hsk vardep varindep
            o   hsk -o vardep varindep

Calcula estimaciones corregidas de heterocedasticidad y sus estadísticos 
asociados. Se ajusta una regresión auxiliar para el logaritmo de los 
cuadrados de los residuos (utilizando los cuadrados de las variables 
independientes, pero no sus productos cruzados) y a partir de esta 
estimación se obtienen los estimadores de mínimos cuadrados ponderados del
modelo inicial. Si se usa la opción -o, se mostrará también la matriz de
varianzas y covarianzas estimada de los coeficientes de la regresión.


#
if
@Programming
Uso:            if condición_boolena
                  instrucción1
                  instrucción2 ...
                endif

Las instrucciones gretl que están dentro del bloque "if ... endif" se 
ejecutan si y sólo si la condición booleana se evalúa como cierta (no cero).
Para conocer la sintaxis de las condiciones booleanas en gretl, ver la 
instrucción 'genr'. Opcionalmente, la orden "endif" puede ir precedida de
una orden "else" (en una línea aparte para ella sola, como las instrucciones 
"if" y "endif"), seguida de un bloque de instrucciones a ejecutar si
la condición booleana original se evalúa como falsa (cero). Los bloques
entre "if", "else" y "endif" pueden contener tantas órdenes como Vd
quiera y estas condiciones pueden estar anidadas. Cada orden "if" debe 
estar emparejada con una orden "endif".


#
import
@Dataset
Uso:            import archivo_csv
                import -o archivo_box

Sin la opción -o, importa datos desde un archivo que tenga formato de
'valores separados por comas' (CSV), como por ejemplo los que se pueden
escribir fácilmente desde cualquier programa de hoja de cálculo. El
archivo debería tener en la primera línea nombres de variables y, en
el resto, una matriz de datos rectangular. Las variables deberían estar
alineadas "por observación" (una columna para cada variable; cada fila
representa una observación).

Con la opción -o, lee un fichero de datos en formato BOX1, como los que
se obtienen utilizando el servicio de extracción de datos del 'US Bureau
of the Census'.

#
info
@Dataset
Uso:          info

Muestra la información contenida en el archivo de cabecera correspondiente
al fichero de datos actual. Esta información debe estar limitada entre los 
caracteres "(*" y "*)", estando situados estos marcadores en líneas separadas.

#
labels
@Dataset
Uso:          labels

Muestra las etiquetas informativas de las variables que se hayan definido
utilizando la instrucción 'genr'.

#
lad
@Estimación
(Estimador de mínima desviación absoluta)
Uso:          lad vardep varindeps

Calcula una regresión que minimiza la suma de las desviaciones absolutas entre 
los valores observados y ajustados de la variable dependiente. Las estimaciones 
de los coeficientes se obtienen utilizando el algoritmo simplex de 
Barrodale-Roberts; surge un mensaje de aviso si la solución no es única. Las 
desviaciones típicas se obtienen por un método 'bootstrap' con 500 
extracciones.

#
lags
@Transformations
Uso:          lags listavar

Crea variables nuevas que son valores retardados de cada una de
las variables que haya en 'listavar'. El número de variables retardadas
que se crean es igual a la periodicidad. Por ejemplo, si la periodicidad
fuera 4 (datos trimestrales), la orden 'lags x y' creará x_1  = x(-1),
x_2 = x(-2), x_3 = x(-3) y x_4 = x(-4); y de igual forma para y. Estas
variables deben referenciarse de forma exacta, es decir, con el carácter de
subrayado.


#
ldiff
@Transformations
Uso:          ldiff listavar

Se calcula la primera diferencia del logaritmo natural de cada 
variable de 'listavar' y el resultado se guarda en una variable nueva
que lleva el prefijo "ld_". Así por ejemplo, "ldiff x y" crea las 
nuevas variables

               ld_x = ln[x(t)] - ln[x(t-1)]
               ld_y = ln[y(t)] - ln[y(t-1)].

#
lmtest
@Tests
Uso:            lmtest
                lmtest -o

Debe utilizarse justo después de una instrucción 'ols'. Calcula el 
estadístico de contraste del multiplicador de Lagrange (LM) para las 
hipótesis alternativas de no linealidad y de heterocedasticidad 
(Contraste de White) o, si se utiliza la opción -o, para correlación de 
orden hasta la periodicidad. También se muestran los coeficientes
de la regresión auxiliar correspondiente. (Ver capítulos 7, 8 y 9 del
libro de Ramanathan para más detalles).

Sólo se usan los cuadrados de las variables independientes, y no sus
productos cruzados. No se pueden obtener los estadísticos de
contraste LM si la generación interna de los cuadrados causa
multicolinealidad exacta.


#
logit
@Estimation
Uso:          logit vardep varindeps

Regresión Logit: la variable dependiente debería ser una variable binaria.
Se obtienen las estimaciones por máxima verosimilitud de los coeficientes de 
'varindeps' por medio de mínimos cuadrados iterados (Método EM, o de 
expectativa-maximización). Como el modelo no es lineal, las pendientes dependen
de los valores de las variables independientes: las pendientes que se muestran 
se evalúan en la media de dichas variables. El estadístico Chi-cuadrado 
contrasta la hipótesis nula de que todos los coeficientes, excepto la 
constante, son cero.


#
logs
@Transformations
Uso:          logs listavar

Se obtiene el logaritmo natural de cada variable en 'listavar' y el
resultado se guarda en una nueva variable con prefijo "l_". Así por ejemplo,
"logs x y" crea las nuevas variables l_x = ln(x) y l_y = ln(y).


#
loop
@Programming
Uso:            loop número_de_veces
                loop while condición
		loop for i=principio..final
Ejemplos:       loop 1000
		loop while essdiff > .00001
		loop for i=1991..2000

Esta instrucción (de guión) da acceso a un modo especial, en el cual el
programa acepta órdenes para repetirlas o un número de veces especificado,
o mientras se satisfaga una condición, o para valores sucesivos de la
variable índice i (interna). Dentro de un bucle, sólo se pueden utilizar
seis instrucciones: genr, ols, print, smpl, store y summary (store no puede
usarse en un bucle de 'while'). Con genr y ols se pueden hacer muchas cosas.
Se sale del modo de introducción de órdenes de bucle con la instrucción
"endloop": en este punto se ejecutarán las órdenes de todo el bloque. Los
bucles construidos mediante "loop" no pueden estar anidados.

La instrucción ols proporciona un cuadro de resultados especial, dependiendo
del tipo de bucle. Si se especifica un "número_de_veces" no se muestran los
resultados de cada regresión particular, en su lugar se ofrece
(a) el valor medio de cada coeficiente estimado a lo largo de todas las
    iteraciones.
(b) las desviaciones típicas de esos coeficientes estimados.
(c) el valor medio de las desviaciones típicas estimadas de cada coeficiente.
(d) las desviaciones típicas de las desviaciones típicas estimadas.
Esto tiene sentido sólo si hay alguna entrada aleatoria en cada iteración.
Esta instrucción está diseñada para el análisis de Monte Carlo. Si se da una
condición "while" se muestran los resultados del modelo especificado
obtenidos en la última estimación del bucle: esto está diseñado para mínimos
cuadrados iterativos.

La instrucción "print" también se comporta de modo diferente en el contexto
de un bucle con "número_de_veces". Muestra la media y desviación típica de
la variable a lo largo de las repeticiones del bucle. Está pensado para
usarlo con variables que toman un solo valor en cada iteración, por ejemplo
la scr (suma de cuadrados de los residuos) de una regresión. La
instrucción "print" funciona de la forma usual para los demás tipos de bucles.

La instrucción "store" (sólo se puede usar una por bucle y sólo en los bucles
de tipo "número_de_veces") escribe los valores de las variables especificadas,
en cada iteración del bucle, al fichero especificado. De este modo, guarda una
copia completa de las variables. Luego puede leerse mediante gretl este fichero
de datos para analizarlo.

Ejemplo de programa de bucle (Monte Carlo):

   genr x = uniform()
   loop 100
   genr u = normal()
   genr y = (10*x) + (20*u)
   ols y const x
   genr r2 = $rsq
   print r2
   genr a = coeff(const)
   genr b = coeff(x)
   store foo.gdt a b
   endloop


#
meantest
@Tests
Uso:            meantest x1 x2
                meantest x1 x2 -o

Calcula el estadístico t para contrastar la hipótesis nula de que las medias
poblacionales de las variables x1 y x2 son iguales y muestra su valor p. Sin
la opción -o, el estadístico se calcula bajo el supuesto de que las varianzas
son iguales para las dos variables; con la opción -o se supone que las
varianzas son distintas. (La opción implicara diferencia sólo si hay
diferentes números de observaciones no perdidas para las dos variables)


#
mpols
@Estimation
(Mínimos cuadrados ordinarios de alta precisión)
Uso:            mpols vardep varindep
Ejemplos:       ols 1 0 2 4 6 7
                ols y 0 x1 x2 x3

Calcula los estimadores de mínimos cuadrados ordinarios con 'vardep' como
variable dependiente y 'varindep' como lista de variables independientes,
utilizando para ello operaciones aritméticas con precisión múltiple. Las 
variables pueden especificarse por su nombre o por su número; utilice el número 
cero para el término constante. Esta instrucción sólo está disponible si gretl 
se configura con soporte para GMP, la biblioteca GNU de precisión múltiple.

Para estimar un ajuste polinómico utilizando aritmética de precisión múltiple al 
generar las potencias necesarias de la variable independiente utilice, por 
ejemplo, la forma

mpols y 0 x ; 2 3 4

Esto hace una regresión de y sobre x, x al cuadrado, x al cubo y x a la cuarta 
potencia. Es decir, los números a la derecha del punto y coma(que deben ser
enteros y positivos) especifican las potencias de x a utilizar. Si se especifica 
más de una variable independiente, la última variable antes del punto y coma es 
la que será elevada a las potencias que se indican.

#
multiply
@Transformations
Uso:            multiply x sufijo vars
Ejemplos:       multiply invpop pc 3 4 5 6
                multiply 1000 big x1 x2 x3

Las variables de la lista "vars" (referenciadas por nombre o número) se
multiplican por x, que puede ser o un valor numérico o el nombre de una
variable previamente definida. Los productos se nombran con el prefijo
que se suministre (máximo tres caracteres). Si hace falta, se truncan los
nombres de las variables originales. Por ejemplo, supongamos que Vd desea
crear las versiones "per cápita" de algunas variables y Vd tiene la variable
"pob" (población). Las instrucciones adecuadas son

  genr invpob = 1/pob
  multiply invpob pc renta gasto

que crearán las variables

rentapc = renta * invpob, gastopc = gasto * invpop.

#
noecho
@Programming
Uso:          noecho

Suprime el eco normal al introducir las instrucciones, cuando se ejecuta
un guión gretl. Ver también la orden "print" (variante de cadena literal).


#
nulldata
@Dataset
Uso:            nulldata tamaño_serie
Ejemplo:        nulldata 100

Establece un conjunto de datos vacío, que contiene sólo una constante, con
periodicidad 1 y el numero de observaciones especificado. Este puede usarse
al hacer simulación: algunas de las órdenes genr (p.ej. genr uniform(),
genr normal(), genr time) generarán datos ficticios desde cero para llenar
el conjunto de datos. La orden "nulldata" también puede ser útil al utilizarla
conjuntamente con "loop".


#
ols
@Estimation
Uso:            ols vardep varindep       o      ols -o vardep varindep
Ejemplos:       ols 1 0 2 4 6 7                  ols -o 1 0 2 4 6 7
                ols y 0 x1 x2 x3                 ols -o y 0 x1 x2 x3

Calcula las estimaciones de mínimos cuadrados ordinarios con 'vardep'
como variable dependiente y siendo 'varindep' una lista de variables 
independientes. Con la opción -o se mostrará la matriz de covarianzas 
de los coeficientes de regresión. Las variables pueden introducirse 
mediante nombres o números. Utilice el número cero para incluir un 
término constante. El programa muestra también los valores p para los 
estadísticos t (a dos colas) y F. Un valor p inferior a 0.01 indica 
significatividad al nivel del 1 por ciento y se denota mediante tres *.
Dos * indican significatividad a niveles entre el 1 y el 5 por ciento.
También se muestran los estadísticos de selección de modelos descritos
en el libro de Ramanathan, sección 4.4.

 
#
omit
@Tests
Uso:            omit varlist
Ejemplos:       omit 5 7 9
		omit xx yy zz

Se ejecuta después de estimar un modelo. Las variables de la lista 'varlist'
se omitirán del modelo anterior y se estimará el nuevo modelo. Si se omite más
de una variable, se muestra el estadístico F de Wald para las variable omitidas
junto con su valor p (sólo para el método MCO). Un valor p inferior a 0.05 
implica que los coeficientes son conjuntamente significativos al nivel del 
5 por ciento.


#
omitfrom
@Tests
Uso:            omitfrom ID_de_modelo varlist
Ejemplo:        omitfrom 2 5 7 9

Funciona como la orden "omit", pero aquí se puede especificar un modelo previo
(usando su número de ID, que se muestra al principio de los resultados del 
modelo) para tomarlo como base al omitir las variables. En el ejemplo de 
arriba se omiten las variables con números 5, 7 y 9 del modelo 2.


#
open
@Dataset
Uso:          open datafile

Abre un fichero de datos. Si ya hay un fichero de datos abierto, se reemplaza 
por el nuevo. El programa intentará detectar el formato del fichero de datos
("nativo", CSV o BOX1) y lo tratará como corresponda.


#
panel
@Dataset
Uso:            panel
                panel -s
	        panel -c

Propone que el conjunto de datos actual sea tratado como un panel (combinando
datos de sección cruzada y de series temporales). Sin ninguna opción o con 
la opción -s, los datos se consideran formados por series temporales apiladas
(bloques sucesivos de datos contienen series temporales para cada unidad de 
sección cruzada). Con la opción -c, los datos se leen como datos de sección 
cruzada apilados (bloques sucesivos contienen secciones cruzadas para cada 
periodo temporal). Ver también la instrucción "setobs".


#
pergm
@Statistics
Uso:            pergm varname
                pergm varname -o

Calcula y muestra (y, en modo interactivo, representa) el espectro de la 
variable especificada. Sin la opción -o se obtiene el periodograma muestral;
con la opción -o se usa una ventana de retardos de Bartlett de tamaño
2*sqrt(tamaño muestral) para estimar el espectro (ver Capítulo 18 del libro
"Análisis Econométrico" de Greene). Cuando se muestra el periodograma 
muestral, se ofrece también un contraste t de integración fraccional: la
hipótesis nula es que el orden de integración de la serie es cero.


#
plot
@Graphs
Uso:            plot x1       plot x1 x2
                plot 3 7      plot -o x1 x2

Representa ( en un gráfico de tipo texto) los valores de los datos de 
las variables especificadas, para el rango de observación actualmente en 
efecto. Cada línea se refiere a una observación y los valores se 
representan horizontalmente. Si se utiliza la opción -o, x1 y x2 se
representan en la misma escala, en otro caso, cada una de ellas se escala
adecuadamente. La opción -o sólo debería usarse si las variables tienen 
aproximadamente el mismo rango de valores (p.ej. la variable dependiente 
observada y su predicción)


#
pooled
@Estimation
Uso:         pooled vardep varindeps

Estima un modelo mediante MCO (ver la instrucción "ols" para detalles sobre 
su sintaxis) y lo marca como 'modelo de panel' (o 'combinado') de manera que
la opción de contraste "diagnósticos de panel" esté disponible. Para consulta 
sobre dichos diagnósticos ver la orden "hausman".

#
print
@Printing

Imprime (=muestra en pantalla) los valores de las variables especificadas para
el rango primobs-ultobs actual, o muestra una cadena literal.

print              escribe el conjunto de datos actual en forma tabular
print 3 6          escribe los valores de las variables números 3 y 6
print x y z        escribe los valores de las variables denominadas x, y y z

Si se utiliza la opción -o las variables se muestran en columnas, en caso
contrario, se muestran en bloques consecutivos. Si se indica la opción "-t", los 
datos se muestran con 10 valores significativos.

print "Cadena literal" escribe la cadena especificada.  Las comillas finales no
son necesarias, pero se necesitan las iniciales para obtener este resultado.

#
probit
@Estimation
Uso:          probit vardep varindeps

Regresión probit: la variable dependiente debe ser una variable binaria. Se
calculan, mediante mínimos cuadrados iterativos (método EM ó de
expectativa-maximización), los estimadores de máxima verosimilitud de los
coeficientes de 'varindeps'. Como el modelo es no lineal, las pendientes
dependen de los valores de las variables independientes: las pendientes que
se muestran se evalúan en las medias de estas variables. El estadístico
Chi-cuadrado contrasta la hipótesis nula de que todos los coeficientes, excepto
el término constante, son cero.


#
pvalue
@Utilities
Uso en modo interactivo:      pvalue
Uso en modo batch (o en modo consola gretl):
     Distribución normal:     pvalue 1 valor_x
     Distribución t:          pvalue 2 g.l. valor_x
     Chi-cuadrado:            pvalue 3 g.l. valor_x
     Distribución F:          pvalue 4 gln gld valor_x
     Distribución Gamma:      pvalue 5 media varianza valor_x

Calcula el área a la derecha del 'valor_x' en la distribución especificada.
g.l. son los grados de libertad, gln son los grados de libertad del numerador,
gld son los grados de libertad del denominador.

#
quit
@Utilities

Salir de gretl. ('q' es un atajo; 'x' sale sin preguntar si se deben guardar
los resultados.)

#
reset
@Tests
Uso:          reset

Debe utilizarse inmediatamente después de una instrucción ols. Realiza el
contraste RESET de especificación de modelos (no linealidad) de Ramsey añadiendo
a la regresión el cuadrado y el cubo de los valores ajustados y calculando el
estadístico F para la hipótesis nula de que los parámetros de los dos términos
añadidos son cero.

#
rhodiff
@Transformations
Uso:           rhodiff listarho ; listavar
Ejemplos:      rhodiff .65 ; 2 3 4
               rhodiff r1 r2 ; x1 x2 x3

Crea las transformaciones rho-diferenciadas de las variables contenidas en
'listvar'(referenciadas mediante número o nombre) y las añade al conjunto de
datos utilizando el prefijo # para las nuevas variables. Dada la variable
v1 en listavar y r1 y r2 en listarho, se crea

v1# = v1(t) - r1*v1(t-1) - r2*v1(t-2)

Las entradas de listarho pueden darse como valores numéricos o mediante
los nombres de variables previamente definidas.

#
rmplot
@Graphs
(Gráfico Rango-Media)
Uso:           rmplot nombrevar

Esta instrucción crea un gráfico simple para ayudar a decidir si una serie
temporal, y(t), tiene o no varianza constante. Se toma la muestra completa
t=1,...,T y se divide en pequeñas submuestras de tamaño arbitrario k [gretl
elige k=sqrt(T)]. La primera submuestra se forma con y(1),...,y(k), la segunda
con y(k+1),...,y(2k), y así sucesivamente. Para cada submuestra se calcula la
media muestral y el rango (=máximo-mínimo) y se construye un gráfico con las
medias en el eje horizontal y los rangos en el vertical. De esta forma, cada
submuestra está representada por un punto en este plano. Si la varianza de la
serie fuera constante los rangos de las submuestras no deberían depender de sus
medias; si se observa que los puntos se aproximan a una recta con pendiente
creciente, esto sugiere que la varianza de la serie aumenta cuando la media
aumenta; si los puntos se aproximan a una recta con pendiente decreciente, esto
sugiere que la varianza está disminuyendo cuando la media aumenta.

Además del gráfico, gretl presenta las medias y los rangos para cada submuestra,
el coeficiente estimado para la pendiente en una regresión MCO de los rangos
sobre las medias y el valor p para el contraste de la hipótesis nula de que esta
pendiente es cero. Si el coeficiente de pendiente es significativo al nivel de
significación del 10 por ciento, en el gráfico se muestra también la recta
ajustada en la regresión de los rangos sobre las medias.

#
run
@Programming
Uso:          run fichero

Si el fichero "fichero" contiene instrucciones gretl, esta orden (que se
invoca desde dentro de gretl) las ejecutará de una en una. Esta es una forma
muy útil de ejecutar instrucciones 'batch' desde una sesión interactiva.

#
runs
@Tests
Uso:            runs nombre_var

Realiza el contraste de rachas no paramétrico sobre aleatoriedad de la
variable especificada. Si Vd desea contrastar la aleatoriedad de las
desviaciones respecto a la mediana de una variable x1 que tiene una mediana
distinta de cero, lo puede hacer de la siguiente forma

genr signx1 = x1 - median(x1)
runs signx1

#
scatters
@Graphs
Uso:            scatters vary ; listavarx    o
		scatters listavary ; varx
Ejemplos:       scatters 1 ; 2 3 4 5
                scatters 1 2 3 4 5 6 ; time

Representa (mediante gnuplot) gráficos bivariantes (scatters) por parejas de 
variables, de vary con respecto a todas las variables de listavarx, o de 
todas las variables de listavary con respecto a varx. En el primer ejemplo de
arriba se sitúa la variable 1 en el eje y y se dibujan cuatro gráficos, el 
primero con la variable 2 en el eje x, el segundo con la variable 3 en el 
eje x, y así sucesivamente. Revisar un conjunto de gráficos como éste puede 
ser interesante al hacer un análisis de datos exploratorio. El número máximo 
de gráficos es seis; cualquier variable extra en la lista será ignorada.


#
seed
@Programming
Uso:            seed entero

Establece la semilla para el generador de números pseudo-aleatorios para
la distribución uniforme y normal (ver la instrucción 'genr'). Por defecto
la semilla se establece cuando se inicia el programa, dependiendo de la hora
del sistema. Si se desean obtener secuencias de números pseudo-aleatorios
repetibles será necesario establecer la semilla de forma manual.


#
setobs
@Dataset
Uso:            setobs periodicidad primobs
Ejemplos:       setobs 4 1990.1
                setobs 12 1978.03
                setobs 20 1.01

Utilice esta orden para forzar al programa a interpretar el conjunto de datos 
actual como de series temporales o de panel, cuando los datos se han leído 
inicialmente como series simples sin fecha. La "periodicidad" debe ser un 
número entero; "primobs" es una cadena que representa la fecha o 
identificación de panel de la primera observación. Utilice un dígito después 
del punto en "primobs" para los datos con periodicidad menor que 10, dos 
dígitos (con un primer cero si es necesario) para periodicidad entre 10 y 99.

En caso de utilizar datos diarios se requiere una forma especial de la cadena
"primobs", concretamente la fecha ha de seguir el patrón YY/MM/DD, por ejemplo
si los datos comienzan el 15 de febrero de 1955 sería "55/02/15" . Si la 
parte YY es menor que 50 se supone que el año pertenece al siglo XXI, en caso
contrario se supone que está en el siglo XX. (Con datos diarios, se aceptan
las dos periodicidades, 5 y 7)


#
setmiss
@Dataset
Uso:           setmiss -1
               setmiss 100 varx

Esta orden se utiliza para hacer que el programa interprete un valor numérico
específico (el primer parámetro de la instrucción) como código de "valor
perdido" al usar datos importados. Si este valor es el único parámetro, como
en el primer ejemplo de arriba, la interpretación se considerará para todas
las series del conjunto de datos. Si se usa, como segundo parámetro, una lista
de variables (por nombre o número) la interpretación se restringe a las
variables especificadas. Así, en el segundo ejemplo se interpreta "100" como
"valor perdido" pero sólo para la variable denominada varx.


#
shell
@Utilities
Uso:		! [shell command]

Un "!" al comienzo de la línea de instrucciones gretl se interpreta como 
una salida al "shell" del usuario. Así se pueden ejecutar instrucciones
del shell desde dentro de gretl.

#
sim
@Dataset
Uso:            sim primobs ultobs y a0 a1 a2 ...

Simula valores para y para los periodos desde "primobs" hasta "ultobs". La 
variable y debe haber sido definida antes con los valores iniciales 
apropiados. "primobs" y "ultobs" deben ser consistentes con la periodicidad.
La fórmula que se usa es:

     y(t) = a0(t) + a1(t)*y(t-1) + a2(t)*y(t-2) + ...

ai(t) pueden ser constantes o los nombres de variables definidas previamente.

Ejemplos:

  sim 1979.2 1983.1 y 0 0.9  [genera y(t) = 0.9*y(t-1)]
  sim 15 25 y 10 0.8 x       [genera y(t) = 10 + 0.8*y(t-1) + x(t)*y(t-2)]

#
smpl
@Dataset
Uso:           smpl primobs ultobs
               smpl -o var_ficticia
               smpl -o
               smpl -r <expresión booleana>
	       smpl full

Restablece el rango muestral. En la primera forma, "primobs" y "ultobs"
deben ser consistentes con la periodicidad de los datos. En la segunda forma,
"var_ficticia" debe ser una variable indicador con valores 0 ó 1: la
muestra se restringirá a aquellas observaciones en las que el valor
indicador sea 1. La tercera forma, smpl -o, quita todas las observaciones 
para las cuales los valores de una o más variables estén 'perdidos'.
La cuarta forma, usando la opción -r, restringe la muestra a los casos 
que satisfagan la condición dada. La última forma, "smpl full",
recupera el rango completo de los datos.



    smpl 3 10                para datos con periodicidad 1
    smpl 1950 1990           para datos anuales con periodicidad 1
    smpl 1960.2 1982.4       para datos trimestrales
    smpl 1960.04 1985.10     para datos mensuales
    smpl 1960.2 ;            para dejar la observación final sin cambiar
    smpl ; 1984.3            para dejar la observación inicial sin cambiar
    smpl -o dum1             para crear una muestra basada en "dum1"
    smpl -r sqft>1400        para restringir la muestra a los casos en los que 
                             la variable sqft tenga un valor mayor que 1400


Hay que señalar un punto especial sobre las formas "-o" y "-r" de smpl:
cualquier información "estructural" en el fichero de cabecera de datos (que
concierna a la naturaleza de series temporales o de panel de los datos) se
pierde al ejecutar esta orden. Se puede reimponer de nuevo la estructura
utilizando la instrucción "setobs".


#
spearman
@Statistics
Uso:            spearman x y
                spearman x y -o

Calcula el coeficiente de correlación por rangos de Spearman para las dos 
variables x e y. No es necesario ordenar las variables y asignar los 
rangos manualmente; la instrucción ya tiene en cuenta esto. Si se 
proporciona la opción -o, se muestran los datos originales junto a los
ordenados.

La ordenación automática se hace de mayor a menor (es decir, al dato mayor se 
le asigna rango 1). Si Vd necesita invertir ese orden, puede crear una nueva 
variable cuyos valores sean los de la variable original cambiados de signo.
Por ejemplo:

  genr altx = -x
  spearman altx y

#
square
@Transformations
Uso:          square x y       o     square -o x y

Genera nuevas variables que son los cuadrados y productos cruzados de las 
variables seleccionadas (-o crea los productos cruzados). En el ejemplo de
arriba las variables creadas serán sq_x = x^2, sq_y = y^2 y x_y = x * y. 
Si una de las variables es una variable ficticia no se tomará su cuadrado,
ya que se obtendría la misma variable.


#
store
@Dataset
Usos:            store nombre_fichero opción
                 store nombre_fichero opción lista_var
Ejemplos:        store misdatos.gdt
                 store misdatos.csv -c
                 store misdatosbin.gdt -o 2 3 4

"nombre_fichero" es el nombre del fichero en el que se guardarán las 
variables. Si no se proporciona "lista_var" se guardarán los valores de 
todas las variables, en caso contrario sólo se grabarán al fichero las
variables especificadas.

Los valores posibles de "opción" son:

ninguno: los datos se guardan en formato xml
  -z  : como el anterior, pero usando compresión de datos tipo gzip
  -o  : los datos se guardan como binarios, en doble precisión
  -s  : los datos se guardan como binarios, en precisión simple
  -c  : los datos se guardan en formato 'valores separados por comas' (CSV), 
        que pueden leerse directamente mediante cualquier programa de 
	hoja de cálculo.
  -r  : los datos se guardan en el formato nativo de GNU R. Así se pueden
        cargar utilizando la instrucción de R 'source()'.
  -m  : los datos se guardan en el formato nativo de GNU Octave. La primera
  	variable citada se toma como variable dependiente y se escribe como
	vector columna; los datos restantes se escriben como una matriz,
	denominada 'X', con una variable por columna.

  -t  : los datos se guardan en el formato "tradicional" de gretl, como en el
  	programa ESL de Ramanathan, con un fichero de datos ascii y un fichero
	de "cabecera" (header .hdr).

#
summary
@Statistics
summary          muestra los estadísticos principales para todas las variables
summary 3 7 9    muestra los estadísticos principales para las variables 
                 número 3, 7, y 9
summary x y z    muestra los estadísticos principales para las variables 
                 x, y y z

Como resultado se ofrecen los siguientes estadísticos: media, desviación
típica (dt), coeficiente de variación (CV= CURTOSIS), mediana, mínimo,
máximo, coeficiente de asimetría y exceso de curtosis.


#
tabprint
@Printing
Uso:            tabprint
                tabprint -o

Debe ejecutarse después de la estimación de una modelo por medio de 'ols'. 
Copia el modelo estimado en forma de entorno tabular de LaTeX, a un fichero
con nombre "model_N.tex", donde N es el número de modelos estimados hasta el
momento en la sesión actual. Esto puede incorporarse en un documento LaTeX.
Ver también la orden 'eqnprint'.

Si se proporciona la opción -o el fichero que se guarda es un documento 
completo LaTeX listo para ser procesado; en caso contrario debe incluirse 
dentro de un documento.


#
testuhat
@Tests
Uso:          testuhat

Debe seguir a una orden de estimación de modelos. Da la distribución de 
frecuencias de los residuos del modelo y un contraste Chi-cuadrado de
normalidad.

#
tsls
@Estimation
Uso:            tsls vardep listavar1 ; listavar2       [-o es opcional]
Ejemplo:        tsls y1 0 y2 y3 x1 x2 ; 0 x1 x2 x3 x4 x5 x6

Calcula las estimaciones de los parámetros de mínimos cuadrados en dos
etapas (MC2E). "vardep" es la variable dependiente, "listavar1" es la lista de
variables independientes (incluyendo variables endógenas del lado derecho
de la ecuación) en la ecuación estructural para las cuales se necesitan las
estimaciones MC2E. "listavar2" es la lista combinada de variables exógenas y
predeterminadas en todas las ecuaciones. Si "listavar2" no es al menos tan
larga como "listavar1", el modelo no está identificado. La opción -o mostrará
la matriz de covarianzas de los coeficientes. En el ejemplo de arriba, las ys
son las variables endógenas y las xs son las variables exógenas y
predeterminadas.


#
var
@Estimation
Uso:            var orden vardep varindep
Ejemplos:       var 4 x1 const time x2 x3
                var 3 1 0 2 3 4

Organiza y estima (vía MCO) una autorregresión vectorial. El primer 
argumento especifica el orden del retardo, después se proporciona la 
estructura para la primera ecuación, de igual forma que en la instrucción 
'ols'. No hay que incluir retardos entre los elementos de la lista 
"varindep" -- se añadirán automáticamente. Se ejecutará una regresión 
para cada variable de la lista, excluyendo la constante, la tendencia 
temporal y las posibles variables ficticias. Los resultados de cada 
ecuación incluyen los contrastes F para restricciones cero de todos
los retardos de cada variable y un contraste F para el máximo retardo.


#
varlist
@Dataset
Uso:          varlist

Muestra una lista de las variables definidas actualmente. "list" y "ls"
son sinónimos.


#
vartest
@Tests
Uso:          vartest x1 x2

Calcula el estadístico F para la hipótesis nula de que las varianzas 
poblacionales de las variables x1 y x2 son iguales y muestra su valor p.


#
wls
@Estimation
Uso:          wls varpesos vardep varindep            [-o opcional]

Se calculan los estimadores de mínimos cuadrados ponderados siendo
"varpesos" la variable de ponderaciones, "vardep" la variable
dependiente y "varindep" la lista de variables independientes. Más
concretamente, se ejecuta una regresión MCO de varpesos*vardep
con respecto a varpesos*varindep. 

Si la variable de ponderaciones es una variable ficticia, esto es 
equivalente a eliminar todas las observaciones que tengan valor cero 
para "varpesos".

Con la opción -o se mostrará la matriz de covarianzas de los 
coeficientes. Se pueden recuperar algunas variables internas 
utilizando la instrucción 'genr'. Para ello es necesario utilizar 
la orden 'genr' inmediatamente después de esta instrucción. Escriba 
"help genr" para ver más detalles sobre esto.

  

