#
add
@Contrastes
(Añadir)
Añade variables a un modelo y contrasta su significatividad.

Se añaden las variables seleccionadas al modelo anterior y se
estima el nuevo modelo. Si se añade más de una variable, se presenta el
estadístico F para el contraste conjunto de significación de las
variables añadidas (sólo para el método MCO) junto con su valor p. Un
valor p inferior a 0.05 indica que los coeficientes son conjuntamente
significativos al nivel del 5 por ciento.

#
adf
@Contrastes
Contraste aumentado de Dickey-Fuller

Esta orden requiere un orden de retardos entero

Calcula los estadísticos para dos contrastes de Dickey-Fuller. En ambos
casos la hipótesis nula es que la variable en cuestión presenta una raíz
unitaria.

El primero es un contraste t basado en el modelo

(1 - L)x(t) = m + g * x(t-1) + e(t).

La hipótesis nula es que g = 0.

En el segundo contraste (aumentado) se estima una regresión no
restringida (teniendo como regresores una constante, una tendencia
temporal, el primer retardo de la variable y "orden" retardos de la
primera diferencia) y una regresión restringida (quitando la tendencia
temporal y el primer retardo). El estadístico de contraste es F
calculado como

[(SCRr - SCRnr)/2]/[SCRnr/(T - k)]

donde T es el tamaño muestral y k el número de parámetros del modelo no
restringido. Es necesario tener en cuenta que los valores críticos para
estos estadísticos no son los usuales.


#
ar
@Estimación
Estimación generalizada (autorregresiva) de Cochrane-Orcutt

Calcula las estimaciones de un modelo utilizando el procedimiento
iterativo generalizado de Cochrane-Orcutt. Las iteraciones se terminan
cuando las sucesivas sumas de cuadrados residuales no varían más del
0.005 por ciento o cuando han transcurrido 20 iteraciones.

En la 'lista de retardos AR' se especifica la estructura del proceso de
error. Por ejemplo, la entrada "1 3 4" corresponde a

     u(t) = rho1*u(t-1) + rho3*u(t-3) + rho4*u(t-4) + et


#
arch
@Contrastes
Realiza un contraste de ARCH (Autoregressive Conditional
Heteroskedasticity)

Esta instrucción requiere un orden de retardos entero.

Contrasta la existencia de un ARCH del orden especificado, en el
modelo. Si el estadístico LM tiene un valor p inferior a 0.10, entonces
se realiza la estimación del ARCH. Si, en la regresión auxiliar,  la
varianza estimada de cualquier observación no es positiva, entonces se
usa el correspondiente residuo al cuadrado. Después se estima el modelo
original utilizando mínimos cuadrados ponderados.



#
boxplots
@Gráficos
(Gráficos de caja)
Análisis exploratorio de los datos.

Estos gráficos (debidos a Tukey y Chambers) presentan la distribución de
una variable. La caja central incluye el 50 por ciento de los datos que
están en el medio de la distribución, es decir, está delimitada por el
primer y el tercer cuartiles. Las "patillas" se extienden desde esta
caja hasta los valores mínimo y máximo. Se dibuja una línea a lo largo de
la caja en la situación de la mediana.

En el caso de cajas recortadas, el recorte muestra los límites de un
intervalo de confianza aproximada del 90 por ciento para la mediana.
Este se obtiene mediante el método bootstrap, que puede tardar un rato
si la serie es muy larga.

Haciendo "click" con el ratón en ventana del gráfico de caja aparece un
menú que permite guardar el gráfico en formato postscript encapsulado
(EPS) o como fichero postscript de página completa. Bajo el sistema X
windows (linux) también se puede guardar el gráfico en formato XPM; bajo
MS Windows se puede copiar al porta-papeles como mapa de bits.

El menú también ofrece la posibilidad de abrir una ventana de texto que
muestra un "resumen de 5 números" (mínimo, primer cuartil, mediana,
tercer cuartil, máximo) y, si se ha elegido la opción de gráfico
"recortado", un intervalo de confianza para la mediana.

Se pueden controlar algunos detalles de los "gráficos de caja" de gretl
por medio de un fichero (de texto plano) de nombre .boxplotrc que se
busca sucesivamente en el directorio de trabajo actual, el directorio
home del usuario (que corresponde a la variable de entorno HOME) y el
directorio gretl del usuario (que se puede observar y cambiar en
Archivo, Preferencias, Menú general). Las opciones que pueden cambiarse
de esta manera son: la fuente a usar al producir una salida postscript
(debe ser un nombre de fuente genérica postscript válido; por defecto es
Helvética), el tamaño de la fuente en puntos (también para la salida
postscript; por defecto es 12), el mínimo y el máximo para el rango
del eje y, la anchura y altura del gráfico en pixels (por defecto, 560 x
448), si se deben imprimir los valores numéricos para los cuartiles y la
mediana (por defecto, no imprimirlos), y si los 'outliers' (puntos que
quedan más allá de 1.5 veces el rango intercuartílico desde la caja
central) deberían indicarse de forma separada (por defecto, no).
Por ejemplo:

font = Times-Roman
fontsize = 16
max = 4.0
min = 0
width = 400
height = 448
numbers = %3.2f
outliers = true

En la segunda línea del final, el valor asociado a "numbers" es una
cadena de formato "printf" como en el lenguaje de programación C;
cuando se especifica, esto controla la inscripción de la mediana y los
cuartiles cerca del gráfico de caja, si no se da una opción "numbers",
esos valores no se escriben. En el ejemplo, los valores se escribirán
a tres dígitos, con dos dígitos después del punto decimal.

No es necesario establecer todas las opciones, y el orden no importa.

Las líneas que no sigan el modelo "clave=valor" son ignoradas, como las
líneas que comiencen con la almohadilla, #.

Después de cada variable especificada en la orden 'gráfico de
caja' puede añadirse una expresión booleana entre paréntesis, para
delimitar la muestra para la variable en cuestión. Debe insertarse un
espacio entre el nombre o número de la variable y la expresión.
Supongamos que tenemos datos de 'salarios' para hombres y mujeres y
tenemos una variable ficticia 'GÉNERO' que toma valor 1 para los hombres
y 0 para las mujeres. En este caso, podemos representar gráficos de caja
comparativos, con las siguientes líneas en el diálogo de la orden
'gráfico de caja':

  salarios (GÉNERO=1) salarios (GÉNERO=0)

#
chow
@Contrastes
Contraste de Chow de homogeneidad estructural

Esta orden necesita un número de observaciones (o fechas, si los datos
tienen fecha).

Debe ejecutarse después de una regresión MCO. Crea una variable ficticia


que es igual a 1 desde el punto de ruptura que se especifique hasta el
final de la muestra y 0 en el resto. También crea términos de
interacción entre esta variable ficticia y las variables independientes
originales. Se ejecuta una regresión aumentada incluyendo esos términos
y se calcula un estadístico F tomando la regresión aumentada como 'no
restringida' y la original como 'restringida'. Este estadístico es
adecuado para contrastar la hipótesis nula de que no hay cambio
estructural en el citado punto de ruptura.


#
coint
@Contrastes
Contraste de cointegración

Con esta orden (que necesita un orden de retardos entero) se realizan
los contrastes de Dickey-Fuller de la hipótesis nula de que cada una de
las variables seleccionadas tiene una raíz unitaria, considerando el
orden de retardos dado. Se estima la regresión cointegrante y se realiza


un contraste ADF sobre los residuos de la regresión.

También se ofrece el estadístico de Durbin-Watson de la regresión
cointegrante.

(Hay que señalar que para ninguno de estos estadísticos pueden
utilizarse las tablas estadísticas usuales)


#
compact
@Conjunto de datos
(Compactar)
Escribiendo datos a una frecuencia inferior.

Cuando se añade a un conjunto de datos una serie que es de una
frecuencia superior, es necesario "compactar" esa nueva serie. Por
ejemplo, una serie mensual tendrá que ser "compactada" para introducirla


en un conjunto de datos trimestral. Se ofrecen tres opciones para el
"compactado":

1. Promediado: el valor escrito en el conjunto de datos será la media
aritmética de los valores de la serie en cuestión. Por ejemplo, el valor


introducido para el primer trimestre de 1990 será la media de los
valores de enero, febrero y marzo de 1990.

2. Valores de 'final de periodo': el valor escrito en el conjunto de
datos es el último valor del periodo correspondiente en los datos
de más alta frecuencia. Por ejemplo, en el primer trimestre de 1990 se
introduciría el valor de marzo de 1990.

3. Valores de 'principio de periodo': el valor escrito en el conjunto de
datos es el primer valor del periodo correspondiente en los datos de
más alta frecuencia. Por ejemplo, en el primer trimestre de 1990 se
introduciría el valor de enero de 1990.



#
corc
@Estimación
Modelo de Cochrane-Orcutt

Esta orden calcula las estimaciones de un modelo usando el procedimiento
iterativo de Cochrane-Orcutt. Las iteraciones acaban cuando dos valores
sucesivos de rho no difieren en más de 0.001 o cuando se han realizado
ya 20 iteraciones.

La regresión transformada final se estima para el rango de observación
primobs+1 ultobs que esté actualmente en efecto.


#
dialog box for models
@Estimación
(cuadro de diálogo para los modelos)
Para elegir la variable dependiente, seleccione una variable en la lista


de la izquierda y presione el botón "Elegir->" que señala a la caja
correspondiente a "variable dependiente". Si se activa el cuadro
"Establecer por defecto" la variable que se elija será preseleccionada
como variable dependiente la próxima vez que se abra el cuadro de
diálogo de modelos. Atajo: haciendo doble "click" sobre una variable
del lado izquierdo se selecciona y se establece como variable
dependiente por defecto.

Para elegir las variables independientes, selecciónelas a la izquierda y


presione el botón "Añadir->" (o haga "click" con el botón derecho del
ratón). Se pueden seleccionar varias variables contiguas arrastrándolas
con el ratón. Se puede seleccionar un grupo de variables no
contiguas haciendo click sobre ellas mientras se pulsa la tecla Ctrl.


#
diff
@Transformaciones
(Diferencia regular)

Se calcula la primera diferencia de cada variable de la lista
dada y el resultado se guarda en una nueva variable con prefijo "d_".
Así por ejemplo, la nueva variable
d_x = x(t) - x(t-1).

#
export
@Conjunto de datos
Exporta datos de gretl a otros formatos.

Se pueden exportar datos en formato CSV (Valores separados por comas):
estos datos pueden abrirse desde hojas de cálculo y muchos otros
programas.

También se pueden exportar datos a los formatos nativos de GNU R y GNU
Octave. Para más información sobre estos programas (ambos desarrollan
análisis estadísticos avanzados) por favor visite sus propias
páginas web, http://www.r-project.org/ y http://www.octave.org/


#
factorized plot
@Gráficos
(Gráfico con factor de separación)

Esta orden necesita que Vd elija tres variables y la última de ellas
ha de ser una variable ficticia (con valor 1 ó 0). La variable Y se
representa con respecto a la variable X, con los puntos coloreados de
forma diferente dependiendo del valor de la tercera variable.

Por ejemplo: si vd tiene datos de salarios y nivel de educación para una
muestra de varias personas; y Vd también tiene una variable ficticia con
valor 1 para los hombres y 0 para las mujeres (como en el fichero que se
suministra con gretl, data7-2). Un "gráfico factorizado" de SALARIOS
con respecto a EDUCACIÓN usando la variable SEXO como factor mostrará
los puntos para los hombres en un color y para las mujeres en otro (con
una leyenda para identificarlos).


#
genr
@Transformaciones
Genera una nueva variable
Uso:             nuevo_nombre_variable = transformación

Crea nuevas variables, normalmente por medio de transformaciones de
variables ya existentes. Ver también diff, logs, lags, ldiff, multiply
y square como atajos.

Los operadores matemáticos que se soportan son, en orden de precedencia:
^ (exponenciación); *, / y % (módulo o resto); + y -.

Los operadores booleanos (de nuevo en orden de precedencia) son ! (NO
lógico), & (Y lógico), | (O lógico), >, <, = y != (No igual). Los
operadores booleanos pueden usarse al construir variables ficticias: por
ejemplo (x>10) devuelve 1 si x(t)>10 y en caso contrario 0.

Las funciones que se soportan pertenecen a estos grupos:

- Funciones matemáticas standard: abs, cos, exp, int (parte entera), ln
(logaritmo natural: log es un sinónimo), sin (seno), sqrt (raíz
cuadrada).

- Funciones estadísticas: mean (media aritmética), median (mediana), var
(varianza), sd (desviación típica o estandard), sum, cov (covarianza),
corr (coeficiente de correlación, min (mínimo), max (máximo).

- Funciones de series temporales: lag (retardo), lead (adelanto), diff
(primera diferencia), ldiff (log-diferencia, o primera diferencia del
logaritmo natural).

- Misceláneas: cum (acumulación), sort (ordenación), uniform
(distribución uniforme), normal (distribución normal), missing (devuelve
1 si la variable tiene la observación perdida, en caso contrario 0),
misszero (reemplaza el código de observación perdida por un 0), zeromiss
(operación inversa de misszero).

Todas las funciones anteriores, a excepción de cov, corr, uniform y
normal, toman como único argumento o el nombre de una variable (hay
que notar que, en una orden genr,  no es posible referirse a las
variables utilizando su número de ID) o una expresión compuesta que se
evalúa en una variable (p.ej. ln((x1+x2)/2)). cov y corr requieren dos
argumentos (dos variables) y devuelven, respectivamente, la covarianza y
el coeficiente de correlación entre las dos variables mencionadas.
uniform() y normal() no tienen argumentos y devuelven, respectivamente,
series pseudoaleatorias obtenidas a partir de las distribuciones
uniforme (0-100) y normal standard (ver también la instrucción seed).

Hay varias variables que se definen internamente al ejecutar una
regresión, que pueden usarse también en transformaciones, como son:

  $ess         suma de cuadrados de los residuos
  $rsq         R-cuadrado no corregido
  $T           número de observaciones utilizado por el modelo
  $df          grados de libertad
  $trsq        TR^2 (T veces el R-cuadrado, siendo T el tamaño
               muestral)
  $sigma       desviación típica de los residuos
  $lnl         log-verosimilitud (en modelos logit y probit)
  coeff(var)   coeficiente estimado de var
  stderr(var)  desviación típica estimada del estimador de var
  rho(i)       coeficiente autorregresivo de iésimo orden de los
               residuos
  vcv(xi,xj)   covarianza entre los coeficientes de las variables
               xi y xj

La variable interna $nobs contiene el número de observaciones del
dominio muestral actual, que puede ser o puede no ser igual al $T del
último modelo.

La variable interna $pd contiene la periodicidad o frecuencia de los datos (por 
ejemplo, 4 para datos trimestrales, 12 para mensuales).

La variable interna t hace referencia a las observaciones, comenzando en
1. Así, es posible hacer "genr dum15 = (t=15)" para generar una variable
ficticia con valor 1 para la observación 15 y 0 para el resto.

Ejemplos de fórmulas válidas:

   y = x1^3          [x1 al cubo]
   y=ln((x1+x2)/x3)  [argumento compuesto a función ln]
   z=x>y             [hace z(t) igual a 1 si x(t) > y(t) en caso
                      contrario 0]
   y=x(-2)           [x retardada 2 periodos]
   y=x(2)            [x adelantada 2 periodos]
   y = mean(x)       [media aritmética]
   y = diff(x)       [y(t) = x(t) - x(t-1)]
   y = ldiff(x)      [y = ln(x(t)) - ln(x(t-1))]
                      ldiff(x) es la tasa instantánea de
                      crecimiento de x.
   y = sort(x)       [ordena x en orden creciente y la guarda en y]
   y = -sort(-x)     [ordena x en orden decreciente]
   y = int(x)        [trunca x y guarda su valor entero como y]
   y = abs(x)        [guarda los valores absolutos de x]
   y = sum(x)        [suma los valores de x excluyendo las entradas
                      de los valores perdidos -999]
   y = cum(x)        [acumula x: y(t) es la suma de x hasta t]
   aa = $ess         [aa = suma de cuadrados de los residuos de la
                      última regresión]
   x = coeff(sqft)   [recoge el coeficiente de sqft del
                      último modelo]
   rho4 = rho(4)     [recoge el coeficiente autorregresivo de
                      4º orden del último modelo (se supone un modelo
                      ar)]
   cv=vcv(x1, x2)    [covarianza de los coeficientes de x1 y x2
                      en el último modelo]
   x=uniform()/100   [variable pseudoaleatoria uniforme, rango 0 a 1]
   x=3*normal()      [variable pseudoaleatoria normal, con media 0
                      y desviación típica 3]

Sugerencias sobre variables ficticias:

* Supongamos que x se codifica con valores 1, 2 o 3 y Vd desea tres
variables ficticias, d1 si x=1, 0 en caso contrario, d2=1 si x=2, 0 en
caso contrario y así sucesivamente. Para crearlas, use las fórmulas
d1 = (x=1), d2 = (x=2), y d3 = (x=3).
* Para obtener z = máx(x,y) genere d=x>y y después
z=(x*d)+(y*(1-d))

#
graphing
@Gráficos
(Gráficos)
generando gráficos de varios tipos

Gretl llama a un programa aparte, gnuplot, para generar los gráficos.
Gnuplot es un programa gráfico de múltiples características con
miles de opciones. Gretl le proporciona a Vd acceso directo, vía una
interface gráfica, a sólo un pequeño subconjunto de esas opciones e
intenta elegir para vd los valores adecuados; también permite que vd
tome completamente el control sobre los detalles del gráfico si así lo
desea.

Bajo MS Windows vd puede hacer click en la esquina de arriba a la
izquierda de la ventana del gráfico, obteniendo así un menú porta-papeles
le permite elegir varias cosas (incluyendo copiar el gráfico al
porta-papeles de Windows y enviarlo a la impresora).

Para tener un control completo sobre el gráfico, siga este
procedimiento:

- Cierre la ventana del gráfico.
- Desde el menú de sesión, seleccione "Añadir el último gráfico".
- En la ventana de iconos de sesión, haga click con el botón
derecho del ratón sobre el icono del nuevo gráfico y seleccione o
"Editar usando GUI" o "Editar las órdenes de gráfico". La entrada de
" Editar usando GUI" abre un controlador gráfico para gnuplot que le
permite refinar varios aspectos del gráfico. La entrada de "Editar las
órdenes de gráfico" abre una ventana de editor que contiene el fichero
actual de instrucciones de Gnuplot para generar el gráfico: esto le
proporciona a vd un control completo sobre los detalles del gráfico --si
vd conoce algo sobre gnuplot. Para más detalles,ver
http://ricardo.ecn.wfu.edu/gnuplot.html or www.gnuplot.org.

#
hccm
@Estimación
(mcch)
Matriz de covarianzas consistente ante heterocedasticidad

Esta orden ejecuta una regresión donde los coeficientes se estiman
por medio del procedimiento MCO standard, pero las desviaciones típicas
de los estimadores de los coeficientes se calculan de una manera que es
robusta ante la heterocedasticidad. Concretamente, se usa el
procedimiento "jacknife" de MacKinnon-White


#
hilu
@Estimación
Método de Hildreth-Lu

Calcula las estimaciones de un modelo utilizando el procedimiento de
búsqueda de Hildreth-Lu (refinado mediante el método de CORC
[Cochrane-Orcutt]). Se representa la suma de cuadrados de los residuos
del modelo transformado con respecto a los valores de rho
desde -0.99 hasta 0.99. La regresión final transformada  se calcula para
el rango de observación primobs+1 ultobs que esté actualmente en efecto.


#
hsk
@Estimación
(corrección de heterocedasticidad)
Estimaciones corregidas de Heterocedasticidad

Se ejecuta una regresión MCO y se guardan los residuos. El logaritmo
del cuadrado de dichos residuos se constituye como variable
dependiente en una regresión auxiliar, en cuyo lado derecho de la
ecuación están las variables independientes originales y sus cuadrados.
Los valores ajustados en la regresión auxiliar se usan entonces para
construir una serie de ponderaciones y el modelo original se reestima
utilizando mínimos cuadrados ponderados. Este resultado final es el que
aparece en el cuadro de resultados.

La serie de ponderaciones se forma como 1/sqrt(exp(fit)), donde "fit"
representa a los valores ajustados obtenidos de la regresión auxiliar.

#
lad
@Estimación
(Estimador de mínima desviación absoluta)
Uso:          lad vardep varindeps

Calcula una regresión que minimiza la suma de las desviaciones absolutas
entre los valores observados y los valores ajustados de la variable dependiente. 
Las estimaciones de los coeficientes se obtienen utilizando el algoritmo simplex 
de Barrodale-Roberts; se muestra un aviso si la solución no es única. Las 
desviaciones típicas se obtienen utilizando un método 'bootstrap' con 500
iteraciones.

#
lags
@Transformaciones
(retardos)

Crea nuevas variables que son valores retardados de cada una de las
variables de la lista que se suministra. El número de contrapartes
retardadas para cada una de las variables listadas es igual a la
periodicidad de los datos. Por ejemplo, si la periodicidad es 4 (datos
trimestrales), se crearán cuatro términos retardados; si en la lista que
se ha suministrado está la variable "x", la orden crea x_1 = x(t-1),
x_2 = x(t-2), x_3 = x(t-3) y x_4 = x(t-4).


#
ldiff
@Transformaciones

Se obtiene la primera diferencia del logaritmo natural de cada variable
de la lista suministrada y el resultado se guarda en una nueva variable
con el prefijo "ld_".

Así por ejemplo, la nueva variable ld_x = ln[x(t)] - ln[x(t-1)].

#
logit
@Estimación
Regresión Logit

La variable dependiente debería ser una variable binaria. Se utiliza
el método de mínimos cuadrados iterativos (el método EM o de
expectativa-maximización) para obtener las estimaciones
máximo-verosímiles de los coeficientes de las variables independientes.
Como el modelo es no lineal, las pendientes dependen de los valores de
las variables independientes: las pendientes que gretl muestra se
evalúan en las medias de dichas variables. El estadístico Chi-cuadrado
contrasta la hipótesis nula de que todos los coeficientes, excepto la
constante, son cero.


#
logs
@Transformaciones

Se calcula el logaritmo natural de cada una de las variables de la lista


que se suministra y el resultado se guarda en una nueva variable con 
prefijo "l_". Así por ejemplo la nueva variable l_x = ln(x).


#
loop
@Programación
(Bucle)
instrucciones repetidas

Uso:            loop número_de_veces
                loop while condición
		loop for i=principio..final
Ejemplos:       loop 1000
		loop while essdiff > .00001
		loop for i=1991..2000

Esta instrucción (de guión) abre un modo especial en el cual el programa
acepta instrucciones a repetir o un número de veces específico, o
mientras una condición se satisfaga, o para los sucesivos valores
enteros de una variable índice (interna) i.

Dentro de un bucle sólo se pueden utilizar 7 instrucciones:
genr, ols, print, sim, smpl, store y summary (store no puede usarse en
un bucle "while"). Con genr y ols es posible hacer bastantes cosas.
Se sale de este modo especial de introducir instrucciones de bucle
mediante la orden "endloop": en este momento se ejecutan las órdenes
que estén en la "pila".

Los bucles no pueden estar anidados. La instrucción ols muestra un
resultado especial dependiendo del tipo de bucle. Si se especifica un
"numero_de_veces" no se muestran los resultados de cada regresión
individual, sino que se obtiene una salida con (a) la media de
cada coeficiente estimado a lo largo de todas las repeticiones, (b) las
desviaciones típicas de estos coeficientes estimados, (c) la media
de las desviaciones típicas estimadas para cada coeficiente y (d) la
desviación típica de las desviaciones típicas estimadas. Esto sólo
tiene sentido si hay algún input aleatorio en cada paso. La instrucción
está diseñada para el análisis de Monte Carlo.

Si se da una condición "while", se muestran los resultados del modelo
especificado a partir de la última vuelta del bucle: esto está diseñado
para mínimos cuadrados iterativos.

La instrucción "print" también se comporta de forma diferente en el
contexto de un bucle con "número_de_veces". En concreto esta instrucción
muestra la media y la desviación típica de la variable a lo largo de las
repeticiones del bucle. Esto está diseñado para variables que toman un
solo valor en cada iteración, por ejemplo la scr (suma de cuadrados
residual $ess ) de una regresión. La instrucción "print" se comporta de
la forma usual con las otras construcciones de bucle.

La instrucción "store" (se usa sólo una de ellas por bucle y sólo en un
bucle con "número_de_veces") escribe los valores de las variables
especificadas desde cada iteración del bucle a un fichero. De esta
manera se guarda una grabación completa de las variables. Este fichero
de datos puede después ser leído y analizado mediante gretl.

Ejemplo de código de bucle (Monte Carlo):

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
lmtest
@Contrastes
Contraste de Multiplicador de Lagrange

Bajo este encabezamiento se encuentran varios contrastes de hipótesis. 
Lo que tienen en común es que el contraste incluye la estimación de una
regresión auxiliar, en la que la variable dependiente es el residuo de 
alguna regresión "original". Entre las variables del lado derecho
se incluyen las de la regresión original y algunas adicionales. El
estadístico de contraste se calcula como (tamaño muestral x
R-cuadrado) de la regresión auxiliar: este se distribuye como una
Chi-cuadrado con grados de libertad iguales al número de variables 
adicionales, bajo la hipótesis nula de que las variables adicionales no 
tienen poder explicativo sobre el residuo. Un valor muy alto de este 
estadístico (valor p pequeño) sugiere que esta hipótesis nula debería 
ser rechazada.


#
markers
@Conjunto de datos
(marcadores)
Añade marcadores de caja al conjunto de datos.

Esta instrucción necesita el nombre del fichero que contenga los
"marcadores de caja", es decir, pequeñas etiquetas que
identifican a las observaciones individuales en el conjunto de datos
(por ejemplo, nombres o códigos de países o de ciudades). Estas
etiquetas no deberían tener más de 8 caracteres. El fichero debería
tener un marcador por línea y debería haber tantos marcadores como 
observaciones en el conjunto de datos. Si se satisfacen estas 
condiciones y se encuentra el fichero especificado, se añadirán los 
"marcadores de caja"; estos se podrán ver cuando Vd elija "Mostrar
valores" en el menú "Datos" de gretl.


#
meantest
@Contrastes
(Contraste de medias)

Calcula el estadístico t para el contraste de la hipótesis nula de que
las medias poblacionales, de las dos variables elegidas, son iguales. 
También muestra su valor p. La instrucción puede ejecutarse con o sin el 
supuesto de que las varianzas de las dos variables son iguales (aunque 
esto supondrá una diferencia en el estadístico de contraste sólo si hay 
un número diferente de valores no-perdidos para las dos variables).


#
missing values
@Conjunto de datos
(Valores perdidos)

Establece un valor numérico que será interpretado como "valor perdido" o
"no disponible", o para una serie de datos particular (bajo el
menú de "Variable") o globalmente para el conjunto de datos completo
(bajo el menú "Muestra").

Gretl tiene su propio código interno para los valores perdidos, pero a
veces los datos importados pueden emplear un código diferente. Por
ejemplo, si una serie determinada está codificada de forma que el valor
-1 indica "no disponible", se puede seleccionar "Establecer código de
'valor perdido'" bajo el menú de "Variable" y escribir el valor "-1"
(sin las comillas). Gretl entonces leerá los -1 como observaciones
perdidas.


#
mpols
@Estimación
(MCO de precisión múltiple)

Calcula las estimaciones de mínimos cuadrados ordinarios utilizando 
operaciones aritméticas con precisión múltiple. Esta orden sólo está 
disponible si gretl está configurado con soporte para GMP, la biblioteca 
de precisión múltiple GNU. Hay que señalar que la precisión de los 
resultados de la regresión puede verse limitada por (a) la precisión de 
los datos que se leen del fichero y (b) cualquier transformación 
realizada usando la instrucción genr, que trabaja utilizando operaciones 
aritméticas ordinarias de punto flotante y doble precisión.

#
nulldata
@Conjunto de datos

Establece un conjunto de datos "vacío", que contiene sólo una constante, 
con periodicidad 1 y el número de observaciones que se especifique. Esto 
puede utilizarse por motivo de simulación: algunas instrucciones genr
(p.ej. genr uniform(), genr normal(), genr time) generarán datos 
artificiales desde cero para rellenar el conjunto de datos. La 
instrucción "nulldata" también puede ser útil combinada con "loop".

#
ols
@Estimación
Método de mínimos cuadrados ordinarios

Calcula las estimaciones de mínimos cuadrados ordinarios de los
coeficientes del modelo especificado. Muestra los valores p para los 
estadísticos t (a dos colas) y F. Un valor p inferior a 0.01 indica 
significatividad al nivel del 1 por ciento. También se muestran una 
serie de estadísticos de selección de modelos.

Ver "/Temas/Estimación/Cuadro de diálogo" para ayuda sobre el uso del 
cuadro de diálogo.

#
omit
@Contrastes
Omite variables de un modelo y contrasta su significatividad conjunta

Las variables elegidas se sustraen del modelo anterior y se estima el 
nuevo modelo. Si se omite más de una variable, se mostrará el 
estadístico F de Wald para las variables omitidas junto con su valor p 
(sólo para el método MCO). Un valor p inferior a 0.05 indica que los 
coeficientes son conjuntamente significativos al nivel del 5 por ciento.


#
online databases
@Conjunto de datos
(Bases de datos en línea)
Acceso a bases de datos vía internet

Gretl puede acceder a las bases de datos del sitio web de gretl, en
Wake Forest University (su ordenador debe de estar conectado a 
internet para que esto funcione).

Bajo el menú "Archivo, Revisar bases de datos" seleccione la entrada
"sobre servidor". Ahora debería aparecer una ventana mostrando las bases 
de datos gretl disponibles en Wake Forest (dependiendo del lugar en que
Vd se encuentre y de su conexión a internet, esto puede tardar
unos segundos). Junto al nombre de la base de datos y a una pequeña 
descripción aparecerá una entrada de "Estado local": esto indica si Vd 
ha instalado la base de datos de forma local (sobre el disco duro de su 
ordenador) y si es así, si se encuentra actualizada con la versión
que actualmente hay en el servidor. Si Vd ha instalado localmente una
determinada base de datos y ésta se encuentra actualizada, no hay 
ninguna ventaja por acceder a ella mediante el servidor. Pero para una 
base de datos que no esté instalada y/o actualizada, Vd puede desear un 
listado de las series de datos: haga "click" sobre "Obtener listado de
series". Esto hace aparecer una nueva ventana desde la cual se pueden 
visualizar los valores de la serie de datos que se elija, representar 
esos valores o importarlos al espacio de trabajo de gretl. Estas tareas 
pueden realizarse usando el menú "Series" o por medio del menú 
contextual que aparece el hacer "click" con el botón derecho del ratón 
sobre una serie dada. También es posible buscar entre el listado una 
variable determinada (con el menú "Buscar").

Si se desea un acceso a los datos más rápido, o un acceso a los datos 
"off line", es posible seleccionar la línea que muestra la base de datos 
que interese, en la ventana inicial de bases de datos, y presionar el 
botón "Instalar". Esto hará que la base de datos se descargue en formato 
comprimido, después se puede descomprimir e instalar en el disco
duro. Más adelante podremos encontrarla bajo el menú "Archivo, Revisar 
bases de datos, Nativa gretl". (Esta característica de gretl depende de
otros proyectos de software de 'fuente abierta': la biblioteca de 
compresión de datos zlib y el programa descargador GNU "wget", de los 
cuales gretl ha tomado prestados algunos trozos de código).

#
panel
@Conjunto de datos
Establece estructura de datos de panel

Las dos opciones disponibles aquí son "series temporales apiladas" y 
"secciones cruzadas apiladas". Si se quiere hacer uso de la
instrucción "MCO combinados" y sus diagnósticos de panel asociados, 
gretl debe saber en qué forma están organizados los datos. "Series 
temporales apiladas" significa que los bloques del fichero de datos son 
series temporales para cada una de las unidades de sección cruzada. Por 
ejemplo, las primeras 10 filas de datos podrían representar los valores 
de ciertas variables para el país A durante 10 periodos, las siguientes 
10 filas los valores para el país B durante los mismos 10 periodos, y 
así sucesivamente. "Secciones cruzadas apiladas" significa que los 
bloques del fichero de datos son secciones cruzadas para cada uno de los 
periodos. Por ejemplo, las primeras 6 filas de datos podrían representar 
los valores de ciertas variables para los países A a F para el año 1970,
las siguientes 6 filas los valores para los mismos países en 1971, y así
sucesivamente.

Si se guarda el fichero de datos después de establecer este atributo, la
información se grabará en el fichero de datos y no será necesario
establecerlo de nuevo la próxima vez que se usen estos datos.


#
pooled
@Estimación
Estimación de MCO combinados

Esta instrucción se utiliza con datos de panel. Para sacar provecho de
ella, se debería especificar un modelo sin ninguna variable ficticia
representando a las unidades de sección cruzada. La rutina presenta
estimaciones de MCO combinados directamente, las cuales tratan las
variaciones de sección cruzada y de series temporales de igual forma.
Este modelo puede ser o puede no ser apropiado. Bajo el menú de
"Contrastes" en la ventana del modelo se puede encontrar una entrada
"diagnósticos de panel", en la cual se contrastan MCO combinados contra
las principales alternativas: los modelos de efectos fijos y de efectos
aleatorios.

En el modelo de efectos fijos se añade una variable ficticia para todas
las unidades de sección cruzada excepto una, permitiendo así al término
constante de la regresión variar a través de las unidades (individuos).
Se presenta un estadístico F para contrastar la significatividad
conjunta de dichas variables ficticias: si el valor p de este contraste
es pequeño, esto es una indicación en contra de la hipótesis nula (de
que el modelo simple combinado es el adecuado) y en favor del modelo de
efectos fijos.

Por otro lado, el modelo de efectos aleatorios descompone la varianza
residual en dos partes, una parte específica de la unidad de sección
cruzada o "grupo" y la otra específica de la observación particular.
(Este estimador sólo puede calcularse si el panel es suficientemente
"ancho", es decir, si el número de unidades de sección cruzada que hay
en el conjunto de datos es superior al número de parámetros a estimar).
El estadístico LM de Breusch-Pagan contrasta la hipótesis nula (de
nuevo, de que el estimador de MCO combinados es el adecuado) contra la
alternativa de efectos aleatorios.

Es muy posible que el modelo de MCO Combinados sea rechazado contra
ambas alternativas (efectos fijos y efectos aleatorios). ¿Cómo se puede
entonces determinar cuál de los dos estimadores es más apropiado? El
contraste de Hausman (que también se presenta, dado que se puede estimar
el modelo de efectos fijos) da una indicación en este sentido. Si
que el error específico de unidad --o grupo-- está incorrelacionado con
las variables independientes, el estimador de efectos aleatorios es más
eficiente que el de efectos fijos; en caso contrario el estimador de
efectos aleatorios es inconsistente y entonces será preferible el
estimador de efectos fijos. La hipótesis nula para el contraste de
Hausman es que el error específico de grupo no está muy correlacionado
(y por tanto, es preferible el estimador de efectos fijos). Entonces,
para este contraste, un valor p pequeño es una indicación en contra del
modelo de efectos aleatorios y a favor del de efectos fijos.

Para un desarrollo riguroso de este tema ver  "Análisis
Econométrico" de Greene (4ª edición), capítulo 14.

#
probit
@Estimación
Regresión Probit

La variable dependiente debería ser una variable binaria. Se utiliza
el método de mínimos cuadrados iterativos (el método EM o de
expectativa-maximización) para obtener las estimaciones
máximo-verosímiles de los coeficientes de las variables independientes.
Como el modelo es no lineal, las pendientes dependen de los valores de
las variables independientes: las pendientes que gretl muestra se
evalúan en las medias de dichas variables. El estadístico Chi-cuadrado
contrasta la hipótesis nula de que todos los coeficientes, excepto la
constante, son cero.

#
range-mean
@Gráficos
Gráfico Rango-Media

Este es un gráfico simple para ayudar a decidir si una serie temporal, y(t), 
tiene o no varianza constante. Se toma la muestra completa t=1,...,T y se divide 
en pequeñas submuestras de tamaño arbitrario k [gretl elige k=sqrt(T)]. La 
primera submuestra se forma con y(1),...,y(k), la segunda con y(k+1),...,y(2k), 
y así sucesivamente. Para cada submuestra se calcula la media muestral y el 
rango (=máximo-mínimo) y se construye un gráfico con las medias en el eje 
horizontal y los rangos en el vertical. De esta forma, cada submuestra está 
representada por un punto en este plano. Si la varianza de la serie fuera 
constante los rangos de las submuestras no debería depender de sus medias; si se 
observa que los puntos se aproximan a una recta con pendiente creciente, esto 
sugiere que la varianza de la serie aumenta cuando la media aumenta; si los 
puntos se aproximan a una recta con pendiente decreciente, esto sugiere que la 
varianza está disminuyendo cuando la media aumenta.

Además del gráfico, gretl presenta una ventana de resultados que muestra las 
medias y los rangos para cada submuestra, el coeficiente estimado para la
pendiente en una regresión MCO de los rangos sobre las medias y el valor p para 
el contraste de la hipótesis nula de que esta pendiente es cero. Si el 
coeficiente de pendiente es significativo al nivel de significación del 10 por 
ciento, en el gráfico se muestra también la recta ajustada en la regresión de 
los rangos sobre las medias.

#
rhodiff
@Transformaciones
Uso:            rhodiff rho listavar
Ejemplo:        rhodiff .65 2 3 4

Crea las correspondientes variable rho-diferenciadas de las variables
(dadas por número o nombre) de listavar y las añade al conjunto de
datos.  Sea la variable v1 de la lista entonces se crea
rd_v1 = v1(t) - rho*v1(t-1).

#
scatters
@Gráficos
Múltiples gráficos cruzados por parejas

Representa un conjunto de gráficos cruzados de la "variable del eje Y"
elegida con respecto a las "variables del eje X" elegidas
consecutivamente. Puede ser útil echar un vistazo a estos gráficos al 
hacer un análisis exploratorio de datos. El número máximo de gráficos es 
seis; cualquier variable extra en el eje X será ignorada.

#
seed
@Programación
Pone en marcha el generador de números aleatorios

Requiere como entrada un entero. Establece la semilla para el generador 
de números pseudoaleatorios utilizado por las opciones "Aleatoria 
uniforme" y "Aleatoria normal" bajo el menú de "Datos, Añadir 
variables". Por defecto, la semilla se establece, utilizando la hora del 
sistema, cuando se inicia el programa. Si se desea obtener secuencias de
números pseudoaleatorios repetibles es necesario establecer la semilla
de forma manual.


#
setobs
@Conjunto de datos
Establece la frecuencia de los datos y la observación inicial

Utilice esta orden para forzar al programa a interpretar el conjunto de 
datos actual como "de series temporales" o "de panel" cuando los datos 
se han leído inicialmente como series simples sin fecha. Se necesitan 
dos parámetros: una frecuencia entera y una observación inicial 
(normalmente una fecha).

Ejemplos de entradas válidas:

  4 1990.1       Interpretar los datos como trimestrales, 
                 comenzando en 1990, trimestre 1.
  12 1978.03     Interpretar los datos como mensuales, comenzando en 
                 marzo de 1978
  20 1.01        Frecuencia de datos 20, comenzando en la observación
                 1.01 (datos de panel)
  5 72/01/10     Datos diarios (5 días por semana), desde 10 de enero
                 de 1972
  7 02/01/10     Datos diarios (7 días por semana), desde 10 de enero
                 de 2002

#
sim
@Conjunto de datos
(Simulación)
Introducir datos simulados en una variable

Esta instrucción requiere una observación inicial, una observación 
final, el nombre de una variable (ya existente en el conjunto de datos) 
en la cual introducir los valores y una lista de coeficientes
autorregresivos, que pueden ser constantes numéricas o nombres 
de variables. Por ejemplo, si en el diálogo de simulación se introduce
    
    1979.2 1983.1 y 0 0.9

esto rellenará y, desde 1979.2 hasta 1983.1 con los valores

    y(t) = 0 + 0.9 y(t-1)

De forma similar

    15 25 y 10 0.8 x

generará desde la observación 15 a la 25:

    y(t) = 10 + 0.8 y(t-1) + x(t) y(t-2)

#
sampling
@Conjunto de datos
(Muestreo)
Seleccionar una submuestra del conjunto de datos actual.

Si se elige "Muestra/Definir a partir de v.ficticia..." es necesario 
suministrar el nombre de una variable ficticia (indicador) que 
debería tener los valores 0 o 1 en cada observación. La muestra se 
restringirá a aquellas observaciones para las cuales la variable 
ficticia tome valor 1. (Haciendo "click" en la línea de una variable en 
la ventana principal de los datos se insertará el nombre de esa variable 
en la caja de diálogo).

Si se elige "Muestra/Restringir a partir de criterio..." es necesario 
proporcionar una expresión booleana (lógica), del mismo tipo que 
se utilizaría para definir una variable ficticia. Por ejemplo, la 
expresión "sqft > 1400" seleccionará sólo los casos para los que la 
variable sqft tenga un valor mayor que 1400. Las condiciones pueden 
estar concatenadas utilizando los operadores lógicos "&" (Y lógico) y
"|" (O lógico).

La entrada de menú "Muestra/Quitar todas las obs. con valores 
perdidos..." redefine la muestra excluyendo todas las observaciones para 
las que los valores de una o más variables están "perdidos" (dejando 
sólo las observaciones completas, en las que hay un valor numérico para 
todas las variables).

Al definir una muestra a partir de una variable ficticia, una expresión 
booleana o por el criterio de los valores perdidos hay que señalar que
cualquier información "estructural" en el fichero de encabezamiento de 
datos (con relación a la naturaleza de 'series temporales' o 'de panel' 
de los datos) se perderá. Es posible reimponer de nuevo la estructura 
mediante "Muestra/Establecer frecuencia, observación inicial..."

Para sólo reorganizar la muestra especificando una observación inicial y 
una final ver "smpl".


#
smpl
@Conjunto de datos
(Establecer el rango muestral)

Restablece el rango muestral especificando una observación inicial y
una observación final (Muestra/Establecer rango...). Este mecanismo se
utiliza para establecer una submuestra de una serie temporal de datos.
Las observaciones inicial y final dadas deberían estar en un formato
consistente con la frecuencia de los datos, p.ej. "1985.1" para datos
trimestrales o ""1996.03" para datos mensuales (marzo de 1996).

#
spearman
@Estadísticos

Calcula el coeficiente de correlación por rangos de Spearman para un par 
de variables especificado. No es necesario antes ordenar y hacer el
ranking de las variables, la función se encarga de ello.

El ranking automático es de mayor a menor (es decir, al dato mayor se le 
asigna rango 1). Si se necesita invertir este ranking, esto se puede 
hacer creando una nueva variable que sea el negativo de la original. Por 
ejemplo:

  genr altx = -x
  spearman altx y

#
square
@Transformaciones
(Cuadrados)

Genera nuevas variables que son los cuadrados de las variables de la
lista dada. Las nuevas variables se nombran con el prefijo "sq_", así
por ejemplo, la nueva variable sq_x = x^2 (al cuadrado).

#
store
@Conjunto de datos
(Guardar los datos)

Guarda un conjunto de datos gretl. Hay dos opciones para el formato de
los datos guardados.

(1) "Formato Standard": los datos se guardan en el formato xml de gretl.
(2) "Comprimido gzip": como arriba, pero utilizando compresión de tipo
     gzip. Esto ahorra espacio de disco, puede ser útil para conjuntos
     de datos grandes.
     
Nótese que si Vd. desea guardar el valor de cualquier escalar generado 
en una sesión de gretl (en lugar de una serie de datos), Vd. debería usar
la instrucción "store" en la ventana de consola de gretl o en un guión
de gretl y especificar la lista de variables a guardar.

#
tsls
@Estimación
(mc2e)
Mínimos cuadrados en dos etapas

Esta orden necesita la selección de dos listas de variables: las 
variables independientes que aparecerán en el modelo dado y un conjunto 
de "instrumentos". Este último comprende las variables exógenas y/o 
predeterminadas que pueden usarse como regresores para obtener valores 
estimados de las variables endógenas del lado derecho e la ecuación del 
modelo. Si alguna de las variables del lado derecho del modelo son
exógenas, éstas deberían aparecer en ambas listas.


#
var
@Estimación
Autorregresión vectorial

Esta instrucción necesita un orden de retardos entero. Hay que 
seleccionar las variables dependiente e independientes para la primera 
ecuación del sistema; estas se permutarán para obtener las
ecuaciones restantes. NO INCLUIR NINGUNA VARIABLE RETARDADA en la lista 
de variables independientes --se añadirán automáticamente.

En general, se realizará una regresión para cada variable de la lista, 
excepto la constante, la tendencia temporal y cualquier variable 
ficticia. La salida de cada ecuación incluye los contrastes F para 
restricciones cero sobre todos los retardos de cada una de las 
variables y un contraste F para el máximo retardo.

#
vartest
@Contrastes
(Contraste de varianzas)

Calcula el estadístico F para la hipótesis nula de que las varianzas
poblacionales para las dos variables elegidas son iguales y muestra su
valor p.

#
wls
@Estimación
(mcp)
Método de mínimos cuadrados ponderados

Sea "varpond" la variable elegida en la caja de "variable de
ponderaciones". Se ejecuta una regresión MCO, donde la variable
dependiente es el producto de varpond por la variable dependiente
elegida y las variables independientes  también se multiplican por
varpond. Si varpond es una variable ficticia, esto es equivalente a 
eliminar todas las observaciones que tengan el valor cero para varpond.

