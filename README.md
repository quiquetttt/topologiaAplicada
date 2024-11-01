# topologiaAplicada
Practicas de topologia con jiiiiin

En este entrega debes entregar un archivo de Python que contenga el código relativo a las prácticas propuestas en clase y los ejemplos utilizados para probar las funciones. El archivo debe llamarse "Nombredelalumno.py"

Cada práctica tiene un peso de 2 puntos sobre la nota final. Los contenidos de cada práctica son:

Práctica 1: Complejos simpliciales
-----------------------------------

Clase Complejos simpliciales que permita almacenar complejos simpliciales introduciendo sus símplices maximales y un flotante (filtración) y que contenga métodos que permitan:
Calcular la dimensión del complejo simplicial.
Calcular el conjunto de todas las caras del complejo simplicial.
Calcular el conjunto de todas las caras de una dimensión dada.
Calcular la estrella de un símplice.
Calcular el link de un símplice.
Calcular la característica de Euler.
Calcular el número de componentes de conexas.
Añadir símplices nuevos (con su flotante).
Almacenar la lista de todos los símplices ordenados según el valor del flotante y en caso de que dos tengan el mismo flotante ordenar por dimensión (las caras aparecen primero)
Calcular el complejo simplicial formado por todos los símplices cuyo flotante asociado sea menor o igual que un flotante dado.

Práctica 2: Alfa Complejos
--------------------------

Definir una función que calcule la filtración de alfa complejos asociada a un conjunto de puntos del plano.
Definir una función que represente gráficamente dicho alfa complejo.
Definir una función que calcule la filtración de complejos de Vietoris-Rips asociada a un conjunto de puntos.

Práctica 3: Homología Simplicial
--------------------------------

Definir una función que calcule la forma normal de Smith de una matriz con coeficientes en Z2
Definir dentro de la clase Complejo simplicial métodos que permitan:
 Calcular la matriz borde para cada dimensión.
Calcular los números de Betti.
Calcular los números de Betti de los siguientes complejos simpliciales:
El tetraedro.
El borde del tetraedro.
El toro con las dos triangulaciones vistas en clase.
El plano proyectivo.
La botella de Klein.
El anillo.
El sombrero del asno.
Del complejo simplicial de la transparencia 4 del documento Homología Simplicial II.
Del doble toro.
De algunos alfa complejos.
Crear una función que calcule los números de Betti b_0 y b_1 de un complejo simplicial contenido en el plano utilizando el algoritmo incremental. 
Calcular los números de Betti de algunos alfa complejos del plano utilizando el algoritmo incremental.

Practica 4: Homología persistente
---------------------------------

Definir un método en la clase complejos simpliciales que calcule la matriz borde generalizado de un complejo simplicial filtrado.
Definir una función que calcule el low de una columna de una matriz.
Definir una función que reduzca por columnas una matriz cuadrada según el algoritmo matricula de cálculo de persistencia.
Definir una función que tenga como datos de entrara un conjunto finito de puntos en el plano y que calcule los puntos del diagrama de persistencia de la filtración de alfa complejos asociada a dicho conjunto de puntos.
Definir una función que dibuje el diagrama de persistencia de una filtración de alfa complejos.
Definir una función que dibujo los códigos de barras de una filtración de alfa complejos.
Probar las funciones con distintos conjuntos finitos de puntos en el plano que aproximen conjuntos conocidos (curvas como la circunferencia, la figura ocho, la elipse...)

