# Proyecto de Análisis de Conexiones Neuronales

Este proyecto se centra en el análisis de conexiones neuronales mediante la aplicación de métodos de correlación cruzada. A continuación, se presentan los puntos más destacados relacionados con la metodología utilizada.

## Requisitos

El código está escrito en **C#**, por lo que es necesario instalar **.NET** y sus dependencias para poder ejecutarlo. Asegúrate de tener la versión adecuada de .NET instalada en tu sistema.

## Método de Correlación

El trabajo se basa en la **correlación cruzada normalizada (NCC)**, que cuantifica las interacciones entre pares de neuronas. Se reconocen dos algoritmos:

1. **Correlación Cruzada Normalizada (NCC)**:
   - Cuantifica la actividad de un neurón objetivo en relación con un neurón de referencia.
   - La fórmula es:

   $C_{xy}(\tau) = \frac{1}{\sqrt{N_x N_y}} \sum_{s=1}^{N_x} x(t) y(t - \tau)$

   - Este método es eficaz para detectar conexiones excitatorias, pero tiene limitaciones en la detección de conexiones inhibitorias.

2. **Filtrado de Correlación Cruzada Normalizada (FNCCH)**:
   - Se introduce un segundo algoritmo para detectar conexiones inhibitorias.
   - La fórmula es:

   $FNCCH = \arg \max_{t} \left| C_{xy}(t) - \frac{1}{W} \sum_{v=-W/2}^{W/2} C_{xy}(v) \right|$

   - Este algoritmo permite distinguir entre conexiones excitatorias e inhibitorias al considerar la media de la NCC y detectar los valores positivos y negativos de disparo.

## Generación de Datos Sustitutos

Además, este código puede generar **datos sustitutos** para evaluar la importancia o significancia de las correlaciones encontradas. El método utilizado para esta generación es el de **shuffling**, que implica reorganizar los datos de forma aleatoria para crear un conjunto de datos que preserve la distribución original pero rompa las correlaciones entre las variables. Esto permite comparar las correlaciones observadas con las que se esperarían por azar.

## Resultados

Los resultados de este análisis permitirán una comprensión más profunda de las interacciones neuronales, especialmente en términos de cómo los diferentes tipos de conexiones influyen en la actividad neuronal. Se presentan ejemplos visuales de cómo se detectan las conexiones excitatorias e inhibitorias en los correlogramas cruzados.

---

Este resumen proporciona una visión clara de los métodos utilizados en el proyecto. Puedes ajustarlo según sea necesario para que se adapte mejor a tu estilo y al contenido del resto del documento.
