# Metodología de Análisis de Actividad Neuronal

## Introducción
En esta investigación, analizamos la actividad neuronal a partir de datos electrofisiológicos. El conjunto de datos consiste en tres horas de registros electrofisiológicos, en los cuales calculamos correlaciones para entender mejor las interacciones neuronales.

## Segmentación de Datos
Para minimizar los efectos negativos causados por la no estacionariedad de las series temporales, dividimos los datos en segmentos de 250 segundos. Dentro de cada segmento, se calculan los coeficientes de variación cada 10 segundos, promediándose luego para obtener un único valor representativo.

## Cálculo de Correlaciones
### Desafíos
El cálculo de todas las correlaciones cruzadas es una de las tareas más intensivas en este trabajo, debido a la gran cantidad de correlaciones a calcular y la necesidad de generar simultáneamente datos surrogados para aplicar un umbral que elimine correlaciones espúreas.

### Herramienta Utilizada
Para abordar estos desafíos, utilizamos el software **Spicodyn**, diseñado para ofrecer alta eficiencia computacional en el cálculo de correlaciones y la generación de datos surrogados. Sin embargo, Spicodyn presenta limitaciones al estar diseñado principalmente para experimentos *in vitro*, lo que requiere una adaptación para datos obtenidos en experimentos *in vivo*.

## Proceso de Construcción de Matrices de Conectividad
El proceso para construir matrices de conectividad funcional implica varias etapas:

1. **Cálculo del Coeficiente de Variación (CV)**: Se calculan los CV en intervalos de 10 segundos, organizándose y agrupándose en conjuntos de CVs adyacentes.
2. **Extracción de Datos**: Para cada serie de 250 segundos, se extraen los datos por neurón, permitiendo el cálculo de todas las correlaciones entre pares de neuronas.
3. **Definición de Parámetros**: Se definen variables locales necesarias antes de iniciar el cálculo, como la ventana de correlación, el número de datos surrogados, y la frecuencia de muestreo.

### Importancia de la Clasificación por Disparos
La clasificación por disparos es crucial para garantizar que los datos obtenidos se puedan identificar con neuronas individuales, lo cual es fundamental para el análisis.

## Resultados
- **Matrices de Conectividad Funcional**: Se construyeron matrices de conectividad funcional a partir de las correlaciones cruzadas calculadas. Estas matrices permiten visualizar las interacciones neuronales significativas.
- **Visualizaciones**: Se presentan gráficos que muestran los coeficientes de variación y las matrices de conectividad, facilitando la interpretación de los resultados.

## Conclusión
Este proceso metodológico permite un análisis detallado de la actividad neuronal, proporcionando insights sobre la conectividad funcional entre neuronas. En futuras investigaciones, se continuarán explorando métricas de red para profundizar en la comprensión de las dinámicas neuronales.

## Referencias
- Spicodyn: Herramienta utilizada para el cálculo de correlaciones y generación de datos surrogados.
- Neuroimaging Tools and Resources Collaboratory (NITRC): Entidad que respalda el desarrollo de herramientas en neurociencia.

---

**Fuente**: Elaborada por el autor (2023)



