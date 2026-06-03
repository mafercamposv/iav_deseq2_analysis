## Influenza A Virus (IAV) · Python para Bioinformática



# Objetivo general

Construir una primera versión de un programa en Python que permita analizar resultados de expresión diferencial obtenidos con DESeq2 para identificar genes significativamente sobreexpresados y subexpresados durante una infección por Influenza A Virus (IAV).

El objetivo no es únicamente programar, sino aprender a:

- pensar la lógica antes de escribir código
- diseñar una solución paso a paso
- trabajar con archivos tabulares reales de bioinformática
- documentar el problema, el diseño y los casos de prueba
- usar Python como herramienta para resolver problemas científicos
- utilizar IA como apoyo razonado y no como sustituto del análisis

En esta guía construiremos primero una versión mínima funcional del programa.

La versión mínima hará lo siguiente:

1. Leer un archivo TSV con resultados de DESeq2.
2. Identificar genes diferencialmente expresados usando `padj` y `log2FoldChange`.
3. Clasificarlos como sobreexpresados o subexpresados.
4. Guardar una tabla de resultados.
5. Mostrar un resumen simple.

La integración con el archivo GFF se trabajará después como una extensión biológica opcional.

Durante la guía construiremos el programa de forma incremental. No intentaremos resolver todo al mismo tiempo: primero desarrollaremos una versión mínima funcional y después discutiremos posibles mejoras.



# 1. Contexto biológico

## ¿Qué es DESeq2?

DESeq2 es una herramienta ampliamente utilizada en bioinformática para analizar datos de RNA-seq.

Su objetivo principal es identificar genes que cambian significativamente su expresión entre dos condiciones biológicas.

Por ejemplo:

- células sanas vs células infectadas
- tratamiento vs control
- cáncer vs tejido normal

DESeq2 genera una tabla donde cada fila representa un gen y cada columna contiene información estadística sobre el cambio de expresión.

Ejemplo simplificado:

| gene | log2FoldChange | padj |
|---|---:|---:|
| MX1 | 4.2 | 0.0001 |
| IFIT1 | 5.1 | 0.00001 |
| GAPDH | 0.1 | 0.89 |



## ¿Qué son las células control?

Las células control son muestras utilizadas como referencia.

Representan la condición basal, normal o sin tratamiento.

En este ejercicio:

- las células control corresponden a células no infectadas
- las células experimentales corresponden a células infectadas con Influenza A Virus

El programa analizará qué genes cambian su expresión al comparar células infectadas contra células control.



## ¿Qué es `log2FoldChange`?

La columna `log2FoldChange` indica cuánto cambia la expresión de un gen entre dos condiciones biológicas.

En este ejercicio:

- condición control → células no infectadas
- condición experimental → células infectadas con Influenza A Virus

Interpretación:

- valores positivos → el gen aumenta su expresión en células infectadas
- valores negativos → el gen disminuye su expresión en células infectadas
- valores cercanos a cero → casi no cambia

Ejemplo:

| log2FoldChange | Interpretación |
|---:|---|
| 5.0 | fuerte aumento |
| 2.0 | aumento moderado |
| 0.1 | casi sin cambio |
| -2.5 | disminución importante |



## ¿Qué es `padj`?

La columna `padj` representa el valor p ajustado (*adjusted p-value*).

Se utiliza para estimar qué tan confiable es el cambio observado.

En RNA-seq se analizan miles de genes al mismo tiempo, por lo que DESeq2 corrige estadísticamente los valores p para reducir falsos positivos.

Interpretación general:

- valores pequeños → mayor confianza en que el cambio es real
- valores grandes → el cambio podría deberse al azar

Ejemplo:

| padj | Interpretación |
|---:|---|
| 0.0001 | muy significativo |
| 0.01 | significativo |
| 0.2 | poco confiable |
| 0.8 | probablemente sin cambio real |



## ¿Qué significa que la expresión de un gen suba o baje?

La combinación de `log2FoldChange` y `padj` permite decidir:

- si el gen cambia significativamente
- y hacia qué dirección cambia

Por ejemplo:

| gene | log2FoldChange | padj | Interpretación |
|---|---:|---:|---|
| MX1 | 4.2 | 0.0001 | aumento fuerte y significativo |
| IFIT1 | 5.1 | 0.00001 | aumento fuerte y significativo |
| GAPDH | 0.1 | 0.89 | sin cambio relevante |

Por eso, si un gen como `MX1` tiene un `log2FoldChange` positivo y un `padj` pequeño, podemos decir que aumenta significativamente su expresión en células infectadas respecto al control.


## Umbrales clásicos de expresión diferencial

En muchos análisis de expresión diferencial se usan criterios como:

- `padj < 0.05`
- `abs(log2FoldChange) >= 1`

Interpretación:

- `padj < 0.05` sugiere que el cambio es estadísticamente confiable
- `abs(log2FoldChange) >= 1` sugiere que el cambio tiene una magnitud relevante

En este ejercicio usaremos esos umbrales como criterio básico.



## ¿Qué significa `abs(log2FoldChange)`?

La función `abs()` obtiene el valor absoluto de un número.

Esto significa que elimina el signo y permite medir la magnitud del cambio sin importar si el gen sube o baja.

Ejemplos:

| Valor | `abs(valor)` |
|---:|---:|
| 2.5 | 2.5 |
| -2.5 | 2.5 |
| 0.3 | 0.3 |

En este ejercicio usamos:

```text
abs(log2FoldChange) >= 1
```

Esto significa que nos interesan genes cuyo cambio sea suficientemente grande, ya sea positivo o negativo.



## ¿Cómo se traduce esto a una decisión computacional?

Para que el programa pueda tomar una decisión, convertiremos el criterio biológico en una regla lógica.

Un gen será considerado diferencialmente expresado si cumple dos condiciones al mismo tiempo:

1. `padj < 0.05`
2. `abs(log2FoldChange) >= 1`

Después, usaremos el signo de `log2FoldChange` para clasificarlo:

- si `log2FoldChange > 0`, el gen se clasificará como sobreexpresado
- si `log2FoldChange < 0`, el gen se clasificará como subexpresado

Esta será la regla central del programa.
