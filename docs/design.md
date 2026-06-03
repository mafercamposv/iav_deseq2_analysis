# Diseño del programa

## Idea general

El programa leerá un archivo de resultados de DESeq2 línea por línea. Para cada gen, extraerá el nombre, log2FoldChange y padj. Después evaluará si el gen cumple los criterios de significancia y lo clasificará como sobreexpresado o subexpresado.

## Estructuras de datos propuestas

Se usará una lista para almacenar genes significativos.

Cada gen leído puede representarse como una tupla:

(gene, log2FoldChange, padj)

Cada gen significativo puede representarse como una tupla:

(gene, log2FoldChange, padj, status)

Donde status puede ser:

- upregulated
- downregulated
```