# Contexto del proyecto

## Problema

Se desea analizar un archivo de resultados de DESeq2 para identificar genes diferencialmente expresados durante una infección por Influenza A Virus.

## Objetivo del programa

Construir una primera versión de un programa en Python que lea un archivo TSV, filtre genes significativos y clasifique genes sobreexpresados y subexpresados.


## Archivos de entrada

- data/iav_deseq2_results.tsv

## Archivo de salida esperado

- results/iav_significant_genes.tsv

## Ejemplo de salida esperada

```text
gene	log2FoldChange	padj	status
MX1	4.2	0.0001	upregulated
IFIT1	5.1	0.00001	upregulated
```

En esta primera versión trabajaremos únicamente con la información del archivo DESeq2.

Más adelante exploraremos cómo integrar información adicional de genes usando archivos GFF.
