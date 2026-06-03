# Análisis de resultados DESeq2

Este proyecto contiene un script en Python para filtrar y clasificar genes a partir de resultados de DESeq2.

## Qué hace el programa

El script `analyze_iav.py`:

- Lee un archivo TSV con resultados de DESeq2.
- Extrae el nombre del gen, el valor de `log2FoldChange` y el valor ajustado de p (`padj`).
- Filtra genes que cumplen criterios de significancia.
- Clasifica los genes significativos como `upregulated` o `downregulated`.
- Escribe el resultado en un archivo de salida TSV.
- Muestra un resumen en pantalla con el total de genes significativos y el número de genes sobreexpresados y subexpresados.

## Formato de entrada esperado

El archivo de entrada debe ser un TSV con un encabezado que incluya, como mínimo, las columnas:

- `gene_id` o similar en la primera columna
- `log2FoldChange`
- `padj`

Ejemplo de encabezado válido:

```text
gene_id	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
```

El script detecta automáticamente las columnas de `log2FoldChange` y `padj` a partir del encabezado. Si el encabezado no contiene esas etiquetas, asumirá que `log2FoldChange` está en la segunda columna y `padj` en la tercera.

El programa ignora:

- líneas vacías
- líneas con menos columnas de las necesarias
- valores no numéricos en `log2FoldChange` o `padj`

## Formato de salida

El archivo de salida es un TSV con las columnas:

- `Gene`
- `Log2FoldChange`
- `padj`
- `Classification`

Ejemplo de salida:

```text
Gene	Log2FoldChange	padj	Classification
MX1	8.69	0.006693	upregulated
IFIT1	8.69	0.006693	upregulated
```

## Criterios de filtrado

Por defecto, el programa considera un gen significativo cuando:

- `padj < 0.05`
- `abs(log2FoldChange) >= 1`

Puedes cambiar estos umbrales desde la línea de comandos.

## Uso

Ejemplo de ejecución:

```bash
python analyze_iav.py data/iav_deseq2_results.tsv results/iav_significant_genes.tsv
```

Con umbrales personalizados:

```bash
python analyze_iav.py data/iav_deseq2_results.tsv results/iav_significant_genes.tsv \
  --lfc_threshold 2.0 --padj_threshold 0.01
```

## Notas

- El script actual solo procesa los resultados de DESeq2 y no agrega información adicional de otros archivos.
- Si el archivo de entrada no existe o no se puede leer, el script muestra un mensaje de error genérico.
