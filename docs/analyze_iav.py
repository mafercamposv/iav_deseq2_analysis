def load_deseq2_results(filename):
    valid_genes = []
    with open(filename, "r") as file:
        header = file.readline()  # Leer el encabezado
        for line in file:
            if line.strip() == "":
                continue  # Ignorar líneas vacías
            columns = line.split("\t")
            if len(columns) < 3:
                continue  # Ignorar líneas incompletas
            gene_name = columns[0]
            try:
                log2_fold_change = float(columns[1])
                p_value = float(columns[2])
                valid_genes.append((gene_name, log2_fold_change, p_value))
            except ValueError:
                continue  # Ignorar líneas con valores no numéricos
    return valid_genes


# Responsabilidad: Cargar los resultados de DESeq2 desde un archivo y devolver una lista de genes válidos con sus estadísticas.


def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True
    else:
        return False


# Responsabilidad: Determinar si un gen es significativamente diferencialmente expresado según los umbrales de log2 fold change y p-valor ajustado.


def classify_gene(log2_fold_change):
    if log2_fold_change > 0:
        return "upregulated"
    elif log2_fold_change < 0:
        return "downregulated"
    else:
        return "not significant"


# Responsabilidad: Clasificar un gen como "upregulated", "downregulated" o "not significant" según su log2 fold change.


def filter_genes(results, lfc_threshold, padj_threshold):
    filtered_genes = []
    for gene_name, log2_fold_change, p_value in results:
        if is_significant(log2_fold_change, p_value, lfc_threshold, padj_threshold):
            classification = classify_gene(log2_fold_change)
            filtered_genes.append(
                (gene_name, log2_fold_change, p_value, classification)
            )
    return filtered_genes


# Responsabilidad: Filtrar los genes significativos y clasificarlos según su log2 fold change, devolviendo una lista de genes con su clasificación.


def write_results(output_file, filtered_genes):
    with open(output_file, "w") as file:
        file.write(
            "Gene\tLog2FoldChange\tPValue\tClassification\n"
        )  # Escribir encabezado
        for gene_name, log2_fold_change, p_value, classification in filtered_genes:
            file.write(
                f"{gene_name}\t{log2_fold_change}\t{p_value}\t{classification}\n"
            )  # Escribir genes filtrados


# Responsabilidad: Escribir los genes filtrados y clasificados en un archivo de salida con un formato específico.
