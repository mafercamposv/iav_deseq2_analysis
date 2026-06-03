def load_deseq2_results(filename):
    """Read a DESeq2 result file and return valid gene records.

    Args:
        filename (str): Path to the DESeq2 TSV file.

    Returns:
        list[tuple[str, float, float]]: A list of tuples with gene name,
            log2 fold change, and adjusted p-value (padj).

    Notes:
        The function detects the columns `log2FoldChange` and `padj` from the
        header. If the header does not include those columns, it uses the
        second and third columns by default.
    """
    valid_genes = []
    with open(filename, "r") as file:
        header = file.readline().strip()  # Leer el encabezado
        header_columns = header.split("\t")
        header_lower = [column.lower() for column in header_columns]

        try:
            lfc_index = header_lower.index("log2foldchange")
        except ValueError:
            lfc_index = 1

        if "padj" in header_lower:
            padj_index = header_lower.index("padj")
        elif "pvalue" in header_lower:
            padj_index = header_lower.index("pvalue")
        else:
            padj_index = 2

        for line in file:
            if line.strip() == "":
                continue  # Ignorar líneas vacías
            columns = line.strip().split("\t")
            if len(columns) <= max(0, lfc_index, padj_index):
                continue  # Ignorar líneas incompletas
            gene_name = columns[0]
            try:
                log2_fold_change = float(columns[lfc_index])
                padj = float(columns[padj_index])
                valid_genes.append((gene_name, log2_fold_change, padj))
            except ValueError:
                continue  # Ignorar líneas con valores no numéricos
    return valid_genes


# Responsabilidad: Cargar los resultados de DESeq2 desde un archivo y devolver una lista de genes válidos con sus estadísticas. Detecta la columna de log2FoldChange y padj a partir del encabezado.


def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    """Determine whether a gene passes significance thresholds.

    Args:
        log2_fold_change (float): Log2 fold change for the gene.
        padj (float): Adjusted p-value for the gene.
        lfc_threshold (float): Minimum absolute log2 fold change required.
        padj_threshold (float): Maximum adjusted p-value allowed.

    Returns:
        bool: True if the gene is significant, otherwise False.
    """
    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True
    else:
        return False


# Responsabilidad: Determinar si un gen es significativamente diferencialmente expresado según los umbrales de log2 fold change y p-valor ajustado.


def classify_gene(log2_fold_change):
    """Classify a gene as upregulated or downregulated.

    Args:
        log2_fold_change (float): Log2 fold change for the gene.

    Returns:
        str: "upregulated" if positive, "downregulated" if negative,
            or "not significant" if zero.
    """
    if log2_fold_change > 0:
        return "upregulated"
    elif log2_fold_change < 0:
        return "downregulated"
    else:
        return "not significant"


# Responsabilidad: Clasificar un gen como "upregulated", "downregulated" o "not significant" según su log2 fold change.


def filter_genes(results, lfc_threshold, padj_threshold):
    """Filter DESeq2 gene records by significance and classify them.

    Args:
        results (list[tuple[str, float, float]]): Gene records with name,
            log2 fold change, and adjusted p-value.
        lfc_threshold (float): Minimum absolute log2 fold change required.
        padj_threshold (float): Maximum adjusted p-value allowed.

    Returns:
        list[tuple[str, float, float, str]]: Filtered gene records including
            classification labels.
    """
    filtered_genes = []
    for gene_name, log2_fold_change, padj in results:
        if is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
            classification = classify_gene(log2_fold_change)
            filtered_genes.append((gene_name, log2_fold_change, padj, classification))
    return filtered_genes


# Responsabilidad: Filtrar los genes significativos y clasificarlos según su log2 fold change, devolviendo una lista de genes con su clasificación.


def write_results(output_file, filtered_genes):
    """Write filtered gene records to an output TSV file.

    Args:
        output_file (str): Path to the output TSV file.
        filtered_genes (list[tuple[str, float, float, str]]): Filtered gene
            records with classification.
    """
    with open(output_file, "w") as file:
        file.write(
            "Gene\tLog2FoldChange\tpadj\tClassification\n"
        )  # Escribir encabezado
        for gene_name, log2_fold_change, padj, classification in filtered_genes:
            file.write(
                f"{gene_name}\t{log2_fold_change}\t{padj}\t{classification}\n"
            )  # Escribir genes filtrados


# Responsabilidad: Escribir los genes filtrados y clasificados en un archivo de salida con un formato específico.


def print_summary(filtered_genes):
    """Print a summary of filtered gene classifications.

    Args:
        filtered_genes (list[tuple[str, float, float, str]]): Filtered gene
            records with classification.
    """
    total_genes = len(filtered_genes)
    upregulated_genes = sum(1 for gene in filtered_genes if gene[3] == "upregulated")
    downregulated_genes = sum(
        1 for gene in filtered_genes if gene[3] == "downregulated"
    )
    print(f"Total significant genes: {total_genes}")
    print(f"Upregulated genes: {upregulated_genes}")
    print(f"Downregulated genes: {downregulated_genes}")


# Responsabilidad: Imprimir un resumen con el conteo de genes significativos y su clasificación.


def main():
    """Parse arguments and run the DESeq2 analysis pipeline.

    This is the program entry point. It loads input data, filters genes,
    writes results to a file, and prints a summary.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Analyze DESeq2 results")
    parser.add_argument("input_file", help="Path to the DESeq2 results file")
    parser.add_argument(
        "output_file", help="Path to the output file for filtered results"
    )
    parser.add_argument(
        "--lfc_threshold", type=float, default=1.0, help="Log2 fold change threshold"
    )
    parser.add_argument(
        "--padj_threshold", type=float, default=0.05, help="Adjusted p-value threshold"
    )
    args = parser.parse_args()

    try:
        results = load_deseq2_results(args.input_file)
        filtered_genes = filter_genes(results, args.lfc_threshold, args.padj_threshold)
        write_results(args.output_file, filtered_genes)
        print_summary(filtered_genes)
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
