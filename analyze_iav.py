def load_deseq2_results(filename):
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
    for gene_name, log2_fold_change, padj in results:
        if is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
            classification = classify_gene(log2_fold_change)
            filtered_genes.append((gene_name, log2_fold_change, padj, classification))
    return filtered_genes


# Responsabilidad: Filtrar los genes significativos y clasificarlos según su log2 fold change, devolviendo una lista de genes con su clasificación.


def write_results(output_file, filtered_genes):
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
