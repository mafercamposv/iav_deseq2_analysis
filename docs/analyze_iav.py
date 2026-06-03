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
