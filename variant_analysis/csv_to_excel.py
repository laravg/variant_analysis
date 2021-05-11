import pandas as pd
import os

path = 'results/resultsPrimersArchaea'

for root, subdirs, files in os.walk(path):
    for filename in files:
        if '_all.csv' in filename:
            file_path = os.path.join(root, filename)
            read_file = pd.read_csv(file_path, sep=';')
            if 'species' in filename:
                cols = ['Genomes_mean_size', 'Genes_mean_size', 'Num_genes',
                        'Num_variant_genes', 'Num_genes_+strand', 'Num_genes_-strand']
                read_file[cols] = read_file[cols].replace({'_': 0})
                read_file[cols] = read_file[cols].apply(pd.to_numeric)
            excel_file_path = file_path.replace('.csv', '.xlsx')
            read_file.to_excel(excel_file_path, index=None, header=True)
