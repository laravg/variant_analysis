
import pandas as pd
import os
from pathlib import Path

groups = {714: ["CP012959", "CP016553", "CP012958"],
          1597: ["CP002618,CP002619", "FM177140", "CP002616,CP002617", "CP005486,CP005487", "CP001084,CP000935"],
          76857: ["CP013121", "LN831027"],
          76859: ["CP012713", "CP012715"],
          712710: ["CP017038", "CP028365"]
          }


path = 'results/resultsPrimersBacteriaRERUNstrains'
pathpath = Path(path)
for dir in os.listdir(path):
    dirpath = pathpath / dir
    print(dirpath)
    if dirpath.is_dir() and 'genomes' not in str(dirpath):
        for group in groups:
            new_df = pd.DataFrame()
            for id in groups[group]:
                name = f'{id}.txt_{dir}'
                species_path = dirpath / name / '7_species_info' / 'species_all.csv'
                species_df = pd.read_csv(species_path, sep=';')
                new_df = pd.concat([new_df, species_df])
            new_df.to_csv(dirpath / f'{group}.csv', sep=';', index=False)
        # for file in os.listdir(dirpath):
        #     if file
        #     print(file)
    # for filename in files:
    #     if '_all.csv' in filename:
    #         file_path = os.path.join(root, filename)
    #         read_file = pd.read_csv(file_path, sep=';')
    #         if 'species' in filename:
    #             cols = ['Genomes_mean_size', 'Genes_mean_size', 'Num_genes',
    #                     'Num_unique_genes', 'Num_genes_+strand', 'Num_genes_-strand']
    #             read_file[cols] = read_file[cols].replace({'_': 0})
    #             read_file[cols] = read_file[cols].apply(pd.to_numeric)
    #         excel_file_path = file_path.replace('.csv', '.xlsx')
    #         read_file.to_excel(excel_file_path, index=None, header=True)
