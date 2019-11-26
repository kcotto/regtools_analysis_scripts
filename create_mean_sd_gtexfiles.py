import glob
import pandas as pd

files = glob.glob('/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices/*')

for file in files:
    outputfile_beginning = '/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices_avg'
    file_name = file.split('/')[-1]
    new_file = f'{outputfile_beginning}/{file_name}'
    data = pd.read_csv(file, sep='\t')
    data['GTEx_mean'] = data.mean(axis=1)
    data['GTEx_sd'] = data.std(axis=1)
    columns_to_save = ['Name', 'Description', 'GTEx_mean', 'GTEx_sd']
    data.to_csv(new_file, sep='\t', columns=columns_to_save, index=False)
    print(f'Finished writing {new_file}')
