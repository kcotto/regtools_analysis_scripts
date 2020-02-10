import csv

metadata_file = '/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
matrix_file = '/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct'

sample_to_tissuetype = {}
tissuetype_to_sample = {}

with open(metadata_file, 'r') as meta:
    metadata_reader = csv.reader(meta, delimiter='\t')
    next(metadata_reader)
    for line in metadata_reader:
        sample = line[0]
        tissue_type = line[6]
        sample_to_tissuetype[sample] = tissue_type
        if tissue_type not in tissuetype_to_sample:
            tissuetype_to_sample[tissue_type] = []
        tissuetype_to_sample[tissue_type].append(sample)

set_of_tissues = set(sample_to_tissuetype.values())

output_files_names = []
output_files = {}

for tissue in set_of_tissues:
    tissue_nospace = tissue.replace(" ", "")
    output_files_names.append(f'/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices/{tissue_nospace}.tsv')
    output_files[tissue] = open(f'/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices/{tissue_nospace}.tsv', 'w')
    samples = '\t'.join(tissuetype_to_sample[tissue])
    output_files[tissue].write(f'Name\tDescription\t{samples}')

with open(matrix_file, 'r') as master_matrix:
    line_count = 0
    next(master_matrix)
    next(master_matrix)
    matrix_reader = csv.DictReader(master_matrix, delimiter='\t')
    for line in matrix_reader:
        line_count += 1
        if line_count % 1000 == 0:
            print(line_count)
        name_desc = f'{line["Name"]}\t{line["Description"]}'
        for output_file in output_files.values():
            output_file.write(f'\n{name_desc}')
        for sampletype in sample_to_tissuetype.keys():
            tissue = sample_to_tissuetype[sampletype]
            output_file = output_files[tissue]
            if sampletype in line:
                output_file.write(f'\t{line[sampletype]}')
            else:
                output_file.write('\t0')
