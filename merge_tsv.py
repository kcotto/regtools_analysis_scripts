import sys
import csv

def file_to_dict(filename: str) -> tuple:
    result = dict()
    header = []
    keys_in_order = []
    with open(filename, 'r') as infile:
        header = infile.readline().strip().split('\t')
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            key = f'{row[0]}-{row[1]}'
            keys_in_order.append(key)
            result[key] = row
    return (result, header, keys_in_order)

if len(sys.argv) == 4:
    input_file_1 = sys.argv[1]
    input_file_2 = sys.argv[2]
    output_file = sys.argv[3]
else:
    input_file_1 = '/Users/dave/Downloads/Cervix-Ectocervix.tsv'
    input_file_2 = '/Users/dave/Downloads/Cervix-Endocervix.tsv'
    output_file = './cervix.tsv'
file1, header1, keys1_in_order = file_to_dict(input_file_1)
file2, header2, keys2_in_order = file_to_dict(input_file_2)
output_header = header1 + header2[2:]
keys1 = set(file1.keys())
keys2 = set(file2.keys())
keys_in_file1_not_in_file2 = set()
keys_in_file2_not_in_file1 = set()
for key1 in keys1:
    if key1 not in keys2:
        print(f'Warning: Key from file 1: {key1} NOT IN file 2')
        keys_in_file1_not_in_file2.add(key1)
for key2 in keys2:
    if key2 not in keys1:
        print(f'Warning: Key from file 2: {key2} NOT IN file 1')
        keys_in_file2_not_in_file1.add(key2)

data_columns_in_file1 = len(next(iter(file1.items()))[1]) - 2
data_columns_in_file2 = len(next(iter(file2.items()))[1]) - 2
file1_zeros = ['0' for _ in range(data_columns_in_file1)]
file2_zeros = ['0' for _ in range(data_columns_in_file2)]
with open(output_file, 'w') as output:
    output.write('\t'.join(output_header))
    output.write('\n')
    for key1 in keys1_in_order:
        output.write('\t'.join(file1[key1]))
        output.write('\t')
        if key1 in keys_in_file1_not_in_file2:
            output.write('\t'.join(file2_zeros))
        else:
            output.write('\t'.join(file2[key1][2:]))
        output.write('\n')
    for key2 in keys_in_file2_not_in_file1:
        output.write('\t'.join(file2[key2][:2]))
        output.write('\t'.join(file1_zeros))
        output.write('\t')
        output.write('\t'.join(file2[key2][2:]))
        output.write('\n')
