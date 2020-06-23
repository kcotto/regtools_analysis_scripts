import csv
import subprocess
import sys
import os

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

# with open('/home/ec2-user/old_files.txt', 'r') as old, open('/home/ec2-user/fixed_files.txt', 'r') as new:
#     old_files = old.readlines()
#     new_files = new.readlines()
#     old_file_keys = set()
#     new_file_keys = set()
#     for file in old_files:
#         fields = file.strip().split('/')
#         key = f'{fields[0]}|{fields[-1]}'
#         old_file_keys.add(key)
#     for file in new_files:
#         fields = file.strip().split('/')
#         key = f'{fields[0]}|{fields[-1]}'
#         new_file_keys.add(key)
    # keys_to_copy = old_file_keys - new_file_keys
    # for key in keys_to_copy:
    #     fields = key.split('|')
    #     run(f'aws s3 cp s3://regtools-results-unstranded/{fields[0]}/compare_junctions2/hist/{fields[1]} s3://regtools-results-unstranded/{fields[0]}/compare_junctions/hist/{fields[1]}')


with open('/home/ec2-user/fixed_files.txt', 'r') as new:
    new_files = new.readlines()
    for new_file in new_files:
        fields = new_file.strip().split('/')
        run(f'aws s3 cp s3://regtools-results-unstranded/{fields[0]}/compare_junctions3/hist/{fields[-1]} new_file.tsv')
        run(f'aws s3 cp s3://regtools-results-unstranded/{fields[0]}/compare_junctions2/hist/{fields[-1]} old_file.tsv')
        with open('old_file.tsv', 'r') as old, open('new_file.tsv', 'r') as new, open('merged.tsv', 'w') as outputfile:
            old_reader = csv.reader(old, delimiter='\t')
            new_reader = csv.reader(new, delimiter='\t')
            keys = set()
            header = next(new_reader)
            header_line = '\t'.join(header)
            outputfile.write(header_line + '\n')
            next(old_reader)
            for line in new_reader:
                keys.add(line[9])
                joined_line = '\t'.join(line)
                outputfile.write(joined_line + '\n')
            for line in old_reader:
                if line[9] not in keys:
                    joined_line = '\t'.join(line)
                    outputfile.write(joined_line + '\n')
            run(f'aws s3 cp merged.tsv s3://regtools-results-unstranded/{fields[0]}/compare_junctions/hist/{fields[-1]}')
        os.remove('old_file.tsv')
        os.remove('new_file.tsv')
        os.remove('merged.tsv')
