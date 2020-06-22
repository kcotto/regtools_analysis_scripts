import csv

with open('/Users/kcotto/Downloads/junction_pvalues_default_old.tsv', 'r') as old, open('/Users/kcotto/Downloads/junction_pvalues_default_new.tsv', 'r') as new, open('/Users/kcotto/Downloads/junction_pvalues_default_merged.tsv', 'w') as outputfile:
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
