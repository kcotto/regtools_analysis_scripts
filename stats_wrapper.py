import glob
import subprocess
import os

tags = ['default', 'i50e5', 'E', 'I']
cwd = os.getcwd()


for tag in tags:
    target_lines_per_file = 75000
    lines_per_file = 0
    input_file = f'all_splicing_variants_{tag}.bed'
    lines = open(input_file).readlines()
    count = len(lines)
    if count <= lines_per_file:
        subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/data/compare_junctions_hist_v2.R {tag} {input_file}')
    else:
        header = lines[0]
        lines.pop(0)
        lines.sort()
        filenum = 1
        small_filename = f'small_file_{filenum}.txt'
        smallfile = open(small_filename, "w")
        smallfile.write(header)
        lines_per_file += target_lines_per_file
        for lineno, line in enumerate(lines):
            smallfile.write(line)
            if lineno >= lines_per_file:
                fields1 = line.split('\t')
                variant1 = f'{fields1[0]}_{fields1[1]}_{fields1[2]}'
                fields2 = lines[lineno+1].split('\t')
                variant2 = f'{fields2[0]}_{fields2[1]}_{fields2[2]}'
                if variant1 != variant2:
                    smallfile.close()
                    filenum += 1
                    small_filename = f'small_file_{filenum}.txt'
                    smallfile = open(small_filename, "w")
                    smallfile.write(header)
                    lines_per_file += target_lines_per_file
    # get chunks
    files = glob.glob('small_file_*')
    for file in files:
        subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/data/compare_junctions_hist_v2.R {tag} {file}')
    subprocess.run(f"awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print} ' small_file*.txt > junction_pvalues_{tag}.tsv")
    os.remove('small_file*')
