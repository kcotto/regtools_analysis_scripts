import glob
import subprocess
import os

tags = ['default', 'I']
cwd = os.getcwd()


for tag in tags:
    print(f'Working on {tag} tag')
    lines_per_file = 50000
    smallfile = None
    input_file = f'all_splicing_variants_{tag}.bed'
    with open(input_file, 'r') as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'small_file_{}.txt'.format(lineno + lines_per_file)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    #get chunks
    files = glob.glob('small_file_*')
    print(f'{len(files)} chunk(s) created')
    count = 1
    for file in files:
        print(f'Starting to run file chunk {count}')
        subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/data/compare_junctions_hist_v2.R {tag} {file}', shell=True, check=True)
        count += 1
    output_file = f'compare_junctions/hist/junctions_pvalues_{tag}.tsv'
    if len(files) > 1:
        subprocess.run(f'head -n1 small_file_{lines_per_file}.txt_out.tsv > junctions_pvalues.tsv', shell=True, check=True)
        subprocess.run('for fname in *.txt_out.tsv; do tail -n+2 $fname >> junctions_pvalues.tsv; done', shell=True, check=True)
        if os.path.exists(output_file):
            os.remove(output_file)
        os.rename('junctions_pvalues.tsv', output_file)
    else:
        tmp_file = glob.glob('*.txt_out.tsv')[0]
        if os.path.exists(output_file):
            os.remove(output_file)
        os.rename(tmp_file, output_file)
    subprocess.run('rm -rf small_file*', shell=True, check=True)
