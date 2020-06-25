import subprocess
import os
import glob
import shutil
import sys

def run(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

cohorts = ['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA', 
          'OSCC', 'SCLC', 'HCC']

tags = ['default', 'i50e5', 'E', 'I']

for cohort in cohorts:
    os.makedirs(f'{cohort}/compare_junctions/hist/', exist_ok=True)
    os.chdir(f'{cohort}/compare_junctions/hist/')
    print(f'Copying files for {cohort}')
    run(f'aws s3 cp s3://regtools-results-unstranded/{cohort}/compare_junctions/hist/ . --recursive')
    os.chdir('../..')
    for tag in tags:
        print(f'Running R script for {cohort}_{tag}')
        run(f'Rscript --vanilla /home/ec2-user/regtools/scripts/filter_and_BH.R {tag}')
    files_to_copy = glob.glob('compare_junctions/hist/junction_pvalues_*')
    print(os.getcwd())
    for file in files_to_copy:
        print(file)
        run(f'aws s3 cp ~/{file} cp s3://regtools-results-unstranded/{cohort}/compare_junctions/hist/')
    

