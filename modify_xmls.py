import subprocess
import glob
import csv
import sys
import os
import shutil
from pathlib import Path
import argparse

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

def edit_spliceai(filename, cohort):
    results_files = 's3://regtools-igv-sessions'
    with open(filename, 'r') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        for line in reader:
            spliceai = line['SpliceAI_raw']
            if spliceai != 'NA':
                variant_junction = line['variant_junction_info']
                variant_junction = variant_junction.replace(':', '_')
                samples_field = line['junction_samples']
                samples_field = samples_field.replace(',','_')
                xml_file = '_'.join([samples_field, variant_junction]) + '.xml'
                tag = filename.split('_')[7]
                run(
                f'aws s3 cp {results_files}/{cohort}/{tag}/{xml_file} .')
                with open(xml_file, 'r') as starting_file, open(xml_file + '_new', 'w') as temp_file:
                    lines = starting_file.readlines()
                    for line in lines:
                        temp_file.write(line)
                        if '<Resources>' in line:
                            temp_file.write(f'<Resource path="https://regtools-igv-files.s3.amazonaws.com/spliceAI/{cohort}_{variant_junction}.bed"/>\n')
                os.remove(xml_file)
                os.rename(f'{xml_file}_new', xml_file)
                run(
                f'aws s3 cp {xml_file} {results_files}/{cohort}/{tag}/{xml_file}')


cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']
results_files = 's3://regtools-results-unstranded'

for cohort in cohorts:
    if os.path.exists(f'{cohort}_igv_session'):
        shutil.rmtree(f'{cohort}_igv_session')
    os.mkdir(f'{cohort}_igv_session')
    os.chdir(f'{cohort}_igv_session')

    run(
        f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/{cohort}_junction_pvalues_significant_0.05_filtered_BH_default_gtex_spliceai_w_IGVsessions.tsv .')
    run(
        f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/{cohort}_junction_pvalues_significant_0.05_filtered_BH_i50e5_gtex_spliceai_w_IGVsessions.tsv .')
    run(
        f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/{cohort}_junction_pvalues_significant_0.05_filtered_BH_E_gtex_spliceai_w_IGVsessions.tsv .')
    run(
        f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/{cohort}_junction_pvalues_significant_0.05_filtered_BH_I_gtex_spliceai_w_IGVsessions.tsv .')

    files = glob.glob('*junction_pvalues_significant_0.05_filtered_BH*.tsv')
    for file in files:
        edit_spliceai(file, cohort)
    os.chdir('..')
    if os.path.exists(f'{cohort}_igv_session'):
        shutil.rmtree(f'{cohort}_igv_session')


