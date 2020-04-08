import csv
import subprocess
import sys
import os
import shutil
def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)
def call(cmd):
    subprocess.call(cmd, shell=True)
results_files = 's3://regtools-results-unstranded'
cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']

def get_vep_variants(cohort, filename):
    output_file = f'{cohort}_vep_splicingvariants.txt'
    f = open(output_file, 'a')
    # write_header = not os.path.exists(output_file)
    # if write_header:
    #     f.write(f'Sample\tCohort\tDA\tD\tA\tNDA\tN\tNovel Junctions\tTotal Junctions\n')
    with open(filename, 'r') as input_samples:
        samples = input_samples.readlines()
        for item in samples:
            fields = item.split('\t')
            if fields[1] == cohort:
                sample = fields[0] 
                vep_vcf = f'{sample}_master.vep.vcf'
                run(f'aws s3 cp {results_files}/VEP_vcfs/{vep_vcf} .')
                with open(vep_vcf, 'r') as vepfile:
                    reader = csv.reader(vepfile, delimiter='\t')
                    for line in reader: 
                        if line[0].startswith('#'):
                            continue
                        else:
                            chrom = line[0]
                            variant_end = line[1]
                            variant_start = int(variant_end) - 1
                            info = line[7]
                            if 'splice' in info:
                                key = f'{cohort}:{chrom}:{variant_start}-{variant_end}'
                                f.write(f'{key}\n')
                os.remove(vep_vcf)
    f.close()
    run(f'aws s3 cp {output_file} {results_files}/VEP_vcfs/splice_variants/')

# for cohort in cohorts:
#     get_vep_variants(cohort, 'primarysolidtumors_metadata.tsv')

if __name__ == '__main__':
    if len(sys.argv) == 1:
        for cohort in cohorts:
            call(f'./.local/bin/pipenv run python get_vep_splicingvariants.py {cohort} &')
    else:
        get_vep_variants(sys.argv[1], 'primarysolidtumors_metadata.tsv')