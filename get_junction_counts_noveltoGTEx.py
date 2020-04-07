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
# cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
#           'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
#           'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
#           'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']
cohorts = ['CHOL']

def create_gtex_set(filenname):
    gtex = set()
    with open(filenname, 'r') as gtex_file:
        next(gtex_file)
        junctions = gtex_file.readlines()
        for junction in junctions:
            junction = junction.split('\t')[0]
            gtex.add(junction)
    return gtex

def get_junction_counts(cohort, filename, gtex):
    output_file = f'{cohort}_novelcounts.txt'
    f = open(output_file, 'a')
    # write_header = not os.path.exists(output_file)
    # if write_header:
    #     f.write(f'Sample\tCohort\tDA\tD\tA\tNDA\tN\tNovel Junctions\tTotal Junctions\n')
    with open(filename, 'r') as input_samples:
        samples = input_samples.readlines()
        for item in samples:
            counts = {'DA':0, 'D':0, 'A':0, 'NDA':0, 'N':0}
            fields = item.split('\t')
            if fields[1] == cohort:
                sample = fields[0] 
                run(f'aws s3 cp {results_files}/{cohort}/{sample}.tar.gz .')
                annotated_bed = f'{sample}/{sample}_annotated.bed'
                run(f'tar xzf {sample}.tar.gz {annotated_bed}')
                with open(annotated_bed, 'r') as bedfile:
                    reader = csv.DictReader(bedfile, delimiter='\t')
                    for line in reader: 
                        chrom = line['chrom']
                        junction_start = str(int(line['start']) + 1)
                        junction_end = str(int(line['end']) - 1)
                        regtools_key = f'{chrom}_{junction_start}_{junction_end}'
                        if regtools_key not in gtex:  
                            counts[line['anchor']] += 1
                total_junctions = sum(counts.values())
                novel_junctions = total_junctions - counts['DA']
                f.write(f"{sample}\t{cohort}\t{counts['DA']}\t{counts['D']}\t{counts['A']}\t{counts['NDA']}\t{counts['N']}\t{novel_junctions}\t{total_junctions}\n")
                shutil.rmtree(sample, ignore_errors=True)
                os.remove(f'{sample}.tar.gz')
    f.close()
    run(f'aws s3 cp {output_file} {results_files}/junction_counts/')

# for cohort in cohorts:
#     gtex = create_gtex_set('master_GTEx_b38_mean_sd.tsv')
#     get_junction_counts(cohort, 'primarysolidtumors_metadata.tsv', gtex)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        for cohort in cohorts:
            call(f'./.local/bin/pipenv run python get_junction_counts.py {cohort} &')
    else:
        gtex = create_gtex_set('master_GTEx_b38_mean_sd.tsv')
        get_junction_counts(sys.argv[1], 'primarysolidtumors_metadata.tsv', gtex)