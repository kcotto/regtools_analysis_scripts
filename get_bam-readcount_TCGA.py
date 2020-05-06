import subprocess
import glob
import csv
import sys
import os
import requests
import shutil
import json
from pathlib import Path


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)


def call(cmd):
    subprocess.call(cmd, shell=True)


def get_bam(sample: str, token_file: str, chrom: str, bam_window_start: int, bam_window_end: int) -> None:
    print('Obtaining bam manifest file')
    bam_url = 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22BAM%22%5D%7D%7D%2C%7B%22op%22%3A%22AND%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%2C%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.submitter_id%22%2C%22value%22%3A%5B%22{0}%22%5D%7D%7D%5D%7D%5D%7D%5D%7D%5D%7D&query=files.data_format%20in%20%5B%22BAM%22%5D%20and%20files.experimental_strategy%20in%20%5B%22RNA-Seq%22%5D%20AND%20cases.project.program.name%20IN%20%5BTCGA%5D%20and%20cases.samples.submitter_id%20IN%20%5B%22P{0}&return_type=manifest'.format(sample)
    response = requests.get(bam_url)
    if response.status_code != 200:
        print(sample, f'Failed to download bam manifest file: HTTP Status Code: {response.status_code}')
    manifest_contents = response.content.decode('utf-8').strip()
    if not manifest_contents.startswith('id\t'):
        print(sample, f'bam manifest did not have expected format')


    bam_manifest_path = Path(f'{sample}_bam_manifest.txt')
    bam_manifest_path.write_text(manifest_contents)
    run(f'aws s3 cp {token_file} .')
    with open(bam_manifest_path, 'r') as manifest_file:
        reader = csv.reader(manifest_file, delimiter='\t')
        next(manifest_file)
        for line in reader:
            bam_id = line[0]
            bam_file = f'{sample}_{chrom}_{bam_window_start}_{bam_window_end}.bam'
            token_file = "gdc-user-token.txt"

            file_id = f'{bam_id}'

            data_endpt = "https://api.gdc.cancer.gov/slicing/view/{}".format(file_id)

            with open(token_file, "r") as token:
                token_string = str(token.read().strip())

            params = {"region": [f'{chrom}:{bam_window_start}-{bam_window_end}']}

            response = requests.post(data_endpt,
                                     data=json.dumps(params),
                                     headers={
                                         "Content-Type": "application/json",
                                         "X-Auth-Token": token_string
                                     })
            print(f'HTTP POST returned {response.status_code}')
            if response.status_code != 200:
                # If failure, queue the item back to the error queue
                error_json = response.json()
                error_message = error_json['error']
                # write error_message to status queue
                sqs.report_status(job_info, f'Error message: {error_message}')
                sqs.queue_error_item(job_info)
            else:            
                # else do this
                with open(bam_file, "wb") as output_file:
                    output_file.write(response.content)
                run(f'samtools index {bam_file}')


def run_bam_readcount(cohort, filename):
    token_file = 's3://regtools-cwl-sharedfiles/gdc-user-token.txt'
    output_file = f'{cohort}_RNF145_variant_counts.txt'
    f = open(output_file, 'a')
    chrom =  'chr5'
    position = 159169058
    bam_start = int(position) - 1000
    bam_stop = int(position) + 1000
    with open(filename, 'r') as input_samples:
        samples = input_samples.readlines()
        for item in samples:
            fields = item.split('\t')
            file_cohort = fields[1].strip()
            if file_cohort == cohort:
                sample = fields[0] 
                bam_rdcount_file = f'{sample}.out'
                get_bam(sample, token_file, chrom, bam_start, bam_stop)
                bam_file = f'{sample}_{chrom}_{bam_start}_{bam_stop}.bam'
                run(f'/home/ec2-user/bam-readcount/bin -b 20 -q 20 -w 1 -l RNF145.bed -f GRCh38.d1.vd1.fa {bam_file} > {bam_rdcount_file}')
                # run(f'/Users/kcotto/Git/bin/bam-readcount -b 20 -q 20 -w 1 -l /Users/kcotto/Git/bin/RNF145.bed -f GRCh38.d1.vd1.fa {bam_file} > {bam_rdcount_file}')
                with open(bam_rdcount_file, 'r') as inputfile:
                    reader = csv.reader(inputfile, delimiter='\t')
                    for line in reader:
                        if line[0] != chrom:
                            continue
                        chr = line[0]
                        pos = line[1]
                        ref = line[2]
                        depth = line[3]
                        variants = {}
                        for item in line[4:]:
                            fields = item.split(':')
                            if int(fields[1]) > 0 and fields[0] != ref:
                                variants[fields[0]] = fields[1]
                                f.write(f'{sample}\t{chr}\t{pos}\t{ref}\t{depth}\t{fields[0]}\t{variants[fields[0]]}\n')
                                f.flush()

cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']

# for cohort in cohorts:
#     run_bam_readcount(cohort, 'all_TCGA_samples_wannotatedbeds.tsv')

if __name__ == '__main__':
    if len(sys.argv) == 1:
        for cohort in cohorts:
            call(f'./.local/bin/pipenv run python get_bam-readcount_TCGA.py {cohort} &')
    else:
        run_bam_readcount(sys.argv[1], 'all_TCGA_samples_wannotatedbeds.tsv')