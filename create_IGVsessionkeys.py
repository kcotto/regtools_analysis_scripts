import yaml
import subprocess
import glob
import csv
import sys
import os
import requests
import shutil
import json
from pathlib import Path
import argparse

uid = os.getuid()
gid = os.getgid()
cwd = os.getcwd()


# input_parser = argparse.ArgumentParser(
#     description="Create IGV sessions for TCGA data",
# )
# input_parser.add_argument(
#     'yaml_file',
#     help="Yaml file with inputs",
# )

# args = input_parser.parse_args()


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)


# def read_file(filename):
#     with open(filename) as f:
#         return f.read()


# def run_docker_image_as_current_user(image_name_and_all_args: str, stdout=None):
#     docker_args = image_name_and_all_args.split()
#     print(f'Running docker container command "{docker_args[0]} {docker_args[1]}"')
#     run(f'docker run --user {uid}:{gid} -v {cwd}:/manifests {image_name_and_all_args}')


# def get_bam(sample, token_file, chrom, bam_window_start, bam_window_end, variant_junction):
#     print('Obtaining bam manifest file')
#     bam_url = 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22BAM%22%5D%7D%7D%2C%7B%22op%22%3A%22AND%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%2C%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.submitter_id%22%2C%22value%22%3A%5B%22{0}%22%5D%7D%7D%5D%7D%5D%7D%5D%7D%5D%7D&query=files.data_format%20in%20%5B%22BAM%22%5D%20and%20files.experimental_strategy%20in%20%5B%22RNA-Seq%22%5D%20AND%20cases.project.program.name%20IN%20%5BTCGA%5D%20and%20cases.samples.submitter_id%20IN%20%5B%22P{0}&return_type=manifest'.format(sample)
#     response = requests.get(bam_url)
#     if response.status_code != 200:
#         raise Exception(sample, f'Failed to download bam manifest file: HTTP Status Code: {response.status_code}')
#     manifest_contents = response.content.decode('utf-8').strip()
#     if not manifest_contents.startswith('id\t'):
#         raise Exception(sample, f'bam manifest did not have expected format')

#     bam_manifest_path = Path(f'{sample}_bam_manifest.txt')
#     bam_manifest_path.write_text(manifest_contents)
#     # run_docker_image_as_current_user(
#     #     f'mgibio/gdc-client gdc-client download -m /manifests/{bam_manifest_path} -t /manifests/gdc-user-token.txt -d /manifests/igv_session')
#     run(f'aws s3 cp {token_file} .')
#     with open(bam_manifest_path, 'r') as manifest_file:
#         reader = csv.reader(manifest_file, delimiter='\t')
#         next(manifest_file)
#         for line in reader:
#             bam_id = line[0]
#             bam_file = f'{sample}_{variant_junction}.bam'
#             token_file = "gdc-user-token.txt"

#             file_id = f'{bam_id}'

#             data_endpt = "https://api.gdc.cancer.gov/slicing/view/{}".format(file_id)

#             with open(token_file, "r") as token:
#                 token_string = str(token.read().strip())

#             params = {"region": [f'{chrom}:{bam_window_start}-{bam_window_end}']}

#             response = requests.post(data_endpt,
#                                      data=json.dumps(params),
#                                      headers={
#                                          "Content-Type": "application/json",
#                                          "X-Auth-Token": token_string
#                                      })


#             with open(bam_file, "wb") as output_file:
#                 output_file.write(response.content)
#             run(f'samtools index {bam_file}')


# yamlfile = args.yaml_file

# with open(yamlfile) as f:
#     data = yaml.load(f, Loader=yaml.SafeLoader)

# file_location = data['file_location']
# cohorts = data['cohort']
# token_file = data['token_file']
# results_files = data['results_files']
# igv_session_location = data['igv_session_location']
# igv_file_location = data['igv_file_location']
# cluster_location = data['cluster_location']
# bam_location = data['bam_location']
results_files = 's3://regtools-results-unstranded'
cohorts=['SKCM',
          'GBM', 'READ', 'ESCA', 'PAAD', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']

with open('output.txt', 'w') as output:
    for cohort in cohorts:
        if os.path.exists(f'{cohort}_igv_session'):
            shutil.rmtree(f'{cohort}_igv_session')
        os.mkdir(f'{cohort}_igv_session')
        os.chdir(f'{cohort}_igv_session')

        run(
            f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/junction_pvalues_significant_0.05_filtered_BH_default.tsv .')
        run(
            f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/junction_pvalues_significant_0.05_filtered_BH_i50e5.tsv .')
        run(
            f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/junction_pvalues_significant_0.05_filtered_BH_E.tsv .')
        run(
            f'aws s3 cp {results_files}/{cohort}/compare_junctions2/hist/junction_pvalues_significant_0.05_filtered_BH_I.tsv .')

        # {cohort}|{tag(default, etc)}|{junction_samples}|{chrom}|{start}|{end}|{variant_junction_info}
        files = glob.glob('*junction_pvalues_significant_0.05_filtered_BH*.tsv')

        for i in files:
            tag = i.split('_')[-1].split('.')[0]
            with open(i, 'r') as result_file:
                reader = csv.DictReader(result_file, delimiter='\t')
                for line in reader:
                    key = f"{cohort}|{tag}|{line['junction_samples']}|{line['chrom']}|{line['start']}|{line['end']}|{line['variant_junction_info']}"
                    output.write(f'{key}\n')
        os.chdir('..')
        if os.path.exists(f'{cohort}_igv_session'):
            shutil.rmtree(f'{cohort}_igv_session')