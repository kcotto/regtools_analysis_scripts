# pylint: disable=missing-module-docstring,missing-function-docstring,invalid-name

import sys
import re
from typing import List

sample_counts = {
    'CHOL':	36,
    'DLBC':	36,
    'UCS':	56,
    'KICH':	65,
    'MESO':	77,
    'UVM':	77,
    'ACC':	79,
    'SKCM':	103,
    'THYM':	118,
    'GBM':	140,
    'READ':	149,
    'TGCT':	149,
    'ESCA':	151,
    'PAAD':	175,
    'PCPG':	175,
    'SARC':	252,
    'OV':	265,
    'KIRP':	284,
    'CESC':	294,
    'KIRC':	332,
    'LIHC':	368,
    'STAD':	372,
    'BLCA':	405,
    'COAD':	425,
    'PRAD':	482,
    'THCA':	485,
    'LUSC':	489,
    'HNSC':	491,
    'LGG':	498,
    'LUAD':	508,
    'UCEC':	533,
    'BRCA':	1022
}
def get_cohorts(file) -> List[str]:
    while True:
        match = re.match(r'^\+\ cohorts=\((.*?)\)', file.readline())
        if match:
            return match.group(1).split(' ')

def get_variant_counts(file) -> dict:
    while True:
        match = re.match(r'^\+\ wc\ -l\ (.*?)$', file.readline())
        if match:
            counts = {}
            while True:
                match = re.match(r'^\s*(\d+)\ all_splicing_variants_(\w+)\.bed', file.readline())
                if match:
                    counts[match.group(2)] = int(match.group(1))
                else:
                    return counts

def get_times(file) -> dict:
    times = {}
    while True:
        real_time = re.match(r'^real\s*(\d+)m(\d+.\d+)s$', file.readline())
        if real_time:
            user_time = re.match(r'^user\s*(\d+)m(\d+.\d+)s$', file.readline())
            sys_time = re.match(r'^sys\s*(\d+)m(\d+.\d+)s$', file.readline())
            times['real'] = int(real_time.group(1)) * 60 + float(real_time.group(2))
            times['user'] = int(user_time.group(1)) * 60 + float(user_time.group(2))
            times['sys'] = int(sys_time.group(1)) * 60 + float(sys_time.group(2))
            times['cpu'] = times['user'] + times['sys']
            return times

def get_row_count(file) -> int:
    while True:
        match = re.match(r'^\"Number\ of\ rows\ in\ data\.table\"$', file.readline())
        if match:
            row_count = re.match(r'^(\d+)$', file.readline())
            if row_count:
                return int(row_count.group(1))

variant_order = ['default', 'i50e5', 'E', 'I']

with open(sys.argv[1], 'r') as f:
    print('cohort,samples,variant,variant_count,real,user,sys,cpu')
    cohorts = get_cohorts(f)
    for cohort in cohorts:
        variants = get_variant_counts(f)
        for variant in variant_order:
            execution_times = get_times(f)
            print(f'{cohort},{sample_counts[cohort]},{variant},{variants[variant]},{execution_times["real"]},{execution_times["user"]},{execution_times["sys"]},{execution_times["cpu"]}')
