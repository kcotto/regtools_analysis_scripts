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

def read_venn(filename):
    ding_ver = set()
    ding = set()
    ver = set()
    cur_set = None
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        reader.__next__()
        for line in reader:
            if line[0]:
                if line[0] == 'ding us ver':
                    cur_set = ding_ver
                elif line[0] == 'us ver':
                    cur_set = ver
                elif line[0] == 'ding us':
                    cur_set = ding
                else:
                    cur_set = None
            if cur_set is not None:
                cur_set.add(line[2])
    return ding_ver,ding,ver

def read_cancergenes(filename):
    result = set()
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        reader.__next__()
        for line in reader:
            result.add(line[0])
    return result

def read_gencode(filename):
    result = {}
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        for line in reader:
            if line[0].startswith('##'):
                continue
            if line[2] == 'exon':
                info = line[8]
                # if info.find('transcript_type "protein_coding"') == -1 or info.find('tag "basic"') == -1:
                if 'transcript_type "protein_coding"' not in info or \
                        'transcript_support_level "5"' in info or \
                        'transcript_support_level "4"' in info or \
                        'transcript_support_level "NA"' in info:     
                    continue
                chrom = line[0]
                strand = line[6]
                if strand == '+':
                    acceptor = line[3]
                    donor = line[4]
                else:
                    acceptor = line[4]
                    donor = line[3]
                acceptor_key = f'{chrom}_{acceptor}'
                donor_key = f'{chrom}_{donor}'
                if acceptor_key not in result:
                    result[acceptor_key] = 1   
                else:
                    result[acceptor_key] += 1
                if donor_key not in result:
                    result[donor_key] = -1   
                else:
                    result[donor_key] -= 1
    return result

def make_spliceai_bed(filename, gtf_dict, cancer_genes, DV, D, V):
    with open(filename, 'r') as input_file, open(f'{filename}_edited', 'w') as outfile:
        reader = csv.DictReader(input_file, delimiter='\t')
        old_header = reader.fieldnames.copy()
        link_header = old_header[-1]
        old_header.extend(['cancer_gene?','misplice,veridical match?', link_header])
        outfile.write('\t'.join(old_header) + '\n')
        for line in reader:
            spliceai = line['SpliceAI_raw']
            variant_junction = line['variant_junction_info']
            samples_field = line['variant_samples']
            strand = line['strand']
            genes = line['genes'].split(',')
            full_variant = line['variant_info']
            new_spliceai_field = 'NA'
            chrom = variant_junction.split('_')[0]
            if spliceai != 'NA':
                if '.' in spliceai:
                    spliceai_fields = spliceai.split('|')
                    DS_AG = spliceai_fields[2]
                    DS_AL = spliceai_fields[3]
                    DS_DG = spliceai_fields[4]
                    DS_DL = spliceai_fields[5]
                    DP_AG = int(spliceai_fields[6])
                    DP_AL = int(spliceai_fields[7])
                    DP_DG = int(spliceai_fields[8])
                    DP_DL = int(spliceai_fields[9])
                    samples_field = samples_field.replace(',','_')
                    variant_junction = variant_junction.replace(':', '_')
                    variant = int(variant_junction.split('-')[-1])
                    junc_end_1 = int(variant_junction.split('_')[1])
                    junc_end_2 = int(variant_junction.split('_')[2])
                    new_file = f'{samples_field}_{variant_junction}.bed'
                    P_AG = variant + DP_AG
                    P_AL = variant + DP_AL
                    P_DG = variant + DP_DG
                    P_DL = variant + DP_DL
                    red = '255,0,0'
                    blue = '0,0,255'
                    header = 'track name="SpliceAI sites" description="SpliceAI acceptor/donor positions" itemRgb="On"'
                    with open(new_file, 'w') as output_file:
                        fw = lambda a,b,c,d,e,f: output_file.write(f'{a}\t{b}\t{b}\t{c}\t{d}\t{e}\t{b}\t{b}\t{f}\n')
                        output_file.write(header + '\n')
                        fw(chrom, P_AG, 'acceptor_gain', DS_AG, strand, red)
                        fw(chrom, P_AL, 'acceptor_loss', DS_AL, strand, red) 
                        fw(chrom, P_DG, 'donor_gain', DS_DG, strand, blue) 
                        fw(chrom, P_DL, 'donor_loss', DS_DL, strand, blue)
                    new_spliceAI_tags = []
                    if strand == '+':
                        if junc_end_1 == P_DG or junc_end_1 == P_DL:
                            new_spliceAI_tags.append('novel donor match')
                        if junc_end_2 == P_AG or junc_end_2 == P_AL:
                            new_spliceAI_tags.append('novel acceptor match')
                    else:
                        if junc_end_2 == P_DG or junc_end_2 == P_DL:
                            new_spliceAI_tags.append('novel donor match')
                        if junc_end_1 == P_AG or junc_end_1 == P_AL:
                            new_spliceAI_tags.append('novel acceptor match')
                    AG_key = f'{chrom}_{P_AG}'
                    AL_key = f'{chrom}_{P_AL}'
                    DG_key = f'{chrom}_{P_DG}'
                    DL_key = f'{chrom}_{P_DL}'
                    keys = [AG_key, AL_key, DG_key, DL_key]
                    is_acceptor = False
                    is_donor = False
                    for key in keys:
                        if key in gtf_dict:
                            if gtf_dict[key] > 0:
                                is_acceptor = True
                            if gtf_dict[key] < 0:
                                is_donor = True
                    if is_acceptor:
                        new_spliceAI_tags.append('canonical acceptor match')
                    if is_donor:
                        new_spliceAI_tags.append('canonical donor match')
                    if new_spliceAI_tags:
                        new_spliceai_field = ','.join(new_spliceAI_tags)
            cancer_genes_results = []
            for gene in genes:
                if gene in cancer_genes:
                    cancer_genes_results.append('yes')
                else:
                    cancer_genes_results.append('no')
            new_cancergenes_field = ','.join(cancer_genes_results)
            file = filename.split('/')[-1]
            cancer_type = file.split('_')[0]
            paper_key = f'{cancer_type}:{chrom}:{full_variant}'
            if paper_key in DV:
                paper_result = 'misplice,veridical'
            elif paper_key in D:
                paper_result = 'misplice'
            elif paper_key in V:
                paper_result = 'veridical'
            else:
                paper_result = 'NA'
            new_line = [x for x in line.values() if x is not None]
            new_line[27] = new_spliceai_field
            link = line['IGV_link']
            del new_line[-1]
            new_line.append(new_cancergenes_field)
            new_line.append(paper_result)
            new_line.append(link)
            out_line = '\t'.join(new_line)
            outfile.write(out_line + '\n')


cancer_genes = read_cancergenes('Census_all.tsv')
DV,D,V = read_venn('venn_result22030.txt')
gtf = read_gencode('gencode.v29.annotation.gtf')

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
        make_spliceai_bed(file, gtf, cancer_genes, DV, D, V)
        bedfiles = glob.glob('*.bed')
        for bedfile in bedfiles: 
            run(
            f'aws s3 cp {bedfile} s3://regtools-igv-files/spliceAI/')
    os.chdir('..')
    if os.path.exists(f'{cohort}_igv_session'):
        shutil.rmtree(f'{cohort}_igv_session')
