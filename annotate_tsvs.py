import fnmatch
import subprocess
import glob
import csv
import sys
import os
import shutil
from pathlib import Path
import argparse

### created from annotate_spliceai_gtex, modify_tsvs, variant_in_junction, and create_mutually_exclusive

def create_mut_exclusive(cohort):
    default = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_default.tsv'
    i50e5 = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_i50e5.tsv'
    E = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_E.tsv'
    I = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_I.tsv'

    default_set = set()
    i50e5_set = set()

    new_default_file = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_mutex_default.tsv'
    new_i50e5_file = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_mutex_i50e5.tsv'
    new_E_file = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_mutex_E.tsv'
    new_I_file = f'{cohort}_junction_pvalues_significant_0.05_filtered_BH_mutex_I.tsv'

    with open(default, 'r') as d, open(new_default_file, 'w') as outfile:
        line_count = 0
        reader = csv.reader(d, delimiter='\t')
        for line in reader:
            line_count += 1
            if line_count == 1:
                outfile.write('\t'.join(line) + '\n')
            if line_count != 1:
                variant_junction = line[9]
                default_set.add(variant_junction)
                outfile.write('\t'.join(line) + '\n')

    with open(i50e5, 'r') as i50e5_file, open(new_i50e5_file, 'w') as outfile:
        line_count = 0
        reader = csv.reader(i50e5_file, delimiter='\t')
        for line in reader:
            line_count += 1
            if line_count == 1:
                outfile.write('\t'.join(line) + '\n')
            if line_count != 1:
                variant_junction = line[9]
                if variant_junction not in default_set:
                    outfile.write('\t'.join(line) + '\n')
                    i50e5_set.add(variant_junction)
                else:
                    print(f'This record ({variant_junction}) exists in the default results file.')

    combined_set = default_set.union(i50e5_set)

    with open(E, 'r') as E_file, open(new_E_file, 'w') as outfile:
        reader = csv.reader(E_file, delimiter='\t')
        for line in reader:
            line_count += 1
            if line_count == 1:
                outfile.write('\t'.join(line) + '\n')
            if line_count != 1:
                variant_junction = line[9]
                if variant_junction not in combined_set:
                    outfile.write('\t'.join(line) + '\n')
                else:
                    print(f'This record ({variant_junction}) exists in either the default or i50e5 results file.')

    with open(I, 'r') as I_file, open(new_I_file, 'w') as outfile:
        reader = csv.reader(I_file, delimiter='\t')
        for line in reader:
            line_count += 1
            if line_count == 1:
                outfile.write('\t'.join(line) + '\n')
            if line_count != 1:
                variant_junction = line[9]
                if variant_junction not in combined_set:
                    outfile.write('\t'.join(line) + '\n')
                else:
                    print(f'This record ({variant_junction}) exists in either the default or i50e5 results file.')

def matched_annotate_GTEx(regtools_input_file, gtex_file, outputfile):
    with open(regtools_input_file, 'r') as regtools, open(gtex_file, 'r') as gtex:
        gtex_values = {}
        gtex_reader = csv.reader(gtex, delimiter='\t')
        regtools_reader = csv.reader(regtools, delimiter='\t')
        next(gtex_reader)
        header = next(regtools_reader)
        match_count = 0
        with open(outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\tmatched_GTEx_mean\tmatched_GTEx_sd\n')
            for line in gtex_reader:
                mean = line[2]
                stdev = line[3]
                gtex_key = line[0]
                gtex_values[gtex_key] = (mean, stdev)
            for line in regtools_reader:
                chrom = line[4]
                junction_start = str(int(line[5]) + 1)
                junction_end = str(int(line[6]) - 1)
                regtools_key = f'{chrom}_{junction_start}_{junction_end}'
                output_line = '\t'.join(line)
                if regtools_key in gtex_values:
                    match_count += 1
                    print(f'match for {output_line}')
                    outfile.write(f'{output_line}\t{gtex_values[regtools_key][0]}\t{gtex_values[regtools_key][1]}\n')
                else:
                    outfile.write(f'{output_line}\tNA\tNA\n')
    outfile.close()
    print(match_count)


def master_annotate_GTEx(regtools_input_file, gtex_file, outputfile):
    with open(regtools_input_file, 'r') as regtools, open(gtex_file, 'r') as gtex:
        gtex_values = {}
        gtex_reader = csv.reader(gtex, delimiter='\t')
        regtools_reader = csv.reader(regtools, delimiter='\t')
        next(gtex_reader)
        header = next(regtools_reader)
        match_count = 0
        with open(outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\ttotal_GTEx_mean\ttotal_GTEx_sd\n')
            for line in gtex_reader:
                mean = line[2]
                stdev = line[3]
                gtex_key = line[0]
                gtex_values[gtex_key] = (mean, stdev)
            for line in regtools_reader:
                chrom = line[4]
                junction_start = str(int(line[5]) + 1)
                junction_end = str(int(line[6]) - 1)
                regtools_key = f'{chrom}_{junction_start}_{junction_end}'
                output_line = '\t'.join(line)
                if regtools_key in gtex_values:
                    match_count += 1
                    print(f'match for {output_line}')
                    outfile.write(f'{output_line}\t{gtex_values[regtools_key][0]}\t{gtex_values[regtools_key][1]}\n')
                else:
                    outfile.write(f'{output_line}\tNA\tNA\n')
    print(match_count)

def annotate_spliceAI(regtools_input_file, spliceAI_annotationfile, outputfile):
    with open(regtools_input_file, 'r') as regtools, open(spliceAI_annotationfile, 'r') as spliceAI:
        spliceAI_reader = csv.reader(spliceAI, delimiter='\t')
        regtools_reader = csv.reader(regtools, delimiter='\t')
        spliceAI_values = {}
        for line in spliceAI_reader:
            li = str(line[0]).strip()
            if not li.startswith('#'):
                info_field = line[7]
                if 'SpliceAI' in info_field:
                        pattern = 'SpliceAI*'
                        extract_SF = line[7].split(';')
                        matching = fnmatch.filter(extract_SF, pattern)
                        spliceAI_variant_key = f'{line[0]}:{line[1]}'
                        spliceAI_string = matching[0]
                        if spliceAI_string.__contains__(','):
                            split_values = spliceAI_string.split(',')
                            spliceAI_values[spliceAI_variant_key] = split_values[0]
                        else:
                            spliceAI_values[spliceAI_variant_key] = spliceAI_string
        header = next(regtools_reader)
        with open (outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\tSpliceAI_raw\tSpliceAI_match\n')
            for line in regtools_reader:
                chrom = line[1].split(':')[0]
                variant = line[1].split('-')[-1]
                regtools_variant_key = f'{chrom}:{variant}'
                output_line = '\t'.join(line)
                if regtools_variant_key in spliceAI_values:
                    junction_start_match = False
                    junction_end_match = False
                    junction_start = int(line[5])
                    junction_end = int(line[6])
                    spliceAI_info = spliceAI_values[regtools_variant_key]
                    spliceAI_split = spliceAI_info.split('|')
                    positions = spliceAI_split[-4:]
                    for position in positions:
                        if position is not '.':
                            genomic_location = int(variant) + int(position)
                            if junction_start == genomic_location:
                                junction_start_match = True
                            if junction_end == genomic_location:
                                junction_end_match = True
                    if junction_start_match and junction_end_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction start and end match\n')
                    elif junction_start_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction start match\n')
                    elif junction_end_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction end match\n')
                    else:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tNA\n')
                else:
                    outfile.write(f'{output_line}\tNA\tNA\n')


def variant_in_junction(regtools_input_file, outputfile):
    with open(regtools_input_file, 'r') as regtools:
        regtools_reader = csv.reader(regtools, delimiter='\t')
        header = next(regtools_reader)
        with open(outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\tvariant_in_junction?\n')
            for line in regtools_reader:
                output_line = '\t'.join(line)
                variant = int(line[1].split('-')[-1])
                junction_start = int(line[5]) + 1
                junction_end = int(line[6]) - 1
                if junction_start <= variant <= junction_end:
                    outfile.write(f'{output_line}\tyes\n')
                else:
                    outfile.write(f'{output_line}\tno\n')
        outfile.close()



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

def read_vep(filename):
    vep = set()
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=',')
        for line in reader:
            vep.add(line[0])
    return(vep)

            
def read_gencode(filename, gtf_filter=False):
    result = {}
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter='\t')
        for line in reader:
            if line[0].startswith('#'):
                continue
            if line[2] == 'exon':
                info = line[8]
                if gtf_filter:
                    if info.find('transcript_type "protein_coding"') == -1 or info.find('tag "basic"') == -1:
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

def make_spliceai_bed(filename, cohort, gtf_dict, cancer_genes, DV, D, V, vep):
    with open(filename, 'r') as input_file, open(f'{filename}_edited', 'w') as outfile:
        reader = csv.DictReader(input_file, delimiter='\t')
        old_header = reader.fieldnames.copy()
        # link_header = old_header[-1]
        # old_header = old_header[:-1]
        # old_header.extend(['cancer_gene?','misplice,veridical match?', 'VEP splicing?', link_header])
        old_header.extend(['cancer_gene?','misplice,veridical match?', 'VEP splicing?'])
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
                if '|.' in spliceai:
                    continue
                else:
                    spliceai_fields = spliceai.split('|')
                    # DS_AG = spliceai_fields[2]
                    # DS_AL = spliceai_fields[3]
                    # DS_DG = spliceai_fields[4]
                    # DS_DL = spliceai_fields[5]
                    DP_AG = int(spliceai_fields[6])
                    DP_AL = int(spliceai_fields[7])
                    DP_DG = int(spliceai_fields[8])
                    DP_DL = int(spliceai_fields[9])
                    samples_field = samples_field.replace(',','_')
                    variant_junction = variant_junction.replace(':', '_')
                    variant = int(variant_junction.split('-')[-1])
                    junc_end_1 = int(variant_junction.split('_')[1])
                    junc_end_2 = int(variant_junction.split('_')[2])
                    # new_file = f'{cohort}_{variant_junction}.bed'
                    P_AG = variant + DP_AG
                    P_AL = variant + DP_AL
                    P_DG = variant + DP_DG
                    P_DL = variant + DP_DL
                    # red = '255,0,0'
                    # blue = '0,0,255'
                    # header = 'track name="SpliceAI sites" description="SpliceAI acceptor/donor positions" itemRgb="On"'
                    # with open(new_file, 'w') as output_file:
                    #     fw = lambda a,b,c,d,e,f: output_file.write(f'{a}\t{b}\t{b}\t{c}\t{d}\t{e}\t{b}\t{b}\t{f}\n')
                    #     output_file.write(header + '\n')
                    #     fw(chrom, P_AG, 'acceptor_gain', DS_AG, strand, red)
                    #     fw(chrom, P_AL, 'acceptor_loss', DS_AL, strand, red) 
                    #     fw(chrom, P_DG, 'donor_gain', DS_DG, strand, blue) 
                    #     fw(chrom, P_DL, 'donor_loss', DS_DL, strand, blue)
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
            paper_key = f'{cancer_type}:{full_variant}'
            if paper_key in DV:
                paper_result = 'misplice,veridical'
            elif paper_key in D:
                paper_result = 'misplice'
            elif paper_key in V:
                paper_result = 'veridical'
            else:
                paper_result = 'NA'
            if paper_key in vep:
                vep_result = 'yes'
            else:
                vep_result = 'NA'
            new_line = [x for x in line.values() if x is not None]
            new_line[27] = new_spliceai_field
            #link = line['IGV_link']
            # del new_line[-1]
            new_line.append(new_cancergenes_field)
            new_line.append(paper_result)
            new_line.append(vep_result)
            # new_line.append(link)
            out_line = '\t'.join(new_line)
            outfile.write(out_line + '\n')

cohort_mappings = {}

with open('/Users/kcotto/Projects/regtools/TCGA_paper/cohort_build_gtex_mapping.txt', 'r') as mapping_file:
    reader = csv.reader(mapping_file, delimiter='\t')
    next(reader)
    for line in reader:
        cohort_mappings[line[0]] = line[1], line[2]


cohorts = ['HCC', 'SCLC']

# cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
#          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
#          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
#          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'BRCA']

# cohorts = ['OSCC']


significant = 'junction_pvalues_significant_0.05_filtered_BH_mutex_'

tags = ['default', 'i50e5', 'E', 'I']

print('Reading in Cancer genes')
cancer_genes = read_cancergenes('Census_all.tsv')
print('Reading in paper data')
DV,D,V = read_venn('venn_result22030.txt')
print('Reading in GTF')
gtf = read_gencode('/Users/kcotto/Desktop/regtools_filtered_files/Homo_sapiens.GRCh37.87.gtf', gtf_filter=True)
print('Reading in VEP data')
vep = read_vep('/Users/kcotto/Desktop/TCGA/VEP/all_vep_splice_variants_sampleremoved.txt')



for cohort in cohorts:
    os.chdir('/Users/kcotto/Desktop/regtools_filtered_files/')
    create_mut_exclusive(cohort)
    for tag in tags:
        regtools_file = f'/Users/kcotto/Desktop/regtools_filtered_files/{cohort}_{significant}{tag}.tsv'
        spliceAI_file = f'/Users/kcotto/Projects/regtools/TCGA_paper/spliceAI/VCFs_b{cohort_mappings[cohort][0]}/{cohort}_sig_variants.vcf.out'
        matched_gtex_file = f'/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices_b{cohort_mappings[cohort][0]}_avg/{cohort_mappings[cohort][1]}'
        master_gtex_file = f'/Users/kcotto/Projects/regtools/TCGA_paper/master_GTEx_b{cohort_mappings[cohort][0]}_mean_sd.tsv'
        base_file = '.'.join(regtools_file.split('.')[0:2])
        tmp_file = f'{base_file}_tmp.tsv'
        tmp_file2 = f'{base_file}_tmp2.tsv'
        final_output_file = f'{base_file}_gtex_spliceai.tsv'
        tmp_file3 = f'{base_file}_tmp3.tsv'
        matched_annotate_GTEx(regtools_file, matched_gtex_file, tmp_file)
        master_annotate_GTEx(tmp_file, master_gtex_file, tmp_file2)
        annotate_spliceAI(tmp_file2, spliceAI_file, final_output_file)
        variant_in_junction(final_output_file, tmp_file3)
        make_spliceai_bed(tmp_file3, cohort, gtf, cancer_genes, DV, D, V, vep)
        patterns = ('*tmp*.tsv', '*gtex_spliceai.tsv')
        files_to_rm = []
        for pattern in patterns:
            files_to_rm.extend(glob.glob(pattern))
        for file in files_to_rm:
            os.remove(file)
