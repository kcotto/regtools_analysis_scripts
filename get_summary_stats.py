import csv
import subprocess
import sys
import os
import shutil
import glob
import vcf
from statistics import mean

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

def call(cmd):
    subprocess.call(cmd, shell=True)

def create_gtex_set(filenname):
    gtex = set()
    with open(filenname, 'r') as gtex_file:
        next(gtex_file)
        junctions = gtex_file.readlines()
        for junction in junctions:
            junction = junction.split('\t')[0]
            gtex.add(junction)
    return gtex

def get_junction_counts(cohort, sample, gtex):
    output_file = f'{cohort}_junctioncounts_greaterthan5.tsv'
    output_file_gtex = f'{cohort}_junctioncounts_greaterthan5_gtexfiltered.tsv'
    total_junction_output_file = f'{cohort}_totaljunctions.tsv'
    total_counts = {'DA':0, 'D':0, 'A':0, 'NDA':0, 'N':0}
    gtex_counts = {'DA':0, 'D':0, 'A':0, 'NDA':0, 'N':0}
    annotated_bed = f'{sample}/{sample}_annotated.bed'
    write_header_junction_stats = not os.path.exists(output_file)
    if write_header_junction_stats:
        junction_stats = open(output_file, 'w')
        junction_stats.write(f'Sample\tCohort\tDA\tD\tA\tNDA\tN\tTotal Unique Junctions\n')
    else:
          junction_stats = open(output_file, 'a')
    write_header_junction_stats = not os.path.exists(output_file_gtex)
    if write_header_junction_stats:
        junction_stats_gtex = open(output_file_gtex, 'w')
        junction_stats_gtex.write(f'Sample\tCohort\tDA\tD\tA\tNDA\tN\tTotal Unique Junctions\n')
    else:
          junction_stats_gtex = open(output_file_gtex, 'a')
    write_header_total_junction = not os.path.exists(total_junction_output_file)
    if write_header_total_junction:
        total_count = open(total_junction_output_file, 'w')
        total_count.write(f'Sample\tTotal number of junctions\n')
    else:
        total_count = open(total_junction_output_file, 'a')
    with open(annotated_bed, 'r') as bedfile:
        total_junction_count = 0
        reader = csv.DictReader(bedfile, delimiter='\t')
        for line in reader:
            score = int(line['score'])
            chrom = line['chrom']
            junction_start = str(int(line['start']) + 1)
            junction_end = str(int(line['end']) - 1)
            regtools_key = f'{chrom}_{junction_start}_{junction_end}'
            total_junction_count += score
            if regtools_key not in gtex and score >= 5:  
                gtex_counts[line['anchor']] += 1
            if score >= 5:
                total_counts[line['anchor']] += 1
        total_unique_junctions_gtex = sum(gtex_counts.values())
        total_unique_junctions = sum(total_counts.values())
        junction_stats_gtex.write(f"{sample}\t{cohort}\t{gtex_counts['DA']}\t{gtex_counts['D']}\t{gtex_counts['A']}\t{gtex_counts['NDA']}\t{gtex_counts['N']}\t{total_unique_junctions_gtex}\n")
        junction_stats.write(f"{sample}\t{cohort}\t{total_counts['DA']}\t{total_counts['D']}\t{total_counts['A']}\t{total_counts['NDA']}\t{total_counts['N']}\t{total_unique_junctions}\n")
        total_count.write(f'{sample}\t{cohort}\t{total_junction_count}\n')
        junction_stats.close()
        total_count.close()


def get_VAF(cohort, sample_name, vcf_file):
    output_file = f'{cohort}_VAFs.tsv'
    write_header = not os.path.exists(output_file)
    if write_header:
        vaf_summary_file = open(output_file, 'w')
        vaf_summary_file.write(f'chr\tposition\tref\talt\tmean_VAF\tnumber of callers\tsample\n')
    else:
        vaf_summary_file = open(output_file, 'a')
    with open(f'{cohort}_VAFs.tsv', 'a') as vaf_summary_file:
        order = {'A':0, 'C':1, 'G':2, 'T':3}
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        for record in vcf_reader:
            alt_allele = record.ALT[0].sequence
            ref_allele = record.REF
            chr = record.CHROM
            pos = record.POS
            num_variant_callers = len(record.INFO['SF'])
            VAF_values = []
            for sample in record.samples:
                if 'TUMOR' in sample.sample:
                    print(sample.sample)
                    if hasattr(sample.data, 'AF') and sample.data.AF is not None:
                        AF = sample.data.AF * 100.0
                        VAF_values.append(AF)
                    elif hasattr(sample.data, 'BCOUNT') and sample.data.BCOUNT is not None:
                        allele_count = sample.data.BCOUNT[order[alt_allele]]
                        total_count = sample.data.DP
                        VAF = (allele_count/total_count) * 100.0
                        VAF_values.append(VAF)
                    elif hasattr(sample.data, 'FREQ') and sample.data.FREQ is not None:
                        VAF = float(sample.data.FREQ[:-1])
                        VAF_values.append(VAF)
                    elif hasattr(sample.data, 'AD') and sample.data.AD is not None and len(sample.data.AD) == 2:
                        VAF = (sample.data.AD[1]/sample.data.DP) * 100.0
                        VAF_values.append(VAF)
            if len(VAF_values) != num_variant_callers:
                print(f'Error occured. Expected {num_variant_callers} values but only have {len(VAF_values)}.')
            mean_VAF = mean(VAF_values)
            vaf_summary_file.write(f'{chr}\t{pos}\t{ref_allele}\t{alt_allele}\t{mean_VAF}\t{num_variant_callers}\t{sample_name}\n')

def get_regtools_variant_numbers(cohort, sample):
    output_file = f'{cohort}_regtools_variantcounts.tsv'
    write_header = not os.path.exists(output_file)
    if write_header:
        outfile = open(output_file, 'w')
        outfile.write(f'cohort\tsample\ttag\tvariant count\n')
    else:
        outfile = open(output_file, 'a')
    files = glob.glob(f'{sample}/*.tsv')
    for file in files:
        tag = file.split('_')[-1].split('.')[0]
        with open(file, 'r') as inputfile:
            reader = csv.DictReader(inputfile, delimiter = '\t')
            variant_list = set([])
            for line in reader:
                variant_list.add(line['variant_info'])
            variant_count = len(variant_list)
            outfile.write(f'{cohort}\t{sample}\t{tag}\t{variant_count}\n')


def get_sample_files(cohort, filename, gtex):
    with open(filename, 'r') as inputfile:
        samples = inputfile.readlines()
        for item in samples:
            fields = item.split('\t')
            if fields[1] == cohort:
                sample = fields[0] 
                run(f'aws s3 cp {results_files}/{cohort}/{sample}.tar.gz .')
                annotated_bed = f'{sample}/{sample}_annotated.bed'
                vcf_file = f'{sample}/{sample}_master.vcf'
                default_regtools = f'{sample}/cse_identify_filtered_default.tsv'
                i50e5_regtools = f'{sample}/cse_identify_filtered_i50e5.tsv'
                E_regtools = f'{sample}/cse_identify_filtered_E.tsv'
                I_regtools = f'{sample}/cse_identify_filtered_I.tsv'
                files = [annotated_bed, vcf_file, default_regtools, i50e5_regtools, E_regtools, I_regtools]
                for file in files:
                    run(f'tar xzf {sample}.tar.gz {file}')
                # get_VAF(cohort, sample, vcf_file)
                get_junction_counts(cohort, sample, gtex)
                get_regtools_variant_numbers(cohort, sample)
                shutil.rmtree(sample, ignore_errors=True)
                os.remove(f'{sample}.tar.gz')
        files_to_save = glob.glob(f'{cohort}*.tsv')
        for file in files_to_save:
            run(f'aws s3 cp {file} {results_files}/junction_counts/')

results_files = 's3://regtools-results-unstranded'
cohorts=['CHOL', 'DLBC', 'UCS', 'KICH', 'MESO', 'UVM', 'ACC', 'SKCM',
          'THYM', 'GBM', 'READ', 'TGCT', 'ESCA', 'PAAD', 'PCPG', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA', 'COAD',
          'PRAD', 'THCA', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']

# for cohort in cohorts:
#     gtex = create_gtex_set('master_GTEx_b38_mean_sd.tsv')
#     get_sample_files(cohort, 'primarysolidtumors_metadata.tsv', gtex)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        for cohort in cohorts:
            call(f'python3 get_summary_stats.py {cohort} &')
    else:
        gtex = create_gtex_set('master_GTEx_b38_mean_sd.tsv')
        get_sample_files(sys.argv[1], 'primarysolidtumors_metadata.tsv', gtex)