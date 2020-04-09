import glob
import csv
import os


tag = 'i50e5'
default_files = glob.glob(f'/Users/kcotto/Desktop/TCGA/regtools/*{tag}*')
root = f'/Users/kcotto/Desktop/TCGA/variant_analysis/{tag}'

results = {}

summary_file = f'{root}/variantsummary.tsv'
with open(summary_file, 'w') as summary:
    summary.write(f'Cohort\tTotal variants\tVariants w/ single junction\tVariants w/ multiple junctions of same type\tVariants w/ multiple junctions of different types\n')


    for file in default_files:
        cohort = file.split('/')[-1].split('_')[0]
        variant_w_onejunctiontype_outfile = f'{root}/{cohort}_{tag}_variants_w_mult_of_onejunctype.txt'
        variant_w_multjunctiontype_outfile = f'{root}/{cohort}_{tag}_variants_w_mult_of_multjunctype.txt'
        with open(file, 'r') as inputfile, \
            open(variant_w_onejunctiontype_outfile, 'w') as one_junc_outfile, \
            open(variant_w_multjunctiontype_outfile, 'w') as mult_junc_outfile:
            reader = csv.DictReader(inputfile, delimiter='\t')
            one_junc_outfile.write(f'Variants with multiple junctions of the same type\n')
            mult_junc_outfile.write(f'Variants with multiple junctions of different types\n')
            for line in reader:
                variant = line['variant_info']
                anchor = line['anchor']
                if variant not in results:
                    results[variant] = {anchor:1}
                else:
                    if anchor not in results[variant]:
                        results[variant][anchor] = 1
                    else:
                        results[variant][anchor] += 1
        
            number_of_unique_variants = len(results.keys())
            variants_w_singlejunction = 0
            variants_w_mult_onetypeofjunction = 0
            variants_w_mult_multtypeofjunction = 0
            for item in results.keys():
                anchors = list(results[item].items())
                if len(results[item].keys()) == 1 and anchors[0][1] == 1:
                    variants_w_singlejunction += 1
                elif len(results[item].keys()) == 1 and anchors[0][1] > 1:
                    variants_w_mult_onetypeofjunction += 1
                    one_junc_outfile.write(f'{item}\n\t{anchors[0][0]}\t{anchors[0][1]}\n')
                else:
                    variants_w_mult_multtypeofjunction += 1
                    mult_junc_outfile.write(f'{item}\n')
                    for anchor in anchors:
                        mult_junc_outfile.write(f'\t{anchor[0]}\t{anchor[1]}\n')

            
            summary.write(f'{cohort}\t{number_of_unique_variants}\t{variants_w_singlejunction}\t{variants_w_mult_onetypeofjunction}\t{variants_w_mult_multtypeofjunction}\n')

