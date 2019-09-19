#!/bin/bash

export PERL5LIB=~/vcftools_0.1.13/perl

tags=(E I default i50e5)
# cohorts=(CHOL DLBC UCS)
cohorts=(CHOL DLBC UCS KICH MESO UVM ACC SKCM THYM GBM READ TGCT ESCA PAAD PCPG SARC OV KIRP CESC KIRC LIHC STAD BLCA COAD PRAD THCA LUSC HNSC LGG LUAD UCEC BRCA)

for c in ${cohorts[@]}; do 
	mkdir ${c}
	cd ${c}
	aws s3 cp s3://regtools-results-unstranded/${c}/ . --recursive
	for t in *.gz; do
		tar xzf ${t}
		rm ${t}
	done
	for i in TCGA*; do
		for k in ${tags[@]}; do
			$HOME/variants.sh ${i}/cse_identify_filtered_${k}.tsv ${i}/variants_${k}.bed
			#echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_${k}.bed
			#j=${i##samples/}
			#uniq ${i}/variants_${k}.bed | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_${k}.bed
		done
		#cut ${i}/*master.vcf -f 1-9 > ${i}/master_trimmed.vcf
		tar -czf ${i}.tar.gz ${i}
		aws s3 cp ${i}.tar.gz s3://regtools-results-unstranded/${c}/
	done

	#~/vcftools_0.1.13/bin/vcf-concat */*master_trimmed.vcf | ~/vcftools_0.1.13/bin/vcf-sort > all_variants_sorted.vcf
	cd ..
	rm -rf ${c}*
done

