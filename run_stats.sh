#!/bin/bash

export PERL5LIB=~/vcftools_0.1.13/perl

tags=(E I default i50e5)
cohorts=(CHOL DLBC UCS)
# cohorts=(CHOL DLBC UCS KICH MESO UVM ACC SKCM THYM GBM READ TGCT ESCA PAAD PCPG SARC OV KIRP CESC KIRC LIHC STAD BLCA COAD PRAD THCA LUSC HNSC LGG LUAD UCEC BRCA)

for c in ${cohorts[@]}; do 
	mkdir -p ${c}/samples/
	mkdir -p ${c}/compare_junctions/hist/
	cd ${c}/samples/
	aws s3 cp s3://regtools-results-unstranded/${c}/ . --recursive
	for t in *.gz; do
		tar xzf ${t}
		rm ${t}
		rm -rf all*
		rm -rf compare_junctions
	done
	cd ..
	ls samples/ > dir_names.tsv
	for k in ${tags[@]}; do
		echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_${k}.bed
		for i in samples/TCGA*/; do
			j=${i##samples/}
			uniq ${i}/variants_${k}.bed | awk '{split($0, a, ","); if (length(a[2]) != 0) print a[1]"\n"a[2]; else print a[1]}' | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_${k}.bed
			#uniq ${i}/variants_${k}.bed | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_${k}.bed
			#Rscript --vanilla ~/compare_junctions_hist.R ${k}
		done
		# aws s3 cp compare_junctions/ s3://regtools-results-unstranded/${c}/compare_junctions/ --recursive
		# rm -rf compare_junctions/hist/junction_pvalues_$tag.tsv
	done
	Rscript --vanilla ~/compare_junctions_hist_default.R
	Rscript --vanilla ~/compare_junctions_hist_i50e5.R
	Rscript --vanilla ~/compare_junctions_hist_E.R
	Rscript --vanilla ~/compare_junctions_hist_I.R
	for i in all_*; do
		aws s3 cp ${i} s3://regtools-results-unstranded/${c}/
	done
	aws s3 cp compare_junctions/ s3://regtools-results-unstranded/${c}/compare_junctions/ --recursive
	cd ..
	rm -rf ${c}*
done
