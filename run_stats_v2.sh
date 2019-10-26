
tags=(E I default i50e5)
cohorts=(CHOL)
#cohorts=(CHOL DLBC UCS KICH MESO UVM ACC SKCM THYM GBM READ TGCT ESCA PAAD PCPG SARC OV KIRP CESC KIRC LIHC STAD BLCA COAD PRAD THCA LUSC HNSC LGG LUAD UCEC BRCA)
#cohorts=(SKCM THYM GBM READ TGCT ESCA PAAD OV CESC LIHC STAD BLCA COAD PRAD THCA LUSC HNSC LGG LUAD UCEC BRCA)

for c in ${cohorts[@]}; do
	mkdir -p ${c}/samples/
	mkdir -p ${c}/compare_junctions/hist/
	cd ${c}/samples/
	aws s3 cp s3://regtools-results-unstranded/${c}/ . --recursive
	for t in *.gz; do
		tar xzf ${t}
		rm ${t}
	done
	cd ..
	ls samples/ > dir_names.tsv		
	for k in ${tags[@]}; do
		for i in samples/TCGA*; do
			bash /data/variants.sh ${i}/cse_identify_filtered_${k}.tsv ${i}/variants_${k}.bed
		done
		echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_${k}.bed
		for i in samples/TCGA*; do
			j=${i##samples/}
			uniq ${i}/variants_${k}.bed | awk '{split($0, a, ","); if (length(a[2]) != 0) print a[1]"\n"a[2]; else print a[1]}' | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_${k}.bed
		done
	done
	Rscript --vanilla /home/ec2-user/regtools/scripts/compare_junctions_hist.R default
	Rscript --vanilla /home/ec2-user/regtools/scripts/compare_junctions_hist.R i50e5
	Rscript --vanilla /home/ec2-user/regtools/scripts/compare_junctions_hist.R E
	Rscript --vanilla /home/ec2-user/regtools/scripts/compare_junctions_hist.R I
	Rscript --vanilla /home/ec2-user/regtools/scripts/filter_and_BH.R default
	Rscript --vanilla /home/ec2-user/regtools/scripts/filter_and_BH.R i50e5
	Rscript --vanilla /home/ec2-user/regtools/scripts/filter_and_BH.R E
	Rscript --vanilla /home/ec2-user/regtools/scripts/filter_and_BH.R I	
	for i in all_*; do
		aws s3 cp ${i} s3://regtools-results-unstranded/${c}/
	done
	for i in samples/*; do
		tar -czf ${i}.tar.gz ${i}
		printf 'aws s3 cp ${i}.tar.gz s3://regtools-results-unstranded/${c}/'
	done
	aws s3 cp compare_junctions/ s3://regtools-results-unstranded/${c}/compare_junctions/ --recursive
	cd ..
	rm -rf ${c}*
done
