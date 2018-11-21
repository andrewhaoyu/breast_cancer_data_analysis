file="aaa.impute2.gen.gz"


snp effectallele beta
rsd:adf:a:b A	0.93

awk '{if(NR>1) print $1}' pheno_include_all_relatives.def.txt > aaa

snplist="target.snplist"


for (( i=1; i<=22; i=i+1 ))
do
        echo "Chr_$i";
        file="${MGS_imputation_results_dir}/chr.${i}.segment.${j}.gen.gz"
        if [ -f $file ]
        then
                zcat $file | awk -v C=2 -v TF="MGS_${info_low_thr}_${info_high_thr}_file_for_clumping_snp_list" -f ~/bin/Element.Matching.awk | awk -v chr=$i '{printf chr"\t"$2"\t"$3"\t"$4"\t"$5; for(i=6; i<NF; i+=3) {if($(i+0) == 0 && $(i+1) == 0 && $(i+2) == 0) printf "\tNA"; else printf "\t"$(i+0)*2+$(i+1)*1+$(i+2)*0}; printf "\n"}' | gzip >> ${extracted_dosage_file}
        fi
        printf "\n"
done




plink.cmd = paste0('plink --score ', prs.model.file, ' no-sum no-mean-imputation --dosage ', genotype_file, ' noheader skip0=1 skip1=1 format=1 --fam ', fam.file)




zcat chr.1.dosage.gz | gzip > chr.all.dosage.gz
for i in {2..22}
do
	zcat chr.${i}.dosage.gz | gzip >> chr.all.dosage.gz
done

zcat chr.tesgt.gz | awk '{print $1"\t"$4}'