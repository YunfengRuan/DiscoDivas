python $PRScs \
 --ref_dir=${ref_dir} \
 --bim_prefix=${bim_prefix} \
 --sst_file=${work_dir}/formated.gwas/${dat}.chr${chr} \
 --n_gwas=${n_trn} \
 --chrom=${chr} \
 --phi=${phi} \
 --out_dir=${work_dir2}/posterior/PRS-CS.${dat}
