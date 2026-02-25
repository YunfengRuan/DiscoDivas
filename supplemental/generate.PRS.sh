# To generate PRS based on genotype data and PRS model
# genotype: $dosage
# PRS model: $beta

n=$(zcat $beta | head -n 1 | awk '{print NF}')
echo $n

plink2=/medpop/esp2/yruan/tools/plink2

$plink2 \
--pfile ${dosage} \
--memory 6000 \
--score <(zcat $beta | sed 's|NA|0|g' | sed 's|NaN|0|g' | sed 's|Inf|0|g') \
1 2 header-read  list-variants-zs \
cols='fid,nallele,dosagesum,scoresums' \
--score-col-nums 4-$n \
--rm-dup  force-first \
--out ${out}.chr${chr}

