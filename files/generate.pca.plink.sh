beta=./g1000.pca2/g1k_hm3_maf5_woamb_wolr.pca.weight
n=$(zless $beta | head -n 1 | awk '{print NF}')
echo $n

pfile=./pfile # or the path to your own plink format genotype data and its prefix
out=./plink.score # or the path to your score (PCA) data and its prefex

plink2 \
--score $beta 2 6 header-read zs no-mean-imputation cols='fid,nallele,dosagesum,scoresums' \
--score-col-nums 7-${n} \
--pfile  \
--rm-dup exclude-mismatch \
--out $out \
--memory 5000
