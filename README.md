# DiscoDivas

## Disambiguation
This page is for "a Genetic ***Dis***tance-assisted PRS ***Co***mbination Pipeline for ***Div***erse Genetic ***A***ncestrie***s***, an R script for calculating PRS for diverse populations, especially admixed populations. 

For the album of disco music, please refer to this [spotify page](https://open.spotify.com/album/511A9pZXN4GhblH021gdDS)

## Required R packages
- data.table
- dplyr
- stringr
- rio
- optparse

## Command line to run DiscoDivas.R
Here is an example command line to run DiscoDivas.R

```
Rscript DiscoDivas.R -m ukbb.pca.med.1000G.tsv \
-p ukbb.pca.txt.gz \
--prs.list afr.prs.tsv.gz,${wd}/sas.prs.tsv.gz,eas.prs.tsv.gz,eur.prs.tsv.gz \
-s IID,PRS \
-A 1,0.9,1,1 \
--regress.PCA T \
-o combine.linear.prs
```
- `-m ukbb.pca.med.1000G.tsv`
  
  `ukbb.pca.med.1000G.tsv`can be found in the `files`. **Please notice that the order of the list after `-A` and `--prs.list` should has the same the order as the order of cohorts/ populations in this file.** User should replace the file by other file that contains the median value of top PC of the validation cohorts actually being used in the analysis.
  
  The format should be: 
  - must have a header; does't contain row name or row index;
  - the first column is the name of cohorts or populations; the seconds column and other columns after that are PCs
```
  POP	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
  AFR	-8054.5169	855.90292	-956.67184	-13544.445	4926.9673	3771.2851	812.396557	-995.90085	-1181.98158	-1144.22429
  ASN	14494.305	-6036.1673	4714.5286	-16009.741	5103.6832	4492.8499	1036.56551	-2046.52376	-1127.46504	-1249.9723
  CHN	18201.946	10558.105	-1619.620231	-13121.892	6257.718	4450.8678	1118.970433	-2033.12265	-797.63404	-1553.94603
  EUR	14634.677	-11959.621	-3110.3896	-11696.802	5096.9267	4264.5863	1074.55081	-1968.98792	-1316.456079	-1089.85971
  ```
    
- `-p ukbb.pca.txt.gz`

  `ukbb.pca.txt.gz` should be replaced by user's file contains the PC of the testing individual. It can be plain text file or `.gz` file. **Please notice that in all the list, items should be separated by `,` and no space should be inserted.**

  The format is:
  - must have header; does't contain row name or row index;
  - the first column is the individual ID; the seconds column and other columns after that are PCs

- `--prs.list afr.prs.tsv.gz,${wd}/sas.prs.tsv.gz,eas.prs.tsv.gz,eur.prs.tsv.gz`

  The list of PRS files of the testing cohort that are based on the weight fine-tuned in each of the validation cohorts. They can be plain text file or `.gz` file. The format is:
  - must have header; does't contain row name or row index;
  - must have a column of individual ID and a column of PRS
  - The individual ID in the PRS files should be the same as the individual ID in the PCA file

- `-s IID,PRS`
  The column names of individual ID and PRS in the PRS files. If the column names are the different in different files, please list them in the same order of `prs.list` : `ID1,PRS1,ID2,PPRS2,ID3,PRS3,ID4,PRS4`

- `-A 1,0.9,1,1`

  The list of shrinkage parameter for each PRS based on their quality. The default is a list of 1

- `--regress.PCA T`

  Choose whether combining PRS after correcting the ancestry information or not. Default and recommanded value is TRUE/T

- `-o`

  The prefix of output file. The format of output file is `.tsv.gz`


## Supporting files
Please visit `files` for the suppporting files







