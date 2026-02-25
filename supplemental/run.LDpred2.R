library(bigsnpr)
library(data.table)
library(dplyr)
library(parallel)
library(parallelly)
library(optparse)
library(rio)

# measure running time
start_time <- Sys.time()

# parameters
args <- commandArgs(trailingOnly = TRUE)
pop <- tolower(args[1])
sst_file <- args[2]
out <- args[3]
trg_bim <- args[4]

# file directories
ref_bim <- paste0("/stanley/huang_lab/home/yruan/ldpred2.ref/g1000_", pop, "_hm3.bim")
tmp <- tempfile(tmpdir = paste0(out, "-tmp"))
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

pst_file <-paste0(out, ".tsv.gz")

# read summary stats, which has been formatted 

sumstats <- bigreadr::fread2(sst_file)
if (dim(sumstats)[2] == 9) {
colnames(sumstats) <- c("rsid", "a1", "a0", "beta", "p", "beta_se", "n_eff", "chr", "pos") 
} else if (dim(sumstats)[2] == 10) { 
colnames(sumstats) <- c("rsid", "a1", "a0", "beta", "p", "beta_se", "n_eff", "af", "chr", "pos") } else {
print("ERROR: sumstats doesn't provide enough info!")
stop()}

# process summary stats
sumstats <- na.omit(sumstats)
sumstats <- sumstats[sumstats$beta_se != 0, ]
print(dim(sumstats))
sumstats$n_eff <- as.integer(sumstats$n_eff)

# read reference bim
snp_ref <- bigreadr::fread2(ref_bim)
snp_ref <- snp_ref[-3]
names(snp_ref) <- c("chr", "rsid", "pos", "a0", "a1")

# read target bim
snp_trg <- bigreadr::fread2(trg_bim)
snp_trg <- snp_trg[-3]
names(snp_trg) <- c("chr", "rsid", "pos", "a0", "a1")

# match summary stats, target and reference
sumstats <- snp_match(sumstats, snp_trg, join_by_pos = FALSE, match.min.prop = 0.05)[c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")]
sumstats <- snp_match(sumstats, snp_ref, join_by_pos = FALSE, match.min.prop = 0.05)[c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")]
sumstats$idx_trg <- match(sumstats$rsid, snp_trg$rsid)
sumstats$idx_ref <- match(sumstats$rsid, snp_ref$rsid)

# calculate ld scores and matrices
for (chr in 1:22) {
    ind.sst <- which(sumstats$chr == chr)
    ind.ref <- sumstats$idx_ref[ind.sst]
    ind.ld <- match(ind.ref, which(snp_ref$chr == chr))

    corr0 <- readRDS(paste0("/medpop/esp2/yruan/raw.data/ldpred2_ref/ld_1kg_", pop, "/ld_1kg_chr", chr, ".rds"))[ind.ld,ind.ld]

    if (chr == 1) {
        df_beta <- sumstats[ind.sst, c("beta", "beta_se", "n_eff", "idx_trg")]
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        df_beta <- rbind(df_beta, sumstats[ind.sst, c("beta", "beta_se", "n_eff", "idx_trg")])
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# ld score regression
## origin 
ldsc <- snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
## adjust in case beta and beta_se are all 0:
## ldsc <- snp_ldsc(ld, length(ld), chi2 = ifelse(df_beta$beta==0, 0, (df_beta$beta / df_beta$beta_se)^2), sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
print(paste0("h2_est =",h2_est))

# measure running time
ldsc_time <- Sys.time()
print(ldsc_time - start_time)
cat(" for LDSC regression \n")

# ldpred2-grid
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
grid.param.name <- expand.grid(p = paste0('p:',p_seq), h2 = paste0('h2:',c(0.7, 1, 1.4)), sparse = c('F', 'T'))
gpn = do.call("paste", c(grid.param.name, sep = "_"))

## originally using detectCores
## no.cores <- detectCores(all.tests = FALSE, logical = TRUE)
## detectCores seems to give unrealistic result sometimes. if this happens, shift to availableCores()
no.cores <- availableCores()
print(paste(no.cores, "cores are available for LDpred2 parallel running"))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = no.cores)

new<-sumstats %>% select(rsid, a1, a0)
new <- cbind(new, beta_grid)
names(new)<-c(c("SNP", "A1", "A0"), gpn)
new <- new[,colSums(is.na(new))<nrow(new)]
export(new, pst_file, quote=F)

file.remove(paste0(tmp, ".sbk"))

# measure running time
end_time <- Sys.time()
print(end_time - ldsc_time)
cat(" for LDpred2-grid \n")

