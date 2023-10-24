library(data.table)
library(dplyr)
library(stringr)
library(rio)
library(optparse)
source("./function.r")

option_list = list(
  make_option(c("-m", "--med.file"), type="character", default=NULL,
              help="the file that contains the median of PC of the fine-tuning cohort. Have header. 1st row is cohort name",
              metavar="character"),
  make_option(c("-p", "--pca.file"), type="character", default=NULL,
              help="the file that contains PC of the testing cohort. Have header. 1st row is individual ID. Other rows are top PCs in ascending order.",
              metavar="character"),
  make_option(c("--prs.list"), type="character", default=NULL,
              help="the list of prs files. In the same order as in the med.file",
              metavar="character"),
  make_option(c("-A", "--a.list"), type="character", default=NULL,
              help="the list of weight. Separate items with ,", metavar="character"),
  make_option(c("-s", "--select.col"), type="character", default=NULL,
              help="the IID and PRS column in the prs files. Separate items with , ", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="the output name prefix", metavar="character")
);


# med.file <- "/medpop/esp2/yruan/projects/pca.adjusted.effect/mgbb.pca.med.1000G.tsv"
# pca.file <- "/medpop/esp2/yruan/projects/pca.adjusted.effect/g1000.pca/ukbb.g1000_5pop.hm3.pruned.scores.txt.gz"
# prslist <- "/medpop/esp2/yruan/projects/zyu.qt.traits/data/prs.ukbb.from.mgbb/CAD.AFR.scores.txt.gz,/medpop/esp2/yruan/projects/zyu.qt.traits/data/prs.ukbb.from.mgbb/CAD.EAS.scores.txt.gz,/medpop/esp2/yruan/projects/zyu.qt.traits/data/prs.ukbb.from.mgbb/CAD.EUR.scores.txt.gz,/medpop/esp2/yruan/projects/zyu.qt.traits/data/prs.ukbb.from.mgbb/CAD.SAS.scores.txt.gz"

cat("\ncheck the input:\n")
Good.Input=TRUE

cat("\ncheck the number of PCA:\n")
Med <- data.frame(fread(med.file), header = T)
npca1 = ncol(Med) - 1
target.pca <- data.frame(fread(pca.file), header = T )
npca2 = ncol(target.pca) - 1

if(min(npca1) < 5) {
print("Fewer than 5 top PCs in med.file. May lead to inefficient genetic distance calculation! ")}
if(min(npca2) < 5) {
print("Fewer than 5 top PCs in prs.file. May lead to inefficient genetic distance calculation! ")}
npca = min(npca1, npca2)
if ( npca < 3) {
print("Too few top PCs. Please check the input")
Good.Input=F}


cat("\ncheck the number of fine-tuning cohorts\n")
med <- Med[, 2:ncol(Med)]
base.list <- as.matrix(Med)[, 1]
N <- length(base.list)
print(paste0("The combine PRS based on ",N, " cohorts."))
print(base.list)
if (N<=1) {
cat("\nPlease provide 2 or more PRS to combine\n")
Good.Input=F}

if (is.null(prs.list)) {
cat("\nPlease provide list of PRS to combine\n")
Good.Input=F } else {
prs.list <- unlist(strsplit(prs.list, ","))}

if (is.null(A)) { A <- rep(1, N)} else {
A <- as.numeric(unlist(strsplit(a.list, ",")))}
if (any(is.na(A)) | any(A < 0) {
cat("Error: Unexpected value in the a.list. Please check input \n")
Good.Input=F}

if (length(prs.list)!=N | length(a.list)!=N) {
cat("Error: Inconsistent numbers of fine tuning cohorts. Please check input\n")
Good.Input=F}

cat("\ncheck the header to select from the PRS files\n")
if (length(select.col) == 2*n) {select.col <- select.col
} else if (length(select.col) == 2) {select.col <- rep(select.col, N)
} else {print("Please check the --selectcol. Incorrect input")
Good.Input=F}

##### Finish checking the Input #####

if (Good.Input==F) {
stop("The input is incorrect and the program will be stopped. Please check!")
} else {
cat("\nStart to calculate!")

# generate the shrinkage parameter depends on the genetic distance of fine-tuning cohorts
a <- shrinkage.vector(med)

# 
cat("\nCalculate genetic distance between fine-tuning cohorts and individuals in testing cohort\n")
colnames(target.pca) <- c("IID", paste0("PC", seq(npca1)))
mat <- target.pca[, 2:ncol(target.pca)]

dis.matrix <- data.frame(matrix(,nrow=nrow(target.pca),ncol=0))
for (i in seq(N)) {
dis <- Euclidean.matrix(mat, med[i,])
dis.matrix <- cbind(dis.matrix, dis)}
colnames(dis.matrix) <- c(base.list)

cat("\nCalculate weight to combine PRS for individuals in testing cohort\n") 
w.matrix <- data.frame(matrix(,nrow=nrow(target.pca),ncol=0))
for (i in seq(N)) {
w <- a[i]*A[i]/dis.matrix[,i]
w.matrix <- cbind(w.matrix, w)}
w.sum <- rowSums(w.matrix)

IID <- target.pca$IID
a.matrix <- cbind(IID)
for (i in seq(N)) {
a <- w.matrix[,i]/w.sum
a.matrix <- cbind(a.matrix, a)}
colnames(a.matrix) <- c("IID", paste0(base.list, ".a"))

cat("\nGather PRS to combine\n")
prs.matrix <- fread(prs.list[1], select=select.col[c(1,2)])
colnames(prs.matrix) <- c("IID", base.list[1])
for (i in seq(2, N)) {
score <- fread(prs.list[i], select=select.col[c(2*i-1, 2*i)])
colnames(score) <- c("IID", base.list[i])
prs.matrix <- merge(prs.matrix, score)}
Nprs <- nrow(prs.matrix)
print(paste("Read PRS for", Nprs, "individuals"))
Npca <- nrow(target.pca)
print(paste("Read PCA for", Npca, "individuals"))

dat <- merge(prs.matrix, target.pca, by="IID")
Ndat <- nrow(dat)
print(paste(Ndat, "individuals have PRS and PCA information"))
if (Ndat < 0.4*min(Nprs, Npca)){
cat("\nWarning: the overlap between PRS and PCA information is too small. Please check input")}
if (Ndat==0) {
stop("No overlap between PRS and PCA information")} else {
cat("\n Calculate combined PRS\n")
IID <- dat$IID
prs.res.matrix <- cbind(IID)
for (i in seq(N)) {
prs.res <- resid(glm(lm(formula = as.formula(paste(base.list[i], " ~", paste(paste0("PC", seq(10), collapse = "+")))),data = dat)))
prs.res.matrix <- cbind(prs.res.matrix, prs.res)}
colnames(prs.res.matrix) <- c("IID", paste0(base.list, ".res"))

dat <- merge(prs.res.matrix, a.matrix, by= "IID")

prs <- rep(0, nrow(dat))
for (i in seq(N)) {
print(paste(colnames(dat)[1+i], "x", colnames(dat)[1+N+i]))
prs <- prs + dat[, 1+i] * dat[, 1+N+i]
}
out.dat <- cbind(IID, prs)
colnames(out.dat) <- c("IID", "PRS")
export(out.dat, paste0(out, ".tsv.gz"), quote = F)
cat("\nPRS combined! \n  Y(^ W ^)Y  \n")}
}
