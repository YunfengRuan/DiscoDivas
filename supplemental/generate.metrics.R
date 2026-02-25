library(data.table)
library(dplyr)
library(stringr)
library(rcompanion)
library(optparse)

# The function to calculate metrics of PRS performance
# omit the code of reading data and printing the output into files

calculate <- function(cov.formula, PHENO, dat) {
line <- data.table(BETA=0.0, SE=1.0, P=1.0, Adj.R2=0, N=0)
try({Full<-glm(formula = as.formula(paste(PHENO, " ~ PRS.SUM + ", cov.formula)), family = reg.type, data = dat)
Null<-glm(formula = as.formula(paste(PHENO, " ~ ", cov.formula)), family = reg.type, data = dat)
coeff <- summary(Full)$coefficients
BETA <- coeff[2,1]
SE <- coeff[2,2]
P <- coeff[2,4]
if (reg.type == "binomial" ) {
Adj.R2<-nagelkerke(Full, null=Null)$Pseudo.R.squared.for.model.vs.null[3]
} else if (reg.type == "gaussian") {
SSE.null<-sum(Null$residuals**2)
SSE.full<-sum(Full$residuals**2)
Adj.R2 <-(SSE.null-SSE.full)/SSE.null
}
line <- data.table(BETA, SE, P, Adj.R2, N=dim(dat)[1])})
print(line)
return(line)}
