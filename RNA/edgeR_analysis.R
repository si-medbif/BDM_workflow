library(edgeR)

infile <- ""
outfile <- ""

x <- read.delim(infile,row.names="geneid")
group <- factor(c())
block <- factor(c())

y <- DGEList(counts=x,group=group)
#There is no purpose in analysing genes that are not expressed in either experimental condition, so genes are first filtered on expression levels.
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
# The TMM normalization is applied to account for the compositional biases:
y <- calcNormFactors(y)
# Create design matrix
design <- model.matrix(~block+group)
y <- estimateDisp(y,design, robust=TRUE)
# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design, robust = TRUE)
# Check for batch effect
qlf <- glmQLFTest(fit,coef=2:3)
topTags(qlf)
# Differential expression for treatment effect
qlf <- glmQLFTest(fit)
res = topTags(qlf,  n=nrow(qlf))
write.table(res, file=outfile, sep="\t", row.names = TRUE, quote=FALSE)
