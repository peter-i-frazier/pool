# Author: Tom Fei
# Created: 06.10.2015
# Test GenOneRandomPeptide

rm(list = ls())
set.seed(1)

root <- "~/repos/peptide-catalysis/src/"
source(paste0(root, "core/dependency.R"))
ResolveDependency(root)

class.vec <- QuerySql("sfp_AcpS", "SELECT * FROM reduced_AA")
class.vec <- unique(class.vec[, 'class'])
print(class.vec)

GEN.ITER <- 1
P <- c()
for (i in 1:GEN.ITER) {
  peptide <- GenOneRandomPeptide(minL, maxL, minR, maxR,
                                 MAXL, MAXR, class.vec)
  P <- rbind(P, peptide)
}

print(P)
