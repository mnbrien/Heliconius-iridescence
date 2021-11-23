# QTL mapping code  used for Brien et al.
# M.N.Brien 2020
library(qtl)
library(qtl2)
library(qtl2convert)
library(ggplot2)
library(dplyr)

# H. erato data files
a<-read.cross(format="csv", file="erato_f2.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="erato_bc.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"), crosstype="bc", estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="erato_all.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="erato_saxs.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
# H. melpomene data files
a<-read.cross(format="csv", file="melp_f2.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="melp_ec70.csv", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="melp_all.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)
a<-read.cross(format="csv", file="melp_saxs.csv", genotypes = c("AA", "AB","BB"), alleles = c("A", "B"),  estimate.map = FALSE, convertXdata = T)

################################################

#### Genome scans ####

# calculate genotype probabilities
a<- calc.genoprob(a, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane", stepwidth="fixed")

# extract covariates - for f2 families
cross<-as.numeric(pull.pheno(a, "family"))
sex<-as.numeric(pull.pheno(a, "sex"))
crossX<-cbind(cross, sex, sex*cross)
# run scans and calculate significance level
a_scan1<- scanone(a, pheno.col=4, model="normal", method="hk", addcovar=crossX, intcovar=cross)
plot(a_scan1)
perm.a <-scanone(a, pheno.col=4, n.perm=1000, perm.Xsp=T, addcovar=crossX, intcovar=cross, method="hk")
summary(perm.a) 

# scan for single families (e.g. EC70. covariate for family not needed here)
b_scan1<- scanone(a, pheno.col=4, model="normal", method="hk", addcovar=sex)
plot(b_scan1)
perm.b <-scanone(a, pheno.col=4, n.perm=1000, perm.Xsp=T, addcovar=sex, method="hk")
summary(perm.b)



# Combining multiple scans
plot(a_scan1+b_scan1)
# manually add lod scores
all_df <- merge(data.frame(a_scan1), data.frame(b_scan1), by="row.names")
all_df <- all_df %>% mutate(lod = lod.x +lod.y)
all_scan <- select(all_df, c("Row.names", "chr.x", "pos.x", "lod")) %>% rename(chr=chr.x, pos=pos.x)
row.names(all_scan)<-all_scan$Row.names
all_scan<-select(all_scan, c("chr", "pos", "lod"))
all_scan <-all_scan[order(all_scan$chr, all_scan$pos),]
class(all_scan)<-c("scanone", "data.frame")

#to plot combined scan in rqtl2
all_qtl1 <- scan_qtl_to_qtl2(all_scan)
plot(all_qtl1$scan1, all_qtl1$map, ylim=c(0,8))
add_threshold(all_qtl1$map, thresholdA = 5.02, thresholdX = 5.05, lty=2)

###############################################################

## significant markers
summary(a_scan1, perms=perm.a, lodcolumn=1, alpha=0.1, pvalues=TRUE) 
# intervals
bayesint(a_scan1, chr=3, prob=0.95)

# calculate % variance explained
sim <- sim.geno(a, n.draws=128, step=2, err=0.001)
qtl<- makeqtl(sim, chr=3, pos=17.97) 
out.fq <- fitqtl(sim, pheno.col=4, qtl=qtl, formula=y~Q1, forceXcovar =T) 
summary(out.fq)


# effect plots in qtl2 
a2<-convert2cross2(a)
map_a<- insert_pseudomarkers(a2$gmap, step=1)
a_gp<-calc_genoprob(a2, a2$gmap, error_prob = 0.001)
g<- maxmarg(a_gp, a2$gmap, chr=3, pos=16.14, return_char=T)
plot_pxg(g, a2$pheno[,"br"], ylab="BR", SEmult=2, sort=F, swap_axes = F) 


