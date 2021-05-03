     Show Dotfiles Show Owner/Mode
/projectnb/bs859/spr21/users/tpillars/final_project/
library(GMMAT)
pheno<-read.table("/projectnb/bs859/spr21/users/tpillars/final_project/narac_cleaned_m.fam",header=F)
colnames(pheno)<-c("FID","IID","fa","mo","sex","case")
pcs<-read.table("/projectnb/bs859/spr21/users/tpillars/final_project/narac_cleaned_pcs.txt",header=T,as.is=T)


# merge the PC data with the fam file data.
#We need to make sure we keep all the individuals in the original phenotype file (fam file) (all.x=TRUE),
#and in the same order, because this is how plink and GMMAT know the order in the genotype file (bed file).
#
#here, I create a new variable, "order" with the original order of the phenotype(fam) file.
#I can then use it to sort the phenotype file after merging with the PC file:
pheno$order<-1:length(pheno$IID)

#Now we merge with the pc data, making sure to keep ALL the individuals in the fam file:
pheno1<-merge(pheno,pcs,by=c("FID","IID"),all.x=TRUE)
#and re-order the file back to the original order of the file prior to the merge:
pheno1<-pheno1[order(pheno1$order),]

##Read in the GRM:
grm<-as.matrix(read.table("grm_m.rel",header=F))

#grm must also have IDs in the same order as the fam file.  Read in the grm id file:

grm.ids<-read.table("grm_m.rel.id",header=F)

dimnames(grm)[[1]]<-dimnames(grm)[[2]]<-grm.ids[,2]

#compare the grm ids to the ids in pheno1 (our fam file), row by row.  If the
# pheno1 ids equal the grm.ids in all rows,  you will get all "TRUE" values here:

table(pheno1[,1:2]==grm.ids)

## these two commands create the Null models (no SNPs) for the score tests:

model1.0<-glmmkin(case-1~1,data=pheno1,id="IID",kins=grm,family=binomial("logit")) #no covariate
model2.0<-glmmkin(case-1~PC1+PC2+PC3+PC5+PC6+PC7+PC8,data=pheno1,id="IID",kins=grm,family=binomial("logit")) # PC1 and PC2 covariates

## these two commands perform the score test for model 1 and model 2 for all of the SNPs in the plink genotype file wgas3.bed.
## the model results are put into the file specified by "outfile="
glmm.score(model1.0,infile="/projectnb/bs859/spr21/users/tpillars/final_project/narac_cleaned_m",outfile="/projectnb/bs859/spr21/users/tpillars/final_project/finproj_m.glmm.score.nocov")



glmm.score(model2.0,infile="/projectnb/bs859/spr21/users/tpillars/final_project/narac_cleaned_m",outfile="/projectnb/bs859/spr21/users/tpillars/final_project/finproj_m.glmm.score.yescov")
