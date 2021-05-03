module load plink/1.90b6.21
module load eigensoft
module load R
module load metal
module load python2
module load ldsc/2018-08-10
DATADIR=/projectnb/bs859/data/RheumatoidArthritis/final_project
LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'

##################################################### 1 ###########################################################################################################
################################################################################################################################################################
# Perform genetic data cleaning of the NARAC GWAS data. 
# Then, perform PCA on the data to identify study outliers, and 
# create a set of PCs that can be used in association analyses.

# Examine number of variants and individuals
wc -l $DATADIR/narac_hg19.bim
# 544276 variants
wc -l $DATADIR/narac_hg19.fam
# 2062 individuals

# create summary file of allele frequencies for all variants
plink --bfile $DATADIR/narac_hg19 --freqx --allow-no-sex
# this produces plink.frqx has genetype vounts, col 5 is the A1 homozygote count

# look into plink.frqx file
head plink.frqx


#Create a filtered data set where all variants have minor allele frequencies >0.01, fewer
#than 5% missing genotypes, and the HWE p-value in controls is greater than 0.000001 and
#where all individuals have fewer than 5% missing genotypes using a single PLINK
#command.
plink --bfile $DATADIR/narac_hg19 --maf 0.01 --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed \
--out narac_hg19_filtered --allow-no-sex
#Options in effect:
#  --allow-no-sex
#  --bfile /projectnb/bs859/data/RheumatoidArthritis/final_project/narac_hg19
#  --geno 0.05          18402 variants removed due to missing genotype data (--geno).<<<<<<<<<<<<<<<<< is this too high?
#  --hwe 1e-6           --hwe: 663 variants removed due to Hardy-Weinberg exact test.
#  --maf 0.01           22907 variants removed due to minor allele threshold(s)
#  --make-bed
#  --mind 0.05          0 people removed due to missing genotype data (--mind).
#  --out narac_hg19_filtered
#502304 variants and 2062 people pass filters and QC.
#Among remaining phenotypes, 868 are cases and 1194 are controls.
"""5) How does the data set in part 3 differ from the dataset created in part 4? Are there any
advantages or disadvantages to using the default order of operations in PLINK when
creating the cleaned dataset?
The filtered dataset from Question 3 had 308,695 variants and 1364 people who passed the
filters while the filtered dataset from Question 4 had 307,031 variants and 1369 people who
passed the filters. If maximizing sample size is important, then applying the filter on individuals
last is likely to yield a bigger sample, as it did here. However, this may result in keeping
individuals whose data are not optimal. The default order is simpler to implement – one
command versus 2."""

# Perform the test of differential missingness by case status.
plink --bfile narac_hg19_filtered --test-missing --allow-no-sex --out narac_missing
awk '$5<0.0001{print $0}' narac_missing.missing> narac_missing.txt
wc narac_missing.txt
# 6503  32515 364168 narac_missing.txt
# 6503 y SNPs are significant at p<0.0001
sort -gk 5 narac_missing.txt |head -n 1
# CHR         SNP     F_MISS_A     F_MISS_U            P     >>from:  head narac_missing.missing
# 12   rs1316952       0.1002     0.002513    3.923e-30 
# 0.0977 e difference in missingness between the cases and controls for the SNP with the lowest p-value
c) do you have concerns about spurious (false) associations due to non-random
missingness by genotype and phenotype? Explain your answer.
There are a fairly large number of SNPs with significant difference in missingness
by case status. If the missingness is non-random, and is also associated with
genotype, then this could result in spurious (false) associations 


# PRUNING HW2/3
plink --bfile narac_hg19_filtered --allow-no-sex --chr 1-22 --indep-pairwise 10000 5000 0.2 --out narac_5000prune
# Pruning complete.  388732 of 490241 variants removed.
wc narac_prune.prune.in
# 101509 variants remain

# GRM
# not necessary? I don't do anything with this file
plink --bfile narac_hg19_filtered --allow-no-sex --extract narac_prune.prune.in --make-bed --make-rel square --out grm


#Run smartpca on the narac hg19 filtered  data.  Save 10 principal components.
nano narac_pca_in.par
#####Parameter file:
#genotypename: narac_hg19_filtered.bed
#snpname: narac_hg19_filtered.bim
#indivname: narac_hg19_filtered.fam
#evecoutname: narac.evec
#evaloutname: narac.eval
#altnormstyle: NO
#numoutevec: 10
#numoutlieriter: 2
#numoutlierevec: 8
#outliersigmathresh: 8
#outlieroutname: outliers.removed
smartpca -p narac_pca_in.par > narac_10_pc.out

wc outliers.removed
# 0 0 0 outliers.removed
#we know we’ve removed everyone with a PC that is
#>8 SD from the mean in the sample

nano narac_10_pc.out # scroll down 
# PCs 1,2,3,5,6?,8 and 9 all differ at p<0.05, but in order of lowest pval:
# PCs 1,2,5,3,8,9 (6 is 0 so this might be first)
## Statistical significance of differences beween populations:
#                                pop1                  pop2	chisq          p-value   |pop1|   |pop2|
#popdifference:               Control                  Case	 689.136  1.34871e-141    1194     868

#Create pairwise PC plots of PCs that are significantly different at p<= 0.001 by case status.
cp /projectnb/bs859/spr21/users/tpillars/class3/plotPCs.R .
R --vanilla --args narac.evec 1 2 10 <plotPCs.R
R --vanilla --args narac.evec 1 3 10 <plotPCs.R
R --vanilla --args narac.evec 1 5 10 <plotPCs.R
R --vanilla --args narac.evec 1 6 10 <plotPCs.R
R --vanilla --args narac.evec 1 8 10 <plotPCs.R

R --vanilla --args narac.evec 2 3 10 <plotPCs.R
R --vanilla --args narac.evec 2 5 10 <plotPCs.R
R --vanilla --args narac.evec 2 6 10 <plotPCs.R
R --vanilla --args narac.evec 2 8 10 <plotPCs.R

R --vanilla --args narac.evec 3 5 10 <plotPCs.R
R --vanilla --args narac.evec 3 6 10 <plotPCs.R
R --vanilla --args narac.evec 3 8 10 <plotPCs.R

R --vanilla --args narac.evec 5 6 10 <plotPCs.R
R --vanilla --args narac.evec 5 8 10 <plotPCs.R

R --vanilla --args narac.evec 6 8 10 <plotPCs.R

#get pc headers from old file
awk 'NR==1{print $0}' /projectnb/bs859/data/tgen/cleaned/TGEN_cleaned_pcs.txt > pc_header.txt
#split first col and save first col
awk -F '[:]' 'NR>1{ print $1 }' narac.evec > narac.evec.col1
#split first col and save all other cols
awk -F '[,:]' 'NR>1{ print $2 $3 }' narac.evec > narac.evec.colrest
# add col to begining
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' narac.evec.col1 narac.evec.colrest > narac.evec_new
# cat header onto new evec file
cat pc_header.txt narac.evec_new > narac_cleaned_pcs.txt
# retry plink command that didn't lead anywhere from line 190
plink --bfile narac_hg19_filtered --covar narac_cleaned_pcs.txt --covar-name PC1-PC10 --out checkPCs --allow-no-sex --logistic no-snp beta

cat checkPCs.assoc.logistic
TEST    NMISS       BETA         STAT            P
       PC1     2062     -33.48       -9.013    2.009e-19
       PC2     2062      37.32        12.86    7.488e-38
       PC3     2062       12.7        4.488    7.199e-06
       PC4     2062    -0.8863      -0.2737       0.7843
       PC5     2062      25.33        9.446    3.531e-21
       PC6     2062      44.54        15.77     4.82e-56
       PC7     2062      3.787        1.543       0.1228
       PC8     2062      -11.3       -4.465    8.018e-06
       PC9     2062      5.579        2.235      0.02539
      PC10     2062    -0.4734      -0.1886       0.8504
# PCs 1,2,3,5,6,7,& 8 are associates with p<0.05.


###################################################### 2 a ##########################################################################################################
################################################################################################################################################################

#HW GUIDE
#Perform a GWAS using logistic regression, adjusting for the PCs associated with case status at p<0.05 that
#you found in question 1. Use a plink command that hides the covariates so that your plink output has only
#the association results for the SNPs, and that gives the effect estimates as regression coefficients rather
#than odds ratios.
##FINAL PROJECT REQUEST:
#Perform two genome-wide association analyses for rheumatoid arthritis: one using only female
#subjects, and one using only male subjects.

# make fam files for just female and males
awk '($5==2){print$0}' narac_hg19_filtered.fam > narac_cleaned_f.fam
awk '($5==1){print$0}' narac_hg19_filtered.fam > narac_cleaned_m.fam

# create new bim/bed/fam files based on female status
plink --bed narac_hg19_filtered.bed \
--bim narac_hg19_filtered.bim \
--fam narac_hg19_filtered.fam \
--keep narac_cleaned_f.fam \
--make-bed --out narac_cleaned_f
# create new bim/bed/fam files based on male status
plink --bed narac_hg19_filtered.bed \
--bim narac_hg19_filtered.bim \
--fam narac_hg19_filtered.fam \
--keep narac_cleaned_m.fam \
--make-bed --out narac_cleaned_m
####################### gwas with logistic regression #########################################################################################################################################
# perform gwas with logistic regression, PCs with case status at p<0.05 are listed
plink --bfile narac_cleaned_f \
--logistic beta hide-covar \
--covar narac_cleaned_pcs.txt \
--covar-name PC1,PC2,PC3,PC5,PC6,PC7,PC8 --out gwas_f
# perform gwas with logistic regression, PCs with case status at p<0.05 are listed
plink --bfile narac_cleaned_m \
--logistic beta hide-covar \
--covar narac_cleaned_pcs.txt \
--covar-name PC1,PC2,PC3,PC5,PC6,PC7,PC8 --out gwas_m
####################### PRUNE #########################################################################################################################################
# prune f
plink --bfile narac_cleaned_f --indep-pairwise 10000 1000 0.2 --out prune_f # since I'm already cutting the data in half lets to 1000
# prune m
plink --bfile narac_cleaned_m --indep-pairwise 10000 1000 0.2 --out prune_m
####################### GRM #########################################################################################################################################
# make grm for f
plink --bfile narac_cleaned_f \
--exclude prune_f.prune.out --make-rel square --out grm_f
# make grm for m
plink --bfile narac_cleaned_m \
--exclude prune_m.prune.out --make-rel square --out grm_m

cp /projectnb/bs859/spr21/users/tpillars/assignment4/hw4_GMMAT.R . # rename GMMAT_f.R or GMMAT_m.R and modify
R --vanilla <GMMAT_f.R>grm_gmmat_f.log
R --vanilla <GMMAT_m.R>grm_gmmat_m.log

#How many SNPs in female GWAS have p-value < 0.0001?
awk '$9<0.0001 {print $0}' gwas_f.assoc.logistic>logistPCslowpval_f
wc logistPCslowpval_f
# 199  1791 17910 logistPCslowpval_f

#How many SNPs in male GWAS have p-value < 0.0001?
awk '$9<0.0001 {print $0}' gwas_m.assoc.logistic>logistPCslowpval_m
wc logistPCslowpval_m
# 129  1161 11610 logistPCslowpval_m

#show the plink results for the most significant SNP. Show the 3 genotypes for the SNP, and how they are
#coded in the plink logistic regression. Compute the odds ratio. Which allele increases the risk of our Disease?
## female
awk 'NR==1||$9<0.000001{print $0, "OR="exp($7)}' logistPCslowpval_f
# visually look for lowest
#       6   rs6457617   32663851    G        ADD     1467     -1.509       -11.55    7.682e-31 OR=0.221131
grep rs6457617 narac_cleaned_f.bim
# 6 	rs6457617	0	32663851	G	A
### AA is the genotype with code 0, G/A, and GG is 2
### G allele decreases risk of RA because beta is < 0. The AA genotype has highest risk for disease. 
## male
awk 'NR==1||$9<0.000001{print $0, "OR="exp($7)}' logistPCslowpval_m
#       6    rs660895   32577380    G        ADD      569      1.596         8.34    7.426e-17 OR=4.93326
grep rs660895 narac_cleaned_m.bim
# 6	rs660895	0	32577380	G	A
### AA is the genotype with code 0, G/A, and GG is 2
### G allele increases the risk of RA because beta is positive (greater than 0). The AA genotype has lowest risk for disease

## female more p-vals under 0.0001 in no covariate file
awk 'NR==1||$11<0.0001{print $0}' finproj_f.glmm.score.yescov |wc
#    191    2101   13170
awk 'NR==1||$11<0.0001{print $0}' finproj_f.glmm.score.nocov |wc
#    409    4499   28170
## female lowest pval variant from before is lower in nocov
grep rs6457617 finproj_f.glmm.score.yescov
#6	rs6457617	0	32663851	G	A	1467	0.639059	107.379	80.5045	5.24929e-33
grep rs6457617 finproj_f.glmm.score.nocov
#6	rs6457617	0	32663851	G	A	1467	0.639059	171.018	137.43	3.33965e-48

## male more p-vals under 0.0001 in no covariate file
awk 'NR==1||$11<0.0001{print $0}' finproj_m.glmm.score.yescov |wc
#    140    1540    9528
awk 'NR==1||$11<0.0001{print $0}' finproj_m.glmm.score.nocov |wc
#    200    2200   13610
## male lowest pval variant from before is lower in nocov
grep rs660895 finproj_m.glmm.score.yescov
#6	rs660895	0	32577380	G	A	569	0.668717	-55.1507	38.4615	5.95934e-19
grep rs660895 finproj_m.glmm.score.nocov
#6	rs660895	0	32577380	G	A	569	0.668717	-87.9388	59.1656	2.87269e-30

####################### PRODUCE QQ PLOTS HERE#########################################################################################################################################
 #NEED QQ PLOTS TO CHOOSE WHICH YES/NO COV files to move forward with
cp /projectnb/bs859/spr21/users/tpillars/assignment5/qqplot.R .
R --vanilla --args gwas_f.assoc.logistic "gwas_f_asso_logi" ADD <qqplot.R >gwas_f_asso_logi.qq.log
R --vanilla --args gwas_m.assoc.logistic "gwas_m_asso_logi" ADD <qqplot.R >gwas_m_asso_logi.qq.log
cp /projectnb/bs859/spr21/users/tpillars/assignment5/qqplot.R . # modify file to read pval instead of p
mv qqplot.R qqplot_pvalname.R
R --vanilla --args finproj_f.glmm.score.nocov "gwas_f_asso_logi_glmm_nocov" ADD <qqplot_pvalname.R >gwas_f_asso_logi_glmm_nocov.qq.log
R --vanilla --args finproj_f.glmm.score.yescov "gwas_f_asso_logi_glmm_yescov" ADD <qqplot_pvalname.R >gwas_f_asso_logi_glmm_yescov.qq.log
R --vanilla --args finproj_m.glmm.score.nocov "gwas_m_asso_logi_glmm_nocov" ADD <qqplot_pvalname.R >gwas_m_asso_logi_glmm_nocov.qq.log
R --vanilla --args finproj_m.glmm.score.yescov "gwas_m_asso_logi_glmm_yescov" ADD <qqplot_pvalname.R >gwas_m_asso_logi_glmm_yescov.qq.log
# Overall the lambda's are all near 1, the male glmm lambda are basically the same, the female glmm nocov performs better.
# When we look at the most statistically significant variant the logistic glmm_no cov options have lower pvalues.
# Will use the logistic GLMM no covariates to do the METAL (meta) analysis

#######################  M E T A L  #########################################################################################################################################
# begin meta analysis here
# create metal.txt and fill like so:

metal metal.txt > metal.log
cat metal.log
cat METAANALYSIS1.TBL.info

# need  chr:pos:ref:alt format
awk '{print $1":"$4":"$5":"$6}' finproj_f.glmm.score.nocov > new_col_fglmm.txt
awk '{print $1":"$4":"$5":"$6}' finproj_m.glmm.score.nocov > new_col_mglmm.txt
# add col to begining
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' finproj_f.glmm.score.nocov new_col_fglmm.txt > finproj_f1_newcol.glmm.score.nocov
awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],$0}' finproj_m.glmm.score.nocov new_col_mglmm.txt > finproj_m1_newcol.glmm.score.nocov
# rename the new col to NEWCOL
awk 'NR==1 {gsub("CHR:POS:A1:A2", "NEWCOL", $0); quit};1' finproj_f1_newcol.glmm.score.nocov > finproj_f_newcol.glmm.score.nocov
awk 'NR==1 {gsub("CHR:POS:A1:A2", "NEWCOL", $0); quit};1' finproj_m1_newcol.glmm.score.nocov > finproj_m_newcol.glmm.score.nocov

#run metal with reformatted files, use _reformat in output, and use last col as MARKER
metal meta_reformat.txt > metal_reformat.log
#######################  P L O T   M E T A  #########################################################################################################################################
# already have qqplot.R need gwaplot.R
cp /projectnb/bs859/spr21/users/tpillars/class4/gwaplot.R .
# get info from the table needed for plots
awk 'NR==1{print "CHR","BP","P"}; NR>1{split($1,a,":"); print a[1],a[2],$10}' METAANALYSIS1.TBL > toplot.txt
# qqplot for men and women
R --vanilla --args toplot.txt crude ADD < qqplot.R > qqplot_crude.out

# get info from previous text file to plot manhatten
awk '(NR==1||($3<0.05)){print $0}' toplot.txt > toplotman_p05.txt
awk '(NR==1||($3<0.2)){print $0}' toplot.txt > toplotman_p2.txt
awk '(NR==1||($1!=6&&$3<0.005)){print $0}' toplot.txt > toplotman_no_chr6_p005.txt
awk '(NR==1||($1!=6&&$3<0.0005)){print $0}' toplot.txt > toplotman_no_chr6_p0005.txt
awk '(NR==1||($1!=6&&$3<0.00005)){print $0}' toplot.txt > toplotman_no_chr6_p00005.txt
awk '(NR==1||($1==1&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr1_p00005.txt
awk '(NR==1||($1==2&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr2_p00005.txt
awk '(NR==1||($1==3&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr3_p00005.txt
awk '(NR==1||($1==4&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr4_p00005.txt
awk '(NR==1||($1==5&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr5_p00005.txt
awk '(NR==1||($1==7&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr7_p00005.txt
awk '(NR==1||($1==8&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr8_p00005.txt
awk '(NR==1||($1==9&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr9_p00005.txt
awk '(NR==1||($1==10&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr10_p00005.txt
awk '(NR==1||($1==11&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr11_p00005.txt
awk '(NR==1||($1==12&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr12_p00005.txt
awk '(NR==1||($1==13&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr13_p00005.txt
awk '(NR==1||($1==14&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr14_p00005.txt
awk '(NR==1||($1==15&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr15_p00005.txt
awk '(NR==1||($1==16&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr16_p00005.txt
awk '(NR==1||($1==17&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr17_p00005.txt
awk '(NR==1||($1==18&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr18_p00005.txt
awk '(NR==1||($1==19&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr19_p00005.txt
awk '(NR==1||($1==20&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr20_p00005.txt
awk '(NR==1||($1==21&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr21_p00005.txt
awk '(NR==1||($1==22&&$3<0.00005)){print $0}' toplot.txt > toplotman_chr22_p00005.txt

[tpillars@scc2 final_project]$ cat toplotman_chr1_p00005.txt
CHR BP P
1 152391567 4.963e-05
[tpillars@scc2 final_project]$ cat toplotman_chr2_p00005.txt
CHR BP P
2 207981241 2.947e-05
2 36220105 4.41e-05
2 172518151 2.354e-05
2 101749923 4.884e-05
[tpillars@scc2 final_project]$ cat toplotman_chr3_p00005.txt
CHR BP P
3 58811057 1.702e-05
3 172562534 8.079e-06
[tpillars@scc2 final_project]$ cat toplotman_chr4_p00005.txt
CHR BP P
4 157402522 2.962e-06
4 157399244 1.13e-05
4 157399443 1.835e-06
[tpillars@scc2 final_project]$ cat toplotman_chr5_p00005.txt
CHR BP P
5 133047775 1.993e-07
5 9551198 1.078e-05
5 180633194 3.57e-05
[tpillars@scc2 final_project]$ cat toplotman_chr7_p00005.txt
CHR BP P
7 147195865 3.214e-05
7 147192154 2.885e-05
7 89382342 4.259e-05
[tpillars@scc2 final_project]$ cat toplotman_chr8_p00005.txt
CHR BP P
8 14082230 4.482e-05
8 106518817 2.919e-05
8 62636162 1.055e-05
8 20358618 2.306e-06
8 11349186 4.303e-05
8 20322808 8.656e-06
8 20343900 6.487e-06
8 18778837 1.291e-05
8 20351911 6.295e-06
8 14062280 3.706e-05
8 20294939 2.708e-05
8 20774324 2.645e-05
8 20353398 5.832e-06
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 14082230
# CHR	SNP	cM	POS	A1	A2	N	AF	SCORE	VAR	PVAL
8	rs981905	0	14082230	A	G	1488	0.609543	39.6699	136.258	0.000677702
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 106518817
8	rs7016378	0	106518817	C	A	1482	0.812753	27.7068	94.1726	0.00430203
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 62636162
8	rs10094729	0	62636162	G	A	1490	0.858725	26.2137	68.6577	0.0015582
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20358618
8	rs9785133	0	20358618	A	G	1466	0.888131	29.966	58.9201	9.46597e-05
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 11349186
8	rs13277113	0	11349186	A	G	1492	0.739611	-38.1734	116.297	0.000400454
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20322808
8	rs17092514	0	20322808	A	G	1492	0.887735	28.9804	59.988	0.000182755
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20343900
8	rs10106243	0	20343900	G	A	1493	0.8858	29.7557	60.8627	0.000136668
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 18778837
8	rs1038848	0	18778837	A	G	1472	0.598505	53.7927	152.137	1.29352e-05
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20351911
8	rs11204117	0	20351911	G	A	1491	0.882294	30.9267	61.9505	8.52081e-05
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 14062280
8	rs9325690	0	14062280	A	G	1491	0.609658	41.1051	136.319	0.000430564
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20294939
8	rs999746	0	20294939	A	G	1447	0.959917	14.9512	22.8554	0.00176371
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20774324
8	rs17423243	0	20774324	A	C	1449	0.969634	12.9278	16.442	0.00143147
awk '(NR==1||($1==8)){print $0}' finproj_f.glmm.score.nocov | grep 20353398
8	rs4921720	0	20353398	A	G	1491	0.882629	30.734	61.9108	9.38224e-05

[tpillars@scc2 final_project]$ cat toplotman_chr9_p00005.txt
CHR BP P
9 123706382 2.288e-07
9 123643855 2.524e-06
9 123640500 7.273e-07
9 7104617 4.549e-05
9 123701990 1.218e-06
9 123652898 6.011e-07
9 123690239 4.698e-07
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123706382
9	rs2900180	0	123706382	A	G	1493	0.659076	-57.4164	134.546	7.42314e-07
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123643855
9	rs10985073	0	123643855	A	G	1493	0.577026	-52.8022	146.106	1.25187e-05
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123640500
9	rs1953126	0	123640500	A	G	1493	0.660415	-54.9366	134.491	2.16768e-06
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 7104617
9	rs10122120	0	7104617	A	G	1493	0.56564	-39.1129	145.995	0.00120766
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123701990
9	rs10760130	0	123701990	G	A	1493	0.576021	-54.8214	146.836	6.06502e-06
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123652898
9	rs881375	0	123652898	A	G	1489	0.659167	-54.4291	133.608	2.49107e-06
awk '(NR==1||($1==9)){print $0}' finproj_f.glmm.score.nocov | grep 123690239
9	rs3761847	0	123690239	G	A	1489	0.584621	-55.558	145.457	4.09346e-06

[tpillars@scc2 final_project]$ cat toplotman_chr10_p00005.txt
CHR BP P
10 50097819 1.662e-06
10 50119054 1.332e-05
10 50028106 1.208e-05
[tpillars@scc2 final_project]$ cat toplotman_chr11_p00005.txt
CHR BP P
[tpillars@scc2 final_project]$ cat toplotman_chr12_p00005.txt
CHR BP P
12 55916900 4.379e-05
12 54172723 1.59e-05
12 54945065 1.603e-05
[tpillars@scc2 final_project]$ cat toplotman_chr13_p00005.txt
CHR BP P
[tpillars@scc2 final_project]$ cat toplotman_chr14_p00005.txt
CHR BP P
14 105168428 1.859e-05
14 93125282 8.813e-06
14 97087097 1.711e-05
[tpillars@scc2 final_project]$ cat toplotman_chr15_p00005.txt
CHR BP P
15 73637772 1.362e-05
[tpillars@scc2 final_project]$ cat toplotman_chr16_p00005.txt
CHR BP P
16 10059007 1.875e-05
16 26888353 1.119e-05
[tpillars@scc2 final_project]$ cat toplotman_chr17_p00005.txt
CHR BP P
[tpillars@scc2 final_project]$ cat toplotman_chr18_p00005.txt
CHR BP P
18 28619795 5.335e-06
18 28628004 7.315e-06
[tpillars@scc2 final_project]$ cat toplotman_chr19_p00005.txt
CHR BP P
19 32159870 2.117e-05
19 32166846 6.29e-06
[tpillars@scc2 final_project]$ cat toplotman_chr20_p00005.txt
CHR BP P
20 58393002 2.735e-07
20 57704379 4.91e-05
[tpillars@scc2 final_project]$ cat toplotman_chr21_p00005.txt
CHR BP P
[tpillars@scc2 final_project]$ cat toplotman_chr22_p00005.txt
CHR BP P
# plot manhatten at two different pval thresholds
R --vanilla --args  toplotman_p05.txt "p<05" first_manhattan < gwaplot.R  > p05_gwa.log ### more impressive and very consise 
R --vanilla --args  toplotman_p2.txt "p<.2" second_manhattan < gwaplot.R  > p2_gwa.log
R --vanilla --args  toplotman_no_chr6_p005.txt "p<0.005 & no chr6" third_manhattan < gwaplot.R > no_chr6_p005.log
#R --vanilla --args  toplotman_no_chr6_p0005.txt "p<0.0005 & no chr6" third_manhattan < gwaplot.R > no_chr6_p0005.log
#R --vanilla --args  toplotman_no_chr6_p00005.txt "p<0.00005 & no chr6" third_manhattan < gwaplot.R > no_chr6_p00005.log

# Extract the results for the top variants
awk '(NR==1||($10<10e-08)){print $0}' METAANALYSIS1.TBL > bestpval.txt
# ALL MARKERS WITH HIGHLY SIGNIFICANT PVALUES HAVE THE SAME DIRECITON OF EFFECT
# there doesn't seem to be evidence of heterogeneity, all HetPVals are very high

#######################  LD score regression  #########################################################################################################################################
# Use LD score regression to 
#1) estimate the heritability of RA, and 
#2) compare the heritability in the Asian and European populations. 
#Describe your methods and present and explain your results.

RA_GWASmeta_European_v2.txt
RA_GWASmeta_Asian_v2.txt.gz
RA_GWASmeta_TransEthnic_v2.txt.gz
LDSCORES_DIR='/projectnb/bs859/data/ldscore_files'
zcat /projectnb/bs859/data/RheumatoidArthritis/final_project/RA_GWASmeta_Asian_v2.txt.gz | head
###   4873 cases, 17642 controls.
head /projectnb/bs859/data/RheumatoidArthritis/final_project/RA_GWASmeta_European_v2.txt
###  14,361 cases and 43,923 controls.
zcat /projectnb/bs859/data/RheumatoidArthritis/final_project/RA_GWASmeta_TransEthnic_v2.txt.gz | head
### combination

# Format
#SNPID   Chr     Position_hg19  A1      A2      OR_A1  OR_95CIlow     OR_95CIup      P-val
#rs61770171	1	750138	A	G	0.99	0.90	1.08	0.82


#######################  EUROPEAN POPULATION LC SCORE REGRESSION  #######################  

# First we need to get summary statistics in a gzip format
"""munge to reformat and check for problems>>>>>>>>class 12 slide about pre-processing"""
# making sure mean is around 1 and it is at 1.00358
awk '{ total += $6 } END { print total/NR }' $DATADIR/RA_GWASmeta_European_v2.txt
munge_sumstats.py \
--sumstats $DATADIR/RA_GWASmeta_European_v2.txt \
--snp SNPID \
--N-cas 14361 \
--N-con 43923 \
--a1 A1 \
--a2 A2 \
--signed-sumstats OR,1 \
--out munge_Euro_out

# no MAF to remove based on allele freq
# all snps accepted based on pval bounds
# Removed 1344225 variants that were not SNPs (like indels) or were strand-ambiguous. (like AT's and GC's)
# no duplicates

# output munge_Euro_out.sumstats.gz
#Metadata:
#Mean chi^2 = 1.564
#Lambda GC = 1.0
#Max chi^2 = 1143.797
#19981 Genome-wide significant SNPs (some may have been removed by filtering).
# 7403737 SNPs remain.
# zstat is made from the p-val in the original file (back transformed)
# effect is from the OR


# now we can do the actual LD score regression
ldsc.py \
--h2 munge_Euro_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out RA_euro_h2_orig_scale
#Total Observed scale h2: 0.1406 (0.0176) ## saying about 14% of the variation on the observed scale(0,1) of RA is due to the genetic effects
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076) #### no inflation, yay
#Ratio < 0 (usually indicates GC correction).

# want to convert to liability scale
# need pop prevelance and prevelance in sample
# assuming disease prevalence of 0.5% 

###################### CAN DO IN R OR CAN RE-RUN LDSC ##########################
############# this is where we will get the heritibility
# Euro
#14,361 cases and 43,923 controls.
# Sample prev
# 14361/ (14361+43923) = 0.246

ldsc.py \
--h2 munge_Euro_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.0075 \
--samp-prev 0.246 \
--out RA_euro_h2_NEW_scale.0075
# try pop-prev at 0.003, 0.005, 0.0075, 0.01, 0.015, 0.1

# --pop-prev 0.003
#Total Liability scale h2: 0.081 (0.0101)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.005
#Total Liability scale h2: 0.0897 (0.0112)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.0075
# Total Liability scale h2: 0.0979 (0.0122)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.01
#Total Liability scale h2: 0.1046 (0.0131)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.015
#Total Liability scale h2: 0.1154 (0.0144)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.1
#Total Liability scale h2: 0.1993 (0.0249)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1168
#Intercept: 0.9511 (0.0076)
#Ratio < 0 (usually indicates GC correction).

#######################  ASIAN POPULATION LC SCORE REGRESSION  ###########

###   4873 cases, 17642 controls.

# First we need to get summary statistics in a gzip format

# making sure mean is around 1 and it is at 1.00358
# this doesn't work >>> zcat awk '{ total += $6 } END { print total/NR }' $DATADIR/RA_GWASmeta_Asian_v2.txt.gz | head

munge_sumstats.py \
--sumstats $DATADIR/RA_GWASmeta_Asian_v2.txt.gz \
--snp SNPID \
--N-cas 4873 \
--N-con 17642 \
--a1 A1 \
--a2 A2 \
--signed-sumstats OR_A1,1 \
--out munge_Asia_out

#Removed 0 SNPs with missing values.
#Removed 0 SNPs with INFO <= 0.9.
#Removed 0 SNPs with MAF <= 0.01.
#Removed 0 SNPs with out-of-bounds p-values.
#Removed 1026158 variants that were not SNPs or were strand-ambiguous.
#5593713 SNPs remain.

# output munge_Asia_out.sumstats.gz.
#Metadata:
#Mean chi^2 = 1.212
#Lambda GC = 1.0
#Max chi^2 = 675.584
#9741 Genome-wide significant SNPs (some may have been removed by filtering).


# now we can do the actual LD score regression

#ldsc.py \
#--h2 munge_Asia_out.sumstats.gz \
#--ref-ld-chr $LDSCORES_DIR/1000G_Phase3_baselineLD_ldscores/ \ * error 
#--w-ld-chr $LDSCORES_DIR/1000G_Phase3_baselineLD_ldscores/ \
#--out RA_asia_h2_orig_scale

ldsc.py \
--h2 munge_Asia_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out RA_asia_h2_orig_scale

Total Observed scale h2: 0.1363 (0.026)
Lambda GC: 1.0466
Mean Chi^2: 1.0431
Intercept: 0.9802 (0.0094)
Ratio < 0 (usually indicates GC correction).

# Sample prev
# 4873/ (4873+17642) = 0.216

ldsc.py \
--h2 munge_Asia_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.01 \
--samp-prev 0.216 \
--out RA_asia_h2_NEW_scale.01
# try pop-prev at 0.005, 0.0075, 0.01

#--pop-prev 0.005
#Total Liability scale h2: 0.0953 (0.0182)
#Lambda GC: 1.0466
#Mean Chi^2: 1.0431
#Intercept: 0.9802 (0.0094)
#Ratio < 0 (usually indicates GC correction).

#--pop-prev 0.0075
#Total Liability scale h2: 0.104 (0.0199)
#Lambda GC: 1.0466
#Mean Chi^2: 1.0431
#Intercept: 0.9802 (0.0094)
#Ratio < 0 (usually indicates GC correction).

#--pop-prev 0.01
#Total Liability scale h2: 0.1111 (0.0212)
#Lambda GC: 1.0466
#Mean Chi^2: 1.0431
#Intercept: 0.9802 (0.0094)
#Ratio < 0 (usually indicates GC correction).

#######################  TRANS-ETHNIC POPULATION LC SCORE REGRESSION  #######################  

zcat /projectnb/bs859/data/RheumatoidArthritis/final_project/RA_GWASmeta_TransEthnic_v2.txt.gz | head
### combination 19,234cases / 80,799 = 0.238

munge_sumstats.py \
--sumstats $DATADIR/RA_GWASmeta_TransEthnic_v2.txt.gz \
--snp SNPID \
--N-cas 19234 \
--N-con 61565 \
--a1 A1 \
--a2 A2 \
--signed-sumstats OR_A1,1 \
--out munge_TransEthnic_out
#Read 9739303 SNPs from --sumstats file.
#Removed 0 SNPs with missing values.
#Removed 0 SNPs with INFO <= 0.9.
#Removed 0 SNPs with MAF <= 0.01.
#Removed 0 SNPs with out-of-bounds p-values.
#Removed 1498996 variants that were not SNPs or were strand-ambiguous.
#8240307 SNPs remain.

ldsc.py \
--h2 munge_TransEthnic_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--out RA_trans_ethnic_h2_orig_scale

ldsc.py \
--h2 munge_TransEthnic_out.sumstats.gz \
--ref-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--w-ld-chr $LDSCORES_DIR/eur_w_ld_chr/ \
--pop-prev 0.01 \
--samp-prev 0.238 \
--out RA_trans_ethnic_h2_NEW_scale.01
# try pop-prev at 0.003, 0.005, 0.0075, 0.01,

# --pop-prev 0.003
#Total Liability scale h2: 0.0702 (0.007)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1484
#Intercept: 0.9556 (0.0084)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.005
#Total Liability scale h2: 0.0778 (0.0078)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1484
#Intercept: 0.9556 (0.0084)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.0075
#Total Liability scale h2: 0.0849 (0.0085)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1484
#Intercept: 0.9556 (0.0084)
#Ratio < 0 (usually indicates GC correction).

# --pop-prev 0.01
#Total Liability scale h2: 0.0907 (0.0091)
#Lambda GC: 1.0466
#Mean Chi^2: 1.1484
#Intercept: 0.9556 (0.0084)
#Ratio < 0 (usually indicates GC correction).





