# create a plink2 object from the vcf file
plink2 --vcf out/filt_candidate_EOPC_variants.vcf --make-pgen --out out/temp_out

# add sex and phenotype information
# per the PLINK2 guidelines (https://www.cog-genomics.org/plink/2.0/assoc#glm), it is recommended to use a minor allele count of 20 (i.e. keep variants which have a minor allele count >= 20) because the GLM models are not well calibrated when the minor allele count is very small
plink2 --pfile out/temp_out --update-sex out/sex.txt --make-pgen --pheno out/pheno.txt --out out/withPhenoSex_candidate_EOPC_variants --mac 20

# run a logistic regression, which is recommended when dealing with binary phenotypes (case/control)
# include sex, age, smoking history, BMI and the first 5 PCs as covariates
# ancestry was excluded as a covariate because it correlated strongly with the PCs, rendering it redundant
# Per Price et al. (2006) Principal components analysis corrects for stratification in genome-wide association studies, including the top PCs as covariates is necessary to correct for population stratification
# since early-onset pancreatic cancer is age-related, it is necessary to include age as a covariate in the analysis
# furthermore, including sex helps prevent potential sex-related bias, while smoking history and BMI account for lifestyle-related confounders
plink2 \
  --pfile out/withPhenoSex_candidate_EOPC_variants \
  --glm \
  --covar out/covar.txt --covar-name sex age smoking_history bmi PC1 PC2 PC3 PC4 PC5 \
  --covar-variance-standardize \
  --out out/EOPC_assoc \
  --ci 0.95