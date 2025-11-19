import pandas as pd

df = pd.read_csv("phenotypes.tsv",sep='\t')

# generate metadata for sex
sex = df[['IID','sex']]
sex['FID'] = 0
cleaned_sex = sex[['FID','IID','sex']]
cleaned_sex.to_csv("out/sex.txt",sep='\t',index=False)

# generate metadata for case/control phenotype
pheno = df[['IID','pancreas_cancer']]
pheno['FID'] = 0
cleaned_pheno = pheno[['FID','IID','pancreas_cancer']]
# plink2 recognizes 1 as case and 2 as control; reformat the case/control to match
cleaned_pheno.loc[df['pancreas_cancer'] == 1, 'pancreas_cancer'] = 2
cleaned_pheno.loc[df['pancreas_cancer'] == 0, 'pancreas_cancer'] = 1
cleaned_pheno.to_csv("out/pheno.txt",sep='\t',index=False)

# generate covariate metadata
covar = df[['IID','sex','age','ancestry','smoking_history',"bmi",'PC1','PC2','PC3','PC4','PC5']]
covar['FID'] = 0
cleaned_covar = covar[['FID','IID','sex','age','ancestry','smoking_history',"bmi",'PC1','PC2','PC3','PC4','PC5']]
# fill any missing values as "NA" to match plink2's convenction for addressing missing values
cleaned_covar = cleaned_covar.fillna("NA")
cleaned_covar.to_csv("out/covar.txt",sep='\t',index=False)
